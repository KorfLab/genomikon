#include "genomikon.h"
#include "isoform.h"

static int NOISY = 0;

static void chatter(const char* format, ...) {
	if (NOISY == 0) return;
	va_list args;
	fflush(stdout);
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
}


/***************\
 Models Section
\***************/

static gkn_pwm getpwm(const char *filename) {
	gkn_pipe io = gkn_pipe_open(filename, "r");
	gkn_pwm pwm = gkn_pwm_read(io);
	gkn_pipe_close(io);
	return pwm;
}

static gkn_mm getmm(const char *filename) {
	gkn_pipe io = gkn_pipe_open(filename, "r");
	gkn_mm   mm = gkn_mm_read(io);
	gkn_pipe_close(io);
	return mm;
}

static gkn_len getlen(const char *filename) {
	gkn_pipe io = gkn_pipe_open(filename, "r");
	gkn_len len = gkn_len_read(io);
	gkn_pipe_close(io);
	return len;
}

static double score_apwm(const gkn_pwm pwm, const gkn_mRNA tx) {
	double score = 0;
	for (int i = 0; i < tx->introns->size; i++) {
		gkn_feat f = tx->introns->elem[i];
		double s = gkn_pwm_score(pwm, f->seq, f->end -pwm->size +1);
		score += s;
	}
	return score;
}

static double score_dpwm(const gkn_pwm pwm, const gkn_mRNA tx) {
	double score = 0;
	for (int i = 0; i < tx->introns->size; i++) {
		gkn_feat f = tx->introns->elem[i];
		double s = gkn_pwm_score(pwm, f->seq, f->beg);
		score += s;
	}
	return score;
}

static double score_elen(const gkn_len model, const gkn_mRNA tx) {
	double score = 0;
	for (int i = 0; i < tx->exons->size; i++) {
		gkn_feat f = tx->exons->elem[i];
		int len = f->end - f->beg + 1;
		double s = gkn_len_score(model, len);
		score += s;
	}
	return score;
}

static double score_ilen(const gkn_len model, const gkn_mRNA tx) {
	double score = 0;
	for (int i = 0; i < tx->introns->size; i++) {
		gkn_feat f = tx->introns->elem[i];
		int len = f->end - f->beg + 1;
		double s = gkn_len_score(model, len);
		score += s;
	}
	return score;
}

static double score_emm(const gkn_mm mm, const double *emem, const gkn_mRNA tx) {
	double score = 0;
	for (int i = 0; i < tx->exons->size; i++) {
		gkn_feat f = tx->exons->elem[i];
		double s = gkn_mm_score_cache(mm, emem, f->beg, f->end);
		score += s;
	}
	return score;
}

static double score_imm(
	const gkn_mm mm,
	const double *imem, const gkn_mRNA tx,
	const gkn_pwm dpwm, const gkn_pwm apwm)
{
	double score = 0;
	for (int i = 0; i < tx->introns->size; i++) {
		gkn_feat f = tx->introns->elem[i];
		int beg = f->beg + dpwm->size;
		int end = f->end - apwm->size;
		double s = gkn_mm_score_cache(mm, imem, beg, end);
		score += s;
	}
	return score;
}

static double complexity(const gkn_vec txs) {
	gkn_fvec p = gkn_fvec_new();
	double total = 0;
	for (int i = 0; i < txs->size; i++) {
		gkn_mRNA tx = txs->elem[i];
		double w = pow(2, tx->score);
		total += w;
		gkn_fvec_push(p, w);
	}
	for (int i = 0; i < p->size; i++) p->elem[i] /= total;

	double h = 0;
	for (int i = 0; i < p->size; i++) {
		if (p->elem[i] > 0) h -= p->elem[i] * log2(p->elem[i]);
	}

	gkn_fvec_free(p);
	return h;
}

static int cmptx(const gkn_mRNA a, const gkn_mRNA b) {
	if      (a->score < b->score) return -1;
	else if (a->score > b->score) return  1;
	else                          return  0;
}

static int txscoresort(const void *a, const void *b) {
	return cmptx( *(gkn_mRNA *)b, *(gkn_mRNA *)a );
}

static char *usage = "\
isoformer - generate all possible isoforms from sequences\n\n\
usage: isoformer <fasta file> [options]\n\
options:\n\
  --min_exon   <int>   minimum exon length [25]\n\
  --min_intron <int>   minimum intron length [35]\n\
  --max_splice <int>   maximum splices [3]\n\
  --flank      <int>   genomic flank lengths [99]\n\
  --introns    <file>  use introns from GFF file\n\
  --limit      <int>   limit report to this many isoforms [20]\n\
  --apwm       <file>  use acceptor pwm\n\
  --dpwm       <file>  use donor pwm\n\
  --emm        <file>  use exon Markov model\n\
  --imm        <file>  use intron Markov model (requires -dpwm & -apwm)\n\
  --elen       <file>  use exon length model\n\
  --ilen       <file>  use intron length model\n\
  --wdpwm      <float> weight of donor pwm [1.0]\n\
  --wapwm      <float> weight of acceptor pwm [1.0]\n\
  --wemm       <float> weight of exon Markov model [1.0]\n\
  --wimm       <float> weight of intron Markov model [1.0]\n\
  --welen      <float> weight of exon length model [1.0]\n\
  --wilen      <float> weight of intron length model [1.0]\n\
  --icost      <float> cost for each intron [0.0]\n\
";

int main(int argc, char **argv) {
	char *file = NULL; // path to fasta file
	int   emin = 25;   // min exon size
	int   imin = 35;   // min intron size
	int   smax = 3;    // max splices
	int   gen  = 99;   // genomic flank (promoter, downstream)
	char *gff  = NULL; // use gff file for splice site coordinates
	int   head = 20;   // limit output to this many

	gkn_pwm apwm = NULL; // acceptor pwm
	gkn_pwm dpwm = NULL; // donor pwm
	gkn_mm  emm  = NULL; // exon Markov model
	gkn_mm  imm  = NULL; // intron Markov model
	gkn_len elen = NULL; // exon length model
	gkn_len ilen = NULL; // intron length model

	double wap = 1;
	double wdp = 1;
	double wem = 1;
	double wim = 1;
	double wel = 1;
	double wil = 1;
	double ic  = 0;

	gkn_pipe  io = NULL; // for reading fasta files
	gkn_fasta ff = NULL; // for individual fasta entries

	// CLI - setup
	gkn_set_program_name(argv[0]);
	gkn_register_option("--min_exon", 1);
	gkn_register_option("--min_intron", 1);
	gkn_register_option("--max_splice", 1);
	gkn_register_option("--flank",  1);
	gkn_register_option("--dpwm", 1);
	gkn_register_option("--apwm", 1);
	gkn_register_option("--emm", 1);
	gkn_register_option("--imm", 1);
	gkn_register_option("--elen", 1);
	gkn_register_option("--ilen", 1);
	gkn_register_option("--introns", 1);
	gkn_register_option("--limit", 1);
	gkn_register_option("--wdpwm", 1);
	gkn_register_option("--wapwm", 1);
	gkn_register_option("--wemm", 1);
	gkn_register_option("--wimm", 1);
	gkn_register_option("--welen", 1);
	gkn_register_option("--wilen", 1);
	gkn_register_option("--icost", 1);
	gkn_register_option("--noisy", 0);
	gkn_parse_options(&argc, argv);
	if (argc == 1) gkn_exit("%s", usage);

	// CLI - harvest
	file = argv[1];
	if (gkn_option("--noisy"))      NOISY = 1;
	if (gkn_option("--min_exon"))   emin  = atoi(gkn_option("--min_exon"));
	if (gkn_option("--min_intron")) imin  = atoi(gkn_option("--min_intron"));
	if (gkn_option("--max_splice")) smax  = atoi(gkn_option("--max_splice"));
	if (gkn_option("--flank"))      gen   = atoi(gkn_option("--flank"));
	if (gkn_option("--apwm"))       apwm  = getpwm(gkn_option("--apwm"));
	if (gkn_option("--dpwm"))       dpwm  = getpwm(gkn_option("--dpwm"));
	if (gkn_option("--emm"))        emm   = getmm(gkn_option("--emm"));
	if (gkn_option("--imm"))        imm   = getmm(gkn_option("--imm"));
	if (gkn_option("--elen"))       elen  = getlen(gkn_option("--elen"));
	if (gkn_option("--ilen"))       ilen  = getlen(gkn_option("--ilen"));
	if (gkn_option("--introns"))    gff   = gkn_option("--introns");
	if (gkn_option("--limit"))      head  = atoi(gkn_option("--limit"));
	if (gkn_option("--wdpwm"))      wdp   = atof(gkn_option("--wdpwm"));
	if (gkn_option("--wapwm"))      wap   = atof(gkn_option("--wapwm"));
	if (gkn_option("--wemm"))       wem   = atof(gkn_option("--wemm"));
	if (gkn_option("--wimm"))       wim   = atof(gkn_option("--wimm"));
	if (gkn_option("--welen"))      wel   = atof(gkn_option("--welen"));
	if (gkn_option("--wilen"))      wil   = atof(gkn_option("--wilen"));
	if (gkn_option("--icost"))      ic    = atof(gkn_option("--icost"));
	if (imm && !(dpwm && apwm)) gkn_exit("-imm requires -dpwm & -apwm");

	// find isoforms (just one fasta file entry)
	chatter("reading fasta file %s... ", file);
	io = gkn_pipe_open(file, "r");
	ff = gkn_fasta_read(io);
	chatter("done\n");
	chatter("finding isoforms... ");
	isozone iso = isoforms(ff->seq, emin, imin, smax, gen, gff, 0);
	if (head > iso->mRNAs->size) head = iso->mRNAs->size;
	chatter("found %d\n", iso->mRNAs->size);

	// score isoforms
	chatter("scoring isoforms... ");
	double *emem = NULL;
	double *imem = NULL;
	if (emm) emem = gkn_mm_cache(emm, ff->seq);
	if (imm) imem = gkn_mm_cache(imm, ff->seq);

	for (int i = 0; i < iso->mRNAs->size; i++) {
		gkn_mRNA tx = iso->mRNAs->elem[i];
		double score = 0;
		if (apwm) score += score_apwm(apwm, tx) * wap;
		if (dpwm) score += score_dpwm(dpwm, tx) * wdp;
		if (elen) score += score_elen(elen, tx) * wel;
		if (ilen) score += score_ilen(ilen, tx) * wil;
		if (emm)  score += score_emm(emm, emem, tx) * wem;
		if (imm)  score += score_imm(imm, imem, tx, dpwm, apwm) * wim;
		score -= tx->introns->size * ic;
		tx->score = score;
	}
	qsort(iso->mRNAs->elem, iso->mRNAs->size, sizeof(gkn_mRNA), txscoresort);
	chatter("done\n");

	// output
	chatter("writing output to stdout\n");
	gkn_vec txs = iso->mRNAs;
	printf("# name: %s\n", ff->def+1);
	printf("# length: %d\n", (int)strlen(ff->seq));
	printf("# donors: %d\n", iso->dons);
	printf("# acceptors: %d\n", iso->accs);
	printf("# trials: %d\n", iso->trials);
	printf("# isoforms: %d\n", iso->mRNAs->size);
	printf("# complexity: %.4f\n", complexity(txs));

	// calculate probability of each isoform from score
	double max_score = ((gkn_mRNA)(iso->mRNAs->elem[0]))->score;
	gkn_fvec p = gkn_fvec_new();
	double total = 0;
	for (int i = 0; i < head; i++) {
		gkn_mRNA tx = iso->mRNAs->elem[i];
		double w = pow(2, tx->score - max_score); // prevent overflow
		total += w;
		gkn_fvec_push(p, w);
	}
	for (int i = 0; i < head; i++) p->elem[i] /= total;

	// gff gene
	char chrom[64];
	strcpy(chrom, ff->def+1);
	for (int i = 0; i < strlen(chrom); i++) {
		if (isspace(chrom[i])) {
			chrom[i] = '\0';
			break;
		}
	}

	gkn_mRNA tx = iso->mRNAs->elem[0];
	gkn_feat e0 = tx->exons->elem[0];
	gkn_feat en = tx->exons->elem[tx->exons->size -1];
	int gbeg = e0->beg + 1;
	int gend = en->end + 1;
	char gene[128];
	sprintf(gene, "Gene-%s", chrom);
	printf("%s\tisoformer\tgene\t%d\t%d\t.\t+\t.\tID=%s\n\n",
		chrom, gbeg, gend, gene);

	// gff transcripts
	char tid[256];
	for (int i = 0; i < head; i++) {
		gkn_mRNA tx = iso->mRNAs->elem[i];
		sprintf(tid, "tx-%s-%d", chrom, i + 1);
		printf("%s\tisoformer\tmRNA\t%d\t%d\t%.4g\t+\t.\tID=%s;Parent=%s\n",
			chrom, tx->beg +1, tx->end +1, p->elem[i],
			tid, gene);

		// exons
		for (int j = 0; j < tx->exons->size; j++) {
			gkn_feat exon = tx->exons->elem[j];
			printf("%s\tisoformer\texon\t%d\t%d\t%.4g\t+\t.\tParent=%s\n",
				chrom, exon->beg +1, exon->end +1,
				p->elem[i], tid);
		}

		// introns
		for (int j = 0; j < tx->introns->size; j++) {
			gkn_feat intron = tx->introns->elem[j];
			printf("%s\tisoformer\tintron\t%d\t%d\t%.4g\t+\t.\tParent=%s\n",
				chrom, intron->beg +1, intron->end +1,
				p->elem[i], tid);
		}
		printf("\n");
	}

	return 0;
}
