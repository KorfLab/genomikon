/*****************************************************************************\
 geney.c
\*****************************************************************************/

#include "genomikon.h"

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

static gkn_len getlen(const char *filename, int size) {
	gkn_pipe io = gkn_pipe_open(filename, "r");
	gkn_len len = gkn_len_read(io, size);
	gkn_pipe_close(io);
	return len;
}

static char *usage = "\
geney - gene scoring demo\n\n\
usage: geney <fasta file> <gff file> [options]\n\
options:\n\
  -dpwm <file> donor model position weight matrix\n\
  -apwm <file> acceptor model position weight matrix\n\
  -emm  <file> exon Markov model\n\
  -imm  <file> intron Markov model (requires -dpwm & -apwm)\n\
  -elen <file> exon length model\n\
  -ilen <file> intron length model";

int main(int argc, char **argv) {
	char  *ffile = NULL; // path to fasta file
	char  *gfile = NULL; // path to gff file
	char   buff[64];
	gkn_pwm dpwm  = NULL;
	gkn_pwm apwm  = NULL;
	gkn_mm  emm   = NULL;
	gkn_mm  imm   = NULL;
	gkn_len ilen  = NULL;
	gkn_len elen  = NULL;

	// CLI
	gkn_set_program_name(argv[0]);
	gkn_register_option("-dpwm", 1);
	gkn_register_option("-apwm", 1);
	gkn_register_option("-emm", 1);
	gkn_register_option("-imm", 1);
	gkn_register_option("-elen", 1);
	gkn_register_option("-ilen", 1);
	gkn_parse_options(&argc, argv);
	if (argc != 3) gkn_exit("%s", usage);

	ffile = argv[1];
	gfile = argv[2];
	if (gkn_option("-dpwm")) dpwm = getpwm(gkn_option("-dpwm"));
	if (gkn_option("-apwm")) apwm = getpwm(gkn_option("-apwm"));
	if (gkn_option("-emm"))  emm  = getmm(gkn_option("-emm"));
	if (gkn_option("-imm"))  imm  = getmm(gkn_option("-imm"));
	if (gkn_option("-elen")) elen = getlen(gkn_option("-elen"), 1000);
	if (gkn_option("-ilen")) ilen = getlen(gkn_option("-ilen"), 1000);

	// construct mRNA
	gkn_pipe  io = gkn_pipe_open(ffile, "r");
	gkn_fasta ff = gkn_fasta_read(io);
	gkn_mRNA  tx = gkn_mRNA_read(gfile, ff->seq);

	// test scoring functions
	double score = 0;

	if (dpwm) {
		for (int i = 0; i < tx->introns->size; i++) {
			gkn_feat intron = tx->introns->elem[i];
			double s = gkn_pwm_score(dpwm, ff->seq, intron->beg);
			score += s;
			strncpy(buff, ff->seq + intron->beg, dpwm->size);
			buff[dpwm->size] = '\0';
			printf("donor: %f %s\n", s, buff);
		}
	}

	if (apwm) {
		for (int i = 0; i < tx->introns->size; i++) {
			gkn_feat intron = tx->introns->elem[i];
			double s = gkn_pwm_score(apwm, ff->seq, intron->end - apwm->size+1);
			score += s;
			strncpy(buff, ff->seq + intron->end - apwm->size +1, apwm->size);
			buff[apwm->size] = '\0';
			printf("acceptor: %f %s\n", s, buff);
		}
	}

	if (emm) {
		for (int i = 0; i < tx->exons->size; i++) {
			gkn_feat exon = tx->exons->elem[i];
			double s = gkn_mm_score(emm, ff->seq, exon->beg, exon->end);
			score += s;
			char *seq = gkn_feat_seq(exon);
			printf("exon %f %s\n", s, seq);
			free(seq);
		}
	}

	if (imm) {
		for (int i = 0; i < tx->introns->size; i++) {
			gkn_feat intron = tx->introns->elem[i];
			double s = gkn_mm_score(imm, ff->seq, intron->beg + dpwm->size,
				intron->end - apwm->size);
			score += s;
			char *seq = gkn_feat_seq(intron);
			printf("intron %f %s\n", s, seq);
			free(seq);
		}
	}

	if (elen) {
		for (int i = 0; i < tx->exons->size; i++) {
			gkn_feat exon = tx->exons->elem[i];
			int len = exon->end - exon->beg + 1;
			double s = gkn_len_score(elen, len);
			score += s;
			printf("exon %f %d\n", s, len);
		}
	}

	if (ilen) {
		for (int i = 0; i < tx->introns->size; i++) {
			gkn_feat intron = tx->introns->elem[i];
			int len = intron->end - intron->beg + 1;
			double s = gkn_len_score(ilen, len);
			score += s;
			printf("intron %f %d\n", s, len);
		}
	}

	printf("total score: %f\n", score);

	// clean up
	gkn_fasta_free(ff);
	gkn_pipe_close(io);

	return 0;
}
