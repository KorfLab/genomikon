
#include "genomikon.h"
#include "dmg.h"
#include "pwm.h"

struct max_motif {
	double score;
	double entropy;
	char seq[32];
};

static struct max_motif maximum_motif(const gkn_pwm o, const gkn_pwm e) {
	char nts[4] = "ACGT";
	int count[4] = {0, 0, 0, 0};
	struct max_motif mm;
	mm.score = 0;

	for (int i = 0; i < o->size; i++) {
		double maxs = -1e300;
		char maxc = '\0';
		for (int j = 0; j < 4; j++) {
			double s = o->score[i][j] - e->score[i][j];
			if (s > maxs) {
				maxs = s;
				maxc = nts[j];
			}
		}
		mm.score += maxs;
		mm.seq[i] = maxc;
	}
	mm.seq[o->size] = '\0';

	for (int i = 0; i < o->size; i++) {
		switch (mm.seq[i]) {
			case 'A': count[0]++; break;
			case 'C': count[1]++; break;
			case 'G': count[2]++; break;
			case 'T': count[3]++; break;
			default: gkn_exit("no, not possible");
		}
	}

	double h = 0;
	for (int i = 0; i < 4; i++) {
		if (count[i] != 0) {
			double p = (double)count[i]/o->size;
			h -= p * log(p);
		}
	}
	mm.entropy = h;

	return mm;
}

static char *usage = "\
usage: motifamatic <fasta> [options]\n\
options:\n\
  --len <int>    motif length [6]\n\
  --aid <int>    alphabet id [3]\n\
  --mod <int>    model id [0]\n\
  --mms <double> minimum motif score per position [1.0]\n\
  --mmh <double> minimum motif entropy [1.0]\n\
  --neg <fasta>  file of negative sequences\n\
  --P1  <float>  uppercase single letter probability [0.97]\n\
  --P2  <float>  uppercase double letter probability [0.49]\n\
  --P3  <float>  uppercase triple letter probability [0.33]\n\
  --p1  <float>  lowercase single letter probability [0.70]\n\
  --p2  <float>  lowercase double letter probability [0.40]\n\
  --p2  <float>  lowercase double letter probability [0.30]\n\
alphabets:\n\
  0 ACGT\n\
  1 ACGTN\n\
  2 ACGTacgtN\n\
  3 ACGTRYMKWSN\n\
  4 ACGTRYMKWSNacgt\n\
  5 ACGTRYMKWSBDHVN\n\
  6 ACGTRYMKWSBDHVacgtN\n\
  7 ACGTRYMKWSBDHVacgtrymkwsN\n\
  8 ACGTRYMKWSBDHVacgtrymkwsbdhvN\n\
models:\n\
";

int main(int argc, char **argv) {
	char  *file = NULL;
	char  *neg = NULL;
	int    len = 6;
	int    aid = 3;
	int    mod = 0;
	double mms = 1.0; // bits
	double mmh = 1.0; // bits
	double P1 = 0.97;
	double P2 = 0.49;
	double P3 = 0.33;
	double p1 = 0.70;
	double p2 = 0.40;
	double p3 = 0.30;

	// CLI - setup
	gkn_set_program_name(argv[0]);
	gkn_register_option("--len", 1);
	gkn_register_option("--aid", 1);
	gkn_register_option("--mod", 1);
	gkn_register_option("--mms", 1);
	gkn_register_option("--mmh", 1);
	gkn_register_option("--neg", 1);
	gkn_register_option("--P1",  1);
	gkn_register_option("--P2",  1);
	gkn_register_option("--P3" , 1);
	gkn_register_option("--p1",  1);
	gkn_register_option("--p2",  1);
	gkn_register_option("--p3",  1);
	gkn_parse_options(&argc, argv);
	if (argc == 1) gkn_exit("%s", usage);

	// CLI - harvest
	file = argv[1];
	if (gkn_option("--len")) len = atoi(gkn_option("--len"));
	if (gkn_option("--aid")) aid = atoi(gkn_option("--aid"));
	if (gkn_option("--mod")) mod = atoi(gkn_option("--mod"));
	if (gkn_option("--mms")) mms = atof(gkn_option("--mms"));
	if (gkn_option("--mmh")) mmh = atof(gkn_option("--mmh"));
	if (gkn_option("--neg")) neg = gkn_option("--neg");
	if (gkn_option("--P1"))  P1  = atof(gkn_option("--P1"));
	if (gkn_option("--P2"))  P2  = atof(gkn_option("--P2"));
	if (gkn_option("--P3"))  P3  = atof(gkn_option("--P3"));
	if (gkn_option("--p1"))  p1  = atof(gkn_option("--p1"));
	if (gkn_option("--p2"))  p2  = atof(gkn_option("--p2"));
	if (gkn_option("--p3"))  p3  = atof(gkn_option("--p3"));

	// read all sequences
	gkn_vec seqs = gkn_vec_new();
	gkn_pipe io = gkn_pipe_open(file, "r");
	gkn_fasta ff;
	while ((ff = gkn_fasta_read(io)) != NULL) {
		gkn_vec_push(seqs, ff->seq);
		gkn_vec_push(seqs, gkn_revcomp(ff->seq));
	}
	gkn_vec negs = gkn_vec_new();
	if (neg) {
		gkn_pipe nio = gkn_pipe_open(neg, "r");
		while ((ff = gkn_fasta_read(nio)) != NULL) {
			gkn_vec_push(negs, ff->seq);
			gkn_vec_push(negs, gkn_revcomp(ff->seq));
		}
	}

	// create background model
	gkn_pwm bkgd = background_model(seqs, len);
	//gkn_pwm_write(bkgd, stdout); // maybe an option to save it?

	// externally log2, internally loge
	mms *= log(2.718281828);
	mmh *= log(2.718281828);

	// initialize discretized alphabet and probabilities
	set_DNTP(P1, P2, P3, p1, p2, p3);
	char *alph = get_alphabet(aid);

	// scoring loop
	int limit = pow(strlen(alph), len);
	for (int n = 0; n < limit; n++) {
		gkn_pwm motif = num2pwm(alph, n, len);
		struct max_motif max = maximum_motif(motif, bkgd);
		if (max.score > mms * len && max.entropy > mmh) {
			double s = score_motif(seqs, motif, mod);
			printf("%s %s %f %f %f\n", motif->name,
				max.seq, max.score/log(2), max.entropy/log(2), s);

		}
		gkn_pwm_free(motif);
	}

	return 0;
}
