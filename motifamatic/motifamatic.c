
#include "genomikon.h"
#include "dmg.h"

static char *usage = "\
motifamatic - find motifs by enumerating digitized representations\n\n\
usage: motifamatic <fasta file> [options]\n\
options:\n\
  --len <int>   motif length [6]\n\
  --sig <int>   motif signature number: 4, 5, 9, 11, 15, 19, 25, 31 [4]\n\
  --P1  <float> uppercase single letter probability [0.970]\n\
  --P2  <float> uppercase double letter probability [0.480]\n\
  --32  <float> uppercase triple letter probability [0.333]\n\
  --p1  <float> lowercase single letter probability [0.900]\n\
  --p2  <float> lowercase double letter probability [0.400]\n\
  --p2  <float> lowercase double letter probability [0.300]\n\
";

// p-value cutoffs for finding motifs
// background model types...
// probability model like oops, zoops, anr
// complexity filters for shit motifs

int main(int argc, char **argv) {
	char *file = NULL;
	int   len = 6;
	int   sig = 4;
	double P1 = 0.970;
	double P2 = 0.480;
	double P3 = 0.333;
	double p1 = 0.900;
	double p2 = 0.400;
	double p3 = 0.300;

	gkn_pipe  io = NULL; // for reading fasta files
	gkn_fasta ff = NULL; // for individual fasta entries

	// CLI - setup
	gkn_set_program_name(argv[0]);
	gkn_register_option("--len", 1);
	gkn_register_option("--sig", 1);
	gkn_register_option("--P1", 1);
	gkn_register_option("--P2", 1);
	gkn_register_option("--P3", 1);
	gkn_register_option("--p1", 1);
	gkn_register_option("--p2", 1);
	gkn_register_option("--p3", 1);
	gkn_parse_options(&argc, argv);
	if (argc == 1) gkn_exit("%s", usage);

	// CLI - harvest
	file = argv[1];
	if (gkn_option("--len")) len = atoi(gkn_option("--len"));
	if (gkn_option("--sig")) sig = atoi(gkn_option("--sig"));
	if (gkn_option("--P1"))  P1  = atof(gkn_option("--P1"));
	if (gkn_option("--P2"))  P2  = atof(gkn_option("--P2"));
	if (gkn_option("--P3"))  P3  = atof(gkn_option("--P3"));
	if (gkn_option("--p1"))  p1  = atof(gkn_option("--p1"));
	if (gkn_option("--p2"))  p2  = atof(gkn_option("--p2"));
	if (gkn_option("--p3"))  p3  = atof(gkn_option("--p3"));

	dmgen dmg = dmgen_new_custom(sig, P1, P2, P3, p1, p2, p3);
	printf("%d %d %s %s\n", len, dmg->sig, dmg->alph, num2str(dmg, 25, 4));



	// main loop
	io = gkn_pipe_open(file, "r");
	while ((ff = gkn_fasta_read(io)) != NULL) {
		// read in all seqs
	}

	// find motifs in all seqs
	// report aggregate p-values

/*
	printf("%s\n", num2str(500, len, mes));
	double **pwm = num2pwm(500, len, mes);
	for (int i = 0; i < len; i++) {
		printf("%d", i);
		for (int j = 0; j < 4; j++) {
			printf(" %.3f", pwm[i][j]);
		}
		printf("\n");
	}
*/

	return 0;
}
