
#include "genomikon.h"
#include "dmg.h"

static char *usage = "\
usage: motifamatic <fasta> [options]\n\
options:\n\
  --len <int>   motif length [6]\n\
  --sig <int>   motif signature: 4, 5, 9, 11, 15, 19, 25, 29 [4]\n\
  --mod <...>   not sure yet...\n\
  --P1  <float> uppercase single letter probability [0.997]\n\
  --P2  <float> uppercase double letter probability [0.498]\n\
  --P3  <float> uppercase triple letter probability [0.333]\n\
  --p1  <float> lowercase single letter probability [0.800]\n\
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
	int   mod = 0;
	double P1 = 0.997;
	double P2 = 0.498;
	double P3 = 0.333;
	double p1 = 0.800;
	double p2 = 0.400;
	double p3 = 0.300;

	// CLI - setup
	gkn_set_program_name(argv[0]);
	gkn_register_option("--len", 1);
	gkn_register_option("--sig", 1);
	gkn_register_option("--mod", 1);
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
	if (gkn_option("--sig")) sig = atoi(gkn_option("--sig"));
	if (gkn_option("--mod")) mod = atoi(gkn_option("--mod"));
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
	while ((ff = gkn_fasta_read(io)) != NULL) gkn_vec_push(seqs, ff);
	
	// create background model
	gkn_pwm bkgd = background_model(seqs, len);
	gkn_pwm_write(bkgd, stdout); // maybe an option to save it?
	
	// scoring loop
	dmgen dmg = dmgen_new_custom(sig, P1, P2, P3, p1, p2, p3);
	int limit = pow(sig, len);
	for (int n = 0; n < limit; n++) {
		gkn_pwm motif = num2pwm(dmg, n, len);
		printf("%s\n", motif->name);
		for (int i = 0; i < seqs->size; i++) {
			ff = seqs->elem[i];
			double s = score_motif(ff->seq, motif, mod);
			printf("%d %f\n", i, s);
		}
		gkn_pwm_free(motif);
	}
	

	return 0;
}
