#include "genomikon.h"
#include "isoform.h"

static char * random_seq(int length) {
	char *seq = malloc(length + 1);
	for (int i = 0; i < length; i++) {
		switch (rand() % 4) {
			case 0: seq[i] = 'A'; break;
			case 1: seq[i] = 'C'; break;
			case 2: seq[i] = 'G'; break;
			case 3: seq[i] = 'T'; break;
		}
	}
	seq[length] = '\0';
	return seq;
}

static char *usage = "\
isorandom - count all possible isoforms from random sequences\n\n\
usage: isorandom <length> <count> [options]\n\
options:\n\
  --min_exon   <int>   minimum exon length [25]\n\
  --min_intron <int>   minimum intron length [35]\n\
  --max_splice <int>   maximum splices [3]\n\
  --flank      <int>   genomic flank lengths [99]\n\
  --seed       <int>   set random seed\n\
";

int main(int argc, char **argv) {
	int   count;       // number of sequences to generate
	int   length;      // length of sequences to generate
	int   emin = 25;   // min exon size
	int   imin = 35;   // min intron size
	int   smax = 3;    // max splices
	int   gen  = 99;   // genomic flank (promoter, downstream)

	// CLI - setup
	gkn_set_program_name(argv[0]);
	gkn_register_option("--min_exon", 1);
	gkn_register_option("--min_intron", 1);
	gkn_register_option("--max_splice", 1);
	gkn_register_option("--flank",  1);
	gkn_parse_options(&argc, argv);
	if (argc == 1) gkn_exit("%s", usage);

	// CLI - harvest
	length = atoi(argv[1]);
	count  = atoi(argv[2]);
	if (gkn_option("--min_exon"))   emin  = atoi(gkn_option("--min_exon"));
	if (gkn_option("--min_intron")) imin  = atoi(gkn_option("--min_intron"));
	if (gkn_option("--max_splice")) smax  = atoi(gkn_option("--max_splice"));
	if (gkn_option("--flank"))      gen   = atoi(gkn_option("--flank"));
	if (gkn_option("--seed"))       srand(atoi(gkn_option("--seed")));

	// main loop
	for (int i = 0; i < count; i++) {
		char *seq = random_seq(length);
		isozone iso = isoforms(seq, emin, imin, smax, gen, NULL, 1);
		printf("%d\t%d\t%d\t%d\t%d\n", (int)strlen(seq), iso->dons, iso->accs,
			iso->trials, iso->isoforms);
		isozone_free(iso);
		free(seq);
	}

	return 0;
}
