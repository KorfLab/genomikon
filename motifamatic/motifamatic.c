#include "genomikon.h"

static char *usage = "\
motifamatic - find motifs by enumerating digitized representations\n\n\
usage: motifamatic <fasta file> [options]\n\
options:\n\
  --len <int>   motif length [6]\n\
  --mod <int>   motif model: 4, 5, 6, 7, 11, 15, 21 [4]\n\
  --P1  <float> uppercase single letter probability [0.97]\n\
  --P2  <float> uppercase double letter probability [0.48]\n\
  --p1  <float> lowercase single letter probability [0.85]\n\
  --p2  <float> lowercase double letter probability [0.40]\n\
";

char *S4  = "ACGT";
char *S5  = "ACGTN";
char *S6  = "ACGTRY";
char *S7  = "ACGTRYN";
char *S11 = "ACGTRYMKWSN";
char *S15 = "ACGTRYMKWSNacgt";
char *S21 = "ACGTRYMKWSNacgtrymkws";

double P1 = 0.97;
double P2 = 0.48;
double p1 = 0.85;
double p2 = 0.40;

static char *num2mod(int num) {
	switch (num) {
		case 4:  return S4;  break;
		case 5:  return S5;  break;
		case 6:  return S6;  break;
		case 7:  return S7;  break;
		case 11: return S11; break;
		case 15: return S15; break;
		case 21: return S21; break;
		default: gkn_exit("impossible");
	}
	return NULL;
}

static char * num2str(int num, int len, int mod) {
	char *str = malloc(len + 1);
	char *src = num2mod(mod);

	str[len] = '\0';
	for (int i = 0; i < len; i++) {
		int max = pow(mod, (len-i-1));
		int r = 0;
		if (num > max -1) {
			r = num / max;
			num -= r * max;
		}
		str[i] = src[r];
	}
	return str;
}

/*
static double ** num2pwm(int num, int len, int mod) {

}
*/

int main(int argc, char **argv) {
	char *file = NULL; // path to fasta file
	int   len = 6;   // motif length
	int   mod = 4;   // alphabet type

	gkn_pipe  io = NULL; // for reading fasta files
	gkn_fasta ff = NULL; // for individual fasta entries

	// CLI - setup
	gkn_set_program_name(argv[0]);
	gkn_register_option("--len", 1);
	gkn_register_option("--mod", 1);
	gkn_register_option("--P1", 1);
	gkn_register_option("--P2", 1);
	gkn_register_option("--p1", 1);
	gkn_register_option("--p2", 1);
	gkn_parse_options(&argc, argv);
	if (argc == 1) gkn_exit("%s", usage);

	// CLI - harvest
	file = argv[1];
	if (gkn_option("--len")) len = atoi(gkn_option("--len"));
	if (gkn_option("--mod")) mod = atoi(gkn_option("--mod"));
	if (gkn_option("--P1"))  P1  = atof(gkn_option("--P1"));
	if (gkn_option("--P2"))  P1  = atof(gkn_option("--P1"));
	if (gkn_option("--p1"))  p1  = atof(gkn_option("--p1"));
	if (gkn_option("--p2"))  p2  = atof(gkn_option("--p2"));

	fprintf(stderr, "%.0f possible motifs of length %d in alphabet %s\n",
		pow(mod, len), len, num2mod(mod));

	// find isoforms (just one fasta file entry)
	io = gkn_pipe_open(file, "r");
	while ((ff = gkn_fasta_read(io)) != NULL) {
		printf("%s\n", ff->def);
	}

	printf("%d %d\n", len, mod);
	printf("%d %s\n", 5, num2str(55, len, mod));

	return 0;
}
