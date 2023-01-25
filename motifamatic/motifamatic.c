#include "genomikon.h"

static char *usage = "\
motifamatic - find motifs by enumerating digitized representations\n\n\
usage: motifamatic <fasta file> [options]\n\
options:\n\
  --len <int>   motif length [6]\n\
  --mes <int>   motif enumeration scheme: 4, 5, 8, 9, 15, 19, 25 [4]\n\
  --P1  <float> uppercase single letter probability [0.97]\n\
  --P2  <float> uppercase double letter probability [0.48]\n\
  --p1  <float> lowercase single letter probability [0.85]\n\
  --p2  <float> lowercase double letter probability [0.40]\n\
";

// p-value cutoffs for finding motifs
// background model types...
// probability model like oops, zoops, anr
// complexity filters for shit motifs

char *S4  = "ACGT";
char *S5  = "ACGTN";
char *S8  = "ACGTacgt";
char *S9  = "ACGTacgtN";
char *S15 = "ACGTRYMKWSBDHVN";
char *S19 = "ACGTRYMKWSBDHVNacgt";
char *S25 = "ACGTRYMKWSBDHVNacgtrymkws";

double P1 = 0.970;
double P2 = 0.480;
double P3 = 0.333;
double p1 = 0.850;
double p2 = 0.400;
double p3 = 0.300;
double N4 = 0.250;
double Q1;
double Q2;
double Q3;
double q1;
double q2;
double q3;

static char *num2mes(int num) {
	switch (num) {
		case 4:  return S4;  break;
		case 5:  return S5;  break;
		case 8:  return S8;  break;
		case 9:  return S9;  break;
		case 15: return S15; break;
		case 19: return S19; break;
		case 25: return S25; break;
		default: gkn_exit("mes %d not supported", num);
	}
	return NULL;
}

static char * num2str(int num, int len, int mes) {
	char *str = malloc(len + 1);
	char *src = num2mes(mes);

	str[len] = '\0';
	for (int i = 0; i < len; i++) {
		int max = pow(mes, (len-i-1));
		int r = 0;
		if (num > max -1) {
			r = num / max;
			num -= r * max;
		}
		str[i] = src[r];
	}
	return str;
}

static double ** num2pwm(int num, int len, int mes) {
	double ** pwm = malloc(sizeof(double*) * len);
	for (int i = 0; i < len; i++) {
		pwm[i] = malloc(sizeof(double) * 4);
	}

	char *str = num2str(num, len, mes);
	for (int i = 0; i < len; i++) {
		switch (str[i]) {
		case 'A': pwm[i][0]=P1; pwm[i][1]=Q1; pwm[i][2]=Q1; pwm[i][3]=Q1; break;
		case 'C': pwm[i][0]=Q1; pwm[i][1]=P1; pwm[i][2]=Q1; pwm[i][3]=Q1; break;
		case 'G': pwm[i][0]=Q1; pwm[i][1]=Q1; pwm[i][2]=P1; pwm[i][3]=Q1; break;
		case 'T': pwm[i][0]=Q1; pwm[i][1]=Q1; pwm[i][2]=Q1; pwm[i][3]=P1; break;
		case 'R': pwm[i][0]=P2; pwm[i][1]=Q2; pwm[i][2]=P2; pwm[i][3]=Q2; break;
		case 'Y': pwm[i][0]=Q2; pwm[i][1]=P2; pwm[i][2]=Q2; pwm[i][3]=P2; break;
		case 'M': pwm[i][0]=P2; pwm[i][1]=P2; pwm[i][2]=Q2; pwm[i][3]=Q2; break;
		case 'K': pwm[i][0]=Q2; pwm[i][1]=Q2; pwm[i][2]=P2; pwm[i][3]=P2; break;
		case 'W': pwm[i][0]=P2; pwm[i][1]=Q2; pwm[i][2]=Q2; pwm[i][3]=P2; break;
		case 'S': pwm[i][0]=Q2; pwm[i][1]=P2; pwm[i][2]=P2; pwm[i][3]=Q2; break;
		case 'B': pwm[i][0]=Q3; pwm[i][1]=P3; pwm[i][2]=P3; pwm[i][3]=P3; break;
		case 'D': pwm[i][0]=P3; pwm[i][1]=Q3; pwm[i][2]=P3; pwm[i][3]=P3; break;
		case 'H': pwm[i][0]=P3; pwm[i][1]=P3; pwm[i][2]=Q3; pwm[i][3]=P3; break;
		case 'V': pwm[i][0]=P3; pwm[i][1]=P3; pwm[i][2]=P3; pwm[i][3]=Q3; break;
		case 'a': pwm[i][0]=p1; pwm[i][1]=q1; pwm[i][2]=q1; pwm[i][3]=q1; break;
		case 'c': pwm[i][0]=q1; pwm[i][1]=p1; pwm[i][2]=q1; pwm[i][3]=q1; break;
		case 'g': pwm[i][0]=q1; pwm[i][1]=q1; pwm[i][2]=p1; pwm[i][3]=q1; break;
		case 't': pwm[i][0]=q1; pwm[i][1]=q1; pwm[i][2]=q1; pwm[i][3]=p1; break;
		case 'r': pwm[i][0]=p2; pwm[i][1]=q2; pwm[i][2]=p2; pwm[i][3]=q2; break;
		case 'y': pwm[i][0]=q2; pwm[i][1]=p2; pwm[i][2]=q2; pwm[i][3]=p2; break;
		case 'm': pwm[i][0]=p2; pwm[i][1]=p2; pwm[i][2]=q2; pwm[i][3]=q2; break;
		case 'k': pwm[i][0]=q2; pwm[i][1]=q2; pwm[i][2]=p2; pwm[i][3]=p2; break;
		case 'w': pwm[i][0]=p2; pwm[i][1]=q2; pwm[i][2]=q2; pwm[i][3]=p2; break;
		case 's': pwm[i][0]=q2; pwm[i][1]=p2; pwm[i][2]=p2; pwm[i][3]=q2; break;
		case 'b': pwm[i][0]=q3; pwm[i][1]=p3; pwm[i][2]=p3; pwm[i][3]=p3; break;
		case 'd': pwm[i][0]=p3; pwm[i][1]=q3; pwm[i][2]=p3; pwm[i][3]=p3; break;
		case 'h': pwm[i][0]=p3; pwm[i][1]=p3; pwm[i][2]=q3; pwm[i][3]=p3; break;
		case 'v': pwm[i][0]=p3; pwm[i][1]=p3; pwm[i][2]=p3; pwm[i][3]=q3; break;
		case 'N': pwm[i][0]=N4; pwm[i][1]=N4; pwm[i][2]=N4; pwm[i][3]=N4; break;
		default: gkn_exit("impossible");
		}
	}
	return pwm;
}

int main(int argc, char **argv) {
	char *file = NULL; // path to fasta file
	int   len = 6;   // motif length
	int   mes = 4;   // alphabet type

	gkn_pipe  io = NULL; // for reading fasta files
	gkn_fasta ff = NULL; // for individual fasta entries

	// CLI - setup
	gkn_set_program_name(argv[0]);
	gkn_register_option("--len", 1);
	gkn_register_option("--mes", 1);
	gkn_register_option("--P1", 1);
	gkn_register_option("--P2", 1);
	gkn_register_option("--p1", 1);
	gkn_register_option("--p2", 1);
	gkn_parse_options(&argc, argv);
	if (argc == 1) gkn_exit("%s", usage);

	// CLI - harvest
	file = argv[1];
	if (gkn_option("--len")) len = atoi(gkn_option("--len"));
	if (gkn_option("--mes")) mes = atoi(gkn_option("--mes"));
	if (gkn_option("--P1"))  P1  = atof(gkn_option("--P1"));
	if (gkn_option("--P2"))  P1  = atof(gkn_option("--P1"));
	if (gkn_option("--p1"))  p1  = atof(gkn_option("--p1"));
	if (gkn_option("--p2"))  p2  = atof(gkn_option("--p2"));

	// init Q values
	Q1 = (1 - P1) / 3;
	Q2 = (1 - 2*P2) / 2;
	Q3 = (1 - 3*P3) / 3;
	q1 = (1 - p1) / 3;
	q2 = (1 - 2*p2) / 2;
	q3 = (1 - 3*p3) / 3;
	fprintf(stderr, "%.0f possible motifs of length %d in alphabet %s\n",
		pow(mes, len), len, num2mes(mes));

	// main loop
	io = gkn_pipe_open(file, "r");
	while ((ff = gkn_fasta_read(io)) != NULL) {
		// read in all seqs
	}

	// find motifs in all seqs
	// report aggregate p-values


	printf("%s\n", num2str(500, len, mes));
	double **pwm = num2pwm(500, len, mes);
	for (int i = 0; i < len; i++) {
		printf("%d", i);
		for (int j = 0; j < 4; j++) {
			printf(" %.3f", pwm[i][j]);
		}
		printf("\n");
	}

	return 0;
}
