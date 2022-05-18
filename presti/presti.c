#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <float.h>

#include "genomikon.h"

#define CMAX 1000

static double ln_factorial(double n) {
	double f;
	if (n == 0) return 0;
	f = (0.5 * log(2 * 3.1415926))
		+ ((n + 0.5) * log((double)n))
		- n
		+ 1 / (12 * n)
		- 1 / (360 * pow((double)n, (double)3));
	return f;
}

static double memo_poisson[CMAX][CMAX];

static void init_memoizer (void) {
	for (int m = 0; m < CMAX; m++) {
		for (int n = 0; n < CMAX; n++) {
			memo_poisson[m][n] = -1; // value for unassigned
		}
	}
}

static double poisson(int m, int n) {
	if (memo_poisson[m][n] == -1) {
		memo_poisson[m][n] = (n * log((double)m)) - m - ln_factorial((int)n);
	}
	return memo_poisson[m][n];
}

static void yikes(const char *msg) {
	fprintf(stderr, "presti ERROR: %s\n", msg);
	exit(1);
}

static char *usage = "\
presti - convert streams of numbers into a few integer classes\n\n\
usage: presti <file-o-numbers> <n1> <n2> [<n...>] [options]\n\
options:\n\
  -x <float>  switch cost [5]\n\
  -y <float>  scaling factor [1.0]\n\
  -z <int>    maximum value [999]\n\
";

int main (int argc, char ** argv) {

	// Defaults
	double x = 5.0;
	double y = 1.0;
	int    z = 999;

	// CLI - setup
	gkn_set_program_name(argv[0]);
	gkn_register_option("-x", 1);
	gkn_register_option("-y", 1);
	gkn_register_option("-z", 1);
	gkn_register_option("-debug", 0);
	gkn_parse_options(&argc, argv);
	if (argc == 1) gkn_exit("%s", usage);

	// CLI - harvest
	if (gkn_option("-x")) x = atof(gkn_option("-x"));
	if (gkn_option("-y")) y = atof(gkn_option("-y"));
	if (gkn_option("-z")) z = atoi(gkn_option("-z"));
	if (z >= CMAX) yikes("z too high (999 max)");

	// cmdline - define HMM classes
	char *filename = argv[1];
	if (argc < 4) yikes(usage);
	gkn_ivec class = gkn_ivec_new(); // HMM states
	for (int i = 2; i < argc; i++) {
		int n = atoi(argv[i]);
		if (n < 1) yikes("minimum class value = 1");
		if (n >= CMAX) yikes("exceeded maximim class value");
		for (int j = 0; j < class->size; j++) {
			if (n == class->elem[j]) yikes("classes must be unique");
		}
		gkn_ivec_push(class, n);
	}

	// read file - define sequence of numbers
	gkn_ivec value = gkn_ivec_new(); // sequence of ints
	FILE *fp;
	char line[32];
	if ((fp = fopen(filename, "r")) == NULL) yikes("file open error");
	while (fgets(line, sizeof(line), fp) != NULL) {
		if (line[0] == '#') continue;
		float f;
		if (sscanf(line, "%f", &f) != 1) yikes("file read error");
		f *= y;
		if (f > z) f = z;
		gkn_ivec_push(value, (int)f);
	}

	// simplifications
	int states = class->size;
	int *state = class->elem;
	int len    = value->size;
	int *seq   = value->elem;

	// create data structures
	double **score = malloc(states * sizeof(double*));
	int    **trace = malloc(states * sizeof(int*));
	for (int i = 0; i < states; i++) {
		score[i] = calloc(sizeof(double), len);
		trace[i] = calloc(sizeof(int)   , len);
	}

	// init
	init_memoizer();
	for (int i = 0; i < states; i++) {
		score[i][0] = poisson(state[i], seq[0]);
		trace[i][0] = i;
	}

	// fill
	for (int i = 1; i < len; i++) {
		for (int j = 0; j < states; j++) { // current state
			int max_state = -1;
			double max_score = -DBL_MAX;
			double emit = poisson(state[j], seq[i]);
			for (int k = 0; k < states; k++) { // previous state
				double s = emit + score[k][i-1];
				if (j != k) s -= x; // switch cost
				if (s > max_score) {
					max_score = s;
					max_state = k;
				}
			}
			score[j][i] = max_score;
			trace[j][i] = max_state;
		}
	}

	// traceback
	double max_score = -DBL_MAX;
	int max_state = -1;
	for (int i = 0; i < states; i++) {
		if (score[i][len-1] > max_score) {
			max_score = score[i][len-1];
			max_state = i;
		}
	}

	int beg = len;
	int end = len;
	gkn_tvec out = gkn_tvec_new();
	char text[64];
	for (int i = len-1; i > 0; i--) {
		if (trace[max_state][i] == max_state) {
			beg = i;
		} else {
			sprintf(text, "%d\t%d\t%d\n", beg, end, state[max_state]);
			gkn_tvec_push(out, text);
			beg = i;
			end = i;
			max_state = trace[max_state][i];
		}
	}
	sprintf(text, "%d\t%d\t%d\n", beg, end, state[max_state]);
	gkn_tvec_push(out, text);

	// output
	printf("# beg\tend\tclass\n");
	for (int i = out->size -1; i >= 0; i--) {
		printf("%s", out->elem[i]);
	}

	// debugging - not advertised, what's the end of the trellis look like?
	if (gkn_option("-debug")) {
		for (int i = 0; i < states; i++) {
			fprintf(stderr, "%d %f\n", i, score[i][len-1]);
		}
	}

	return 0;
}
