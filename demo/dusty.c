/*****************************************************************************\
 dusty.c
\*****************************************************************************/

#include "genomikon.h"

static double entropy(int a, int c, int g, int t) {
	int total = a + c + g + t;
	if (total == 0) return 0;

	double pa = a / (double)total;
	double pc = c / (double)total;
	double pg = g / (double)total;
	double pt = t / (double)total;

	double h = 0;
	if (a > 0) h -= pa * log(pa);
	if (c > 0) h -= pc * log(pc);
	if (g > 0) h -= pg * log(pg);
	if (t > 0) h -= pt * log(pt);

	return h/log(2);
}

static char *dust1(const gkn_fasta ff, int w, int t, int lc) {
	int A, C, G, T;

	char *mask = malloc(ff->length+1);
	strcpy(mask, ff->seq);
	for (int i = 0; i < ff->length; i++) mask[i] = toupper(mask[i]);

	for (int i = 0; i < ff->length -w + 1; i++) {
		A = 0; C = 0; G = 0; T = 0;
		for (int j = 0; j < w; j++) {
			switch (mask[i+j]) {
				case 'A': A++; break;
				case 'C': C++; break;
				case 'G': G++; break;
				case 'T': T++; break;
			}
		}
		double h = entropy(A, C, G, T);
		if (h < t) {
			int pos = i + w/2;
			if (lc) mask[pos] = tolower(mask[pos]);
			else    mask[pos] = 'N';
		}
	}
	return mask;
}

static char *dust2(const gkn_fasta ff, int w, int t, int lc) {
	int A, C, G, T;
	char *mask = malloc(ff->length+1);
	strcpy(mask, ff->seq);
	for (int i = 0; i < ff->length; i++) mask[i] = toupper(mask[i]);


	// first window
	A = 0; C = 0; G = 0; T = 0;
	for (int i = 0; i < w; i++) {
		switch (mask[i]) {
			case 'A': A++; break;
			case 'C': C++; break;
			case 'G': G++; break;
			case 'T': T++; break;
		}
	}

	double h = entropy(A, C, G, T);
	if (h < t) {
		int pos = w/2;
		if (lc) mask[pos] = tolower(mask[pos]);
		else    mask[pos] = 'N';
	}

	// remaining windows
	for (int i = 1; i < ff->length -w + 1; i++) {
		char on  = mask[i+w-1];
		char off = mask[i-1];

		switch (on) {
			case 'A': A++; break;
			case 'C': C++; break;
			case 'G': G++; break;
			case 'T': T++; break;
		}

		switch (off) {
			case 'A': A--; break;
			case 'C': C--; break;
			case 'G': G--; break;
			case 'T': T--; break;
		}

		double h = entropy(A, C, G, T);
		if (h < t) {
			int pos = i + w/2;
			if (lc) mask[pos] = tolower(mask[pos]);
			else    mask[pos] = 'N';
		}
	}

	return mask;
}

static char *usage = "\
dusty - low complexity demo\n\n\
usage: dusty <file> [options]\n\
options:\n\
  -w <int>    window size [11]\n\
  -h <float>  entropy threshold [1.1]\n\
  -n          use Ns for masking (default lowercase)\n\
  -a <int>    algorithm [1]";

int main(int argc, char **argv) {
	char *file = NULL; // path to fasta file
	int   w = 11;      // window size
	float h = 1.1;     // entropy
	int   lc = 1;      // lowercase?
	int   alg = 1;     // algorithm
	gkn_fasta in, out;

	// Command Line Interface
	gkn_set_program_name(argv[0]);
	gkn_register_option("-i", 1);
	gkn_register_option("-w", 1);
	gkn_register_option("-h", 1);
	gkn_register_option("-n", 0);
	gkn_register_option("-a", 1);
	gkn_parse_options(&argc, argv);

	if (argc == 1) gkn_exit("%s", usage);

	file = argv[1];
	if (gkn_option("-w")) w = atoi(gkn_option("-w"));
	if (gkn_option("-h")) h = atof(gkn_option("-h"));
	if (gkn_option("-n")) lc = 0;
	if (gkn_option("-a")) alg = atoi(gkn_option("-a"));

	// main loop
	gkn_pipe io = gkn_pipe_open(file, "r");
	while ((in = gkn_fasta_read(io)) != NULL) {
		char *mask = NULL;
		switch (alg) {
			case 1: mask = dust1(in, w, h, lc); break;
			case 2: mask = dust2(in, w, h, lc); break;
			default: gkn_exit("algorithm out of range");
		}

		out = gkn_fasta_new(in->def, mask);
		gkn_fasta_write(stdout, out);

		free(mask);
		gkn_fasta_free(in);
		gkn_fasta_free(out);
	}

	return 0;
}
