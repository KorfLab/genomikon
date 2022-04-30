#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <float.h>

#include "genomikon.h"

double LOW_SCORE = -999;

static void squeek(const char *msg) {
	fprintf(stderr, "hmmstar ERROR: %s\n", msg);
	exit(1);
}

static gkn_map read_kmer_file(const char *filename) {
	gkn_map map = gkn_map_new();
	FILE *fp;
	char line[64];
	if ((fp = fopen(filename, "r")) == NULL) squeek("file open error");
	while (fgets(line, sizeof(line), fp) != NULL) {
		if (line[0] == '#') continue;
		float f;
		char kmer[16];
		if (sscanf(line, "%s %f", kmer, &f) != 2) squeek("file read error");
		double *dp = malloc(sizeof(double));
		*dp = f;
		gkn_map_set(map, kmer, dp);
	}
	return map;
}

struct Path {
	int    beg;    // 1-based for output
	int    end;    // 1-based for output
	int    sid;    // state index
};
typedef struct Path *path;

static path new_path(int beg, int end, int sid, int k) {
	path p = malloc(sizeof(struct Path));
	p->beg = beg + k / 2;
	p->end = end + k / 2;
	p->sid = sid;
	return p;
}

static gkn_vec viterbi(const char *seq, int K, gkn_vec mods, double **tm) {

	// allocate memory
	double **score = malloc(sizeof(double*) * mods->size);
	int    **trace = malloc(sizeof(int*)    * mods->size);
	for (int j = 0; j < mods->size; j++) {
		score[j] = malloc(sizeof(double) * (strlen(seq) +1));
		trace[j] = malloc(sizeof(int)    * (strlen(seq) +1));
	}

	// init
	for (int j = 0; j < mods->size; j++) {
		score[j][0] = 0;
		trace[j][0] = -1;
	}

	// fill
	char kmer[16];
	for (int i = 1; i < strlen(seq) -K +2; i++) {
		strncpy(kmer, seq+(i-1), K);
		kmer[K] = '\0';

		for (int j = 0; j < mods->size; j++) {
			int max_state = -1;
			double max_score = -DBL_MAX;
			gkn_map table = mods->elem[j];
			double emit = *(double*)gkn_map_get(table, kmer);

			for (int k = 0; k < mods->size; k++) {
				double prev  = score[k][i-1];
				double trans = tm[k][j];
				double score = prev + emit + trans;
				if (score > max_score) {
					max_score = score;
					max_state = k;
				}
			}
			score[j][i] = max_score;
			trace[j][i] = max_state;
		}
	}

	// trace
	int last = strlen(seq);
	double max_score = score[0][last];
	int    max_state = trace[0][last];
	for (int j = 1; j < mods->size; j++) {
		if (score[j][last] > max_score) {
			max_score = score[j][last];
			max_state = trace[j][last];
		}
	}

	int beg = last -1;
	int end = last -1;
	gkn_vec paths = gkn_vec_new();
	for (int i = last -1; i > 0; i--) {
		if (trace[max_state][i] == max_state) {
			beg = i;
		} else {
			gkn_vec_push(paths, (void*)new_path(beg, end, max_state, K));
			beg = i;
			end = i;
			max_state = trace[max_state][i];
		}
	}
	gkn_vec_push(paths, (void*)new_path(beg, end, max_state, K));

	// free memory
	for (int  j = 0; j < mods->size; j++) {
		free(score[j]);
		free(trace[j]);
	}
	free(score);
	free(trace);

	// reverse paths
	for (int i = 0; i < paths->size / 2; i++) {
		void *temp = paths->elem[i];
		paths->elem[i] = paths->elem[paths->size -i -1];
		paths->elem[paths->size -i -1] = temp;
	}

	return paths;
}

static char *usage = "\
hmmstar - nucleotide viterbi decoder with k-mer tables\n\n\
usage: hmmstar <fast file> <table1-> <table-2> [<table-n...>] [options]\n\
options:\n\
  -x <float>  switch probability X [0.001]\n\
  -y <float>  switch probability Y [0.01]\n\
  -f          use full model (default star)\n\
";

int main (int argc, char ** argv) {

	// Defaults
	double X = 0.001;
	double Y = 0.01;
	int STAR = 1;

	// CLI - setup
	gkn_set_program_name(argv[0]);
	gkn_register_option("-x", 1);
	gkn_register_option("-y", 1);
	gkn_register_option("-f", 0);
	gkn_parse_options(&argc, argv);
	if (argc <= 3) gkn_exit("%s", usage);

	// CLI - harvest
	char *fasta = argv[1];
	gkn_tvec files = gkn_tvec_new();
	for (int i = 2; i < argc; i++) {
		gkn_tvec_push(files, argv[i]);
	}
	if (gkn_option("-x")) X = atof(gkn_option("-x"));
	if (gkn_option("-y")) Y = atof(gkn_option("-y"));
	if (gkn_option("-f")) STAR = 0;

	// read K-mer files
	gkn_vec mods = gkn_vec_new();
	for (int i = 0; i < files->size; i++) {
		gkn_vec_push(mods, (void*)read_kmer_file(files->elem[i]));
	}
	int K = strlen((char*)((gkn_map)mods->elem[0])->keys->elem[0]); // sheesh

	// create null model - could be marginal probs...
	gkn_tvec kmers = gkn_map_keys(mods->elem[0]);
	gkn_map null_model = gkn_map_new();
	for (int i = 0; i < kmers->size; i++) {
		char *kmer = kmers->elem[i];
		double *null_value = malloc(sizeof(double *));
		*null_value = pow(4, -K);
		gkn_map_set(null_model, kmer, null_value);
	}

	// star shape adds an additional state
	if (STAR) {
		gkn_vec_push(mods, (void*) null_model);
		gkn_tvec_push(files, "Null");
	}

	// convert models to log-odds with null expected
	for (int i = 0; i < mods->size; i++) {
		gkn_map table = mods->elem[i];
		for (int j = 0; j < table->keys->size; j++) {
			char *kmer = table->keys->elem[j];
			void *v1 = gkn_map_get(table, kmer);
			void *v2 = gkn_map_get(null_model, kmer);
			double lod = log(*(double*)v1 / *(double*)v2);
			*(double*)v1 = lod;
		}
	}

	// transition matrix
	double **trans = malloc(mods->size * sizeof(double*));
	for (int i = 0; i < mods->size; i++) {
		trans[i] = malloc(mods->size * sizeof(double*));
	}

	if (STAR) {
		int null = mods->size -1; // index of null model
		for (int i = 0; i < mods->size; i++) {
			for (int j = 0; j < mods->size; j++) {
				if (i == null && j == null) trans[i][j] = 1-(mods->size-1) * X;
				else if (j == null)         trans[i][j] = Y;
				else if (i == null)         trans[i][j] = X;
				else if (i == j)            trans[i][j] = 1 - Y;
				else                        trans[i][j] = 0;
			}
		}
	} else {
		for (int i = 0; i < mods->size; i++) {
			for (int j = 0; j < mods->size; j++) {
				if (i == j) trans[i][j] = 1 - (mods->size -1) * X;
				else        trans[i][j] = X;
			}
		}
	}

	// convert all transition probabilities to log space
	for (int i = 0; i < mods->size; i++) {
		for (int j = 0; j < mods->size; j++) {
			if (trans[i][j] == 0) trans[i][j] = LOW_SCORE;
			else                  trans[i][j] = log(trans[i][j]);
		}
	}

	// main loop
	gkn_pipe io = gkn_pipe_open(fasta, "r");
	gkn_fasta in;
	while ((in = gkn_fasta_read(io)) != NULL) {
		gkn_vec paths = viterbi(in->seq, K, mods, trans);
		printf("%s\n", in->def);
		for (int i = 0; i < paths->size; i++) {
			path p = paths->elem[i];
			printf("%d\t%d\t%s\n", p->beg, p->end, files->elem[p->sid]);
		}
	}

	return 0;
}
