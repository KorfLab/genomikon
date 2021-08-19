/******************************************************************************\
 model.c
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef GENOMIKON_MODEL_C
#define GENOMIKON_MODEL_C

#include "model.h"

// utilities

int dna2dec(const char *seq, int off, int k) {
	int idx = 0;
	for (int i = 0; i < k; i++) {
		switch (seq[off+i]) {
			case 'A': case 'a': idx += pow(4, (k -i -1)) * 0; break;
			case 'C': case 'c': idx += pow(4, (k -i -1)) * 1; break;
			case 'G': case 'g': idx += pow(4, (k -i -1)) * 2; break;
			case 'T': case 't': idx += pow(4, (k -i -1)) * 3; break;
			default: return -1;
		}
	}
	return idx;
}

double prob2score(double p) {
	if (p == 0) return -100; // umm...
	return log(p/0.25) / log(2);
}

// PWM

void gkn_pwm_free(gkn_pwm pwm) {
	free(pwm->name);
	for (int i = 0; i < pwm->size; i++) {
		free(pwm->score[i]);
	}
	free(pwm->score);
	free(pwm);
}

gkn_pwm gkn_pwm_read(const char *filename) {
	char    *line = NULL;
	size_t   len = 0;
	ssize_t  read;
	gkn_pipe io  = gkn_pipe_open(filename, "r");
	char     blah[256];
	int      size;
	double **score = NULL;
	double   a, c, g, t;
	int      row = 0;

	while ((read = getline(&line, &len, io->stream)) != -1) {
		if (line[0] == '#') {
			assert(sscanf(line, "# PWM %s %d", blah, &size) == 2);
			score = malloc(sizeof(double*) * size);
			for (int i = 0; i < size; i++) {
				score[i] = malloc(sizeof(double) * 4);
			}
		} else if (sscanf(line, "%lf %lf %lf %lf", &a, &c, &g, &t) == 4) {
			score[row][0] = prob2score(a);
			score[row][1] = prob2score(c);
			score[row][2] = prob2score(g);
			score[row][3] = prob2score(t);
			row++;
		}
	}
	gkn_pipe_close(io);
	if (line) free(line);

	gkn_pwm model = malloc(sizeof(struct gkn_PWM));
	model->name = malloc(strlen(filename)+1);
	strcpy(model->name, filename);
	model->size = size;
	model->score = score;

	return model;
}

double gkn_pwm_score(const gkn_pwm pwm, const char *seq, int pos) {
	double p = 0;
	for (int i = 0; i < pwm->size; i++) {
		switch (seq[i+pos]) {
			case 'A': case 'a': p += pwm->score[i][0]; break;
			case 'C': case 'c': p += pwm->score[i][1]; break;
			case 'G': case 'g': p += pwm->score[i][2]; break;
			case 'T': case 't': p += pwm->score[i][3]; break;
		}
	}
	return p;
}

// Markov model

void gkn_mm_free(gkn_mm mm) {
	free(mm->name);
	free(mm->score);
	free(mm);
}

gkn_mm gkn_mm_read(const char *filename) {
	char    *line = NULL;
	size_t   len = 0;
	ssize_t  read;
	gkn_pipe io  = gkn_pipe_open(filename, "r");
	double  *score = NULL;
	char     kmer[16];
	char     blah[256];
	int      size;
	double   p;

	while ((read = getline(&line, &len, io->stream)) != -1) {
		if (line[0] == '#') {
			assert(sscanf(line, "# MM %s %d", blah, &size) == 2);
			score = malloc(sizeof(double) * size);
		} else if (sscanf(line, "%s %lf", kmer, &p) == 2) {
			int idx = dna2dec(kmer, 0, strlen(kmer));
			if (idx == -1) gkn_exit("alphabet error in: %s", kmer);
			score[idx] = prob2score(p);
		}
	}
	gkn_pipe_close(io);
	if (line) free(line);

	gkn_mm model = malloc(sizeof(struct gkn_MM));
	model->name = malloc(strlen(filename)+1);
	strcpy(model->name, filename);
	model->k = strlen(kmer);
	model->size = size;
	model->score = score;

	return model;
}

double gkn_mm_score(const gkn_mm mm, const char *seq, int pos, int end) {
	double p = 0;
	if (pos < mm->k) pos = mm->k;
	for (int i = pos; i <= end; i++) {
		int idx = dna2dec(seq, i, mm->k);
		if (idx != -1) p += mm->score[idx];
	}
	return p;
}

double * gkn_mm_cache(const gkn_mm mm, const char *seq) {
	int len = strlen(seq);
	double *score = malloc(sizeof(double) * len);
	for (int i = 0; i < mm->k; i++) score[i] = 0;
	for (int i = mm->k; i < len; i++) {
		int idx = dna2dec(seq, i, mm->k);
		if (idx == -1) score[i] = score[i-1];
		else           score[i] = score[i-1] + mm->score[idx];
	}
	return score;
}

double gkn_mm_score_cache(const double *cache, int beg, int end) {
	return cache[end] - cache[beg -1];
}

// Length model

static double find_tail(double val, int x) {
	double lo = 0;
	double hi = 1000;
	double m;

	while (hi - lo > 1) {
		m = (hi + lo) / 2;
		double p = 1 / m;
		double f = pow(1-p, x-1) * p;
		//printf("%f %f %f %f\n", lo, hi, m, f);
		if (f < val) lo += (m - lo) / 2;
		else         hi -= (hi - m) / 2;
	}

	return m;
}

void gkn_len_free(gkn_len model) {
	free(model->name);
	free(model->score);
	free(model);
}

gkn_len gkn_len_read(const char *filename, int limit) {
	char    *line = NULL;
	size_t   len = 0;
	ssize_t  read;
	gkn_pipe io = gkn_pipe_open(filename, "r");
	double  *score = NULL;
	double   p;
	int      idx = 0;
	int      size;
	char     blah[64];

	// read probabilities
	while ((read = getline(&line, &len, io->stream)) != -1) {
		if (line[0] == '#') {
			assert(sscanf(line, "# LEN %s %d", blah, &size) == 2);
			score = malloc(sizeof(double) * size);
			
		} else if (sscanf(line, "%lf", &p) == 1) {
			score[idx] = p;
			idx++;
		}
	}

	gkn_pipe_close(io);
	if (line) free(line);
	
	gkn_len model = malloc(sizeof(struct gkn_LEN));
	model->name = malloc(strlen(filename)+1);
	strcpy(model->name, filename);
	model->score = score;
	model->size = size;
	model->limit = limit;
	model->tail = find_tail(score[size-1], size);
	
	// convert probabilities to scores
	for (int i = 0; i < size; i++) {
		double expect = (double) 1 / limit;
		score[i] = log(score[i]/expect) / log(2);
	}

	return model;
}

double gkn_len_score(const gkn_len len, int x) {
	assert(x >= 0);
	if (x >= len->size) {
		double p = 1 / len->tail;
		double q = pow(1-p, x-1) * p;
		double expect = (double)1 / len->limit;
		double s = log(q/expect) / log(2);
		return s;
	} else {
		return len->score[x];
	}
}

#endif
