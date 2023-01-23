/******************************************************************************\
 model.c
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef GENOMIKON_MODEL_C
#define GENOMIKON_MODEL_C

#include "model.h"

double gkn_p2s(double p) {
	assert(p >= 0 && p <= 1);
	if (p == 0) return -100; // umm...
	return log(p/0.25);
}

double gkn_sum2(double a, double b) {
	if (fabs(a - b) > 20) return (a > b) ? a : b;
	return (log(1 + pow(2.7182818, b-a)) + a);
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

gkn_pwm gkn_pwm_read(gkn_pipe io) {
	char    *line;
	char     name[256];
	int      size;
	double **score = NULL;
	double   a, c, g, t;
	int      row = 0;

	while ((line = gkn_readline(io)) != NULL) {
		if (line[0] == '%') {
			assert(sscanf(line, "%% PWM %s %d", name, &size) == 2);
			score = malloc(sizeof(double*) * size);
			for (int i = 0; i < size; i++) {
				score[i] = malloc(sizeof(double) * 4);
			}
			free(line);
		} else if (sscanf(line, "%lf %lf %lf %lf", &a, &c, &g, &t) == 4) {
			score[row][0] = gkn_p2s(a);
			score[row][1] = gkn_p2s(c);
			score[row][2] = gkn_p2s(g);
			score[row][3] = gkn_p2s(t);
			row++;
			free(line);
		} else {
			free(line);
		}
	}

	gkn_pwm model = malloc(sizeof(struct gkn_PWM));
	model->name = malloc(strlen(name)+1);
	strcpy(model->name, name);
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

gkn_mm gkn_mm_read(gkn_pipe io) {
	char    *line = NULL;
	double  *score = NULL;
	char     kmer[16];
	char     name[256];
	int      size;
	double   p;

	while ((line = gkn_readline(io)) != NULL) {
		if (line[0] == '%') {
			assert(sscanf(line, "%% MM %s %d", name, &size) == 2);
			score = malloc(sizeof(double) * size);
			free(line);
		} else if (sscanf(line, "%s %lf", kmer, &p) == 2) {
			int idx = gkn_str2idx(kmer);
			if (idx == -1) gkn_exit("alphabet error in: %s", kmer);
			score[idx] = gkn_p2s(p);
			free(line);
		} else {
			free(line);
		}
	}

	gkn_mm model = malloc(sizeof(struct gkn_MM));
	model->name = malloc(strlen(name)+1);
	strcpy(model->name, name);
	model->k = strlen(kmer);
	model->size = size;
	model->score = score;

	return model;
}

double gkn_mm_score(const gkn_mm mm, const char *seq, int pos, int end) {
	double p = 0;
	if (pos < mm->k) pos = mm->k;
	for (int i = pos; i < end - mm->k +2; i++) {
		int idx = gkn_mem2idx(seq, i, mm->k);
		if (idx != -1) p += mm->score[idx];
	}
	return p;
}

double * gkn_mm_cache(const gkn_mm mm, const char *seq) {
	int len = strlen(seq);
	double *score = malloc(sizeof(double) * len);
	for (int i = 0; i < mm->k; i++) score[i] = 0;
	for (int i = mm->k; i < len; i++) {
		int idx = gkn_mem2idx(seq, i, mm->k);
		if (idx == -1) score[i] = score[i-1];
		else           score[i] = score[i-1] + mm->score[idx];
	}
	return score;
}

double gkn_mm_score_cache(const gkn_mm mm, const double *v, int beg, int end) {
	return v[end - mm->k +1] - v[beg -1];
}

// Length model

static double find_tail(double val, int x) {
	double lo = 0;
	double hi = 1000; // maybe param
	double m;

	while (hi - lo > 1) {
		m = (hi + lo) / 2;
		double p = 1 / m;
		double f = pow(1-p, x-1) * p;
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

gkn_len gkn_len_read(gkn_pipe io) {
	char    *line = NULL;
	double  *score = NULL;
	double   p;
	int      idx = 0;
	int      size;
	char     name[64];

	// read probabilities
	while ((line = gkn_readline(io)) != NULL) {
		if (line[0] == '%') {
			assert(sscanf(line, "%% LEN %s %d", name, &size) == 2);
			score = malloc(sizeof(double) * size);
			free(line);
		} else if (sscanf(line, "%lf", &p) == 1) {
			score[idx] = p;
			idx++;
			free(line);
		} else {
			free(line);
		}
	}

	gkn_len model = malloc(sizeof(struct gkn_LEN));
	model->name = malloc(strlen(name)+1);
	strcpy(model->name, name);
	model->score = score;
	model->size = size;
	model->tail = find_tail(score[size-1], size);

	// convert probabilities to scores
	double expect = (double) 1 / model->size;
	for (int i = 0; i < size; i++) {
		score[i] = log(score[i]/expect) / log(2); // divide by zero?
	}

	return model;
}

double gkn_len_score(const gkn_len len, int x) {
	assert(x >= 0);
	if (x >= len->size) {
		double p = 1 / len->tail;
		double q = pow(1-p, x-1) * p;
		double expect = (double)1 / len->size;
		double s = log(q/expect) / log(2);
		return s;
	} else {
		return len->score[x];
	}
}

#endif
