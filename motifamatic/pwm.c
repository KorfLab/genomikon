/******************************************************************************\
 pwm.c
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef PWM_C
#define PWM_C

#include "genomikon.h"
#include "pwm.h"

double pwm_distance(gkn_pwm m1, gkn_pwm m2) {
	return 3.14;
}


// motif comparison

/*

manhattan, euclidean, kl distances

sliding motif comparison

gapped motif comparision?



*/



// motif-finding functions

static int count_motifs(const char *seq, gkn_pwm pwm, double t) {
	int n = 0;
	for (int i = 0; i < strlen(seq) - pwm->size + 1; i++) {
		if (gkn_pwm_score(pwm, seq, i) > t) n++;
	}
	return n;
}

double score_motif(gkn_vec seqs, gkn_pwm pwm, int mod) {
	int swm = 0; // sequences with motifs
	int tot = 0; // total motifs found
	for (int i = 0; i < seqs->size; i++) {
		int n = count_motifs(seqs->elem[i], pwm, 5);
		if (n > 0) swm++;
		tot += n;
	}

	// convert swm, tot to prob

	return swm;
}

gkn_pwm background_model(gkn_vec seqs, int len) {
	int a = 0;
	int c = 0;
	int g = 0;
	int t = 0;
	int total = 0;
	for (int i = 0; i < seqs->size; i++) {
		char *seq = seqs->elem[i];
		int seqlen = strlen(seq);
		total += seqlen;
		for (int j = 0; j < seqlen; j++) {
			switch (seq[j]) {
				case 'A': a++; break;
				case 'C': c++; break;
				case 'G': g++; break;
				case 'T': t++; break;
			}
		}
	}

	double ** pwm = malloc(sizeof(double*) * len);
	for (int i = 0; i < len; i++) {
		pwm[i] = malloc(sizeof(double) * 4);
	}
	for (int i = 0; i < len; i++) {
		pwm[i][0] = gkn_p2s((double)a/total);
		pwm[i][1] = gkn_p2s((double)c/total);
		pwm[i][2] = gkn_p2s((double)g/total);
		pwm[i][3] = gkn_p2s((double)t/total);
	}

	gkn_pwm model = malloc(sizeof(struct gkn_PWM));
	char *name = "background";
	model->name = malloc(strlen(name)+1);
	strcpy(model->name, name);
	model->size = len;
	model->score = pwm;

	return model;
}

#endif
