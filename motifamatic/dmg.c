/******************************************************************************\
 dmg.c
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef DMG_C
#define DMG_C

#include "genomikon.h"
#include "dmg.h"

// discrete nucleotide
static double DNTP[128][4]; // probabilities
static const char *DNTS = "ACGTRYMKWSBDHVNacgtrymkwsbdhvn"; // symbols
static const char *DNTA[9] = {
	"ACGT",
	"ACGTN",
	"ACGTacgtN",
	"ACGTRYMKWSN",
	"ACGTRYMKWSNacgt",
	"ACGTRYMKWSBDHVN",
	"ACGTRYMKWSBDHVacgtN",
	"ACGTRYMKWSBDHVacgtrymkwsN",
	"ACGTRYMKWSBDHVacgtrymkwsbdhvN"
};

void init_DNTP(void) {
	set_DNTP(0.97, 0.49, 0.33, 0.70, 0.40, 0.30);
}

void set_DNTP(double P1, double P2, double P3, double p1, double p2, double p3) {
	double Q1 = (1 - P1) / 3;
	double Q2 = (1 - 2 * P2) / 2;
	double Q3 = (1 - 3 * P3) / 3;
	double q1 = (1 - p1) / 3;
	double q2 = (1 - 2 * p2) / 2;
	double q3 = (1 - 3 * p3) / 3;
	double N4 = 0.25;
	double n4 = 0.25;

	for (int i = 0; i < 128; i++) {
		for (int j = 0; j < 4; j++) {
			DNTP[i][j] = -1;
		}
	}

	double prob[30][4] = {
		{P1, Q1, Q1, Q1}, // A
		{Q1, P1, Q1, Q1}, // C
		{Q1, Q1, P1, Q1}, // G
		{Q1, Q1, Q1, P1}, // T
		{P2, Q2, P2, Q2}, // R
		{Q2, P2, Q2, P2}, // Y
		{P2, P2, Q2, Q2}, // M
		{Q2, Q2, P2, P2}, // K
		{P2, Q2, Q2, P2}, // W
		{Q2, P2, P2, Q2}, // S
		{Q3, P3, P3, P3}, // B
		{P3, Q3, P3, P3}, // D
		{P3, P3, Q3, P3}, // H
		{P3, P3, P3, Q3}, // V
		{N4, N4, N4, N4}, // N
		{p1, q1, q1, q1}, // a
		{q1, p1, q1, q1}, // c
		{q1, q1, p1, q1}, // g
		{q1, q1, q1, p1}, // t
		{p2, q2, p2, q2}, // r
		{q2, p2, q2, p2}, // y
		{p2, p2, q2, q2}, // m
		{q2, q2, p2, p2}, // k
		{p2, q2, q2, p2}, // w
		{q2, p2, p2, q2}, // s
		{q3, p3, p3, p3}, // b
		{p3, q3, p3, p3}, // d
		{p3, p3, q3, p3}, // h
		{p3, p3, p3, q3}, // v
		{n4, n4, n4, n4}  // n
	};

	for (int i = 0; i < 30; i++) {
		for (int j = 0; j < 4; j++) {
			DNTP[(int)DNTS[i]][j] = prob[i][j];
		}
	}
}

char * get_alphabet(int n) {
	assert(n >= 0 && n <= 8);
	char *str = malloc(strlen(DNTA[n]) + 1);
	strcpy(str, DNTA[n]);
	return str;
}

char * num2str(const char *alph, int num, int len) {
	char *str = malloc(len + 1);
	//char *src = alph;
	int n = strlen(alph);

	str[len] = '\0';
	for (int i = 0; i < len; i++) {
		int max = pow(n, (len-i-1));
		int r = 0;
		if (num > max -1) {
			r = num / max;
			num -= r * max;
		}
		str[i] = alph[r];
	}
	return str;
}

gkn_pwm num2pwm(const char *alph, int num, int len) {
	double ** pwm = malloc(sizeof(double*) * len);
	for (int i = 0; i < len; i++) {
		pwm[i] = malloc(sizeof(double) * 4);
	}

	char *str = num2str(alph, num, len);
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < 4; j++) {
			pwm[i][j] = gkn_p2s(DNTP[(int)str[i]][j]);
		}
	}

	gkn_pwm model = malloc(sizeof(struct gkn_PWM));
	char name[64];
	sprintf(name, "%s-%d-%s", alph, num, str);
	model->name = malloc(strlen(name)+1);
	strcpy(model->name, name);
	model->size = len;
	model->score = pwm;

	free(str);

	return model;
}

char * pwm2str(gkn_pwm pwm) {
	char *s = malloc(pwm->size + 1);
	return s;
}


gkn_pwm str2pwm(const char *s, const char *name) {
	int len = strlen(s);

	double ** pwm = malloc(sizeof(double*) * len);
	for (int i = 0; i < len; i++) {
		pwm[i] = malloc(sizeof(double) * 4);
	}

	// unfinished - need the actual values

	gkn_pwm model = malloc(sizeof(struct gkn_PWM));
	model->name = malloc(strlen(name)+1);
	strcpy(model->name, name);
	model->size = len;
	model->score = pwm;

	return model;
}


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
