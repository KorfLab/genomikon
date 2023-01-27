/******************************************************************************\
 dmg.c
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef DMG_C
#define DMG_C

#include "genomikon.h"
#include "dmg.h"

static const char *S4  = "ACGT";
static const char *S5  = "ACGTN";
static const char *S9  = "ACGTacgtN";
static const char *S11 = "ACGTRYMKWSN";
static const char *S15 = "ACGTRYMKWSBDHVN";
static const char *S19 = "ACGTRYMKWSBDHVacgtN";
static const char *S25 = "ACGTRYMKWSBDHVacgtrymkwsN";
static const char *S29 = "ACGTRYMKWSBDHVacgtrymkwsbdhvN";

void dmgen_free(dmgen dmg) {
	free(dmg->alph);
	free(dmg);
}

dmgen dmgen_new(int sig) {
	return dmgen_new_custom(sig, 0.997, 0.498, 0.333, 0.8, 0.4, 0.3);
}

dmgen dmgen_new_custom(int sig,
		double P1, double P2, double P3,
		double p1, double p2, double p3) {
	dmgen dmg = malloc(sizeof(struct discrete_motif_generator));
	dmg->sig = sig;
	dmg->alph = malloc(sig+1);
	switch (sig) {
		case 4:  strcpy(dmg->alph, S4);  break;
		case 5:  strcpy(dmg->alph, S5);  break;
		case 9:  strcpy(dmg->alph, S9);  break;
		case 11: strcpy(dmg->alph, S11); break;
		case 15: strcpy(dmg->alph, S15); break;
		case 19: strcpy(dmg->alph, S19); break;
		case 25: strcpy(dmg->alph, S25); break;
		case 29: strcpy(dmg->alph, S29); break;
		default: gkn_exit("illegal sig: %d\n", sig);
	}

	dmg->P1 = P1;
	dmg->P2 = P2;
	dmg->P3 = P3;
	dmg->p1 = p1;
	dmg->p2 = p2;
	dmg->p3 = p3;
	dmg->Q1 = (1 - P1) / 3;
	dmg->Q2 = (1 - 2*P2) / 2;
	dmg->Q3 = (1 - 3*P3) / 3;
	dmg->q1 = (1 - p1) / 3;
	dmg->q2 = (1 - 2*p2) / 2;
	dmg->q3 = (1 - 3*p3) / 3;
	dmg->N4 = 0.25;

	return dmg;
}

char * num2str(dmgen dmg, int num, int len) {
	char *str = malloc(len + 1);
	char *src = dmg->alph;

	str[len] = '\0';
	for (int i = 0; i < len; i++) {
		int max = pow(dmg->sig, (len-i-1));
		int r = 0;
		if (num > max -1) {
			r = num / max;
			num -= r * max;
		}
		str[i] = src[r];
	}
	return str;
}

gkn_pwm num2pwm(dmgen dmg, int num, int len) {
	double ** pwm = malloc(sizeof(double*) * len);
	for (int i = 0; i < len; i++) {
		pwm[i] = malloc(sizeof(double) * 4);
	}
	double P1 = dmg->P1;
	double P2 = dmg->P2;
	double P3 = dmg->P3;
	double p1 = dmg->p1;
	double p2 = dmg->p2;
	double p3 = dmg->p3;
	double Q1 = dmg->Q1;
	double Q2 = dmg->Q2;
	double Q3 = dmg->Q3;
	double q1 = dmg->q1;
	double q2 = dmg->q2;
	double q3 = dmg->q3;
	double N4 = dmg->N4;

	char *str = num2str(dmg, num, len);
	for (int i = 0; i < len; i++) {
		switch (str[i]) {
		//                  A             C             G             T
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

	gkn_pwm model = malloc(sizeof(struct gkn_PWM));
	char name[64];
	sprintf(name, "%s-%d-%s", dmg->alph, num, str);
	model->name = malloc(strlen(name)+1);
	strcpy(model->name, name);
	model->size = len;
	model->score = pwm;

	free(str);

	return model;
}

#endif
