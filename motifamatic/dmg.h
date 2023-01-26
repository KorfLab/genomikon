/******************************************************************************\
 dmg.h
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef DMG_H
#define DMG_H

#include "genomikon.h"

struct discrete_motif_generator {
	int    sig;  // alphabet signature: 4, 5, 8, 9, 15, 19, 25
	char  *alph; // alphabet (copied, needs freeing)
	double P1;   // uppercase single value, eg. A = 0.970
	double P2;   // uppercase double value, eg. R = 0.480
	double P3;   // uppercase triple value, eg. B = 0.333
	double p1;   // lowercase single value, eg. a = 0.900
	double p2;   // lowercase double value, eg. r = 0.450
	double p3;   // lowercase triple value, eg. b = 0.300
	double N4;   // 0.25 each
	double Q1;   // remaining values e.g. C = 0.01, G = 0.01, T = 0.01
	double Q2;   // etc
	double Q3;
	double q1;
	double q2;
	double q3;
};
typedef struct discrete_motif_generator * dmgen;
void dmgen_free(dmgen);
dmgen dmgen_new(int);
dmgen dmgen_new_custom(int, double, double, double, double, double, double);

char * num2str(dmgen, int, int);
gkn_pwm num2pwm(dmgen, int, int);

#endif
