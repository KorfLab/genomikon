/******************************************************************************\
 dmg.h
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef DMG_H
#define DMG_H

#include "genomikon.h"

struct discrete_motif_generator {
	int    aid;      // alphabet identifier
	char   alph[32]; // alphabet
	double P1;       // uppercase single value
	double P2;       // uppercase double value
	double P3;       // uppercase triple value
	double p1;       // lowercase single value
	double p2;       // lowercase double value
	double p3;       // lowercase triple value
	double N4;       // 0.25 each
	double Q1;       // remaining single values (1-P)/3
	double Q2;       // reamining double values (1-2P)/2
	double Q3;       // reamining triple value (1-3p)
	double q1;       // etc for lowercase
	double q2;
	double q3;
};
typedef struct discrete_motif_generator * dmgen;
void dmgen_free(dmgen);
dmgen dmgen_new(int);
dmgen dmgen_new_custom(int, double, double, double, double, double, double);

char * num2str(dmgen, int, int);
gkn_pwm num2pwm(dmgen, int, int);

struct motif_score_result {
	int seqs;
	int seqs_found;
	int motifs_found;
};

double score_motif(gkn_vec, gkn_pwm, int);
gkn_pwm background_model(gkn_vec, int);

#endif
