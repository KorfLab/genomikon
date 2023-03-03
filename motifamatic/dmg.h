/******************************************************************************\
 dmg.h
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef DMG_H
#define DMG_H

#include "genomikon.h"

void init_DNTP();
void set_DNTP(double, double, double, double, double, double);

char * get_alphabet(int);

/*
struct discrete_motif_generator {
	int    aid;      // alphabet identifier
	char   alph[32]; // alphabet
};
typedef struct discrete_motif_generator * dmgen;
void dmgen_free(dmgen);
dmgen dmgen_new(int);
*/

char*   num2str(const char*, int, int);
gkn_pwm num2pwm(const char*, int, int);
char*   pwm2str(gkn_pwm);
gkn_pwm str2pwm(const char*, const char*);

struct motif_score_result {
	int seqs;
	int seqs_found;
	int motifs_found;
};

double score_motif(gkn_vec, gkn_pwm, int);
gkn_pwm background_model(gkn_vec, int);

#endif
