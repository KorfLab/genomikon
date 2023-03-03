/******************************************************************************\
 pwm.h
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef PWM_H
#define PWM_H

#include "genomikon.h"


double pwm_distance(gkn_pwm, gkn_pwm);

struct motif_score_result {
	int seqs;
	int seqs_found;
	int motifs_found;
};

double score_motif(gkn_vec, gkn_pwm, int);
gkn_pwm background_model(gkn_vec, int);


#endif
