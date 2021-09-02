/******************************************************************************\
 model.h
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef GENOMIKON_MODEL_H
#define GENOMIKON_MODEL_H

#include "sequence.h"
#include "toolbox.h"

// Utilities

double gkn_p2s(double);

// Position Weight Matrix

struct gkn_PWM {
	char    *name;   // acceptor, donor
	int      size;   // eg. 6
	double **score;  // score[pos][nt]
};
typedef struct gkn_PWM * gkn_pwm;
void    gkn_pwm_free(gkn_pwm);
gkn_pwm gkn_pwm_read(gkn_pipe);
double  gkn_pwm_score(const gkn_pwm, const char *, int);

// Markov model

struct gkn_MM {
	char   *name;   // exon, intron
	int     k;      // kmer size
	int     size;   // size of array
	double *score;  // score[dna2dec()] = value
};
typedef struct gkn_MM * gkn_mm;
void     gkn_mm_free(gkn_mm);
gkn_mm   gkn_mm_read(gkn_pipe);
double   gkn_mm_score(const gkn_mm, const char *, int, int);
double * gkn_mm_cache(const gkn_mm, const char *);
double   gkn_mm_score_cache(const double *, int, int);

// Length model

struct gkn_LEN {
	char   *name;   // exon, intron, actually unused
	int     size;   // length of defined region
//	int     limit;  // maximum length for scoring
	double *score;  // values for the defined region
	double  tail;   // mean of geometric tail
};
typedef struct gkn_LEN * gkn_len;
void    gkn_len_free(gkn_len);
gkn_len gkn_len_read(gkn_pipe/*, int*/);
double  gkn_len_score(const gkn_len, int);

#endif
