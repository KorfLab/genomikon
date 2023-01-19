/******************************************************************************\
 model.h
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef GENOMIKON_MODEL_H
#define GENOMIKON_MODEL_H

#include "sequence.h"
#include "toolbox.h"

// Utilities

double gkn_p2s(double);          // probability to score
double gkn_sum2(double, double); // sum 2 values in log space

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
double   gkn_mm_score_cache(const gkn_mm, const double *, int, int);

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

// Hidden Markov model

struct gkn_HMM {
	char   *name;    // something descriptive
	gkn_vec states;  // number of states in the model
};
typedef struct gkn_HMM * gkn_hmm;
void    gkn_hmm_free(gkn_hmm);
gkn_hmm gkn_hmm_read(gkn_pipe);

struct gkn_STATE {
	char   *name;        // exon, intron, whatever
	double  init;        // initial probability
	double  term;        // terminal probability
	int     transitions; // number of states this connects out to
	int     emissions;   // number of emission probs (4**n)
	int     durations;   // number of optional duration probabilities

	gkn_map  trans; // map of transition probabilities to other states
	gkn_mm   emits; // emission probabilities for Nth order MM
	gkn_len  durs;  // optional length model
};
typedef struct gkn_STATE * gkn_state;
void      gkn_state_free(gkn_state);
gkn_state gkn_state_read(FILE *);
void      gkn_state_write(FILE *, const gkn_state);

struct gkn_DECODER {
	char    *seq;        // sequence as read
	gkn_hmm  hmm;        // hmm as read
	double **transition; // transition probability matrix
	double **vscore;     // viterbi score
	int    **vtrace;     // viterbi trace
	int    **vjumps;     // viterbi jumps
	double **fscore;     // forward score
	double **bscore;     // backward score
};
typedef struct gkn_DECODER * gkn_decoder;
void        gkn_decoder_free(gkn_decoder);
gkn_decoder gkn_decoder_new(const gkn_hmm, const char*);
void        gkn_decoder_viteri(gkn_decoder);
void        gkn_decoder_forward(gkn_decoder);
void        gkn_decoder_backward(gkn_decoder);

#endif
