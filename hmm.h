/******************************************************************************\
 hmm.h
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef GENOMIKON_HMM_H
#define GENOMIKON_HMM_H

#include "toolbox.h"
#include "sequence.h"
#include "model.h"

// HMM state

struct gkn_STATE {
	char   *name;
	
	// ends
	double  init;  // initial probability
	double  term;  // terminal probability
	
	// transitions
	gkn_tvec adjs; // adjacent state
	gkn_fvec adjp; // adjacent probability
	
	// emission model
	gkn_mm emit;
	
	// length model
	gkn_len len;
};
typedef struct gkn_STATE * gkn_state;
void      gkn_state_free(gkn_state);
gkn_state gkn_state_read(FILE *);
void      gkn_state_write(FILE *, const gkn_state);

// HMM

struct gkn_HMM {
	char     *name;
	gkn_vec   states; // vector of gkn_state
	double  **trans;  // transition probabilities
};
typedef struct gkn_HMM * gkn_hmm;
void    gkn_hmm_free(gkn_hmm);
gkn_hmm gkn_hmm_read(FILE *);
void    gkn_hmm_write(FILE *, const gkn_hmm);

// HMM decoding trellis

struct gkn_TRELLIS {
	const char    *dna;
	const gkn_hmm  hmm;
	double       **vscore;
	int          **vtrace;
	int          **vjumps;
};
typedef struct gkn_TRELLIS * gkn_trellis;
void        gkn_trellis_free(gkn_trellis);
gkn_trellis gkn_trellis_new(const gkn_hmm, const char *);
gkn_ivec    gkn_trellis_viterbi(gkn_trellis);
gkn_vec     gkn_trellis_sviterbi(gkn_trellis);

#endif
