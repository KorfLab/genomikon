/******************************************************************************\
 sequence.h
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef GENOMIKON_SEQUENCE_H
#define GENOMIKON_SEQUENCE_H

#include "model.h"
#include "toolbox.h"

// Utilities
int    gkn_ntindex(const char *, int, int);
int    gkn_mem2idx(const char *, int, int);
int    gkn_str2idx(const char *);
char * gkn_idx2str(int, int);
char * gkn_revcomp(const char*);

// FASTA file
struct gkn_FASTA {
	int    length;
	char * def;
	char * seq;
};
typedef struct gkn_FASTA * gkn_fasta;
void	  gkn_fasta_free(gkn_fasta);
gkn_fasta gkn_fasta_new(const char *, const char *);
gkn_fasta gkn_fasta_read(gkn_pipe);
void	  gkn_fasta_write(FILE *, const gkn_fasta);
void	  gkn_fasta_set_line_length(int);

#endif
