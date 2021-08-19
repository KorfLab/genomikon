/******************************************************************************\
 sequence.h
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef GENOMIKON_SEQUENCE_H
#define GENOMIKON_SEQUENCE_H

#include "toolbox.h"

// Utilities
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
gkn_fasta gkn_fasta_read(FILE *);
void	  gkn_fasta_write(FILE *, const gkn_fasta);
void	  gkn_fasta_set_line_length(int);

#endif
