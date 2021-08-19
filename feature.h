/******************************************************************************\
 feature.h
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef GENOMIKON_FEATURE_H
#define GENOMIKON_FEATURE_H

#include "toolbox.h"

// gff
struct gkn_GFF {
	char  *name;    // chromosome or sequence name
	char  *source;  // whatever (ik for stuff here?)
	char  *type;    // should be SO compliant but... 
	char  *group;   // possibly not used
	int    beg;     // 0-based internally
	int    end;     // 0-based internally
	double score;   // . or double
	char   strand;  // {.+-}
	char   phase;   // {.012}
};
typedef struct gkn_GFF * gkn_gff;
gkn_gff gkn_gff_read(FILE *);
gkn_gff gkn_gff_new(void);
void    gkn_gff_free(gkn_gff);

// features
struct gkn_FEAT {
	const char *seq;    // parent sequence
	int         beg;    // 0-based
	int         end;    // 0-based
	double      score;  // defaults to 0, set manually
};
typedef struct gkn_FEAT * gkn_feat;
gkn_feat gkn_feat_new(const char *, int, int);
void     gkn_feat_free(gkn_feat);
char *   gkn_feat_seq(const gkn_feat);

// mRNA
struct gkn_MRNA {
	const char *seq;     // parent sequence
	int         beg;     // 0-based, start of 5'UTR
	int         end;     // 0-based, end at poly-A tail
	gkn_vec     exons;   // gkn_feat
	gkn_vec     introns; // gkn_feat
	int         atg;     // CDS start, manually set, default -1
	double      score;   // defaults to 0, set manually
};
typedef struct gkn_MRNA * gkn_mRNA;
gkn_mRNA gkn_mRNA_new(const char *, int, int);
gkn_mRNA gkn_mRNA_read(const char *, const char *);
void     gkn_mRNA_free(gkn_mRNA);

#endif
