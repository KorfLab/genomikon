/******************************************************************************\
 align.h
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef GENOMIKON_ALIGN_H
#define GENOMIKON_ALIGN_H

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "toolbox.h"

// scoring matrices

struct gkn_SMAT {
	char name[64];
	int  gap;
	int  score[25][25];
};
typedef struct gkn_SMAT * gkn_smat;
void     gkn_smat_free(gkn_smat);
gkn_smat gkn_smat_blosum(int);
gkn_smat gkn_smat_mng(int, int, int);

// high scoring pairs

struct gkn_HSP {
	int   beg1;
	int   end1;
	int   beg2;
	int   end2;
	int   length;
	int   score;
	char *s1;
	char *s2;
	char *as;
};
typedef struct gkn_HSP * gkn_hsp;
void    gkn_hsp_free (gkn_hsp);
gkn_hsp gkn_hsp_new (void);

// alignment algorithms

gkn_hsp gkn_sw(const char *, const char *, gkn_smat);

#endif
