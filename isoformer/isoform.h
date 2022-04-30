/******************************************************************************\
 isoform.h
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef ISOFORM_H
#define ISOFORM_H

#include "genomikon.h"

struct isoform_zone {
	int     dons;
	int     accs;
	int     trials;
	int     isoforms;
	gkn_vec mRNAs;
};
typedef struct isoform_zone * isozone;
void    isozone_free(isozone);
isozone isoforms(const char*, int, int, int, int, const char*, int);

#endif
