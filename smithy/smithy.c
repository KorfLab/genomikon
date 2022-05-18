/*****************************************************************************\
 smithy.c
\*****************************************************************************/

#include "genomikon.h"

static char *usage = "\
smithy - smith-waterman\n\n\
usage: smithy <string1> <string2>\n\
";

int main(int argc, char **argv) {
	if (argc == 1) gkn_exit("%s", usage);

	gkn_smat b62 = gkn_smat_blosum(62);
	gkn_hsp hsp = gkn_sw(argv[1], argv[2], b62);
	printf("%d %d..%d %d..%d\n", hsp->score, hsp->beg1+1, hsp->end1+1,
		hsp->beg2+1, hsp->end2+1);
	printf("%s\n%s\n%s\n", hsp->s1, hsp->as, hsp->s2);

	return 0;
}
