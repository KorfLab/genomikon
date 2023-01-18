/*****************************************************************************\
 viterby.c
\*****************************************************************************/

#include "genomikon.h"

static char *usage = "\
viterby - a typical HMM Viterbi decoder\n\n\
usage: viterby <hmm> <fasta>\n\
";




int main(int argc, char **argv) {
	if (argc !=3) gkn_exit("%s", usage);
	
	// read the hmm
	
	// read and decode sequences
	gkn_pipe io = gkn_pipe_open(argv[2], "r");
	gkn_fasta in;
	while ((in = gkn_fasta_read(io)) != NULL) {
		printf("name: %s\n", in->def);
		printf("seq: %s\n", in->seq);
		gkn_fasta_free(in);
	}
	gkn_pipe_close(io);

}

