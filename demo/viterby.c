/*****************************************************************************\
 viterby.c
\*****************************************************************************/

#include "genomikon.h"

static char *usage = "\
viterby - gene prediction demo\n\n\
usage: viterby <hmm file> <fasta file> [options]\n\
options:\n\
  -verbose\n\
";

int main (int argc, char *argv[]) {

	// CLI
	gkn_set_program_name(argv[0]);
	gkn_register_option("-verbose", 0);
	gkn_parse_options(&argc, argv);

	if (argc != 3) gkn_exit(usage);

	gkn_pipe hmm_file = gkn_pipe_open(argv[1], "r");
	//gkn_hmm hmm = NULL;
	gkn_pipe_close(hmm_file);

	// process files
	gkn_pipe fasta_file = gkn_pipe_open(argv[2], "r");
	gkn_fasta fasta = NULL;
	while ((fasta = gkn_fasta_read(fasta_file)) != NULL) {


		// clean up

	}

	/* clean up */

	return 0;
}
