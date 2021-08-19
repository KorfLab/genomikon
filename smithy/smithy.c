/*****************************************************************************\
 smithy.c
\*****************************************************************************/

#include "toolbox.h"
#include "sequence.h"
#include "align.h"

void stuff (int m, int n, int g, int a) {

}

static char *usage = "\
sw - smith-waterman for demonstration purposes\n\n\
usage: sw <file1> <file2> [options]\n\
options:\n\
  -m <int>  match score [1]\n\
  -n <int>  mismatch score [-1]\n\
  -g <int>  gap score [-2]\n\
  -a <int>    algorithm [1]";

int main(int argc, char **argv) {
	char *file1 = NULL; // path to fasta file
	char *file2 = NULL; // path to fasta file
	int   m =  1;       // match score, positive
	float n = -1;       // mismatch score, negative
	int   g = -2;       // gap score, negative
	int   a =  1;       // algorithm
	gkn_pipe io1, io2;
	gkn_fasta ff1, ff2;
	
	// Command Line Interface
	gkn_set_program_name(argv[0]);
	gkn_register_option("-m", 1);
	gkn_register_option("-n", 1);
	gkn_register_option("-h", 1);
	gkn_register_option("-a", 1);
	gkn_parse_options(&argc, argv);
	
	if (argc == 1) gkn_exit("%s", usage);
	
	file1 = argv[1];
	file2 = argv[2];
	if (gkn_option("-m")) m = atoi(gkn_option("-m"));
	if (gkn_option("-n")) n = atoi(gkn_option("-n"));
	if (gkn_option("-g")) g = atoi(gkn_option("-g"));
	if (gkn_option("-a")) a = atoi(gkn_option("-a"));
	
	// main loop
	io1 = gkn_pipe_open(file1, "r");
	io2 = gkn_pipe_open(file2, "r");
	while ((ff1 = gkn_fasta_read(io1->stream)) != NULL) {
		while ((ff2 = gkn_fasta_read(io2->stream)) != NULL) {
		
			stuff(m, n, g, a);
		
			gkn_fasta_free(ff1);
			gkn_fasta_free(ff2);
		}
	}

	return 0;
}
