#include "genomikon.h"
#include "isoform.h"

static char *usage = "\
isocounter - count all possible isoforms from sequences\n\n\
usage: iscounter <fasta file> [options]\n\
options:\n\
  --min_exon   <int>   minimum exon length [25]\n\
  --min_intron <int>   minimum intron length [35]\n\
  --max_splice <int>   maximum splices [3]\n\
  --flank      <int>   genomic flank lengths [99]\n\
  --introns    <file>  use introns from GFF file\n\
";

int main(int argc, char **argv) {
	char *file = NULL; // path to fasta file
	int   emin = 25;   // min exon size
	int   imin = 35;   // min intron size
	int   smax = 3;    // max splices
	int   gen  = 99;   // genomic flank (promoter, downstream)
	char *gff  = NULL; // use gff file for splice site coordinates

	gkn_pipe  io = NULL; // for reading fasta files
	gkn_fasta ff = NULL; // for individual fasta entries

	// CLI - setup
	gkn_set_program_name(argv[0]);
	gkn_register_option("--min_exon", 1);
	gkn_register_option("--min_intron", 1);
	gkn_register_option("--max_splice", 1);
	gkn_register_option("--flank",  1);
	gkn_register_option("--introns", 1);
	gkn_parse_options(&argc, argv);
	if (argc == 1) gkn_exit("%s", usage);

	// CLI - harvest
	file = argv[1];
	if (gkn_option("--min_exon"))   emin  = atoi(gkn_option("--min_exon"));
	if (gkn_option("--min_intron")) imin  = atoi(gkn_option("--min_intron"));
	if (gkn_option("--max_splice")) smax  = atoi(gkn_option("--max_splice"));
	if (gkn_option("--flank"))      gen   = atoi(gkn_option("--flank"));
	if (gkn_option("--introns"))    gff   = gkn_option("--introns");

	// find isoforms
	io = gkn_pipe_open(file, "r");
	while ((ff = gkn_fasta_read(io)) != NULL) {
		isozone iso = isoforms(ff->seq, emin, imin, smax, gen, gff, 1);
		for (int i = 0; i < strlen(ff->def); i++) {
			if (ff->def[i] == 32) ff->def[i] = 0;
		}
		printf("%s\t%d\t%d\t%d\t%d\t%d\n", ff->def+1, (int)strlen(ff->seq),
			iso->dons, iso->accs, iso->trials, iso->isoforms);
		gkn_fasta_free(ff);
		isozone_free(iso);
	}

	return 0;
}
