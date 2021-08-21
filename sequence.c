/******************************************************************************\
 sequence.c
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef GENOMIKON_SEQUENCE_C
#define GENOMIKON_SEQUENCE_C

#include "sequence.h"

// Utilities

int gkn_ntindex(const char *seq, int off, int k) {
	int idx = 0;
	for (int i = 0; i < k; i++) {
		switch (seq[off+i]) {
			case 'A': case 'a': idx += pow(4, (k -i -1)) * 0; break;
			case 'C': case 'c': idx += pow(4, (k -i -1)) * 1; break;
			case 'G': case 'g': idx += pow(4, (k -i -1)) * 2; break;
			case 'T': case 't': idx += pow(4, (k -i -1)) * 3; break;
			default: return -1;
		}
	}
	return idx;
}

char * gkn_revcomp (const char *seq) {
	int length = strlen(seq);
	char *str = gkn_malloc(length +1);
	str[strlen(seq)] = '\0';

	for (int i = 1; i <= length; i++) {
		switch (seq[i-1]) {
			case 'A': str[length -i] = 'T'; break;
			case 'a': str[length -i] = 't'; break;
			case 'C': str[length -i] = 'G'; break;
			case 'c': str[length -i] = 'g'; break;
			case 'G': str[length -i] = 'C'; break;
			case 'g': str[length -i] = 'c'; break;
			case 'T': str[length -i] = 'A'; break;
			case 't': str[length -i] = 'a'; break;
			case 'N': str[length -i] = 'N'; break;
			case 'n': str[length -i] = 'n'; break;
			case 'R': str[length -i] = 'Y'; break;
			case 'r': str[length -i] = 'y'; break;
			case 'Y': str[length -i] = 'R'; break;
			case 'y': str[length -i] = 'r'; break;
			case 'W': str[length -i] = 'S'; break;
			case 'w': str[length -i] = 's'; break;
			case 'S': str[length -i] = 'W'; break;
			case 's': str[length -i] = 'w'; break;
			case 'K': str[length -i] = 'M'; break;
			case 'k': str[length -i] = 'm'; break;
			case 'M': str[length -i] = 'K'; break;
			case 'm': str[length -i] = 'k'; break;
			case 'B': str[length -i] = 'V'; break;
			case 'b': str[length -i] = 'v'; break;
			case 'D': str[length -i] = 'H'; break;
			case 'd': str[length -i] = 'h'; break;
			case 'H': str[length -i] = 'D'; break;
			case 'h': str[length -i] = 'd'; break;
			case 'V': str[length -i] = 'B'; break;
			case 'v': str[length -i] = 'b'; break;
			default:  gkn_exit("alphabet error %c", seq[i-1]);
		}
	}

	return str;
}

// FASTA file

void gkn_fasta_free(gkn_fasta ff) {
	free(ff->def);
	free(ff->seq);
	free(ff);
}

gkn_fasta gkn_fasta_new (const char *def, const char *seq) {
	gkn_fasta ff = gkn_malloc(sizeof(struct gkn_FASTA));
	ff->def	= gkn_malloc(strlen(def) +1);
	ff->seq	= gkn_malloc(strlen(seq) +1);
	ff->length = strlen(seq);
	strcpy(ff->def, def);
	strcpy(ff->seq, seq);
	return ff;
}

gkn_fasta gkn_fasta_read(gkn_pipe io) {

	// check for fasta header
	char c = fgetc(io->stream);
	if (c == EOF || (unsigned char)c == 255) return NULL;
	if (c != '>') gkn_exit("fasta? %c %d", c, (int)c);
	ungetc(c, io->stream);

	// def
	char *def = gkn_readline(io);

	// seq
	gkn_vec lines = gkn_vec_new();
	while (1) {
		char c = fgetc(io->stream);
		if (c == EOF) break;
		if (c == '>') {
			ungetc(c, io->stream);
			break;
		}
		char *line = gkn_readline(io);
		if (line == NULL) break;
		gkn_vec_push(lines, line);
	}
	
	int letters = 0;
	for (int i = 0; i < lines->size; i++) {
		char *line = lines->elem[i];
		letters += strlen(line);
	}
	
	char *seq = malloc(letters + 1);
	int off = 0;
	for (int i = 0; i < lines->size; i++) {
		char *line = lines->elem[i];
		strcpy(seq + off, line);
		off += strlen(line);
	}
	
	// clean up
	for (int i = 0; i < lines->size; i++) {
		free(lines->elem[i]);
	}
	gkn_vec_free(lines);

	// return object
	gkn_fasta ff = malloc(sizeof(struct gkn_FASTA));
	ff->def = def;
	ff->seq = seq;
	ff->length = 1; //strlen(seq);
	
	return ff;
}

static int FASTA_LINE_LENGTH = 80;

void gkn_fasta_write(FILE *stream, const gkn_fasta ff) {
	if (ff->def[0] != '>') fprintf(stream, ">");
	fprintf(stream, "%s", ff->def);
	if (ff->def[strlen(ff->def) -1] != '\n') fprintf(stream, "\n");

	for (int i = 0; i < ff->length; i++) {
		fputc(ff->seq[i], stream);
		if ((i+1) % FASTA_LINE_LENGTH == 0) fprintf(stream, "\n");
	}

	fprintf(stream, "\n");
}

void gkn_fasta_set_line_length (int length) {
	FASTA_LINE_LENGTH = length;
}

#endif
