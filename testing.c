/*****************************************************************************\
 testing.c
\*****************************************************************************/

#include <stdio.h>
#include <assert.h>

#include "align.h"
#include "feature.h"
#include "model.h"
#include "sequence.h"
#include "toolbox.h"

static int COUNT = 100;
void test_math(void);
void test_seq(void);
void test_vec(int);
void test_ivec(int);
void test_fvec(int);
void test_tvec(int);
void test_vec(int);
void test_tmap(int);
void test_xtree(int);
void test_map(int);
void test_feat(int);
void test_smat(int);
void test_sw(int);
void test_pipe(int, const char *);
void test_read(int, const char *);
void test_fasta(int, const char *);
void test_gff(int, const char *);
void test_pwm(int, const char *);
void test_mm(int, const char *);
void test_len(int, const char*);

void test_this(void);

static char usage[] = "\
usage: testing [options]\n\
options:\n\
  -math\n\
  -vec -ivec -fvec -tvec -map -tmap -xtree\n\
  -feat -smat -sw\n\
  -pipe  <file>\n\
  -fasta <file>\n\
  -gff   <file>\n\
  -pwm   <file>\n\
  -mm    <file>\n\
  -len   <file>\n\
  -count <int> [100]\n\
";

int main(int argc, char ** argv) {

	/* preamble */
	if (argc == 1) {
		fprintf(stderr, "%s", usage);
		exit(1);
	}

	/* options */
	gkn_set_program_name(argv[0]);
	gkn_register_option("-count", 1);
	gkn_register_option("-math",  0);
	gkn_register_option("-seq",   0);
	gkn_register_option("-vec",   0);
	gkn_register_option("-ivec",  0);
	gkn_register_option("-fvec",  0);
	gkn_register_option("-tvec",  0);
	gkn_register_option("-vec",   0);
	gkn_register_option("-tmap",  0);
	gkn_register_option("-map",   0);
	gkn_register_option("-xtree", 0);
	gkn_register_option("-feat",  0);
	gkn_register_option("-smat",  0);
	gkn_register_option("-sw",    0);
	gkn_register_option("-pipe",  1);
	gkn_register_option("-read",  1);
	gkn_register_option("-fasta", 1);
	gkn_register_option("-gff",   1);
	gkn_register_option("-pwm",   1);
	gkn_register_option("-mm",    1);
	gkn_register_option("-len",   1);
	gkn_register_option("-this",  0);
	gkn_parse_options(&argc, argv);

	/* control */
	if (gkn_option("-count")) COUNT = atoi(gkn_option("-count"));
	int update = (int)((double)COUNT/50);
	if (update < 2) update = 2;

	if (gkn_option("-math"))  test_math();
	if (gkn_option("-seq"))   test_seq();
	if (gkn_option("-vec"))   test_vec(update);
	if (gkn_option("-ivec"))  test_ivec(update);
	if (gkn_option("-fvec"))  test_fvec(update);
	if (gkn_option("-tvec"))  test_tvec(update);
	if (gkn_option("-tmap"))  test_tmap(update);
	if (gkn_option("-map"))	  test_map(update);
	if (gkn_option("-xtree")) test_xtree(update);
	if (gkn_option("-feat"))  test_feat(update);
	if (gkn_option("-smat"))  test_smat(update);
	if (gkn_option("-sw"))    test_sw(update);
	if (gkn_option("-pipe"))  test_pipe(update,  gkn_option("-pipe"));
	if (gkn_option("-read"))  test_read(update,  gkn_option("-read"));
	if (gkn_option("-fasta")) test_fasta(update, gkn_option("-fasta"));
	if (gkn_option("-gff"))   test_gff(update,   gkn_option("-gff"));
	if (gkn_option("-pwm"))   test_pwm(update,   gkn_option("-pwm"));
	if (gkn_option("-mm"))    test_mm(update,    gkn_option("-mm"));
	if (gkn_option("-len"))   test_len(update,   gkn_option("-len"));
	if (gkn_option("-this"))  test_this();

	return 0;
}

void test_math() {
	printf("math ");
	printf("(%.2f %.2f %.2f) ", gkn_p2s(0), gkn_p2s(0.25), gkn_p2s(1));
	printf("(%.2f)", pow(2.71828, gkn_sum2(log(0.375), log(0.125))));
	printf(" done\n");
}

void test_seq() {
	printf("seq ");
	printf(" %d", gkn_ntindex("ACGT", 0, 3));
	nt_index is backwards and should be gkn_idx2seq and gkn_seq2idx
	printf(" done\n");
}

void test_vec(int update) {
	printf("vec ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}
		gkn_vec vec = gkn_vec_new();
		for (int j = 0; j < 99999; j++) {
			gkn_vec_push(vec, NULL);
		}
		for (int j = 0; j < 99999; j++) {
			if (vec->elem[j] != NULL) gkn_exit("ivec elem failure");
		}
		for (int j = 99998; j >= 0; j--) {
			void *v = gkn_vec_pop(vec);
			if (v != NULL) gkn_exit("ivec pop failure");
		}
		gkn_vec_free(vec);

	}
	printf(" done \n");
}

void test_ivec(int update) {
	printf("ivec ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}
		gkn_ivec vec = gkn_ivec_new();
		for (int j = 0; j < 99999; j++) {
			gkn_ivec_push(vec, j);
		}
		for (int j = 0; j < 99999; j++) {
			if (vec->elem[j] != j) gkn_exit("ivec elem failure");
		}
		for (int j = 99998; j >= 0; j--) {
			int v = gkn_ivec_pop(vec);
			if (v != j) gkn_exit("ivec pop failure");
		}
		gkn_ivec_free(vec);

	}
	printf(" done \n");
}

void test_fvec(int update) {
	printf("fvec ");

	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}
		gkn_fvec vec = gkn_fvec_new();
		for (int j = 0; j < 99999; j++) {
			gkn_fvec_push(vec, j);
		}
		for (int j = 0; j < 99999; j++) {
			if (vec->elem[j] != j) gkn_exit("ivec elem failure");
		}
		for (int j = 99998; j >= 0; j--) {
			int v = gkn_fvec_pop(vec);
			if (v != j) gkn_exit("ivec pop failure");
		}
		gkn_fvec_free(vec);

	}
	printf(" done\n");
}

void test_tvec(int update) {
	char text[16];
	sprintf(text, "hello world");

	printf("tvec ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}
		gkn_tvec tvec = gkn_tvec_new();
		for (int j = 0; j < 9999; j++) {
			gkn_tvec_push(tvec, text);
		}
		for (int j = 0; j < 9999; j++) {
			if (strcmp(tvec->elem[j], text) != 0) gkn_exit("tvec elem failure");
		}
		for (int j = 9998; j >= 0; j--) {
			char  *v = gkn_tvec_pop(tvec);
			if (strcmp(v, text) != 0) gkn_exit("tvec pop failure");
			free(v);
		}
		gkn_tvec_free(tvec);
	}
	printf(" done\n");
}

void test_map(int update) {
	char text[32];

	printf("map ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}
		gkn_map map = gkn_map_new();
		for (int j = 0; j < 999; j++) {
			sprintf(text, "key %d", j);
			gkn_map_set(map, text, text);
		}
		gkn_map_free(map);
	}
	printf(" done\n");
}

void test_tmap(int update) {
	char text[32];

	printf("tmap ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}
		gkn_tmap tmap = gkn_tmap_new();
		for (int j = 0; j < 555; j++) {
			sprintf(text, "key %d", j);
			gkn_tmap_set(tmap, text, text);
		}

		gkn_tvec keys = gkn_tmap_keys(tmap);
		for (int k = 0; k < keys->size; k++) {
			char* v = gkn_tmap_get(tmap, keys->elem[k]);
			if (strcmp(v, keys->elem[k]) != 0) {
				printf("key %s has value %s\n", keys->elem[k], v);
				gkn_exit("tmap integrity failure");
			}
		}
		gkn_tvec_free(keys);
		gkn_tmap_free(tmap);
	}
	printf(" done\n");
}

void test_xtree(int update) {
	char text[32];

	printf("xtree ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}
		gkn_xtree xt = gkn_xtree_new();

		for (int j = 0; j < 1111; j++) {
			sprintf(text, "key %d", j);
			gkn_xtree_set(xt, text, NULL);
		}

		gkn_tvec keys = gkn_xtree_keys(xt);
		for (int k = 0; k < keys->size; k++) {
			void *v = gkn_xtree_get(xt, keys->elem[k]);
			if (v != NULL) {
				printf("key %s has non NULL value\n", keys->elem[k]);
				gkn_exit("xtree integrity failure");
			}
		}
		gkn_tvec_free(keys);
		gkn_xtree_free(xt);

	}
	printf(" done\n");
}

void test_feat(int update) {
	char *seq = "NNNNNAAAAAGTAAGTTTTTTTTCAGAAAAANNNNN";
	printf("feat ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}
		for (int j = 0; j < 5555; j++) {
			gkn_feat f = gkn_feat_new(seq, 10, 25);
			char *s = gkn_feat_seq(f);
			if (strcmp(s, "GTAAGTTTTTTTTCAG") != 0) gkn_exit("feat bad");
			free(s);
			gkn_feat_free(f);
		}
	}
	printf(" done\n");
}

void test_pipe(int update, const char *filename) {
	printf("pipe ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}

		for (int j = 0; j < 100; j ++) {
			gkn_pipe io = gkn_pipe_open(filename, "r");
			gkn_pipe_close(io);
		}
	}
	printf(" done\n");
}

void test_read(int update, const char *filename) {
	printf("read ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}

		for (int j = 0; j < 100; j ++) {
			gkn_pipe io = gkn_pipe_open(filename, "r");
			char *line = gkn_readline(io);
			free(line);
			gkn_pipe_close(io);
		}
	}
	printf(" done\n");
}

void test_fasta(int update, const char *filename) {
	gkn_fasta in;
	printf("fasta ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}

		for (int j = 0; j < 20; j ++) {
			gkn_pipe io = gkn_pipe_open(filename, "r");
			while ((in = gkn_fasta_read(io)) != NULL) {
				gkn_fasta_free(in);
			}
			gkn_pipe_close(io);
		}
	}
	printf(" done\n");
}

void test_gff(int update, const char *filename) {
	gkn_gff gff;
	printf("gff ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}

		for (int j = 0; j < 10; j ++) {
			gkn_pipe io = gkn_pipe_open(filename, "r");
			while ((gff = gkn_gff_read(io)) != NULL) {
				gkn_gff_free(gff);
			}
			gkn_pipe_close(io);
		}
	}
	printf(" done\n");
}

void test_pwm(int update, const char *filename) {
	printf("pwm ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}
		for (int j = 0; j < 50; j++) {
			gkn_pipe io = gkn_pipe_open(filename, "r");
			gkn_pwm pwm = gkn_pwm_read(io);
			gkn_pwm_free(pwm);
			gkn_pipe_close(io);
		}
	}
	printf(" done\n");
}

void test_mm(int update, const char *filename) {
	printf("mm ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}
		for (int j = 0; j < 10; j++) {
			gkn_pipe io = gkn_pipe_open(filename, "r");
			gkn_mm mm = gkn_mm_read(io);
			gkn_mm_free(mm);
			gkn_pipe_close(io);
		}
	}
	printf(" done\n");
}

void test_len(int update, const char *filename) {
	printf("len ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}
		for (int j = 0; j < 10; j++) {
			gkn_pipe io = gkn_pipe_open(filename, "r");
			gkn_len len = gkn_len_read(io);
			gkn_len_free(len);
			gkn_pipe_close(io);
		}
	}
	printf(" done\n");
}

void test_smat(int update) {
	printf("smat ");
	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}
		for (int j = 0; j < 999; j++) {
			gkn_smat blosum = gkn_smat_blosum(62);
			gkn_smat_free(blosum);
			gkn_smat nt = gkn_smat_mng(1, -1, -2);
			gkn_smat_free(nt);
		}
	}
	printf(" done\n");
}

void test_sw(int update) {
	printf("sw ");
	gkn_smat b62 = gkn_smat_blosum(62);
	char *s1 = "AAAAACDEF";
	char *s2 = "ACDEFFFFF";

	for (int i = 0; i < COUNT; i++) {
		if (i % update == 0) {
			printf(".");
			fflush(stdout);
		}
		// stuff here
		for (int j = 0; j < 200; j++) {
			gkn_hsp msp = gkn_sw(s1, s2, b62);
			gkn_hsp_free(msp);
		}
	}

	gkn_smat_free(b62);
	printf(" done\n");
}

void test_this(void) {

	// current private test
	printf("hello world\n");

}
