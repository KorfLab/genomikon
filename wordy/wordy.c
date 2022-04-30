/*****************************************************************************\
 wordy.c
\*****************************************************************************/

#include <time.h>
#include "genomikon.h"

#define WORD_LENGTH 64
#define MAX_BOARD 8

struct DICTIONARY {
	gkn_xtree tree;
	int       count[256];
	float     prob[256];
	float     ramp[256];
};
typedef struct DICTIONARY dictionary;

dictionary read_dictionary (const char *filename) {
	char word[1024];

	// blank dictionary
	dictionary d;
	d.tree = gkn_xtree_new();
	for (int i = 0; i < 256; i++) {
		d.count[i] = 0;
		d.prob[i] = 0;
	}

	// read words
	int total = 0;
	gkn_pipe pipe = gkn_pipe_open(filename, "r");
	while (fscanf(pipe->stream, "%s", word) != EOF) {
		int length = strlen(word);
		if (length >= WORD_LENGTH) gkn_exit("max word length exceeded");
		for (int i = 0; i < length; i++) word[i] = tolower(word[i]);
		gkn_xtree_set(d.tree, word, (void*)1);
		for (int i = 0; i < length; i++) {
			d.count[(int)word[i]]++;
			total++;
		}
	}
	gkn_pipe_close(pipe);

	// symbol frequencies
	for (int i = 0; i < 256; i++) {
		d.prob[i] = (float)d.count[i] / (float)total;
	}

	double sum = 0;
	for (int i = 0; i < 256; i++) {
		d.ramp[i] = d.prob[i] + sum;
		sum += d.prob[i];
	}

	return d;
}

char random_char(dictionary d) {
	double r = random() / (double)RAND_MAX;
	for (int i = 0; i < 256; i++) {
		if (r < d.ramp[i]) {
			return i;
		}
	}
	return -1;
}

struct BOARD {
	int size;
	int matrix[MAX_BOARD][MAX_BOARD];
};
typedef struct BOARD * board;

void free_board(board b) {
	free(b);
}

board new_board(int s) {
	assert(s <= MAX_BOARD);
	board b = malloc(sizeof(struct BOARD));
	b->size = s;
	for (int i = 0; i < b->size; i++) {
		for (int j = 0; j < b->size; j++) {
			b->matrix[i][j] = 0;
		}
	}
	return b;
}

board random_board(int s, dictionary d) {
	board b = new_board(s);
	for (int i = 0; i < b->size; i++) {
		for (int j = 0; j < b->size; j++) {
			char c = random_char(d);
			b->matrix[i][j] = c;
		}
	}
	return b;
}

board copy_board (const board parent) {
	board child = new_board(parent->size);
	for (int i = 0; i < parent->size; i++) {
		for (int j = 0; j < parent->size; j++) {
			child->matrix[i][j] = parent->matrix[i][j];
		}
	}
	return child;
}

void show_board (const board b) {
	for (int x = 0; x < b->size +2; x++) printf("-");
	printf("\n");

	for (int y = 0; y < b->size; y++) {
		printf("|");
		for (int x = 0; x < b->size; x++) {
			if (b->matrix[x][y] < 'a') printf("%d", b->matrix[x][y]);
			else                       printf("%c", b->matrix[x][y]);
		}
		printf("|\n");
	}

	for (int x = 0; x < b->size +2; x++) printf("-");
	printf("\n");
}

static const int MOVE[8][2] = {
	{0,1}, {1,1}, {1,0}, {1,-1}, {0,-1}, {-1,-1}, {-1,0}, {-1,1}
};

static void find_words (const board b,
	const dictionary d,
	const int x0,
	const int y0,
	const board path,
	const char *word,
	int length,
	gkn_tvec words)
{
	char new_word[WORD_LENGTH];
	if (gkn_xtree_get(d.tree, word)) {
		gkn_tvec_push(words, word);
	}

	for (int i = 0; i < 8; i++) {
		int x = MOVE[i][0] + x0;
		int y = MOVE[i][1] + y0;

		if (x < 0 || x == b->size) continue; // out of bounds
		if (y < 0 || y == b->size) continue; // out of bounds
		if (path->matrix[x][y])    continue; // already visited

		strcpy(new_word, word);
		new_word[length] = b->matrix[x][y];
		new_word[length+1] = '\0';
		if (!gkn_xtree_check(d.tree, new_word)) continue; // no words extend

		board new_path = copy_board(path);
		new_path->matrix[x][y] = 1;
		find_words(b, d, x, y, new_path, new_word, length +1, words);
		free_board(new_path);
	}
}

gkn_tvec solve_board (const board b, const dictionary d, int min_word) {
	gkn_tvec words = gkn_tvec_new();
	for (int x = 0; x < b->size; x++) {
		for (int y = 0; y < b->size; y++) {
			board path = new_board(b->size);
			path->matrix[x][y] = 1;
			char word[2];
			word[0] = b->matrix[x][y];
			word[1] = '\0';
			find_words(b, d, x, y, path, word, 1, words);
			free_board(path);
		}
	}

	// filter for length and redundancy
	gkn_map unique = gkn_map_new();
	for (int i = 0; i < words->size; i++) {
		if (strlen(words->elem[i]) >= min_word) {
			gkn_map_set(unique, words->elem[i], (void*)1);
		}
	}
	gkn_tvec final = gkn_map_keys(unique);

	// clean up
	gkn_tvec_free(words);
	gkn_map_free(unique);

	return final;
}

struct BEING {
	board genotype;
	int   fitness;
};
typedef struct BEING being;

static int fitness(board b, dictionary d, int min) {
	gkn_tvec words = solve_board(b, d, min);
	int fit = words->size;
	gkn_tvec_free(words);
	return fit;
}

static int wsort (const void *p1, const void *p2) {
	char *a = *(char **) p1;
	char *b = *(char **) p2;
	return strcmp(a, b);
}

static int bsort (const void *p1, const void *p2) {
	being *a = (being*) p1;
	being *b = (being*) p2;
	if (a->fitness > b->fitness) return -1;
	if (b->fitness > a->fitness) return 1;
	return 0;
}

static void mate(
	const being p1,
	const being p2,
	being *child,
	dictionary d,
	int min,
	float mut)
{

	// recombination & mutation
	for (int i = 0; i < p1.genotype->size; i++) {
		for (int j = 0; j < p1.genotype->size; j++) {
			if (rand()/RAND_MAX < 0.5)
				child->genotype->matrix[i][j] = p1.genotype->matrix[i][j];
			else
				child->genotype->matrix[i][j] = p2.genotype->matrix[i][j];

			if (rand() / RAND_MAX < mut)
				child->genotype->matrix[i][j] = random_char(d);
		}
	}

	// finalization
	child->fitness = fitness(child->genotype, d, min);
}

static char *usage = "\
dusty - boggle solver demo\n\n\
usage: wordy <dictionary> [options]\n\
options:\n\
  -size <int>   size of board [4]\n\
  -min <int>    minimum word size [3]\n\
  -pop <int>    population size [1000]\n\
  -gen <int>    generations [100]\n\
  -die <float>  fraction that die each generation [0.5]\n\
  -mut <float>  mutation frequency [0.1]\n\
  -seed <int>   random seed [not set]\n\
  -words        show all words\n\
  -verbose      verbose progress reporting\n\
";

int main(int argc, char **argv) {
	char *file = NULL;       // path to dictionary file
	int   size = 4;          // board size
	int   min = 3;           // window size
	int   pop = 1000;        // population size
	int   gen = 100;         // generations
	float die = 0.5;         // fraction that die each gen
	float mut = 0.1;         // mutation frequency
	int   seed = 0;          // random seed, 0 not used
	int   show_words = 0;    // final report
	int   show_progress = 0; // intermediate report

	// Command Line Interface
	gkn_set_program_name(argv[0]);
	gkn_register_option("-size", 1);
	gkn_register_option("-min", 1);
	gkn_register_option("-pop", 1);
	gkn_register_option("-gen", 1);
	gkn_register_option("-die", 1);
	gkn_register_option("-mut", 1);
	gkn_register_option("-seed", 1);
	gkn_register_option("-words", 0);
	gkn_register_option("-verbose", 0);
	gkn_parse_options(&argc, argv);

	if (argc == 1) gkn_exit("%s", usage);

	file = argv[1];
	if (gkn_option("-size"))    size = atoi(gkn_option("-size"));
	if (gkn_option("-min"))     min = atoi(gkn_option("-min"));
	if (gkn_option("-pop"))     pop = atoi(gkn_option("-pop"));
	if (gkn_option("-gen"))     gen = atoi(gkn_option("-gen"));
	if (gkn_option("-die"))     die = atof(gkn_option("-die"));
	if (gkn_option("-mut"))     mut = atof(gkn_option("-mut"));
	if (gkn_option("-seed"))    seed = atoi(gkn_option("-seed"));
	if (gkn_option("-words"))   show_words = 1;
	if (gkn_option("-verbose")) show_progress = 1;

	// setup
	dictionary d = read_dictionary(file);
	if (seed) srand(seed);
	else srand(time(NULL));
	int half = pop * (1 - die);

	// initial population
	being *population = malloc(pop * sizeof(being));
	for (int i = 0; i < pop; i++) {
		population[i].genotype = random_board(size, d);
		population[i].fitness = fitness(population[i].genotype, d, min);
	}

	// evolve
	for (int g = 0; g < gen; g++) {
		qsort(population, pop, sizeof(being), bsort);
		if (show_progress)
			printf("gen: %d, max: %d\n", g, population[0].fitness);
		for (int i = half; i < pop; i++) {
			int p1 = rand() / (double)RAND_MAX * half;
			int p2 = rand() / (double)RAND_MAX * half;
			mate(population[p1], population[p2], &population[i], d, min, mut);
		}
	}

	// report
	qsort(population, pop, sizeof(being), bsort);
	board b = population[0].genotype;
	gkn_tvec words = solve_board(b, d, min);
	show_board(b);
	printf("words: %d\n", words->size);
	if (show_words) {
		qsort(words->elem, words->size, sizeof(char*), wsort);
		for (int i = 0; i < words->size; i++) {
			printf("%s\n", words->elem[i]);
		}
	}

	return 0;
}
