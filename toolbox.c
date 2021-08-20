/******************************************************************************\
 toolbox.c
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef GENOMIKON_TOOLBOX_C
#define GENOMIKON_TOOLBOX_C

#include "toolbox.h"

static char gkn_version_number[] = "genomikon-2021";
static char gkn_program_name[256] = "name not set";

char * gkn_get_version_number (void) {return gkn_version_number;}
void   gkn_set_program_name (const char *s) {strcpy(gkn_program_name, s);}
char * gkn_get_program_name (void) {return gkn_program_name;}

void * gkn_malloc(size_t size) {
	void *mem = malloc(size);
	if (mem == NULL) gkn_exit("gkn_malloc %d", size);
	return mem;
}

void * gkn_calloc(size_t count, size_t size) {
	void *mem = calloc(count, size);
	if (mem == NULL) gkn_exit("gkn_calloc %d %d", count, size);
	return mem;
}

void * gkn_realloc(void *p, size_t size) {
	void *mem = realloc(p, size);
	if (mem == NULL) gkn_exit("gkn_realloc %d", size);
	return mem;
}


void gkn_ivec_free(gkn_ivec vec) {
	if (vec == NULL) return;
	if (vec->elem) free(vec->elem);
	free(vec);
}

gkn_ivec gkn_ivec_new(void) {
	gkn_ivec vec = gkn_malloc(sizeof(struct gkn_IVEC));
	vec->size  = 0;
	vec->limit = 0;
	vec->elem  = NULL;
	return vec;
}

void gkn_ivec_push(gkn_ivec vec, int val) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit	 = 1;
		else				 vec->limit *= 2;
		vec->elem = gkn_realloc(vec->elem, vec->limit * sizeof(int));
	}
	vec->elem[vec->size] = val;
	vec->last = val;
	vec->size++;
}

int gkn_ivec_pop(gkn_ivec vec) {
	if (vec->size == 0) gkn_exit("can't pop a zero-length vector");
	vec->size--;
	return vec->elem[vec->size];
}

void gkn_fvec_free(gkn_fvec vec) {
	if (vec == NULL) return;
	if (vec->elem) free(vec->elem);
	free(vec);
}

gkn_fvec gkn_fvec_new(void) {
	gkn_fvec vec = gkn_malloc(sizeof(struct gkn_FVEC));
	vec->size  = 0;
	vec->limit = 0;
	vec->elem  = NULL;
	return vec;
}

void gkn_fvec_push(gkn_fvec vec, float val) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit	 = 1;
		else				 vec->limit *= 2;
		vec->elem = gkn_realloc(vec->elem, vec->limit * sizeof(float));
	}
	vec->elem[vec->size] = val;
	vec->last = val;
	vec->size++;
}

float gkn_fvec_pop(gkn_fvec vec) {
	if (vec->size == 0) gkn_exit("can't pop a zero-length vector");
	vec->size--;
	return vec->elem[vec->size];
}

void gkn_tvec_free(gkn_tvec vec) {
	if (vec == NULL) return;
	if (vec->elem) {
		for (int i = 0; i < vec->size; i++) free(vec->elem[i]);
		free(vec->elem);
	}
	free(vec);
}

gkn_tvec gkn_tvec_new(void) {
	gkn_tvec vec = gkn_malloc(sizeof(struct gkn_TVEC));
	vec->size  = 0;
	vec->limit = 0;
	vec->elem  = NULL;
	return vec;
}

void gkn_tvec_push(gkn_tvec vec, const char *text) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit	 = 1;
		else				 vec->limit *= 2;
		vec->elem = gkn_realloc(vec->elem, vec->limit * sizeof(char *));
	}
	vec->elem[vec->size] = gkn_malloc(strlen(text) + 1);
	strcpy(vec->elem[vec->size], text);
	vec->last = vec->elem[vec->size];
	vec->size++;
}

char * gkn_tvec_pop(gkn_tvec vec) {
	if (vec->size == 0) gkn_exit("can't pop a zero-length vector");
	vec->size--;
	return vec->elem[vec->size];
}

void gkn_vec_free(gkn_vec vec) {
	if (vec == NULL) return;
	if (vec->elem) free(vec->elem);
	free(vec);
}

gkn_vec gkn_vec_new(void) {
	gkn_vec vec = gkn_malloc(sizeof(struct gkn_VEC));
	vec->size  = 0;
	vec->limit = 0;
	vec->elem  = NULL;
	return vec;
}

void gkn_vec_push(gkn_vec vec, void *p) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit	 = 1;
		else				 vec->limit *= 2;
		vec->elem = gkn_realloc(vec->elem, vec->limit * sizeof(void *));
	}
	vec->elem[vec->size] = p;
	vec->last = vec->elem[vec->size];
	vec->size++;
}

void * gkn_vec_pop(gkn_vec vec) {
	if (vec->size == 0) gkn_exit("can't pop a zero-length vector");
	vec->size--;
	return vec->elem[vec->size];
}

// hashing materials
static double HASH_MULTIPLIER[7] = {
	3.1415926536, // PI
	2.7182818285, // e
	1.6180339887, // golden mean
	1.7320508076, // square root of 3
	2.2360679775, // square root of 5
	2.6457513111, // square root of 7
	3.3166247904  // square root of 11
};
static float MAX_HASH_DEPTH = 2.0;
static int HashLevelToSlots(int level) {return pow(4, level);}
static int HashFunc(const gkn_map hash, const char *key) {
	double sum = 0;
	for (int i = 0; i < strlen(key); i++)
		sum += key[i] * HASH_MULTIPLIER[i % 7];
	return (int) (hash->slots * (sum - floor(sum)));
}

static void ExpandHash(gkn_map hash) {
	int		 oldslots = hash->slots;
	gkn_vec *oldkey = hash->key;
	gkn_vec *oldval = hash->val;
	gkn_vec	 kvec;
	gkn_vec	 vvec;
	gkn_tvec keys;

	// create the new hash
	hash->level = hash->level +1;
	hash->slots = HashLevelToSlots(hash->level);
	hash->key	= gkn_malloc(hash->slots * sizeof(struct gkn_VEC));
	hash->val	= gkn_malloc(hash->slots * sizeof(struct gkn_VEC));
	for (int i = 0; i < hash->slots; i++) {
		hash->key[i] = gkn_vec_new();
		hash->val[i] = gkn_vec_new();
	}

	// brand new hash?
	if (hash->keys->size == 0) return;

	keys = hash->keys;
	hash->keys = gkn_tvec_new();

	// transfer old stuff to new hash
	for (int i = 0; i < oldslots; i++) {
		kvec = oldkey[i];
		vvec = oldval[i];
		for (int j = 0; j < kvec->size; j++) {
			char *key = kvec->elem[j];
			char *val = vvec->elem[j];
			gkn_map_set(hash, key, val);
		}
	}

	// free old stuff
	for (int i = 0; i < oldslots; i++) {
		kvec = oldkey[i];
		vvec = oldval[i];
		gkn_vec_free(kvec);
		gkn_vec_free(vvec);
	}
	free(oldkey);
	free(oldval);
	gkn_tvec_free(keys);
}

void gkn_map_free(gkn_map hash) {
	if (hash == NULL) return;
	for (int i = 0; i < hash->slots; i++) {
		if (hash->key[i]) gkn_vec_free(hash->key[i]);
		if (hash->val[i]) gkn_vec_free(hash->val[i]);
	}
	gkn_tvec_free(hash->keys);
	gkn_vec_free(hash->vals);
	free(hash->key);
	free(hash->val);
	free(hash);
}

gkn_map gkn_map_new(void) {
	gkn_map hash = gkn_malloc(sizeof(struct gkn_MAP));
	hash->level = 0;
	hash->slots = 0;
	hash->keys	= gkn_tvec_new();
	hash->vals	= gkn_vec_new();
	hash->key	= NULL;
	hash->val	= NULL;
	ExpandHash(hash);
	return hash;
}

void * gkn_map_get(const gkn_map hash, const char *key) {
	int index = HashFunc(hash, key);

	// resolve collisions
	for (int i = 0; i < hash->key[index]->size; i++) {
		char *string = hash->key[index]->elem[i];
		if (strcmp(key, string) == 0) {
			return hash->val[index]->elem[i];
		}
	}
	return NULL; // return is NULL if not found
}

void gkn_map_set(gkn_map hash, const char *key, void *val) {
	int	new_key = 1;
	int index = HashFunc(hash, key);

	// reassign unless new key
	for (int i = 0; i < hash->key[index]->size; i++) {
		char *string = hash->key[index]->elem[i];
		if (strcmp(key, string) == 0) {
			hash->val[index]->elem[i] = val;
			new_key = 0;
			return;
		}
	}

	if (new_key) {
		gkn_tvec_push(hash->keys, key);
		gkn_vec_push(hash->key[index], hash->keys->last);
		gkn_vec_push(hash->vals, val);
		gkn_vec_push(hash->val[index], hash->vals->last);
	}

	// check if we have to expand the hash
	if ((float)hash->keys->size / (float)hash->slots >= MAX_HASH_DEPTH) {
		ExpandHash(hash);
	}
}

gkn_tvec gkn_map_keys(const gkn_map hash) {
	gkn_tvec vec = gkn_tvec_new();
	for (int i = 0; i < hash->keys->size; i++) gkn_tvec_push(vec, hash->keys->elem[i]);
	return vec;
}

gkn_vec gkn_map_vals(const gkn_map hash) {
	gkn_vec vec = gkn_vec_new();
	for (int i = 0; i < hash->vals->size; i++) gkn_vec_push(vec, hash->vals->elem[i]);
	return vec;
}

void gkn_map_stat(const gkn_map hash) {
	int max = 0;
	int min = INT_MAX;
	int total = 0;
	for (int i = 0; i < hash->slots; i++) {
		int count = hash->val[i]->size;
		total += count;
		if (count > max) max = count;
		if (count < min) min = count;
	}
	fprintf(stdout, "HashStats: level=%d slots=%d keys=%d min=%d max=%d ave=%f\n",
		 hash->level, hash->slots, hash->keys->size, min, max,
		 (float)total / (float)hash->slots);
}

// text map

void gkn_tmap_free(gkn_tmap t) {
	gkn_map_free(t->hash);
	gkn_tvec_free(t->tvec);
	free(t);
}

gkn_tmap gkn_tmap_new(void) {
	gkn_tmap t = gkn_malloc(sizeof(struct gkn_TMAP));
	t->hash = gkn_map_new();
	t->tvec = gkn_tvec_new();
	return t;
}

void gkn_tmap_set(gkn_tmap t, const char *key, const char *val) {
	gkn_tvec_push(t->tvec, val);
	gkn_map_set(t->hash, key, t->tvec->last);
}

int gkn_tmap_exists(const gkn_tmap t, const char *key) {
	void *ref = gkn_map_get(t->hash, key);
	if (ref == NULL) return 0;
	return 1;
}

char * gkn_tmap_get(const gkn_tmap t, const char *key) {
	void *ref = gkn_map_get(t->hash, key);
	assert(ref != NULL);
	return ref;
}

gkn_tvec gkn_tmap_keys(const gkn_tmap t) {
	return gkn_map_keys(t->hash);
}

// command line options

static gkn_tvec COMMAND_LINE = NULL;
static gkn_map CL_REGISTER	= NULL;
static gkn_map CL_OPTIONS	= NULL;

void gkn_register_option(const char *name, int flag) {
	if (COMMAND_LINE == NULL) {
		COMMAND_LINE = gkn_tvec_new();
		CL_REGISTER	 = gkn_map_new();
		CL_OPTIONS	 = gkn_map_new();
	}

	switch (flag) {
		case 0: gkn_map_set(CL_REGISTER, name, (void *)1); break;
		case 1: gkn_map_set(CL_REGISTER, name, (void *)2); break;
		default: gkn_exit("gkn_register_option: flag 0 or 1");
	}
}

void gkn_parse_options(int *argc, char **argv) {
	for (int i = 0; i < *argc; i++) {
		char *token = argv[i];
		if (token[0] == '-' && strlen(token) > 1) {
			switch ((size_t)gkn_map_get(CL_REGISTER, token)) {
				case 0:
					gkn_exit("unknown option (%s)", token);
					break;
				case 1:
					gkn_map_set(CL_OPTIONS, token, token);
					break;
				case 2:
					gkn_map_set(CL_OPTIONS, token, argv[i+1]);
					i++;
					break;
				default:
					gkn_exit("not possible");
			}
		} else {
			gkn_tvec_push(COMMAND_LINE, argv[i]);
		}
	}

	*argc = COMMAND_LINE->size;
	for (int i = 0; i < COMMAND_LINE->size; i++) {
		argv[i] = COMMAND_LINE->elem[i];
	}
}

char * gkn_option(const char *tag) {
	return gkn_map_get(CL_OPTIONS, tag);
}

// pipe

void gkn_pipe_close(gkn_pipe pipe) {
	pipe->mode = 0;
	free(pipe->name);
	if (pipe->gzip) pclose(pipe->stream);
	else			fclose(pipe->stream);
	pipe->gzip = 0;
	free(pipe);
}

gkn_pipe gkn_pipe_open(const char *name, const char *mode) {
	char	command[1024];
	int     length = strlen(name);
	gkn_pipe pipe = gkn_malloc(sizeof(struct gkn_PIPE));

	if		(strcmp(mode, "r") == 0)  pipe->mode = 0;
	else if (strcmp(mode, "w") == 0)  pipe->mode = 1;
	else if (strcmp(mode, "r+") == 0) pipe->mode = 2;
	else gkn_exit("r, w, or r+ only in gkn_pipe");

	pipe->name = gkn_malloc(length + 1);
	strcpy(pipe->name, name);

	pipe->gzip = 0;

	if (name[length -3] == '.' &&
		name[length -2] == 'g' &&
		name[length -1] == 'z') pipe->gzip = 1; // .gz
	if (name[length -2] == '.' &&
		name[length -1] == 'z') pipe->gzip = 1; // .z
	if (name[length -2] == '.' &&
		name[length -1] == 'Z') pipe->gzip = 1; // .Z

	if (pipe->gzip) {
		if (pipe->mode != 0) gkn_exit("compressed pipes are read only");
		sprintf(command, "gunzip -c %s", name);
		pipe->stream = popen(command, "r");
	} else {
		if (strcmp(name, "-") == 0) pipe->stream = stdin;
		else						pipe->stream = fopen(name, mode);
	}

	if (pipe->stream == NULL) {
		gkn_exit("failed to open %s\n", name);
	}

	return pipe;
}

char * gkn_readline(gkn_pipe io) {
	char line[4096];
	int  read = 0;
	while (fgets(line, sizeof(line), io->stream) != NULL) {
		if (line[0] == '#') continue;
		if (strlen(line) == 0) continue;
		read = 1;
		break;
	}
	if (read == 0) return NULL;
	char *out = malloc(strlen(line) + 1);
	strcpy(out, line);
	return out;
}

void gkn_exit(const char* format, ...) {
	va_list args;
	fflush(stdout);
	fprintf(stderr, "ERROR from program %s, libarary %s\n",
		gkn_get_program_name(),
		gkn_get_version_number());
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
	fprintf(stderr, "\n");
	exit(1);
}

#endif
