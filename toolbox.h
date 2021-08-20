/******************************************************************************\
 toolbox.h
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef GENOMIKON_TOOLBOX_H
#define GENOMIKON_TOOLBOX_H

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// library and program info
char * gkn_get_version_number(void);
void   gkn_set_program_name(const char *);
char * gkn_get_program_name(void);

// memory
void * gkn_malloc(size_t);
void * gkn_calloc(size_t, size_t);
void * gkn_realloc(void *, size_t);

// integer vector
struct gkn_IVEC {
	int * elem;
	int   size;
	int   limit;
	int   last;
};
typedef struct gkn_IVEC * gkn_ivec;
void	 gkn_ivec_free(gkn_ivec);
gkn_ivec gkn_ivec_new(void);
void	 gkn_ivec_push(gkn_ivec, int);
int      gkn_ivec_pop(gkn_ivec);

// float vector
struct gkn_FVEC {
	float * elem;
	int     size;
	int     limit;
	int     last;
};
typedef struct gkn_FVEC * gkn_fvec;
void	 gkn_fvec_free(gkn_fvec);
gkn_fvec gkn_fvec_new(void);
void	 gkn_fvec_push(gkn_fvec, float);
float    gkn_fvec_pop(gkn_fvec);

// text vector
struct gkn_TVEC {
	char ** elem;
	int     size;
	int     limit;
	char  * last;
};
typedef struct gkn_TVEC * gkn_tvec;
void	 gkn_tvec_free(gkn_tvec);
gkn_tvec gkn_tvec_new(void);
void	 gkn_tvec_push(gkn_tvec, const char *);
char *   gkn_tvec_pop(gkn_tvec);

// generic void * vector
struct gkn_VEC {
	void ** elem;
	int     size;
	int     limit;
	void  * last;
};
typedef struct gkn_VEC * gkn_vec;
void    gkn_vec_free(gkn_vec);
gkn_vec gkn_vec_new(void);
void    gkn_vec_push(gkn_vec, void *);
void *  gkn_vec_pop(gkn_vec);

// generic map (text key, void * value)
struct gkn_MAP {
	int      level;
	int      slots;
	gkn_tvec keys;
	gkn_vec	 vals;
	gkn_vec * key;
	gkn_vec * val;
};
typedef struct gkn_MAP * gkn_map;
void	 gkn_map_free(gkn_map);
gkn_map	 gkn_map_new(void);
void	 gkn_map_set(gkn_map, const char *, void *);
void *	 gkn_map_get(const gkn_map, const char *);
gkn_tvec gkn_map_keys(const gkn_map);
gkn_vec	 gkn_map_vals(const gkn_map);
void	 gkn_map_stat(const gkn_map);

// text map
struct gkn_TMAP {
	gkn_map  hash;
	gkn_tvec tvec;
};
typedef struct gkn_TMAP * gkn_tmap;
void     gkn_tmap_free(gkn_tmap);
gkn_tmap gkn_tmap_new(void);
void     gkn_tmap_set(gkn_tmap, const char *, const char *);
char *   gkn_tmap_get(const gkn_tmap, const char *);
int      gkn_tmap_exists(const gkn_tmap, const char *);
gkn_tvec gkn_tmap_keys(const gkn_tmap);

// command line processing
void   gkn_register_option(const char *, int);
void   gkn_parse_options(int *, char **);
char * gkn_option(const char *);

// pipe
struct gkn_PIPE {
	int    mode; // 0 = read, 1 = write, 2 = r+
	char * name;
	int    gzip;
	FILE * stream;
};
typedef struct gkn_PIPE * gkn_pipe;
gkn_pipe gkn_pipe_open(const char *, const char *);
void     gkn_pipe_close(gkn_pipe);

// input/output
char * gkn_readline(gkn_pipe);
void   gkn_exit(const char *, ...);

#endif
