/******************************************************************************\
 hmm.c
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef GENOMIKON_HMM_C
#define GENOMIKON_HMM_C

#include "hmm.h"

void gkn_state_free (gkn_state s) {
	free(s->name);
	gkn_tvec_free(s->adjs);
	gkn_fvec_free(s->adjp);
	gkn_mm_free(s->emit);
	gkn_len_free(s->len);
}

/*

dnState dnState_read (FILE * stream) {
	int		i;
	char	tag[64], name[64];
	double	value;
	dnState s = ik_malloc(sizeof(struct dnState));
		
	if (fscanf(stream, "%s %s %lf %lf %d %d %d", tag, name, &s->init, &s->term,
		&s->transitions, &s->order, &s->durations) != 7)
		ik_exit (1, "denada-State parse error 1");
	s->init = p2s(s->init);
	s->term = p2s(s->term);
			
	// name 
	s->name = ik_malloc(strlen(name) +1);
	strcpy(s->name, name);
	
	// transitions 
	s->state = ik_tvec_new();
	s->score = ik_fvec_new();
	for (i = 0; i < s->transitions; i++) {
		if (fscanf(stream, "%s", name) != 1)
			ik_exit(2, "denada-State parse error 2");
		ik_tvec_push(s->state, name);
		if (fscanf(stream, "%lf", &value) != 1)
			ik_exit(3, "denada-State parse error 3");
		ik_fvec_push(s->score, p2s(value));
	}
	
	// emissions
	s->count = pow(4, (double)s->order +1);
	s->emit = ik_malloc(sizeof(double) * s->count);
	for (i = 0; i < s->count; i++) {
		if (fscanf(stream, "%lf", &value) != 1)
			ik_exit(4, "denada-State parse error 4");
		s->emit[i] = p2s(value);
	}
	
	// durations
	if (s->durations) s->duration = ik_malloc(sizeof(double) * s->durations);
	for (i = 0; i < s->durations; i++) {
		if (fscanf(stream, "%lf", &value) != 1) 
			ik_exit(5, "denada-State parse error 5");
		s->duration[i] = p2s(value);
	}
	
	// explicit length states may not start or end
	if (s->durations) {
		if (LOGARITHMS) {
			if (s->init != MIN_SCORE) ik_exit(1, "explicit states can't init");
			if (s->term != MIN_SCORE) ik_exit(1, "explicit states can't term");
		} else {
			if (s->init != 0) ik_exit(1, "explict states can't init");
			if (s->term != 0) ik_exit(1, "explict states can't term");
		}
	}
	
	return s;
}

void dnState_write (FILE * stream, const dnState s) {
	int i;
	
	// header
	fprintf(stream, "\tdenada-State %s %.3f %.3f %d %d\n\n", 
			s->name, s2p(s->init), s2p(s->term), s->transitions, s->durations);
	
	// transitions
	for (i = 0; i < s->transitions; i++) {
		fprintf(stream, "\t\t%s %.3f\n", s->state->elem[i], s2p(s->score->elem[i]));
	}
	
	// emissions
	for (i = 0; i < s->count; i++) {
		if (i % 4 == 0) fprintf(stream, "\n\t\t");
		fprintf(stream, "%.3f ", s2p(s->emit[i]));
	}
	if (s->durations) fprintf(stream, "\n");
	
	// durations
	for (i = 0; i < s->durations; i++) {
		if (i % 5 == 0) fprintf(stream, "\n\t\t");
		fprintf(stream, "%.3f ", s2p(s->duration[i]));
	}
	fprintf(stream, "\n\n");
}



void dnHMM_free (dnHMM hmm) {
	int i;
	ik_free(hmm->name);
	for (i = 0; i < hmm->states; i++) dnState_free(hmm->state[i]);
}

dnHMM dnHMM_read (FILE * stream) {
	int		 i, j, from, to, count;
	double	 score;
	char	 line[1024], tag[256], name[256];
	dnHMM	 hmm = ik_malloc(sizeof(struct dnHMM));
	ik_map	 name_to_number;
	
	// clear out
	hmm->name = NULL;
	hmm->states = 0;
	for (i = 0; i < dnSTATE_LIMIT; i++) hmm->state[i] = NULL;
	for (i = 0; i < dnSTATE_LIMIT; i++)
		for (j = 0; j < dnSTATE_LIMIT; j++) hmm->transition[i][j] = MIN_SCORE;
	
	// find denada-HMM header, skip comments (lines starting with #)
	while (fgets(line, sizeof(line), stream) != NULL) {
		if (line[0] == '#') continue;
		
		if (sscanf(line, "%s %s %d", tag, name, &hmm->states) == 3) {
			if (strcmp(tag, "denada-HMM") != 0) {
				ik_exit(1, "unrecognized HMM type: %s", line);
			} else {
				break;
			}
		}
	}
	
	// name
	hmm->name = ik_malloc(strlen(name) +1);
	strcpy(hmm->name, name);
	
	// states
	for (i = 0; i < hmm->states; i++) hmm->state[i] = dnState_read(stream);
	
	// state names and numbers
	name_to_number = ik_map_new();
	for (i = 0; i < hmm->states; i++) {
		ik_map_set(name_to_number, hmm->state[i]->name, (void*)i);
	}

	// transition
	for (i = 0; i < hmm->states; i++) {
		from = i;
		for (j = 0; j < hmm->state[i]->transitions; j++) {
			to = (size_t)ik_map_get(name_to_number, hmm->state[i]->state->elem[j]);
			score = hmm->state[i]->score->elem[j];
			hmm->transition[from][to] = score;
		}
	}


	for (i = 0; i < hmm->states; i++) {
		if (hmm->state[i]->durations == 0) continue;
		count = 0;
		for (j = 0; j < hmm->states; j++) {
			if (hmm->transition[i][j] != MIN_SCORE) count++;
		}
		if (count > 1) {
			ik_exit(1, "explicit duration states allowed only 1 link");
		}
	}
				
	for (i = 0; i < hmm->states; i++) {
		for (j = 0; j < hmm->states; j++) {
			if (i == j) continue;
			if (hmm->transition[i][j] == MIN_SCORE) continue;
			if (hmm->state[i]->durations && hmm->state[j]->durations)
				ik_exit(1, "explicit duration states must not be adjacent");
		}
	}
	
	
	// clean up
	ik_map_free(name_to_number);
	
	return hmm;
}

void dnHMM_write (FILE *stream, const dnHMM hmm) {
	int i;
	
	fprintf(stream, "denada-HMM %s %d\n\n", hmm->name, hmm->states);
	for (i = 0; i < hmm->states; i++) {
		dnState_write(stream, hmm->state[i]);
	}
}

static double emission (const dnState state, const char * seq, int pos) {
	int i, p, c, index;
	
	if (LOGARITHMS) {
		if (pos < state->order) return -2;
		if (seq[pos] == 4)		return -2;
	} else {
		if (pos < state->order) return 0.25;
		if (seq[pos] == 4)		return 0.25;
	}
	
	index = 0;
	for (i = 0; i <= state->order; i++) {
		c = seq[pos -i];
		if (c == 4) {
			if (LOGARITHMS) return -2;
			else			return 0.25;
		}
		p = ik_POWER[4][i];
		index += p * c;
	}
	
	return state->emit[index];
}

struct dnMax {
	double score;
	int	   state;
	int	   jumps;
};

static struct dnMax viterbi_calc (dnMatrix m, int this_coor, int this_state) {
	int		prev_state, prev_coor = this_coor -1;
	struct	dnMax max;
	double	emit, trans, prev, total_score, max_score, max_state;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
	
	max.score = MIN_SCORE;
	max.state = this_state;
	max.jumps = 0; // normal states do not jump
	
	emit = emission(hmm->state[this_state], dna->num, this_coor);
	if (emit == MIN_SCORE) return max;
	
	max_score = MIN_SCORE;
	max_state = -1;
	for (prev_state = 0; prev_state < hmm->states; prev_state++) {
		trans = hmm->transition[prev_state][this_state];
		prev = m->vscore[prev_state][prev_coor];
		
		if (trans == MIN_SCORE) continue;
		if (prev == MIN_SCORE) continue;
		
		total_score = emit + trans + prev;
		if (total_score > max_score) {
			max_score = total_score;
			max_state = prev_state;
		}
	}
	
	if (max_score == MIN_SCORE) return max;
	max.score = max_score;
	max.state = max_state;
	
	return max;
}


static struct dnMax viterbiX_calc (
								  dnMatrix m,
								  const int this_coor,
								  const int this_state)
{
	int		i, prev_state, limit;
	double	emit_score[dnDURATION_LIMIT], state_score[dnDURATION_LIMIT];
	double	prev_score, trans_score, this_score, total_score;
	struct	dnMax max;
	dnHMM	hmm = m->hmm;
	ik_dna	dna = m->dna;
	dnState state = m->hmm->state[this_state];
	
	max.score = MIN_SCORE;
	max.state = -1;
	max.jumps = 0;
	
	// pre-calculate emit_score for all distances as a cumulative score
	emit_score[0] = emission(state, dna->num, this_coor);	
	limit = state->durations;
	if (limit > this_coor) limit = this_coor;
	for (i = 1; i < limit; i++) {
		this_score = emission(state, dna->num, this_coor - i);
		prev_score = emit_score[i-1];
		if (this_score == MIN_SCORE) emit_score[i] = MIN_SCORE;
		else if (prev_score == MIN_SCORE) emit_score[i] = MIN_SCORE;
		else emit_score[i] = this_score + prev_score;
	}
	
	// pre-calculate state_score (emission + duration) 
	for (i = 0; i < limit; i++) {
		if (state->duration[i] == MIN_SCORE) state_score[i] = MIN_SCORE;
		else if (emit_score[i] == MIN_SCORE) state_score[i] = MIN_SCORE;
		else state_score[i] = emit_score[i] + state->duration[i];
	}
	
	// calculate max for each state and duration 
	for (prev_state = 0; prev_state < hmm->states; prev_state++) {
		if (prev_state == this_state) continue;
		trans_score = hmm->transition[prev_state][this_state];
		if (trans_score == MIN_SCORE) continue;
		for (i = 0; i < limit; i++) {
			if (state_score[i] == MIN_SCORE) continue;
			prev_score = m->vscore[prev_state][this_coor -i];
			if (prev_score == MIN_SCORE) continue;
			total_score = prev_score + state_score[i] + trans_score;
			if (total_score > max.score) {
				max.score = total_score;
				max.state = prev_state;
				max.jumps = i;
			}
		}
	}
	
	return max;
}


static struct dnMax viterbiP_calc (
								   dnMatrix m, 
								   int this_coor, 
								   int this_state) 
{
	int		prev_state, prev_coor = this_coor -1;
	struct	dnMax max;
	double	emit, total_score, max_score, max_state;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
		
	max.score = 0;
	max.state = this_state;
	max.jumps = 0; // normal states do not jump 
	
	emit = emission(hmm->state[this_state], dna->num, this_coor);
	
	max_score = 0;
	max_state = -1;
	for (prev_state = 0; prev_state < hmm->states; prev_state++) {
		total_score = emit * hmm->transition[prev_state][this_state]
			* m->vscore[prev_state][prev_coor];
		if (total_score > max_score) {
			max_score = total_score;
			max_state = prev_state;
		}
	}
	
	if (max_score == 0) return max;
	max.score = max_score;
	max.state = max_state;
	
	return max;
}

static struct dnMax viterbiXP_calc (
								   dnMatrix m, 
								   const int this_coor, 
								   const int this_state) 
{
	int		i, prev_state, limit;
	double	emit_score[dnDURATION_LIMIT], state_score[dnDURATION_LIMIT];
	double	prev_score, this_score, total_score;
	struct	dnMax max;
	dnHMM	hmm = m->hmm;
	ik_dna	dna = m->dna;
	dnState state = m->hmm->state[this_state];
		
	max.score = 0;
	max.state = -1;
	max.jumps = 0;
	
	// pre-calculate emit_score for all distances as a cumulative score
	emit_score[0] = emission(state, dna->num, this_coor);	
	limit = state->durations;
	if (limit > this_coor) limit = this_coor;
	for (i = 1; i < limit; i++) {
		this_score = emission(state, dna->num, this_coor - i);
		prev_score = emit_score[i-1];
		emit_score[i] = this_score * prev_score;
	}
	
	// pre-calculate state_score (emission * duration)
	for (i = 0; i < limit; i++) {
		state_score[i] = emit_score[i] * state->duration[i];
	}
	
	// calculate max for each state and duration
	for (prev_state = 0; prev_state < hmm->states; prev_state++) {
		for (i = 0; i < limit; i++) {
			prev_score = m->vscore[prev_state][this_coor -i];
			total_score = 
				prev_score 
				* state_score[i] 
				* hmm->transition[prev_state][this_state];
			if (total_score > max.score) {
				max.score = total_score;
				max.state = prev_state;
				max.jumps = i;
			}
		}
	}
	
	return max;
}

static void dnMatrix_viterbi (dnMatrix m) {
	int		i, j, k, max_state;
	double	max_score;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
	struct	dnMax max;
		
	// initialization
	i = 0;
	for (j = 0; j < hmm->states; j++) {
		if (hmm->state[j]->init == MIN_SCORE) {
			m->vscore[j][i] = MIN_SCORE;
		} else if (emission(hmm->state[j], dna->num, i) == MIN_SCORE) {
			m->vscore[j][i] = MIN_SCORE;
		} else {
			m->vscore[j][i] = hmm->state[j]->init 
				+ emission(hmm->state[j], dna->num, 0);
		}
		m->vtrace[j][i] = j;
		m->vjumps[j][i] = 0;
	}
	
	// induction
	for (i = 1; i < dna->length; i++) {
		for (j = 0; j < hmm->states; j++) {
			if (hmm->state[j]->durations) max = viterbiX_calc(m, i, j);
			else						  max = viterbi_calc(m, i, j);
			m->vscore[j][i] = max.score;
			m->vtrace[j][i] = max.state;
			m->vjumps[j][i] = max.jumps;
		}
	}
	
	// termination
	i = dna->length -1;
	max_state = -1;
	max_score = MIN_SCORE;
	for (j = 0; j < hmm->states; j++) {
		if (hmm->state[j]->term == MIN_SCORE) {
			m->vscore[j][i] = MIN_SCORE;
		} else {
			m->vscore[j][i] += hmm->state[j]->term;
		}
		if (m->vscore[j][i] > max_score) {
			max_score = m->vscore[j][i];
			max_state = j;
		}
	}
	m->vmax = max_score;
		
	// trace
	i = dna->length -1;
	j = max_state;			   
	ik_ivec_push(m->vpath, j); 
	while (i > 0) {
		if (m->vtrace[j][i] == -1) {
			ik_exit(1, "fatal traceback error");
		} else if (m->vtrace[j][i] == j) {
			ik_ivec_push(m->vpath, j);
			i--;
		} else {
			if (m->vjumps[j][i] == 0) {
				ik_ivec_push(m->vpath, m->vtrace[j][i]);
				j = m->vtrace[j][i];
				i--;
			} else {
				max_state = m->vtrace[j][i];
				for (k = 0; k < m->vjumps[j][i]; k++) {
					ik_ivec_push(m->vpath, j);
				}
				i -= m->vjumps[j][i];
				j = max_state;
			}
		}
	}
}

static void dnMatrix_viterbiP (dnMatrix m) {
	int		i, j, k, max_state;
	double	max_score;
	ik_dna	dna = m->dna;
	dnHMM	hmm = m->hmm;
	struct	dnMax max;
	
	// initialization
	i = 0;
	for (j = 0; j < hmm->states; j++) {
		if (hmm->state[j]->init == 0) {
			m->vscore[j][i] = 0;
		} else if (emission(hmm->state[j], dna->num, i) == 0) {
			m->vscore[j][i] = 0;
		} else {
			m->vscore[j][i] = 
			hmm->state[j]->init * emission(hmm->state[j], dna->num, 0);
		}
		m->vtrace[j][i] = j;
		m->vjumps[j][i] = 0;
	}
	
	// induction
	for (i = 1; i < dna->length; i++) {
		for (j = 0; j < hmm->states; j++) {
			if (hmm->state[j]->durations) max = viterbiXP_calc(m, i, j);
			else						  max = viterbiP_calc(m, i, j);
			m->vscore[j][i] = max.score;
			m->vtrace[j][i] = max.state;
			m->vjumps[j][i] = max.jumps;
		}
	}
	
	// termination
	i = dna->length -1;
	max_state = -1;
	max_score = 0;
	for (j = 0; j < hmm->states; j++) {
		if (hmm->state[j]->term == 0) {
			m->vscore[j][i] = 0;
		} else {
			m->vscore[j][i] *= hmm->state[j]->term;
		}
		if (m->vscore[j][i] > max_score) {
			max_score = m->vscore[j][i];
			max_state = j;
		}
	}
	m->vmax = max_score;
	
	if (max_state == -1) {
		ik_exit(1, "Viterbi underflowed, -no-logs was a bad idea");
	}
	
	// trace
	i = dna->length -1;
	j = max_state;			   
	ik_ivec_push(m->vpath, j);
   
	while (i > 0) {
		if (m->vtrace[j][i] == -1) {
			ik_exit(1, "fatal traceback error");
		} else if (m->vtrace[j][i] == j) {
			ik_ivec_push(m->vpath, j);
			i--;
		} else {
			if (m->vjumps[j][i] == 0) {
				ik_ivec_push(m->vpath, m->vtrace[j][i]);
				j = m->vtrace[j][i];
				i--;
			} else {
				max_state = m->vtrace[j][i];
				for (k = 0; k < m->vjumps[j][i]; k++) {
					ik_ivec_push(m->vpath, j);
				}
				i -= m->vjumps[j][i];
				j = max_state;
			}
		}
	}
}


void dnMatrix_free (dnMatrix m) {
	int i;
	
	for (i = 0; i < m->hmm->states; i++) {
		ik_free(m->vscore[i]);
		ik_free(m->vtrace[i]);
		ik_free(m->vjumps[i]);
		ik_free(m->bscore[i]);
		ik_free(m->fscore[i]);
	}
	
	ik_free(m->vscore);
	ik_free(m->vtrace);
	ik_free(m->vjumps);
	ik_free(m->bscore);
	ik_free(m->fscore);
	ik_ivec_free(m->vpath);
	ik_ivec_free(m->ppath);
}

dnMatrix dnMatrix_new (const dnHMM hmm, const ik_dna dna) {
	int		 i, j;
	dnMatrix m = ik_malloc(sizeof(struct dnMatrix));
	
	m->dna = dna;
	m->hmm = hmm;
	
	m->vscore = ik_malloc(sizeof(double) * hmm->states);
	m->vtrace = ik_malloc(sizeof(int)	* hmm->states);
	m->vjumps = ik_malloc(sizeof(int)	* hmm->states);
	m->bscore = ik_malloc(sizeof(double) * hmm->states);
	m->fscore = ik_malloc(sizeof(double) * hmm->states);
	for (i = 0; i < hmm->states; i++) {
		m->vscore[i] = ik_malloc(sizeof(double) * dna->length);
		m->vtrace[i] = ik_malloc(sizeof(int)   * dna->length);
		m->vjumps[i] = ik_malloc(sizeof(int)   * dna->length);
		m->bscore[i] = ik_malloc(sizeof(double) * dna->length);
		m->fscore[i] = ik_malloc(sizeof(double) * dna->length);
	}
	
	for (i = 0; i < hmm->states; i++) {
		for (j = 0; j < dna->length; j++) {
			m->vscore[i][j] = 0;
			m->vtrace[i][j] = -1;
			m->bscore[i][j] = 0;
			m->fscore[i][j] = 0;
		}
	}
	
	m->vpath = ik_ivec_new();
	m->vmax = 0;
	m->ppath = ik_ivec_new();
	
	return m;
}

*/

#endif
