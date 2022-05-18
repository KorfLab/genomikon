/******************************************************************************\
 isoform.c
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef ISOFORM_C
#define ISOFORM_C

#include "genomikon.h"
#include "isoform.h"

/*************\
 Private Area
\*************/

static void combo(gkn_vec ans, gkn_ivec tmp, int n, int left, int k) {
	if (k == 0) {
		gkn_ivec keep = gkn_ivec_new();
		for (int i = 0; i < tmp->size; i++) gkn_ivec_push(keep, tmp->elem[i]);
		gkn_vec_push(ans, (void*)keep);
		return;
	}

	for (int i = left; i <= n; i++) {
		gkn_ivec_push(tmp, i - 1);
		combo(ans, tmp, n, i + 1, k - 1);
		gkn_ivec_pop(tmp);
	}
}

static gkn_vec get_combinations(const gkn_ivec sites, int k) {
	gkn_vec  indexes = gkn_vec_new();
	gkn_ivec tmp = gkn_ivec_new();

	combo(indexes, tmp, sites->size, 1, k);
	for (int i = 0; i < indexes->size; i++) {
		gkn_ivec iv = indexes->elem[i];
		for (int j = 0; j < iv->size; j++) {
			int idx = iv->elem[j];
			int val = sites->elem[idx];
			iv->elem[j] = val;
		}
	}
	gkn_ivec_free(tmp);

	return indexes;
}

static gkn_mRNA build_mRNA(
	const char *seq, int beg, int end,
	const gkn_ivec dons, const gkn_ivec accs)
{
	assert(beg <= end);
	gkn_mRNA tx = gkn_mRNA_new(seq, beg, end);
	tx->exons   = gkn_vec_new();
	tx->introns = gkn_vec_new();
	tx->score   = 0;

	assert(dons->size == accs->size);

	if (dons->size == 0) {
		gkn_vec_push(tx->exons, gkn_feat_new(seq, beg, end));
		return tx;
	}

	// introns
	for (int i = 0; i < dons->size; i++) {
		int beg = dons->elem[i];
		int end = accs->elem[i];
		gkn_feat f = gkn_feat_new(seq, beg, end);
		gkn_vec_push(tx->introns, f);
	}

	// exons
	gkn_feat ei = gkn_feat_new(seq, beg, dons->elem[0] -1);
	gkn_vec_push(tx->exons, ei);
	for (int i = 1; i < dons->size; i++) {
		int a = accs->elem[i-1] +1;
		int b = dons->elem[i] -1;
		gkn_feat ex = gkn_feat_new(seq, a, b);
		gkn_vec_push(tx->exons, ex);
	}
	gkn_feat et = gkn_feat_new(seq, accs->elem[accs->size -1] + 1, end);
	gkn_vec_push(tx->exons, et);

	return tx;
}

static int short_intron(const gkn_ivec dons, const gkn_ivec accs, int min) {
	for (int i = 0; i < dons->size; i++) {
		int len = accs->elem[i] - dons->elem[i] + 1;
		if (len < min) return 1;
	}
	return 0;
}

static int short_exon(
	const gkn_ivec dons, const gkn_ivec accs,
	int seqlen, int flank, int min)
{

	// first exon
	int exon_beg = flank + 1;
	int exon_end = dons->elem[0] -1;
	int exon_len = exon_end - exon_beg + 1;
	if (exon_len < min) return 1;

	// last exon
	exon_beg = accs->elem[accs->size -1] + 1;
	exon_end = seqlen - flank -1;
	exon_len = exon_end - exon_beg + 1;
	if (exon_len < min) return 1;

	// internal exons
	for (int i = 1; i < dons->size; i++) {
		int beg = accs->elem[i-1] + 1;
		int end = dons->elem[i] -1;
		int len = end - beg;
		if (len < min) return 1;
	}
	return 0;
}

static int canonical(const gkn_gff gff, const char *seq) {
	if (seq[gff->beg]   != 'G') return 0;
	if (seq[gff->beg+1] != 'T') return 0;
	if (seq[gff->end-1] != 'A') return 0;
	if (seq[gff->end]   != 'G') return 0;
	return 1;
}

static int intsort(const void *a, const void *b) {
	return ( *(int*)a - *(int*)b );
}

static int not_found(gkn_ivec vec, int val) {
	for (int i = 0; i < vec->size; i++) {
		if (vec->elem[i] == val) return 0;
	}
	return 1;
}

static void gff_sites(const char *file, gkn_ivec dons, gkn_ivec accs,
		const char *seq, int gtag) {
	gkn_gff gff;
	gkn_pipe io = gkn_pipe_open(file, "r");
	while ((gff = gkn_gff_read(io)) != NULL) {
		if (strcmp(gff->type, "intron") != 0) continue;
		if (gff->strand != '+') continue;
		if (gtag && !canonical(gff, seq)) continue;
		if (not_found(dons, gff->beg)) gkn_ivec_push(dons, gff->beg);
		if (not_found(accs, gff->end)) gkn_ivec_push(accs, gff->end);
		gkn_gff_free(gff);
	}

	qsort(dons->elem, dons->size, sizeof(int), intsort);
	qsort(accs->elem, accs->size, sizeof(int), intsort);

	gkn_pipe_close(io);
}

static void gtag_sites(
	const char *seq,
	int beg, int end,
	gkn_ivec dons, gkn_ivec accs)
{
	for (int i = beg; i < end; i++) {
		if (seq[i] == 'G' && seq[i+1] == 'T') gkn_ivec_push(dons, i);
		if (seq[i] == 'A' && seq[i+1] == 'G') gkn_ivec_push(accs, i+1);
	}
}

/************\
 Public Area
\************/

void isozone_free(isozone iz) {
	// descend into mRNA and free all of that shit
	gkn_vec_free(iz->mRNAs);
	free(iz);
}

isozone isoforms(
	const char *seq,
	int emin, int imin, int smax, int gen,
	const char *gff, int count_only)
{

	// inits
	int len = strlen(seq);
	gkn_ivec dons = gkn_ivec_new();
	gkn_ivec accs = gkn_ivec_new();

	if   (gff) gff_sites(gff, dons, accs, seq, 1); // canonical only for now
	else       gtag_sites(seq, gen +emin, len -gen -emin, dons, accs);

	int nsites = dons->size < accs->size ? dons->size : accs->size;
	if (nsites > smax) nsites = smax;

	// main loop
	int trials = 0;
	int forms = 0;
	gkn_vec txs = gkn_vec_new();
	for (int n = 1; n <= nsites; n++) {

		// create combos
		gkn_vec dcombos = get_combinations(dons, n);
		gkn_vec acombos = get_combinations(accs, n);

		// create isoforms
		for (int i = 0; i < dcombos->size; i++) {
			for (int j = 0; j < acombos->size; j++) {
				gkn_ivec dsites = dcombos->elem[i];
				gkn_ivec asites = acombos->elem[j];
				assert(dsites->size == asites->size);

				trials++;
				if (short_intron(dsites, asites, imin)) continue;
				if (short_exon(dsites, asites, len, gen, emin)) continue;
				
				forms++;
				if (count_only) continue;

				// save isoform
				gkn_mRNA tx = build_mRNA(seq, gen, len - gen, dsites, asites);
				gkn_vec_push(txs, tx);
			}
		}

		// free combos
		for (int i = 0; i < dcombos->size; i++) {
			gkn_ivec v = dcombos->elem[i];
			gkn_ivec_free(v);
		}
		gkn_vec_free(dcombos);
		for (int i = 0; i < acombos->size; i++) {
			gkn_ivec v = acombos->elem[i];
			gkn_ivec_free(v);
		}
		gkn_vec_free(acombos);
		
	}

	// return values
	isozone iz = malloc(sizeof(struct isoform_zone));
	iz->dons     = dons->size;
	iz->accs     = accs->size;
	iz->trials   = trials;
	iz->isoforms = forms;
	iz->mRNAs    = txs;

	// clean up
	gkn_ivec_free(dons);
	gkn_ivec_free(accs);

	return iz;
}

#endif
