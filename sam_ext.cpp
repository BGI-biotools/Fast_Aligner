/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2019-07-02 17:43:22
  *Edit History: 
***********************************************************/

#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include <sam.h>
#include <faidx.h>
#include <khash_str2int.h>

#include "var.h"
#include "bedidx.h"
#include "sam_bio.h"
#include "sam_ext.h"

int mplp_get_ref(mplp_aux_t *ma, int tid,  char **ref, int *ref_len) {
    mplp_ref_t *r = ma->ref;

    //printf("get ref %d {%d/%p, %d/%p}\n", tid, r->ref_id[0], r->ref[0], r->ref_id[1], r->ref[1]);

    if (!r || !ma->conf->fai) {
        *ref = NULL;
        return 0;
    }

    // Do we need to reference count this so multiple mplp_aux_t can
    // track which references are in use?
    // For now we just cache the last two. Sufficient?
    if (tid == r->ref_id[0]) {
        *ref = r->ref[0];
        *ref_len = r->ref_len[0];
        return 1;
    }
    if (tid == r->ref_id[1]) {
        // Last, swap over
        int tmp;
        tmp = r->ref_id[0];  r->ref_id[0]  = r->ref_id[1];  r->ref_id[1]  = tmp;
        tmp = r->ref_len[0]; r->ref_len[0] = r->ref_len[1]; r->ref_len[1] = tmp;

        char *tc;
        tc = r->ref[0]; r->ref[0] = r->ref[1]; r->ref[1] = tc;
        *ref = r->ref[0];
        *ref_len = r->ref_len[0];
        return 1;
    }

    // New, so migrate to old and load new
    free(r->ref[1]);
    r->ref[1]     = r->ref[0];
    r->ref_id[1]  = r->ref_id[0];
    r->ref_len[1] = r->ref_len[0];

    r->ref_id[0] = tid;
    r->ref[0] = faidx_fetch_seq(ma->conf->fai,
                                ma->h->target_name[r->ref_id[0]],
                                0,
                                INT_MAX,
                                &r->ref_len[0]);

    if (!r->ref[0]) {
        r->ref[0] = NULL;
        r->ref_id[0] = -1;
        r->ref_len[0] = 0;
        *ref = NULL;
        return 0;
    }

    *ref = r->ref[0];
    *ref_len = r->ref_len[0];
    return 1;
}

int mplp_func(void *data, bam1_t *b)
{
    char *ref;
    mplp_aux_t *ma = (mplp_aux_t*)data;
    int ret, skip = 0, ref_len;
    do {
        int has_ref;
        ret = ma->iter? sam_itr_next(ma->fp, ma->iter, b) : sam_read1(ma->fp, ma->h, b);
        if (ret < 0) break;
        // The 'B' cigar operation is not part of the specification, considering as obsolete.
        //  bam_remove_B(b);
        if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) { // exclude unmapped reads
            skip = 1;
            continue;
        }
        if (ma->conf->rflag_require && !(ma->conf->rflag_require&b->core.flag)) { skip = 1; continue; }
        if (ma->conf->rflag_filter && ma->conf->rflag_filter&b->core.flag) { skip = 1; continue; }
        if (ma->conf->bed && ma->conf->all == 0) { // test overlap
            skip = !bed_overlap(ma->conf->bed, ma->h->target_name[b->core.tid], b->core.pos, bam_endpos(b));
            if (skip) continue;
        }
        if (ma->conf->rghash) { // exclude read groups
            uint8_t *rg = bam_aux_get(b, "RG");
            skip = (rg && khash_str2int_get(ma->conf->rghash, (const char*)(rg+1), NULL)==0);
            if (skip) continue;
        }
        if (ma->conf->flag & MPLP_ILLUMINA13) {
            int i;
            uint8_t *qual = bam_get_qual(b);
            for (i = 0; i < b->core.l_qseq; ++i)
                qual[i] = qual[i] > 31? qual[i] - 31 : 0;
        }

        if (ma->conf->fai && b->core.tid >= 0) {
            has_ref = mplp_get_ref(ma, b->core.tid, &ref, &ref_len);
            if (has_ref && ref_len <= b->core.pos) { // exclude reads outside of the reference sequence
                fprintf(stderr,"[%s] Skipping because %d is outside of %d [ref:%d]\n",
                        __func__, b->core.pos, ref_len, b->core.tid);
                skip = 1;
                continue;
            }
        } else {
            has_ref = 0;
        }

        skip = 0;
        if (has_ref && (ma->conf->flag&MPLP_REALN)) sam_prob_realn(b, ref, ref_len, (ma->conf->flag & MPLP_REDO_BAQ)? 7 : 3);
        if (has_ref && ma->conf->capQ_thres > 10) {
            int q = sam_cap_mapq(b, ref, ref_len, ma->conf->capQ_thres);
            if (q < 0) skip = 1;
            else if (b->core.qual > q) b->core.qual = q;
        }
        if (b->core.qual < ma->conf->min_mq) skip = 1;
        else if ((ma->conf->flag&MPLP_NO_ORPHAN) && (b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) skip = 1;
    } while (skip);
    return ret;
}

// added

void
mplp_conf_init (mplp_conf_t * conf)
{
  memset(conf, 0, sizeof(mplp_conf_t));
  conf->min_mq = 13; // -q 13
  conf->min_baseQ = 13;
  conf->capQ_thres = 0;
  conf->max_depth = MPLP_MAX_DEPTH;
  conf->max_indel_depth = MPLP_MAX_INDEL_DEPTH;
  conf->openQ = 40; conf->extQ = 20; conf->tandemQ = 100;
  conf->min_frac = 0.002; conf->min_support = 1;
  conf->flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_SMART_OVERLAPS;
  conf->flag |= MPLP_PRINT_MAPQ; // -O
  conf->flag |= MPLP_PRINT_POS; // -s
  conf->argc = -1; conf->argv = NULL;
  conf->rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
  conf->output_fname = NULL;
  conf->all = 0;
}

static inline mplp_aux_t *
mpileup_aux_init (char * bam, mplp_conf_t * conf)
{
	mplp_aux_t * aux;

  aux = (mplp_aux_t *) ckalloc (1, sizeof(mplp_aux_t));
  aux->fp = xsam_open (bam, "rb");

  if (hts_set_opt(aux->fp, CRAM_OPT_DECODE_MD, 0))
    err_mesg ("failed to set CRAM_OPT_DECODE_MD value");
  if (conf->fai_fname && hts_set_fai_filename(aux->fp,conf->fai_fname)!=0)
    err_mesg ("failed to set fai file name");

  aux->conf = conf;

  aux->ref = (mplp_ref_t *) ckalloc (1, sizeof(mplp_ref_t));
  aux->ref->ref[0] = aux->ref->ref[1] = NULL;
  aux->ref->ref_id[0] = aux->ref->ref_id[1] = -1;
  aux->ref->ref_len[0] = aux->ref->ref_len[1] = 0;

	return aux;
}

static inline void
mpileup_aux_free (mplp_aux_t * aux)
{
	sam_close (aux->fp);
	free (aux->ref);
	free (aux);
}

mplp_aux_t **
mplp_aux_init (char * ref_file, int n_bams, char ** bams)
{
  int i;
  bam_hdr_t * bam_hdr;
  mplp_conf_t * conf;
  mplp_aux_t ** aux;

  conf = (mplp_conf_t *) ckmalloc (sizeof(mplp_conf_t));
  mplp_conf_init (conf);
  conf->fai = fai_load (ref_file);
  if (conf->fai == NULL)
    err_mesg ("[%s] fail to load ref.fa.fai!", __func__);
  conf->fai_fname = ref_file;

  aux = (mplp_aux_t **) ckalloc (n_bams, sizeof(mplp_aux_t *));
  for (i=0; i<n_bams; ++i)
    aux[i] = mpileup_aux_init (bams[i], conf);
  bam_hdr = xsam_hdr_read (aux[0]->fp);
  for (i=0; i<n_bams; ++i)
    aux[i]->h = bam_hdr;

  return aux;
}

void
mplp_aux_free (mplp_aux_t ** aux, int n_bams)
{
	int i;

	bam_hdr_destroy (aux[0]->h);
	free ((void*)aux[0]->conf);
	fai_destroy (aux[0]->conf->fai);

	for (i=0; i<n_bams; ++i)
		mpileup_aux_free (aux[i]);
	free (aux);
}

int
next2useful_softclip (const bam_pileup1_t * p, int min_base_qual)
{
  uint8_t * bq;
  int32_t l_seq;

  l_seq = p->b->core.l_qseq;
  if (p->qpos >= l_seq)
    return 0;

  bq = bam_get_qual (p->b);
  if (bq[p->qpos] <= min_base_qual)
    return 0;

  if (p->is_head && p->qpos>0
      && bq[p->qpos-1] > min_base_qual)
    return 1;

  if (p->is_tail && p->qpos<l_seq-1
      && bq[p->qpos+1] > min_base_qual)
    return 1;

  return 0;
}
