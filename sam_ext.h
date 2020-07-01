/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2019-07-02 18:04:22
  *Edit History: 
***********************************************************/

#ifndef XDK_SAM_EXT_H
#define XDK_SAM_EXT_H

#include <sam.h>
#include <faidx.h>

#include "mp.h"

/*---------------------------------------------------------------------------*/
/*-------------------------- Read Piles Function ----------------------------*/
/*---------------------------------------------------------------------------*/

#define MPLP_BCF        1
#define MPLP_VCF        (1<<1)
#define MPLP_NO_COMP    (1<<2)
#define MPLP_NO_ORPHAN  (1<<3)
#define MPLP_REALN      (1<<4)
#define MPLP_NO_INDEL   (1<<5)
#define MPLP_REDO_BAQ   (1<<6)
#define MPLP_ILLUMINA13 (1<<7)
#define MPLP_IGNORE_RG  (1<<8)
#define MPLP_PRINT_POS  (1<<9)
#define MPLP_PRINT_MAPQ (1<<10)
#define MPLP_PER_SAMPLE (1<<11)
#define MPLP_SMART_OVERLAPS (1<<12)
#define MPLP_PRINT_QNAME (1<<13)

//#define MPLP_MAX_DEPTH 8000
//#define MPLP_MAX_INDEL_DEPTH 250
// modified
#define MPLP_MAX_DEPTH 200000
#define MPLP_MAX_INDEL_DEPTH 10000

typedef struct {
    int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag, all;
    int rflag_require, rflag_filter;
    int openQ, extQ, tandemQ, min_support; // for indels
    double min_frac; // for indels
    char *reg, *pl_list, *fai_fname, *output_fname;
    faidx_t *fai;
    void *bed, *rghash;
    int argc;
    char **argv;
} mplp_conf_t;

typedef struct {
    char *ref[2];
    int ref_id[2];
    int ref_len[2];
} mplp_ref_t;

#define MPLP_REF_INIT {{NULL,NULL},{-1,-1},{0,0}}

typedef struct {
    samFile *fp;
    hts_itr_t *iter;
    bam_hdr_t *h;
    mplp_ref_t *ref;
    const mplp_conf_t *conf;
} mplp_aux_t;

typedef struct {
    int n;
    int *n_plp, *m_plp;
    bam_pileup1_t **plp;
} mplp_pileup_t;

#ifdef __cplusplus
extern "C" {
#endif

  int mplp_get_ref(mplp_aux_t *ma, int tid,  char **ref, int *ref_len);
  int mplp_func(void *data, bam1_t *b);

  // added
  void mplp_conf_init (mplp_conf_t * conf);
  mplp_aux_t ** mplp_aux_init (char * ref_file, int n_bams, char ** bams);
	void mplp_aux_free (mplp_aux_t ** aux, int n_bams);
  int next2useful_softclip (const bam_pileup1_t * p, int min_base_qual);

#ifdef __cplusplus
}
#endif

#endif
