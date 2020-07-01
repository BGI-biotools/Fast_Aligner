/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2020-02-17 17:52:53
  *Edit History: 
***********************************************************/

#ifndef X_EBAM_H
#define X_EBAM_H

#include <stdint.h>

#include <sam.h>
#include <hts.h>

#include "mp.h"
#include "bmp.h"
#include "tmp.h"
#include "array.h"
#include "utils.h"
#include "sample.h"
#include "mem2/bwa.h"
#include "mem2/bwamem.h"
#include "mem2/bntseq.h"

/*---------------------------------------------------------------------------*/
/*----------------------------- ESAM Functions ------------------------------*/
/*---------------------------------------------------------------------------*/

struct ebam_s;
typedef struct ebam_s ebam_t;

struct aln_res_s;
typedef struct aln_res_s aln_res_t;

struct ebam_file_s;
typedef struct ebam_file_s ebam_file_t;

struct ebam_set_s;
typedef struct ebam_set_s ebam_set_t;

struct ebam_s {
  bam1_t b;

	// for tmp file output
	FILE * out;

  // for alignment sort
  uint64_t key;

  // for marking duplicates
	// end[0], end[1], and frag_end = (tid << 32) | end_pos
  int64_t end[2];
	int64_t frag_end;
	int32_t tile;
	int32_t x, y;
	int16_t score;
	int8_t lib_id;
	uint8_t pair_orit:4; // pair orientation
	uint8_t frag_orit:4; // frag orientation

	uint8_t flag;
};

#define EBAM_IS_PAIR 0x1
#define EBAM_IS_SEC  0x2

MP_DEF (erd, ebam_t);
BMP_DEF (erd, ebam_t);
ARR_DEF (erp, ebam_t*);

struct aln_res_s {
	arr_t(erp) * rds;
};

MP_DEF (ar, aln_res_t);

struct ebam_file_s {
  char * path;
  FILE * fp;
  int64_t n_reads;
};

#define N_PIPE_THREADS 2

/* inline functions */

static inline void
ebam_init2 (ebam_t * eb)
{
  memset (&eb->b, 0, sizeof(bam1_t));
  eb->end[0] = eb->end[1] = -1;
	eb->score = -1;
	eb->out = NULL;
}

static inline void
ebam_clear (ebam_t * eb)
{
  eb->b.l_data = 0;
  eb->end[0] = eb->end[1] = -1;
	eb->score = -1;
	eb->out = NULL;
}

static inline void
ebam_free2 (ebam_t * eb)
{
  if (eb->b.data != NULL)
    free (eb->b.data);
}

static inline void
aln_res_init2 (aln_res_t * res)
{
	res->rds = arr_init (erp);
}

static inline void
aln_res_clear (aln_res_t * res)
{
	arr_clear (erp, res->rds);
}

static inline void
aln_res_free2 (aln_res_t * res)
{
	arr_free (erp, res->rds);
}

#ifdef __cplusplus
extern "C" {
#endif

	ebam_set_t * ebam_set_init (void);

  int mem_aln2ebam (const mem_opt_t *opt, const bntseq_t *bns,
      bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_, int pid);

	int ebam_reset (bseq1_t * seqs, int n_seqs);

	int ebam_set_dump (bseq1_t * seqs, int n_seqs);

	void ebam_mkd_score_calc (bseq1_t * r1, bseq1_t * r2);

  int load_ebam_files (void);

	int ebam_list_dump (void);

  int ebam_list_remove (void);

	void ebam_dump (ebam_t * r);

	int load_reads (const char * file, ebam_t * reads, int32_t n_reads);

#ifdef __cplusplus
}
#endif

#endif
