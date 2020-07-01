/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-07-13 15:48:40
  *Edit History: 
***********************************************************/

#ifndef SAM_UTILS_H
#define SAM_UTILS_H

#include <sam.h>

#include "mp.h"

/*----------------------------------------------------------------------------*/
/*----------------------------[ Pool for bam1_t ]-----------------------------*/
/*----------------------------------------------------------------------------*/

MP_DEF (hts_sam, bam1_t);
typedef mp_t(hts_sam) hts_sam_pool_t;

static inline void hts_sam_init2 (bam1_t * b, void * data)
{
	memset (b, 0, sizeof(bam1_t));
}

static inline void hts_sam_free2 (bam1_t * b)
{
	free (b->data);
}

#define hts_sam_pool_cnt(pool) (mp_cnt(pool))

#define hts_sam_pool_init() (mp_init(hts_sam, hts_sam_init2, NULL))
#define hts_sam_pool_free(pool) (mp_free(hts_sam, (pool), hts_sam_free2))
#define hts_sam_pool_clear(pool) (mp_clear(hts_sam, (pool), NULL))
#define hts_sam_pool_resize(pool,size) (mp_resize(hts_sam, (pool), (size)))
#define hts_sam_pool_alloc(pool) (mp_alloc(hts_sam, (pool)))
#define hts_sam_pool_at(pool,idx) (mp_at(hts_sam, (pool), (idx)))

#endif
