/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2019-12-31 11:31:59
  *Edit History: 
***********************************************************/

#ifndef XDK_SAMPLE_H
#define XDK_SAMPLE_H

#include "mp.h"
#include "str.h"
#include "utils.h"

struct bio_fq_s;
typedef struct bio_fq_s bio_fq_t;

struct bio_lib_s;
typedef struct bio_lib_s bio_lib_t;

struct bio_spl_s;
typedef struct bio_spl_s bio_spl_t;

struct bio_fq_s {
  uint64_t flag;
  int32_t rd_len[2];
	str_t fq[2];
	str_t ad[2];
  str_t ad_str[2];
};
MP_DEF (bfq, bio_fq_t);

struct bio_lib_s {
	mp_t(bfq) * fqs;
	str_t * name;
};
MP_DEF (blib, bio_lib_t);

struct bio_spl_s {
	mp_t(blib) * libs;
	str_t * name;
};

/**********************************************************
 ******************** inline functions ********************
 **********************************************************/

// bio_fq_t

static inline void
bio_fq_init2 (bio_fq_t * fq)
{
  str_init2 (fq->fq);
  str_init2 (fq->fq+1);
  str_init2 (fq->ad);
  str_init2 (fq->ad+1);
  str_init2 (fq->ad_str);
  str_init2 (fq->ad_str+1);
  fq->flag = 0;
}

static inline void
bio_fq_free2 (bio_fq_t * fq)
{
	str_free2 (fq->fq);
	str_free2 (fq->fq+1);
	str_free2 (fq->ad);
	str_free2 (fq->ad+1);
	str_free2 (fq->ad_str);
	str_free2 (fq->ad_str+1);
}

// bio_lib_t

static inline void
bio_lib_init2 (bio_lib_t * lib)
{
	lib->fqs = mp_init (bfq, bio_fq_init2);
	lib->name = str_init ();
}

static inline void
bio_lib_free2 (bio_lib_t * lib)
{
	mp_free (bfq, lib->fqs, bio_fq_free2);
	str_free (lib->name);
}

// bio_spl_t

static inline void
bio_spl_init2 (bio_spl_t * spl)
{
	spl->libs = mp_init (blib, bio_lib_init2);
	spl->name = str_init ();
}

static inline void
bio_spl_free2 (bio_spl_t * spl)
{
	mp_free (blib, spl->libs, bio_lib_free2);
	str_free (spl->name);
}

/**********************************************************
 *********************** functions ************************
 **********************************************************/

#ifdef __cplusplus
extern "C" {
#endif

	bio_spl_t * bio_spl_load (const char * spl_cfg);

	int bio_spl_check (bio_spl_t * spl);

  void bio_spl_free (bio_spl_t * spl);

	void bio_spl_adapter_parse (str_t * ad, str_t * ad_str, char * s);

	void bio_spl_dump (FILE * out, bio_spl_t * spl);

#ifdef __cplusplus
}
#endif

#endif
