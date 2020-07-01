/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2017-05-16 15:07:45
  *Edit History: 
***********************************************************/

#ifndef XDK_MP_H
#define XDK_MP_H

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include "utils.h"

#define MP_INIT_SIZE 2

#define MP_UNIQSORT 1
#define MP_SORTED   2

#define MP_DEF2(SCOPE, name, type_t) \
	struct mp_##name##_s { \
		int64_t n, m; \
    int64_t mn; \
		type_t * pool; \
		void (*init_type_f) (type_t*); \
    int32_t * multi; \
    int64_t multi_max; \
    uint64_t flag; \
	}; \
	typedef struct mp_##name##_s mp_##name##_t; \
	\
	SCOPE mp_##name##_t * mp_init_##name ( \
			void (*init_type_f)(type_t*)) { \
		int64_t i; \
		mp_##name##_t * mp; \
		mp = (mp_##name##_t *) ckalloc (1, sizeof(mp_##name##_t)); \
		mp->n = 0; \
    mp->mn = 0; \
		mp->m = MP_INIT_SIZE; \
    mp->multi_max = -1; \
    mp->flag = 0; \
		mp->init_type_f = init_type_f; \
		mp->pool = (type_t *) ckalloc (mp->m, sizeof(type_t)); \
		if (init_type_f != NULL) \
			for (i=0; i<mp->m; i++) \
				init_type_f (mp->pool+i); \
		return mp; \
	} \
	\
	SCOPE void mp_free_##name (mp_##name##_t * mp, void (*free_type_f)(type_t*)) { \
		int64_t i; \
    int64_t max; \
    if (mp->init_type_f == NULL) \
      max = mp->n>mp->mn ? mp->n : mp->mn; \
    else \
      max = mp->m; \
		if (free_type_f != NULL) { \
			for (i=0; i<max; i++) \
				free_type_f (mp->pool + i); \
		} \
    if (mp->multi_max > 0) \
      free (mp->multi); \
		free (mp->pool); \
		free (mp); \
	} \
	SCOPE void mp_clear_##name (mp_##name##_t * mp, void (*clear_type_f)(type_t*)) { \
		int64_t i; \
		if (clear_type_f != NULL) { \
			for (i=0; i<mp->m; i++) \
				clear_type_f (mp->pool + i); \
		} else { \
      memset (mp->pool, 0, mp->m*sizeof(type_t)); \
    } \
    if (mp->n > mp->mn) \
      mp->mn = mp->n; \
		mp->n = 0; \
    mp->flag = 0; \
	} \
	SCOPE void mp_resize_##name (mp_##name##_t * mp, int64_t new_size) { \
		int64_t i, old_max; \
		if (mp->m >= new_size) \
			return; \
		old_max = mp->m; \
		if (mp->m <= 0) \
			mp->m = MP_INIT_SIZE; \
		while (mp->m < new_size) { \
			if (mp->m < 0x100000) \
				mp->m <<= 1; \
			else \
				mp->m += 0x100000; \
		} \
		mp->pool = (type_t *) ckrealloc (mp->pool, mp->m*sizeof(type_t)); \
		if (mp->init_type_f == NULL) \
			memset (mp->pool+old_max, 0, (mp->m-old_max)*sizeof(type_t)); \
		else { \
			for (i=old_max; i<mp->m; i++) \
				mp->init_type_f (mp->pool+i); \
		} \
	} \
	SCOPE type_t * mp_alloc_##name (mp_##name##_t * mp) { \
		mp_resize_##name (mp, mp->n+1); \
    mp->flag = 0; \
		return mp->pool + mp->n++; \
	} \
	SCOPE int64_t mp_add_##name (mp_##name##_t * mp, type_t * src, \
			void (*copy_type_f)(type_t*,type_t*)) { \
		type_t * dst; \
		dst = mp_alloc (name, mp); \
		if (src == NULL) \
			err_mesg ("src == NULL"); \
		if (copy_type_f == NULL) \
			memcpy (dst, src, sizeof(type_t)); \
		else \
			copy_type_f (dst, src); \
    mp->flag = 0; \
		return mp->n - 1; \
	} \
	SCOPE type_t * mp_at_##name (mp_##name##_t * mp, int64_t idx) { \
		if (idx<0 || idx>=mp->n) \
			return NULL; \
		return mp->pool + idx; \
	} \
	SCOPE void mp_dump_##name (mp_##name##_t * mp, FILE * out, \
			void(*type_dump_f)(FILE*,type_t*,void*), void * data) { \
		int64_t i; \
    assert (out != NULL); \
    assert (type_dump_f != NULL); \
		for (i=0; i<mp->n; i++) \
      type_dump_f (out, mp->pool+i, data); \
	} \
  SCOPE void mp_dump_gz_##name (mp_##name##_t * mp, gzFile out, \
      void(*type_dump_f)(gzFile,type_t*,void*), void * data) { \
    int64_t i; \
    assert (out != NULL); \
    assert (type_dump_f != NULL); \
    for (i=0; i<mp->n; ++i) \
      type_dump_f (out, mp->pool+i, data); \
  } \
	SCOPE void mp_copy_##name (mp_##name##_t * dst, mp_##name##_t * src, \
			void (*copy_type_f)(type_t*,type_t*)) { \
		int64_t i; \
		mp_resize_##name (dst, src->n); \
		dst->n = src->n; \
		if (copy_type_f == NULL) { \
			memcpy (dst->pool, src->pool, src->n*sizeof(type_t)); \
		} else { \
			for (i=0; i<src->n; ++i) \
				copy_type_f (dst->pool+i, src->pool+i); \
		} \
    dst->flag = src->flag; \
	} \
  SCOPE void mp_append_##name (mp_##name##_t * dst, mp_##name##_t * src, \
      void (*copy_type_f)(type_t*,type_t*)) { \
    int64_t i; \
    int64_t n_old; \
    type_t * dst_ptr; \
    type_t * src_ptr; \
    n_old = dst->n; \
    dst->n += src->n; \
    mp_resize_##name (dst, dst->n); \
    dst_ptr = dst->pool + n_old; \
    src_ptr = src->pool; \
    if (copy_type_f == NULL) { \
      memcpy (dst_ptr, src_ptr, src->n*sizeof(type_t)); \
    } else { \
      for (i=0; i<src->n; ++i,++dst_ptr,++src_ptr) \
        copy_type_f (dst_ptr, src_ptr); \
    } \
    dst->flag = 0; \
  } \
  \
  SCOPE type_t * mp_last_##name ( \
      mp_##name##_t * mp) { \
    if (mp->n <= 0) \
      err_mesg ("memory pool is empty!"); \
    return mp->pool + mp->n-1; \
  } \
  \
  SCOPE void mp_std_sort_##name ( \
      mp_##name##_t * mp, \
      CompFunc comp_func) { \
    if (mp->flag & MP_SORTED) \
      return; \
    qsort (mp->pool, mp->n, sizeof(type_t), comp_func); \
    mp->flag |= MP_SORTED; \
  } \
  \
  SCOPE void mp_multi_mem_resize_##name ( \
      mp_##name##_t * mp) { \
    assert (mp->n >= 0); \
    if (mp->multi_max < 0) { \
      mp->multi_max = ((mp->n>>3) + 1) << 3; \
      mp->multi = (int32_t *) ckalloc (mp->multi_max, sizeof(int32_t)); \
    } else if (mp->multi_max < mp->n) { \
      mp->multi_max = ((mp->n>>3) + 1) << 3; \
      mp->multi = (int32_t *) ckrealloc (mp->multi, mp->multi_max*sizeof(int32_t)); \
    } \
  } \
  \
  SCOPE void mp_uniqsort_##name ( \
      mp_##name##_t * mp, \
      CompFunc comp_func, \
      IsEqualFunc equal_func) { \
    if (mp->n <= 1) { \
      mp_multi_mem_resize_##name (mp); \
      mp->multi[0] = 1; \
      return; \
    } \
    int64_t i, j, k; \
    type_t tmp; \
    type_t * dst; \
    type_t * src; \
    if (mp->flag & MP_UNIQSORT) \
      return; \
    qsort (mp->pool, mp->n, sizeof(type_t), comp_func); \
    mp_multi_mem_resize_##name (mp); \
    i = k = 0; \
    dst = mp->pool + i; \
    for (j=1; j<mp->n; ++j) { \
      src = mp->pool + j; \
      if (!equal_func(dst,src)) { \
        mp->multi[i]=j-k; k=j; \
        ++i, ++dst; \
        if (i != j) { \
          memcpy (&tmp, dst, sizeof(type_t)); \
          memcpy (dst, src, sizeof(type_t)); \
          memcpy (src, &tmp, sizeof(type_t)); \
        } \
      } \
    } \
    mp->multi[i] = j-k; \
    mp->n = i + 1; \
    mp->flag |= MP_UNIQSORT; \
  } \
  SCOPE int64_t mp_bisearch_##name ( \
      mp_##name##_t * mp, \
      CompFunc comp_func, \
      type_t * target, \
      int64_t beg, \
      int64_t end) { \
    int ret; \
    int64_t high, low, mid;\
    type_t * t; \
    low = beg>=0 ? beg : 0; \
    high = end>=0 ? end-1 : mp->n-1; \
    while (low <= high) { \
      mid = (low + high) / 2; \
      t = mp->pool + mid; \
      ret = comp_func (target, t); \
      if (ret == 0) \
        return mid; \
      if (ret < 0) \
        high = mid; \
      else \
        low = mid; \
    } \
    return -1; \
  } \
  \
  static inline int mp_is_sorted_##name ( \
      mp_##name##_t * mp) { \
    return mp->flag & MP_SORTED; \
  } \
  \
  static inline int mp_is_uniqsorted_##name ( \
      mp_##name##_t * mp) { \
    return mp->flag & MP_UNIQSORT; \
  } \
	static int __mp_##name##_def_end__ = 1

#define MP_DEF(name, type_t) \
	MP_DEF2(static inline, name, type_t)

#define mp_t(name) mp_##name##_t
#define mp_cnt(mp) ((mp)->n)

#define mp_init(name,func) mp_init_##name((func))
#define mp_free(name,mp,f) mp_free_##name((mp), (f))
#define mp_clear(name,mp,f) mp_clear_##name((mp), (f))
#define mp_resize(name,mp,size) mp_resize_##name((mp), (size))
#define mp_alloc(name,mp) mp_alloc_##name(mp)
#define mp_add(name,mp,ptr,copy_f) mp_add_##name((mp), (ptr), (copy_f))
#define mp_at(name,mp,idx) mp_at_##name((mp), (idx))
#define mp_dump(name,mp,fp,f,d) mp_dump_##name((mp), (fp), (f), (d))
#define mp_dump_gz(name,mp,fp,f,d) mp_dump_gz_##name((mp), (fp), (f), (d))
#define mp_copy(name,dst,src,f) mp_copy_##name((dst), (src), (f))
#define mp_append(name,dst,src,f) mp_append_##name((dst), (src), (f))
#define mp_std_sort(name,mp,f) mp_std_sort_##name((mp),(f))
#define mp_uniqsort(name,mp,cmp_f,equ_f) mp_uniqsort_##name((mp),(cmp_f),(equ_f))
#define mp_bisearch(name,mp,cmp_f,t,beg,end) mp_bisearch_##name((mp),(cmp_f),(t),(beg),(end))
#define mp_is_sorted(name,mp) mp_is_sorted_##name(mp)
#define mp_is_uniqsorted(name,mp) mp_is_uniqsorted_##name(mp)
#define mp_last(name,mp) mp_last_##name(mp)

#endif
