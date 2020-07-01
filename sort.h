/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2017-08-11 16:06:40
  *Edit History: 
***********************************************************/

#ifndef XDK_SORT
#define XDK_SORT

#include <stdint.h>
#include <stdlib.h>

#include "mp.h"
#include "array.h"
#include "utils.h"

#define NUMBASES 256
#define STEP     8

/*
 * merge sort
 *
 * radix sort
 *   requirement:
 *     1. type_t must contain an segment 'key' with the type 'uint64_t'
 */

#define XSORT_DEF(name, type_t, comp_type_f) \
	MP_DEF (xsort##name, type_t); \
	ARR_DEF (xsort_ptr##name, type_t*); \
	struct xsort_##name##_s { \
		mp_t(xsort##name) * pool; \
		arr_t(xsort_ptr##name) * ptrs[2]; \
		void (*init_type_f) (type_t*); \
		void (*copy_type_f) (type_t*,type_t*); \
    size_t * remainders; \
	}; \
	typedef struct xsort_##name##_s xsort_##name##_t; \
	\
	static xsort_##name##_t * xsort_init_##name ( \
			void (*init_type_f)(type_t*), \
			void (*copy_type_f)(type_t*,type_t*)) { \
		xsort_##name##_t * xs; \
		xs = (xsort_##name##_t *) ckalloc (1, sizeof(xsort_##name##_t)); \
		xs->pool = mp_init (xsort##name, init_type_f); \
		xs->ptrs[0] = arr_init (xsort_ptr##name); \
		xs->ptrs[1] = arr_init (xsort_ptr##name); \
		xs->init_type_f = init_type_f; \
		xs->copy_type_f = copy_type_f; \
    xs->remainders = (size_t *) ckmalloc (NUMBASES * sizeof(size_t)); \
		return xs; \
	} \
	\
	static void xsort_free_##name ( \
			xsort_##name##_t * xs, \
			void (*free_type_f) (type_t*)) { \
		arr_free (xsort_ptr##name, xs->ptrs[0]); \
		arr_free (xsort_ptr##name, xs->ptrs[1]); \
		mp_free (xsort##name, xs->pool, free_type_f); \
		free (xs); \
	} \
	\
	static int xmerge_sort_##name (xsort_##name##_t * xs, type_t * array, int32_t n) { \
		int ret, curr=0; \
		int32_t l, m, r; \
		int32_t i, j, k, w, dw; \
		type_t ** a; \
    type_t ** b; \
    type_t * buf; \
		if (comp_type_f == NULL) \
			err_mesg ("compare function must be set!"); \
		if (n <= 1) \
			return 0; \
		mp_resize (xsort##name, xs->pool, n); \
		arr_resize (xsort_ptr##name, xs->ptrs[0], n); \
		arr_resize (xsort_ptr##name, xs->ptrs[1], n); \
		for (i=0; i<n; ++i) \
			xs->ptrs[0]->arr[i] = array + i; \
		for (w=1; w<n; w<<=1) { \
		  a = xs->ptrs[curr]->arr; \
		  b = xs->ptrs[1-curr]->arr; \
			dw = w << 1; \
			for (l=0; l+w<n; l+=dw) { \
				m = l + w; \
				r = m + w; \
				if (r > n) \
					r = n; \
				k=l, i=l, j=m; \
				while (i<m && j<r) { \
					if (comp_type_f(a[i],a[j]) <= 0) \
						b[k++] = a[i++]; \
					else \
						b[k++] = a[j++]; \
				} \
				while (i < m) \
					b[k++] = a[i++]; \
				while (j < r) \
					b[k++] = a[j++]; \
			} \
			for ( ; l<n; ++l) \
				b[l] = a[l]; \
      curr = 1 - curr; \
		} \
		buf = xs->pool->pool; \
		if (xs->copy_type_f == NULL) { \
			for (i=0; i<n; ++i) \
				memcpy (buf+i, b[i], sizeof(type_t)); \
			memcpy (array, buf, n*sizeof(type_t));  \
		} else { \
			for (i=0; i<n; ++i) \
				xs->copy_type_f (buf+i, b[i]); \
			for (i=0; i<n; ++i) \
				xs->copy_type_f (array+i, buf+i); \
		} \
    return 0; \
	} \
  \
  static int xradix_sort_##name ( \
			xsort_##name##_t * xs, \
			type_t * array, \
			int32_t n, \
			uint64_t (*get_key)(type_t*)) { \
    int shift=0, curr=0, max_digit=0; \
    size_t * r = xs->remainders; \
    int32_t i, j; \
		uint64_t key; \
		uint64_t key_max=0; \
    type_t ** a; \
    type_t ** b; \
    type_t * buf; \
		if (n <= 1) \
			return 0; \
    mp_resize (xsort##name, xs->pool, n); \
    arr_resize (xsort_ptr##name, xs->ptrs[0], n); \
    arr_resize (xsort_ptr##name, xs->ptrs[1], n); \
    for (i=0; i<n; ++i) { \
      xs->ptrs[0]->arr[i] = array + i; \
			key = get_key (array + i); \
			if (key > key_max) \
				key_max = key; \
		} \
		while (key_max) { \
			++max_digit; \
			key_max >>= 1; \
		} \
    while (shift < max_digit) { \
      memset (r, 0, NUMBASES*sizeof(size_t)); \
      a = xs->ptrs[curr]->arr; \
      b = xs->ptrs[1-curr]->arr; \
      for (i=0; i<n; ++i) \
        ++ (r[(get_key(a[i])>>shift)%NUMBASES]); \
      for (i=1; i<NUMBASES; ++i) \
        r[i] += r[i-1]; \
      for (i=n-1; i>=0; --i) { \
        j = --(r[(get_key(a[i])>>shift)%NUMBASES]); \
        b[j] = a[i]; \
      } \
      shift += STEP; \
      curr = 1 - curr; \
    } \
    buf = xs->pool->pool; \
    if (xs->copy_type_f == NULL) { \
      for (i=0; i<n; ++i) \
        memcpy (buf+i, b[i], sizeof(type_t)); \
      memcpy (array, buf, n*sizeof(type_t)); \
    } else { \
      for (i=0; i<n; ++i) \
        xs->copy_type_f (buf+i, b[i]); \
      for (i=0; i<n; ++i) \
        xs->copy_type_f (array+i, buf+i); \
    } \
    return 0; \
  } \
  \
	static int __xsort_def_end_##name##__ = 1

#define XSORT_LITE_DEF(name,type_t,comp_type_f) \
  static int xinsert_sort_lite_##name (type_t * array, int32_t n) { \
    int32_t i, j, k; \
    type_t* p; \
    type_t* q; \
    type_t t; \
    for (i=1; i<n; ++i) { \
      p = array + i; \
      for (j=i-1; j>=0; --j) { \
        q = array + j; \
        if (comp_type_f(p,q) >= 0) \
          break; \
      } \
      for (k=i-1; k>j; --k) { \
        t = array[k+1]; \
        array[k+1] = array[k]; \
        array[k] = t; \
      } \
    } \
    return 0; \
  } \
  \
  static int xheap_adjust_lite_##name (type_t * array, int32_t n, int32_t i) { \
    int l, r, m; \
    type_t* p; \
    type_t t; \
    l = (i << 1) + 1; \
    r = l + 1; \
    while (l < n) { \
      if (r<n && comp_type_f(array+r,array+i)>0) \
        m = r; \
      else \
        m = i; \
      if (comp_type_f(array+l,array+m) > 0) \
        m = l; \
      if (m == i) \
        break; \
      t = array[i]; \
      array[i] = array[m]; \
      array[m] = t; \
      i = m; \
      l = (i << 1) + 1; \
      r = l + 1; \
    } \
    return 0; \
  } \
  \
  static int xheap_make_lite_##name (type_t * array, int32_t n) { \
    int32_t i; \
    for (i=(n>>1)-1; i>=0; --i) \
      xheap_adjust_lite_##name (array, n, i); \
    return 0; \
  } \
  \
	static int __xsort_lite_def_end_##name##__ = 1

#define xsort_t(name) xsort_##name##_t

#define xsort_init(name,init_f,copy_f) xsort_init_##name(init_f,copy_f)
#define xsort_free(name,xs,free_f) xsort_free_##name((xs),(free_f))
#define xmerge_sort(name,xs,a,n) xmerge_sort_##name((xs),(a),(n))
#define xradix_sort(name,xs,a,n,key_f) xradix_sort_##name((xs),(a),(n),(key_f))

#define xinsert_sort_lite(name,arr,n) xinsert_sort_lite_##name((arr),(n))
#define xheap_make_lite(name,arr,n) xheap_make_lite_##name((arr),(n))
#define xheap_adjust_lite(name,arr,n,i) xheap_adjust_lite_##name((arr),(n),(i))

#endif
