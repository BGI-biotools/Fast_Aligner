/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-11-02 16:16:44
  *Edit History: 
***********************************************************/

#ifndef _XDK_HASH_H
#define _XDK_HASH_H

/*
 * Warning: These hash functions are not thread safe
 */

#include <stdint.h>

#include "bmp.h"
#include "utils.h"

#define XH_FAIL 0
#define XH_EXIST 1
#define XH_NEW 2

typedef uint64_t (*HashFunc) (const void*);

struct xh_item_s;
typedef struct xh_item_s xh_item_t;

struct xh_s;
typedef struct xh_s xh_t;

struct xh_item_s {
  void * key, * val;
  uint64_t hash_val;
  int64_t id;
  uint32_t multi:31;
  uint32_t deleted:1;
  xh_item_t * next;
};

struct xh_s {
  uint64_t size;
  uint64_t max;
  uint64_t cnt;
  uint64_t del_c;
  double load_factor;
  HashFunc hash_func;
  IsEqualFunc is_equal_func;
  xh_item_t * pool;
  xh_item_t ** slots;
};

/**********************************************************
 ********************* Shared Functions *******************
 **********************************************************/

xh_t * _xh_init (int64_t size, double load_factor, HashFunc hash_f, IsEqualFunc is_equal_f);
void _xh_clear (xh_t * h);
void _xh_free (xh_t * h);

/**********************************************************
 ************************ Hash Set ************************
 **********************************************************/

#ifdef __cplusplus
extern "C" {
#endif

/*
 * if key is new added to set:
 *   return XH_NEW
 * if key already exist in set:
 *   return XH_EXIST
 * if key fail to add
 *   return XH_FAIL
 */
int _xh_set_add (xh_t * h, void * key);

/*
 * if key is new added to set:
 *   *old_key==NULL and return XH_NEW
 * if key already exist in set:
 *   *old_key==old_key_in_hash and return XH_EXIST
 * if key fail to add
 *   return XH_FAIL
 */
int _xh_set_add2 (xh_t * h, void * new_key, void ** old_key);

/*
 * if key is in the set:
 *   return XH_EXIST
 * if key is not found in the set:
 *   return XH_FAIL
 */
int _xh_set_search (xh_t * h, void * key);

/*
 * if key is in the set:
 *   return the existed original key
 * if key is not found in the set:
 *   return NULL
 */
void * _xh_set_search2 (xh_t * h, void * key);

/*
 * if key is in the set:
 *   return the existed original xh_item_t
 * if key is not found in the set:
 *   return NULL
 */
xh_item_t * _xh_set_search3 (xh_t * h, void * key);

#ifdef __cplusplus
}
#endif

#define HASH_SET_DEF(name, key_t) \
  BMP_DEF (xh_set_key##name, key_t); \
  struct xh_set_##name##_s { \
    xh_t * hash; \
    bmp_t(xh_set_key##name) * key_bmp; \
    void (*key_copy_f) (key_t*, key_t*); \
  }; \
  typedef struct xh_set_##name##_s xh_set_##name##_t; \
  \
  static inline xh_set_##name##_t * xh_set_init_##name ( \
      int64_t size, double load_factor, \
      void (*key_init_f)(key_t*), \
      void (*key_copy_f)(key_t*,key_t*), \
      HashFunc hash_f, \
      IsEqualFunc is_equal_f) { \
    xh_set_##name##_t * xh; \
    xh = (xh_set_##name##_t *) ckmalloc (sizeof(xh_set_##name##_t)); \
    xh->hash = _xh_init (size, load_factor, hash_f, is_equal_f); \
    xh->key_bmp = bmp_init (xh_set_key##name, size, key_init_f); \
    xh->key_copy_f = key_copy_f; \
    return xh; \
  } \
  \
  static inline void xh_set_clear_##name ( \
      xh_set_##name##_t * xh, \
      void (*key_clear_f) (key_t*)) { \
    bmp_clear (xh_set_key##name, xh->key_bmp, key_clear_f); \
    _xh_clear (xh->hash); \
  } \
  \
  static inline void xh_set_free_##name ( \
      xh_set_##name##_t * xh, \
      void (*key_free_f) (key_t*)) { \
    bmp_free (xh_set_key##name, xh->key_bmp, key_free_f); \
    _xh_free (xh->hash); \
    free (xh); \
  } \
  \
  static inline int xh_set_add_##name ( \
      xh_set_##name##_t * xh, \
      key_t * key) { \
    int ret; \
    key_t * new_key; \
    new_key = bmp_alloc (xh_set_key##name, xh->key_bmp); \
    if (xh->key_copy_f == NULL) \
      memcpy (new_key, key, sizeof(key_t)); \
    else \
      xh->key_copy_f (new_key, key); \
    ret = _xh_set_add (xh->hash, (void*)new_key); \
    if (ret != XH_NEW) \
      bmp_pop (xh_set_key##name, xh->key_bmp); \
    return ret; \
  } \
  static inline key_t * xh_set_add2_##name ( \
      xh_set_##name##_t * xh, \
      key_t * key) { \
    int ret; \
    key_t * new_key; \
    key_t * old_key; \
    new_key = bmp_alloc (xh_set_key##name, xh->key_bmp); \
    if (xh->key_copy_f == NULL) \
      memcpy (new_key, key, sizeof(key_t)); \
    else \
      xh->key_copy_f (new_key, key); \
    ret = _xh_set_add2 (xh->hash, (void*)new_key, (void**)(&old_key)); \
    if (ret != XH_NEW) { \
      bmp_pop (xh_set_key##name, xh->key_bmp); \
      new_key = old_key; \
    } \
    return new_key; \
  } \
  \
  static inline int xh_set_search_##name ( \
      xh_set_##name##_t * xh, \
      key_t * key) { \
    return _xh_set_search (xh->hash, (void*)key); \
  } \
  static inline key_t * xh_set_search2_##name ( \
      xh_set_##name##_t * xh, \
      key_t * key) { \
    return (key_t*) _xh_set_search2 (xh->hash, (void*)key); \
  } \
  static inline xh_item_t * xh_set_search3_##name ( \
      xh_set_##name##_t * xh, \
      key_t * key) { \
    return _xh_set_search3 (xh->hash, (void*)key); \
  } \
  static inline int xh_set_key_iter_init_##name ( \
      xh_set_##name##_t * xh) { \
    return bmp_iter_init_xh_set_key##name (xh->key_bmp); \
  } \
  static inline key_t * xh_set_key_iter_next_##name ( \
      xh_set_##name##_t * xh) { \
    return bmp_iter_next_xh_set_key##name (xh->key_bmp); \
  } \
  static int __xh_set_##name##_def_end__ = 1

#define xh_set_t(name) xh_set_##name##_t
#define xh_set_cnt(xh) ((xh)->hash->cnt)

#define xh_set_init(name, size, load_factor, key_init_f, key_copy_f, key_hash_f, key_equal_f) \
  xh_set_init_##name (size, load_factor, key_init_f, key_copy_f, key_hash_f, key_equal_f)

#define xh_set_clear(name, xh, key_clear_f) \
  xh_set_clear_##name (xh, key_clear_f)

#define xh_set_free(name, xh, key_free_f) \
  xh_set_free_##name (xh, key_free_f)

#define xh_set_add(name, xh, key) \
  xh_set_add_##name (xh, key)

#define xh_set_add2(name, xh, key) \
  xh_set_add2_##name (xh, key)

#define xh_set_search(name, xh, key) \
  xh_set_search_##name (xh, key)

#define xh_set_search2(name, xh, key) \
  xh_set_search2_##name (xh, key)

#define xh_set_search3(name, xh, key) \
  xh_set_search3_##name (xh, key)

#define xh_set_key_iter_init(name, xh) \
  xh_set_key_iter_init_##name (xh)

#define xh_set_key_iter_next(name, xh) \
  xh_set_key_iter_next_##name (xh)

static inline uint64_t
itg_hash_f (const void * key)
{
  return (uint64_t) (*(int32_t*)key);
}

static inline int
itg_equal_f (const void * a, const void * b)
{
  int32_t va = * (int32_t *) a;
  int32_t vb = * (int32_t *) b;

  return va == vb;
}

HASH_SET_DEF (itg, int32_t);

#define xh_itg_set_init() xh_set_init(itg,128,0.75,NULL,NULL,itg_hash_f,itg_equal_f)

/**********************************************************
 ************************ Hash Maps ***********************
 **********************************************************/

#ifdef __cplusplus
extern "C" {
#endif

/*
 * if key is new added to set:
 *   return XH_NEW
 * if key already exist in set:
 *   return XH_EXIST
 * if key fail to add
 *   return XH_FAIL
 */
int _xh_map_add (xh_t * h, void * key, void * val);

/*
 * if key is in the set:
 *   return val
 * if key is not found in the set:
 *   return NULL
 */
void * _xh_map_search (xh_t * h, void * key);

#ifdef __cplusplus
}
#endif

#define HASH_MAP_DEF(name, key_t, val_t) \
  BMP_DEF (xh_map_key##name, key_t); \
  BMP_DEF (xh_map_val##name, val_t); \
  struct xh_map_##name##_s { \
    xh_t * hash; \
    bmp_t(xh_map_key##name) * key_bmp; \
    bmp_t(xh_map_val##name) * val_bmp; \
    void (*key_copy_f) (key_t*,key_t*); \
    void (*val_copy_f) (val_t*,val_t*); \
  }; \
  typedef struct xh_map_##name##_s xh_map_##name##_t; \
  \
  static inline xh_map_##name##_t * xh_map_init_##name ( \
      int64_t size, double load_factor, \
      void (*key_init_f)(key_t*), \
      void (*val_init_f)(val_t*), \
      void (*key_copy_f)(key_t*,key_t*), \
      void (*val_copy_f)(val_t*,val_t*), \
      HashFunc hash_f, \
      IsEqualFunc is_equal_f) { \
    xh_map_##name##_t * xh; \
    xh = (xh_map_##name##_t *) ckmalloc (sizeof(xh_map_##name##_t)); \
    xh->hash = _xh_init (size, load_factor, hash_f, is_equal_f); \
    xh->key_bmp = bmp_init (xh_map_key##name, size, key_init_f); \
    xh->val_bmp = bmp_init (xh_map_val##name, size, val_init_f); \
    xh->key_copy_f = key_copy_f; \
    xh->val_copy_f = val_copy_f; \
    return xh; \
  } \
  \
  static inline void xh_map_clear_##name ( \
      xh_map_##name##_t * xh, \
      void (*key_clear_f) (key_t*), \
      void (*val_clear_f) (val_t*)) { \
    bmp_clear (xh_map_key##name, xh->key_bmp, key_clear_f); \
    bmp_clear (xh_map_val##name, xh->val_bmp, val_clear_f); \
    _xh_clear (xh->hash); \
  } \
  \
  static inline void xh_map_free_##name ( \
      xh_map_##name##_t * xh, \
      void (*key_free_f) (key_t*), \
      void (*val_free_f) (val_t*)) { \
    bmp_free (xh_map_key##name, xh->key_bmp, key_free_f); \
    bmp_free (xh_map_val##name, xh->val_bmp, val_free_f); \
    _xh_free (xh->hash); \
    free (xh); \
  } \
  \
  static inline int xh_map_add_##name ( \
      xh_map_##name##_t * xh, \
      key_t * key, \
      val_t * val) { \
    int ret; \
    key_t * new_key; \
    val_t * new_val; \
    new_key = bmp_alloc (xh_map_key##name, xh->key_bmp); \
    if (xh->key_copy_f == NULL) \
      memcpy (new_key, key, sizeof(key_t)); \
    else \
      xh->key_copy_f (new_key, key); \
    new_val = bmp_alloc (xh_map_val##name, xh->val_bmp); \
    if (xh->val_copy_f == NULL) \
      memcpy (new_val, val, sizeof(val_t)); \
    else \
      xh->val_copy_f (new_val, val); \
    ret = _xh_map_add (xh->hash, (void*)new_key, (void*)new_val); \
    if (ret != XH_NEW) { \
      new_key = bmp_pop (xh_map_key##name, xh->key_bmp); \
      new_val = bmp_pop (xh_map_val##name, xh->val_bmp); \
    } \
    return ret; \
  } \
  \
  static inline val_t * xh_map_search_##name ( \
      xh_map_##name##_t * xh, \
      key_t * key) { \
    return (val_t*) _xh_map_search (xh->hash, key); \
  } \
  \
  static inline int xh_map_key_iter_init_##name ( \
      xh_map_##name##_t * xh) { \
    return bmp_iter_init_xh_map_key##name (xh->key_bmp); \
  } \
  \
  static inline key_t * xh_map_key_iter_next_##name ( \
      xh_map_##name##_t * xh) { \
    return bmp_iter_next_xh_map_key##name (xh->key_bmp); \
  } \
  \
  static int __xh_map_##name##_def_end__ = 1

#define xh_map_t(name) xh_map_##name##_t

#define xh_map_init(name, size, load_factor, key_init_f, val_init_f, key_copy_f, val_copy_f, key_hash_f, key_equal_f) \
  xh_map_init_##name (size, load_factor, key_init_f, val_init_f, key_copy_f, val_copy_f, key_hash_f, key_equal_f)

#define xh_map_clear(name, xh, key_clear_f, val_clear_f) \
  xh_map_clear_##name (xh, key_clear_f, val_clear_f)

#define xh_map_free(name, xh, key_free_f, val_free_f) \
  xh_map_free_##name (xh, key_free_f, val_free_f)

#define xh_map_add(name, xh, key, val) \
  xh_map_add_##name (xh, key, val)

#define xh_map_search(name, xh, key) \
  xh_map_search_##name (xh, key)

#define xh_map_key_iter_init(name,xh) \
  xh_map_key_iter_init_##name (xh)

#define xh_map_key_iter_next(name,xh) \
  xh_map_key_iter_next_##name (xh)

HASH_MAP_DEF (i2i, int32_t, int32_t);

#endif
