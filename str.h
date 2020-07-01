/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2017-05-16 16:38:00
  *Edit History: 
***********************************************************/

#ifndef XDK_STR_H
#define XDK_STR_H

#include <stdint.h>

#include "mp.h"
#include "hash.h"
#include "utils.h"
#include "hash_func.h"

struct str_s;
typedef struct str_s str_t;

struct str_s {
	char * s;
	int32_t l, m;
};

#ifdef __cplusplus
extern "C" {
#endif

	str_t * str_init (void);
	void str_init2 (str_t * str);

	void str_free (str_t * str);
	void str_free2 (str_t * str);

	void str_clear (str_t * str);

	int str_resize (str_t * str, int32_t l);

	int str_dump (FILE * fp, str_t * str);
	int str_write (FILE * fp, str_t * str);
	int str_read (FILE * fp, str_t * str);

	void str_copy (str_t * dst, str_t * src);

	str_t * str_dup (str_t * str);

	int str_assign (str_t * str, const char * s);
	int str_assign2 (str_t * str, const char * s, int32_t l);

	int str_append (str_t * str, const char * s, int32_t l);

  int str_add (str_t * str, char ch);
  void str_pop (str_t * str);

	// added on May 30, 2017
	int str_cmp (str_t * s1, str_t * s2);
  int str_cmp2 (const void * s1, const void * s2);
  int str_ncmp (str_t * s1, const char * s2, int32_t l);

	int str_equal (str_t * s1, str_t * s2);
  int str_equal2 (const void * s1, const void * s2);
  int str_nequal (str_t * s1, const char * s2, int32_t l);

  int str_is_empty (str_t * s);

  // update s->l according to strlen(s->s)
  // please be sure there is a '\0' !!!
  void str_update (str_t * s);

#ifdef __cplusplus
}
#endif

MP_DEF (xstr, str_t);
typedef mp_t(xstr) str_set_t;

#define str_set_init() mp_init(xstr, str_init2)
#define str_set_clear(set) mp_clear(xstr, (set), str_clear)
#define str_set_free(set) mp_free(xstr, set, str_free2)
#define str_set_alloc(set) mp_alloc(xstr, set)
#define str_set_resize(set,size) mp_resize(xstr, (set), (size))
#define str_set_at(set,idx) mp_at(xstr, (set), (idx))
#define str_set_cnt(set) mp_cnt(set)
#define str_set_std_sort(set) mp_std_sort(xstr, (set), str_cmp2)
#define str_set_uniqsort(set) mp_uniqsort(xstr, (set), str_cmp2, str_equal2)

HASH_SET_DEF (xstr, str_t);
HASH_MAP_DEF (s2i, str_t, int32_t);

static inline uint64_t
str_hash_func (const void * key)
{
  str_t * s = (str_t *) key;
  return blizzard_hash_func (s->s, s->l, 1);
}

#ifdef __cplusplus
extern "C" {
#endif

  void str_set_add (str_set_t * set, const char * seq);

	void str_set_copy (str_set_t * dst, str_set_t * src);

  str_set_t * load_file_list (const char * file_list);

#ifdef __cplusplus
}
#endif

#endif
