/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2017-08-29 14:52:59
  *Edit History: 
***********************************************************/

#ifndef XDK_TMP_H
#define XDK_TMP_H

struct tmp_dir_s;
typedef struct tmp_dir_s tmp_dir_t;

#ifdef __cplusplus
extern "C" {
#endif

	tmp_dir_t * tmp_dir_init (const char * dir);

	void tmp_dir_free (tmp_dir_t * tmp_dir);

	char * tmp_dir_add1file (tmp_dir_t * tmp_dir, const char * surfix);

#ifdef __cplusplus
}
#endif

#endif
