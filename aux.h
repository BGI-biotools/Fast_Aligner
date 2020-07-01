/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2017-05-30 12:48:42
  *Edit History: 
***********************************************************/

#ifndef XDK_AUX_H
#define XDK_AUX_H

#include <stdint.h>

#include "str.h"

#define AUX_VAL_TYPE_CNT 11

#define AUX_INT8 0
#define AUX_UINT8 1
#define AUX_INT16 2
#define AUX_UINT16 3
#define AUX_INT32 4
#define AUX_UINT32 5
#define AUX_INT64 6
#define AUX_UINT64 7
#define AUX_FLOAT 8
#define AUX_DOUBLE 9
#define AUX_STR 10

struct aux_item_s;
typedef struct aux_item_s aux_item_t;

struct aux_s;
typedef struct aux_s aux_t;

struct mp_aux_s;
typedef struct mp_aux_s aux_mp_t;

struct aux_item_s {
	char key[4];
	char tag;
	uint8_t type;
	uint16_t deleted;
	str_t * s;
};

struct aux_s {
	aux_mp_t * pool;
};

#ifdef __cplusplus
extern "C" {
#endif

	aux_t * aux_init (void);
	
	void aux_clear (aux_t * aux);
	
	void aux_free (aux_t * aux);
	
	int32_t aux_cnt (aux_t * aux);
	
	aux_item_t * aux_at (aux_t * aux, int32_t idx);
	
	aux_item_t * aux_find (aux_t * aux, char * key);
	
	int aux_add (aux_t * aux, char * key, uint8_t * val, uint8_t type, char tag);
	
	int aux_del (aux_t * aux, char * key);
	
	int aux_edit (aux_t * aux, char * key, uint8_t * val);
	
	int aux_read (FILE * fp, aux_t * aux);
	
	int aux_write (FILE * fp, aux_t * aux);
	
	int aux_copy (aux_t * dst, aux_t * src);
	
	int aux2hts_aux (aux_item_t * item, str_t * s);

#ifdef __cplusplus
}
#endif

#endif
