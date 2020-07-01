/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2017-05-22 11:22:52
  *Edit History: 
***********************************************************/

#ifndef XDK_BIO_H
#define XDK_BIO_H

#include <stdint.h>

#include "str.h"
#include "hash_func.h"

/**********************************************************
 ***************** Bio Macro Definitions ******************
 **********************************************************/

#define BASEQ_PHRED_OFFSET 33
#define MAX_READ_BASEQ_SUM 16383

#define int_comp(i) ((i)^0x02)
#define int2base(i) ("ACTG"[(i)])
#define base2int(b) ((b)>>1&0x03)

extern const char base_rc_tbl[128];
extern const char base2int_tbl[128];
extern const uint32_t cigar_bwa2samtools[5];

/**********************************************************
 ***************** Seq Address Functions ******************
 **********************************************************/

#ifdef __cplusplus
extern "C" {
#endif

str_t * reverse_complement1seq_i (str_t * seq, str_t * rv_seq);
str_t * reverse_complement1seq_b (str_t * seq, str_t * rv_seq);

#ifdef __cplusplus
}
#endif

/**********************************************************
 **************** Genome Address Functions ****************
 **********************************************************/

static inline uint64_t
gchr_key_hash_func (const void * key)
{
  str_t * s = (str_t *) key;
	return blizzard_hash_func (s->s, s->l, 1);
}

static inline int
gchr_key_equal_func (const void * key1, const void * key2)
{
  str_t * s1 = (str_t *) key1;
  str_t * s2 = (str_t *) key2;

  return str_equal (s1, s2);
}

struct genome_info_s;
typedef struct genome_info_s genome_info_t;

struct xh_map_gchr_s;
typedef struct xh_map_gchr_s chr_hash_t;

struct genome_info_s {
	str_t * target_name;
	int64_t * chr_idx_array;
	int64_t * target_len;
	int32_t n_targets;
	int64_t max_chr_len;
	int64_t genome_len;
	char * gseq;
	chr_hash_t * chr_hash;
};

#ifdef __cplusplus
extern "C" {
#endif

	genome_info_t * load_ref_genome_info (const char * ref);
	void free_genome_info (genome_info_t * ginfo);
	
	int check_ref_idx (const char * ref, const char * samtools);
	
	char * load_genome_seq (const char * ref, const genome_info_t * ginfo);
	char * get_ref_seq (genome_info_t * genone_info, int32_t tid, int64_t beg);
	int get_ref_seq2 (genome_info_t * genome_info, int32_t tid, int64_t beg, int64_t len, str_t * s);
	int get_ref_seq3 (genome_info_t * genome_info, int64_t gbeg, int64_t len, str_t * s);
	
	chr_hash_t * hash_ref_chr (const char * ref, const char * samtools);
	chr_hash_t * hash_info_chr (const genome_info_t * ginfo);
	int32_t get_chr_idx (chr_hash_t * chr_hash, str_t * chr_name);
	int32_t get_chr_idx2 (chr_hash_t * chr_hash, char * chr);
	int tandom_duplicate (char * ref_seq, int32_t left_limit, int32_t right_limit,
	    int32_t * dup_left_limit, int32_t * dup_right_limit);

	int32_t read_baseq_sum (char * qual, int32_t l_seq, int32_t offset);
	void read_name_parse (const char * name, int32_t l_name, int32_t * tile, int32_t * x, int32_t * y);

#ifdef __cplusplus
}
#endif

#endif
