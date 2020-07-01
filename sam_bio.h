/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-10-30 17:50:26
  *Edit History: 
***********************************************************/

#ifndef XDK_SAM_BIO_H
#define XDK_SAM_BIO_H

#include <stdint.h>

#include <sam.h>
#include <hts.h>

#include "mp.h"
#include "bio.h"
#include "hash.h"
#include "cigar.h"

#ifdef __cplusplus
extern "C" {
#endif

	genome_info_t * load_bam_genome_info (const char * bam);
	chr_hash_t * hash_bam_chr (const char * bam);
	int validate_ref_and_bam (const char * ref, const char * bam);
	
	samFile * xsam_open (const char * bam, const char * mode);
	bam_hdr_t * xsam_hdr_read (samFile * fp);
	hts_idx_t * xbam_index_load (const char * bam);
	
	int multi_mapping (bam1_t * b);
	int32_t mismatch_count (char * ref, bam1_t * b);
	int32_t idl_len_calc (bam1_t * b);
	int32_t sc_len_calc (bam1_t * b);
  int32_t bam_rpos2qpos (uint32_t * cigars, int32_t n_cigars,
      int32_t rd_beg, int32_t ref_pos, int32_t * cigar_idx, int32_t * cigar_pos);
  int bam_read_equ_ref (bam1_t * b);
  int bam_read_has_softclip (bam1_t * b);
  void bam_read_seq_dump (FILE * out, bam1_t * b);
  int bam_get_uncliped_begin (int32_t align_begin, int32_t n_cigar, uint32_t * cigars);
  int bam_get_uncliped_end (int32_t align_end, int32_t n_cigar, uint32_t * cigars);
  int bwa_get_uncliped_begin (int32_t align_begin, int32_t n_cigar, uint32_t * cigars);
  int bwa_get_uncliped_end (int32_t align_end, int32_t n_cigar, uint32_t * cigars);

	// TODO
	int proper_mapped_mates (bam1_t * b, int32_t max_isize);
	int bam_isize_stat (const char * bam, int32_t * cnter, int32_t n_pairs, int32_t max_isize);
	// TODO
	//int bam_isize_estimate (const char * bam, int32_t * mean, int32_t * sd);
	//int bam_isize_estimate2 (const char * bam, int32_t * min, int32_t * max);

#ifdef __cplusplus
}
#endif

struct dup_cls_s;
typedef struct dup_cls_s dup_cls_t;
struct dup_engine_s;
typedef struct dup_engine_s dup_engine_t;
struct duplex_s;
typedef struct duplex_s duplex_t;

struct dup_cls_s {
  uint64_t barcode;
  uint64_t seg_range; // beg<<32 | end
  int32_t n_rds[2]; // idx 0: support; idx 1: not support
  int32_t direct[2]; // idx 0: foorward support; idx 1: backward support
  int32_t sup[4];
  int32_t qpos;
  int32_t ref_start;
  str_t * dup_info;
  cigar_t * cigar;
};

HASH_SET_DEF (dup_cls, dup_cls_t);

struct dup_engine_s {
  xh_set_t(dup_cls) * clusters;
  char base;
  int8_t has_umi;
  int16_t lbc;
  int32_t var_pos;
  bam1_t * b;
  str_t * str_buf;
};

struct duplex_s {
  int32_t single_support;
  int32_t multi_ambiguous_support;
  int32_t multi_consistent_support;
  int32_t duplex_support;
};

static inline void
dup_cls_init2 (dup_cls_t * cls)
{
  cls->dup_info = str_init ();
  cls->cigar = cigar_init ();
}

static inline void
dup_cls_clear (dup_cls_t * cls)
{
  cls->sup[0] = cls->sup[1] = 0;
  cls->sup[2] = cls->sup[3] = 0;
  cls->n_rds[0] = cls->n_rds[1] = 0;
  cls->direct[0] = cls->direct[1] = 0;
  str_clear (cls->dup_info);
  cigar_clear (cls->cigar);
}

static inline void
dup_cls_free2 (dup_cls_t * cls)
{
  str_free (cls->dup_info);
  cigar_free (cls->cigar);
}

static inline void
dup_cls_copy (dup_cls_t * dst, dup_cls_t * src)
{
  dst->barcode = src->barcode;
  dst->seg_range = src->seg_range;
  cigar_copy (dst->cigar, src->cigar);
  dst->qpos = src->qpos;
  dst->ref_start = src->ref_start;
}

#ifdef __cplusplus
extern "C" {
#endif

	dup_engine_t * dup_engine_init (int has_umi, int lbc);
	void dup_engine_free (dup_engine_t * engine);
	void dup_engine_clear (dup_engine_t * engine);
	int dup_engine_add1candidate (dup_engine_t * engine, bam1_t * b, int32_t ref_pos);
	int duplex_calculate (samFile * in, hts_itr_t * iter,
	    int32_t tid, int32_t var_pos, char base, duplex_t * d, dup_engine_t * engine);
  int duplex_calculate_dump (const char * out_file, samFile * in, hts_itr_t * iter,
      int32_t tid, int32_t ref_pos, char ref_base, char alt_base,
      duplex_t * d, dup_engine_t * engine, genome_info_t * ginfo);
  bam_hdr_t * xsam_hdr_create (genome_info_t * ginfo);
  void xsam_hdr_add_line (bam_hdr_t * h, const char * mesg);
	void bam_serialize (str_t * s, const bam1_t * b, int is_be);

#ifdef __cplusplus
}
#endif

#endif
