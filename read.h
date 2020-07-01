/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-07-02 17:34:32
  *Edit History: 
***********************************************************/

#ifndef XDK_READ_H
#define XDK_READ_H

#include <sam.h>

#include "mp.h"
#include "bmp.h"
#include "aux.h"
#include "str.h"
#include "cigar.h"
#include "array.h"

/***********************************************************
 *                                                         *
 *                    Read Definition                      *
 *                                                         *
 ***********************************************************/

struct sam_core_s;
typedef struct sam_core_s sam_core_t;

struct sam_s;
typedef struct sam_s sam_t;

struct xread_s;
typedef struct xread_s xread_t;

struct sam_core_s {
  int16_t tid;
  int16_t mtid;
  int32_t pos;
  int32_t end;
  int32_t mpos;
  int32_t isize;
  int16_t qual;
  uint16_t flag;
};

struct sam_s {
  sam_core_t core;  
  cigar_t * cigar;
  aux_t * aux;
  char * iq;
  char * dq;
  char * gq;
};

struct xread_s {
	char * bases;
  char * quals;
	int32_t l, m;
  str_t * name;
  sam_t * aln;
};

MP_DEF (xrd, xread_t);

#ifdef __cplusplus
extern "C" {
#endif

xread_t * read_init (void);
void read_init2 (xread_t * r);
void read_resize (xread_t * r, int32_t new_size, int resize_all);
void read_clear (xread_t * r);
void read_free (xread_t * r);
void read_free2 (xread_t * r);
sam_t * read_aln_init (int init_all, int32_t rd_len);
void read_aln_copy (sam_t * dst, sam_t * src);
int read_aln_assign (sam_t * sc, bam1_t * b, int copy_aux);
void read_assign (xread_t * r, const char * name, const char * anno,
    const char * bases, const char * quals);

#ifdef __cplusplus
}
#endif

/***********************************************************
 *                                                         *
 *                 Read Common Functions                   *
 *                                                         *
 ***********************************************************/

#ifdef __cplusplus
extern "C" {
#endif

// if read quality system follows phred 33 system, return 1; else, return 0
int read_has_phred33_quals (xread_t * r);

// estimate read length from fastq file
// if first 100 reads are the same length, return the size; else, return -1
int read_length_estimate (const char * fastq_file);

int read_is_low_qual (char * quals, int l, int low_qual, int max_num_low_qual);
int read_has_too_many_Ns (char * bases, int l, int max_num_Ns);

/*
 * if the adatper can be aligned to the read, return the start of alignment(>=0); else, return -1
 */
int read_adapter_align (char * basess, int l, const char * adapter, int l_adapter, int min_num_match, int max_num_mismatch);

struct xh_set_adapter_s;
typedef struct xh_set_adapter_s adapter_set_t;

adapter_set_t * read_adapter_hash_init (void);
void read_adapter_hash_clear (adapter_set_t * set);
void read_adapter_hash_free (adapter_set_t * set);
int read_load_adapter_list (const char * adapter_list, adapter_set_t * adapter_hash);
int read_is_in_adapter_list (xread_t * r, adapter_set_t * adapter_hash);
int read_is_in_adapter_list2 (char * rd_name, int l_name, adapter_set_t * adapter_hash);

int read_fastq_dump (FILE * fp, xread_t * r, str_t * buf);

struct xrd_filt_s;
typedef struct xrd_filt_s xrd_filt_t;

struct xrd_filt_s {
  uint32_t flag;

  // low qual
  int32_t low_qual;
  int32_t max_num_low_qual[2];
  int32_t max_num_Ns[2];

  // adapter list
  adapter_set_t * ad_set[2];

  // adapter seq
  int32_t min_num_match[2];
  int32_t max_num_mismatch[2];
  str_t * ad_seq[2]; // backward adapter sequence
};

xrd_filt_t * xrd_filter_init (void);

int xrd_filter_set (xrd_filt_t * filter, uint32_t flag, int low_qual,
    int max_num_low_qual4read1, int max_num_low_qual4read2,
    int max_num_Ns4read1, int max_num_Ns4read2,
    int min_num_match4read1, int min_num_match4read2,
    int max_num_mismatch4read1, int max_num_mismatch4read2,
    const char * forw_ad_list, const char * back_ad_list,
    const char * forw_ad_seq, const char * back_ad_seq);

/*
 * if this pair/read is low quality, return 1; else, return 0
 */
int xrd_filter_low_qual (xrd_filt_t * filter,
    char * rd_name, int l_name,
    char * rd1_bases, char * rd1_quals, int l1,
    char * rd2_bases, char * rd2_quals, int l2);

void xrd_filter_free (xrd_filt_t * filter);

#ifdef __cplusplus
}
#endif

/***********************************************************
 *                                                         *
 *                  PE Read Definition                     *
 *                                                         *
 ***********************************************************/

#define PE_FLAG_LOW_QUAL 0x1
#define PE_FLAG_R1_SET   0x2
#define PE_FLAG_R2_SET   0x4

struct _pe_mate_s;
typedef struct _pe_mate_s _pe_mate_t;

struct pe_read_s;
typedef struct pe_read_s pe_read_t;

struct _pe_mate_s {
	char * bases;
	char * quals;
	int32_t l, m;
	sam_t * aln;
};

struct pe_read_s {
	_pe_mate_t r[2];
	str_t * name;
  uint64_t flag;
};

MP_DEF (xpe, pe_read_t);
BMP_DEF (xpe, pe_read_t);
HASH_SET_DEF (xpe, pe_read_t);
ARR_DEF (xpp, pe_read_t *);

static inline uint64_t
pe_read_hash_func (const void * key)
{
	uint32_t val;
	pe_read_t * pe = (pe_read_t *) key;
	val = blizzard_hash_func (pe->name->s, pe->name->l, 1);
	return val;
}

static inline int
pe_read_same_func (const void * key1, const void * key2)
{
	int val;
	pe_read_t * pe1 = (pe_read_t *) key1;
	pe_read_t * pe2 = (pe_read_t *) key2;
	val = str_equal (pe1->name, pe2->name);
	return val;
}

#ifdef __cplusplus
extern "C" {
#endif

pe_read_t * pe_read_init (void);
void pe_read_init2 (pe_read_t * pe);
void pe_read_resize (pe_read_t * pe, int32_t r1_new_size, int32_t r2_new_size);
void pe_read_clear (pe_read_t * pe);
void pe_read_free (pe_read_t * pe);
void pe_read_free2 (pe_read_t * pe);
void pe_read_aln_init (pe_read_t * pe);
void pe_read_copy (pe_read_t * dst, pe_read_t * src);
int pe_read_assign1mate (pe_read_t * pe, bam1_t * b, int mate_idx);

void pe_mate_copy (_pe_mate_t * dst, _pe_mate_t * src);
void pe_mate_resize (_pe_mate_t * mate, int32_t new_size);
void pe_mate_assign (_pe_mate_t * mate, bam1_t * b);

#ifdef __cplusplus
}
#endif

/***********************************************************
 *                                                         *
 *                        Read IO                          *
 *                                                         *
 ***********************************************************/

struct read_file_s;
typedef struct read_file_s read_file_t;

typedef hts_itr_t read_itr_t;

struct read_file_s {
	samFile * fp;
	bam_hdr_t * header;

	int32_t l_data;
	int32_t m_data;
	uint8_t * data;
};

#ifdef __cplusplus
extern "C" {
#endif

read_file_t * read_open (const char * file, const char * mode);

int read_close (read_file_t * fp);

int read_hdr_write (read_file_t * out, read_file_t * in);

int read_read (read_file_t * in, xread_t * r);

int read_dump (read_file_t * out, xread_t * r);

hts_itr_t * read_itr_query (const hts_idx_t * idx, int32_t tid, int32_t beg, int32_t end);

int read_itr_next (read_file_t * fp, read_itr_t * itr, xread_t * r);

#ifdef __cplusplus
}
#endif

#endif
