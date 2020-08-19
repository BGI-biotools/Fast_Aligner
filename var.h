/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2020-01-29 20:56:36
  *Edit History: 
***********************************************************/

#ifndef X_VAR_H
#define X_VAR_H

#include <sam.h>
#include <hts.h>

#include "bio.h"
#include "tmp.h"
#include "str.h"
#include "ebam.h"
#include "read.h"

#define DEBUG 0

#define BIN_SHIFT 23
#define BIN_SIZE  8388608

// shared variables
extern char * out_dir;
extern char * tmp_path;
extern char * spl_name;
extern int is_be;
extern int n_chr;
extern int n_libs;
extern int n_threads;
extern int need_filt;
extern int cal_mc;
extern int * n_chr_bins;
extern int64_t total_reads;
extern int64_t mapped_reads;
extern int64_t unmapped_reads;
extern genome_info_t * ginfo;
extern bam_hdr_t * bam_hdr;
extern bio_spl_t * spl;
extern ebam_file_t * unmapped_ebam_file;
extern ebam_file_t *** mapped_ebam_files; // lib/chr/bin

// only used in alignment step
extern xrd_filt_t * fastq_filter;
extern str_t * pth_str_buf;
extern mp_t(ar) * aln_res_set;
extern bmp_t(erd) ** raw_reads;
extern ebam_file_t ** batch_files;

// debug
extern FILE * dfp;

#endif
