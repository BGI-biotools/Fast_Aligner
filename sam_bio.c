/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-10-30 17:56:34
  *Edit History: 
***********************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sam.h>

#include "bio.h"
#include "bmp.h"
#include "hash.h"
#include "str.h"
#include "utils.h"
#include "sam_bio.h"
#include "sam_utils.h"
#include "hts_endian.h"

HASH_MAP_DEF (gchr, str_t, int32_t);

genome_info_t *
load_bam_genome_info (const char * bam)
{
	int32_t i;
	int32_t l;
	str_t * s;
	samFile * in;
	bam_hdr_t * head;
	genome_info_t * ginfo;

	ginfo = (genome_info_t *) ckalloc (1, sizeof(genome_info_t));
  ginfo->gseq = NULL;
  ginfo->chr_hash = NULL;

	in = xsam_open (bam, "rb");
	head = xsam_hdr_read (in);
	ginfo->n_targets = head->n_targets;

	ginfo->target_len = (int64_t *) ckalloc (ginfo->n_targets, sizeof(int64_t));
	ginfo->chr_idx_array = (int64_t *) ckalloc (ginfo->n_targets+1, sizeof(int64_t));
	ginfo->target_name = (str_t *) ckalloc (ginfo->n_targets, sizeof(str_t));

	ginfo->chr_idx_array[0] = 0;
	ginfo->max_chr_len = 0;
	for (i=0; i<ginfo->n_targets; i++) {
		str_init2 (ginfo->target_name+i);
		str_assign (ginfo->target_name+i, head->target_name[i]);

		ginfo->target_len[i] = head->target_len[i];
		if (ginfo->target_len[i] > ginfo->max_chr_len)
			ginfo->max_chr_len = ginfo->target_len[i];

		ginfo->chr_idx_array[i+1] = ginfo->chr_idx_array[i] + ginfo->target_len[i];
	}

	ginfo->genome_len = ginfo->chr_idx_array[ginfo->n_targets];
	ginfo->gseq = NULL;
	ginfo->chr_hash = NULL;

	sam_close (in);
	bam_hdr_destroy (head);

	return ginfo;
}

xh_map_t(gchr) *
hash_bam_chr (const char * bam)
{
	int32_t i;
	str_t * s;
	samFile * in;
	bam_hdr_t * head;
	xh_map_t(gchr) * set;

	s = str_init ();
	set = xh_map_init (gchr, 256, 0.75, str_init2, NULL, str_copy, NULL, gchr_key_hash_func, gchr_key_equal_func);
	in = xsam_open (bam, "rb");
	head = xsam_hdr_read (in);

	for (i=0; i<head->n_targets; i++) {
		str_assign (s, head->target_name[i]);
		xh_map_add (gchr, set, s, &i);
	}

	str_free (s);
	sam_close (in);
	bam_hdr_destroy (head);

	return set;
}

int
validate_ref_and_bam(const char * ref, const char * bam)
{
	int32_t i;
	str_t * s1;
	str_t * s2;
	genome_info_t * bam_ginfo;
	genome_info_t * ref_ginfo;

	bam_ginfo = load_bam_genome_info (bam);
	ref_ginfo = load_ref_genome_info (ref);

	if (bam_ginfo->n_targets != ref_ginfo->n_targets)
		return -1;
	if (bam_ginfo->max_chr_len != ref_ginfo->max_chr_len)
		return -1;
	if (bam_ginfo->genome_len != ref_ginfo->genome_len)
		return -1;

	if (memcmp(bam_ginfo->target_len, ref_ginfo->target_len,
				bam_ginfo->n_targets*sizeof(int64_t)) != 0)
		return -1;

	if (memcmp(bam_ginfo->chr_idx_array, ref_ginfo->chr_idx_array,
				(bam_ginfo->n_targets+1)*sizeof(int64_t)) != 0)
		return -1;

	for (i=0; i<bam_ginfo->n_targets; i++) {
		s1 = bam_ginfo->target_name + i;
		s2 = ref_ginfo->target_name + i;
		if (s1->l != s2->l)
			return -1;
		if (memcmp(s1->s, s2->s, s1->l) != 0)
			return -1;
	}

	return 0;
}

samFile *
xsam_open (const char * bam, const char * mode)
{
	samFile * fp;

	if ((fp = sam_open(bam,mode)) == NULL)
		err_mesg ("fail to open bam file '%s'!", bam);

	return fp;
}

bam_hdr_t *
xsam_hdr_read (samFile * fp)
{
	bam_hdr_t * hdr;

	if ((hdr = sam_hdr_read(fp)) == NULL)
		err_mesg ("fail to load bam header!");

	return hdr;
}

hts_idx_t *
xbam_index_load (const char * bam)
{
	hts_idx_t * bam_idx;

	if ((bam_idx = bam_index_load(bam)) == NULL)
		err_mesg ("fail to load index for bam '%s'!", bam);

	return bam_idx;
}

int
multi_mapping (bam1_t * b)
{
  uint8_t * aux;
  int as, xs;

  aux = bam_aux_get (b, "X0");
  if (aux!=NULL && bam_aux2i(aux)>1)
    return 1;
  aux = bam_aux_get (b, "XT");
  if (aux!=NULL && bam_aux2A(aux)=='R')
    return 1;
  aux = bam_aux_get (b, "X1");
  if (aux!=NULL && bam_aux2i(aux)>0)
    return 1;

  aux = bam_aux_get (b, "XA");
  if (aux != NULL)
    return 1;
  as = xs = -1;
  aux = bam_aux_get (b, "AS");
  if (aux != NULL)
    as = bam_aux2i (aux);
  aux = bam_aux_get (b, "XS");
  if (aux != NULL)
    xs = bam_aux2i (aux);
  if (as>=0 && xs>=0 && as<=xs)
    return 1;

  if (b->core.qual == 0)
    return 1;

  return 0;
}

int32_t
mismatch_count (char * ref, bam1_t * b)
{
  char rd_base;
  char ref_base;
  uint8_t * bs;
  int32_t i, j;
  int32_t opr;
  int32_t len;
  int32_t n_cigar;
  int32_t tot_len;
  int32_t n_mis;
  int32_t ref_os;
  int32_t rd_os;
  int32_t opr_type;
  uint32_t * cigar;

  n_cigar = b->core.n_cigar;
  cigar = bam_get_cigar (b);
  if (n_cigar<1 || cigar==NULL)
    return 0;

  bs = bam_get_seq (b);
  rd_os = ref_os = n_mis = 0;
  for (i=0; i<n_cigar; ++i) {
    opr = bam_cigar_op (cigar[i]);
    len = bam_cigar_oplen (cigar[i]);

    if (opr == BAM_CMATCH) {
      for (j=0; j<len; ++j) {
        ref_base = toupper (*(ref+ref_os+j));
        rd_base = seq_nt16_str[bam_seqi(bs,(rd_os+j))];
        if (ref_base != rd_base)
          ++n_mis;
      }
    }

    opr_type = bam_cigar_type (opr);
    if (opr_type & 1)
      rd_os += len;
    if (opr_type & 2)
      ref_os += len;
  }

  return n_mis;
}

int32_t
idl_len_calc (bam1_t * b)
{
  int32_t i;
  int32_t opr;
  int32_t len;
  int32_t rd_os;
  int32_t n_cigar;
  int32_t tot_len;
  uint32_t * cigar;

  n_cigar = b->core.n_cigar;
  cigar = bam_get_cigar (b);
  if (n_cigar<=1 || cigar==NULL)
    return 0;

  tot_len = rd_os = 0;
  for (i=0; i<n_cigar; ++i) {
    opr = bam_cigar_op (cigar[i]);
    len = bam_cigar_oplen (cigar[i]);
    if (opr==BAM_CINS || opr==BAM_CDEL)
      tot_len += len;
  }

  return tot_len;
}

int32_t
sc_len_calc (bam1_t * b)
{
  int32_t opr;
  int32_t len;
  int32_t n_cigar;
  int32_t tot_len;
  uint32_t * cigar;

  n_cigar = b->core.n_cigar;
  cigar = bam_get_cigar (b);
  if (n_cigar<=1 || cigar==NULL)
    return 0;

  tot_len = 0;
  opr = bam_cigar_op (cigar[0]);
  len = bam_cigar_oplen (cigar[0]);
  if (opr == BAM_CSOFT_CLIP)
    tot_len += len;

  opr = bam_cigar_op (cigar[n_cigar-1]);
  len = bam_cigar_oplen (cigar[n_cigar-1]);
  if (opr == BAM_CSOFT_CLIP)
    tot_len += len;

  return tot_len;
}

int32_t
bam_rpos2qpos (uint32_t * cigars, int32_t n_cigars,
    int32_t rd_beg, int32_t point, int32_t * cigar_idx, int32_t * cigar_pos)
{
  int32_t i;
  int32_t opr;
  int32_t len;
  int32_t qpos;
  int32_t rpos;
  int32_t type;

  for (i=qpos=0,rpos=rd_beg; i<n_cigars; ++i) {
    opr = bam_cigar_op (cigars[i]);
    len = bam_cigar_oplen (cigars[i]);

    if (opr == BAM_CMATCH) {
      if (rpos+len <= point) {
        qpos += len;
        rpos += len;
        continue;
      }
      if (rpos > point)
        break;
      *cigar_idx = i;
      *cigar_pos = point - rpos;
      return qpos + (point-rpos);
    }

    type = bam_cigar_type (opr);
    if (type & 1)
      qpos += len;
    if (type & 2)
      rpos += len;
  }

  return -1;
}

int
bam_read_equ_ref (bam1_t * b)
{
  uint8_t * aux;
  int nm;
  int32_t oprlen;
  uint32_t * cigars;
  bam1_core_t * c;

  c = &b->core;
  if (c->flag & 4)
    return 0;
  if (c->n_cigar != 1)
    return 0;
  cigars = bam_get_cigar (b);
  oprlen = bam_cigar_oplen (cigars[0]);
  if (oprlen != c->l_qseq)
    return 0;

  aux = bam_aux_get (b, "NM");
  if (aux!=NULL && bam_aux2i(aux) > 0)
    return 0;

  return 1;
}

int
bam_read_has_softclip (bam1_t * b)
{
  int32_t opr;
  int32_t len;
  int32_t n_cigar;
  uint32_t * cigars;

  if ((n_cigar = b->core.n_cigar) <= 1)
    return 0;

  cigars = bam_get_cigar (b);
  opr = bam_cigar_op (cigars[0]);
  if (opr == BAM_CSOFT_CLIP)
    return 1;
  opr = bam_cigar_op (cigars[n_cigar-1]);
  if (opr == BAM_CSOFT_CLIP)
    return 1;

  return 0;
}

void
bam_read_seq_dump (FILE * out, bam1_t * b)
{
  uint8_t * s;
  int32_t i;

  s = bam_get_seq (b);
  for (i=0; i<b->core.l_qseq; ++i)
    fprintf (out, "%c", seq_nt16_str[bam_seqi(s,i)]);
}

int
bam_get_uncliped_begin (int32_t align_begin, int32_t n_cigar, uint32_t * cigars)
{
  int i;
  int opr;
  int len;

  for (i=0; i<n_cigar; ++i) {
    opr = bam_cigar_op (cigars[i]);
    if (opr==BAM_CSOFT_CLIP || opr==BAM_CHARD_CLIP) {
      len = bam_cigar_oplen (cigars[i]);
      align_begin -= len;
    } else
      break;
  }

  return align_begin;
}

int
bam_get_uncliped_end (int32_t align_end, int32_t n_cigar, uint32_t * cigars)
{
  int i;
  int opr;
  int len;

  for (i=n_cigar-1; i>=0; --i) {
    opr = bam_cigar_op (cigars[i]);
    if (opr==BAM_CSOFT_CLIP || opr==BAM_CHARD_CLIP) {
      len = bam_cigar_oplen (cigars[i]);
      align_end += len;
    } else
      break;
  }

  return align_end;
}

/* BWA Cigar MIDSH */
int
bwa_get_uncliped_begin (int32_t align_begin, int32_t n_cigar, uint32_t * cigars)
{
  int i;
  int opr;
  int len;

  for (i=0; i<n_cigar; ++i) {
    opr = bam_cigar_op (cigars[i]);
    if (opr==3 || opr==4) {
      len = bam_cigar_oplen (cigars[i]);
      align_begin -= len;
    } else
      break;
  }

  return align_begin;
}

int
bwa_get_uncliped_end (int32_t align_end, int32_t n_cigar, uint32_t * cigars)
{
  int i;
  int opr;
  int len;

  for (i=n_cigar-1; i>=0; --i) {
    opr = bam_cigar_op (cigars[i]);
    if (opr==3 || opr==4) {
      len = bam_cigar_oplen (cigars[i]);
      align_end += len;
    } else
      break;
  }

  return align_end;
}

// TODO
int
proper_mapped_mates (bam1_t * b, int32_t max_isize)
{

}

int
bam_isize_stat (const char * bam, int32_t * cnter, int32_t n_pairs, int32_t max_isize)
{
	int32_t n;
	bam1_t * b;
	bam1_core_t * c;
	samFile * in;
	bam_hdr_t * header;

	memset (cnter, 0, n_pairs*sizeof(int32_t));

	n = 0;
	b = bam_init1 ();
	c = &b->core;
	in = xsam_open (bam, "rb");
	header = xsam_hdr_read (in);
	while (sam_read1(in,header,b) >= 0) {
		if (!proper_mapped_mates(b, max_isize))
			continue;
		if (c->isize < 0)
			continue;
		assert (c->isize>0 && c->isize<max_isize);
		++(cnter[c->isize]);
		if (++n > n_pairs)
			break;
	}
	bam_destroy1 (b);
	sam_close (in);
	bam_hdr_destroy (header);

	return 0;
}

#define N_ISIZE_STAT 10000
#define MAX_ISIZE    1000

// TODO
int
bam_isize_estimate (const char * bam, int32_t * mean, int32_t * sd)
{
	int32_t * isize_cnter;

	isize_cnter = (int32_t *) ckmalloc (MAX_ISIZE * sizeof(int32_t));
	bam_isize_stat (bam, isize_cnter, N_ISIZE_STAT, MAX_ISIZE);
}

// TODO
int
bam_isize_estimate2 (const char * bam, int32_t * min, int32_t * sd)
{
	int32_t * isize_cnter;

	isize_cnter = (int32_t *) ckmalloc (MAX_ISIZE * sizeof(int32_t));
	bam_isize_stat (bam, isize_cnter, N_ISIZE_STAT, MAX_ISIZE);
}

char
get_read_base (bam1_t * b, int32_t pos, int32_t * qpos)
{
  int32_t opr;
  int32_t len;
  int32_t type;
  int32_t i, j;
  int32_t rd_pos;
  int32_t ref_pos;
  uint32_t * cigar;
  bam1_core_t * c;

  c = &b->core;
  if (c->n_cigar <= 0)
    return -1;
  if ((cigar=bam_get_cigar(b)) == NULL)
    return -1;

  rd_pos = 0;
  ref_pos = c->pos;
  for (i=0; i<c->n_cigar; ++i) {
    opr = bam_cigar_op (cigar[i]);
    len = bam_cigar_oplen (cigar[i]);
    if (opr == BAM_CMATCH) {
      if (ref_pos+len <= pos) {
        rd_pos += len;
        ref_pos += len;
        continue;
      }
      if (ref_pos > pos)
        break;
      *qpos = rd_pos + (pos-ref_pos);
      return seq_nt16_str[bam_seqi(bam_get_seq(b),(*qpos))];
    }

    type = bam_cigar_type (opr);
    if (type & 1)
      rd_pos += len;
    if (type & 2)
      ref_pos += len;
  }

  return -1;
}

// duplex calculation

static inline uint64_t
dup_cls_hash_f (const void * key)
{
  dup_cls_t * c = (dup_cls_t *) key;
  return c->seg_range;
}

static inline int
dup_cls_equal_f (const void * a, const void * b)
{
  dup_cls_t * ca = (dup_cls_t *) a;
  dup_cls_t * cb = (dup_cls_t *) b;

  return ca->seg_range == cb->seg_range;
}

static inline uint64_t
umi_dup_cls_hash_f (const void * key)
{
  dup_cls_t * c = (dup_cls_t *) key;
  return c->barcode ^ c->seg_range;
}

static inline int
umi_dup_cls_equal_f (const void * a, const void * b)
{
  dup_cls_t * ca = (dup_cls_t *) a;
  dup_cls_t * cb = (dup_cls_t *) b;

  if (ca->barcode == cb->barcode
      && ca->seg_range == cb->seg_range)
    return 1;
  else
    return 0;
}

dup_engine_t *
dup_engine_init (int has_umi, int l_bc)
{
  dup_engine_t * engine;

  engine = (dup_engine_t *) ckmalloc (sizeof(dup_engine_t));
  engine->b = bam_init1 ();
  engine->str_buf = str_init ();
  if (has_umi) {
    engine->has_umi = 1;
    engine->lbc = l_bc;
    engine->clusters = xh_set_init (dup_cls, 4096, 0.75,
        dup_cls_init2, dup_cls_copy, umi_dup_cls_hash_f, umi_dup_cls_equal_f);
  } else {
    engine->has_umi = 0;
    engine->lbc = -1;
    engine->clusters = xh_set_init (dup_cls, 4096, 0.75,
        dup_cls_init2, dup_cls_copy, dup_cls_hash_f, dup_cls_equal_f);
  }

  return engine;
}

void
dup_engine_free (dup_engine_t * engine)
{
  xh_set_free (dup_cls, engine->clusters, dup_cls_free2);
  bam_destroy1 (engine->b);
  str_free (engine->str_buf);
  free (engine);
}

void
dup_engine_clear (dup_engine_t * engine)
{
  xh_set_clear (dup_cls, engine->clusters, dup_cls_clear);
}

static inline int64_t
umi_barcode_calculate (uint8_t * str, int len)
{
  int i;
  int64_t bc;

  for (i=bc=0; i<len; ++i) {
    if (str[i]=='N' || str[i]=='n')
      return -1;
    bc = (bc<<2) | base2int(str[i]);
  }

  return bc;
}

#define MAX_ISIZE4DUPLEX 1000
static int
read_valid4dup (bam1_t * b)
{
  char * mc;
  uint8_t * aux;
  uint32_t * cigar;
  bam1_core_t * c;

  c = &b->core;
  if ((c->flag & 4)
      || (c->flag & 8)
      || (c->flag & 256)
      || (c->flag & 2048))
    return 0;

  if (c->tid != c->mtid
      || c->isize == 0
      || abs(c->isize) > MAX_ISIZE4DUPLEX)
    return 0;

  if (c->n_cigar <= 0)
    return 0;
  if ((cigar=bam_get_cigar(b)) == NULL)
    return 0;

  if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP
      || bam_cigar_op(cigar[c->n_cigar-1]) == BAM_CSOFT_CLIP)
    return 0;

  if ((aux = bam_aux_get(b,"MC")) != NULL
      && (mc = bam_aux2Z(aux)) != NULL
      && strchr(mc,'S') != NULL)
    return 0;

  return 1;
}

#define DUP_FORW 0
#define DUP_BACK 1

static inline int
calculate_barcode_and_range (bam1_t * b, dup_cls_t * cls, int has_umi, int lbc, int * direct)
{
  int l;
  int l_data;
  int64_t left_bc;
  int64_t right_bc;
  uint64_t seg_beg;
  uint64_t seg_end;
  bam1_core_t * c;

  if (has_umi) {
    l = lbc / 2;
    l_data = strlen ((char*)b->data);
    assert (l_data > lbc);
    left_bc = umi_barcode_calculate (b->data+l_data-lbc, l);
    right_bc = umi_barcode_calculate (b->data+l_data-l, l);
    if (left_bc<0 || right_bc<0)
      return -1;
    if (left_bc < right_bc) {
      cls->barcode = ((uint64_t)left_bc<<32) | (uint64_t)right_bc;
      *direct = DUP_FORW;
    } else {
      cls->barcode = ((uint64_t)right_bc<<32) | (uint64_t)left_bc;
      *direct = DUP_BACK;
    }
  }

  c = &b->core;
  if (c->isize > 0) {
    seg_beg = c->pos;
    seg_end = seg_beg + c->isize;
    cls->seg_range = seg_beg<<32 | seg_end;
  } else {
    seg_end = c->pos + bam_cigar2qlen(c->n_cigar,bam_get_cigar(b));
    seg_beg = seg_end + c->isize;
    cls->seg_range = seg_end<<32 | seg_beg;
  }

  if (c->flag & 64)
    *direct = DUP_FORW;
  else
    *direct = DUP_BACK;

  return 0;
}

int
dup_engine_add1candidate (dup_engine_t * engine, bam1_t * b, int32_t ref_pos)
{
  char base;
  int d;
  int32_t qpos;
  uint64_t seg_beg;
  uint64_t seg_end;
  cigar_t cigar;
  dup_cls_t cls;

  if (!read_valid4dup(b))
    return -1;

  if ((base = get_read_base(b,ref_pos,&qpos)) < 0)
    return -1;

  cigar.cigar = bam_get_cigar (b);
  cigar.n = b->core.n_cigar;
  cls.cigar = &cigar;
  cls.qpos = qpos;
  cls.ref_start = b->core.pos;
  if (calculate_barcode_and_range(b, &cls, engine->has_umi, engine->lbc, &d) != 0)
    return -1;
  xh_set_add (dup_cls, engine->clusters, &cls);

  return 0;
}

static inline int
check_read_cigar (cigar_t * cigar, bam1_t * b)
{
  if (cigar->n != b->core.n_cigar)
    return -1;
  if (memcmp(cigar->cigar,bam_get_cigar(b),cigar->n*sizeof(uint32_t)) != 0)
    return -1;

  return 0;
}

static void
add1dup_info (str_t * s, bam1_t * b, int32_t qpos, char * ref_seq, str_t * str_buf)
{
  char ref_base;
  char qry_base;
  char * ch;
  char * buf;
  uint8_t * md;
  uint8_t * bs;
  uint8_t * qs;
  int32_t i, j;
  int32_t opr;
  int32_t len;
  int32_t type;
  int32_t ref_idx;
  int32_t qry_idx;
  uint32_t * cigar;

  ref_idx = qry_idx = 0;
  bs = bam_get_seq (b);
  qs = bam_get_qual (b);
  cigar = bam_get_cigar (b);
  for (i=0; i<b->core.n_cigar; ++i) {
    opr = bam_cigar_op (cigar[i]);
    len = bam_cigar_oplen (cigar[i]);

    if (opr == BAM_CMATCH) {
      for (j=0; j<len; ++j) {
        ref_base = toupper (ref_seq[ref_idx+j]);
        qry_base = seq_nt16_str[bam_seqi(bs,(qry_idx+j))];
        if (ref_base == qry_base)
          str_add (s, '.');
        else
          str_add (s, qry_base);
      }
    } else if (opr == BAM_CINS) {
      for (j=0; j<len; ++j) {
        qry_base = seq_nt16_str[bam_seqi(bs,(qry_idx+j))];
        str_add (s, qry_base);
      }
    } else if (opr == BAM_CDEL) {
      for (j=0; j<len; ++j)
        str_add (s, '*');
    }

    type = bam_cigar_type (opr);
    if (type & 1)
      qry_idx += len;
    if (type & 2)
      ref_idx += len;
  }

  buf = str_buf->s;
  *buf = '\0';

  // segment direction
  ch = buf;
  *ch++ = '\t';
  if (b->core.flag & 64)
    *ch++ = 'F';
  else
    *ch++ = 'R';

  // others
  md = bam_aux_get (b, "MD");
  sprintf (ch, "\t%d\t%d\t%d\t%d\t%s\t%s\n",
      qs[qpos], b->core.flag, b->core.qual, b->core.isize,
      md==NULL?"*":bam_aux2Z(md), b->data);
  str_append (s, buf, strlen(buf));
}

static int
add1read2cluster (dup_engine_t * engine, bam1_t * b,
    int32_t tid, int32_t ref_pos, char alt_base, int dump_info, genome_info_t * ginfo)
{
  char rd_base;
  int direct;
  int32_t qpos;
  dup_cls_t item;
  dup_cls_t * cls;

  if (!read_valid4dup(b))
    return -1;
  if (calculate_barcode_and_range(b,&item,engine->has_umi,engine->lbc,&direct) != 0)
    return -1;
  if ((cls = xh_set_search2(dup_cls,engine->clusters,&item)) == NULL)
    return -1;
  if ((rd_base=get_read_base(b,ref_pos,&qpos)) < 0)
    return -1;
  if (qpos!=cls->qpos || b->core.pos!=cls->ref_start)
    return -1;
  if (check_read_cigar(cls->cigar,b) != 0)
    return -1;
  if (rd_base == alt_base) {
    ++(cls->n_rds[0]);
    ++(cls->direct[direct]);
  }else
    ++(cls->n_rds[1]);
  ++(cls->sup[base2int(rd_base)]);

  if (dump_info)
    add1dup_info (cls->dup_info, b, qpos, ginfo->gseq+ginfo->chr_idx_array[tid]+b->core.pos, engine->str_buf);

  return 0;
}

static void
duplex_stat (dup_engine_t * engine, duplex_t * d)
{
  dup_cls_t * cls;

  memset (d, 0, sizeof(duplex_t));
  xh_set_key_iter_init (dup_cls, engine->clusters);
  while ((cls = xh_set_key_iter_next(dup_cls,engine->clusters)) != NULL) {
    assert (cls->n_rds[0]>=0 && cls->n_rds[1]>=0);
    assert (cls->direct[0]>=0 && cls->direct[1]>=0);
    if (cls->n_rds[0] == 0)
      continue;

    // duplex support
    if (cls->direct[0]>0 && cls->direct[1]>0)
      ++(d->duplex_support);
    // multi consistent support
    else if (cls->n_rds[0]>1 && cls->n_rds[1]==0)
      ++(d->multi_consistent_support);
    // multi ambiguous support
    else if (cls->n_rds[0]>1 && cls->n_rds[0]>cls->n_rds[1])
      ++(d->multi_ambiguous_support);
    // single support
    else {
      if (cls->n_rds[0]==1 && cls->n_rds[1]==0)
        ++(d->single_support);
    }
  }
}

int
duplex_calculate (samFile * in, hts_itr_t * iter,
    int32_t tid, int32_t ref_pos, char base, duplex_t * d, dup_engine_t * engine)
{
  bam1_t * b;

  b = engine->b;
  while (sam_itr_next(in,iter,b) >= 0)
    add1read2cluster (engine, b, tid, ref_pos, base, 0, NULL);

  duplex_stat (engine, d);

  return 0;
}

static void
dump_ref_seq (FILE * out, cigar_t * cigar, int32_t qpos, char * ref_seq)
{
  int is_started;
  int32_t i, j;
  int32_t opr;
  int32_t len;
  int32_t type;
  int32_t offset;
  int32_t ref_idx;
  int32_t qry_idx;

  qry_idx = offset = is_started = 0;
  for (i=0; i<cigar->n; ++i) {
    opr = bam_cigar_op (cigar->cigar[i]);
    len = bam_cigar_oplen (cigar->cigar[i]);

    if (opr == BAM_CMATCH) {
      if (qry_idx+len > qpos) {
        assert (qry_idx <= qpos);
        offset += qpos - qry_idx;
        break;
      }
    }

    type = bam_cigar_type (opr);
    if (type & 2)
      is_started = 1;
    if (is_started && type!=0)
      offset += len;
    if (type & 1)
      qry_idx += len;
  }
  for (i=0; i<offset; ++i)
    fprintf (out, " ");
  fprintf (out, "^\n");

  ref_idx = is_started = 0;
  for (i=0; i<cigar->n; ++i) {
    opr = bam_cigar_op (cigar->cigar[i]);
    len = bam_cigar_oplen (cigar->cigar[i]);
    type = bam_cigar_type (opr);
    if (type & 2)
      is_started = 1;
    if (type == 0)
      continue;
    if (type == 1) {
      for (j=0; j<len; ++j)
        fprintf (out, "*");
    } else {
      for (j=0; j<len; ++j)
        fprintf (out, "%c", ref_seq[ref_idx+j]);
    }
    if (type & 2)
      ref_idx += len;
  }
  fprintf (out, "\n");
}

int
duplex_calculate_dump (const char * out_file, samFile * in, hts_itr_t * iter,
    int32_t tid, int32_t ref_pos, char ref_base, char alt_base,
    duplex_t * d, dup_engine_t * engine, genome_info_t * ginfo)
{
  str_t * buf;
  bam1_t * b;
  dup_cls_t * cls;
  FILE * out;

  b = engine->b;
  str_resize (engine->str_buf, 4096);
  while (sam_itr_next(in,iter,b) >= 0)
    add1read2cluster (engine, b, tid, ref_pos, alt_base, 1, ginfo);

  duplex_stat (engine, d);

  out = ckopen (out_file, "w");
  fprintf (out, ">%s:%d:%c->%c\t%d:%d:%d:%d\n",
      ginfo->target_name[tid].s, ref_pos+1, ref_base, alt_base,
      d->duplex_support, d->multi_consistent_support,
      d->multi_ambiguous_support, d->single_support);

  xh_set_key_iter_init (dup_cls, engine->clusters);
  while ((cls = xh_set_key_iter_next(dup_cls,engine->clusters)) != NULL) {
    fprintf (out, "\n");
    fprintf (out, "%s\t%d\t%c\t%c\t%d\t%d\t%d\t%d\n",
        ginfo->target_name[tid].s, ref_pos+1, ref_base, alt_base,
        cls->sup[0], cls->sup[1], cls->sup[2], cls->sup[3]);
    dump_ref_seq (out, cls->cigar, cls->qpos, ginfo->gseq+ginfo->chr_idx_array[tid]+cls->ref_start);
    fprintf (out, "%s", cls->dup_info->s);
  }
  fclose (out);

  return 0;
}

bam_hdr_t *
xsam_hdr_create (genome_info_t * ginfo)
{
  char * c;
  int i;
  str_t * s;
  bam_hdr_t * h;

  h = bam_hdr_init ();
  h->n_targets = ginfo->n_targets;
  h->target_len = (uint32_t *) ckmalloc (ginfo->n_targets * sizeof(uint32_t));
  h->target_name = (char **) ckmalloc (ginfo->n_targets * sizeof(char *));
  for (i=0; i<ginfo->n_targets; ++i) {
    h->target_len[i] = (uint32_t)(ginfo->target_len[i]);
    h->target_name[i] = strdup (ginfo->target_name[i].s);
  }

  s = str_init ();
  str_resize (s, LONG_LINE_MAX);
  c = s->s;
  sprintf (c, "@HD\tVN:1.6\tSO:coordinate\n");
  for (i=0; i<ginfo->n_targets; ++i) {
    c = s->s + strlen (s->s);
    sprintf (c, "@SQ\tSN:%s\tLN:%u\n", h->target_name[i], h->target_len[i]);
  }
  h->l_text = strlen (s->s);
  h->text = s->s;

  return h;
}

void
xsam_hdr_add_line (bam_hdr_t * h, const char * mesg)
{
  char * c;
  char * t;
  int32_t l_mesg;
  str_t * s;

  s = str_init ();
  l_mesg = strlen (mesg);
  str_resize (s, h->l_text+l_mesg);
  sprintf (s->s, "%s%s\n", h->text, mesg);

  t = h->text;
  h->text = s->s;
  h->l_text += l_mesg + 1;
  s->s = t;
  str_free (s);
}

static void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host)
{
    uint32_t *cigar = (uint32_t*)(data + c->l_qname);
    uint32_t i;
    for (i = 0; i < c->n_cigar; ++i) ed_swap_4p(&cigar[i]);
}

// only used in short CIGAR, n_cigar <= 0xffff
void
bam_serialize (str_t * s, const bam1_t * b, int is_be)
{
	char * ch;
	const bam1_core_t * c;
	int i, ok;
	uint32_t y, l;
	uint32_t block_len;
	uint32_t x[8];

	c = &b->core;

	assert (c->n_cigar <= 0xffff);

	block_len = b->l_data - c->l_extranul + 32;
	str_resize (s, block_len);
	ch = s->s, l = 0;

	x[0] = c->tid;
	x[1] = c->pos;
	x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | (c->l_qname - c->l_extranul);
	x[3] = (uint32_t)c->flag << 16 | (c->n_cigar & 0xffff);
	x[4] = c->l_qseq;
	x[5] = c->mtid;
	x[6] = c->mpos;
	x[7] = c->isize;

	if (is_be) {
		for (i=0; i<8; ++i)
			ed_swap_4p (x+i);
		swap_data (c, b->l_data, b->data, 1);

		y = block_len;
		ed_swap_4p (&y);
		memcpy (ch, &y, 4);
		l += 4; ch += 4;
	}

	memcpy (ch, x, 32);
	l += 32; ch += 32;

	y = c->l_qname - c->l_extranul;
	memcpy (ch, b->data, y);
	l += y; ch += y;

	y = b->l_data - c->l_qname;
	memcpy (ch, b->data + c->l_qname, y);
	l += y; ch += y;

	s->l = l;
}
