/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-07-02 17:41:40
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include <sam.h>
#include <bgzf.h>

#include "aux.h"
#include "hash.h"
#include "cigar.h"
#include "read.h"
#include "utils.h"
#include "xendian.h"

#define INIT_SIZE 128
#define INIT_SIZE_BT_WIDTH 7

/***********************************************************
 *                                                         *
 *                    Read Definition                      *
 *                                                         *
 ***********************************************************/

xread_t *
read_init (void)
{
	xread_t * r;

	r = (xread_t *) ckalloc (1, sizeof(xread_t));
	read_init2 (r);

	return r;
}

void
read_init2 (xread_t * r)
{
	r->l = 0;
	r->m = INIT_SIZE;
	r->bases = (char *) ckalloc (INIT_SIZE, 1);
	r->quals = (char *) ckalloc (INIT_SIZE, 1);
  r->name = str_init ();
  r->aln = NULL;
}

void
read_resize (xread_t * r, int32_t new_size, int resize_all)
{
	if (new_size < r->m)
		return;

	r->m = ((new_size>>INIT_SIZE_BT_WIDTH)+1) << INIT_SIZE_BT_WIDTH;
	free (r->bases);
	free (r->quals);
	r->bases = (char *) ckalloc (r->m, 1);
	r->quals = (char *) ckalloc (r->m, 1);

  if (resize_all) {
    assert (r->aln->iq != NULL);
    free (r->aln->iq);
    free (r->aln->dq);
    free (r->aln->gq);
	  r->aln->iq = (char *) ckalloc (r->m, 1);
	  r->aln->dq = (char *) ckalloc (r->m, 1);
	  r->aln->gq = (char *) ckalloc (r->m, 1);
  }
}

void
read_clear (xread_t * r)
{
	r->l = 0;
  str_clear (r->name);
}

void
read_free (xread_t * r)
{
	read_free2 (r);
	free (r);
}

void
read_free2 (xread_t * r)
{
	free (r->bases);
	free (r->quals);
}

sam_t *
read_aln_init (int init_all, int32_t rd_len)
{
  sam_t * aln;

  aln = (sam_t *) ckalloc (1, sizeof(sam_t));
  aln->cigar = cigar_init ();
  aln->aux = aux_init ();
  if (init_all) {
	  aln->iq = (char *) ckalloc (rd_len, 1);
	  aln->dq = (char *) ckalloc (rd_len, 1);
	  aln->gq = (char *) ckalloc (rd_len, 1);
  } else
    aln->iq = aln->dq = aln->gq = NULL;

	return aln;
}

void
read_aln_copy (sam_t * dst, sam_t * src)
{
	memcpy (&dst->core, &src->core, sizeof(sam_core_t));
	cigar_copy (dst->cigar, src->cigar);
	aux_copy (dst->aux, src->aux);
}

int
read_aln_assign (sam_t * aln, bam1_t * b, int copy_aux)
{
	uint32_t * cigar;
	sam_core_t * sc;
	bam1_core_t * c;

	c = &b->core;
	sc = &aln->core;
	sc->flag  = c->flag;
	if (c->n_cigar <= 0) {
		sc->tid = -1;
		return -1;
	}

	sc->tid   = c->tid;
	sc->mtid  = c->mtid;
	sc->pos   = c->pos;
	sc->mpos  = c->mpos;
	sc->isize = c->isize;
	sc->qual  = c->qual;

	cigar = bam_get_cigar (b);
	cigar_assign (aln->cigar, cigar, c->n_cigar);
	sc->end = sc->pos + bam_cigar2rlen (c->n_cigar, cigar);

	return 0;
}

void
read_assign (xread_t * r, const char * name, const char * anno,
                const char * bases, const char * quals)
{
  int32_t l_seq;

  if (name != NULL)
   str_assign (r->name, name);

  l_seq = strlen (bases);
  if (l_seq != strlen(quals))
    err_mesg ("bases and quals are not same length!");

  read_resize (r, l_seq, 0);
  memcpy (r->bases, bases, l_seq);
  memcpy (r->quals, quals, l_seq);
  r->l = l_seq;
}

/***********************************************************
 *                                                         *
 *                  Fastq IO Definition                    *
 *                                                         *
 ***********************************************************/

static int
fq_buffer_load (fq_io_t * fq, int is_first_batch)
{
	char * ch;
	char * line;
	int32_t l;
	xread_t * s;

	// end of file
	if (!is_first_batch && fq->n<fq->m)
		return 1;

	line = fq->buf;
	fq->idx = fq->n = 0;
	while (gzgets(fq->fp,line,LINE_MAX)) {
		s = fq->seqs + fq->n++;

		// read name
		if ((ch=strrchr(line,'/')) != NULL) {
			*ch = '\0';
			str_assign (s->name, line);
		} else
			str_assign3 (s->name, line, '\n');

		// read bases
		assert (gzgets(fq->fp,line,LINE_MAX) != NULL);
		l = strlen (line);
		if (line[l-1] == '\n')
			--l;
		read_resize (s, l, 0);
		s->l = l;
		memcpy (s->bases, line, l);
		s->bases[l] = '\0';

		// read flag
		assert (gzgets(fq->fp,line,LINE_MAX) != NULL);

		// read quals
		assert (gzgets(fq->fp,line,LINE_MAX) != NULL);
		l = strlen (line);
		if (line[l-1] == '\n')
			--l;
		assert (l == s->l);
		memcpy (s->quals, line, l);
		s->quals[l] = '\0';

		if (fq->n >= fq->m)
			break;
	}

	if (fq->n == 0) // end of file
		return 1;
	else // fq->n > 0
		return 0;
}

fq_io_t *
fq_open (const char * path, const char * mode)
{
	int32_t i;
	fq_io_t * fq;

	fq = (fq_io_t *) ckmalloc (sizeof(fq_io_t));

	fq->n = 0;
	fq->m = FQ_SEQ_BATCH;
	fq->seqs = (xread_t *) ckmalloc (FQ_SEQ_BATCH * sizeof(xread_t));
	for (i=0; i<FQ_SEQ_BATCH; ++i)
		read_init2 (fq->seqs+i);

	fq->flag = 0;
	fq->buf = ALLOC_LINE;
	if (strchr(mode,'r')) {
		fq->flag |= FQ_IN;
		fq->fp = ckgzopen (path, "r");
		if (fq_buffer_load(fq,1) != 0)
			err_mesg ("'%s' is empty", path);
	} else {
		fq->flag |= FQ_OUT;
		fq->fp = ckgzopen (path, "w");
	}

	return fq;
}

xread_t *
fq_next (fq_io_t * fq)
{
	assert (fq->flag & FQ_IN);

	if (fq->idx < fq->n)
		return fq->seqs + fq->idx++;

	if (fq_buffer_load(fq,0) != 0)
		return NULL;

	return fq->seqs + fq->idx++;
}

void
fq_close (fq_io_t * fq)
{
	int32_t i;

	for (i=0; i<fq->m; ++i)
		read_free2 (fq->seqs+i);
	free (fq->seqs);
	free (fq->buf);
	gzclose (fq->fp);
	free (fq);
}

// original qualities, before -33/-64
// if read quality system follows phred 33 system, return 1; else, return 0
int
fq_has_phred33_quals (fq_io_t * fq)
{
	int64_t i;

	assert (fq->n > 0);
	for (i=0; i<fq->n; ++i)
		if (read_has_phred33_quals(fq->seqs+i))
			return 1;

	return 0;
}

// original qualities, before -33/-64
// if read quality system follows phred 33 system, return 1; else, return 0
int
fq_has_phred33_quals2 (const char * fq_file)
{
	int ret;
	fq_io_t * fq;

	fq = fq_open (fq_file, "r");
	ret = fq_has_phred33_quals (fq);
	fq_close (fq);

	return ret;
}

/***********************************************************
 *                                                         *
 *                 Read Common Functions                   *
 *                                                         *
 ***********************************************************/

// original qualities, before -33/-64
// if read quality system follows phred 33 system, return 1; else, return 0
int
read_has_phred33_quals (xread_t * r)
{
  int i;

  for (i=0; i<r->l; ++i)
    if (r->quals[i] < 64)
      return 1;

  return 0;
}

int
read_length_estimate (const char * fastq_file)
{
  char line[LINE_MAX];
  int pre_len;
  int cur_len;
  int cnter;
  gzFile in;

  cnter = 0;
  pre_len = -1;
  in = ckgzopen (fastq_file, "r");
  while (gzgets(in, line, LINE_MAX)) {
    // read sequence
    if (gzgets(in, line, LINE_MAX) == NULL)
      err_mesg ("invalid fastq file!");
    chomp (line);
    cur_len = strlen (line);
    if (pre_len!=-1 && cur_len!=pre_len)
      return -1;
    pre_len = cur_len;

    if (++cnter >= 100)
      break;

    // flag line
    if (gzgets(in, line, LINE_MAX) == NULL)
      err_mesg ("invalid fastq file!");

    // qual sequence
    if (gzgets(in, line, LINE_MAX) == NULL)
      err_mesg ("invalid fastq file!");
  }
  gzclose (in);

  return cur_len;
}

// fixed qualities, after -33/-64
int
read_is_low_qual (char * quals, int l, int low_qual, double max_num_low_qual)
{
  int i;
  double num_low_qual;

	num_low_qual = 0;
  for (i=0; i<l; ++i) {
		if (quals[i] <= low_qual) {
			num_low_qual += 1;
			if (num_low_qual >= max_num_low_qual)
				return 1;
		}
	}

  return 0;
}

int
read_has_too_many_Ns (char * bases, int l, double max_num_Ns)
{
  int i;
  double num_Ns;

	num_Ns = 0;
  for (i=0; i<l; ++i) {
		if (bases[i] == 'N') {
			num_Ns += 1;
			if (num_Ns >= max_num_Ns)
				return 1;
		}
	}

  return 0;
}

int
read_adapter_align (char * bases, int l, const char * adapter, int l_adapter,
		int min_num_match, int max_num_mismatch, float mis_grad, int adapter_edge)
{
	int i, j;
	int mis;
	int num_match;
	int mis_tmp;
	int min_edge5 = l_adapter - 5;

	for (i=0; i<5; ++i) {
		mis = 0;
		num_match = 0;
		mis_tmp = i / mis_grad;
		for (j=0; j<i+min_edge5; ++j) {
			if (adapter[l_adapter-i-min_edge5+j] == bases[j]) {
				if (++num_match >= min_num_match)
					return 0;
			} else {
				num_match = 0;
				if (++mis > mis_tmp)
					break;
			}
		}
		if (mis <= mis_tmp)
			return 0;
	}

	for (i=0; i<=l-l_adapter; ++i) {
		mis = 0;
		num_match = 0;
		for (j=0; j<l_adapter; ++j) {
			if (adapter[j] == bases[i+j]) {
				if (++num_match >= min_num_match)
					return i;
			} else {
				num_match = 0;
				if (++mis > max_num_mismatch)
					break;
			}
		}
		if (mis <= max_num_mismatch)
			return i;
	}

	for (i=0; i<l_adapter-adapter_edge; ++i) {
		mis = 0;
		num_match = 0;
		mis_tmp = i / mis_grad;
		for (j=0; j<i+adapter_edge; ++j) {
			if (adapter[j] == bases[l-i-adapter_edge+j]) {
				if (++num_match >= min_num_match)
					return l - i - adapter_edge;
			} else {
				num_match = 0;
				if (++mis > mis_tmp)
					break;
			}
		}
		if (mis <= mis_tmp)
			return l - i - adapter_edge;
	}

	return -1;
}

/*
int
read_adapter_align (char * bases, int l, const char * adapter, int l_adapter, int min_num_match, int max_num_mismatch)
{
  char * rd_seq;
  int find;
  int ad_idx;
  int rd_idx;
  int len;
  int len1;
  int len2;
  int i;
  int mis;
  int right;
  int max_len_map;
  int map[LINE_MAX];

  find = -1;
  ad_idx = l_adapter - min_num_match;
  right = l - min_num_match;
  rd_seq = bases;

  for (rd_idx=0; rd_idx<=right; ) {
    len1 = l_adapter - ad_idx;
    len2 = l - rd_idx;
    len = (len1<len2) ? len1 : len2;

    mis = map[0] = 0;
    for (i=0; i<len; ++i) {
      if (adapter[ad_idx+i] == rd_seq[rd_idx+i])
        ++map[mis];
      else {
        ++mis;
        map[mis] = 0;
      }
    }

    max_len_map = 0;
    for (i=0; i<mis; ++i)
      if (map[i] > max_len_map)
        max_len_map = map[i];

    if (mis<=max_num_mismatch || max_len_map>=min_num_match) {
      find = rd_idx;
      break;
    }

    if (ad_idx > 0)
      --ad_idx;
    else
      ++rd_idx;
  }

  return find;
}
*/

HASH_SET_DEF (adapter, str_t);

adapter_set_t *
read_adapter_hash_init (void)
{
  return xh_set_init (adapter, 65536, 0.75, str_init2, str_copy, str_hash_func, str_equal2);
}

void
read_adapter_hash_clear (adapter_set_t * set)
{
  xh_set_clear (adapter, set, str_clear);
}

void
read_adapter_hash_free (adapter_set_t * set)
{
  xh_set_free (adapter, set, str_free2);
}

int
read_load_adapter_list (const char * adapter_list, adapter_set_t * adapter_hash)
{
  char * name;
  char * line;
  str_t * s;
  gzFile in;

  s = str_init ();
  name = ALLOC_LINE;
  line = ALLOC_LINE;
  if ((in = gzopen(adapter_list, "r")) == NULL)
    err_mesg ("fail to open file '%s' to read!", adapter_list);

  name[0] = '@';
  while (gzgets(in, line, LINE_MAX)) {
    if (*line == '#') // header
      continue;

    sscanf (line, "%s", name+1);
    str_cut_tail (name, "/1");
    str_cut_tail (name, "/2");

    str_assign (s, name);
    xh_set_add (adapter, adapter_hash, s);
  }

  gzclose (in);
  free (name);
  free (line);
  str_free (s);

  return 0;
}

int
read_is_in_adapter_list (xread_t * r, adapter_set_t * adapter_hash)
{
  if (xh_set_search(adapter, adapter_hash, r->name) == XH_EXIST)
    return 0;
  else
    return 1;
}

int
read_is_in_adapter_list2 (char * rd_name, int l_name, adapter_set_t * adapter_hash)
{
  str_t s;

  s.s = rd_name;
  s.l = l_name;
  if (xh_set_search(adapter, adapter_hash, &s) == XH_EXIST)
    return 0;
  else
    return 1;
}

int
read_fastq_dump (FILE * fp, xread_t * r, str_t * buf)
{
  char * ch;
  int32_t i, l;

  str_clear (buf);

  str_append (buf, r->name->s, r->name->l);
  str_append (buf, "\n", 1);

  str_append (buf, r->bases, r->l);
  str_append (buf, "\n", 1);

  str_append (buf, "+\n", 2);

  str_resize (buf, buf->l+r->l+1);
  ch = buf->s + buf->l;
  for (i=0; i<r->l; ++i,++ch)
    *ch = r->quals[i] + 33;
  *(ch++) = '\n';
  *ch = '\0';

  fprintf (fp, "%s", buf->s);

  return 0;
}

xrd_filt_t *
xrd_filter_init (void)
{
  xrd_filt_t * filter;

  filter = (xrd_filt_t *) ckmalloc (sizeof(xrd_filt_t));
  filter->ad_set[0] = read_adapter_hash_init ();
  filter->ad_set[1] = read_adapter_hash_init ();
  filter->ad_seq[0] = str_init ();
  filter->ad_seq[1] = str_init ();

	return filter;
}

int
xrd_filter_set (xrd_filt_t * filter, uint32_t flag, int low_qual,
    double max_num_low_qual4read1, double max_num_low_qual4read2,
    double max_num_Ns4read1, double max_num_Ns4read2,
    int min_num_match4read1, int min_num_match4read2,
    int max_num_mismatch, const char * forw_ad_list, const char * back_ad_list,
    const char * forw_ad_seq, const char * back_ad_seq)
{
  filter->flag = flag;
  filter->low_qual = low_qual;

  filter->max_num_low_qual[0] = max_num_low_qual4read1;
  filter->max_num_Ns[0] = max_num_Ns4read1;

  if (flag & BIO_FASTQ_AD_LIST) {
    read_adapter_hash_clear (filter->ad_set[0]);
    read_load_adapter_list (forw_ad_list, filter->ad_set[0]);
  } else if (flag & BIO_FASTQ_AD_STR) {
		filter->adapter_edge = 6;
  	filter->max_num_mismatch = max_num_mismatch;

    str_assign (filter->ad_seq[0], forw_ad_seq);
  	filter->min_num_match[0] = min_num_match4read1;
		filter->mis_grad[0] = (filter->ad_seq[0]->l - filter->adapter_edge) / (filter->max_num_mismatch+1);
  }

  if (!(flag & BIO_FASTQ_PE))
    return 0;

  filter->max_num_low_qual[1] = max_num_low_qual4read2;
  filter->max_num_Ns[1] = max_num_Ns4read2;

  if (flag & BIO_FASTQ_AD_LIST) {
    read_adapter_hash_clear (filter->ad_set[1]);
    read_load_adapter_list (back_ad_list, filter->ad_set[1]);
  } else if (flag & BIO_FASTQ_AD_STR) {
    str_assign (filter->ad_seq[1], back_ad_seq);
  	filter->min_num_match[1] = min_num_match4read2;
		filter->mis_grad[1] = (filter->ad_seq[1]->l - filter->adapter_edge) / (filter->max_num_mismatch+1);
  }

  return 0;
}

int
xrd_filter_low_qual (xrd_filt_t * filter,
    char * rd_name, int l_name,
    char * rd1_bases, char * rd1_quals, int l1,
    char * rd2_bases, char * rd2_quals, int l2)
{
  if (read_is_low_qual(rd1_quals,l1,filter->low_qual,filter->max_num_low_qual[0]))
    return 1;
  if (read_has_too_many_Ns(rd1_bases,l1,filter->max_num_Ns[0]))
    return 1;

  if (filter->flag & BIO_FASTQ_PE) {
    if (read_is_low_qual(rd2_quals,l2,filter->low_qual,filter->max_num_low_qual[1]))
      return 1;
    if (read_has_too_many_Ns(rd2_bases,l2,filter->max_num_Ns[1]))
      return 1;
  }

  if (filter->flag & BIO_FASTQ_AD_STR) {
    if (read_adapter_align(rd1_bases, l1, filter->ad_seq[0]->s, filter->ad_seq[0]->l,
          filter->min_num_match[0], filter->max_num_mismatch,
					filter->mis_grad[0], filter->adapter_edge) >= 0) {
      return 1;
		}
    if ((filter->flag & BIO_FASTQ_PE)
        && read_adapter_align(rd2_bases, l2, filter->ad_seq[1]->s, filter->ad_seq[1]->l,
          filter->min_num_match[1], filter->max_num_mismatch,
					filter->mis_grad[1], filter->adapter_edge) >= 0) {
      return 1;
		}
  } else if (filter->flag & BIO_FASTQ_AD_LIST)
    if (read_is_in_adapter_list2(rd_name,l_name,filter->ad_set[0]))
      return 1;

  return 0;
}

void
xrd_filter_free (xrd_filt_t * filter)
{
  read_adapter_hash_free (filter->ad_set[0]);
  read_adapter_hash_free (filter->ad_set[1]);
  str_free (filter->ad_seq[0]);
  str_free (filter->ad_seq[1]);
}
/***********************************************************
 *                                                         *
 *                  PE Read Definition                     *
 *                                                         *
 ***********************************************************/

pe_read_t *
pe_read_init (void)
{
  pe_read_t * pe;

  pe = (pe_read_t *) ckalloc (1, sizeof(pe_read_t));
	pe_read_init2 (pe);

  return pe;
}

void
pe_read_init2 (pe_read_t * pe)
{
	int i;
	_pe_mate_t * r;

	pe->name = str_init ();
	pe->flag = 0;
	for (i=0; i<2; ++i) {
		r = pe->r + i;
		r->l = 0;
		r->m = INIT_SIZE;
		r->bases = (char *) ckmalloc (INIT_SIZE);
		r->quals = (char *) ckmalloc (INIT_SIZE);
		r->aln = NULL;
	}
}

void
pe_read_resize (pe_read_t * pe, int32_t r1_new_size, int32_t r2_new_size)
{
	_pe_mate_t * r;

	// read 1
	r = pe->r;
	if (r1_new_size > r->m) {
		r->m = ((r1_new_size>>INIT_SIZE_BT_WIDTH)+1) << INIT_SIZE_BT_WIDTH;
		free (r->bases);
		free (r->quals);
		r->bases = (char *) ckmalloc (r->m);
		r->quals = (char *) ckmalloc (r->m);
	}

	// read 2
	r = pe->r + 1;
	if (r2_new_size > r->m) {
		r->m = ((r2_new_size>>INIT_SIZE_BT_WIDTH)+1) << INIT_SIZE_BT_WIDTH;
		free (r->bases);
		free (r->quals);
		r->bases = (char *) ckmalloc (r->m);
		r->quals = (char *) ckmalloc (r->m);
	}
}

void
pe_read_clear (pe_read_t * pe)
{
	pe->flag = 0;
	str_clear (pe->name);
	pe->r[0].l = 0;
	pe->r[1].l = 0;
}

void
pe_read_free (pe_read_t * pe)
{
	pe_read_free2 (pe);
  free (pe);
}

void
pe_read_free2 (pe_read_t * pe)
{
	str_free (pe->name);
	free (pe->r[0].bases);
	free (pe->r[0].quals);
	free (pe->r[1].bases);
	free (pe->r[1].quals);
}

void
pe_read_aln_init (pe_read_t * pe)
{
	pe->r[0].aln = read_aln_init (0, pe->r[0].m);
	pe->r[1].aln = read_aln_init (0, pe->r[1].m);
}

void
pe_read_copy (pe_read_t * dst, pe_read_t * src)
{
	// copy flag
	dst->flag = src->flag;

	// copy read name
	str_copy (dst->name, src->name);

	// copy pe mates
	pe_mate_copy (dst->r, src->r);
	pe_mate_copy (dst->r+1, src->r+1);
}

int
pe_read_assign1mate (pe_read_t * pe, bam1_t * b, int mate_idx)
{
	uint8_t * bs;
	uint8_t * qs;
	int32_t i;
	_pe_mate_t * r;
	bam1_core_t * c;

	c = &b->core;
	if (mate_idx == 0) {
		pe->flag |= PE_FLAG_R1_SET;
		pe_read_resize (pe, c->l_qseq, -1);
	} else if (mate_idx == 1) {
		pe->flag |= PE_FLAG_R2_SET;
		pe_read_resize (pe, -1, c->l_qseq);
	} else
		err_mesg ("invalid mate index!");

	pe_mate_assign (pe->r+mate_idx, b);

	return 0;
}

void
pe_mate_copy (_pe_mate_t * dst, _pe_mate_t * src)
{
	pe_mate_resize (dst, src->l);
	dst->l = src->l;
	memcpy (dst->bases,src->bases, src->l);
	memcpy (dst->quals,src->quals, src->l);

	if (src != NULL)
		read_aln_copy (dst->aln, src->aln);
	else
		dst->aln = NULL;
}

void
pe_mate_resize (_pe_mate_t * mate, int32_t new_size)
{
	if (new_size <= mate->m)
		return;

	mate->m = ((new_size>>INIT_SIZE_BT_WIDTH)+1) << INIT_SIZE_BT_WIDTH;
	free (mate->bases);
	free (mate->quals);
	mate->bases = (char *) ckmalloc (mate->m);
	mate->quals = (char *) ckmalloc (mate->m);
}

void
pe_mate_assign (_pe_mate_t * mate, bam1_t * b)
{
	uint8_t * bs;
	uint8_t * qs;
	int32_t i;
	bam1_core_t * c;

	c = &b->core;

	// copy base sequence
	bs = bam_get_seq (b);
	qs = bam_get_qual (b);
	for (i=0; i<c->l_qseq; ++i) {
		mate->bases[i] = seq_nt16_str[bam_seqi(bs,i)];
		mate->quals[i] = qs[i];
	}

	// assign alignment information
	if (mate->aln == NULL)
		mate->aln = read_aln_init (0, 0);
	read_aln_assign (mate->aln, b, 0);
}

/***********************************************************
 *                                                         *
 *                        Read IO                          *
 *                                                         *
 ***********************************************************/

typedef bam_hdr_t read_hdr_t;

static inline int
aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}

static void
swap_bam_data(uint8_t *data, int is_host,
		int32_t l_data, int32_t l_qname, int32_t n_cigar, int32_t l_qseq)
{
    uint8_t *s;
    uint32_t *cigar = (uint32_t*)(data + l_qname);
    uint32_t i, n;
    s = data + n_cigar*4 + l_qname + l_qseq + (l_qseq + 1)/2;
    for (i = 0; i < n_cigar; ++i) ed_swap_4p(&cigar[i]);
    while (s < data + l_data) {
        int size;
        s += 2; // skip key
        size = aux_type2size(*s); ++s; // skip type
        switch (size) {
        case 1: ++s; break;
        case 2: ed_swap_2p(s); s += 2; break;
        case 4: ed_swap_4p(s); s += 4; break;
        case 8: ed_swap_8p(s); s += 8; break;
        case 'Z':
        case 'H':
            while (*s) ++s;
            ++s;
            break;
        case 'B':
            size = aux_type2size(*s); ++s;
            if (is_host) memcpy(&n, s, 4), ed_swap_4p(s);
            else ed_swap_4p(s), memcpy(&n, s, 4);
            s += 4;
            switch (size) {
            case 1: s += n; break;
            case 2: for (i = 0; i < n; ++i, s += 2) ed_swap_2p(s); break;
            case 4: for (i = 0; i < n; ++i, s += 4) ed_swap_4p(s); break;
            case 8: for (i = 0; i < n; ++i, s += 8) ed_swap_8p(s); break;
            }
            break;
        }
    }
}

static void
unpack_data (xread_t * r, uint8_t * data, int32_t l_data,
		int32_t l_qname, int32_t n_cigar, int32_t l_qseq)
{
	char key[2];
	uint8_t * s;
	uint8_t * aux;
	uint8_t * bb;
	uint8_t * bq;
	int type;
	int32_t i;
	int32_t l;
  sam_t * aln;

	s = data;
  aln = r->aln;

	r->name->l = l_qname - 1;
	str_resize (r->name, r->name->l);
	str_assign (r->name, (const char *)data);
	s += l_qname;

	aln->cigar->n = n_cigar;
	cigar_resize (aln->cigar, n_cigar);
	memcpy (aln->cigar->cigar, s, n_cigar<<2);
	s += n_cigar << 2;

	read_resize (r, r->l, 0);
	bb = s; s += (r->l+1) >> 1;
	bq = s; s += r->l;
	for (i=0; i<r->l; ++i) {
		r->bases[i] = seq_nt16_str[bam_seqi(bb,i)];
		r->quals[i] = bq[i];
	}

	aux = s;
	aux_clear (aln->aux);
	while (aux+4 <= data+l_data) {
		key[0] = aux[0]; key[1] = aux[1];
		aux += 2; type = *aux++;
		if (type == 'A') {
			l = 1;
			aux_add (aln->aux, key, aux, AUX_UINT8, type);
		} else if (type == 'C') {
			l = 1;
			aux_add (aln->aux, key, aux, AUX_UINT8, type);
		} else if (type == 'c') {
			l = 1;
			aux_add (aln->aux, key, aux, AUX_INT8, type);
		} else if (type == 'S') {
			l = 2;
			aux_add (aln->aux, key, aux, AUX_UINT16, type);
		} else if (type == 's') {
			l = 2;
			aux_add (aln->aux, key, aux, AUX_INT16, type);
		} else if (type == 'I') {
			l = 4;
			aux_add (aln->aux, key, aux, AUX_UINT32, type);
		} else if (type == 'i') {
			l = 4;
			aux_add (aln->aux, key, aux, AUX_INT32, type);
		} else if (type == 'f') {
			l = 4;
			aux_add (aln->aux, key, aux, AUX_FLOAT, type);
		} else if (type == 'd') {
			l = 8;
			aux_add (aln->aux, key, aux, AUX_DOUBLE, type);
		} else if (type == 'Z' || type == 'H') {
			l = 0;
			while (*(aux+l))
				++l;
			aux_add (aln->aux, key, aux, AUX_STR, type);
			l += 1;
		}
		aux += l;
	}
}

static void
pack_data (xread_t * r, uint8_t ** data, int32_t * l_data, int32_t * m_data)
{
	char * ch;
	int32_t i;
	int32_t l_qname;
	int32_t n_cigar;
	int32_t l_qseq;
	str_t s;
	aux_item_t * item;

	s.l = 0;
	s.m = *m_data;
	s.s = (char *) (*data);
	memset (s.s, 0, s.m);

	l_qname = r->name->l + 1;
	n_cigar = r->aln->cigar->n;
	l_qseq  = r->l;

	str_resize (&s, l_qname + (n_cigar<<2) + l_qseq + (l_qseq+1)>>1);

	// query name
	memcpy (s.s, r->name->s, r->name->l+1); s.l += l_qname;

	// cigar
	memcpy (s.s+s.l, r->aln->cigar->cigar, n_cigar<<2); s.l += n_cigar<<2;

	// sequence
	ch = s.s + s.l;
  memset (ch, 0, (l_qseq+1)>>1);
	for (i=0; i<r->l; ++i)
		ch[i>>1] |= seq_nt16_table[(int)(r->bases[i])] << ((~i&1)<<2);
	s.l += (l_qseq+1)>>1;

	// quality
	memcpy (s.s+s.l, r->quals, l_qseq); s.l += l_qseq;

	// auxiliary
	for (i=0; i<aux_cnt(r->aln->aux); ++i) {
		item = aux_at (r->aln->aux, i);
		aux2hts_aux (item, &s);
	}

	*data   = (uint8_t*)s.s;
	*l_data = s.l;
	*m_data = s.m;
}

void
read_dump_fq (FILE * fp, xread_t * r)
{
	int32_t i;

	if (r == NULL)
		err_mesg ("[%s] r==NULL!", __func__);
}

static int
dump_data (BGZF * fp, int32_t l_data, uint8_t * data, int32_t * x)
{
	int i;
	int ok;
	int32_t y;
	int32_t block_len;

	block_len = l_data + 32;
	ok = (bgzf_flush_try(fp, 4+block_len) >= 0);
	if (fp->is_be) {
		for (i=0; i<8; ++i)
			swap_endian_4p (x+i);
		y = block_len;
	if (ok)
			ok = (bgzf_write(fp, swap_endian_4p(&y), 4) >= 0);
	} else {
		if (ok)
			ok = (bgzf_write(fp, &block_len, 4) >= 0);
	}
	if (ok)
		ok = (bgzf_write(fp, x, 32) >= 0);
	if (ok)
		ok = (bgzf_write(fp, data, l_data) >= 0);

	return ok ? 4+block_len : -1;
}

read_file_t *
read_open (const char * file, const char * mode)
{
	read_file_t * fp;

	fp = (read_file_t *) ckalloc (1, sizeof(read_file_t));

	fp->fp = sam_open (file, mode);

	if (strchr(mode,'r'))
		fp->header = sam_hdr_read (fp->fp);

	fp->l_data = 0;
	fp->m_data = 512;
	fp->data   = (uint8_t *) ckalloc (512, 1);

	return fp;
}

int
read_close (read_file_t * fp)
{
	sam_close (fp->fp);
	free (fp->data);
	free (fp);
}

int
read_hdr_write (read_file_t * out, read_file_t * in)
{
	return sam_hdr_write (out->fp, in->header);
}

int
read_read (read_file_t * in, xread_t * r)
{
	uint8_t * s;
	uint8_t * data;
	int32_t i;
	int32_t ret;
	int32_t l_blk;
	int32_t l_name;
	int32_t n_cigar;
	int32_t l_data;
	uint32_t x[8];
	BGZF * fp;
	sam_core_t * core;

	fp = in->fp->fp.bgzf;
	core = &r->aln->core;

	if ((ret = bgzf_read(fp, &l_blk, 4)) != 4) {
		if (ret == 0)
			return -1; // normal end-of-file
		else
			return -2; // truncated
	}

	if (bgzf_read(fp, x, 32) != 32)
		return -3;
	if (fp->is_be) {
		swap_endian_4p (&l_blk);
		for (i=0; i<8; ++i)
			swap_endian_4p (x+i);
	}

	core->tid   = x[0];
	core->pos   = x[1];
	core->qual  = x[2]>>8 & 0xff;
	l_name      = x[2] & 0xff;
	core->flag  = x[3] >> 16;
	n_cigar     = x[3] & 0xffff;
	r->l        = x[4];
	core->mtid  = x[5];
	core->mpos  = x[6];
	core->isize = x[7];

	l_data = l_blk - 32;
	if (l_data<0 || r->l<0 || l_name<1)
		return -4;
	if (l_data < (n_cigar<<2) + l_name + (r->l+1)>>1 + r->l)
		return -4;

	if (in->m_data < l_data) {
		in->m_data = l_data;
		kroundup32 (in->m_data);
		in->data = (uint8_t *) ckrealloc (in->data, in->m_data);
		if (in->data == NULL)
			return -4;
	}

	data = in->data;

	if (bgzf_read(fp, data, l_data) != l_data)
		return -4;

	if (fp->is_be)
		swap_bam_data (data, 0, l_data, l_name, n_cigar, r->l);

	unpack_data (r, data, l_data, l_name, n_cigar, r->l);

	return 0;
}

int
read_dump (read_file_t * out, xread_t * r)
{
	int32_t l;
	int32_t l_ref_consume;
	uint32_t bin;
	uint32_t x[8];
	sam_core_t * core;
	BGZF * fp;

	core = &r->aln->core;
	fp = out->fp->fp.bgzf;

	if (r->aln->cigar->n <= 0)
		l_ref_consume = 1;
	else
		l_ref_consume = cigar2ref_len (r->aln->cigar->n, r->aln->cigar->cigar);
	bin = hts_reg2bin (core->pos, core->pos+l_ref_consume, 14, 5);

	x[0] = core->tid;
	x[1] = core->pos;
	x[2] = (bin<<16) | (core->qual<<8) | (r->name->l+1);
	x[3] = (core->flag<<16) | r->aln->cigar->n;
	x[4] = r->l;
	x[5] = core->mtid;
	x[6] = core->mpos;
	x[7] = core->isize;

	pack_data (r, &out->data, &out->l_data, &out->m_data);
	if (fp->is_be)
		swap_bam_data (out->data, 1, out->l_data, r->name->l+1, r->aln->cigar->n, r->l);
	l = dump_data (fp, out->l_data, out->data, (int32_t*)x);
	if (fp->is_be)
		swap_bam_data (out->data, 0, out->l_data, r->name->l+1, r->aln->cigar->n, r->l);

	return l;
}

hts_itr_t *
read_itr_query (const hts_idx_t * idx, int32_t tid, int32_t beg, int32_t end)
{
	return hts_itr_query (idx, tid, beg, end, NULL);
}

int
read_itr_next (read_file_t * rfp, read_itr_t * itr, xread_t * r)
{
	int ret;
	int32_t tid;
	int32_t beg;
	int32_t end;
	BGZF * fp;

	if (itr==NULL || itr->finished)
		return -1;

	if (itr->off == 0)
		return -1;

	fp = rfp->fp->fp.bgzf;
	for (;;) {
		if (itr->curr_off==0 || itr->curr_off>=itr->off[itr->i].v) {
			if (itr->i == itr->n_off-1) {
				ret = -1;
				break;
			}
			if (itr->i<0 || itr->off[itr->i].v!=itr->off[itr->i+1].u) {
				ret = bgzf_seek (rfp->fp->fp.bgzf, itr->off[itr->i+1].u, SEEK_SET);
				itr->curr_off = bgzf_tell (rfp->fp->fp.bgzf);
			}
			++itr->i;
		}
		if ((ret = read_read(rfp, r)) >= 0) {
			itr->curr_off = bgzf_tell (rfp->fp->fp.bgzf);
			tid = r->aln->core.tid;
			beg = r->aln->core.pos;

			if (tid!=itr->tid || beg>itr->end) {
				ret = -1;
				break;
			}
			end = beg + cigar2ref_len(r->aln->cigar->n, r->aln->cigar->cigar);
			if (end>itr->beg && itr->end>beg) {
				itr->curr_tid = tid;
				itr->curr_beg = beg;
				itr->curr_end = end;
				return ret;
			}
		} else
			break;
	}

	itr->finished = 1;
	return ret;
}
