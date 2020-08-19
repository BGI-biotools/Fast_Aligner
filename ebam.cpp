#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <hts_endian.h>

#include "mp.h"
#include "xth.h"
#include "var.h"
#include "ebam.h"
#include "sort.h"
#include "cigar.h"
#include "sample.h"
#include "sam_bio.h"

/*---------------------------------------------------------------------------*/
/*---------------------------- Static Functions -----------------------------*/
/*---------------------------------------------------------------------------*/

static inline void
bam1_data_copy (bam1_t * b, str_t * buf)
{
	uint8_t * data;
	int32_t m_data;

	data = b->data;
	m_data = b->m_data;
	if (m_data < buf->l) {
		m_data = buf->l;
		kroundup32 (m_data);
		data = (uint8_t *) ckrealloc (data, m_data);
	}
	memcpy (data, buf->s, buf->l);
	b->l_data = buf->l;
	b->m_data = m_data;
	b->data   = data;
}

static inline void
add_cigar (const mem_opt_t * opt, mem_aln_t * p, str_t * buf, int which)
{
	int i;

	for (i=0; i<p->n_cigar; ++i) {
		uint32_t elem;
		int opr = p->cigar[i] & 0xf;
		int len = p->cigar[i] >> 4;
		if (!(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt && (opr == 3 || opr == 4))
			opr = which? 4 : 3; // use hard clipping for supplementary alignments
		elem = (len << 4) | (cigar_bwa2samtools[opr]);
		memcpy (buf->s+buf->l, &elem, 4);
		buf->l += 4;
	}
}

static inline void
add_cigar_str (const mem_opt_t * opt, mem_aln_t * p, str_t * buf, int which)
{
	char num[16];
	int i;

	if (p->n_cigar <= 0) {
		str_add (buf, '*');
		return;
	}

	for (i=0; i<p->n_cigar; ++i) {
		uint32_t elem;
		int opr = p->cigar[i] & 0xf;
		if (!(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt && (opr == 3 || opr == 4))
			opr = which? 4 : 3; // use hard clipping for supplementary alignments
		sprintf (num, "%d", p->cigar[i] >> 4);
		str_append (buf, num, strlen(num));
		str_add (buf, "MIDSH"[opr]);
	}
}

/*---------------------------------------------------------------------------*/
/*----------------------------- ESAM Functions ------------------------------*/
/*---------------------------------------------------------------------------*/

int
mem_aln2ebam (const mem_opt_t *opt, const bntseq_t *bns,
    bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_, int pid)
{
  char * ch;
  int i, j, l;
  int32_t rlen[2];
  str_t * buf;
  bam1_core_t * c;
  ebam_t * r;

  mem_aln_t ptmp = list[which], *p = &ptmp, mtmp, *m = 0;

  if (m_) mtmp = *m_, m = &mtmp;

  buf = pth_str_buf + pid;
	str_clear (buf);

	r = bmp_alloc (erd, raw_reads[pid]);
	arr_add (erp, s->res->rds, r);
	c = &r->b.core;

  if (p->rid>=0 && p->n_cigar>0)
    rlen[0] = bwa_cigar2ref_len (p->n_cigar, p->cigar);
  else
    rlen[0] = -1;

  if (m->rid>=0 && m->n_cigar>0)
    rlen[1] = bwa_cigar2ref_len (m->n_cigar, m->cigar);
  else
    rlen[1] = -1;

	// set flag
	p->flag |= m? 0x1 : 0; // is paired in sequencing
	p->flag |= p->rid < 0? 0x4 : 0; // is mapped
	p->flag |= m && m->rid < 0? 0x8 : 0; // is mate mapped
	if (p->rid < 0 && m && m->rid >= 0) // copy mate to alignment
		p->rid = m->rid, p->pos = m->pos, p->is_rev = m->is_rev, p->n_cigar = 0;
	if (m && m->rid < 0 && p->rid >= 0) // copy alignment to mate
		m->rid = p->rid, m->pos = p->pos, m->is_rev = p->is_rev, m->n_cigar = 0;
	p->flag |= p->is_rev? 0x10 : 0; // is on the reverse strand
	p->flag |= m && m->is_rev? 0x20 : 0; // is mate on the reverse strand

  // copy basic info
  c->flag  = p->flag;
  c->tid   = p->rid;
  c->pos   = p->pos;
  c->qual  = p->mapq;

  // query name
  c->l_qname = strlen(s->name) + 1;
  memcpy (buf->s, s->name, c->l_qname);
  buf->l += c->l_qname;

  c->l_extranul = (4 - (c->l_qname&3)) & 3;
  memcpy (buf->s+buf->l, "\0\0\0\0", c->l_extranul);
  buf->l += c->l_extranul;
  c->l_qname += c->l_extranul;

  // cigar
  c->n_cigar = p->n_cigar;
  if (rlen[0] >= 0) {
		add_cigar (opt, p, buf, which);

    l = (!(c->flag&BAM_FUNMAP)) ? rlen[0] : 1;
    c->bin = hts_reg2bin (c->pos, c->pos+l, 14, 5);
  }

  // mate information
  c->mtid = m->rid;
  c->mpos = m->pos;
  if (p->rid==m->rid && p->rid!=-1) {
    if (m->n_cigar==0 || p->n_cigar==0) {
      c->isize = 0;
    } else {
      int64_t p0 = p->pos + (p->is_rev? rlen[0] - 1 : 0);
      int64_t p1 = m->pos + (m->is_rev? rlen[1] - 1 : 0);
      c->isize = -(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0));
    }
  } else
    c->isize = 0;

	// print SEQ and QUAL
	if (p->flag & 0x100) { // for secondary alignments, don't write SEQ and QUAL
    c->l_qseq = 0;
	} else if (!p->is_rev) { // the forward strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) { // have cigar && not the primary alignment && not softclip all
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qb += p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qe -= p->cigar[p->n_cigar-1]>>4;
		}
    ch = buf->s + buf->l;
    c->l_qseq = l = qe - qb;
    memset (ch, 0, (l+1)>>1);
    for (i=qb,j=0; i<qe; ++i,++j) {
      unsigned char base = "ACGTN"[(int)s->seq[i]];
      ch[j>>1] |= seq_nt16_table[base] << ((~j&1)<<2);
    }
    buf->l += (l+1) >> 1;
		if (s->qual) { // printf qual
			for (i=qb; i<qe; ++i)
				buf->s[buf->l++] = s->qual[i] - BASEQ_PHRED_OFFSET;
    }
	} else { // the reverse strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) {
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qe -= p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qb += p->cigar[p->n_cigar-1]>>4;
		}
    ch = buf->s + buf->l;
    c->l_qseq = l = qe - qb;
    memset (ch, 0, (l+1)>>1);
    for (i=qe-1,j=0; i>=qb; --i,++j) {
      unsigned char base = "TGCAN"[(int)s->seq[i]];
      ch[j>>1] |= seq_nt16_table[base] << ((~j&1)<<2);
    }
    buf->l += (l+1) >> 1;
		if (s->qual) { // printf qual
      for (i=qe-1; i>=qb; --i)
        buf->s[buf->l++] = s->qual[i] - BASEQ_PHRED_OFFSET;
    }
	}

	// print optional tags
	if (p->n_cigar) {
    buf->s[buf->l++] = 'N';
    buf->s[buf->l++] = 'M';
    buf->s[buf->l++] = 'i';
    i32_to_le (p->NM, (uint8_t*)(buf->s+buf->l));
    buf->l += 4;

    buf->s[buf->l++] = 'M';
    buf->s[buf->l++] = 'D';
    buf->s[buf->l++] = 'Z';
    ch = (char*)(p->cigar + p->n_cigar);
    l = strlen (ch) + 1;
    memcpy (buf->s+buf->l, ch, l);
    buf->l += l;

		if (cal_mc && m && m->n_cigar) {
			buf->s[buf->l++] = 'M';
			buf->s[buf->l++] = 'C';
			buf->s[buf->l++] = 'Z';
			add_cigar_str (opt, m, buf, which);
			str_add (buf, '\0');
		}
	}

  if (p->score >= 0) {
    buf->s[buf->l++] = 'A';
    buf->s[buf->l++] = 'S';
    buf->s[buf->l++] = 'i';
    i32_to_le (p->score, (uint8_t*)(buf->s+buf->l));
    buf->l += 4;
  }

  if (p->sub >= 0) {
    buf->s[buf->l++] = 'X';
    buf->s[buf->l++] = 'S';
    buf->s[buf->l++] = 'i';
    i32_to_le (p->sub, (uint8_t*)(buf->s+buf->l));
    buf->l += 4;
  }

  if (bwa_rg_id[0]) {
    buf->s[buf->l++] = 'R';
    buf->s[buf->l++] = 'G';
    buf->s[buf->l++] = 'Z';
    l = strlen (bwa_rg_id) + 1;
    memcpy (buf->s+buf->l, bwa_rg_id, l);
    buf->l += l;
  }

	if (!(p->flag & 0x100)) { // not multi-hit
		for (i = 0; i < n; ++i)
			if (i != which && !(list[i].flag&0x100)) break;
		if (i < n) { // there are other primary hits; output them
      buf->s[buf->l++] = 'S';
      buf->s[buf->l++] = 'A';
      buf->s[buf->l++] = 'Z';
			for (i = 0; i < n; ++i) {
				const mem_aln_t *r = &list[i];
				int k;
				if (i == which || (r->flag&0x100)) continue; // proceed if: 1) different from the current; 2) not shadowed multi hit

        ch = buf->s + buf->l;
        sprintf (ch, "%s,%ld,%c,", bns->anns[r->rid].name, r->pos+1, "+-"[r->is_rev]);
        buf->l += strlen (ch);

				for (k = 0; k < r->n_cigar; ++k) {
          ch = buf->s + buf->l;
          sprintf (ch, "%d%c", r->cigar[k]>>4, "MIDSH"[r->cigar[k]&0xf]);
          buf->l += strlen (ch);
				}

        ch = buf->s + buf->l;
        sprintf (ch, ",%d,%d;", r->mapq, r->NM);
        buf->l += strlen (ch);
			}
    	++buf->l; // including the tailing NULL
		}

		if (p->alt_sc > 0) {
      buf->s[buf->l++] = 'p';
      buf->s[buf->l++] = 'a';
      buf->s[buf->l++] = 'f';
      float_to_le ((double)p->score / p->alt_sc, (uint8_t*)(buf->s+buf->l));
      buf->l += sizeof(float);
    }
	}

  if (p->XA) {
    ch = buf->s + buf->l;
    sprintf (ch, "XAZ%s", p->XA);
    buf->l += strlen(ch) + 1; // including the tailling NULL
  }
	
	int64_t pos;
	uint64_t tid;

	// create key for sort
	tid = c->tid==-1 ? bns->n_seqs : c->tid;
	r->key = (uint64_t)tid<<32 | (c->pos+1)<<1 | ((c->flag&BAM_FREVERSE)!=0);

  if (c->tid<0 || c->pos<0) {
		r->out = unmapped_ebam_file->fp;
		__sync_fetch_and_add (&(unmapped_ebam_file->n_reads), 1);
  } else {
    i = c->pos >> BIN_SHIFT;    
    assert (i < n_chr_bins[c->tid]);
		r->out = batch_files[c->tid][i].fp;
		__sync_fetch_and_add (&(batch_files[c->tid][i].n_reads), 1);
  }

	// create information for markdup
	if (c->flag & 4
			|| c->flag & 256
			|| c->flag & 2048) {
		r->end[0] = -1;
		bam1_data_copy (&r->b, buf);
		return 0;
	}

	pos = c->flag&BAM_FREVERSE
            ? bwa_get_uncliped_end(c->pos+rlen[0], p->n_cigar, p->cigar)-1
            : bwa_get_uncliped_begin(c->pos, p->n_cigar, p->cigar);
	r->end[0] = c->tid;
	r->end[0] = (r->end[0] << 32) | pos;
	r->frag_end = r->end[0];

  if (c->flag & BAM_FMUNMAP || rlen[1] < 0 || c->mtid < 0 || r->end[0] < 0) {
    r->end[1] = -1;
  } else {
    pos = (c->flag & BAM_FMREVERSE)
              ? bwa_get_uncliped_end(c->mpos+rlen[1], m->n_cigar, m->cigar)-1
              : bwa_get_uncliped_begin(c->mpos, m->n_cigar, m->cigar);
		r->end[1] = c->mtid;
    r->end[1] = (r->end[1] << 32) | pos;
  }

  int64_t tmp, o1, o2;
	r->flag = 0;
	r->pair_orit = 0;
	r->frag_orit = 0;
	r->tile = r->x = r->y = -1;
	read_name_parse (s->name, c->l_qname, &r->tile, &r->x, &r->y);
  if (r->end[0]!=-1 && r->end[1]!=-1) { // read mapped and mate mapped
    if (r->end[0] > r->end[1]) {
      tmp = r->end[0];
      r->end[0] = r->end[1];
      r->end[1] = tmp;

			r->flag |= EBAM_IS_SEC;
	    o1 = c->flag&32 ? 0 : 1;
	    o2 = c->flag&16 ? 0 : 1;
    } else {
	    o1 = c->flag&16 ? 0 : 1;
	    o2 = c->flag&32 ? 0 : 1;
    }

		//if (r->end[0]==r->end[1] && c->flag&128)
		//	r->flag |= EBAM_IS_SEC;

		if (r->end[0] == r->end[1]) {
			uint64_t mate_key = (uint64_t)m->rid<<32 | (m->pos+1)<<1 | ((c->flag&BAM_FMREVERSE)!=0);
			if (r->key > mate_key
					|| (r->key==mate_key && c->flag&128))
				r->flag |= EBAM_IS_SEC;
		}

		r->flag |= EBAM_IS_PAIR;
		r->pair_orit = (o1<<1) | o2;
		r->pair_orit += 1;
  }

	if (r->frag_end != -1) {
		r->frag_orit = c->flag&16 ? 0 : 1;
		r->frag_orit += 1;
  }

	bam1_data_copy (&r->b, buf);

  return 0;
}

int
ebam_reset (bseq1_t * seqs, int n_seqs)
{
	int i;
	bmp_t(erd) * reads;

	for (i=0; i<n_threads; ++i)
		bmp_clear (erd, raw_reads[i], ebam_clear);

	mp_clear (ar, aln_res_set, aln_res_clear);
	mp_resize (ar, aln_res_set, n_seqs);
	aln_res_set->n = n_seqs;
	for (i=0; i<n_seqs; ++i) {
		seqs[i].res = aln_res_set->pool + i;
	}

	return 0;
}

int
ebam_set_dump (bseq1_t * seqs, int n_seqs)
{
	int32_t i, j;
	int32_t n;
	ebam_t ** r;
	bseq1_t * s;

	for (i=0; i<n_seqs; ++i) {
		s = seqs + i;

		n = s->res->rds->n;
		r = s->res->rds->arr;
		for (j=0; j<n; ++j)
			ebam_dump (r[j]);
	}

	return 0;
}

void
ebam_mkd_score_calc (bseq1_t * s1, bseq1_t * s2)
{
	int r1_is_mapped;
	int r2_is_mapped;
	int32_t i;
	int32_t n;
	int64_t r1_score;
	int64_t r2_score;
	ebam_t * r1_prim;
	ebam_t * r2_prim;
	ebam_t ** rds;

	r1_score = read_baseq_sum (s1->qual, s1->l_seq, BASEQ_PHRED_OFFSET);
	if (r1_score > MAX_READ_BASEQ_SUM)
		r1_score = MAX_READ_BASEQ_SUM;

	// single end
	if (s2 == NULL) {
		assert (arr_cnt(s1->res->rds) > 0);
		r1_prim = s1->res->rds->arr[0];
		if (!(r1_prim->b.core.flag & 4))
			r1_prim->score = r1_score;
		return;
	}

	r2_score = read_baseq_sum (s2->qual, s2->l_seq, BASEQ_PHRED_OFFSET);
	if (r2_score > MAX_READ_BASEQ_SUM)
		r2_score = MAX_READ_BASEQ_SUM;

	assert (arr_cnt(s1->res->rds) > 0);
	r1_prim = s1->res->rds->arr[0];
	r1_is_mapped = !(r1_prim->b.core.flag & 4);

	assert (arr_cnt(s2->res->rds) > 0);
	r2_prim = s2->res->rds->arr[0];
	r2_is_mapped = !(r2_prim->b.core.flag & 4);

	if (r1_is_mapped) {
		if (r2_is_mapped)
			r1_prim->score = r2_prim->score = r1_score + r2_score;
		else
			r1_prim->score = r1_score;
	} else {
		if (r2_is_mapped)
			r2_prim->score = r2_score;
	}
}

#define FN_ALIGN_DONE_TAG ".alignment_done.tag"
#define FN_EBAM_FILE_LIST "ebam_file.list"

/*
 * // file for unmapped alignments
 * -1      -1      -1      n_reads  path
 *
 * // file for proper alignments
 * lib_id  chr_id  bin_id  n_reads  path
 */

int
load_ebam_files (void)
{
  char line[4096];
  char file[4096];
	int32_t i, j;
  int32_t lib_id;
  int32_t chr_id;
  int32_t bin_id;
  int32_t n_reads;
	FILE * in;
  ebam_file_t * ef;

  sprintf (file, "%s/%s", out_dir, FN_ALIGN_DONE_TAG);
  if (access(file,0) != 0)
    return -1;

  // init ebam files
  unmapped_ebam_file = (ebam_file_t *) ckmalloc (sizeof(ebam_file_t));
  mapped_ebam_files = (ebam_file_t ***) ckmalloc (n_libs * sizeof(ebam_file_t **));
  for (i=0; i<n_libs; ++i) {
    mapped_ebam_files[i] = (ebam_file_t **) ckmalloc (n_chr * sizeof(ebam_file_t *));
    for (j=0; j<n_chr; ++j)
      mapped_ebam_files[i][j] = (ebam_file_t *) ckmalloc (n_chr_bins[j] * sizeof(ebam_file_t));
  }

	sprintf (file, "%s/%s", out_dir, FN_EBAM_FILE_LIST);
	in = ckopen (file, "r");

	// load information for unmapped reads
	fgets (line, 4096, in);
  sscanf (line, "%*s %*s %*s %d %s", &n_reads, file);
  unmapped_ebam_file->path = strdup (file);
  unmapped_ebam_file->n_reads = n_reads;

  while (fgets(line, 4096, in)) {
    sscanf (line, "%d %d %d %d %s", &lib_id, &chr_id, &bin_id, &n_reads, file);
    ef = mapped_ebam_files[lib_id][chr_id] + bin_id;
    ef->path = get_abs_path (file);
    ef->n_reads = n_reads;
  }
	fclose (in);

  return 0;
}

int
ebam_list_dump (void)
{
	char file[4096];
	int32_t i, j, k;
	FILE * out;
	ebam_file_t * ef;

	sprintf (file, "%s/%s", out_dir, FN_ALIGN_DONE_TAG);
	out = ckopen (file, "w");
	fclose (out);

	sprintf (file, "%s/%s", out_dir, FN_EBAM_FILE_LIST);
	out = ckopen (file, "w");
	fprintf (out, "-1\t-1\t-1\t%ld\t%s\n", unmapped_ebam_file->n_reads, unmapped_ebam_file->path);

	for (i=0; i<n_libs; ++i) {
		for (j=0; j<n_chr; ++j) {
			for (k=0; k<n_chr_bins[j]; ++k) {
				ef = mapped_ebam_files[i][j] + k;
				fprintf (out, "%d\t%d\t%d\t%ld\t%s\n", i, j, k, ef->n_reads, ef->path);
			}
		}
	}
	fclose (out);

	return 0;
}

int
ebam_list_remove (void)
{
  char cmd[4096];

  sprintf (cmd, "rm -f %s/%s", out_dir, FN_ALIGN_DONE_TAG);
  cksystem (cmd);

  sprintf (cmd, "rm -f %s/%s", out_dir, FN_EBAM_FILE_LIST);
  cksystem (cmd);

  return 0;
}

void
ebam_dump (ebam_t * r)
{
	uint8_t orit;
	int32_t l_data;
	FILE * out;

	out = r->out;
	l_data = r->b.l_data;
	ckfwrite (&r->key, 8, 1, out);

	orit = (r->pair_orit<<4) | r->frag_orit;
	ckfwrite (&orit, 1, 1, out);
	ckfwrite (&r->flag, 1, 1, out);

	ckfwrite (&r->tile, 4, 1, out);
	ckfwrite (&r->x, 4, 1, out);
	ckfwrite (&r->y, 4, 1, out);

	ckfwrite (r->end, 8, 2, out);
	ckfwrite (&r->frag_end, 8, 1, out);
	ckfwrite (&r->score, 4, 1, out);
	ckfwrite (&r->b.core, sizeof(bam1_core_t), 1, out);
	ckfwrite (&l_data, 4, 1, out);
	ckfwrite (r->b.data, 1, l_data, out);
}

int
load_reads (const char * file, ebam_t * reads, int32_t n_reads)
{
	uint8_t orit;
  uint8_t * data;
  int32_t cnter;
  int32_t l_data;
  int32_t m_data;
	uint64_t key;
  FILE * in;
  ebam_t * r;  

  cnter = 0;
  in = ckopen (file, "rb");
  for (;;) {
		if (fread(&key,8,1,in)<1 || feof(in))
			break;

    r = reads + cnter++;
    assert (cnter <= n_reads);
		r->key = key;

		ckfread (&orit, 1, 1, in);
		r->pair_orit = orit >> 4;
		r->frag_orit = orit & 0xf;
		ckfread (&r->flag, 1, 1, in);

		ckfread (&r->tile, 4, 1, in);
		ckfread (&r->x, 4, 1, in);
		ckfread (&r->y, 4, 1, in);

		ckfread (r->end, 8, 2, in);
		ckfread (&r->frag_end, 8, 1, in);
		ckfread (&r->score, 4, 1, in);
    ckfread (&r->b.core, sizeof(bam1_core_t), 1, in);
    ckfread (&l_data, 4, 1, in);

    data = r->b.data;
    m_data = r->b.m_data;
    if (m_data < l_data) {
      m_data = l_data;
      kroundup32 (m_data);
      data = (uint8_t*) ckrealloc (data, m_data);
    }
    ckfread (data, 1, l_data, in);
		r->b.l_data = l_data;
    r->b.m_data = m_data;
    r->b.data = data;
  }
  fclose (in);

  return 0;
}
