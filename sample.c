/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2019-12-31 13:21:42
  *Edit History: 
***********************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "read.h"
#include "sample.h"

static void
parse_info (str_t * s, char * line)
{
  char * c;

  if ((c=strchr(line,'=')) == NULL)
    err_mesg ("invalid info: '%s'", line);

  c = skip_blanks (c+1);
  str_assign (s, c);
}

static void
parse_adapter_info (str_t * ad, str_t * ad_str, char * s)
{
	char * c;

	if ((c=strchr(s,'=')) == NULL)
    err_mesg ("invalid info: '%s'", s);

	c = skip_blanks (c+1);
	bio_spl_adapter_parse (ad, ad_str, c);
}

bio_spl_t *
bio_spl_load (const char * spl_cfg)
{
  char * c;
  char * line;
  char * path;
  FILE * in;
	bio_spl_t * spl;
  bio_lib_t * lib;
  bio_fq_t * fq;

	spl = (bio_spl_t *) ckmalloc (sizeof(bio_spl_t));
	bio_spl_init2 (spl);

  lib = NULL;
  fq = NULL;
  line = ALLOC_LINE;
  path = ALLOC_LINE;
  in = ckopen (spl_cfg, "r");
  while (fgets(line,LINE_MAX,in)) {
    chomp (line);
    c = skip_blanks (line);
    if (*c == '[') {
      if (strncmp(c+1,"lib",3) == 0)
        lib = mp_alloc (blib, spl->libs);
      else if (strncmp(c+1,"lane",4) == 0) {
        assert (lib != NULL);
        fq = mp_alloc (bfq, lib->fqs);
      }
    } else {
      if (strncmp(c,"LB",2)==0 || strncmp(c,"lb",2)==0) {
        assert (lib != NULL);
        parse_info (lib->name, c);
      } else if (strncmp(c,"fq1",3) == 0) {
        assert (fq != NULL);
        parse_info (fq->fq, c);
      } else if (strncmp(c,"fq2",3) == 0) {
        assert (fq != NULL);
        parse_info (fq->fq+1, c);
      } else if (strncmp(c,"ad1_str",7) == 0) {
        assert (fq != NULL);
        parse_adapter_info (fq->ad, fq->ad_str, c);
      } else if (strncmp(c,"ad2_str",7) == 0) {
        assert (fq != NULL);
        parse_adapter_info (fq->ad+1, fq->ad_str+1, c);
      } else if (strncmp(c,"ad1",3) == 0) {
        assert (fq != NULL);
        parse_adapter_info (fq->ad, fq->ad_str, c);
      } else if (strncmp(c,"ad2",3) == 0) {
        assert (fq != NULL);
        parse_adapter_info (fq->ad+1, fq->ad_str+1, c);
      }
    }
  }
  fclose (in);
  free (line);
  free (path);

	return spl;
}

int
bio_spl_check (bio_spl_t * spl)
{
	char * line;
	int32_t n_fqs;
  int64_t i, j;
  bio_lib_t * lib;
	bio_fq_t * fq;

	n_fqs = 0;
	line = ALLOC_LINE;
  for (i=0; i<mp_cnt(spl->libs); ++i) {
    lib = mp_at (blib, spl->libs, i);
    if (str_is_empty(lib->name)) {
      sprintf (line, "lib%ld", i);
      str_assign (lib->name, line);
    }
    for (j=0; j<mp_cnt(lib->fqs); ++j,++n_fqs) {
      fq = mp_at (bfq, lib->fqs, j);

      // check fastq files
      if (str_is_empty(fq->fq))
        err_mesg ("fastq file for read1 must be provided!");
      if (!str_is_empty(fq->fq+1))
        fq->flag |= BIO_FASTQ_PE;

      // check information for adapters
      if (str_is_empty(fq->ad)) { // no adapter lists, then check adapter sequence
        if (!str_is_empty(fq->ad_str)) {
          if ((fq->flag & BIO_FASTQ_PE)
              && (str_is_empty(fq->ad_str+1)))
            err_mesg ("backward adapter sequence must be provided!");
          fq->flag |= BIO_FASTQ_AD_STR;
        }
      } else {
        if ((fq->flag & BIO_FASTQ_PE)
            && (str_is_empty(fq->ad+1)))
          err_mesg ("list for adatper 2 must be provided!");
        fq->flag |= BIO_FASTQ_AD_LIST;
      }

      // misc
      if ((fq->rd_len[0] = read_length_estimate(fq->fq[0].s)) < 0)
        err_mesg ("invalid fastq file: '%s'!", fq->fq[0].s);

      if ((fq->flag & BIO_FASTQ_PE)
          && (fq->rd_len[1] = read_length_estimate(fq->fq[1].s)) < 0)
        err_mesg ("invalid fastq file: '%s'!", fq->fq[1].s);
    }
  }
  free (line);

	if (n_fqs == 0)
		err_mesg ("No fastq files are provided!\n");

	return 0;
}

void
bio_spl_free (bio_spl_t * spl)
{
  bio_spl_free2 (spl);
  free (spl);
}

void
bio_spl_adapter_parse (str_t * ad, str_t * ad_str, char * s)
{
	char * c;

	if (access(s,0) == 0) {
		str_assign (ad, s);
		return;
	}

	for (c=s; *c!='\0'; ++c) {
		*c = toupper (*c);
		if (*c != 'A'
				&& *c != 'C'
				&& *c != 'T'
				&& *c != 'G')
			err_mesg ("invalid adapter info: '%s'", s);
	}

	str_assign (ad_str, s);
}

void
bio_spl_dump (FILE * out, bio_spl_t * spl)
{
	int32_t i, j;
	bio_lib_t * lib;
	bio_fq_t * fq;

	fprintf (out, "\n");
  fprintf (out, "total %ld libs\n", mp_cnt(spl->libs));
	for (i=0; i<mp_cnt(spl->libs); ++i) {
		lib = mp_at (blib, spl->libs, i);
		fprintf (out, " %s\n", lib->name->s);
		for (j=0; j<mp_cnt(lib->fqs); ++j) {
			fq = mp_at (bfq, lib->fqs, j);
			if (fq->flag & BIO_FASTQ_PE) {
				fprintf (out, "  pair end reads\n");
				fprintf (out, "  fq1: %s\n", fq->fq[0].s);
				fprintf (out, "  fq2: %s\n", fq->fq[1].s);
				if (fq->flag & BIO_FASTQ_AD_STR) {
					fprintf (out, "  forw ad seq: %s\n", fq->ad_str[0].s);
					fprintf (out, "  back ad seq: %s\n", fq->ad_str[1].s);
				} else if (fq->flag & BIO_FASTQ_AD_LIST) {
					fprintf (out, "  forw ad list: %s\n", fq->ad[0].s);
					fprintf (out, "  back ad list: %s\n", fq->ad[1].s);
				}
			} else {
				fprintf (out, "  single end reads\n");
				fprintf (out, "  fq: %s\n", fq->fq[0].s);
				if (fq->flag & BIO_FASTQ_AD_STR)
					fprintf (out, "  forw ad seq: %s\n", fq->ad_str[0].s);
				else if (fq->flag & BIO_FASTQ_AD_LIST)
					fprintf (out, "  forw ad list: %s\n", fq->ad[0].s);
			}
			fprintf (out, "\n");
		}
	}
}
