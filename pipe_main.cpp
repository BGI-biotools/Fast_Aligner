/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2019-12-27 13:58:58
  *Edit History: 
***********************************************************/

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sam.h>
#include <bgzf.h>

#include "bio.h"
#include "var.h"
#include "ebam.h"
#include "read.h"
#include "sample.h"
#include "sam_bio.h"
#include "sort_and_mark.h"
#include "mem2/bwa.h"
#include "mem2/macro.h"
#include "mem2/bwamem.h"
#include "mem2/fastmap.h"
#include "mem2/FMI_search.h"

extern FMI_search *fmi;
extern uint8_t *ref_string;
extern uint64_t proc_freq, tprof[LIM_R][LIM_C];
extern int readLen;

void update_a(mem_opt_t *opt, const mem_opt_t *opt0);
int process(void *shared, gzFile gfp, gzFile gfp2, int pipe_threads);

// global variables
char * out_dir;
char * tmp_path;
char * spl_name;
int is_be; // big end system
int n_chr;
int n_libs;
int n_threads;
int need_filt;
int cal_mc;
int * n_chr_bins;
int64_t total_reads;
int64_t mapped_reads;
int64_t unmapped_reads;
genome_info_t * ginfo;
bam_hdr_t * bam_hdr;
bio_spl_t * spl;
ebam_file_t * unmapped_ebam_file;
ebam_file_t *** mapped_ebam_files;

xrd_filt_t * fastq_filter;
str_t * pth_str_buf;
mp_t(ar) * aln_res_set;
bmp_t(erd) ** raw_reads;
ebam_file_t ** batch_files;

samFile * debug_fp;

/**********************************************************
 ******************* Parameter Analyze ********************
 **********************************************************/

struct opt_s;
typedef struct opt_s opt_t;

struct opt_s {
	char * in_cfg;
	char * out_dir;
	char * tmp_path;
	char * spl_name;
  char * ref_file;
	char * ad_list;
	char * fq_pair;
  int n_threads;

	int low_qual;
	double lq_rate;
	double n_rate;

  int64_t fixed_chunk_size;
  int ignore_alt;
	int is_forced;
	int need_mkd;
	int need_filt;
	int cal_mc; // calculate mate cigar
  mem_opt_t * mem_opt;
  mem_pestat_t * pes;
  ktp_aux_t * aux;
  gzFile fp, fp2;
};

static opt_t *
opt_init (void)
{
  int i;
  opt_t * opt;

  opt = (opt_t *) ckmalloc (sizeof(opt_t));
  opt->in_cfg = ALLOC_LINE;
  opt->out_dir = ALLOC_LINE;
  opt->tmp_path = ALLOC_LINE;
  opt->spl_name = ALLOC_LINE;
  opt->ref_file = ALLOC_LINE;
  opt->ad_list  = ALLOC_LINE;
  opt->fq_pair  = ALLOC_LINE;

	opt->n_threads = 1;
	opt->low_qual  = 5;
	opt->lq_rate   = 0.5;
	opt->n_rate    = 0.05;
	opt->need_mkd  = 1;
	opt->need_filt = 1;
	opt->cal_mc    = 0;

  opt->mem_opt = mem_opt_init ();
  opt->aux = (ktp_aux_t *) ckalloc (1, sizeof(ktp_aux_t));
  opt->pes = (mem_pestat_t *) ckalloc (4, sizeof(mem_pestat_t));
  for (i=0; i<4; ++i)
    opt->pes[i].failed = 1;

  opt->ignore_alt = 0;
	opt->is_forced = 0;
  opt->fixed_chunk_size = -1;

  return opt;
}

static int
opt_analyze (opt_t * o, int argc, char * argv[])
{
  char * p;
  int c;
  int tmp;
  double dtmp;
  mem_opt_t opt0;
  mem_opt_t * opt;
  ktp_aux_t * aux;
  mem_pestat_t * pes;

  opt = o->mem_opt;
  aux = o->aux;
  aux->opt = opt;
  pes = o->pes;
  memset (&opt0, 0, sizeof(mem_opt_t));
	// remained alphabets: F g H J l q R u z Z
	while ((c = getopt(argc, argv, "F:xei:o:n:1:2:3:b:t:paMCSPVYjk:c:v:s:r:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:G:h:y:K:X:f")) != -1)
	{
    if (c == 'i')
      strcpy (o->in_cfg, optarg);
    else if (c == 'o') {
      strcpy (o->out_dir, optarg);
			chop (o->out_dir, '/');
		} else if (c == 'n')
      strcpy (o->spl_name, optarg);
    else if (c == '1') {
      if ((tmp = atoi(optarg)) < 0) {
        fprintf (stderr, "low quality threshold must >= 0");
        return -1;
      }
      o->low_qual = tmp;
    } else if (c == '2') {
      if ((dtmp = atof(optarg)) < 0) {
        fprintf (stderr, "low quality rate >= 0");
        return -1;
      }
      o->lq_rate = dtmp;
    } else if (c == '3') {
      if ((dtmp = atof(optarg)) < 0) {
        fprintf (stderr, "N rate threshold must >= 0");
        return -1;
      }
      o->n_rate = dtmp;
    } else if (c == 't') {
      o->n_threads = atoi (optarg);
      if (o->n_threads < 1)
        o->n_threads = 1;
      opt->n_threads = o->n_threads;
      assert(opt->n_threads >= INT_MIN && opt->n_threads <= INT_MAX);
		} else if (c == 'f') {
			o->is_forced = 1;
		} else if (c == 'x') {
			o->need_mkd = 0;
		} else if (c == 'b') {
			strcpy (o->ad_list, optarg);
		} else if (c == 'e') {
			o->need_filt = 0;
		} else if (c == 'F') {
			o->cal_mc = 1;
    // original bwa mem2 options:
    } else if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
		else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1, assert(opt->a >= INT_MIN && opt->a <= INT_MAX);
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1, assert(opt->b >= INT_MIN && opt->b <= INT_MAX);
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1, assert(opt->T >= INT_MIN && opt->T <= INT_MAX);
		else if (c == 'U')
			opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1, assert(opt->pen_unpaired >= INT_MIN && opt->pen_unpaired <= INT_MAX);
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
		else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
		else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
		else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
		else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'j') o->ignore_alt = 1;
		else if (c == 'r')
			opt->split_factor = atof(optarg), opt0.split_factor = 1.;
		else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
		else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
		else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
		else if (c == 'G')
			opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
		else if (c == 'N')
			opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
		else if (c == 'W')
			opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
		else if (c == 'y')
			opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;
		else if (c == 'C') aux->copy_comment = 1;
		else if (c == 'K') o->fixed_chunk_size = atoi(optarg);
		else if (c == 'X') opt->mask_level = atof(optarg);
		else if (c == 'h')
		{
			opt0.max_XA_hits = opt0.max_XA_hits_alt = 1;
			opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->max_XA_hits_alt = strtol(p+1, &p, 10);
		}
		else if (c == 'Q')
		{
			opt0.mapQ_coef_len = 1;
			opt->mapQ_coef_len = atoi(optarg);
			opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
		}
		else if (c == 'O')
		{
			opt0.o_del = opt0.o_ins = 1;
			opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->o_ins = strtol(p+1, &p, 10);
		}
		else if (c == 'E')
		{
			opt0.e_del = opt0.e_ins = 1;
			opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->e_ins = strtol(p+1, &p, 10);
		}
		else if (c == 'L')
		{
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
			opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->pen_clip3 = strtol(p+1, &p, 10);
		}
		else if (c == 'I')
		{
			aux->pes0 = pes;
			pes[1].failed = 0;
			pes[1].avg = strtod(optarg, &p);
			pes[1].std = pes[1].avg * .1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].std = strtod(p+1, &p);
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
			if (pes[1].low < 1) pes[1].low = 1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].high = (int)(strtod(p+1, &p) + .499);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].low  = (int)(strtod(p+1, &p) + .499);
		}
	}

  if (optind >= argc)
    err_mesg ("fail to load reference file");
  else
    strcpy (o->ref_file, argv[optind++]);

	if (optind < argc)
		strcpy (o->fq_pair, argv[optind++]);

	if (!*o->in_cfg && !*o->fq_pair)
		return -1;

	if (!*o->spl_name)
		strcpy (o->spl_name, "sample");

	if (!*o->out_dir)
		strcpy (o->out_dir, ".");

  update_a (opt, &opt0);

  /* Matrix for SWA */
  bwa_fill_scmat(opt->a, opt->b, opt->mat);

  /* Chunk Size */
  if (o->fixed_chunk_size > 0)
    aux->task_size = o->fixed_chunk_size;
  else
    aux->task_size = opt->chunk_size * o->n_threads;
  tprof[MISC][1] = opt->chunk_size = aux->actual_chunk_size = aux->task_size;

	// create tmperory directory
  ckcreate_dir (o->out_dir);
  sprintf (o->tmp_path, "%s/tmp", o->out_dir);
  ckcreate_dir (o->tmp_path);

	out_dir  = o->out_dir;
	tmp_path = o->tmp_path;
	spl_name = o->spl_name;

  return 0;
}

static void
opt_destroy (opt_t * opt)
{
  free (opt->in_cfg);
  free (opt->out_dir);
  free (opt->ref_file);
  free (opt->spl_name);
  free (opt->mem_opt);
  free (opt->aux);
  free (opt->pes);
  free (opt);
}

/**********************************************************
 **************** Other Static Functions ******************
 **********************************************************/

static void
read_align_init (bio_fq_t * fq, opt_t * opt)
{
  int i;
  int l_rd;
  int is_pe;
  int max_num_low_qual[2];
  int max_num_Ns[2];
  int min_num_match[2];
  int max_num_mismatch[2];

  is_pe = fq->flag & BIO_FASTQ_PE;

  max_num_low_qual[0] = (int) ceil (opt->lq_rate * (double)fq->rd_len[0]);
  max_num_Ns[0] = (int) ceil (opt->n_rate * (double)fq->rd_len[0]);
  min_num_match[0] = (int) ceil (0.5 * (double)fq->rd_len[0]);
  max_num_mismatch[0] = 1;
  if (is_pe) {
    max_num_low_qual[1] = (int) (opt->lq_rate * (double)fq->rd_len[1]);
    max_num_Ns[1] = (int) (opt->n_rate * (double)fq->rd_len[1]);
    min_num_match[1] = (int) ceil (0.5 * (double)fq->rd_len[1]);
    max_num_mismatch[1] = 1;
  }

  xrd_filter_set (fastq_filter, fq->flag, opt->low_qual,
      max_num_low_qual[0], max_num_low_qual[1],
      max_num_Ns[0], max_num_Ns[1],
      min_num_match[0], min_num_match[1],
      max_num_mismatch[0], max_num_mismatch[1],
      fq->ad[0].s, fq->ad[1].s,
      fq->ad_str[0].s, fq->ad_str[1].s);

  // READS file opterations
  opt->fp = ckgzopen (fq->fq[0].s, "r");
  opt->aux->ks = kseq_init (opt->fp);

  opt->aux->ks2 = 0;
  if (fq->flag & BIO_FASTQ_PE) {
    opt->fp2 = ckgzopen (fq->fq[1].s, "r");
    opt->aux->ks2 = kseq_init (opt->fp2);
    opt->mem_opt->flag |= MEM_F_PE;
    assert (opt->aux->ks2 != 0);
  }

  readLen = fq->rd_len[0];

	total_reads += opt->aux->n_processed;
	opt->aux->n_processed = 0;
}

static void
read_align_reset (opt_t * opt)
{
  kseq_destroy (opt->aux->ks);
  gzclose (opt->fp);

  if (opt->aux->ks2) {
    kseq_destroy (opt->aux->ks2);
    gzclose (opt->fp2);
  }
}

static void
alignment_init (void)
{
  int i, j, k;
	tmp_dir_t * tmp_dir;

	aln_res_set = mp_init (ar, aln_res_init2);

	raw_reads = (bmp_t(erd) **) ckmalloc (n_threads * sizeof(bmp_t(erd) *));
	for (i=0; i<n_threads; ++i)
		raw_reads[i] = bmp_init (erd, 4096, ebam_init2);

	tmp_dir = tmp_dir_init (tmp_path);
  mapped_ebam_files = (ebam_file_t ***) ckmalloc (n_libs * sizeof(ebam_file_t **));
	for (i=0; i<n_libs; ++i) {
		mapped_ebam_files[i] = (ebam_file_t **) ckmalloc (n_chr * sizeof(ebam_file_t *));
		for (j=0; j<n_chr; ++j) {
			mapped_ebam_files[i][j] = (ebam_file_t *) ckmalloc (n_chr_bins[j] * sizeof(ebam_file_t));
			for (k=0; k<n_chr_bins[j]; ++k) {
				mapped_ebam_files[i][j][k].path = get_abs_path (tmp_dir_add1file(tmp_dir,"ebam"));
				mapped_ebam_files[i][j][k].n_reads = 0;
			}
		}
	}

	unmapped_ebam_file = (ebam_file_t *) ckmalloc (sizeof(ebam_file_t));
	unmapped_ebam_file->path = get_abs_path (tmp_dir_add1file(tmp_dir,"ebam"));
	unmapped_ebam_file->fp = ckopen (unmapped_ebam_file->path, "wb");
	unmapped_ebam_file->n_reads = 0;

	tmp_dir_free (tmp_dir);
}

/**********************************************************
 **************** Main Function for pipe ******************
 **********************************************************/

static void
mem2_origin_usage (mem_opt_t * opt)
{
  fprintf(stderr, "Options for bwa_mem2:\n");
  fprintf(stderr, "  Algorithm options:\n");
	fprintf(stderr, "   -k INT        minimum seed length [%d]\n", opt->min_seed_len);
	fprintf(stderr, "   -w INT        band width for banded alignment [%d]\n", opt->w);
	fprintf(stderr, "   -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
	fprintf(stderr, "   -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
	fprintf(stderr, "   -y INT        seed occurrence for the 3rd round seeding [%ld]\n", (long)opt->max_mem_intv);
	fprintf(stderr, "   -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
	fprintf(stderr, "   -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->drop_ratio);
	fprintf(stderr, "   -W INT        discard a chain if seeded bases shorter than INT [0]\n");
	fprintf(stderr, "   -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
	fprintf(stderr, "   -S            skip mate rescue\n");
	fprintf(stderr, "   -o            output file name missing\n");
	fprintf(stderr, "   -P            skip pairing; mate rescue performed unless -S also in use\n");
	fprintf(stderr, "  Scoring options:\n");
	fprintf(stderr, "   -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [%d]\n", opt->a);
	fprintf(stderr, "   -B INT        penalty for a mismatch [%d]\n", opt->b);
	fprintf(stderr, "   -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
	fprintf(stderr, "   -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", opt->e_del, opt->e_ins);
	fprintf(stderr, "   -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
	fprintf(stderr, "   -U INT        penalty for an unpaired read pair [%d]\n", opt->pen_unpaired);
	fprintf(stderr, "  Input/output options:\n");
	fprintf(stderr, "   -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n");
	fprintf(stderr, "   -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
	fprintf(stderr, "   -T INT        minimum score to output [%d]\n", opt->T);
	fprintf(stderr, "   -h INT[,INT]  if there are <INT hits with score >80%% of the max score, output all in XA [%d,%d]\n", opt->max_XA_hits, opt->max_XA_hits_alt);
	fprintf(stderr, "   -a            output all alignments for SE or unpaired PE\n");
	fprintf(stderr, "   -C            append FASTA/FASTQ comment to SAM output\n");
	fprintf(stderr, "   -V            output the reference FASTA header in the XR tag\n");
	fprintf(stderr, "   -Y            use soft clipping for supplementary alignments\n");
	fprintf(stderr, "   -M            mark shorter split hits as secondary\n");
	fprintf(stderr, "   -I FLOAT[,FLOAT[,INT[,INT]]]\n");
	fprintf(stderr, "                 specify the mean, standard deviation (10%% of the mean if absent), max\n");
	fprintf(stderr, "                 (4 sigma from the mean if absent) and min of the insert size distribution.\n");
	fprintf(stderr, "                 FR orientation only. [inferred]\n");
	fprintf (stderr, "\n");
	fprintf(stderr, "  Note: 1. '-o' and '-t' are removed because they are in the list of basic options\n");
	fprintf(stderr, "        2. '-p', '-R', and '-H' are remove because this tool is designed for multiple fastq files\n");
	fprintf (stderr, "\n");
}

static int
usage (void)
{
  mem_opt_t * mem_opt;

	fprintf (stderr, "\n");
	fprintf (stderr, "Usage:  fast_aligner mem2 [options] <reference> [fq1,fq2]\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "Basic Options:\n");
	fprintf (stderr, "   -i   STR   config file for input fastqs\n");
	fprintf (stderr, "   -o   STR   output directory [./]\n");
	fprintf (stderr, "   -n   STR   sample name [sample]\n");
	fprintf (stderr, "   -t   INT   number of threads [1]\n");
	fprintf (stderr, "   -e   NUL   skip the step of filtering reads with low quality\n");
	fprintf (stderr, "   -x   NUL   skip the step of masking duplicates\n");
	fprintf (stderr, "   -b   STR   adapter list file, or adapter sequence, format: ad1,ad2\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "Advanced Options\n");
	fprintf (stderr, "   -1   INT   low quality threshold [5]\n");
	fprintf (stderr, "   -2   FLT   low quality rate [0.5]\n");
	fprintf (stderr, "   -3   FLT   N rate threshold [0.05]\n");
	fprintf (stderr, "   -f   NUL   re-run and overwrite all intermediate files\n");
	fprintf (stderr, "   -F   NUL   calculate the cigar of mate read and add it to the record of read\n");
	fprintf (stderr, "\n");

  mem_opt = mem_opt_init ();
  mem2_origin_usage (mem_opt);
  free (mem_opt);

	fprintf (stderr, "Examples:\n");
	fprintf (stderr, "   1. fast_aligner mem2 -i input.cfg -o out_dir hg19.fa\n");
	fprintf (stderr, "   2. fast_aligner mem2 -b AAGTGACAA,AAGTGAGCCAAGGAGTTG -o out_dir hg19.fa input_r1.fq.gz,input_r2.fq.gz\n");

  fprintf (stderr, "\n");

	return 0;
}

static void
load_ref4align (char * ref_file, int ignore_alt)
{
	/* Load bwt2/FMI index */
	uint64_t tim = __rdtsc();

	fprintf(stderr, "Ref file: %s\n", ref_file);
	fmi = new FMI_search(ref_file);
	tprof[FMI][0] += __rdtsc() - tim;
	
	// reading ref string from the file
	tim = __rdtsc();
	fprintf(stderr, "Reading reference genome..\n");
	
  char binary_seq_file[200];
  sprintf(binary_seq_file, "%s.0123", ref_file);
	
	fprintf(stderr, "Binary seq file = %s\n", binary_seq_file);
	FILE *fr = fopen(binary_seq_file, "r");
	
	if (fr == NULL) {
		fprintf(stderr, "Error: can't open %s input file\n", binary_seq_file);
		exit(0);
	}
	
	int64_t rlen = 0;
	fseek(fr, 0, SEEK_END); 
	rlen = ftell(fr);
	ref_string = (uint8_t*) _mm_malloc(rlen, 64);
	rewind(fr);

	/* Reading ref. sequence */
	fread(ref_string, 1, rlen, fr);

	uint64_t timer  = __rdtsc();
	tprof[REF_IO][0] += timer - tim;
	
	fclose(fr);
	fprintf(stderr, "Reference genome size: %ld bp\n", rlen);
	fprintf(stderr, "Done readng reference genome !!\n\n");

  int i;
  if (ignore_alt)
    for (i=0; i<fmi->idx->bns->n_seqs; ++i)
      fmi->idx->bns->anns[i].is_alt = 0;
}

static void
tmp_file_open (int lib_idx, int n_chr)
{
	int i, j;

  batch_files = mapped_ebam_files[lib_idx];
	for (i=0; i<n_chr; ++i)
		for (j=0; j<n_chr_bins[i]; ++j)
			batch_files[i][j].fp = ckopen (batch_files[i][j].path, "wb");
}

static void
tmp_file_close (int lib_idx)
{
	int i, j;

	for (i=0; i<n_chr; ++i)
		for (j=0; j<n_chr_bins[i]; ++j)
			fclose (batch_files[i][j].fp);
}

static samFile *
bam_file_init (int argc, char ** argv)
{
  char line[4096];
  char * s;
  int i, ret;
  bio_lib_t * lib;
  samFile * out;

  sprintf (line, "%s/%s.markdup.bam", out_dir, spl_name);
	printf ("%s\n", line);
  out = xsam_open (line, "wb");
	if (bgzf_mt(out->fp.bgzf, n_threads, 64) != 0)
		err_mesg ("[bgzf_mt failed!]");
  bam_hdr = xsam_hdr_create (ginfo);

  for (i=0; i<mp_cnt(spl->libs); ++i) {
    lib = mp_at (blib, spl->libs, i);
    sprintf (line, "@RG\tSM:%s\tID:%s\tLB:%s", spl_name, spl_name, lib->name->s);
    xsam_hdr_add_line (bam_hdr, line);
  }

  s = line;
  sprintf (s, "@PG\tID:fast_aln\tVN:0.1.0\tCL:fast_aln");
  for (i=0; i<argc; ++i) {
    s = line + strlen (line);
    sprintf (s, " %s", argv[i]);
  }
  xsam_hdr_add_line (bam_hdr, line);

  ret = sam_hdr_write (out, bam_hdr);

  return out;
}

static bio_spl_t *
bio_spl_load_lite (char * ad_list, char * fq_list)
{
	char * ch;
	bio_spl_t * spl;
	bio_lib_t * lib;
	bio_fq_t * fq;

	spl = (bio_spl_t *) ckmalloc (sizeof(bio_spl_t));
	bio_spl_init2 (spl);

	lib = mp_alloc (blib, spl->libs);
	fq = mp_alloc (bfq, lib->fqs);

	if ((ch = strchr(ad_list,',')) == NULL) {
		bio_spl_adapter_parse (fq->ad, fq->ad_str, ad_list);

		ch = strchr (fq_list, ',');
		if (*ad_list) {
			// single end
			if (ch != NULL)
				err_mesg ("Information for 2 adapters must be provided, if there are a pair of fastq files!");
			str_assign (fq->fq, fq_list);
		} else {
			// no adapters are provided
			if (ch != NULL) {
				*ch = '\0';
				str_assign (fq->fq, fq_list);
				str_assign (fq->fq+1, ch+1);
			} else
				str_assign (fq->fq, fq_list);
		}
	} else {
		// pair end
		*ch = '\0';
		bio_spl_adapter_parse (fq->ad, fq->ad_str, ad_list);
		bio_spl_adapter_parse (fq->ad+1, fq->ad_str+1, ch+1);

		ch = strchr (fq_list, ',');
		if (ch == NULL)
			err_mesg ("A pair of fastq files must be provided, if there are 2 adapters!");
		*ch = '\0';
		str_assign (fq->fq, fq_list);
		str_assign (fq->fq+1, ch+1);
	}

	return spl;
}

static samFile *
prog_init (opt_t * opt, int argc, char ** argv)
{
	int i;
	samFile * out;

	ginfo = load_ref_genome_info (opt->ref_file);

	if (*opt->in_cfg)
		spl = bio_spl_load (opt->in_cfg);
	else
		spl = bio_spl_load_lite (opt->ad_list, opt->fq_pair);
	bio_spl_check (spl);
	bio_spl_dump (stdout, spl);

	n_chr = ginfo->n_targets;
	n_libs = mp_cnt (spl->libs);
	n_threads = opt->n_threads;

  n_chr_bins = (int *) ckmalloc (n_chr * sizeof(int));
  for (i=0; i<ginfo->n_targets; ++i)
    n_chr_bins[i] = (ginfo->target_len[i]>>BIN_SHIFT) + 1;

	need_filt = opt->need_filt;
	cal_mc    = opt->cal_mc;

	out = bam_file_init (argc, argv);
	is_be = out->is_be;

	pth_str_buf = (str_t *) ckmalloc (n_threads * sizeof(str_t));
	for (i=0; i<n_threads; ++i) {
		str_init2 (pth_str_buf+i);
		str_resize (pth_str_buf+i, LONG_LINE_MAX);
	}

	total_reads = 0;
	mapped_reads = 0;
	unmapped_reads = 0;

	return out;
}

static void
align_mem_free (void)
{

}

static void
unmapped_file_close (void)
{
	fclose (unmapped_ebam_file->fp);
}

int
wgs_pipe_main (int argc, char * argv[])
{
	if (argc < 2)
		return usage ();

	int ret;
  int64_t i, j;
	time_t time_tag;
  opt_t * opt;
  bio_lib_t * lib;
  bio_fq_t * fq;
  bam_hdr_t * h;
  samFile * out;

  opt = opt_init ();
  if (opt_analyze(opt,argc,argv) != 0) {
    opt_destroy (opt);
    return usage();
  }

	// program initialization
	time (&time_tag);
	out = prog_init (opt, argc, argv);
	fprintf (stderr, "Program init costs %lds\n", time(NULL)-time_tag);

  // align reads to reference
	if (opt->is_forced || load_ebam_files()!=0) {
		time (&time_tag);
	  load_ref4align (opt->ref_file, opt->ignore_alt);
	  alignment_init ();
	  fastq_filter = xrd_filter_init ();
	  for (i=0; i<n_libs; ++i) {
	    lib = mp_at (blib, spl->libs, i);
			tmp_file_open (i, ginfo->n_targets);
	    for (j=0; j<mp_cnt(lib->fqs); ++j) {
	      fq = mp_at (bfq, lib->fqs, j);
	      read_align_init (fq, opt);
	      process (opt->aux, opt->fp, opt->fp2, 2);
	      read_align_reset (opt);
	    }
			tmp_file_close (i);
	  }
		unmapped_file_close ();
		align_mem_free ();
		ebam_list_dump ();
		fprintf (stderr, "Alignment costs %lds\n", time(NULL)-time_tag);
		delete (fmi);
	}

  // sort and markdup
	time (&time_tag);
	if (opt->need_mkd)
  	ebam_sort_and_markdup (out, spl);
	else
		ebam_sort (out, spl);
	fprintf (stderr, "Sort and markdup costs %lds\n", time(NULL)-time_tag);

  sam_close (out);
  bam_hdr_destroy (bam_hdr);

	// build bam index
	time (&time_tag);
	char file[4096];
	sprintf (file, "%s/%s.markdup.bam", out_dir, spl_name);
	ret = sam_index_build3 (file, NULL, 0, n_threads);
	fprintf (stderr, "Build bam index costs %lds\n", time(NULL)-time_tag);

  // free memory
  opt_destroy (opt);
  bio_spl_free (spl);

#if !DEBUG
  // remove tmperory files
  char cmd[4096];
  sprintf (cmd, "rm -rf %s", tmp_path);
	time (&time_tag);
  cksystem (cmd);
  ebam_list_remove ();
	fprintf (stderr, "Remove temporary files costs %lds\n", time(NULL)-time_tag);
#endif

	return 0;
}
