/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2020-04-02 15:37:28
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mp.h"
#include "xth.h"
#include "var.h"
#include "ebam.h"
#include "sort.h"
#include "sample.h"
#include "sam_bio.h"
#include "sort_and_mark.h"

#define SIG_WAIT       0
#define SIG_WORK       1
#define SIG_STOP       2
#define SIG_LOAD_READS 3
#define SIG_SORT_READS 4
#define SIG_MARKDUP    5

static inline int
ebam_ptr_comp (ebam_t ** pa, ebam_t ** pb)
{
  if ((*pa)->key <= (*pb)->key)
    return 1;
  else
    return -1;
}

static inline int
ebam_ptr_comp4merge (ebam_t ** pa, ebam_t ** pb)
{
  if ((*pa)->key < (*pb)->key)
    return 1;
  else if ((*pa)->key > (*pb)->key)
    return -1;
	else if ((*pa)->lib_id <= (*pb)->lib_id)
		return 1;
	else
		return -1;
}

static inline uint64_t
ebam_ptr_get_key (ebam_t ** ptr)
{
	return (*ptr)->key;
}

XSORT_DEF (erp, ebam_t*, ebam_ptr_comp);
XSORT_LITE_DEF (erp_mrg, ebam_t*, ebam_ptr_comp4merge);

static inline int
erp_mkd_pair_cmp (ebam_t ** a, ebam_t ** b)
{
	int ret;
  ebam_t * pa = *a;
  ebam_t * pb = *b;

	ret = pa->end[0]<pb->end[0]
					? -1
					: (pa->end[0]>pb->end[0] ? 1 : 0);

	if (ret == 0)
		ret = pa->end[1]<pb->end[1]
						? -1
						: (pa->end[1]>pb->end[1] ? 1 : 0);

	if (ret == 0)
		ret = pa->pair_orit - pb->pair_orit;

	if (ret == 0) {
		int pa_is_sec = pa->flag & EBAM_IS_SEC;		
		int pb_is_sec = pb->flag & EBAM_IS_SEC;		
		ret = pa_is_sec - pb_is_sec;
	}

	if (ret == 0)
		ret = pa->tile - pb->tile;

	if (ret == 0)
		ret = pa->x - pb->x;

	if (ret == 0)
		ret = pa->y - pb->y;

	return ret;
}

static inline int
erp_mkd_frag_cmp (ebam_t ** a, ebam_t ** b)
{
	int ret;
	ebam_t * pa = *a;
	ebam_t * pb = *b;

	ret = pa->frag_end<pb->frag_end
					? -1
					: (pa->frag_end>pb->frag_end ? 1 : 0);

	if (ret == 0)
		ret = pa->frag_orit - pb->frag_orit;

	if (ret == 0)
		ret = pa->tile - pb->tile;

	if (ret == 0)
		ret = pa->x - pb->x;

	if (ret == 0)
		ret = pa->y - pb->y;

	return ret;
}

XSORT_DEF (erp_mkd_pair, ebam_t*, erp_mkd_pair_cmp);
XSORT_DEF (erp_mkd_frag, ebam_t*, erp_mkd_frag_cmp);

static ebam_t EMPTY_EBAM;

static void
empty_ebam_init (void)
{
	EMPTY_EBAM.key    = UINT64_MAX;
	EMPTY_EBAM.end[0] = INT64_MAX;
	EMPTY_EBAM.end[1] = INT64_MAX;
	EMPTY_EBAM.score  = -2;
}

typedef struct {
  int tid;
  int bin_idx;
  int lib_idx;
  int64_t bin_beg;
  int64_t bin_end;
  int64_t cur_idx;
  int32_t * pre_bat_left;
  int32_t * cur_bat_left;
	int32_t * n_reads;
  mp_t(erd) *** reads;

  // sort
  arr_t(erp) ** blk_rds;
  arr_t(erp) *** sort_rds;
	xsort_t(erp) ** sort_engine;

  // markdup
  arr_t(erp) ** mkd_rds;

  int32_t * rd_idx;
  int32_t * beg_idxs;
  int32_t * end_idxs;
  ebam_t ** heap_ptrs;
  ebam_t *** read_ptrs;
} shared_data_t;

static int pre_bat_idx;
static int cur_bat_idx;
static int32_t blk_size;

/*---------------------------------------------------------------------------*/
/*---------------------------- Static Functions -----------------------------*/
/*---------------------------------------------------------------------------*/

static shared_data_t *
shared_data_init (bio_spl_t * spl)
{
  int h, i, j;
  shared_data_t * sd;

  sd = (shared_data_t *) ckmalloc (sizeof(shared_data_t));
	sd->n_reads = (int32_t *) ckmalloc (n_libs * 4);;
  sd->reads = (mp_t(erd) ***) ckmalloc (2 * sizeof(mp_t(erd) **));
  sd->sort_rds = (arr_t(erp) ***) ckmalloc (2 * sizeof(arr_t(erp) **));

  for (h=0; h<2; ++h) {
    sd->reads[h] = (mp_t(erd) **) ckmalloc (n_libs * sizeof(mp_t(erd) *));
  	sd->sort_rds[h] = (arr_t(erp) **) ckmalloc (n_libs * sizeof(arr_t(erp) *));
    for (i=0; i<n_libs; ++i) {
      sd->reads[h][i] = mp_init (erd, ebam_init2);
			sd->sort_rds[h][i] = arr_init (erp);
		}
  }

  sd->blk_rds = (arr_t(erp) **) ckmalloc (n_threads * sizeof(arr_t(erp) *));
  sd->sort_engine = (xsort_t(erp) **) ckmalloc (n_threads * sizeof(xsort_t(erp) *));
  for (i=0; i<n_threads; ++i) {
    sd->blk_rds[i] = arr_init (erp);
	  sd->sort_engine[i] = xsort_init (erp, NULL, NULL);
  }

  sd->pre_bat_left = (int32_t *) ckmalloc (n_libs * 4);
  sd->cur_bat_left = (int32_t *) ckmalloc (n_libs * 4);
  sd->mkd_rds = (arr_t(erp) **) ckmalloc (n_libs * sizeof(arr_t(erp) *));
  for (i=0; i<n_libs; ++i)
    sd->mkd_rds[i] = arr_init (erp);

  sd->rd_idx = (int32_t *) ckmalloc (n_libs * sizeof(int32_t));
  sd->beg_idxs = (int32_t *) ckmalloc (n_libs * sizeof(int32_t));
  sd->end_idxs = (int32_t *) ckmalloc (n_libs * sizeof(int32_t));
  sd->heap_ptrs = (ebam_t **) ckmalloc (n_libs * sizeof(ebam_t *));
  sd->read_ptrs = (ebam_t ***) ckmalloc (n_libs * sizeof(ebam_t **));

	empty_ebam_init ();
  blk_size = BIN_SIZE / n_threads + 1;

  return sd;
}

static void
shared_data_resize (shared_data_t * sd, int32_t tid, int32_t bin_idx)
{
  int32_t i, j;

  for (i=0; i<n_libs; ++i) {
		sd->n_reads[i] = mapped_ebam_files[i][tid][bin_idx].n_reads;
		mp_resize (erd, sd->reads[cur_bat_idx][i], sd->n_reads[i]);
  }
}

static inline void
find_block_end (mp_t(erd) * sorted_reads, int32_t * tid, int32_t * pos)
{
  ebam_t * r = mp_last (erd, sorted_reads);
  *tid = r->b.core.tid;
  *pos = r->b.core.pos;
}

static int
mkd_pair_cmp (ebam_t * r1, ebam_t * r2)
{
	if (r1->end[0] != r2->end[0])
		return 1;

	if (r1->end[1] != r2->end[1])
		return 1;

	if (r1->pair_orit != r2->pair_orit)
		return 1;

	int r1_is_sec = r1->flag & EBAM_IS_SEC;
	int r2_is_sec = r2->flag & EBAM_IS_SEC;

	if (r1_is_sec != r2_is_sec)
		return 1;

	return 0;
}

static int
mkd_frag_cmp (ebam_t * r1, ebam_t * r2)
{
	if (r1->frag_end != r2->frag_end)
		return 1;

	if (r1->frag_orit != r2->frag_orit)
		return 1;

	return 0;
}

static inline void
mark_duplicate_core (ebam_t ** reads, int32_t bidx, int32_t eidx)
{
	int32_t i;
	int32_t max_idx;
	int32_t max_score;

	max_idx = bidx;
	max_score = reads[bidx]->score;
	for (i=bidx+1; i<eidx; ++i) {
		if (reads[i]->score > max_score) {
			max_idx = i;
			max_score = reads[i]->score;
		}
	}

	for (i=bidx; i<eidx; ++i)
		if (i != max_idx)
			reads[i]->b.core.flag |= 1024;
}

static void
read_pair_mkd (ebam_t ** reads, int32_t bidx, int32_t eidx)
{
	assert (bidx < eidx);
	if (eidx-bidx == 1)
		return;

	if (reads[bidx]->end[1] < 0)
		return;

	mark_duplicate_core (reads, bidx, eidx);
}

static void
read_frag_mkd (ebam_t ** reads, int32_t bidx, int32_t eidx, int has_pair, int has_frag)
{
	int32_t i;
	ebam_t * r;

	if (!has_frag)
		return;

	if (has_pair) {
		for (i=bidx; i<eidx; ++i)
			if (!(reads[i]->flag & EBAM_IS_PAIR))
				reads[i]->b.core.flag |= 1024;
		return;
	}

	mark_duplicate_core (reads, bidx, eidx);
}

static int
read_pairs_mkd (ebam_t ** reads, int32_t n_reads)
{
	int32_t i;
	int32_t beg;
	ebam_t * r;
	ebam_t * pr;

	beg = 0;
	pr  = reads[0];
	for (i=1; i<n_reads; ++i) {
		r = reads[i];
		if (mkd_pair_cmp(r,pr) != 0) {
			read_pair_mkd (reads, beg, i);
			beg = i;
			pr = r;
		}
	}
	read_pair_mkd (reads, beg, i);

	return 0;
}

static int
read_frags_mkd (ebam_t ** reads, int32_t n_reads)
{
	int has_pair;
	int has_frag;
	int32_t i;
	int32_t beg;
	ebam_t * r;
	ebam_t * pr;

	beg = 0;
	has_pair = 0;
	has_frag = 0;
	pr  = reads[0];

	if (pr->flag & EBAM_IS_PAIR)
		has_pair |= 1;
	else
		has_frag |= 1;

	for (i=1; i<n_reads; ++i) {
		r = reads[i];
		if (mkd_frag_cmp(r,pr) != 0) {
			read_frag_mkd (reads, beg, i, has_pair, has_frag);
			beg = i;
			pr = r;
			has_pair = 0;
			has_frag = 0;
		}

		if (r->flag & EBAM_IS_PAIR)
			has_pair |= 1;
		else
			has_frag |= 1;
	}
	read_frag_mkd (reads, beg, i, has_pair, has_frag);

	return 0;
}

static void
simple_reads_dump_core (samFile * out, ebam_t ** reads, int32_t beg_idx, int32_t end_idx)
{
	int ret;
	int64_t i;
	str_t * s;
	ebam_t * r;

	s = pth_str_buf;
  for (i=beg_idx; i<end_idx; ++i) {
    ret = sam_write1 (out, bam_hdr, &(reads[i]->b));

		if (!(reads[i]->b.core.flag & 256)
				&& !(reads[i]->b.core.flag & 2048))
			++mapped_reads;
	}
}

static void
unmapped_reads_dump (mp_t(erd) * reads, samFile * out)
{
	int ret;
	int64_t i;

	if ((reads->n=unmapped_ebam_file->n_reads) == 0)
		return;

	mp_resize (erd, reads, reads->n);
	load_reads (unmapped_ebam_file->path, reads->pool, reads->n);
	for (i=0; i<mp_cnt(reads); ++i)
		ret = sam_write1 (out, bam_hdr, &(reads->pool[i].b));
}

static void
simple_reads_dump (samFile * out,
		arr_t(erp) * pre_reads, int32_t pre_bat_left,
		arr_t(erp) * cur_reads, int32_t cur_bat_left)
{
  if (pre_bat_left >= 0)
    simple_reads_dump_core (out, pre_reads->arr, pre_bat_left, pre_reads->n);

  simple_reads_dump_core (out, cur_reads->arr, 0, cur_bat_left);
}

static void
merge_and_dump_core (samFile * out, arr_t(erp) ** reads,
    int32_t * idxs, int32_t * ends, int n_libs, ebam_t ** heap_ptrs, ebam_t *** read_ptrs)
{
	int ret;
	int n_done;
  int32_t i;
	int32_t lib_id;
	str_t * s;

	// n_done means the number of libs that have been addressed
	n_done = 0;

  for (i=0; i<n_libs; ++i) {
		if (idxs[i] >= ends[i]) {
			heap_ptrs[i] = &EMPTY_EBAM;
			heap_ptrs[i]->lib_id = i;
			++n_done;
			continue;
		}
    read_ptrs[i] = reads[i]->arr;
    heap_ptrs[i] = read_ptrs[i][idxs[i]++];
		heap_ptrs[i]->lib_id = i;
  }

	xheap_make_lite (erp_mrg, heap_ptrs, n_libs);

	s = pth_str_buf;
	while (n_done < n_libs) {
		// the first read is empty means all heap_ptrs point to empty ebam
		// so we need jump out of the loop
		if (heap_ptrs[0]->score == -2)
			break;
		ret = sam_write1 (out, bam_hdr, &heap_ptrs[0]->b);

		// load the next read
		lib_id = heap_ptrs[0]->lib_id;
		if (idxs[lib_id] >= ends[lib_id]) {
			heap_ptrs[0] = &EMPTY_EBAM;
			++n_done;
		} else
			heap_ptrs[0] = read_ptrs[lib_id][idxs[lib_id]++];
		heap_ptrs[0]->lib_id = lib_id;

		// heap adjust
		xheap_adjust_lite (erp_mrg, heap_ptrs, n_libs, 0);
	}
}

static int
merge_and_dump (samFile * out, shared_data_t * sd, int n_libs, int is_last_bin)
{
  int32_t i;
	arr_t(erp) ** pre_reads;
	arr_t(erp) ** cur_reads;

	pre_reads = sd->sort_rds[pre_bat_idx];
	cur_reads = sd->sort_rds[cur_bat_idx];

  if (sd->pre_bat_left[0] < 0) {
    for (i=1; i<n_libs; ++i)
      assert (sd->pre_bat_left[i] < 0);
  } else {
    for (i=0; i<n_libs; ++i) {
      assert (sd->pre_bat_left[i] >= 0);
      sd->beg_idxs[i] = sd->pre_bat_left[i];
      sd->end_idxs[i] = arr_cnt (pre_reads[i]);
    }
    merge_and_dump_core (out, pre_reads,
        sd->beg_idxs, sd->end_idxs, n_libs, sd->heap_ptrs, sd->read_ptrs);
  }

  for (i=0; i<n_libs; ++i) {
    sd->beg_idxs[i] = 0;
    if (is_last_bin)
      sd->end_idxs[i] = arr_cnt (cur_reads[i]);
    else
      sd->end_idxs[i] = sd->cur_bat_left[i];
  }
  merge_and_dump_core (out, cur_reads,
      sd->beg_idxs, sd->end_idxs, n_libs, sd->heap_ptrs, sd->read_ptrs);

  return 0;
}

static int
is_block_empty (shared_data_t * sd)
{
	int i;
	int64_t n_reads;
	mp_t(erd) * reads;

	n_reads = 0;
	for (i=0; i<n_libs; ++i) {
		reads = sd->reads[cur_bat_idx][i];
		n_reads += reads->n;
	}

	if (n_reads == 0)
		return 1;
	else
		return 0;
}

/*---------------------------------------------------------------------------*/
/*------------------------------ ESAM Worker --------------------------------*/
/*---------------------------------------------------------------------------*/

static void *
ebam_worker (void * data)
{
  char * file;
  int pid;
  int pth_idx;
	int blk_idx;
  int * signal;
  int32_t batch_end;
  int32_t pre_bat_left;
	int64_t i, j;
  int64_t n_reads;
	int64_t bin_beg;
	int64_t beg_idx;
  int64_t * cur_idx;
	int64_t * pth_cnt;
	ebam_t * r;
	arr_t(erp) * rd_ptrs;
	arr_t(erp) * mkd_rds;
	arr_t(erp) * sort_rds;
	arr_t(erp) * pre_sort_rds;
  mp_t(erd) * reads;
  xth_data_t * xd;
  shared_data_t * sd;
	xsort_t(erp_mkd_pair) * mkd_pair_sort;
	xsort_t(erp_mkd_frag) * mkd_frag_sort;

  xd = (xth_data_t *) data;
  pid = xd->thid;
  sd = (shared_data_t *) (xd->data);
  signal = xd->signal;
  cur_idx = &sd->cur_idx;
	pth_cnt = (int64_t *) ckmalloc (n_threads * sizeof(int64_t));
	mkd_pair_sort = xsort_init (erp_mkd_pair, NULL, NULL);
	mkd_frag_sort = xsort_init (erp_mkd_frag, NULL, NULL);

  for (;;) {
    if (*signal == SIG_LOAD_READS) {
      for (;;) {
        pth_idx = __sync_fetch_and_add (cur_idx, 1);
        if (pth_idx >= n_libs)
          break;

        reads = sd->reads[cur_bat_idx][pth_idx];
        n_reads = mapped_ebam_files[pth_idx][sd->tid][sd->bin_idx].n_reads;
        mp_clear (erd, reads, ebam_clear);

				sort_rds = sd->sort_rds[cur_bat_idx][pth_idx];
				arr_resize (erp, sort_rds, n_reads);
        reads->n = n_reads;
				sort_rds->n = n_reads;
        if (n_reads == 0)
          continue;

        file = mapped_ebam_files[pth_idx][sd->tid][sd->bin_idx].path;
        load_reads (file, reads->pool, n_reads);
      }

      *signal = SIG_WAIT;
    } else if (*signal == SIG_SORT_READS) {
			bin_beg = sd->bin_beg;
      reads = sd->reads[cur_bat_idx][sd->lib_idx];
			if ((n_reads = reads->n) == 0) {
				*signal = SIG_WAIT;
				continue;
			}

			rd_ptrs = sd->blk_rds[pid];
			arr_clear (erp, rd_ptrs);
			sort_rds = sd->sort_rds[cur_bat_idx][sd->lib_idx];
			memset (pth_cnt, 0, n_threads*sizeof(int64_t));

			for (i=0; i<mp_cnt(reads); ++i) {
				r = mp_at (erd, reads, i);
				blk_idx = (r->b.core.pos - bin_beg) / blk_size;
				assert (blk_idx>=0 && blk_idx<n_threads);
				++(pth_cnt[blk_idx]);
				if (blk_idx != pid)
					continue;
				arr_add (erp, rd_ptrs, r);
			}

			if (rd_ptrs->n == 0) {
				*signal = SIG_WAIT;
				continue;
			}

			assert (rd_ptrs->n == pth_cnt[pid]);
			xradix_sort (erp, sd->sort_engine[pid], rd_ptrs->arr, rd_ptrs->n, ebam_ptr_get_key);

			for (beg_idx=i=0; i<pid; ++i)
				beg_idx += pth_cnt[i];

			for (i=beg_idx,j=0; j<rd_ptrs->n; ++i,++j) {
				assert (i < n_reads);
				sort_rds->arr[i] = arr_at (erp, rd_ptrs, j);
			}

      *signal = SIG_WAIT;
    } else if (*signal == SIG_MARKDUP) {
			for (;;) {
				pth_idx = __sync_fetch_and_add (cur_idx, 1);
				if (pth_idx >= n_libs)
					break;

				sort_rds = sd->sort_rds[cur_bat_idx][pth_idx];
				pre_sort_rds = sd->sort_rds[pre_bat_idx][pth_idx];
				if (sort_rds->n == 0)
					continue;

        // find batch dump range
        batch_end = sd->bin_end - 100;
        for (i=sort_rds->n-1; i>=0; --i)
          if (sort_rds->arr[i]->b.core.pos < batch_end)
            break;
        sd->cur_bat_left[pth_idx] = i + 1;

				mkd_rds = sd->mkd_rds[pth_idx];
				arr_clear (erp, mkd_rds);

        pre_bat_left = sd->pre_bat_left[pth_idx];
				if (sd->bin_idx == 0)
					assert (pre_bat_left == -1);

        if (pre_bat_left >= 0) {
          for (i=pre_bat_left; i<arr_cnt(pre_sort_rds); ++i) {
            r = arr_at (erp, pre_sort_rds, i);
            if (r->end[0] < 0)
              continue;
            arr_add (erp, mkd_rds, r);
          }
        }

        for (i=0; i<arr_cnt(sort_rds); ++i) {
          r = arr_at (erp, sort_rds, i);
          if (r->end[0] < 0)
            continue;
          arr_add (erp, mkd_rds, r);
        }

				if (mkd_rds->n == 0)
					continue;

				// mark read pairs
        xmerge_sort (erp_mkd_pair, mkd_pair_sort, mkd_rds->arr, mkd_rds->n);
				read_pairs_mkd (mkd_rds->arr, mkd_rds->n);

				// mark read fragments
				xmerge_sort (erp_mkd_frag, mkd_frag_sort, mkd_rds->arr, mkd_rds->n);
				read_frags_mkd (mkd_rds->arr, mkd_rds->n);
			}

      *signal = SIG_WAIT;
    } else if (*signal == SIG_STOP) {
      *signal = SIG_WAIT;
      break;
    }

    usleep (1);
  }

	xsort_free (erp_mkd_pair, mkd_pair_sort, NULL);
	xsort_free (erp_mkd_frag, mkd_frag_sort, NULL);

  return (void *) 0;
}

int
ebam_sort_and_markdup (samFile * out, bio_spl_t * spl)
{
  int i, j, k;
	time_t beg;
	time_t beg2;
  xth_engine_t * engine;
  shared_data_t * sd;

  sd = shared_data_init (spl);
  engine = xth_engine_start (ebam_worker, (void*)sd, n_threads, SIG_WAIT);

  cur_bat_idx = 0;
	pre_bat_idx = 1;
  for (i=0; i<n_chr; ++i) {
    sd->tid = i;
		for (j=0; j<n_libs; ++j)
			sd->pre_bat_left[j] = -1;

    for (j=0; j<n_chr_bins[i]; ++j) {
			time (&beg);
      sd->bin_idx = j;
      sd->bin_beg = j << BIN_SHIFT;
      sd->bin_end = (j+1) << BIN_SHIFT;
      shared_data_resize (sd, i, j);

#if DEBUG
			printf ("\n");
			printf ("tid: %d, bin: %d, begins\n", i, j);
			time (&beg2);
#endif
      sd->cur_idx = 0;
      xth_engine_send_signal (engine, SIG_LOAD_READS, SIG_WAIT);
			if (is_block_empty(sd)) {
#if DEBUG
				printf ("tid: %d, bin: %d, total cost %lds\n", i, j, time(NULL)-beg);
#endif
				continue;
			}
#if DEBUG
			printf ("  load reads cost: %lds\n", time(NULL)-beg2);
#endif

#if DEBUG
			time (&beg2);
#endif
      for (k=0; k<n_libs; ++k) {
        sd->cur_idx = 0;
        sd->lib_idx = k;
        xth_engine_send_signal (engine, SIG_SORT_READS, SIG_WAIT);
      }
#if DEBUG
			printf ("  sort reads cost: %lds\n", time(NULL)-beg2);
#endif

#if DEBUG
			time (&beg2);
#endif
     	sd->cur_idx = 0;
     	xth_engine_send_signal (engine, SIG_MARKDUP, SIG_WAIT);
#if DEBUG
			printf ("  mark duplicates cost: %lds\n", time(NULL)-beg2);
#endif

#if DEBUG
			time (&beg2);
#endif
      if (n_libs == 1)
        simple_reads_dump (out, sd->sort_rds[pre_bat_idx][0], sd->pre_bat_left[0],
						sd->sort_rds[cur_bat_idx][0], sd->cur_bat_left[0]);
      else
        merge_and_dump (out, sd, n_libs, j==n_chr_bins[i]);
#if DEBUG
			printf ("  dump reads cost: %lds\n", time(NULL)-beg2);
#endif

      cur_bat_idx = 1 - cur_bat_idx;
			pre_bat_idx = 1 - pre_bat_idx;

      for (k=0; k<n_libs; ++k)
        sd->pre_bat_left[k] = sd->cur_bat_left[k];
#if DEBUG
			printf ("tid: %d, bin: %d, total cost %lds\n", i, j, time(NULL)-beg);
#endif
    }
  }

	unmapped_reads_dump (sd->reads[0][0], out);

#if DEBUG
	printf ("\n");
#endif

  xth_engine_stop (engine, SIG_STOP, SIG_WAIT);

  return 0;
}

int
ebam_sort (samFile * out, bio_spl_t * spl)
{
  int i, j, k;
	time_t beg;
	time_t beg2;
  xth_engine_t * engine;
  shared_data_t * sd;

  sd = shared_data_init (spl);
  engine = xth_engine_start (ebam_worker, (void*)sd, n_threads, SIG_WAIT);

	cur_bat_idx = 0;
	for (j=0; j<n_libs; ++j)
		sd->pre_bat_left[j] = -1;

  for (i=0; i<n_chr; ++i) {
    sd->tid = i;

    for (j=0; j<n_chr_bins[i]; ++j) {
			time (&beg);
      sd->bin_idx = j;
      sd->bin_beg = j << BIN_SHIFT;
      sd->bin_end = (j+1) << BIN_SHIFT;
      shared_data_resize (sd, i, j);

#if DEBUG
			printf ("\n");
			printf ("tid: %d, bin: %d, begins\n", i, j);
			time (&beg2);
#endif
      sd->cur_idx = 0;
      xth_engine_send_signal (engine, SIG_LOAD_READS, SIG_WAIT);
			if (is_block_empty(sd)) {
#if DEBUG
				printf ("tid: %d, bin: %d, total cost %lds\n", i, j, time(NULL)-beg);
#endif
				continue;
			}
#if DEBUG
			printf ("  load reads cost: %lds\n", time(NULL)-beg2);
#endif

#if DEBUG
			time (&beg2);
#endif
      for (k=0; k<n_libs; ++k) {
        sd->cur_idx = 0;
        sd->lib_idx = k;
        xth_engine_send_signal (engine, SIG_SORT_READS, SIG_WAIT);
      }
#if DEBUG
			printf ("  sort reads cost: %lds\n", time(NULL)-beg2);
#endif

#if DEBUG
			time (&beg2);
#endif
      if (n_libs == 1)
        simple_reads_dump (out, NULL, -1, sd->sort_rds[cur_bat_idx][0], mp_cnt(sd->sort_rds[cur_bat_idx][0]));
      else {
				for (k=0; k<n_libs; ++k)
					sd->cur_bat_left[k] = mp_cnt (sd->sort_rds[cur_bat_idx][k]);
        merge_and_dump (out, sd, n_libs, j==n_chr_bins[i]);
			}
#if DEBUG
			printf ("  dump reads cost: %lds\n", time(NULL)-beg2);
			printf ("tid: %d, bin: %d, total cost %lds\n", i, j, time(NULL)-beg);
#endif

    }
  }

	unmapped_reads_dump (sd->reads[0][0], out);

#if DEBUG
	printf ("\n");
#endif

  xth_engine_stop (engine, SIG_STOP, SIG_WAIT);

  return 0;
}
