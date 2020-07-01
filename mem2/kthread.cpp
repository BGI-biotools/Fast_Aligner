/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Vasimuddin Md, Sanchit Misra, Intel Corporation, Heng Li.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
         Heng Li <hli@jimmy.harvard.edu>.
*****************************************************************************************/

#include "kthread.h"
// #include <omp.h>
#include <stdio.h>

extern uint64_t proc_freq, tprof[LIM_R][LIM_C];
extern int nthreads, affy[256];
int g_itr;

static inline long steal_work(kt_for_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	// return k >= t->n? -1 : k;
	return k*BATCH_SIZE >= t->n? -1 : k;
}

#if 0  // GSHARED
static void *ktf_worker(void *data)
{
	ktf_worker_t *w = (ktf_worker_t*)data;
	long i, val = 0;
	for (;;) {
		// i = __sync_fetch_and_add(&w->i, w->t->n_threads);
		i = __sync_fetch_and_add(&g_itr, 1);
		int st = i * BATCH_SIZE;
		if (st >= w->t->n) break;
		int ed = (i + 1) * BATCH_SIZE < w->t->n? (i + 1) * BATCH_SIZE : w->t->n;
		w->t->func(w->t->data, st, ed-st, w->i);
	}
	pthread_exit(0);
}

#else

static void *ktf_worker(void *data)
{
	ktf_worker_t *w = (ktf_worker_t*)data;
	long i, val = 0;
	int tid = w->i;
#if 0 && (__liunx__)
	fprintf(stderr, "i: %d, CPU: %d\n", tid , sched_getcpu());
#endif
	
	for (;;) {
		i = __sync_fetch_and_add(&w->i, w->t->n_threads);
		// i = __sync_fetch_and_add(&g_itr, 1);
		int st = i * BATCH_SIZE;
		if (st >= w->t->n) break;
		int ed = (i + 1) * BATCH_SIZE < w->t->n? (i + 1) * BATCH_SIZE : w->t->n;
		w->t->func(w->t->data, st, ed-st, tid);
	}

	while ((i = steal_work(w->t)) >= 0) {
		int st = i * BATCH_SIZE;
		int ed = (i + 1) * BATCH_SIZE < w->t->n? (i + 1) * BATCH_SIZE : w->t->n;
		w->t->func(w->t->data, st, ed-st, tid);
	}
	pthread_exit(0);
}
#endif

void kt_for(void (*func)(void*, int, int, int), void *data, int n)
{
	int i;
	kt_for_t t;
	pthread_t *tid;
	t.func = func, t.data = data, t.n_threads = nthreads, t.n = n;
	t.w = (ktf_worker_t*) malloc (nthreads * sizeof(ktf_worker_t));
	tid = (pthread_t*) malloc (nthreads * sizeof(pthread_t));
	for (i = 0; i < nthreads; ++i)
		t.w[i].t = &t, t.w[i].i = i;

	pthread_attr_t attr;
    pthread_attr_init(&attr);
	
	// printf("getcpu: %d\n", sched_getcpu());
	g_itr = 0;
	for (i = 0; i < nthreads; ++i) {
#if 0 && (__linux__)
		cpu_set_t cpus;
		CPU_ZERO(&cpus);
		// CPU_SET(i, &cpus);
		CPU_SET(affy[i], &cpus);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);	
		pthread_create(&tid[i], &attr, ktf_worker, &t.w[i]);
#else
		pthread_create(&tid[i], NULL, ktf_worker, &t.w[i]);
#endif
	}
	for (i = 0; i < nthreads; ++i) pthread_join(tid[i], 0);

    free(t.w);
	free(tid);
}
