/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2019-03-24 22:16:32
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xth.h"
#include "utils.h"

int
xth_parallel (Worker worker, void * data, int nth)
{
  int i;
  pthread_t * pids;
  xth_data_t * xdatas;

  pids = (pthread_t *) ckalloc (nth, sizeof(pthread_t));
  xdatas = (xth_data_t *) ckalloc (nth, sizeof(xth_data_t));
  for (i=0; i<nth; ++i) {
    xdatas[i].data = data;
    xdatas[i].thid = i;
    xdatas[i].nth = nth;
    ckpthread_create (pids+i, NULL, worker, (void*)(xdatas+i));
  }

  for (i=0; i<nth; ++i)
    ckpthread_join (pids[i]);

  free (pids);
  free (xdatas);

  return 0;
}

struct xth_engine_s {
  pthread_t * pids;
  int * signals;
  xth_data_t * datas;
  int nth;
};

xth_engine_t *
xth_engine_start (Worker worker, void * data, int nth, int start_signal)
{
  int i;
  xth_engine_t * engine;
  xth_data_t * xdatas;

  xdatas = (xth_data_t *) ckalloc (nth, sizeof(xth_data_t));
  engine = (xth_engine_t *) ckalloc (1, sizeof(xth_engine_t));
  engine->nth = nth;
  engine->datas = xdatas;
  engine->pids = (pthread_t *) ckalloc (nth, sizeof(pthread_t));
  engine->signals = (int *) ckalloc (nth, sizeof(int));
  for (i=0; i<nth; ++i) {
    engine->signals[i] = start_signal;
    xdatas[i].data = data;
    xdatas[i].thid = i;
    xdatas[i].nth = nth;
    xdatas[i].signal = engine->signals + i;
    ckpthread_create (engine->pids+i, NULL, worker, (void*)(xdatas+i));
  }

  return engine;
}

void
xth_engine_send_signal (xth_engine_t * engine, int work_signal, int wait_signal)
{
  send_work_signal (engine->signals, work_signal, wait_signal, engine->nth);
}

void
xth_engine_stop (xth_engine_t * engine, int stop_signal, int wait_signal)
{
  int i;

  xth_engine_send_signal (engine, stop_signal, wait_signal);
  for (i=0; i<engine->nth; ++i)
    ckpthread_join (engine->pids[i]);

  free (engine->datas);
  free (engine->pids);
  free (engine->signals);
  free (engine);
}

#define XTH_CNTER_WAIT  0
#define XTH_CNTER_WORK  1
#define XTH_CNTER_RESET 2
#define XTH_CNTER_STOP  3

struct xth_cnter_s {
  int64_t * cur_idx;
  int64_t step;
  char * item_name;
  pthread_t pid;
  xth_engine_t * engine;
};

static void *
counter (void * data)
{
  int * signal;
  int64_t pre_idx;
  xth_data_t * xd;
  xth_cnter_t * cnt;

  xd = (xth_data_t *) data;
  cnt = (xth_cnter_t *) xd->data;
  signal = xd->signal;

  for (;;) {
    if (*signal == XTH_CNTER_WORK) {
      pre_idx = *cnt->cur_idx;
      for (;;) {
        if (*cnt->cur_idx >= pre_idx+cnt->step) {
          fprintf (stderr, "[parallel counter] %ld %ss have been addressed\n", pre_idx+cnt->step, cnt->item_name);
          pre_idx += cnt->step;
        }

        if (*signal == XTH_CNTER_RESET) {
          fprintf (stderr, "[parallel counter] %ld %ss have been addressed\n", *cnt->cur_idx, cnt->item_name);
          break;
        }

        if (*signal == XTH_CNTER_STOP) {
          fprintf (stderr, "[parallel counter] %ld %ss have been addressed\n", *cnt->cur_idx, cnt->item_name);
          break;
        }
      }

      *signal = XTH_CNTER_WAIT;
    } else if (*signal == XTH_CNTER_RESET) {
      *signal = XTH_CNTER_WAIT;
    } else if (*signal == XTH_CNTER_STOP) {
      *signal = XTH_CNTER_WAIT;
      break;
    }

    usleep (100);
  }

  return (void *) 0;
}

xth_cnter_t *
xth_counter_init (int64_t * cur_idx, int64_t step, const char * item_name)
{
  xth_cnter_t * cnt;

  cnt = (xth_cnter_t *) ckalloc (1, sizeof(xth_cnter_t));
  cnt->cur_idx = cur_idx;
  cnt->step = step;
  cnt->item_name = strdup (item_name);
  cnt->engine = xth_engine_start (counter, (void*)cnt, 1, XTH_CNTER_WAIT);

  return cnt;
}

void
xth_counter_start (xth_cnter_t * cnter)
{
  cnter->engine->signals[0] = XTH_CNTER_WORK;
}

void
xth_counter_reset (xth_cnter_t * cnter)
{
  cnter->engine->signals[0] = XTH_CNTER_RESET;
}

void
xth_counter_stop (xth_cnter_t * cnter)
{
  xth_engine_stop (cnter->engine, XTH_CNTER_STOP, XTH_CNTER_WAIT);
  free (cnter);
}

#if 0
struct xth_pipeline_s;
typedef struct xth_pipeline_s xth_pipeline_t;

typedef struct {
  int data_idx;
  int step_idx;
  int pid;
  PipelineWorker worker;
  xth_pipeline_t * xp;
} xp_worker_t;

struct xth_pipeline_s {
  xp_worker_t * wokers;
  void * shared_data;
  void ** run_data;
  pthread_t * pids;
  int n_step;
  int n_data;
  int nth;
};

static void *
xp_core (void * data)
{
  int ret;
  xp_worker_t * w;
  xth_pipeline_t * xp;

  w = (xp_worker_t *) data;
  xp = w->xp;
  while (w->step_idx < xp->n_step) {
    pthread_mutex_lock (&xp->mutex);
    for (;;) {
      for (i=0; i<xp->nth; ++i) {
        if (i == w->pid)
          continue;
        if (xp->workers[i].step_idx <= w->step_idx
            && xp->workers[i].data_idx < w->data_idx)
          break;
      }
      if (i == xp->nth)
        break;
      pthread_cond_wait (&xp->cv, &p->mutex);
    }
    pthread_mutect_unlock (&xp->mutex);

    ret = w->worker (xp->shared_data, xp->run_data[w->data_idx], w->step_idx);

    pthread_mutex_lock (&xp->mutex);
  }
}

int
xth_pipeline (int nth, int n_steps, void * shared_data,
    DataInit3 data_init3, DataFree3 data_free3, PipelineWorker pipeline_worker)
{
  int i;
  pthread_t * pids;
  xp_worker_t * worker;
  xth_pipeline_t * xp;

  assert (data_init3 != NULL);
  assert (data_free3 != NULL);
  assert (shared_data != NULL);

  xp = (xth_pipeline_t *) ckalloc (1, sizeof(xth_pipeline_t));
  xp->n_step = n_step;
  xp->n_data = nth;
  xp->nth = nth;

  xp->shared_data = shared_data;
  xp->run_data = (void **) ckalloc (xp->n_data, sizeof(void *));
  for (i=0; i<xp->n_data; ++i)
    data_init3 (xp->run_data+i);

  xp->pids = (pthread_t *) ckalloc (xp->nth, sizeof(pthread_t));
  xp->workers = (xp_worker_t *) ckalloc (xp->nth, sizeof(xp_worker_t));
  for (i=0; i<xp->nth; ++i) {
    worker = xp->workers + i;
    worker->pid = i;
    worker->data_idx = 0;
    worker->step_idx = 0;
    worker->worker = pipeline_worker;
    ckpthread_create (xp->pids+i, NULL, xp_core, (void*)(worker));
  }

  for (i=0; i<xp->nth; ++i)
    ckpthread_join (xp->pids[i]);

  for (i=0; i<xp->n_data; ++i)
    data_free3 (xp->run_data[i]);
  free (xp->run_data);
  free (xp->pids);
  free (xp->workers);
  free (xp);

  return 0;
}
#endif
