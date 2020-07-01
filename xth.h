/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2019-03-24 22:12:48
  *Edit History: 
***********************************************************/

#ifndef XDK_XTH_H
#define XDK_XTH_H

#include <stdint.h>

typedef void* (*Worker) (void*);
typedef void (*DataInit3) (void**);
typedef void (*DataFree3) (void*);
typedef int (*PipelineWorker) (void*,void*,int);

struct xth_data_s;
typedef struct xth_data_s xth_data_t;

struct xth_engine_s;
typedef struct xth_engine_s xth_engine_t;

struct xth_data_s {
  void * data;
  int * signal;
  int thid;
  int nth;
};

#ifdef __cplusplus
extern "C" {
#endif

	int xth_parallel (Worker worker, void * data, int nth);
	
	xth_engine_t * xth_engine_start (Worker worker, void * data, int nth, int start_signal);
	void xth_engine_send_signal (xth_engine_t * engine, int work_signal, int wait_signal);
	void xth_engine_stop (xth_engine_t * engine, int stop_signal, int wait_signal);
	
	struct xth_cnter_s;
	typedef struct xth_cnter_s xth_cnter_t;
	
	xth_cnter_t * xth_counter_init (int64_t * cur_idx, int64_t step, const char * item_name);
	void xth_counter_start (xth_cnter_t * cnter);
	void xth_counter_reset (xth_cnter_t * cnter);
	void xth_counter_stop (xth_cnter_t * cnter);

#ifdef __cplusplus
}
#endif

#if 0
int xth_pipeline_init (int nth, int n_steps, void ** datas,
    DataInit3 data_init3, DataFree3 data_free3, PipelineWorker worker);
#endif

#endif
