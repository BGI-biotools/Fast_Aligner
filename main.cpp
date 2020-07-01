/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2019-12-27 11:35:12
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "simd_check.h"
#include "mem2/macro.h"
#include "mem2/profiling.h"

// ----------------------------------
uint64_t proc_freq, tprof[LIM_R][LIM_C], prof[LIM_R];
int nthreads;
int num_ranks = 1, myrank = 0;
int64_t reference_seq_len;
// ----------------------------------

static int
usage (void)
{
	fprintf (stderr, "\n");
	fprintf (stderr, "Program: Efficient Alignment Tool for Short Reads\n");
	fprintf (stderr, "Contact: chenxi1@genomics.cn\n");
	fprintf (stderr, "Version: v2.0.0\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "Usage:   fast_aligner <command> <arguments>\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "Command: index  prepare some permanent files for analysis\n");
	//fprintf (stderr, "         mem    filter-bwa_mem-sort-markdup pipeline\n");
	fprintf (stderr, "         mem2   filter-bwa_mem2-sort-markdup pipeline\n");
	fprintf (stderr, "\n");

	return 1;
}

extern int pipe_db_build_main (int argc, char * argv[]);
extern int mem_pipe_main (int argc, char * argv[]);
extern int wgs_pipe_main (int argc, char * argv[]);

int
main (int argc, char * argv[])
{
	if (argc < 2)
		return usage ();

	int ret;
	uint32_t flag_check;

	simd_info_collect ();

		fprintf(stderr, "-----------------------------\n");
#if __AVX512BW__
		fprintf(stderr, "Executing in AVX512 mode!!\n");
		flag_check = CPU_AVX512_SUP | OS_AVX512_SUP;
		if (!simd_check(flag_check))
			err_mesg ("AVX512 is not supported!\n");
#endif

#if ((!__AVX512BW__) & (__AVX2__))
		fprintf(stderr, "Executing in AVX2 mode!!\n");
		flag_check = CPU_AVX2_SUP | OS_AVX_SUP;
		if (!simd_check(flag_check))
			err_mesg ("AVX2 is not supported!\n");
#endif

#if ((!__AVX512BW__) && (!__AVX2__) && (__SSE2__))
		fprintf(stderr, "Executing in SSE4.1 mode!!\n");
#endif

#if ((!__AVX512BW__) && (!__AVX2__) && (!__SSE2__))
		fprintf(stderr, "Executing in Scalar mode!!\n");
#endif
		fprintf(stderr, "-----------------------------\n");

	if (strcmp(argv[1],"index") == 0) {
		ret = pipe_db_build_main (argc-1, argv+1);
	} else if (strcmp(argv[1],"mem") == 0) {
		ret = mem_pipe_main (argc-1, argv+1);
	} else if (strcmp(argv[1],"mem2") == 0) {
		ret = wgs_pipe_main (argc-1, argv+1);
	} else
		err_mesg ("invalid command!");

	return 0;
}
