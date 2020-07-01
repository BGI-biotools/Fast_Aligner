/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2020-04-02 15:35:15
  *Edit History: 
***********************************************************/

#ifndef SORT_AND_MKD_H
#define SORT_AND_MKD_H

#ifdef __cplusplus
extern "C" {
#endif

	int ebam_sort_and_markdup (samFile * out, bio_spl_t * spl);

	int ebam_sort (samFile * out, bio_spl_t * spl);

#ifdef __cplusplus
}
#endif

#endif
