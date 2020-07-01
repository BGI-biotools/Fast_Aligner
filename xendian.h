/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-01-26 11:39:36
  *Edit History: 
***********************************************************/

#ifndef XDH_ENDIAN_H
#define XDH_ENDIAN_H

#include <stdint.h>

static inline int is_big_endian (void)
{
	int64_t n = 1;
	return !(*((char*)(&n)));
}

static inline uint16_t swap_endian_2 (uint16_t v)
{
	return (uint16_t) (((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}

static inline void * swap_endian_2p (void * x)
{
	*(uint16_t*) x = swap_endian_2 (*(uint16_t*)x);
	return x;
}

static inline uint32_t swap_endian_4 (uint32_t v)
{
	v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
	return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}

static inline void * swap_endian_4p (void * x)
{
	*(uint32_t*) x = swap_endian_4 (*(uint32_t*)x);
}

static inline uint64_t swap_endian_8 (uint64_t v)
{
	v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
	v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
	return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}

static inline void * swap_endian_8p (void * x)
{
	*(uint64_t*) x = swap_endian_8 (*(uint64_t*)x);
	return x;
}

#endif
