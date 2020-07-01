/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2020-04-22 17:02:44
  *Edit History: 
***********************************************************/

#ifndef XDK_SIMD_CHECK_H
#define XDK_SIMD_CHECK_H

#define CPUID_EAX 0
#define CPUID_EBX 1
#define CPUID_ECX 2
#define CPUID_EDX 3

#define OCR0_EAX  100

#define CPU_MMX_SUP     0x1U
#define CPU_SSE_SUP     0x2U
#define CPU_SSE2_SUP    0x4U
#define CPU_SSE3_SUP    0x8U
#define CPU_SSSE3_SUP   0x10U
#define CPU_FMA_SUP     0x20U
#define CPU_SSE4_1_SUP  0x40U
#define CPU_SSE4_2_SUP  0x80U
#define CPU_AES_SUP     0x100U
#define CPU_AVX_SUP     0x200U
#define CPU_XSAVE_SUP   0x400U
#define CPU_OSXSAVE_SUP 0x800U
#define CPU_AVX2_SUP    0x1000U
#define CPU_AVX512_SUP  0x2000U

#define OS_MMX_SUP      0x4000U

// XMM, SSE series
#define OS_SSE_SUP      0x8000U

// YMM, AVX/AVX2
#define OS_AVX_SUP      0x10000U

// ZMM0 - ZMM15, AVX512
// ZMM16 - ZMM31, AVX512
// OPMASK, AVX512
#define OS_AVX512_SUP   0x20000U

#ifdef __cplusplus
extern "C" {
#endif

	int simd_info_collect (void);

	int simd_check (int flag);

	int simd_info_dump (void);

#ifdef __cplusplus
}
#endif

#endif
