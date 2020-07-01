/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2020-04-22 17:41:22
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "simd_check.h"

#define N_SIMD_INFO 18

static uint32_t simd_info;

static const char * simd_mesg[N_SIMD_INFO] = {
											"CPU Support MMX",
											"CPU Support SSE",
											"CPU Support SSE2",
											"CPU Support SSE3",
											"CPU Support SSSE3",
											"CPU Support FMA",
											"CPU Support SSE4_1",
											"CPU Support SSE4_2",
											"CPU Support AES",
											"CPU Support AVX",
											"CPU Support XSAVE",
											"CPU Support OSXSAVE",
											"CPU Support AVX2",
											"CPU Support AVX512",

											"OS Support MMX",
											"OS Support SSE",
											"OS Support AVX/AVX2",
											"OS Support AVX512",
                     };

int
simd_info_collect (void)
{
	uint32_t a, b, c, d;
	uint32_t eax, ecx;

	simd_info = 0;

	// MMX
	eax = 1, ecx = 0;
	__asm__ ("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(eax), "c"(ecx));
	if (d & (1 << 23))
		simd_info |= CPU_MMX_SUP;

	// SSE
	eax = 1, ecx = 0;
	__asm__ ("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(eax), "c"(ecx));
	if (d & (1 << 25))
		simd_info |= CPU_SSE_SUP;

	// SSE2
	eax = 1, ecx = 0;
	__asm__ ("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(eax), "c"(ecx));
	if (d & (1 << 26))
		simd_info |= CPU_SSE2_SUP;

	// SSE3
	eax = 1, ecx = 0;
	__asm__ ("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(eax), "c"(ecx));
	if (c & (1 << 0))
		simd_info |= CPU_SSE3_SUP;

	// SSSE3
	eax = 1, ecx = 0;
	__asm__ ("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(eax), "c"(ecx));
	if (c & (1 << 9))
		simd_info |= CPU_SSSE3_SUP;

	// FMA
	eax = 1, ecx = 0;
	__asm__ ("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(eax), "c"(ecx));
	if (c & (1 << 12))
		simd_info |= CPU_FMA_SUP;

	// SSE4_1
	eax = 1, ecx = 0;
	__asm__ ("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(eax), "c"(ecx));
	if (c & (1 << 19))
		simd_info |= CPU_SSE4_1_SUP;

	// SSE4_2
	eax = 1, ecx = 0;
	__asm__ ("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(eax), "c"(ecx));
	if (c & (1 << 20))
		simd_info |= CPU_SSE4_2_SUP;

	// AES
	eax = 1, ecx = 0;
	__asm__ ("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(eax), "c"(ecx));
	if (c & (1 << 25))
		simd_info |= CPU_AES_SUP;

	// AVX
	eax = 1, ecx = 0;
	__asm__ ("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(eax), "c"(ecx));
	if (c & (1 << 28))
		simd_info |= CPU_AVX_SUP;

	// XSAVE
	eax = 1, ecx = 0;
	__asm__ ("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(eax), "c"(ecx));
	if (c & (1 << 26))
		simd_info |= CPU_XSAVE_SUP;

	// OSXSAVE
	eax = 1, ecx = 0;
	__asm__ ("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(eax), "c"(ecx));
	if (c & (1 << 27))
		simd_info |= CPU_OSXSAVE_SUP;

	// AVX2
	eax = 7, ecx = 0;
	__asm__ ("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(eax), "c"(ecx));
	if (b & (1 << 5))
		simd_info |= CPU_AVX2_SUP;

	// AVX512
	eax = 7, ecx = 0;
	__asm__ ("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(eax), "c"(ecx));
	if (b & (1 << 16))
		simd_info |= CPU_AVX512_SUP;

	// ------------------

	uint32_t xcr0_eax;

	__asm__ ("xgetbv" : "=a"(xcr0_eax) : "c"(0));

	// MMX
	if (xcr0_eax & (1 << 0))
		simd_info |= OS_MMX_SUP;

	// XMM
	if (xcr0_eax & (1 << 1))
		simd_info |= OS_SSE_SUP;

	// YMM
	if (xcr0_eax & (1 << 2))
		simd_info |= OS_AVX_SUP;

	// ZMM0-15, ZMM16-31, OPMASK
	if ((xcr0_eax & (1 << 6)) // ZMM0-15
			|| (xcr0_eax & (1 << 7)) // ZMM16-31
			|| (xcr0_eax & (1 << 5))) // OPMASK
		simd_info |= OS_AVX512_SUP;

	return 0;
}

int
simd_check (int flag)
{
	return (simd_info & flag) == flag;
}

int
simd_info_dump (void)
{
	int i;

	for (i=0; i<N_SIMD_INFO; ++i)
		if (simd_info & (1 << i))
			printf ("%s\n", simd_mesg[i]);

	return 0;
}
