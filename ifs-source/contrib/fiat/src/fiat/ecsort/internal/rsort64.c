/*
 * (C) Copyright 2005- ECMWF.
 * (C) Copyright 2013- Meteo-France.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <signal.h>
#include "raise.h"

/* rsort64_() : 64-bit Fortran-callable RADIX-sort */

/* 
   by Sami Saarinen, ECMWF, 22/02/2005 : Initial version derived from rsort32.c
           - " -            23/02/2005 : Fixes and some optimizations
           - " -            07/07/2005 : Mods in index_adj & bitsum
           - " -            07/02/2007 : Intercepting alloc (IBM & NEC SX) + NEC SX vectorization
      

   Thanks to Mike Fisher, ECMWF
   and Cray SCILIB ORDERS()-function developers
*/

/* 
   Methods:

   2 :          64-bit doubles (IEEE) : signbit + 11-bit exp + 52-bits mantissa
   4 :   Signed 64-bit ints
   5 : Unsigned 64-bit ints

*/

typedef unsigned int            Uint32;
typedef unsigned long long int  Uint64;
typedef unsigned char           Uchar;

/* scalar prozezzorz */
static int SpeedUp = 1;

#define SORT_R64  2
#define SORT_I64  4
#define SORT_U64  5

typedef long long int ll_t;

#define  ALLOC(x,size)    \
 { ll_t bytes = (ll_t)sizeof(*x) * (size); \
   bytes = (bytes < 1) ? 1 : bytes; \
   x = malloc(bytes); \
   if (!x) { fprintf(stderr, \
		     "malloc() of %s (%lld bytes) failed in file=%s, line=%d\n", \
		     #x, bytes, __FILE__, __LINE__); RAISE(SIGABRT); } }

#define FREE(x)           if (x) { free(x); x = NULL; }

#define BITSUM(x) bitsum[x] += ((item >> x) & 1ull)

#define SIGNBIT64   0x8000000000000000ull
#define MASKALL64   0xFFFFFFFFFFFFFFFFull
#define ZEROALL64   0x0000000000000000ull

#define CVMGM(a,b,c) ( ((c) & SIGNBIT64) ? (a) : (b) )

#define N64BITS 64

void 
rsort64_(const    int *Mode,       /* if < 10, then index[] needs to be initialized ; method = modulo 10 */
	 const    int *N,          /* no. of 64-bit elements */
	 const    int *Inc,        /* stride in terms of 64-bit elements */
	 const    int *Start_addr, /* Fortran start address i.e. normally == 1 */
	       Uint64  Data[],     /* 64-bit elements to be sorted */
	          int  index[],    /* sorting index */
	 const    int *Index_adj,  /* 0=index[] is a C-index, 1=index[] is a Fortran-index (the usual case) */
	          int *retc)
{
  int mode = *Mode;
  int method = mode%10;
  int n = *N;
  int rc = n;
  int inc = *Inc;
  int index_adj = *Index_adj;
  int addr = (*Start_addr) - 1; /* Fortran to C */
  int i, j, jj;
  Uchar xorit = 0;
  Uchar copytmp = 0;
  Uchar alloc_data = 0;
  Uint64 *data = NULL;
  int *tmp = NULL;
  Uint32 bitsum[N64BITS];

  if (method != SORT_R64    &&
      method != SORT_I64    &&
      method != SORT_U64 ) {
    rc = -1;
    goto finish;
  }

  if (n <= 0) {
    if (n < 0) rc = -2;
    goto finish;
  }

  if (inc < 1) {
    rc = -3;
    goto finish;
  }

  if (mode < 10) {
    /* index[] needs to be initialized */
    for (i=0; i<n; i++) index[i] = i + index_adj; /* Fortran-index */
  }

  j = addr;
  data = &Data[j];

  alloc_data = ((inc > 1) || (method == SORT_R64));
  if (alloc_data) ALLOC(data, n);

  if (method == SORT_R64) {
    for (i=0; i<n; i++) {
      Uint64 mask = CVMGM(MASKALL64, SIGNBIT64, Data[j]);
      data[i] = Data[j] ^ mask;
      j += inc;
    }
  }
  else if (method == SORT_I64) {
    if (inc == 1) { /* optimization */
      for (i=0; i<n; i++) {
	data[i] ^= SIGNBIT64;
      }
    }
    else {
      for (i=0; i<n; i++) {
	data[i] = Data[j] ^ SIGNBIT64;
	j += inc;
      }
    }
    xorit = 1;
  }
  else if (inc > 1) {
    for (i=0; i<n; i++) {
      data[i] = Data[j];
      j += inc;
    }
  }

  /* Check whether particular "bit-columns" are all zero or one */

  for (j=0; j<N64BITS; j++) bitsum[j] = 0;

  for (i=0; i<n; i++) {
    Uint64 item = data[i];
    /* Unrolled, full vector */
    BITSUM(0) ; BITSUM(1) ; BITSUM(2) ; BITSUM(3) ;
    BITSUM(4) ; BITSUM(5) ; BITSUM(6) ; BITSUM(7) ;
    BITSUM(8) ; BITSUM(9) ; BITSUM(10); BITSUM(11);
    BITSUM(12); BITSUM(13); BITSUM(14); BITSUM(15);
    BITSUM(16); BITSUM(17); BITSUM(18); BITSUM(19);
    BITSUM(20); BITSUM(21); BITSUM(22); BITSUM(23);
    BITSUM(24); BITSUM(25); BITSUM(26); BITSUM(27);
    BITSUM(28); BITSUM(29); BITSUM(30); BITSUM(31);
    BITSUM(32); BITSUM(33); BITSUM(34); BITSUM(35); 
    BITSUM(36); BITSUM(37); BITSUM(38); BITSUM(39); 
    BITSUM(40); BITSUM(41); BITSUM(42); BITSUM(43);
    BITSUM(44); BITSUM(45); BITSUM(46); BITSUM(47);
    BITSUM(48); BITSUM(49); BITSUM(50); BITSUM(51);
    BITSUM(52); BITSUM(53); BITSUM(54); BITSUM(55);
    BITSUM(56); BITSUM(57); BITSUM(58); BITSUM(59);
    BITSUM(60); BITSUM(61); BITSUM(62); BITSUM(63);
  }

  ALLOC(tmp, n);

  jj = 0;
  for (j=0; j<N64BITS; j++) {
    int sum = bitsum[j];
    if (sum > 0 && sum < n) { /* if 0 or n, then the whole column of bits#j 0's or 1's */
      Uint64 mask = (1ull << j);
      int *i1, *i2;
      
      if (jj%2 == 0) {
	i1 = index;
	i2 = tmp;
	copytmp = 1;
      }
      else {
	i1 = tmp;
	i2 = index;
	copytmp = 0;
      }
      
      if (SpeedUp == 0) {
	int k = 0;
	for (i=0; i<n; i++) /* Gather zero bits */
	  if ( (data[i1[i]-index_adj] & mask) ==    0 ) i2[k++] = i1[i];
	
	for (i=0; i<n; i++) /* Gather one bits */
	  if ( (data[i1[i]-index_adj] & mask) == mask ) i2[k++] = i1[i];
      }
      else
      {
	int k1 = 0, k2 = n-sum;
	for (i=0; i<n; i++) { /* Gather zero & one bits in a single sweep */
	  Uint64 value = data[i1[i]-index_adj] & mask;
	  i2[value == 0 ? k1++ : k2++] = i1[i];
	} /* for (i=0; i<n; i++) */
	if (k1 + sum != n || k2 != n) {
	  fprintf(stderr,
		  "***Programming error in rsort64_(): k1 + sum != n || k2 != n; k1=%d,k2=%d,sum=%d,n=%d\n",
		  k1,k2,sum,n);
	  RAISE(SIGABRT);
	}
      }
      
      jj++;
    } /* if (sum > 0 && sum < n) */
  }

  if (copytmp) for (i=0; i<n; i++) index[i] = tmp[i];

  FREE(tmp);

  if (!alloc_data && xorit && inc == 1) {
    /* 64-bit signed ints : backward */
    for (i=0; i<n; i++) data[i] ^= SIGNBIT64;
  }

  if (alloc_data) FREE(data);

 finish:

  *retc = rc;
}
