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

/* rsort32_() : 32-bit Fortran-callable RADIX-sort */

/* 
   by Sami Saarinen, ECMWF, 3/2/1998 
         - " -              1/2/2000 : BIG_ENDIAN & LITTLE_ENDIAN labels renamed to *_INDIAN
                                       since they may conflict with the ones in <sys/endian.h>
         - " -              3/1/2001 : reference to valloc() removed; ALLOC() modified
	 - " -	           25/1/2001 : BIG_INDIAN removed (as label)
			               LITTLE_INDIAN called as LITTLE
         - " -            ??/9?/2001 : Speedup in rsort32
         - " -             14/3/2002 : rsort32_func implemeted to enable to run alternative sorting
                                       routine than rsort32
         - " -            18/02/2005 : Handle 64-bit (signed) ints
                                       IBM malloc() may call __malloc()/__free() [see getcurheap.c]
         - " -            21/02/2005 : Some optimization & endian detection on-the-fly
         - " -            07/07/2005 : Mods in index_adj & bitsum
	                               Added support for 64-bit unsigned ints
           - " -          07/02/2007 : Intercepting alloc (IBM & NEC SX) + NEC SX vectorization

   Thanks to Mike Fisher, ECMWF
   and Cray SCILIB ORDERS()-function developers
*/

/* 
   Methods:

   0 : Unsigned 32-bit ints
   1 :   Signed 32-bit ints
   2 :          64-bit doubles (IEEE) : signbit + 11-bit exp + 52-bits mantissa
   3 :          32-bit floats  (IEEE) : signbit +  8-bit exp + 23-bits mantissa
   4 :   Signed 64-bit ints
   5 : Unsigned 64-bit ints

*/

typedef unsigned int  Uint32;
typedef unsigned char Uchar;

/* scalar prozezzorz */
static int SpeedUp = 1;

#define SORT_UINT 0
#define SORT_INT  1
#define SORT_R64  2
#define SORT_R32  3
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

#define BITSUM(x) bitsum[x] += ((item >> x) & 1U)

#define SIGNBIT32   0x80000000
#define MASKALL32   0xFFFFFFFF
#define ZEROALL32   0x00000000

#define CVMGM(a,b,c) ( ((c) & SIGNBIT32) ? (a) : (b) )

#define N32BITS 32

void 
rsort32_(const    int *Mode,
	 const    int *N,
	 const    int *Inc,
	 const    int *Start_addr,
	       Uint32  Data[],
	          int  index[],
	 const    int *Index_adj,
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
  Uint32 *data = NULL;
  int *tmp = NULL;
  Uint32 bitsum[N32BITS];
  int lsw, msw;

  if (method != SORT_UINT   &&
      method != SORT_INT    &&
      method != SORT_R64    &&
      method != SORT_R32    &&
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

  { /* Little/big-endian selection */
    extern int ec_is_little_endian();
    int i_am_little = ec_is_little_endian();

    if (i_am_little) { 
      /* We are on little-endian machine */
      lsw   =  0;
      msw   =  1;
    }
    else { 
      /* We are on big-endian machine */
      lsw   =  1;
      msw   =  0;
    }
  }

  if (method == SORT_R64    || 
      method == SORT_I64    ||
      method == SORT_U64) {
    inc  *= 2;
    addr *= 2;
  }

  if (mode < 10) {
    /* index[] needs to be initialized */
    for (i=0; i<n; i++) index[i] = i + index_adj;
  }

  alloc_data = ((inc > 1) 
		|| (method == SORT_R32) 
		|| (method == SORT_R64)
		|| (method == SORT_I64)
		|| (method == SORT_U64)
		);
  if (alloc_data) ALLOC(data, n);

  if (method == SORT_R32) {
    j = addr;
    for (i=0; i<n; i++) {
      Uint32 mask = CVMGM(MASKALL32, SIGNBIT32, Data[j]);
      data[i] = Data[j] ^ mask;
      j += inc;
    }

    method = SORT_UINT;
  }
  else if (method == SORT_R64) {
    int Method = 10 + SORT_UINT;
    const int aStart_addr = 1;
    int aN = n;
    const int aInc = 1;
    int aIndex_adj = index_adj;

    /* Least significant word */
    j  = addr + lsw;
    jj = addr + msw;
    for (i=0; i<n; i++) {
      Uint32 mask = CVMGM(MASKALL32, ZEROALL32, Data[jj]);
      data[i] = Data[j] ^ mask;
      j += inc;
      jj += inc;
    }

    rsort32_(&Method, &aN, &aInc, &aStart_addr, data, index, &aIndex_adj, &rc);

    if (rc != n) goto finish;

    /* Most significant word */
    jj = addr + msw;
    for (i=0; i<n; i++) {
      Uint32 mask = CVMGM(MASKALL32, SIGNBIT32, Data[jj]);
      data[i] = Data[jj] ^ mask;
      jj += inc;
    }

    method = SORT_UINT;
  }
  else if (method == SORT_I64 || method == SORT_U64) {
    int Method = 10 + SORT_UINT;
    const int aStart_addr = 1;
    int aN = n;
    const int aInc = 1;
    int aIndex_adj = index_adj;

    /* Least significant word */
    jj = addr + lsw;
    for (i=0; i<n; i++) {
      data[i] = Data[jj];
      jj += inc;
    }

    rsort32_(&Method, &aN, &aInc, &aStart_addr, data, index, &aIndex_adj, &rc);

    if (rc != n) goto finish;

    /* Most significant word */
    if (method == SORT_I64) {
      jj = addr + msw;
      for (i=0; i<n; i++) {
	data[i] = Data[jj] ^ SIGNBIT32;
	jj += inc;
      }
    }
    else { /* unsigned 64-bit ints i.e. method == SORT_U64 */
      jj = addr + msw;
      for (i=0; i<n; i++) {
	data[i] = Data[jj];
	jj += inc;
      }
    }

    method = SORT_UINT;
  }
  else if (inc > 1) {
    j = addr;
    for (i=0; i<n; i++) {
      data[i] = Data[j];
      j += inc;
    }
  }
  else {
    data = &Data[addr];
  }

  xorit = (method == SORT_INT);

  /* Check whether particular "bit-columns" are all zero or one */

  for (j=0; j<N32BITS; j++) bitsum[j] = 0;

  for (i=0; i<n; i++) {
    Uint32 item;
    if (xorit) data[i] ^= SIGNBIT32;
    item = data[i];
    /* Unrolled, full vector */
    BITSUM(0) ; BITSUM(1) ; BITSUM(2) ; BITSUM(3) ;
    BITSUM(4) ; BITSUM(5) ; BITSUM(6) ; BITSUM(7) ;
    BITSUM(8) ; BITSUM(9) ; BITSUM(10); BITSUM(11);
    BITSUM(12); BITSUM(13); BITSUM(14); BITSUM(15);
    BITSUM(16); BITSUM(17); BITSUM(18); BITSUM(19);
    BITSUM(20); BITSUM(21); BITSUM(22); BITSUM(23);
    BITSUM(24); BITSUM(25); BITSUM(26); BITSUM(27);
    BITSUM(28); BITSUM(29); BITSUM(30); BITSUM(31);
  }

  ALLOC(tmp, n);

  jj = 0;
  for (j=0; j<N32BITS; j++) {
    int sum = bitsum[j];
    if (sum > 0 && sum < n) { /* if 0 or n, then the whole column of bits#j 0's or 1's */
      Uint32 mask = (1U << j);
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
	  Uint32 value = data[i1[i]-index_adj] & mask;
	  i2[value == 0 ? k1++ : k2++] = i1[i];
	} /* for (i=0; i<n; i++) */
	if (k1 + sum != n || k2 != n) {
	  fprintf(stderr,
		  "***Programming error in rsort32_(): k1 + sum != n || k2 != n; k1=%d,k2=%d,sum=%d,n=%d\n",
		  k1,k2,sum,n);
	  RAISE(SIGABRT);
	}
      }
      
      jj++;
    }
  }

  if (copytmp) {
    for (i=0; i<n; i++) index[i] = tmp[i];
  }

  FREE(tmp);

  if (!alloc_data && xorit && inc == 1) {
    /* 32-bit signed ints : backward */
    for (i=0; i<n; i++) data[i] ^= SIGNBIT32;
  }

  if (alloc_data) FREE(data);

 finish:

  *retc = rc;
}

static 
void (*default_rsort32_func)(const    int *Mode,
			     const    int *N,
			     const    int *Inc,
			     const    int *Start_addr,
			           Uint32  Data[],
			              int  index[],
			     const    int *Index_adj,
			              int *retc) =
     rsort32_
;

void 
rsort32_func_(const    int *Mode,
	      const    int *N,
	      const    int *Inc,
	      const    int *Start_addr,
	      Uint32        Data[],
	      int           index[],
	      const    int *Index_adj,
	      int          *retc)
{
  default_rsort32_func(Mode,
		       N,
		       Inc,
		       Start_addr,
		       Data,
		       index,
		       Index_adj,
		       retc);
}

void 
rsort32_setup_(void (*func)(const    int *Mode,
			    const    int *N,
			    const    int *Inc,
			    const    int *Start_addr,
			    Uint32        Data[],
			    int           index[],
			    const    int *Index_adj,
			    int          *retc),
	       int *speedup)
{
  if (func) default_rsort32_func = func;
  if (speedup) SpeedUp = *speedup;
}
