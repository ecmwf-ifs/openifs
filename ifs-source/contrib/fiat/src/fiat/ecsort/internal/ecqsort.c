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

/* ecqsort_() : Fortran-callable quick-sort */

/* 
   by Sami Saarinen, ECMWF, 7/07/2005 : Interface derived from rsort32.c & rsort64.c
           - " -            4/09/2006 : Dr.Hook call for kwiksort_u64_index
           - " -            7/02/2007 : Intercepting alloc (IBM & NEC SX) + NEC SX vectorization
           - " -            3/07/2007 : Rewritten to use qsort() standard library routine
           - " -           15/10/2007 : Fast qsort() added for simple 1-dim cases (see ../include/ecsort_shared.h)
           - " -           16/10/2007 : Reverse-flag added to avoid explicit negation of the array (cheaper)
           - " -           03/12/2007 : Dr.Hook calls removed; disturbs when in OMP-region on IBM ; a compiler bug ?
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

typedef unsigned long long int  Uint64;
typedef          long long int  Sint64;
typedef unsigned int            Uint32;
typedef          int            Sint32;

typedef    short int            Sint16;
typedef   signed char           Sint8;

typedef long long int ll_t;

#define  ALLOC(x,size)    \
 { ll_t bytes = (ll_t)sizeof(*x) * (size); \
   bytes = (bytes < 1) ? 1 : bytes; \
   x = malloc(bytes); \
   if (!x) { fprintf(stderr, \
		     "malloc() of %s (%lld bytes) failed in file=%s, line=%d\n", \
		     #x, bytes, __FILE__, __LINE__); RAISE(SIGABRT); } }

#define FREE(x)           if (x) { free(x); x = NULL; }

#if defined(NO_TRUNC)
/* For systems without trunc() -function [an extension of ANSI-C, but usually available] */
#define trunc(x) ((x) - fmod((x),1))
#else
extern double trunc(double d);
#endif

#define MakeKwikSort(T) \
typedef struct { \
  const T *valueptr; \
  int j; \
  int idx; \
} T##Str_t; \
\
static int \
T##cmp(const T##Str_t *a, const T##Str_t *b) { \
  if      ( *a->valueptr > *b->valueptr ) return  1; \
  else if ( *a->valueptr < *b->valueptr ) return -1; \
  else { /* *a->valueptr == *b->valueptr */ \
    /* the next line is essential for the stable qsort() */ \
    return (a->j > b->j) ? 1 : -1; \
  } \
} \
static int \
T##cmp_rev(const T##Str_t *a, const T##Str_t *b) { \
  if      ( *a->valueptr < *b->valueptr ) return  1; \
  else if ( *a->valueptr > *b->valueptr ) return -1; \
  else { /* a->valueptr == b->valueptr */ \
    /* the next line is essential for the stable qsort() */ \
    return (a->j > b->j) ? 1 : -1; \
  } \
} \
\
static void \
kwiksort_##T(const T v[], int n, int index[], int inc, \
             int index_adj, int mode, int irev) \
{ \
  int j; \
  T##Str_t *x = NULL; \
  ALLOC(x, n); \
  if (mode < 10) { \
    /* index[] needs to be initialized */ \
    if (inc == 1) { \
      for (j=0; j<n; j++) { x[j].valueptr = &v[j]; \
                            x[j].j = j; \
                            x[j].idx = j + index_adj; /* C -> Fortran */ \
                          } \
    } \
    else { \
      for (j=0; j<n; j++) { x[j].valueptr = &v[j * inc]; \
                            x[j].j = j; \
                            x[j].idx = j + index_adj; /* C -> Fortran */ \
                          } \
    } \
  } \
  else { \
    if (inc == 1) { \
      for (j=0; j<n; j++) { \
        int tmpidx = index[j] - index_adj; /* Fortran -> C */ \
        x[j].valueptr = &v[tmpidx]; \
        x[j].j = j; \
        x[j].idx = index[j]; \
      } \
    } \
    else { \
      for (j=0; j<n; j++) { \
        int tmpidx = index[j] - index_adj; /* Fortran -> C */ \
        x[j].valueptr = &v[tmpidx * inc]; \
        x[j].j = j; \
        x[j].idx = index[j]; \
      } \
    } \
  } \
  qsort(x, n, sizeof(*x), \
	irev ? \
	(int (*)(const void *, const void *))T##cmp_rev : \
	(int (*)(const void *, const void *))T##cmp); \
  for (j=0; j<n; j++) index[j] = x[j].idx; /* Re-arranged indices */ \
  FREE(x); \
}

#define MakeFastSort(T) \
static int \
T##fcmp(const T *a, const T *b) { \
  if      ( *a > *b ) return  1; \
  else if ( *a < *b ) return -1; \
  else return 0; \
} \
static int \
T##fcmp_rev(const T *a, const T *b) { \
  if      ( *a < *b ) return  1; \
  else if ( *a > *b ) return -1; \
  else return 0; \
} \
\
static void \
FastSort_##T(T v[], int n, int irev) \
{ \
  qsort(v, n, sizeof(*v), \
	irev ? \
	(int (*)(const void *, const void *))T##fcmp_rev : \
	(int (*)(const void *, const void *))T##fcmp); \
}

#define kwiksort(T) \
MakeKwikSort(T) \
MakeFastSort(T)

#define SORT_UINT 0
kwiksort(Uint32)

#define SORT_INT  1
kwiksort(Sint32)

#define SORT_R64  2
kwiksort(double)

#define SORT_R32  3
kwiksort(float)

#define SORT_I64  4
kwiksort(Sint64)

#define SORT_U64  5
kwiksort(Uint64)

#define DoSort(T) { \
  T *data = Data; \
  { \
    kwiksort_##T(&data[addr], n, index, inc, index_adj, mode, irev); \
  } \
}

#define DoFastSort(T) { \
  T *data = Data; \
  { \
    FastSort_##T(data, n, irev); \
  } \
}

void 
ecqsort_(const    int *Mode,
	 const    int *N,
	 const    int *Inc,
	 const    int *Start_addr,
	         void *Data,
	          int  index[],
	 const    int *Index_adj,
	 const    int *Reverse,
	          int *retc)
{
  int mode = *Mode;
  int method = mode%10;
  int n = *N;
  int rc = n;
  int inc = *Inc;
  int index_adj = *Index_adj;
  int irev = *Reverse;
  int addr = (*Start_addr) - 1; /* Fortran to C */

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

  switch (method) {
  case SORT_UINT:
    DoSort(Uint32);
    break;
  case SORT_INT:
    DoSort(Sint32);
    break;
  case SORT_R64:
    DoSort(double);
    break;
  case SORT_R32:
    DoSort(float);
    break;
  case SORT_I64:
    DoSort(Sint64);
    break;
  case SORT_U64:
    DoSort(Uint64);
    break;
  }

 finish:

  *retc = rc;
}

void 
ecqsortfast_(const    int *Mode,
	     const    int *N,
	             void *Data,
	     const    int *Reverse,
	              int *retc)
{
  int mode = *Mode;
  int method = mode%10;
  int n = *N;
  int rc = n;
  int irev = *Reverse;

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

  switch (method) {
  case SORT_UINT:
    DoFastSort(Uint32);
    break;
  case SORT_INT:
    DoFastSort(Sint32);
    break;
  case SORT_R64:
    DoFastSort(double);
    break;
  case SORT_R32:
    DoFastSort(float);
    break;
  case SORT_I64:
    DoFastSort(Sint64);
    break;
  case SORT_U64:
    DoFastSort(Uint64);
    break;
  }

 finish:

  *retc = rc;
}

#define MakeMergeFuncs(T)  \
static int \
T##_Merge(T data[], int amax, int bmax) \
{ \
  int i, j, N = amax + bmax; \
  T *a = data; \
  T *b = &data[amax]; \
  i=amax-1; j=N-1; \
  if (a[i] > b[0]) { \
    int k; \
    T *c = NULL; \
    ALLOC(c, bmax); \
    memcpy(c, b, bmax * sizeof(T)); \
    k=bmax-1; \
    while ((i >= 0) && (k >= 0)) { \
      if (a[i] >= c[k]) data[j--] = a[i--]; else data[j--] = c[k--]; \
    } \
    while (k >= 0) data[j--] = c[k--]; \
    FREE(c); \
  } \
  return N; \
} \
static int \
T##_MergeIdx(const T data[], int amax, int bmax, int index[], const int rank[]) \
{ \
  int i, j, N = amax + bmax; \
  int *a = index; \
  int *b = &index[amax]; \
  i=amax-1; j=N-1; \
  if (data[a[i]] > data[b[0]]) { \
    int k; \
    int *c = NULL; \
    ALLOC(c, bmax); \
    memcpy(c, b, bmax * sizeof(int)); \
    k=bmax-1; \
    while ((i >= 0) && (k >= 0)) { \
      T dai = data[a[i]]; \
      T dck = data[c[k]]; \
      if (dai > dck) index[j--] = a[i--]; \
      else if (dai < dck) index[j--] = c[k--]; \
      else { /* dai == dck : the rank[] decides */ \
        if (rank[a[i]] > rank[c[k]]) index[j--] = a[i--]; else index[j--] = c[k--]; \
      } \
    } \
    while (k >= 0) index[j--] = c[k--]; \
    FREE(c); \
  } \
  return N; \
} \
static int \
T##_Merge_rev(T data[], int amax, int bmax) \
{ \
  int i, j, N = amax + bmax; \
  T *a = data; \
  T *b = &data[amax]; \
  i=amax-1; j=N-1; \
  if (a[i] < b[0]) { \
    int k; \
    T *c = NULL; \
    ALLOC(c, bmax); \
    memcpy(c, b, bmax * sizeof(T)); \
    k=bmax-1; \
    while ((i >= 0) && (k >= 0)) { \
      if (a[i] <= c[k]) data[j--] = a[i--]; else data[j--] = c[k--]; \
    } \
    while (k >= 0) data[j--] = c[k--]; \
    FREE(c); \
  } \
  return N; \
} \
static int \
T##_MergeIdx_rev(const T data[], int amax, int bmax, int index[], const int rank[]) \
{ \
  int i, j, N = amax + bmax; \
  int *a = index; \
  int *b = &index[amax]; \
  i=amax-1; j=N-1; \
  if (data[a[i]] < data[b[0]]) { \
    int k; \
    int *c = NULL; \
    ALLOC(c, bmax); \
    memcpy(c, b, bmax * sizeof(int)); \
    k=bmax-1; \
    while ((i >= 0) && (k >= 0)) { \
      T dai = data[a[i]]; \
      T dck = data[c[k]]; \
      if (dai < dck) index[j--] = a[i--]; \
      else if (dai > dck) index[j--] = c[k--]; \
      else { /* dai == dck : the rank[] decides */ \
        if (rank[a[i]] > rank[c[k]]) index[j--] = a[i--]; else index[j--] = c[k--]; \
      } \
    } \
    while (k >= 0) index[j--] = c[k--]; \
    FREE(c); \
  } \
  return N; \
}

MakeMergeFuncs(Uint32)
MakeMergeFuncs(Sint32)
MakeMergeFuncs(double)
MakeMergeFuncs(float)
MakeMergeFuncs(Sint64)
MakeMergeFuncs(Uint64)

#define DoMerge(T) { \
  T *X = Data; \
  X += addr; \
  if (index && nidx >= n) { \
    int j; \
    if (index_adj) for (j=0; j<n; j++) index[j] -= index_adj; /* Fortran -> C */ \
    rc = irev ? \
      T##_MergeIdx_rev(X, amax, bmax, index, rank) : \
      T##_MergeIdx(X, amax, bmax, index, rank); \
    if (index_adj) for (j=0; j<n; j++) index[j] += index_adj; /* C -> Fortran */ \
  } \
  else { \
    rc = irev ? T##_Merge_rev(X, amax, bmax) : T##_Merge(X, amax, bmax); \
  } \
}

void ecmerge2_(const int *Mode, 
	       const int *Start_addr,
	       const int *Amax,
	       const int *Bmax,
	            void *Data,
	             int *index,
	       const int *Nidx,
	       const int *Index_adj,
	       const int *Reverse,
	       const int  rank[],
	             int *retc)
{
  int mode = *Mode;
  int method = mode%10;
  int addr = (*Start_addr) - 1; /* Fortran to C */
  int amax = *Amax;
  int bmax = *Bmax;
  int n = amax + bmax;
  int nidx = *Nidx;
  int index_adj = *Index_adj;
  int irev = *Reverse;
  int rc = 0;

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

  switch (method) {
  case SORT_UINT:
    DoMerge(Uint32);
    break;
  case SORT_INT:
    DoMerge(Sint32);
    break;
  case SORT_R64:
    DoMerge(double);
    break;
  case SORT_R32:
    DoMerge(float);
    break;
  case SORT_I64:
    DoMerge(Sint64);
    break;
  case SORT_U64:
    DoMerge(Uint64);
    break;
  }

 finish:
  *retc = rc;
}
