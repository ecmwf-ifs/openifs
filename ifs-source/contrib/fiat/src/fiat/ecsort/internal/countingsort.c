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

/* ec_countingsort_() : Fortran-callable counting-sort */

/*
   by Sami Saarinen, ECMWF, 01/11/2007 : 1st working version

   Algorithm derived from C++ implementation in Wikipedia
   However, the index[]'ed version added by ourselves.
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

#define SORT_UINT 0

#define SORT_INT  1

#define SORT_R64  2

#define SORT_R32  3

#define SORT_I64  4

#define SORT_U64  5

typedef unsigned long long int  Uint64;
typedef          long long int  Sint64;
typedef unsigned int            Uint32;
typedef          int            Sint32;

typedef          long long int   ll_t;
typedef unsigned long long int u_ll_t;

#define  ALLOC(x,size)    \
 { ll_t bytes = (ll_t)sizeof(*x) * (size); \
   bytes = (bytes < 1) ? 1 : bytes; \
   x = malloc(bytes); \
   if (!x) { fprintf(stderr, \
                     "malloc() of %s (%lld bytes) failed in file=%s, line=%d\n", \
                     #x, bytes, __FILE__, __LINE__); RAISE(SIGABRT); } }

#define  CALLOC(x,size)    \
 { ll_t sz = (ll_t)(size); \
   sz = (sz < 1) ? 1 : sz; \
   x = calloc(sz, sizeof(*x)); \
   if (!x) { ll_t bytes = (ll_t)sizeof(*x) * (sz); \
             fprintf(stderr, \
                     "calloc() of %s (%lld bytes) failed in file=%s, line=%d\n", \
                     #x, bytes, __FILE__, __LINE__); RAISE(SIGABRT); } }

#define FREE(x)           if (x) { free(x); x = NULL; }

#ifdef DEBUG
#define AZZERT(cond) if (cond) ABOR1("Azzertion failed: "#cond)
#else
#define AZZERT(cond) 
#endif

/* Applicable for 32-bits only */

#define SIGNBIT32   0x80000000u
#define MASKALL32   0xFFFFFFFFu

#define CVMGM32(a,b,c) ( ((c) & SIGNBIT32) ? (a) : (b) )

static const int Npasses32 = 2; /* i.e. 2 x 16-bit passes == 32-bits */

/* Applicable for 64-bits only */

#define SIGNBIT64   0x8000000000000000ull
#define MASKALL64   0xFFFFFFFFFFFFFFFFull

#define CVMGM64(a,b,c) ( ((c) & SIGNBIT64) ? (a) : (b) )

static const int Npasses64 = 4; /* i.e. 4 x 16-bit passes == 64-bits */

/* CountingSort */

#define MASKALL16 0xFFFF

#define NCOUNT (MASKALL16+1)

typedef struct {
  int *sorted;
  int counts[NCOUNT];
} cs_shared_t;

#define FOR(i,s,e) for (i = s; i < e; ++i)

#define SHIFTMASK(a,shift) ((int)(((a) >> shift) & mask))

#define CntSort(T,NB,shift,idummy) \
static void \
CSortSM##shift##NB(const T A[], const int n, int local_index[], \
                   const T idummy, cs_shared_t *cs, const int irev) \
{ /* Note: A[] is a contiguous, stride=1, local data -- a shade copy of the original "void *Data"-array */ \
  const T mask = MASKALL16; \
  int i, tmp, nr, min = MASKALL16, max = min; \
  FOR(i,0,n) { \
    tmp = SHIFTMASK(A[i],shift); if (tmp < min) min = tmp; else if (tmp > max) max = tmp; \
  } \
  nr = max - min + 1; \
  AZZERT(nr <= 0 || nr > NCOUNT); \
  if (nr > 1) { /* i.e. max > min */ \
    /* nr == 1 would have meant that all values were equal --> skip */ \
    int j, icnt; \
    int *counts = cs->counts; \
    memset(counts, 0, nr * sizeof(*counts)); \
    if (irev) { /* Reverse, descending order */ \
      FOR(i,0,n) { tmp = max - SHIFTMASK(A[i],shift); AZZERT(tmp < 0 || tmp >= NCOUNT); ++counts[ tmp ]; } \
    } \
    else { /* Ascending order */ \
      FOR(i,0,n) { tmp = SHIFTMASK(A[i],shift) - min; AZZERT(tmp < 0 || tmp >= NCOUNT); ++counts[ tmp ]; } \
    } \
    /* Cascade counts to get cumulative counts */ \
    FOR(j,1,nr) counts[j] += counts[j-1]; \
    { \
      int *sorted = cs->sorted; \
      if (!sorted) { ALLOC(sorted, n); cs->sorted = sorted; } \
      if (irev) { /* Reverse, descending order */ \
	for (i = n-1; i >= 0; --i) { \
          j = local_index[i]; \
	  tmp = max - SHIFTMASK(A[j],shift); AZZERT(tmp < 0 || tmp >= NCOUNT); \
          icnt = --counts[ tmp ]; AZZERT(icnt < 0 || icnt >= n); \
	  sorted[icnt] = j; \
	} \
      } \
      else { /* Ascending order */ \
	for (i = n-1; i >= 0; --i) { \
          j = local_index[i]; \
	  tmp = SHIFTMASK(A[j],shift) - min; AZZERT(tmp < 0 || tmp >= NCOUNT); \
          icnt = --counts[ tmp ]; AZZERT(icnt < 0 || icnt >= n); \
	  sorted[icnt] = j; \
	} \
      } \
      memcpy(local_index, sorted, n * sizeof(int)); \
    } \
  } /* if (nr > 1) */ \
}

#define Helpers(T,NB) \
static T *  \
signmask##NB(const T Data[], int n, int inc, const int *index, int index_adj) \
{ \
  T *A = NULL; \
  int i; \
  ALLOC(A, n); \
  if (index && index_adj == 0) { \
    if (inc == 1) { \
      FOR(i,0,n) { \
        int j = index[i]; \
        T mask = CVMGM##NB(MASKALL##NB, SIGNBIT##NB, Data[j]); \
        A[i] = Data[j] ^ mask; \
      } \
    } else { \
      FOR(i,0,n) { \
        int j = index[i]*inc; \
        T mask = CVMGM##NB(MASKALL##NB, SIGNBIT##NB, Data[j]); \
        A[i] = Data[j] ^ mask; \
      } \
    } \
  } \
  else if (index) { \
    if (inc == 1) { \
      FOR(i,0,n) { \
        int j = (index[i] - index_adj); \
        T mask = CVMGM##NB(MASKALL##NB, SIGNBIT##NB, Data[j]); \
        A[i] = Data[j] ^ mask; \
      } \
    } else { \
      FOR(i,0,n) { \
        int j = (index[i] - index_adj)*inc; \
        T mask = CVMGM##NB(MASKALL##NB, SIGNBIT##NB, Data[j]); \
        A[i] = Data[j] ^ mask; \
      } \
    } \
  } else { \
    if (inc == 1) { \
      FOR(i,0,n) { \
        T mask = CVMGM##NB(MASKALL##NB, SIGNBIT##NB, Data[i]); \
        A[i] = Data[i] ^ mask; \
      } \
    } else { \
      FOR(i,0,n) { \
        int j = i*inc; \
        T mask = CVMGM##NB(MASKALL##NB, SIGNBIT##NB, Data[j]); \
        A[i] = Data[j] ^ mask; \
      } \
    } \
  } \
  return A; \
} \
static T * \
justcopy##NB(const T Data[], int n, int inc, const int *index, int index_adj) \
{ \
  T *A = NULL; \
  int i; \
  ALLOC(A, n); \
  if (index && index_adj == 0) { \
    if (inc == 1) { \
      FOR(i,0,n) A[i] = Data[index[i]]; \
    } else { \
      FOR(i,0,n) A[i] = Data[index[i]*inc]; \
    } \
  } \
  else if (index) { \
    if (inc == 1) { \
      FOR(i,0,n) A[i] = Data[(index[i]-index_adj)]; \
    } else { \
      FOR(i,0,n) A[i] = Data[(index[i]-index_adj)*inc]; \
    } \
  } \
  else { \
    if (inc == 1) { \
      memcpy(A, Data, n * sizeof(T)); \
    } else { \
      FOR(i,0,n) A[i] = Data[i*inc]; \
    } \
  } \
  return A; \
} \
static void \
sorted##NB(T Data[], int n, int inc, const int local_index[], T work[]) \
{ \
  int i; \
  if (inc == 1) { \
    FOR(i,0,n) work[i] = Data[local_index[i]]; \
    memcpy(Data, work, n * sizeof(T)); \
  } \
  else { \
    FOR(i,0,n) work[i] = Data[local_index[i]*inc]; \
    FOR(i,0,n) Data[i * inc] = work[i]; \
  } \
}


static int *
CreateIndex(int n)
{
  int *local_index = NULL;
  int i;
  ALLOC(local_index, n);
  FOR(i,0,n) local_index[i] = i;
  return local_index;
}

static void
Local2GlobalIndex(int n, const int local_index[], int index[], void *work)
{
  int i;
  int *tmpidx = work;
  FOR(i,0,n) {
    int lc = local_index[i];
    AZZERT(lc < 0 || lc >= n);
    tmpidx[i] = index[lc];
  }
  memcpy(index, tmpidx, n * sizeof(*index));
}


CntSort(Uint32,32,0,idummy)
CntSort(Uint32,32,shift,shift)
Helpers(Uint32,32)


CntSort(Uint64,64,0,idummy)
CntSort(Uint64,64,shift,shift)
Helpers(Uint64,64)


#define DoSort(T,copyfun,NB) { \
  int j; \
  int Npasses = Npasses##NB; \
  int *local_index = NULL; \
  T *data = Data; \
  T *dada = copyfun(&data[addr], n, inc, index, index ? index_adj : 0); \
  T shift = 0; \
  cs_shared_t cs; \
  cs.sorted = NULL; \
  local_index = CreateIndex(n); \
  CSortSM0##NB(dada, n, local_index, shift, &cs, irev); \
  for (j=1; j<Npasses; j++) { \
    shift += 16; \
    CSortSMshift##NB(dada, n, local_index, shift, &cs, irev); \
  } \
  if (index) { \
    Local2GlobalIndex(n, local_index, index, dada); \
  } \
  else { /* No index[] supplied */ \
    sorted##NB(&data[addr], n, inc, local_index, dada); \
  } \
  FREE(local_index); \
  FREE(cs.sorted); \
  FREE(dada); \
  rc = n; \
}

void
ec_countingsort_(const    int *Mode,
		 const    int *N,
		 const    int *Inc,
		 const    int *Start_addr,
            	         void *Data,
		          int *index,
		 const    int *Nindex,
		 const    int *Index_adj,
		 const    int *Reverse,
	                  int *retc)
{
  int mode = *Mode;
  int method = mode%10;
  int n = *N;
  int rc = n;
  int inc = *Inc;
  int addr = (*Start_addr) - 1; /* Fortran to C */
  int nidx = *Nindex; /* Must be >= n or otherwise the index[] is disregarded */
  int index_adj = *Index_adj;
  int irev = *Reverse;

  if (n <= 0) goto finish;

  if (nidx < n) index = NULL;
  
  switch (method) {
  case SORT_UINT:
    DoSort(Uint32,justcopy32,32);
    break;
  case SORT_INT:
    DoSort(Uint32,signmask32,32);
    break;
  case SORT_R64:
    DoSort(Uint64,signmask64,64);
    break;
  case SORT_R32:
    DoSort(Uint32,signmask32,32);
    break;
  case SORT_I64:
    DoSort(Uint64,signmask64,64);
    break;
  case SORT_U64:
    DoSort(Uint64,justcopy64,64);
    break;
  default:
    rc = -1;
    break;
  }

 finish:
  *retc = rc;
}

