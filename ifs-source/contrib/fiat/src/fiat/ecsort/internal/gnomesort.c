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

/* 
   gnomesort.c : The easiest sort on Earth ? Yes, it is. 

   This is how a Dutch Garden Gnome sorts a line of flower pots. Basically,
   he looks at the flower pot next to him and the previous one; 
   if they are in the right order he steps one pot forward, otherwise 
   he swaps them and steps one pot backwards. 
   Boundary conditions: if there is no previous pot, he steps forwards; 
   if there is no pot next to him, he is done.

   See more Wikipedia & http://www.cs.vu.nl/~dick/gnomesort.html

   Author: Sami Saarinen, ECMWF, 12-Nov-2007
*/

/* Fortran callable: ecgnomesort_() */

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

#define GnomeSort(T) \
static void GnomeSort_##T(T a[], const int n, const int inc) \
{ \
  int i = 0; \
  if (inc == 1) { \
    while (i < n) { \
      if (i == 0 || a[i-1] <= a[i]) ++i; \
      else {T tmp = a[i]; a[i] = a[i-1]; a[--i] = tmp;} \
    } \
  } else { \
    while (i < n) { \
      if (i == 0 || a[(i-1)*inc] <= a[i*inc]) ++i; \
      else {T tmp = a[i*inc]; a[i*inc] = a[(i-1)*inc]; a[(--i)*inc] = tmp;} \
    } \
  } \
} \
static void GnomeSortIdx_##T(const T a[], const int n, const int inc, int index[], const int index_adj) \
{ \
  int i = 0; \
  if (inc == 1) { \
    if (index_adj) { \
      while (i < n) { \
        if (i == 0 || a[index[i-1]-index_adj] <= a[index[i]-index_adj]) ++i; \
        else {int tmp = index[i]; index[i] = index[i-1]; index[--i] = tmp;} \
      } \
    } else { /* index_adj == 0 */ \
       while (i < n) { \
        if (i == 0 || a[index[i-1]] <= a[index[i]]) i++; \
        else {int tmp = index[i]; index[i] = index[i-1]; index[--i] = tmp;} \
      } \
    } \
  } else { /* inc != 1 */ \
    if (index_adj) { \
      while (i < n) { \
        if (i == 0 || a[(index[i-1]-index_adj)*inc] <= a[(index[i]-index_adj)*inc]) ++i; \
        else {int tmp = index[i]; index[i] = index[i-1]; index[--i] = tmp;} \
      } \
    } else { /* index_adj == 0 */ \
       while (i < n) { \
        if (i == 0 || a[index[i-1]*inc] <= a[index[i]*inc]) ++i; \
        else {int tmp = index[i]; index[i] = index[i-1]; index[--i] = tmp;} \
      } \
    } \
  } \
}

#define DoSort(T) { \
  T *data = Data; \
  if (index && nidx >= n) { \
    GnomeSortIdx_##T(&data[addr], n, inc, index, index_adj); \
  } else { \
    GnomeSort_##T(&data[addr], n, inc); \
  } \
}

#define SORT_UINT 0
GnomeSort(Uint32)

#define SORT_INT  1
GnomeSort(Sint32)

#define SORT_R64  2
GnomeSort(double)

#define SORT_R32  3
GnomeSort(float)

#define SORT_I64  4
GnomeSort(Sint64)

#define SORT_U64  5
GnomeSort(Uint64)

void
ecgnomesort_(const    int *Mode,
	     const    int *N,
	     const    int *Inc,
	     const    int *Start_addr,
	             void *Data,
                      int *index,
	     const    int *Nindex,
	     const    int *Index_adj,
	              int *retc)
{
  int mode = *Mode;
  int method = mode%10;
  int n = *N;
  int rc = n;
  int inc = *Inc;
  int nidx = *Nindex; /* Must be >= n or otherwise the index[] is disregarded */
  int index_adj = *Index_adj;
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
