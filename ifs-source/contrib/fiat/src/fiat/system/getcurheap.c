/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>

typedef  long long int  ll_t;

#define getcurheap           getcurheap_
#define getmaxcurheap        getmaxcurheap_
#define getcurheap_thread    getcurheap_thread_
#define getmaxcurheap_thread getmaxcurheap_thread_
#define getmaxloc            getmaxloc_
#define resetmaxloc          resetmaxloc_
#define profile_heap_get     profile_heap_get_

static ll_t maxloc = 0;    /* For stackcheck */
static ll_t begloc = 0;    /* For stackcheck */
static int heapcheck = 0;  /* Fro heapcheck */

extern ll_t gethwm_();

static ll_t maxcurheap = 0;

void
profile_heap_get(ll_t val[], 
		 const int *Nval, 
		 const int *Icase,
		 int *nret)
     /* Fortran callable */
{
  *nret = 0;
}


ll_t
getcurheap()
{
  // Cray linker: if you intend to link with -hstd_alloc and use Cray C compiler, then compile this file with -DSTD_ALLOC too
#if !defined(STD_ALLOC) && (defined(_CRAYC) || defined(USE_TCMALLOC))
  extern size_t get_tcmalloc_current_allocated_bytes_();
  return get_tcmalloc_current_allocated_bytes_();
#else
  ll_t rc = gethwm_();
  if (rc > maxcurheap) maxcurheap = rc;
  return rc;
#endif
}

ll_t
getcurheap_thread(const int *thread_id)
{
  return getcurheap();
}

/* Maximum (total) current (virtual mem) allocation encountered */

ll_t
getmaxcurheap()
{
  return 0;
}

/* Maximum (total) current (virtual mem) allocation encountered per thread */

ll_t
getmaxcurheap_thread(const int *thread_id) /* ***Note: YOMOML thread id */
{
  return 0;
}

ll_t
getmaxloc()
{
  ll_t z=maxloc-begloc;
  return z;
}

void
resetmaxloc()
{
  maxloc=0;
}

void

setheapcheck_()
{
  heapcheck=1;
}
