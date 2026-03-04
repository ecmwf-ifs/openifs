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
#include <stdlib.h>
#include <unistd.h>

#include "getstatm.h"

typedef  long long int  ll_t;

static ll_t maxhwm = 0;

#define gethwm gethwm_

// Cray linker: if you intend to link with -hstd_alloc and use Cray C compiler, then compile this file with -DSTD_ALLOC too
#if !defined(STD_ALLOC) && (defined(_CRAYC) || defined(USE_TCMALLOC))
ll_t
gethwm()
{
  extern size_t get_tcmalloc_heap_size_();
  return get_tcmalloc_heap_size_();
}
#elif defined(LINUX)
static ll_t basesize = -1;
static size_t pagesize = 4096;
ll_t gethwm()
{
  struct statm sm;
  ll_t rc = 0;
  if (getstatm(&sm) == 0) {
    if (basesize < 0) { /* the very first time */
      basesize = sm.size;
      pagesize = getpagesize();
      if (pagesize <= 0) pagesize = 4096;
    }
    rc = (sm.size - basesize) * pagesize;
    if (rc > maxhwm) maxhwm = rc;
  }
  return rc;
}

#else
ll_t gethwm()
{
  ll_t rc = (ll_t)((char *)sbrk(0) - (char *)0);
  return rc;
}
#endif

ll_t getmaxhwm_()
{
  ll_t rc = gethwm_();
  if (rc > maxhwm) maxhwm = rc;
  return maxhwm;
}
