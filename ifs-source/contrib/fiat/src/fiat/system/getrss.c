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

#define getrss getrss_

#if defined(LINUX) || defined(__APPLE__)
static ll_t basesize = -1;
static size_t pagesize = 4096;
ll_t getrss()
{
  struct statm sm;
  ll_t rc = 0;
  if (getstatm(&sm) == 0) {
    if (basesize < 0) { /* the very first time */
      basesize = sm.resident;
      pagesize = getpagesize();
      if (pagesize <= 0) pagesize = 4096;
    }
    rc = (sm.resident - basesize) * pagesize;
  }
  return rc;
}
#else
ll_t getrss()
{
  ll_t rc = (ll_t)((char *)sbrk(0) - (char *)0);
  return rc;
}
#endif
