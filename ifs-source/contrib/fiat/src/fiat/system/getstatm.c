/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "getstatm.h"

#if defined(LINUX)
#include <stdio.h>
#include <stdlib.h>

int getstatm(struct statm *sm)
{
  static int dont_bother = 0;
  if (!sm || dont_bother) {
    return -2;
  }
  else {
    FILE *statfile = fopen ("/proc/self/statm", "r");
    if (!statfile) {
      dont_bother = 1;
      return -1;
    }
    {
      int nelem = fscanf(statfile, "%d %d %d %d %d %d %d",
                         &(sm->size), &(sm->resident),
                         &(sm->shared), &(sm->trs), &(sm->drs),
                         &(sm->lrs), &(sm->dt));
    }
    fclose(statfile);
  }
  return 0;
}

#else

int getstatm(struct statm *sm)
{
  return -1; /* Not implemented */
}

#endif
