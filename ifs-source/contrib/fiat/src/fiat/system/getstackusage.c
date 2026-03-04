/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

typedef long long int ll_t;
typedef unsigned long long int ull_t;

#if defined(LINUX) && defined(USE_MEMORY_MONITOR)
#include <stdio.h>
#include <stdlib.h>

ll_t getstackusage_()
{
  ll_t rc = 0;
  static int dont_bother = 0;
  if (dont_bother) {
    rc = -2;
  }
  else {
    FILE *statfile = fopen ("/proc/self/stat", "r");
    if (!statfile) {
      dont_bother = 1;
      rc = -1;
    }
    else {
      char dm[80];
      ull_t startstack, kstkesp; /* stack start & ESP, the 28th and 29th columns, respectively */
      /* Maybe not the brightest coding, but has to suffice for now (SS) */
      int nelem = fscanf(statfile, 
		     "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %llu %llu ",
		     dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,&startstack,&kstkesp);
      if (nelem != 29) {
	dont_bother = 1;
	rc = -3;
      }
      else {
	rc = (ll_t)(startstack - kstkesp);
      }
      fclose(statfile);
    }
  }
  return rc;
}

#else

ll_t getstackusage_() { return 0L; }

#endif
