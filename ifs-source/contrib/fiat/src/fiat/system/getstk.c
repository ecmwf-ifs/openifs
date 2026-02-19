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
#include <string.h>
#include <unistd.h>

typedef  long long int  ll_t;

static ll_t maxstack = 0;

#if defined(LINUX)

ll_t getvmpeak_()
{
  ll_t virtmempeak = 0;
  FILE *fp = fopen("/proc/self/status","r");
  if (fp) {
    char in[4096];
    while (fgets(in,sizeof(in),fp) == in) {
      if (strncmp(in,"VmPeak:",7) == 0) {
	ll_t value;
	int nf = sscanf(in,"%*s %lld kB",&value);
	if (nf == 1) virtmempeak = value * (ll_t) 1024;
	break;
      }
    }
    fclose(fp);
  }
  return virtmempeak;
}

ll_t linux_getstackusage_()
{
  ll_t stackused = 0;
  FILE *fp = fopen("/proc/self/status","r");
  if (fp) {
    char in[4096];
    while (fgets(in,sizeof(in),fp) == in) {
      if (strncmp(in,"VmStk:",6) == 0) {
	ll_t value;
	int nf = sscanf(in,"%*s %lld kB",&value);
	if (nf == 1) stackused = value * (ll_t) 1024;
	break;
      }
    }
    fclose(fp);
  }
  return stackused;
}

ll_t getstk_()
{
  /*
  extern ll_t getstackusage_();
  ll_t stackused = getstackusage_();
  */
  extern ll_t linux_getstackusage_();
  ll_t stackused = linux_getstackusage_();
  if (stackused > maxstack) maxstack = stackused;
  return stackused;
}

#else
ll_t
getstk_() 
{ 
  extern ll_t getstackusage_();
  static ll_t init_stack = -1;
  ll_t stackused = 0;
  if (init_stack == -1) init_stack = getstackusage_();
  stackused = getstackusage_() - init_stack;
  if (stackused > maxstack) maxstack = stackused;
  return stackused;
}

#endif

/* Maximum stacksize encountered */

ll_t
getmaxstk_()
{
  ll_t stackused = getstk_();
  if (stackused > maxstack) maxstack = stackused;
  return maxstack;
}

#if !defined(LINUX)
ll_t getvmpeak_() { return 0L; }
#endif

