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
#include <math.h>
#include <ctype.h>

#include "drhook.h"

#ifdef __INTEL_COMPILER
// Otherwise SIGFPE is not triggered
#pragma STDC FENV_ACCESS ON
#endif

static int trig_sigfpe = 0;

void set_trig_sigfpe() {
    char* env = getenv("SIGFPE");
    trig_sigfpe = env ? atoi(env) : 1;
}

double mysin(double x)
{
  DRHOOK_START(mysin);
  {
    double div = 180;
    static int jcnt = 0;
    if (trig_sigfpe && (++jcnt%15000 == 0)) {
      fprintf(stderr, "Trigger div-by-zero\n" );
      div = 0; /* Trigger divide-by-zero  i.e. "My sin" ;-) */
    }
    x = sin(x/div);
  }
  DRHOOK_END(); /* mysin */
  return x;
}

double mycos(double x)
{
  DRHOOK_START(mycos);
  x = 1 - mysin(x);
  DRHOOK_END(); /* mycos */
  return x;
}

double mysqrt(double x)
{
  DRHOOK_START(mysqrt);
  x = sqrt(x);
  DRHOOK_END(); /* mysqrt */
  return x;
}

double sub(int j)
{
  double ans = 0;
  DRHOOK_START(sub);
  {
    int i, i1 = (j-1)*100, i2 = j*100;
    for (i=i1; i<i2; i++) {
      ans += mysin(i) + mysqrt(fabs(mycos(i)));
    }
  }
  DRHOOK_END(); /* sub */
  return ans;
}

int main(int argc, char *argv[])
{
  int rc = 0;
  DRHOOK_START(real_main);

  set_trig_sigfpe();
  fprintf(stderr,"trig_sigfpe = %d\n",trig_sigfpe);
  
  {
    int j, n = 1000;
    DRHOOK_START(LOOP_BLOCK);
    for (j=0; j<n; j++) {
      (void) sub(j);
    }
    DRHOOK_END(); /* LOOP_BLOCK */
  }

  DRHOOK_END(); /* real_main */
  return rc;
}
