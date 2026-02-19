/*
 * (C) Copyright 1989- ECMWF.
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * 
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction
 * 
 * (C) Copyright 1989- Meteo-France.
 * 
*/

#ifdef RS6K
#include <stdlib.h>
#include <stdio.h>
#include <sys/processor.h>
#include <unistd.h>
#include <sys/thread.h>
#include <sys/processor.h>
#include <string.h>

void aff_(int* icpu, int* ncpu)
{
   *icpu=mycpu();
   *ncpu=sysconf(_SC_NPROCESSORS_ONLN);
}

tid_t jbind_(proc)
  int *proc;
{
  tid_t t;
  t=thread_self();
  int ii;
/*
  fprintf(stderr,"proc=%d\n",*proc);
  fprintf(stderr,"thread=%d\n",t);
*/
  if( 0 != bindprocessor(BINDTHREAD,t,*proc) )
    fprintf(stderr,"bindprocessor failed for thread = %d\n",t);
}

int is_smt_on_()
{
  FILE *out;
  char v2[80];
  int i, jj;
  out=popen("mpstat -s | grep cpu | grep -v System | sed s/cpu//g","r");
  jj=-99;
  fscanf(out,"%d",&jj);
  return jj;
}

smtctl_(int *N, int *i0, int *i1)
{
  FILE *out0, *out1;
  int ii;
  char v1[80];
  out0=popen("bindprocessor -s 0","r");
  out1=popen("bindprocessor -s 1","r");
  fscanf(out0,"%s %s %s %s",v1,v1,v1,v1);
  fscanf(out1,"%s %s %s %s",v1,v1,v1,v1);
  for (ii=0; ii<*N; ii++) 
  {
    fscanf(out0,"%d",&i0[ii]);
    fscanf(out1,"%d",&i1[ii]);
  }
}
#endif
