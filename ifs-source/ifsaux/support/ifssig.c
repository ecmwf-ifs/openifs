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

/***********************************************************************
.TITLE    ECMWF Utility
.NAME     IFSSIG
.SECTION  METAPPS
.AUTHOR   Otto Pesonen
.DATE     29-JUL-1996 / 22-JUL-1996 / OP
.VERSION  1.0
.LANGUAGE ANSI-C
.FILE     ifssig.c
.DATE     22-AUG-1996 / 22-AUG-1996 / OP
.VERSION  1.1
*
*  IFS signal handler to produce checkpoint files.
*
*  Compile:
*
*  cc -c ifssig.c
*
***********************************************************************/

#include <stdio.h>
#include <signal.h>
#include <unistd.h>

#include "mpl.h"

#undef SIGNAL_STACK                /* Signals do not stack */

#define CHECKPOINT SIGHUP
#define RESTART    SIGINT

#define IFSSIGMASK (1<<(CHECKPOINT-1)) | (1<<(RESTART-1))

static int *re_start;
static int *go_on;

static void (*sigcheck)();
static void (*sigrestart)();

static int    mask;

static void catch(int sig)
{
  switch(sig)
  {
    case CHECKPOINT:
      *go_on = 0;
      *re_start = 1;

#ifdef SIGNAL_STACK
      if( sigcheck != SIG_DFL && sigcheck != SIG_IGN )
        sigcheck(sig);
#endif
      signal(sig,catch);
      printf("ifssig:got CHECKPOINT signal\n");
      break;

    case RESTART:
      *go_on = 1;
      *re_start = 1;

#ifdef SIGNAL_STACK
      if( sigrestart != SIG_DFL && sigrestart != SIG_IGN )
        sigrestart(sig);
#endif

      signal(sig,catch);
      printf("ifssig:got RESTART signal\n");
      break;
  }
  fflush(stdout);
}

void ifssigb()
/**************************************************************************
?  Unblock IFS signals. This will call the signal-handler if there is any
|  such signals pending.
|  Once the signals are processed, block them until we come in here again
***********************************************************************/
{
  sigset_t set;

  sigemptyset(&set);
  sigaddset(&set, CHECKPOINT);
  sigaddset(&set, RESTART);
  sigprocmask(SIG_UNBLOCK, &set, 0);

  sigemptyset(&set);
  sigaddset(&set, CHECKPOINT);
  sigaddset(&set, RESTART);
  sigprocmask(SIG_BLOCK, &set, 0);
}

void ifssigb_() { ifssigb(); }
void ifssig(int *goon, int *restart)
/***********************************************************************
?  Register for IFS signals.
|  And block them untill we want to process them.
***********************************************************************/
{
  go_on = goon;
  re_start = restart;

  *go_on = 0;
  *re_start = 0;

  if( mpl_myrank() == 1 )
  {
    sigcheck = signal(CHECKPOINT, catch);
    sigrestart = signal(RESTART, catch);
  }
  else
  {
    sigcheck = signal(CHECKPOINT, SIG_IGN);
    sigrestart = signal(RESTART, SIG_IGN);
  }

  ifssigb();
  fflush(stdout);
}

void ifssig_(int *goon, int *restart) { ifssig(goon, restart); }

#ifdef DEBUG

main()
{
  int i;
  int goon    = 0;
  int restart = 0;

  printf("PID is %d\n",getpid());

  ifssig(&goon, &restart);

  for( i=0 ; i<33 && !goon ; i++ )
  {
    sleep(3);

    ifssigb();

    if(goon)
    { 
      printf("i=%d goon %d\n",i,goon);
    }

    if(restart)
    {
      printf("i=%d restart %d\n",i,restart);
      restart=0;
    }
  }
}

#endif

void sigmaster_()
{
  FILE *out;
  char v1[132];
  int i;
  size_t len;
  len=132;
/*   out=popen("hostname","r"); */
  i=gethostname(v1,len);
  i=getpid();
  fprintf(stderr,"MASTER-HOSTNAME-PID %s %d\n",v1,i);
}


