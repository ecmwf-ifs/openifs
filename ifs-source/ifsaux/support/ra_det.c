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

#if defined RS6K
#include <sys/rset.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/thread.h>

#define CHECK_PROC 1
#define CHECK_THREAD 2
#define ra_det ra_det_

/*------ Code supplied by Will Weir------------*/
int ra_det(int *det_type, int *prnt)
{
  int rs,preDetachrs,postDetachrs;
  int thread_num;
  rstype_t rstype;
  rsid_t rsid;
  rsethandle_t rset;
  char *preRSET;
  char *postRSET;
  char *tid_string;
/*  preRSET=malloc(25*sizeof(char));
  postRSET=malloc(25*sizeof(char)); */
  tid_string=malloc(4*sizeof(char));

  pid_t mypid;
  tid_t mytid;
  
  mypid=getpid();
  mytid=thread_self();

  rset=rs_alloc(RS_EMPTY);

  if ( *det_type==CHECK_PROC ) {
    rstype=R_PROCESS;
    rsid.at_pid=RS_MYSELF;
    sprintf(tid_string," %d",mypid);
  }
  else if ( *det_type==CHECK_THREAD ) {
    rstype=R_THREAD;
    rsid.at_tid=RS_MYSELF;
    thread_num=omp_get_thread_num();
    sprintf(tid_string," %d: %d",thread_num,mytid);
  }
  else {
    fprintf(stderr,"Error: rstype not R_PROCESS (1) or R_THREAD (2)\n");
/*    exit(-1); */
  }

/* Check current association */
  preDetachrs=ra_getrset(rstype,rsid,0,rset);

  switch (preDetachrs) {
    case RS_EFFECTIVE_RSET:
                   preRSET="RS_EFFECTIVE_RSET";
                                      break;
    case RS_PARTITION_RSET:
                   preRSET="RS_PARTITION_RSET";
                                      break;
    case RS_DEFAULT_RSET:
                     preRSET="RS_DEFAULT_RSET";
                                      break;
    case RS_THREAD_RSET:
                      preRSET="RS_THREAD_RSET";
                                      break;
    case RS_THREAD_PARTITION_RSET:
            preRSET="RS_THREAD_PARTITION_RSET";
                                      break;
/*----does not work on P6----------
    case RS_ADVISORY_RSET:
                    preRSET="RS_ADVISORY_RSET";
                                      break;
*/
  }

  rs=ra_detachrset(rstype,rsid,0);

  postDetachrs=ra_getrset(rstype,rsid,0,rset);

  switch (postDetachrs) {
    case RS_EFFECTIVE_RSET:
                   postRSET="RS_EFFECTIVE_RSET";
                                      break;
    case RS_PARTITION_RSET:
                   postRSET="RS_PARTITION_RSET";
                                      break;
    case RS_DEFAULT_RSET:
                     postRSET="RS_DEFAULT_RSET";
                                      break;
    case RS_THREAD_RSET:
                      postRSET="RS_THREAD_RSET";
                                      break;
    case RS_THREAD_PARTITION_RSET:
            postRSET="RS_THREAD_PARTITION_RSET";
                                      break;
/*----does not work on P6----------
    case RS_ADVISORY_RSET:
                    postRSET="RS_ADVISORY_RSET";
                                      break;
*/
  }


  if ( rs!=0 ) {
    perror("ra_detachrset error:");
/*    exit(-1); */
  } else {
    if(*prnt==1) fprintf(stderr,"ra_det:%s Successfully detached from rset: %s. Current rset: %s\n",tid_string,preRSET,postRSET);
  }
}
#endif
