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
#include <sys/rset.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/thread.h>

#define CHECK_PROC 1
#define CHECK_THREAD 2
#define ra_check ra_check_

/*------ Code supplied by Will Weir------------*/
int ra_check(const int *check_type, int *prnt)
{
  int rs;
  int thread_num;
  rstype_t rstype;
  rsid_t rsid;
  unsigned int flags;
  rsethandle_t rset;

  pid_t mypid;
  tid_t mytid;
  
  mypid=getpid();
  mytid=thread_self();

  rset=rs_alloc(RS_EMPTY);

  if ( *check_type == CHECK_PROC )
  {
/*    fprintf(stderr,"ra_check: CHECK_PROC\n"); */
    rstype=R_PROCESS;
    rsid.at_pid=RS_MYSELF;
  }
  else if ( *check_type == CHECK_THREAD )
  {
/*    fprintf(stderr,"ra_check: CHECK_THREAD\n"); */
    rstype=R_THREAD;
    rsid.at_tid=RS_MYSELF;
    thread_num=omp_get_thread_num();
  }
  else
  {
    fprintf(stderr,"ra_check: check_type $d not CHECK_PROC or CHECK_THREAD\n",*check_type);
    return(1);
  }
  
/*  rsid.at_tid=RS_MYSELF; */
/*  rsid.at_tid=pid; */
/*  rsid.at_pid=pid; */

  rs=ra_getrset(rstype,rsid,0,rset);
  if ( rs < 0 ) {
    perror("ra_check: ra_getrset error:");
    exit(-1);
  } else if ( rs==RS_EFFECTIVE_RSET ) {
    if(*prnt==1) fprintf(stderr,"pid: %d  thread: %d tid: %d RSET type RS_EFFECTIVE_RSET\n",mypid,thread_num,mytid);
  }
  else if ( rs==RS_PARTITION_RSET ) {
    if(*prnt==1) fprintf(stderr,"pid: %d  thread: %d tid: %d RSET type RS_PARTITION_RSET\n",mypid,thread_num,mytid);
  }
  else if ( rs==RS_DEFAULT_RSET ) {
    if(*prnt==1) fprintf(stderr,"pid: %d  thread: %d tid: %d RSET type RS_DEFAULT_RSET\n",mypid,thread_num,mytid);
  }
  else if ( rs==RS_THREAD_RSET ) {
    if(*prnt==1) fprintf(stderr,"pid: %d  thread: %d tid: %d RSET type RS_THREAD_RSET\n",mypid,thread_num,mytid);
  }
  else if ( rs==RS_THREAD_PARTITION_RSET ) {
    if(*prnt==1) fprintf(stderr,"pid: %d  thread: %d tid: %d RSET type RS_THREAD_PARTITION_RSET\n",mypid,thread_num,mytid);
  }
/*--------does not work on P6----------
  else if ( rs==RS_ADVISORY_RSET ) {
    if(*prnt==1) fprintf(stderr,"pid: %d  thread: %d tid: %d RSET type RS_ADVISORY_RSET\n",mypid,thread_num,mytid);
  }
*/
  if ( rs == RS_DEFAULT_RSET ) {
    if(*prnt==1) fprintf(stderr,"ra_check: RS_DEFAULT_RSET returning 1\n");
    return(1);
  } else {
    if(*prnt==1) fprintf(stderr,"ra_check: NOT RS_DEFAULT_RSET returning 0\n");
    return(0);
  }

}
#endif
