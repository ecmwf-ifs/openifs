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

#ifdef VPP
#include <stdio.h>
#include <sys/param.h>
#include <sys/types.h>
#include <stropts.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/procfs.h>
#include <sys/proc.h>
                                  /*********************************/
void  getpidinfo_ (               /**  To get process pid info    **/
    int   *pid                    /*********************************/
      )
 
  {
    int  fid;
    int  ipid;
    int  ierr;
    int  namelen;
    char *pname = "/proc/xxxxx";
    static prpsinfo_t info;
    static proc_t     proc;
    long  vrtsize;
    int   stksize;
    u_int brksize;
    long  cpu_time;

     if (*pid != 0)
      { ipid = *pid;
        sprintf (pname,"/proc/%5.5d",ipid); }
    else
      { ipid = getpid();
        sprintf (pname,"/proc/%5.5d",ipid); }
 
    fid = open (pname, O_RDONLY, 0);
    ioctl (fid, PIOCPSINFO, (char *)&info);
    ioctl (fid, PIOCGETPR,  (char *)&proc);
    close (fid);
 
    vrtsize = info.pr_size;
    vrtsize = vrtsize*32*1024;
    stksize = proc.p_stksize;
    brksize = proc.p_brksize;
    cpu_time = info.pr_time.tv_sec + info.pr_time.tv_nsec*0.0000000001;
 
    fprintf(stderr,"%-5s  %-12s  ","Pid","StackSize");
    fprintf(stderr,"%-12s  %-12s  ","HeapSize","VirtualSize");
    fprintf(stderr,"%-5s  \n","Time");
    fprintf(stderr,"%-5d  %-12d  ",ipid,stksize);
    fprintf(stderr,"%-12d  %-12d  ",brksize,vrtsize);
    fprintf(stderr,"%-5d  \n",cpu_time);
 
  } /* getpidinfo */
#endif
