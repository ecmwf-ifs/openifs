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
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "drhook.h"
#include "mpl.h"

extern void abor1_(const char msg[], int msglen);

extern int ec_sleep(const int nsec);

#pragma weak abor1_

// Forward declarations
extern void LinuxTraceBack(const char *prefix, const char *timestr, void *sigcontextptr);
extern void ec_microsleep(int usecs); // from ec_env.c
extern void fortran_mpi_abort(int rc);
extern int fortran_mpi_initialized();

static const char tabort_lockfile[] = "tabort_lock";

void tabort_delete_lockfile() {
  if (access(tabort_lockfile, F_OK) != -1) {
    // File is found
    remove(tabort_lockfile);
  }
}

void tabort_delete_lockfile_() {
  tabort_delete_lockfile();
}

void tabort_() {
  const int sig = SIGABRT;
  int rc = 128 + sig;
  static volatile sig_atomic_t irecur = 0;
  if (++irecur == 1) { // only one thread per task ever gets here
    // Only the fastest MPI task calls LinuxTraceBack -- avoids messy outputs
    int mpi_init = fortran_mpi_initialized();
    int nfirst = (mpi_init ? 0 : 1);
    if( nfirst == 0 ) {
      int fd = open(tabort_lockfile,O_CREAT|O_TRUNC|O_EXCL,S_IRUSR|S_IWUSR);
      if (fd >= 0) {
        close(fd);
        nfirst = 1;
      }
    }
    if (nfirst) {
      drhook_calltree();
      LinuxTraceBack(NULL,NULL,NULL);
    }
    else {
      const int nsecs = 100;
      ec_sleep(nsecs);
    }
    if (mpi_init) {
      fortran_mpi_abort(rc); // calls MPI_ABORT with MPI_COMM_WORLD
    }
    else {
      abort();
    }
  }
  // Still here ?? get the hell out of here ... now !!
  _exit(rc);
}

void abort_()
{
  if (abor1_) { // Call only if available
    static volatile sig_atomic_t irecur = 0;
    if (++irecur == 1) {
      const char msg[] = "Fortran ABORT()";
      abor1_(msg,strlen(msg));
    }
  }
  tabort_();
}

void _gfortran_abort()
{
  abort_();
}
