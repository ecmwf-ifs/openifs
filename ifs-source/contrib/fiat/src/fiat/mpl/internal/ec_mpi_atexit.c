/*
 * (C) Copyright 2005- ECMWF.
 * (C) Copyright 2013- Meteo-France.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/* ec_mpi_atexit.c */

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "mpl.h"

/* A routine to be called at the very end in case MPI wasn't finalized */
/* Registered *only* by MPL_INIT */
/* Disable this feature via : export EC_MPI_ATEXIT=0 */

void ec_mpi_atexit_(void)
{
  char *env = getenv("EC_MPI_ATEXIT");
  int do_it = env ? atoi(env) : 1;
  static int callnum = 0;
  ++callnum;
  if (do_it) {
    if (callnum == 1) {
      /* register */
      atexit(ec_mpi_atexit_);
    }
    else if (callnum == 2) {
      /* action : finish MPI via F90 mpl_end (in mpl_bindc.F90) */
      mpl_end();
    }
  }
}

void ec_mpi_atexit(void)
{
  ec_mpi_atexit_();
}

