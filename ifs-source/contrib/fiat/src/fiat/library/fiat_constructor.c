/*
 * (C) Copyright 2022- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdio.h>
#include <stdlib.h>

#ifdef FIAT_ATTRIBUTE_CONSTRUCTOR_SUPPORTED
#define ATTRIBUTE_CONSTRUCTOR __attribute__((constructor))
#else
#define ATTRIBUTE_CONSTRUCTOR
#endif

// Forward declarations
extern int ec_mpirank();
extern void tabort_delete_lockfile();

/*
 * fiat_debug()
 *
 * Access FIAT_DEBUG environment variable. Default treated as if FIAT_DEBUG=0
 */
static int fiat_debug() {
  char* env = getenv("FIAT_DEBUG");
  return env ? atoi(env) : 0;
}

/*
 * fiat_constructor() is called before main() upon loading of this library
 */
void ATTRIBUTE_CONSTRUCTOR fiat_constructor() {
  int mpi_rank = ec_mpirank();
  if (mpi_rank == 0 && fiat_debug()) {
    fprintf(stderr,"FIAT_DEBUG fiat_constructor()\n");
  }
  if( mpi_rank == 0 ) {
    tabort_delete_lockfile();
  }
}
