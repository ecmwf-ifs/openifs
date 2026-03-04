/*
 * (C) Copyright 2014- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "ectrans/transi.h"

#include "transi_test.h"

int main ( int arc, char **argv )
{
  int nloen[] = {18,25,36,40,45,54,60,64,72,72,80,90,96,100,108,120,120,128,135,144,144,150,160,160,180,180,180,192,192,200,200,216,216,216,225,225,240,240,240,256,256,256,256,288,288,288,288,288,288,288,288,288,300,300,300,300,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,320,300,300,300,300,288,288,288,288,288,288,288,288,288,256,256,256,256,240,240,240,225,225,216,216,216,200,200,192,192,180,180,180,160,160,150,144,144,135,128,120,120,108,100,96,90,80,72,72,64,60,54,45,40,36,25,18};

  int start = allocated();

  print_mem("Initially allocated: ",allocated()-start);

  if( strcmp(ectrans_version(),"transi/contrib") == 0 ) // equals
    trans_set_handles_limit(2);
  else /// Older versions cannot reuse existing handles, so allocate enough.
    trans_set_handles_limit(3);

  trans_use_mpi(false);

  trans_init();


  int iter=0;

  //int start_loop = allocated();
  for( iter=0; iter<3; ++iter )
  {
    printf("iteration %d\n",iter+1);
    int start_iter = allocated();
    struct Trans_t trans;

    trans_new(&trans);
    trans_set_resol(&trans,sizeof(nloen)/sizeof(int),nloen);
    trans_set_trunc(&trans,159);

    trans_setup(&trans);

    print_mem("Allocated in iteration:       ",allocated()-start_iter);

    trans_delete(&trans);

    print_mem("Possibly leaked in iteration: ",allocated()-start_iter);
    //print_mem("total leaked in loop:   ",allocated()-start_loop);

    // No memory leaks in subsequent iterations
    if( iter > 0 )
      ASSERT(allocated()-start_iter == 0);

  }

  printf( "trans_finalize()\n" );
  trans_finalize();

  print_mem("Total leaked: ",allocated()-start);


  display_mallinfo();

  return 0;
}

