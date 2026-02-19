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

#include "ectrans/transi.h"
#include "transi_test.h"

// Following dummy functions are implementation details
// that don't contribute to this example. They could be
// replaced with grib_api functionality
void read_rspec( struct Trans_t* trans, double* rspec[], int nfld );

int main ( int argc, char **argv )
{
  trans_use_mpi( test_use_mpi() );

  double begin = transi_test_time();
  struct Trans_t trans;
  trans_new(&trans);
  double start;
  int N = 80;
  int nfld = 10;
  printf( "Grid resolution: N = %d\n",N);
  printf( "Number of fields to transform: %d\n",nfld);
  printf( "\n" );

  trans_init();

  // Read resolution information
  set_standard_rgg(&trans,N,-1);

  // Register resolution in trans library
  start = transi_test_time();
  trans_setup(&trans);


  if( trans.myproc == 1 ) printf( "trans_setup() ... ");
  if( trans.myproc == 1 ) print_time("",transi_test_time()-start);

  //for( nfld = 200; nfld<=250; ++nfld )
  {
    double* rspec;
    read_rspec(&trans,&rspec,nfld);
    double* rgp = malloc( sizeof(double) * nfld*trans.ngptot );

    // Inverse Transform
    struct InvTrans_t invtrans = new_invtrans(&trans);
      invtrans.nscalar   = nfld;
      invtrans.rspscalar = rspec;
      invtrans.rgp       = rgp;

    if( trans.myproc == 1 ) printf( "trans_invtrans(%3d) ...   ",nfld);
    start = transi_test_time();
    trans_invtrans(&invtrans);
    if( trans.myproc == 1 ) print_time("",transi_test_time()-start);

    free(rgp);
    free(rspec);
  }

  trans_delete(&trans);

  trans_finalize();

  if( trans.myproc == 1 ) print_time("END transi_timings   total: ",transi_test_time()-begin);


  return 0;
}

//---------------------------------------------------------------------------
// Dummy functions, used in this example



void read_rspec(struct Trans_t* trans, double* rspec[], int nfld )
{
  int i,j;
  *rspec = malloc( sizeof(double) * nfld*trans->nspec2 );
  for( i=0; i<trans->nspec2; ++i )
  {
    for( j=0; j<nfld; ++j )
    {
      (*rspec)[i*nfld+j] = 0;
    }
  }
}
