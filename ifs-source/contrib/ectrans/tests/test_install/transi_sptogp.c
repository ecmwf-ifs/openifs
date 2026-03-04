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

#include "ectrans/transi.h"

#define TRANS_CHECK( CALL ) do {\
  int errcode = CALL;\
  if( errcode != TRANS_SUCCESS) {\
    printf("ERROR: %s failed @%s:%d:\n%s\n",#CALL,__FILE__,__LINE__,trans_error_msg(errcode));\
    abort();\
  }\
} while(0)

/*! @example transi_sptogp.c
 *
 * Transform spectral to gridpoint
 *
 * This is an example of how to setup and
 * transform global spectral data to
 * global gridpoint data
 */

// Following dummy functions are implementation details
// that don't contribute to this example. They could be
// replaced with grib_api functionality
void read_grid( struct Trans_t* trans );
void read_rspecg( struct Trans_t* trans, double* rspecg[], int* nfrom[], int* nfld );
void write_rgpg( struct Trans_t* trans, double* rgpg[], int nfld );

int main ( int arc, char **argv )
{
  int jfld;
  struct Trans_t trans;
  TRANS_CHECK(trans_new(&trans));

  // Read resolution information
  read_grid(&trans);

  // Register resolution in trans library
  TRANS_CHECK( trans_setup(&trans) );

  // Declare global spectral data
  int nfld;
  double* rspecg = NULL;
  int* nfrom = NULL;

  // Read global spectral data (could be from grib file)
  read_rspecg(&trans,&rspecg,&nfrom,&nfld);

  // Distribute data to all procs
  double* rspec  = malloc( sizeof(double) * nfld *trans.nspec2  );
  struct DistSpec_t distspec = new_distspec(&trans);
    distspec.nfrom  = nfrom;
    distspec.rspecg = rspecg;
    distspec.rspec  = rspec;
    distspec.nfld   = nfld;
  TRANS_CHECK(trans_distspec(&distspec));


  // Transform sp to gp fields
  double* rgp = malloc( sizeof(double) * nfld*trans.ngptot );

  struct InvTrans_t invtrans = new_invtrans(&trans);
    invtrans.nscalar   = nfld;
    invtrans.rspscalar = rspec;
    invtrans.rgp       = rgp;
  TRANS_CHECK( trans_invtrans(&invtrans) );


  // Gather all gridpoint fields
  double* rgpg = NULL;
  if( trans.myproc == 1 )
    rgpg = malloc( sizeof(double) * nfld*trans.ngptotg );

  int* nto =  malloc( sizeof(int) * nfld );
  for( jfld=0; jfld<nfld; ++jfld )
    nto[jfld] = 1;

  struct GathGrid_t gathgrid = new_gathgrid(&trans);
    gathgrid.rgp  = rgp;
    gathgrid.rgpg = rgpg;
    gathgrid.nfld = nfld;
    gathgrid.nto  = nto;
  TRANS_CHECK( trans_gathgrid(&gathgrid) );

  // Write global spectral data (could be to grib file)
  if( trans.myproc == 1 )
    write_rgpg(&trans,&rgpg,nfld);

  // Deallocate and finalize
  free(rgp);
  free(rgpg);
  free(rspec);
  free(rspecg);
  free(nfrom);
  free(nto);

  TRANS_CHECK( trans_delete(&trans) );

  TRANS_CHECK( trans_finalize() );

  return 0;
}

//---------------------------------------------------------------------------
// Dummy functions, used in this example

void read_grid(struct Trans_t* trans)
{
  int i;

  trans->ndgl=160;
  trans->nloen=malloc(sizeof(int)*trans->ndgl);
  for( i=0; i<trans->ndgl; i++)  trans->nloen[i] = trans->ndgl*2;

  // Assume Linear Grid
  trans->nsmax=1279;
      //(2*trans->ndgl-1)/2;
}


void read_rspecg(struct Trans_t* trans, double* rspecg[], int* nfrom[], int* nfld )
{
  int i;
  int jfld;
  if( trans->myproc == 1 ) printf("read_rspecg ...\n");
  *nfld = 2;
  if( trans->myproc == 1 )
  {
    *rspecg = malloc( sizeof(double) * (*nfld)*trans->nspec2g );
    for( i=0; i<trans->nspec2g; ++i )
    {
      (*rspecg)[i*(*nfld) + 0] = (i==0 ? 1. : 0.); // scalar field 1
      (*rspecg)[i*(*nfld) + 1] = (i==0 ? 2. : 0.); // scalar field 2
    }
  }
  *nfrom = malloc( sizeof(int) * (*nfld) );
  for (jfld=0; jfld<(*nfld); ++jfld)
    (*nfrom)[jfld] = 1;
  if( trans->myproc == 1 ) printf("read_rspecg ... done\n");
}

void write_rgpg(struct Trans_t* trans, double* rgpg[], int nfld )
{
  int jfld;
  if( trans->myproc == 1 ) printf("write_rgpg ...\n");
  for( jfld=0; jfld<nfld; ++jfld )
  {
    // output global field rgpg[jfld]
  }
  if( trans->myproc == 1 ) printf("write_rgpg ... done\n");
}

