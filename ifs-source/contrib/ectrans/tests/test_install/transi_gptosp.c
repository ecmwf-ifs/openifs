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

#define TRANS_CHECK( CALL ) do {\
  int errcode = CALL;\
  if( errcode != TRANS_SUCCESS) {\
    printf("ERROR: %s failed @%s:%d:\n%s\n",#CALL,__FILE__,__LINE__,trans_error_msg(errcode));\
    abort();\
  }\
} while(0)

/*! @example transi_gptosp.c
 *
 * Transform gridpoint to spectral
 *
 * This is an example of how to setup and
 * transform global gridpoint data to
 * global spectral data
 */

// Following dummy functions are implementation details
// that don't contribute to this example. They could be
// replaced with grib_api functionality
void read_grid( struct Trans_t* trans );
void read_rgpg( struct Trans_t* trans, double* rgpg[], int* nfrom[], int* nfld );
void write_rspecg( struct Trans_t* trans, double* rspecg[], int nfld );


int main ( int arc, char **argv )
{
  int jfld;
  struct Trans_t trans;
  trans_new(&trans);

  // Read resolution information
  read_grid(&trans);

  // Register resolution in trans library
  trans_setup(&trans);

  // Declare global gridpoint data
  int nfld;
  double* rgpg = NULL;
  int* nfrom = NULL;

  // Read global gridpoint data (could be from grib file)
  read_rgpg(&trans,&rgpg,&nfrom,&nfld);

  // Distribute data to all procs
  double* rgp  = malloc( sizeof(double) * nfld *trans.ngptot  );
  struct DistGrid_t distgrid = new_distgrid(&trans);
    distgrid.nfrom = nfrom;
    distgrid.rgpg  = rgpg;
    distgrid.rgp   = rgp;
    distgrid.nfld  = nfld;
  trans_distgrid(&distgrid);


  // Transform gp to sp fields
  double* rspec = malloc( sizeof(double) * nfld*trans.nspec2 );

  struct DirTrans_t dirtrans = new_dirtrans(&trans);
    dirtrans.nscalar   = nfld;
    dirtrans.rgp       = rgp;
    dirtrans.rspscalar = rspec;
  trans_dirtrans(&dirtrans);


  // Gather all spectral fields
  double* rspecg = NULL;
  if( trans.myproc == 1 )
    rspecg = malloc( sizeof(double) * nfld*trans.nspec2g );

  int* nto =  malloc( sizeof(int) * nfld );
  for( jfld=0; jfld<nfld; ++jfld )
    nto[jfld] = 1;

  double* rnorm = malloc( sizeof(double) * nfld );
  struct SpecNorm_t specnorm = new_specnorm(&trans);
    specnorm.rspec = rspec;
    specnorm.nfld  = nfld;
    specnorm.rnorm = rnorm;
  trans_specnorm(&specnorm);

  for( jfld=0; jfld<nfld; ++jfld ) {
    if( trans.myproc == 1 ) {
      printf("specnorm[%d] = %f\n",jfld,rnorm[jfld]);
    }
  }

  struct GathSpec_t gathspec = new_gathspec(&trans);
    gathspec.rspec  = rspec;
    gathspec.rspecg = rspecg;
    gathspec.nfld   = nfld;
    gathspec.nto    = nto;
  trans_gathspec(&gathspec);

  // Write global spectral data (could be to grib file)
  if( trans.myproc == 1 )
    write_rspecg(&trans,&rspecg,nfld);

  // Deallocate and finalize
  free(rgp);
  free(rgpg);
  free(rspec);
  free(rspecg);
  free(nfrom);
  free(nto);
  free(rnorm);

  trans_delete(&trans);

  trans_finalize();

  return 0;
}

//---------------------------------------------------------------------------
// Dummy functions, used in this example

void read_grid(struct Trans_t* trans)
{
  int i;
  int T159[] = {
     18,  25,  36,  40,  45,  50,  60,  64,  72,  72,
     80,  90,  96, 100, 108, 120, 120, 125, 135, 144,
    144, 150, 160, 160, 180, 180, 180, 192, 192, 200,
    200, 216, 216, 216, 225, 225, 240, 240, 240, 243,
    250, 250, 256, 270, 270, 270, 288, 288, 288, 288,
    288, 288, 300, 300, 300, 300, 320, 320, 320, 320,
    320, 320, 320, 320, 320, 320, 320, 320, 320, 320,
    320, 320, 320, 320, 320, 320, 320, 320, 320, 320,
    320, 320, 320, 320, 320, 320, 320, 320, 320, 320,
    320, 320, 320, 320, 320, 320, 320, 320, 320, 320,
    320, 320, 320, 320, 300, 300, 300, 300, 288, 288,
    288, 288, 288, 288, 270, 270, 270, 256, 250, 250,
    243, 240, 240, 240, 225, 225, 216, 216, 216, 200,
    200, 192, 192, 180, 180, 180, 160, 160, 150, 144,
    144, 135, 125, 120, 120, 108, 100,  96,  90,  80,
     72,  72,  64,  60,  50,  45,  40,  36,  25,  18,
  };
  trans->ndgl  = sizeof(T159)/sizeof(int);
  trans->nloen = malloc( sizeof(T159) );
  for( i = 0; i<trans->ndgl; i++)  trans->nloen[i] = T159[i];

  // Assume Linear Grid
  trans->nsmax=(2*trans->ndgl-1)/2;
}


void read_rgpg(struct Trans_t* trans, double* rgpg[], int* nfrom[], int* nfld )
{
  int i;
  int jfld;
  if( trans->myproc == 1 ) printf("read_rpgp ...\n");
  *nfld = 2;
  if( trans->myproc == 1 )
  {
    *rgpg = malloc( sizeof(double) * (*nfld)*trans->ngptotg );
    for( i=0; i<trans->ngptotg; ++i )
    {
      (*rgpg)[0*trans->ngptotg + i] = 1.; // scalar field 1
      (*rgpg)[1*trans->ngptotg + i] = 2.; // scalar field 2
    }
  }
  *nfrom = malloc( sizeof(int) * (*nfld) );
  for (jfld=0; jfld<(*nfld); ++jfld)
    (*nfrom)[jfld] = 1;
  if( trans->myproc == 1 ) printf("read_rpgp ... done\n");
}

void write_rspecg(struct Trans_t* trans, double* rspecg[], int nfld )
{
  int i;
  if( trans->myproc == 1 ) printf("write_rspecg ...\n");
  for( i=0; i<trans->nspec2g; ++i )
  {
    // output global fields rspecg[i][0:nfld-1]
  }
  if( trans->myproc == 1 ) printf("write_rspecg ... done\n");
}

