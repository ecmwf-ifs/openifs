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

#include "transi_test.h"

static bool check_values = true;

// ----------------------------------------------------------------------------

void test_sptogp(int nlon, int nlat, int nsmax)
{
  double sptogp_tol = 1.e-6;
  printf("test_sptogp( nlon=%d, nlat=%d, nsmax=%d )\n",nlon,nlat,nsmax);

  int nout = 2;

  struct Trans_t trans;
  TRANS_CHECK( trans_new(&trans) );
  TRANS_CHECK( trans_set_resol_lonlat(&trans,nlon,nlat) );
  TRANS_CHECK( trans_set_trunc(&trans,nsmax) );

  trans.fft = TRANS_FFTW;
  trans.flt = 0;

  TRANS_CHECK( trans_setup(&trans) );
  printf("ndgl  = %d\n",trans.ndgl);
  printf("nsmax = %d\n",trans.nsmax);
  printf("ngptotg  = %d\n",trans.ngptotg);

  ASSERT(trans.ngptotg == nlon*nlat);
  if( trans.nproc == 1 ) {
    ASSERT(trans.ngptot  == nlon*nlat + nlon*(nlat%2));
  }
  ASSERT(trans.nspec2 == trans.nspec2g);
  ASSERT( trans.nproc == 1 );

  // In case of odd number of latitudes, there is one latitude extra in distributed gridpoints
  // This can only be checked with nproc==1

  // Allocate gridpoint data
  int nscalar = 1;
  int nvordiv = 1;
  int nfld  = 2*nvordiv+nscalar;
  double* rgp  = malloc( sizeof(double) * nfld * trans.ngptot );

  // Load data on proc 1
  double* rgpg = NULL;
  if( trans.myproc == 1 )
  {
    rgpg = malloc( sizeof(double) * nfld*trans.ngptotg );
  }


  // // Allocate spectral data
  double* rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
  double* rspvor    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspdiv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );

  // Gather spectral field (for fun)
  int* nto = malloc( sizeof(int) * nfld );
  int k;
  for( k=0; k<nfld; ++k )
    nto[k] = 1;

  int* nfrom = nto;

  double* rspscalarg = NULL;
  if( trans.myproc == 1 ) {
    rspscalarg = malloc( sizeof(double) * nscalar*trans.nspec2g );
    int i,j;
    for( j=0; j<nscalar; ++j ) {
      for( i=0; i<trans.nspec2g; ++i ) {
        rspscalarg[i*nscalar+j] = 0.;
      }
    }
    for( j=0; j<nvordiv; ++j ) {
      for( i=0; i<trans.nspec2g; ++i ) {
        rspvor[i*nvordiv+j] = 0.;
        rspdiv[i*nvordiv+j] = 0.;
    }
  }
}


  printf("trans_distspec()\n");
  // Distribute spectral field (for fun)
  struct DistSpec_t distspec = new_distspec(&trans);
    distspec.rspec  = rspscalar;
    distspec.rspecg = rspscalarg;
    distspec.nfld   = nscalar;
    distspec.nfrom  = nto;
  TRANS_CHECK( trans_distspec(&distspec) );

  printf("trans_invtrans()\n");
  // Inverse Transform
  struct InvTrans_t invtrans = new_invtrans(&trans);
    invtrans.nscalar   = nscalar;
    invtrans.nvordiv   = nvordiv;
    invtrans.rspscalar = rspscalar;
    invtrans.rspvor    = rspvor;
    invtrans.rspdiv    = rspdiv;
    invtrans.rgp       = rgp;
  TRANS_CHECK( trans_invtrans(&invtrans) );

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nfld; ++j)
    {
      for( i=0; i<nout; ++i )
      {
        printf("rgp[%d][%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
      }
      if( check_values )
      {
        for( i=0; i<trans.ngptot; ++i )
        {
          //if ( fabs(rgp[j*trans.ngptot+i] - (j+1) )/(double)(j+1) > sptogp_tol )
          //  printf("error --> rgp[fld=%d][pt=%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
          ASSERT( fabs(rgp[j*trans.ngptot+i]) < sptogp_tol );
        }
      }
    }
  }



  printf("trans_gathgrid()\n");
  // Gather gridpoint fields
  struct GathGrid_t gathgrid = new_gathgrid(&trans);
    gathgrid.rgp  = rgp;
    gathgrid.rgpg = rgpg;
    gathgrid.nto  = nfrom;
    gathgrid.nfld = nfld;
  TRANS_CHECK( trans_gathgrid(&gathgrid) );


  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nfld; ++j)
    {
      for( i=0; i<nout; ++i )
      {
        printf("rgpg[fld=%d][pt=%d] : %f\n",j,i,rgpg[j*trans.ngptotg+i]);
      }
      if( check_values )
      {
        for( i=0; i<trans.ngptotg; ++i )
        {
          ASSERT( fabs(rgp[j*trans.ngptot+i]) < sptogp_tol );
        }
      }
    }
  }

  // Deallocate arrays
  free(rgp);
  free(rgpg);
  free(rspscalar);
  free(rspscalarg);
  free(rspvor);
  free(rspdiv);
  free(nto);

  TRANS_CHECK( trans_delete(&trans) );
}

// ----------------------------------------------------------------------------

int main ( int arc, char **argv )
{
  trans_use_mpi( test_use_mpi() );

  setbuf(stdout,NULL); // unbuffered stdout


  // ATLAS-149: not all values are assigned, unless FOUBUF_IN(:) = 0 in ltinv_ctl_mod.F90

  //test_sptogp(360,37,1279); // As in MIR-183 ( 1/5 lonlat grid )
  test_sptogp(6,3,63); //  60/90 lonlat grid


  TRANS_CHECK( trans_finalize() );

  return 0;
}

