/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ectrans/transi.h"
#include "transi_test.h"

// ----------------------------------------------------------------------------

void test_invtrans_adjoint(int nlon, int nlat, int nsmax)
{
double adjoint_tol = 1.e-6;
printf("test_invtrans_adjoint( nlon=%d, nlat=%d, nsmax=%d )\n",nlon,nlat,nsmax);

int nout = 2;

struct Trans_t trans;
TRANS_CHECK( trans_new(&trans) );

int* nloen  = malloc( sizeof(int) * nlat);
{
    int i;
    for( i=0; i<nlat; ++i )
      nloen[i] = nlon;
}
TRANS_CHECK( trans_set_resol(&trans,nlat, nloen) );
TRANS_CHECK( trans_set_trunc(&trans,nsmax) );

trans.fft = TRANS_FFTW;
trans.flt = 0;

TRANS_CHECK( trans_setup(&trans) );
printf("ndgl  = %d\n",trans.ndgl);
printf("nsmax = %d\n",trans.nsmax);
printf("ngptotg  = %d\n",trans.ngptotg);
printf("llatlon  = %d\n",trans.llatlon);
ASSERT(trans.ngptotg == nlon*nlat);

// Allocate gridpoint data
int nscalar = 2;
int nvordiv = 1;
int nfld  = 2*nvordiv+nscalar;
double* rgp  = malloc( sizeof(double) * nfld * trans.ngptot );

// Load data on proc 1
double* rgpg = NULL;
if( trans.myproc == 1 )
{
  rgpg = malloc( sizeof(double) * nfld*trans.ngptotg );
  int i;
  for( i=0; i<trans.ngptotg; ++i )
  {
    rgpg[0*trans.ngptotg+i] = 1.; // U
    rgpg[1*trans.ngptotg+i] = 2.; // V
    rgpg[2*trans.ngptotg+i] = 3.; // scalar 1
    rgpg[3*trans.ngptotg+i] = 4.; // scalar 2
  }
}

// Distribute gridpoint data from proc 1
int* nfrom = malloc( sizeof(int) * nfld );
nfrom[0] = 1; // U
nfrom[1] = 1; // V
nfrom[2] = 1; // scalar 1
nfrom[3] = 1; // scalar 2

printf("trans_distgrid()\n");
struct DistGrid_t distgrid = new_distgrid(&trans);
  distgrid.nfrom = nfrom;
  distgrid.rgpg  = rgpg;
  distgrid.rgp   = rgp;
  distgrid.nfld  = nfld;
TRANS_CHECK( trans_distgrid(&distgrid) );

// Allocate spectral data

double* rspscalar = calloc( nscalar*trans.nspec2, sizeof(double));
double* rspvor    = calloc( nvordiv*trans.nspec2, sizeof(double));
double* rspdiv    = calloc( nvordiv*trans.nspec2, sizeof(double));

printf("trans_invtrans_adj()\n");
struct InvTransAdj_t invtrans_adj = new_invtrans_adj(&trans);
  invtrans_adj.nscalar   = nscalar;
  invtrans_adj.nvordiv   = nvordiv;
  invtrans_adj.rgp       = rgp;
  invtrans_adj.rspscalar = rspscalar;
  invtrans_adj.rspvor    = rspvor;
  invtrans_adj.rspdiv    = rspdiv;
TRANS_CHECK( trans_invtrans_adj(&invtrans_adj) );

// Gather spectral field (for fun)
int* nto = malloc( sizeof(int) * nscalar );
nto[0] = 1;
nto[1] = 1;

double* rspscalarg = NULL;
if( trans.myproc == 1 )
  rspscalarg = malloc( sizeof(double) * nscalar*trans.nspec2g );

printf("trans_gathspec()\n");
struct GathSpec_t gathspec = new_gathspec(&trans);
  gathspec.rspec  = rspscalar;
  gathspec.rspecg = rspscalarg;
  gathspec.nfld   = nscalar;
  gathspec.nto    = nto;
TRANS_CHECK( trans_gathspec(&gathspec) );

double adj_value1 = 0.0;
if( trans.myproc == 1 )
{
  int i,j;

  for( j=0; j<nscalar; ++j)
  {
    for( i=0; i<nout; ++i )
      printf("rspscalarg[%d][%d] : %f\n",j,i,rspscalarg[i*nscalar+j]);

    // The first 2*nlat terms are on m=0
    for( i=0; i<2*nlat; ++i )
      adj_value1 += rspscalarg[i*nscalar+j] * rspscalarg[i*nscalar+j];

    for( i=2*nlat; i<trans.nspec2g; ++i )
      adj_value1 += 2.0 * rspscalarg[i*nscalar+j] * rspscalarg[i*nscalar+j];

    for( i=0; i<trans.nspec2g; ++i )
      printf("error -> rspscalarg[fld=%d][wave=%d] : %f\n",j,i,rspscalarg[i*nscalar+j]);

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

printf("trans_gathgrid()\n");
// Gather gridpoint fields
struct GathGrid_t gathgrid = new_gathgrid(&trans);
  gathgrid.rgp  = rgp;
  gathgrid.rgpg = rgpg;
  gathgrid.nto  = nfrom;
  gathgrid.nfld = nfld;
TRANS_CHECK( trans_gathgrid(&gathgrid) );


double adj_value2 = 0.0;
if( trans.myproc == 1 )
{
  int i,j;
  for( j=0; j<nfld; ++j)
  {
    for( i=0; i<nout; ++i )
    {
      printf("rgpg[fld=%d][pt=%d] : %f\n",j,i,rgpg[j*trans.ngptotg+i]);
    }
    if( j>=2 )
    {
      for( i=0; i<trans.ngptotg; ++i )
        adj_value2 += rgpg[j*trans.ngptotg+i] * (double)(j+1);
    }
  }
}
printf("rgpg[adjval1=%f][adjval2=%f] :\n", adj_value1, adj_value2);
ASSERT( fabs(adj_value1 - adj_value2 )/adj_value1 < adjoint_tol );

// Deallocate arrays
free(nloen);
free(rgp);
free(rgpg);
free(rspscalar);
free(rspscalarg);
free(rspvor);
free(rspdiv);
free(nfrom);
free(nto);

TRANS_CHECK( trans_delete(&trans) );

}
//
// ----------------------------------------------------------------------------

int main ( int arc, char **argv )
{
  trans_use_mpi( test_use_mpi() );

  setbuf(stdout,NULL); // unbuffered stdout

// The adjoint test works for standard gaussian latitude grid
//   with no points on the equator or poles.
// nsmax = nlat - 1

  printf("-----------------------------\n");
  test_invtrans_adjoint(8,4,3);

  TRANS_CHECK( trans_finalize() );

  return 0;
}

