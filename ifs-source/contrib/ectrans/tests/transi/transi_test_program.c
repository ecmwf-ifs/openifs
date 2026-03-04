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

void read_grid(struct Trans_t*);

int main ( int arc, char **argv )
{
  trans_use_mpi( test_use_mpi() );

  printf("ectrans version int = %d\n",ectrans_version_int());
  printf("ectrans version = %s\n",ectrans_version());
  printf("ectrans version str = %s\n",ectrans_version_str());
  printf("ectrans git sha1 [0:7] = %s\n",ectrans_git_sha1_abbrev(7));
  printf("ectrans git sha1 [0:12] = %s\n",ectrans_git_sha1_abbrev(12));
  printf("ectrans git sha1 = %s\n",ectrans_git_sha1());

  //printf("transi started\n");
  int nout = 3;
  struct Trans_t trans;
  trans_new(&trans);

  read_grid(&trans);
  trans_setup(&trans);
  trans_inquire(&trans,"numpp,ngptotl,nmyms,nasm0,npossp,nptrms,nallms,ndim0g,nvalue");
  trans_inquire(&trans,"nfrstlat,nlstlat,nptrlat,nptrfrstlat,nptrlstlat,nsta,nonl,ldsplitlat");
  trans_inquire(&trans,"nultpp,nptrls,nnmeng");
  trans_inquire(&trans,"rmu,rgw,npms,rlapin,ndglu");

  //Check values of numpp
  if( trans.myproc == 1 )
  {
    printf("nprtrw = %d\n",trans.nprtrw);
    int i;
    for( i=0; i<trans.nprtrw; ++i)
      printf("%d : %d\n",i,trans.numpp[i]);
  }

  // Allocate gridpoint data
  int nscalar = 2;
  int nvordiv = 1;
  int nfld  = 2*nvordiv+nscalar;
  double* rgp  = malloc( sizeof(double) * nfld *trans.ngptot  );

  // Load data on proc 1
  double* rgpg = NULL;
  if( trans.myproc == 1 )
  {
    rgpg = malloc( sizeof(double) * 4*trans.ngptotg );
    int i;
    for( i=0; i<trans.ngptotg; ++i )
    {
      rgpg[0*trans.ngptotg+i] = 1.; // U
      rgpg[1*trans.ngptotg+i] = 2.; // V
      rgpg[2*trans.ngptotg+i] = 3.; // scalar 1
      rgpg[3*trans.ngptotg+i] = 4.; // scalar 2
    }
  }
  // if( trans.myproc == 2 )
  // {
  //   rgpg = malloc( sizeof(double) * 1*trans.ngptotg );
  //   int i;
  //   for( i=0; i<trans.ngptotg; ++i )
  //   {
  //     rgpg[i+0*trans.ngptotg] = 4.; // scalar 2
  //   }
  // }

  // Distribute gridpoint data from proc 1
  int* nfrom = malloc( sizeof(int) * nfld );
  nfrom[0] = 1; // U
  nfrom[1] = 1; // V
  nfrom[2] = 1; // scalar 1
  nfrom[3] = 1; // scalar 2

  struct DistGrid_t distgrid = new_distgrid(&trans);
    distgrid.nfrom = nfrom;
    distgrid.rgpg  = rgpg;
    distgrid.rgp   = rgp;
    distgrid.nfld  = nfld;
  trans_distgrid(&distgrid);

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nfld; ++j)
    {
      for( i=0; i<nout; ++i )
      {
        printf("rgp[%d][%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
      }
      for( i=0; i<trans.ngptot; ++i )
      {
        if( fabs(rgp[j*trans.ngptot+i]-(j+1)) > 1.e-5)
          printf("rgp[%d][%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
      }

    }
  }

  // Allocate spectral data

  double* rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
  double* rspvor    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspdiv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );

  // Direct Transform
  struct DirTrans_t dirtrans = new_dirtrans(&trans);
    dirtrans.nscalar   = nscalar;
    dirtrans.nvordiv   = nvordiv;
    dirtrans.rgp       = rgp;
    dirtrans.rspscalar = rspscalar;
    dirtrans.rspvor    = rspvor;
    dirtrans.rspdiv    = rspdiv;
  trans_dirtrans(&dirtrans);

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nscalar; ++j)
    {
      for( i=0; i<nout; ++i )
        printf("rspscalar[%d][%d] : %f\n",j,i,rspscalar[i*nscalar+j]);
      for( i=0; i<trans.nspec2; ++i )
      {
        if( fabs(rspscalar[i*nscalar+j]) > 1.e-5)
          printf("rspscalar[%d][%d] : %f\n",j,i,rspscalar[i*nscalar+j]);
      }
    }
  }

  // Allocate fields for u*cos(theta) and v*cos(theta)
  double* rspu    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );

  // Convert vorticity & divergence to u*cos(theta) & v*cos(theta)
  printf("Converting spectral vorticity-divergence to u*cos(lat)-v*cos(lat)...\n");
  struct VorDivToUV_t vordiv_to_UV = new_vordiv_to_UV();
    vordiv_to_UV.rspvor = rspvor;
    vordiv_to_UV.rspdiv = rspdiv;
    vordiv_to_UV.rspu   = rspu;
    vordiv_to_UV.rspv   = rspv;
    vordiv_to_UV.nfld   = nvordiv;
    vordiv_to_UV.ncoeff = trans.nspec2;
    vordiv_to_UV.nsmax  = trans.nsmax;
  trans_vordiv_to_UV(&vordiv_to_UV);
  printf("Converting spectral vorticity-divergence to u*cos(lat)-v*cos(lat)...done\n");


  // Gather spectral field (for fun)
  int* nto = malloc( sizeof(int) * nscalar );
  nto[0] = 1;
  nto[1] = 1;

  double* rspscalarg = NULL;
  if( trans.myproc == 1 )
    rspscalarg = malloc( sizeof(double) * nscalar*trans.nspec2g );

  struct GathSpec_t gathspec = new_gathspec(&trans);
    gathspec.rspec  = rspscalar;
    gathspec.rspecg = rspscalarg;
    gathspec.nfld   = nscalar;
    gathspec.nto    = nto;
  trans_gathspec(&gathspec);

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nscalar; ++j)
    {
      for( i=0; i<nout; ++i )
        printf("rspscalarg[%d][%d] : %f\n",j,i,rspscalarg[i*nscalar+j]);
      for( i=0; i<trans.nspec2g; ++i )
      {
        if( fabs(rspscalarg[i*nscalar+j]) > 1.e-5 && i > 0)
          printf("rspscalarg[%d][%d] : %f\n",j,i,rspscalarg[i*nscalar+j]);
      }
    }
  }

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nscalar; ++j)
    {
      for( i=0; i<nout; ++i )
        printf("rspscalarg[%d][%d] : %f\n",j,i,rspscalarg[i*nscalar+j]);
      for( i=0; i<trans.nspec2g; ++i )
      {
        if( fabs(rspscalarg[i*nscalar+j]) > 1.e-5 && i > 0)
          printf("rspscalarg[%d][%d] : %f\n",j,i,rspscalarg[i*nscalar+j]);
      }
    }
  }

  // Allocate fields for u*cos(theta) and v*cos(theta)
  double* rspvorg    = malloc( sizeof(double) * nvordiv*trans.nspec2g );
  double* rspdivg    = malloc( sizeof(double) * nvordiv*trans.nspec2g );

  gathspec = new_gathspec(&trans);
    gathspec.rspec  = rspvor;
    gathspec.rspecg = rspvorg;
    gathspec.nfld   = nvordiv;
    gathspec.nto    = nto;
  trans_gathspec(&gathspec);

  gathspec = new_gathspec(&trans);
    gathspec.rspec  = rspdiv;
    gathspec.rspecg = rspdivg;
    gathspec.nfld   = nvordiv;
    gathspec.nto    = nto;
  trans_gathspec(&gathspec);

  // Allocate fields for u*cos(theta) and v*cos(theta)
  double* rspug    = malloc( sizeof(double) * nvordiv*trans.nspec2g );
  double* rspvg    = malloc( sizeof(double) * nvordiv*trans.nspec2g );

  // Convert vorticity & divergence to u*cos(theta) & v*cos(theta)
  printf("Converting spectral vorticity-divergence to U-V globally...\n");
  struct VorDivToUV_t vordiv_to_UV_g = new_vordiv_to_UV();
    vordiv_to_UV_g.rspvor = rspvorg;
    vordiv_to_UV_g.rspdiv = rspdivg;
    vordiv_to_UV_g.rspu   = rspug;
    vordiv_to_UV_g.rspv   = rspvg;
    vordiv_to_UV_g.nfld   = nvordiv;
    vordiv_to_UV_g.ncoeff = trans.nspec2g;
    vordiv_to_UV_g.nsmax  = trans.nsmax;
  trans_vordiv_to_UV(&vordiv_to_UV_g);
  printf("Converting spectral vorticity-divergence to U-V globally...done\n");


  // Distribute spectral field (for fun)
  struct DistSpec_t distspec = new_distspec(&trans);
    distspec.rspec  = rspscalar;
    distspec.rspecg = rspscalarg;
    distspec.nfld   = nscalar;
    distspec.nfrom  = nto;
  trans_distspec(&distspec);

  // Inverse Transform
  struct InvTrans_t invtrans = new_invtrans(&trans);
    invtrans.nscalar   = nscalar;
    invtrans.nvordiv   = nvordiv;
    invtrans.rspscalar = rspscalar;
    invtrans.rspvor    = rspvor;
    invtrans.rspdiv    = rspdiv;
    invtrans.rgp       = rgp;
  trans_invtrans(&invtrans);

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<3; ++j)
    {
      for( i=0; i<nout; ++i )
      {
        printf("rgp[%d][%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
      }
    }
  }

  // Gather gridpoint fields
  struct GathGrid_t gathgrid = new_gathgrid(&trans);
    gathgrid.rgp  = rgp;
    gathgrid.rgpg = rgpg;
    gathgrid.nto  = nfrom;
    gathgrid.nfld = nfld;
  trans_gathgrid(&gathgrid);


  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<3; ++j)
    {
      for( i=0; i<nout; ++i )
      {
        printf("rgpg[%d][%d] : %f\n",j,i,rgpg[j*trans.ngptotg+i]);
      }
    }
  }

  // Deallocate arrays
  free(rgp);
  free(rgpg);
  free(rspscalar);
  free(rspscalarg);
  free(rspvor);
  free(rspvorg);
  free(rspdiv);
  free(rspdivg);
  free(rspu);
  free(rspug);
  free(rspv);
  free(rspvg);
  free(nfrom);
  free(nto);

  printf("cleanup");
  trans_delete(&trans);

  trans_finalize();

  //printf("transi finished\n");
  return 0;
}







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
