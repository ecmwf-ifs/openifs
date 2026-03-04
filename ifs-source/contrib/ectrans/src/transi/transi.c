/*
 * (C) Copyright 2014- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/*
 *  @file  transi/trans.c
 *  @brief C-interface to the IFS trans-library
 *  @author Willem Deconinck (nawd)
 *  @date Jul 2014
 */

#include <stdlib.h>
#include <string.h>
#include "transi.h"

/*
 * These functions are to be used in the fortran part (trans_module.F90) to
 * allocate and deallocate arrays of type C_PTR
 */
void transi_malloc_bool  (void* ptr[], int len) { *ptr = malloc(sizeof(bool  ) * len); }
void transi_malloc_int   (void* ptr[], int len) { *ptr = malloc(sizeof(int   ) * len); }
void transi_malloc_float (void* ptr[], int len) { *ptr = malloc(sizeof(float ) * len); }
void transi_malloc_double(void* ptr[], int len) { *ptr = malloc(sizeof(double) * len); }
void transi_free(void* ptr[]) { free(*ptr); *ptr=NULL; }

#define TRANS_ERROR            -1
#define TRANS_NOTIMPL          -2
#define TRANS_MISSING_ARG      -3
#define TRANS_UNRECOGNIZED_ARG -4
#define TRANS_STALE_ARG        -5

const char* trans_error_msg(int errcode)
{
  switch( errcode )
  {
    case TRANS_SUCCESS:
  return "No error";
    case TRANS_ERROR:
  return "Trans: Error";
    case TRANS_NOTIMPL:
  return "Trans: Not (yet) implemented";
    case TRANS_MISSING_ARG:
  return "Trans: Required member of the argument structure is missing or not allocated";
    case TRANS_UNRECOGNIZED_ARG:
  return "Trans: Unrecognized argument";
    case TRANS_STALE_ARG:
  return "Trans: Passed argument was already used in previous call";
    default:
  return "Trans: Unknown error";
  }
}

int trans_new( struct Trans_t* trans )
{
  trans->llatlon = 0;
  trans->lsplit = true;
  trans->flt = -1;
  trans->fft = TRANS_FFTW;
  trans->nsmax = -1;
  trans->ndgl = -1;
  trans->nlon = -1;
  trans->nloen = NULL;
  trans->readfp = NULL;
  trans->writefp = NULL;
  trans->cache = NULL;
  trans->cachesize = 0;
  return TRANS_SUCCESS;
}

int trans_set_resol( struct Trans_t* trans, int ndgl, const int* nloen )
{
  size_t i;
  trans->ndgl = ndgl;
  trans->nloen = malloc( sizeof(int) * ndgl );
  for ( i = 0; i < ndgl; ++i )
    trans->nloen[i] = nloen[i];
  return TRANS_SUCCESS;
}

int trans_set_resol_lonlat( struct Trans_t* trans, int nlon, int nlat )
{
  size_t i;
  if( nlat%2 == 0 ) // The shifted lonlat grid (excluding poles and equator)
  {
    trans->ndgl = nlat;
    trans->nlon = nlon;
    trans->llatlon = 2;
    if( trans->nloen ) free(trans->nloen);
    trans->nloen = malloc( sizeof(int) * nlat );
    for ( i = 0; i < nlat; ++i )
      trans->nloen[i] = nlon;
  }
  else // The lonlat grid including poles and equator
  {
    trans->ndgl = nlat-1; // Internally coefficients are computed with ndgl+2 (equator duplicated)
    trans->nlon = nlon;
    trans->llatlon = 1;
  }
  return TRANS_SUCCESS;
}

int trans_set_trunc( struct Trans_t* trans, int nsmax )
{
  trans->nsmax = nsmax;
  return TRANS_SUCCESS;
}

int trans_set_read(struct Trans_t* trans, const char* filepath)
{
  trans->readfp = malloc( sizeof(char)*1024 );
  strcpy(trans->readfp, filepath);
  return TRANS_SUCCESS;
}

int trans_set_write(struct Trans_t* trans, const char* filepath)
{
  trans->writefp = malloc( sizeof(char)*1024 );
  strcpy(trans->writefp, filepath);
  return TRANS_SUCCESS;
}

int trans_set_cache(struct Trans_t* trans, const void* cache , size_t cachesize)
{
  trans->cache = cache;
  trans->cachesize = cachesize;
  return TRANS_SUCCESS;
}

struct DirTrans_t new_dirtrans(struct Trans_t* trans)
{
  struct DirTrans_t dirtrans;
  dirtrans.count = 0;
  dirtrans.rgp = NULL;
  dirtrans.rspscalar = NULL;
  dirtrans.rspvor = NULL;
  dirtrans.rspdiv = NULL;
  dirtrans.ngpblks = 1;
  dirtrans.nproma = trans->ngptot;
  dirtrans.nscalar = 0;
  dirtrans.nvordiv = 0;
  dirtrans.lglobal = 0;
  dirtrans.trans = trans;
  return dirtrans;
}

struct InvTrans_t new_invtrans(struct Trans_t* trans)
{
  struct InvTrans_t invtrans;
  invtrans.count = 0;
  invtrans.rspscalar = NULL;
  invtrans.rspvor = NULL;
  invtrans.rspdiv = NULL;
  invtrans.rgp = NULL;
  invtrans.ngpblks = 1;
  invtrans.nproma = trans->ngptot;
  invtrans.nscalar = 0;
  invtrans.nvordiv = 0;
  invtrans.lscalarders = 0;
  invtrans.luvder_EW = 0;
  invtrans.lvordivgp = 0;
  invtrans.lglobal = 0;
  invtrans.trans = trans;
  return invtrans;
}

struct InvTransAdj_t new_invtrans_adj(struct Trans_t* trans)
{
  struct InvTransAdj_t invtrans_adj;
  invtrans_adj.count = 0;
  invtrans_adj.rspscalar = NULL;
  invtrans_adj.rspvor = NULL;
  invtrans_adj.rspdiv = NULL;
  invtrans_adj.rgp = NULL;
  invtrans_adj.ngpblks = 1;
  invtrans_adj.nproma = trans->ngptot;
  invtrans_adj.nscalar = 0;
  invtrans_adj.nvordiv = 0;
  invtrans_adj.lscalarders = 0;
  invtrans_adj.luvder_EW = 0;
  invtrans_adj.lvordivgp = 0;
  invtrans_adj.lglobal = 0;
  invtrans_adj.trans = trans;
  return invtrans_adj;
}


struct DistGrid_t new_distgrid(struct Trans_t* trans)
{
  struct DistGrid_t distgrid;
  distgrid.count = 0;
  distgrid.rgpg = NULL;
  distgrid.rgp = NULL;
  distgrid.nfrom = NULL;
  distgrid.ngpblks = 1;
  distgrid.nproma = trans->ngptot;
  distgrid.nfld = 0;
  distgrid.trans = trans;
  return distgrid;
}

struct GathGrid_t new_gathgrid(struct Trans_t* trans)
{
  struct GathGrid_t gathgrid;
  gathgrid.count = 0;
  gathgrid.rgpg = NULL;
  gathgrid.rgp = NULL;
  gathgrid.nto = NULL;
  gathgrid.ngpblks = 1;
  gathgrid.nproma = trans->ngptot;
  gathgrid.nfld = 0;
  gathgrid.trans = trans;
  return gathgrid;
}

struct DistSpec_t new_distspec(struct Trans_t* trans)
{
  struct DistSpec_t distspec;
  distspec.count = 0;
  distspec.rspecg = NULL;
  distspec.rspec = NULL;
  distspec.nfrom = NULL;
  distspec.trans = trans;
  return distspec;
}

struct GathSpec_t new_gathspec(struct Trans_t* trans)
{
  struct GathSpec_t gathspec;
  gathspec.count = 0;
  gathspec.rspecg = NULL;
  gathspec.rspec = NULL;
  gathspec.nto = NULL;
  gathspec.trans = trans;
  return gathspec;
}

struct VorDivToUV_t new_vordiv_to_UV()
{
  struct VorDivToUV_t vdtouv;
  vdtouv.count = 0;
  vdtouv.rspvor = NULL;
  vdtouv.rspdiv = NULL;
  vdtouv.rspu = NULL;
  vdtouv.rspv = NULL;
  vdtouv.nfld = 0;
  vdtouv.ncoeff = 0;
  vdtouv.nsmax = 0;
  return vdtouv;
}

struct SpecNorm_t new_specnorm(struct Trans_t* trans)
{
  struct SpecNorm_t specnorm;
  specnorm.rspec = NULL;
  specnorm.nmaster = 1;
  specnorm.rmet = NULL;
  specnorm.rnorm = NULL;
  specnorm.nfld = 0;
  specnorm.trans = trans;
  specnorm.count = 0;
  return specnorm;
}
