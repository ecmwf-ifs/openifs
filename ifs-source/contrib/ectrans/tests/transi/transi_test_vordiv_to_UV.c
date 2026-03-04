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

// // Following dummy functions are implementation details
// // that don't contribute to this example. They could be
// // replaced with grib_api functionality
// void read_rspec( struct Trans_t* trans, double* rspec[], int nfld );

int main ( int arc, char **argv )
{
  trans_use_mpi( test_use_mpi() );

  int nfld = 1;
  int nsmax_array[] = {3};//159,160,1279,1280};
  int ni = sizeof(nsmax_array)/sizeof(int);
  int i,j;
  for( i=0; i<ni; ++i)
  {
    int nsmax = nsmax_array[i];
    int nspec2 = (nsmax+1)*(nsmax+2);

    printf( "nsmax = %d\n", nsmax );

    double* field_vor = malloc( sizeof(double) * nfld*nspec2 );
    double* field_div = malloc( sizeof(double) * nfld*nspec2 );
    double* field_U   = malloc( sizeof(double) * nfld*nspec2 );
    double* field_V   = malloc( sizeof(double) * nfld*nspec2 );


    for( j=0; j<nfld*nspec2; ++j ) {
      field_vor[j] = 0;
      field_div[j] = 0;
    }


    struct VorDivToUV_t vod_to_UV = new_vordiv_to_UV();
    vod_to_UV.nfld   = 1;              // number of distributed fields
    vod_to_UV.ncoeff = nspec2;         // number of spectral coefficients (equivalent to NSPEC2 for distributed or NSPEC2G for global)
    vod_to_UV.nsmax  = nsmax;          // spectral resolution (T)
    vod_to_UV.rspvor = field_vor;      // spectral array for vorticity    DIMENSIONS(1:NFLD,1:NSPEC2)
    vod_to_UV.rspdiv = field_div;      // spectral array for divergence   DIMENSIONS(1:NFLD,1:NSPEC2)
    vod_to_UV.rspu   = field_U;        // spectral array for u*cos(theta) DIMENSIONS(1:NFLD,1:NSPEC2)
    vod_to_UV.rspv   = field_V;        // spectral array for v*cos(theta) DIMENSIONS(1:NFLD,1:NSPEC2)

    TRANS_CHECK( trans_vordiv_to_UV(&vod_to_UV) );

    free(field_vor);
    free(field_div);
    free(field_U);
    free(field_V);
  }

  return 0;
}
