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


int read_bytes(const char* fp, void** buffer, size_t* size)
{
  FILE *fileptr;
  fileptr = fopen(fp, "rb");     // Open the file in binary mode
  fseek(fileptr, 0, SEEK_END);   // Jump to the end of the file
  *size = ftell(fileptr);        // Get the current byte offset in the file
  rewind(fileptr);               // Jump back to the beginning of the file

  *buffer = (void *)malloc((*size+1)*sizeof(char)); // Enough memory for file + \0
  fread(*buffer, *size, 1, fileptr); // Read in the entire file
  fclose(fileptr);
  return 0;
}


void test_io()
{
  int N=80;
  int T=-1;
  struct Trans_t trans;
  void* buffer;
  size_t size;
  const char* filepath = "TL159_lp";
  double start,end;
  int mem;

  TRANS_CHECK( trans_new(&trans) );
  set_standard_rgg(&trans,N,T);
  TRANS_CHECK( trans_set_write(&trans,filepath) );
  mem = allocated();
  start = transi_test_time();
  TRANS_CHECK( trans_setup(&trans) );
  end=transi_test_time();
  print_time("Timing rgg compute+write: ",end-start);
  print_mem( "Alloc  rgg compute+write: ",allocated()-mem);
  TRANS_CHECK( trans_delete(&trans) );


  read_bytes(filepath,&buffer,&size);
  print_mem( "Cache size:", size);

  TRANS_CHECK( trans_new(&trans) );
  set_standard_rgg(&trans,N,T);
  TRANS_CHECK( trans_set_cache(&trans,buffer,size) );
  mem = allocated();
  start = transi_test_time();
  TRANS_CHECK( trans_setup(&trans) );
  print_time("Timing rgg use cache: ",transi_test_time()-start);
  print_mem( "Alloc  rgg use cache: ",allocated()-mem);
  TRANS_CHECK( trans_delete(&trans) );


  TRANS_CHECK( trans_new(&trans) );
  set_standard_rgg(&trans,N,T);
  TRANS_CHECK( trans_set_read(&trans,filepath) );
  mem = allocated();
  start = transi_test_time();
  TRANS_CHECK( trans_setup(&trans) );
  print_time("Timing rgg read: ",transi_test_time()-start);
  print_mem( "Alloc  rgg read: ",allocated()-mem);
  TRANS_CHECK( trans_delete(&trans) );

  free(buffer);
}

void test_io_lonlat(int nlon, int nlat, int nsmax, int flt)
{
  printf("\nll.%dx%d --- T%d\n",nlon,nlat,nsmax);
  struct Trans_t trans;
  void* buffer = NULL;
  size_t size;
  char filepath[200];
  sprintf(filepath,"T%d_ll.%dx%d_flt%d",nsmax,nlon,nlat,flt);
  double start;
  int mem;

//  int nlon=1280;
//  int nlat=641;
//  int nsmax=639;

//  int nlon=640;
//  int nlat=321;
//  int nsmax=319;

//  int nlon=320;
//  int nlat=161;
//  int nsmax=159;

  bool lonlat=true;
  int N = (nlat-1)/2;

//  int nlon=160;
//  int nlat=81;
//  int nsmax=79;

  int nscalar = 2;
  int nvordiv = 1;
  int nfld = 2*nvordiv+nscalar;

  double* rspscalar = NULL;
  double* rspvor    = NULL;
  double* rspdiv    = NULL;
  double* rgp = NULL;
  struct InvTrans_t invtrans;

  // ---------------------------------------
  // Writing
  // ---------------------------------------
  int j;
  for( j=0; j<1; ++j ) {
  printf("Writing\n");

  TRANS_CHECK( trans_new(&trans) );
  trans.flt = flt;
  if( lonlat )
  {
    TRANS_CHECK( trans_set_resol_lonlat(&trans,nlon,nlat) );
    TRANS_CHECK( trans_set_trunc(&trans,nsmax) );
  }
  else
  {
    set_standard_rgg(&trans,N,nsmax);
  }
  TRANS_CHECK( trans_set_write(&trans,filepath) );
  mem = allocated();
  start = transi_test_time();
  TRANS_CHECK( trans_setup(&trans) );
  print_time("Timing lonlat compute+write: ",transi_test_time()-start);
  print_mem ("Alloc  lonlat compute+write: ",allocated()-mem);


  if( nscalar && rspscalar == NULL )
    rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
  if( nvordiv && rspvor == NULL && rspdiv == NULL ) {
    rspvor    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
    rspdiv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  }
  if( nfld && rgp == NULL )
    rgp  = malloc( sizeof(double) * nfld * trans.ngptot );

  printf("trans_invtrans()...\n");
  double start_transform = transi_test_time();
  invtrans = new_invtrans(&trans);
    invtrans.nscalar   = nscalar;
    invtrans.nvordiv   = nvordiv;
    invtrans.rspscalar = rspscalar;
    invtrans.rspvor    = rspvor;
    invtrans.rspdiv    = rspdiv;
    invtrans.rgp       = rgp;
  TRANS_CHECK( trans_invtrans(&invtrans) );
  print_time("trans_invtrans()...done in ",transi_test_time()-start_transform);

  TRANS_CHECK( trans_delete(&trans) );
  }

  // ---------------------------------------
  // Use Cache
  // ---------------------------------------
  if(1) {
  printf("Use Cache\n");

  read_bytes(filepath,&buffer,&size);

  print_mem( "Cache size:", size);

  TRANS_CHECK( trans_new(&trans) );
  trans.flt = flt;
  if( lonlat )
  {
    TRANS_CHECK( trans_set_resol_lonlat(&trans,nlon,nlat) );
    TRANS_CHECK( trans_set_trunc(&trans,nsmax) );
  }
  else
  {
    set_standard_rgg(&trans,N,nsmax);
  }

  TRANS_CHECK( trans_set_cache(&trans,buffer,size) );
  mem = allocated();
  start = transi_test_time();
  TRANS_CHECK( trans_setup(&trans) );
  print_time("Timing lonlat use cache: ",transi_test_time()-start);
  print_mem ("Alloc  lonlat use cache: ",allocated()-mem);

  if( nscalar && rspscalar == NULL )
    rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
  if( nvordiv && rspvor == NULL && rspdiv == NULL ) {
    rspvor    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
    rspdiv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  }
  if( nfld && rgp == NULL )
    rgp  = malloc( sizeof(double) * nfld * trans.ngptot );

  printf("trans_invtrans()...\n");
  double start_transform = transi_test_time();
  invtrans = new_invtrans(&trans);
    invtrans.nscalar   = nscalar;
    invtrans.nvordiv   = nvordiv;
    invtrans.rspscalar = rspscalar;
    invtrans.rspvor    = rspvor;
    invtrans.rspdiv    = rspdiv;
    invtrans.rgp       = rgp;
  TRANS_CHECK( trans_invtrans(&invtrans) );
  print_time("trans_invtrans()...done in ",transi_test_time()-start_transform);


  TRANS_CHECK( trans_delete(&trans) );
  }

  // ---------------------------------------
  // Reading
  // ---------------------------------------
  printf("Reading\n");
  TRANS_CHECK( trans_new(&trans) );
  trans.flt = flt;
  if( lonlat )
  {
    TRANS_CHECK( trans_set_resol_lonlat(&trans,nlon,nlat) );
    TRANS_CHECK( trans_set_trunc(&trans,nsmax) );
  }
  else
  {
    set_standard_rgg(&trans,N,nsmax);
  }
  TRANS_CHECK( trans_set_read(&trans,filepath) );
  mem = allocated();
  start = transi_test_time();
  TRANS_CHECK( trans_setup(&trans) );
  print_time("Timing lonlat read: ",transi_test_time()-start);
  print_mem ("Alloc  lonlat read: ",allocated()-mem);

  if( nscalar && rspscalar == NULL )
    rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
  if( nvordiv && rspvor == NULL && rspdiv == NULL ) {
    rspvor    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
    rspdiv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  }
  if( nfld && rgp == NULL )
    rgp  = malloc( sizeof(double) * nfld * trans.ngptot );

  printf("trans_invtrans()...\n");
  invtrans = new_invtrans(&trans);
    invtrans.nscalar   = nscalar;
    invtrans.nvordiv   = nvordiv;
    invtrans.rspscalar = rspscalar;
    invtrans.rspvor    = rspvor;
    invtrans.rspdiv    = rspdiv;
    invtrans.rgp       = rgp;
  TRANS_CHECK( trans_invtrans(&invtrans) );
  printf("trans_invtrans()...done\n");

  TRANS_CHECK( trans_delete(&trans) );

  // ---------------------------------------


  if( buffer ) free(buffer);
  if( rgp ) free(rgp);
  if( rspscalar ) free(rspscalar);
  if( rspvor ) free(rspvor);
  if( rspdiv ) free(rspdiv);

}

int main ( int arc, char **argv )
{
  trans_use_mpi( test_use_mpi() );
  setbuf(stdout, NULL);
  TRANS_CHECK( trans_init() );

  test_io();

  int flt = false;
  test_io_lonlat(320,161,511,flt);
  test_io_lonlat(320,161,159,flt);
//  test_io_lonlat(2400,1201,799,flt);


  TRANS_CHECK( trans_finalize() );

  return 0;
}

