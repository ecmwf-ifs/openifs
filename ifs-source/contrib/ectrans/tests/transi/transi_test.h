/*
 * (C) Copyright 2014- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#ifndef TRANSI_TEST_H
#define TRANSI_TEST_H

#include <stdio.h>
#include <stdlib.h>

#include "ectrans/transi.h"

#define TRANS_CHECK( CALL ) do {\
  int errcode = CALL;\
  if( errcode != TRANS_SUCCESS) {\
    printf("ERROR: %s failed @%s:%d:\n%s\n",#CALL,__FILE__,__LINE__,trans_error_msg(errcode));\
    exit(1);\
  }\
} while(0)

#define TRANS_CHECK_ERROR( CALL, ERR ) do {\
  int errcode = CALL;\
  if( errcode != ERR ) {\
    printf("ERROR: %s should fail with errcode %s(%d) @%s:%d:\n%s\n",#CALL,#ERR,ERR,__FILE__,__LINE__,trans_error_msg(errcode));\
    exit(1);\
  }\
} while(0)

#define ASSERT( assertion ) do {\
  if( !(assertion) ) {\
    printf("ERROR: Assertion `%s' failed @%s:%d\n",#assertion,__FILE__,__LINE__);\
    exit(1);\
  }\
} while(0)

#define TRANS_ERROR            -1
#define TRANS_NOTIMPL          -2
#define TRANS_MISSING_ARG      -3
#define TRANS_UNRECOGNIZED_ARG -4
#define TRANS_STALE_ARG        -5

double transi_test_time();
void print_time(const char* str,double elapsed);
void display_mallinfo(void);
int allocated();
void print_mem(const char* str,const int bytes);
void set_standard_rgg(struct Trans_t* trans, int N, int T);

int test_use_mpi();

#endif
