/*
! (C) Copyright 1989- ECMWF.
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * 
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction
 * 
! (C) Copyright 1989- Meteo-France.
 * 
*/

#ifndef _FI_LIBC_H
#define _FI_LIBC_H

#include <stdio.h>

typedef int fi_integer4;
typedef long long int fi_integer8;

void fi_fopen_ (FILE ** fpp, const char * path, const char * mode, int path_len, int mode_len);
void fi_fclose_ (fi_integer4 * err, FILE ** fpp);
void fi_fread_ (fi_integer8 * err, void * ptr, fi_integer8 * size, fi_integer8 * nmemb, FILE ** fpp);
void fi_fwrite_ (fi_integer8 * err, const void * ptr, fi_integer8 * size, fi_integer8 * nmemb, FILE ** fpp);
void fi_fseek_ (fi_integer4 * err, FILE ** fpp, fi_integer8 * offset, fi_integer4 * whence);
void fi_ftell_ (fi_integer8 * err, FILE ** fpp);
void fi_fileno_ (fi_integer4 * err, FILE ** fpp);
void fi_fstat_ (fi_integer4 *err, fi_integer4 *fd, fi_integer8 buf[13]);
void fi_ftruncate_ (fi_integer4* err, fi_integer4 *fd, fi_integer8 *length);
void fi_mkdir_ (const char * path, int path_len);
void fi_chdir_ (fi_integer4 * ierr, const char * path, int path_len);
void fi_fflush_ (fi_integer4* err, FILE **fpp);

#endif
