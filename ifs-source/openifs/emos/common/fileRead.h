/**
* Copyright 1981- ECMWF.
*
* This software is licensed under the terms of the Apache Licence
* Version 2.0 which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*
* In applying this licence, ECMWF does not waive the privileges and immunities
* granted to it by virtue of its status as an intergovernmental organisation
* nor does it submit to any jurisdiction.
*/
#ifndef FILE_READ_H
#define FILE_READ_H

#include <stdio.h>

#if defined(linux) && !defined(darwin)
# if !defined __off64_t_defined
typedef __off64_t off64_t;
#define __off64_t_defined
#endif
#endif

#ifdef FOPEN64
#define OFF_T off64_t
#else
#define OFF_T off_t
#endif

static fortint fileRead(char * buffer, fortint length, void * file) {
/*
//  buffer = buffer to fill,
//  length = size of buffer in bytes,
//  file   = file pointer from fopen().
//
//  Returns the number of bytes read.
//  On END_OF_FILE, returns negative value for number of bytes read .
*/

fortint nbytes;

  nbytes = (fortint) fread( buffer, 1, length, (FILE *) file);
  if ( feof((FILE *) file ) ) {
     clearerr( (FILE *) file );
     return (-nbytes);
  }
  return nbytes;
}

#ifdef FOPEN64
static OFF_T fileSeek(void * file, OFF_T offset, fortint from) {

  return ( (OFF_T) fseeko64((FILE *) file, offset, from) );
}
#else
static fortint fileSeek(void * file, fortint offset, fortint from) {

  return ( (fortint) fseek((FILE *) file, offset, from) );
}
#endif


#ifdef FOPEN64
static OFF_T fileTell(void * file) {

  return ( (OFF_T) ftello64((FILE *) file) );
}
#else
static fortint fileTell(void * file) {

  return ( (fortint) ftell((FILE *) file) );
}
#endif

#ifdef FOPEN64
fortint readprod( char *, char * , fortint * ,
                  fortint (*fileRead)(char *, fortint, void *),
                  OFF_T (*fileSeek)(void *, OFF_T, fortint),
                  OFF_T (*fileTell)(void *),
                  void * );
#else
fortint readprod( char *, char * , fortint * ,
                  fortint (*fileRead)(char *, fortint, void *),
                  fortint (*fileSeek)(void *, fortint, fortint),
                  fortint (*fileTell)(void *),
                  void * );
#endif

#endif /* end of FILE_READ_H */
