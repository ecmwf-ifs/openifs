/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdio.h>
#include <string.h>

extern int MPI_Get_version(int *version, int *subversion);
extern int MPI_Get_library_version(char *version, int *resultlen);

/* Depending on MPI version these function may not exist -- thus weak definition */
#pragma weak MPI_Get_version
#pragma weak MPI_Get_library_version

void ecmpi_version_(int *version,
		    int *subversion,
		    char *library_version,
		    int *resultlen
		    /* hidden length */
		    ,const int len_library_version)
{
  int slen = 0;
  if (version && subversion) {
    if (MPI_Get_version) {
      (void) MPI_Get_version(version, subversion);
    }
    else {
      *version = 0;
      *subversion = 0;
    }
  }
  if (library_version && len_library_version > 0) {
    if (MPI_Get_library_version) {
      char s[4096];
      (void) MPI_Get_library_version(s,&slen);
      if (slen > len_library_version) slen = len_library_version;
      while (slen > 0 && (s[slen-1] == '\n' || s[slen-1] == '\0')) slen--;
      memset(library_version,' ',len_library_version);
      if (slen > 0) memcpy(library_version,s,slen);
    }
  }
  if (resultlen) *resultlen = slen;
}

void ecmpi_version(int *version,
		   int *subversion,
		   char *library_version,
		   int *resultlen
		   /* hidden length */
		   ,const int len_library_version)
{
  ecmpi_version_(version,subversion,library_version,resultlen,len_library_version);
}

