/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/* fnec.c */

#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef __NEC__
void exit_(const int *exit_code) {
  exit(exit_code ? *exit_code : 0 );
}

void getenv_(const char *s, char *value,
             /* Hidden arguments */
             int slen,
             const int valuelen)
{
  extern void ec_getenv(const char *s, char *value, /* Hidden arguments */ int slen, const int valuelen);
  ec_getenv(s, value, slen, valuelen);
}

void sleep_(const int *nsec) {
	extern unsigned int ec_sleep(const int nsec);
	(void)ec_sleep(*nsec);
}

void flush_(const int *io) {
	/* temporary fix */
}

int hostnm_(char a[], int alen) {
	extern void ec_gethostname(char a[],  /* Hidden argument */ int alen);
	ec_gethostname(a,alen);
	return 0;
}

void system_(const char *name, /* Hidden arguments */ int length) {
  char string[4096];
  int return_code = 0;

  strncpy(string, name, length);
  while(length && string[length-1] == ' ')length --;
  string[length]=0;

  return_code = system(string);
  if ( return_code != 0 ) {
    if ( errno != 0 )
      return_code = -errno;
    else
      return_code = -return_code;
  }

  return(return_code);
}

#endif
