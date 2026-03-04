/*
 * (C) Copyright 2005- ECMWF.
 * (C) Copyright 2013- Meteo-France.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <string.h>
#include <unistd.h>

extern void abor1fl_(const char *filename, const int *linenum, const char *s, 
                     int filenamelen, int slen);
extern void abor1_(const char *s, int slen);

void abor1(const char* filename, const int linenum, const char* s) {
  if( filename ) {
    s ? abor1fl_( filename, &linenum,  s, strlen(filename), strlen(s) )
      : abor1fl_( filename, &linenum, "", strlen(filename), 0 );
  }
  else {
    s ? abor1_( s, strlen(s) )
      : abor1_( "", 0 );
  }
  _exit(1); /* Should never end up here */
}
