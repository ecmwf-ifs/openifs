/*
 * (C) Copyright 1989- ECMWF.
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * 
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction
 * 
 * (C) Copyright 1989- Meteo-France.
 * 
*/

#include <stdio.h>

/*
 * Used to get a string representation of an address from FORTRAN.
 * Useful to generate unique ids.
 */

void paddrs_ (char * s, void * p, int slen)
{
  int i, n = snprintf (s, slen, "%p", p);
  for (i = n; i < slen; i++) s[i] = ' ';
}
