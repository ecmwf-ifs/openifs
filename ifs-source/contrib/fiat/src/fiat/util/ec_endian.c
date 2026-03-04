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

/* endian.c */

/* 
   Please note: the following 2 routines
   cannot be named as "is_little_endian()"
   and "is_big_endian()", since there is a clash
   with the new Magics++ library  (by SS, 21-Mar-2006)

   --> consequently "ec_" prefix was added
*/

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>

int ec_is_little_endian()
{
  /* Little/big-endian runtime auto-detection */
  const unsigned int ulbtest = 0x12345678;
  const unsigned char *clbtest = (const unsigned char *)&ulbtest;
  
  if (*clbtest == 0x78) { 
    /* We are on a little-endian machine */
    return 1;
  }
  else { 
    /* We are on a big-endian machine */
    return 0;
  }
}

int ec_is_big_endian()
{
  return !ec_is_little_endian();
}

/* Fortran interface */

int ec_is_big_endian_()    { return ec_is_big_endian();    }
int ec_is_little_endian_() { return ec_is_little_endian(); }

