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

#include <string.h>

static const char * T = "abcd";
static const int dTle = 0x64636261;
static const int dTbe = 0x61626364;

void iswap_isle_ (int * reqd)
{
  int * dT = (int *)T;

  if (*dT == dTle)
    *reqd = 1;
  else
    *reqd = 0;
}

void iswap (char * a, const char * b, int t, int n, int d)
{
  long i, j;

  if (d)
    {
      for (i = 0; i < n; i++)
        for (j = 0; j < t / 2; j++)
          {
            char c = b[i*t+j]; /* a and b may be the same */
            a[i*t+j] = b[i*t+t-j-1];
            a[i*t+t-j-1] = c;
          }
    }
  else if (a != b)
    {
      size_t nbytes = t * n;
      memcpy (a, b, nbytes);
    }

}

void iswap_ (char * a, const char * b, const int * _t, const int * _n)
{
  int * dT = (int *)T;
  iswap (a, b, *_t, *_n, *dT == dTle);
}

void jswap_ (char * a, const char * b, const int * _t, const int * _n)
{
  iswap (a, b, *_t, *_n, 1);
}


