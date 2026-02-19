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
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


void save_ (const char * p, const char * f, const int * p_len, const unsigned long int f_len)
{
  char _f[f_len+1];
  FILE * fp;

  memset (_f, 0, sizeof (_f));
  memcpy (_f, f, f_len);

  fp = fopen (_f, "w");
  fwrite (p, 1, *p_len, fp);
  fclose (fp);
}

void load_ (char * p, const char * f, int * p_len, const unsigned long int f_len)
{
  char _f[f_len+1];
  FILE * fp;
  struct stat st;
  int size;

  memset (_f, 0, sizeof (_f));
  memcpy (_f, f, f_len);

  stat (_f, &st);
  size = st.st_size;

  fp = fopen (_f, "r");
  fread (p, 1, size, fp);
  fclose (fp);

  *p_len = size;

}
