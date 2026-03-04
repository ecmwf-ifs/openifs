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
#include <time.h>
void cgmtime_(int *gmt )
{
  time_t rawtime;
  struct tm *info;
  time(&rawtime);
  info=gmtime(&rawtime);
  gmt[0]=info->tm_sec;
  gmt[1]=info->tm_min;
  gmt[2]=info->tm_hour;
  gmt[3]=info->tm_mday;
  gmt[4]=info->tm_mon;
  gmt[5]=info->tm_year;
  gmt[6]=info->tm_wday;
  gmt[7]=info->tm_yday;
  gmt[8]=info->tm_isdst;
  return;
}
