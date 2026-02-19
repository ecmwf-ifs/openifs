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
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#ifdef RS6K
#include <dlfcn.h>
#endif
#include <sys/time.h>

long long int jio_time_()
{
  long long int    elapsed;
#ifdef RS6K
  timebasestruct_t t_start;
  int              sec, n_sec;
  read_real_time(&t_start, TIMEBASE_SZ);
  time_base_to_time(&t_start, TIMEBASE_SZ);
  sec   = t_start.tb_high;
  n_sec = t_start.tb_low;
  elapsed = (long long int) sec*1000000000 + (long long int)n_sec;
#else
  elapsed = 0;
#endif
  return elapsed;
}
