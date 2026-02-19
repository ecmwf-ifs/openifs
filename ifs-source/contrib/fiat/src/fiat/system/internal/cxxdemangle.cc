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
#include <stdlib.h>
#include <cxxabi.h>

extern "C" 
char *cxxdemangle(const char *mangled_name, int *status)
{
  int istat = 0;
#ifdef __NEC__
  // Where is libstdc++ on NEC Aurora ??
  char *demangled_name = NULL;
#else
  char *demangled_name = abi::__cxa_demangle(mangled_name, NULL, NULL, &istat);
#endif
  if (status) *status = istat;
  return demangled_name; // this must be free()'d by the user
}
