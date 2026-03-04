/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#if defined(__GNUC__)

#include <unistd.h>

int getpid_() { /* GNU Fortran did not recognize this ? Here it comes then */
  return (int)getpid();
}

#endif /* defined(__GNUC__) */
