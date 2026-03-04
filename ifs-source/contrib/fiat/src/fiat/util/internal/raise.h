/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#ifndef _RAISE_H_
#define _RAISE_H_

/* raise.h */

#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>

#include "abor1.h"

#define RAISE(x) { \
  if ((x) == SIGABRT) { \
    ABOR1("*** Fatal error; aborting (SIGABRT) ..."); \
    _exit(1); /* Should never end up here */ \
  } \
  else raise(x); \
}

#endif /* _RAISE_H_ */
