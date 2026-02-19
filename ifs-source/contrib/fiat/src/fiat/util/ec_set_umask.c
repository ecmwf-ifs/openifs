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

/* ec_set_umask.c */

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>

void ec_set_umask_(void)
{
  char *env = getenv("EC_SET_UMASK");
  if (env) {
    int newmask;
    int n = sscanf(env,"%o",&newmask);
    if (n == 1) {
      int oldmask = umask(newmask);
      fprintf(stderr,
	      "*** EC_SET_UMASK : new/old = %o/%o (oct), %d/%d (dec), %x/%x (hex)\n",
	      newmask,oldmask,
	      newmask,oldmask,
	      newmask,oldmask);
    } /* if (n == 1) */
  } /* if (env) */
}

