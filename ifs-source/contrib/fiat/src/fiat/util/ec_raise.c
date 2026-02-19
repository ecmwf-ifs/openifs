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

/* ec_raise.c */

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>

/* CALL ec_raise(6) == CALL abort() */

void ec_raise_(const int *sig) { raise(*sig); }
void ec_raise(const int *sig) { ec_raise_(sig); }

