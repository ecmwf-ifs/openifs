/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/* ecmwf_transfer.c  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Used in module strhandler (stransfer) */

void
ecmwf_transfer_(void *out, const int *Len_out,
		const void *in, const int *Len_in
		/* Possible hidden argument (not referred) */
		, int Sta_lin)
{
  int len = *Len_out;
  if (*Len_in < len) len = *Len_in;
  if (len > 0) memcpy(out,in,len);
}
