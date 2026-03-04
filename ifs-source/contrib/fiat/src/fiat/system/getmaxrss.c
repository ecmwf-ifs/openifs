/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

typedef  long long int  ll_t;
#include <sys/resource.h>
#include <sys/time.h>
ll_t
getmaxrss_()
{
  const ll_t scaler = 1024; /* in kilobytes */
  ll_t rc = 0;
  struct rusage r;
  rc = getrusage(RUSAGE_SELF, &r);
  rc = (rc == 0) ? (ll_t) r.ru_maxrss * scaler : 0;
  return rc;
}
