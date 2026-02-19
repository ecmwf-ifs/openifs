/**
* Copyright 1981- ECMWF.
*
* This software is licensed under the terms of the Apache Licence
* Version 2.0 which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*
* In applying this licence, ECMWF does not waive the privileges and immunities
* granted to it by virtue of its status as an intergovernmental organisation
* nor does it submit to any jurisdiction.
*/

#ifndef FORTINT_H
#define FORTINT_H

#ifdef INTEGER_IS_LONG
#define fortint long
#else
#ifdef INTEGER_IS_INT
#define fortint int
#else
#if defined hpR64 || defined hpiaR64
#define fortint long long
#else
#define fortint long
#endif
#endif
#endif

#endif /* end of FORTINT_H */
