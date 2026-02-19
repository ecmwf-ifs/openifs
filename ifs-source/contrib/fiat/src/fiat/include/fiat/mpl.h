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

/* mpl.h */
#ifndef _MPL_H_
#define _MPL_H_

#ifdef __cplusplus
extern "C" {
#endif

int  mpl_init();
int  mpl_end();
int  mpl_myrank(); // Note return value is 1-based as opposed to MPI_Rank which is 0-based

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* _MPL_H_ */
