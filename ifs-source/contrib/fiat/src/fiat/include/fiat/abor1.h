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

/* abor1.h */

#ifndef _ABOR1_H_
#define _ABOR1_H_

#ifdef __cplusplus
extern "C" {
#endif
void abor1(const char* filename, const int linenum, const char* s);

void set_abor1_exception_handler();

#define ABOR1(txt) abor1( __FILE__, __LINE__, (txt) ) 

#ifdef __cplusplus
} // extern "C"
#endif
#endif /* _ABOR1_H_ */
