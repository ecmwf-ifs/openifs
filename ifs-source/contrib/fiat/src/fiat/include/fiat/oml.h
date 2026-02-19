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

/* oml.h */
#ifndef _OML_H_
#define _OML_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef long long int oml_lock_t; /* i.e. 64-bit integer */

int  oml_my_thread();         // Equivalent to omp_get_thread_num() + 1
int  oml_get_max_threads();   // Equivalent to omp_get_max_threads()
int  oml_get_num_threads();   // Equivalent to omp_get_num_threads()
int  oml_in_parallel();

void oml_set_debug(int);
int  oml_get_debug();
void oml_init_lock();
void oml_init_lockid(oml_lock_t*);
void oml_init_lockid_with_name(oml_lock_t*, const char* name);
void oml_set_lock();
void oml_set_lockid(oml_lock_t*);
void oml_unset_lock();
void oml_unset_lockid(oml_lock_t*);
int  oml_test_lock();
int  oml_test_lockid(oml_lock_t*);

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* _OML_H_ */
