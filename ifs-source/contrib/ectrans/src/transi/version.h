/*
 * (C) Copyright 2014- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#ifndef ectrans_version_h
#define ectrans_version_h

#ifndef __cplusplus
// C99 header, defines bool as _Bool ( only required for C compiler )
#include <stdbool.h>
#else
extern "C" {
#endif

const char * ectrans_version();

unsigned int ectrans_version_int();

const char * ectrans_version_str();

const char * ectrans_git_sha1();

const char * ectrans_git_sha1_abbrev(unsigned int length);

const char * ectrans_fiat_version();

unsigned int ectrans_fiat_version_int();

const char * ectrans_fiat_version_str();

const char * ectrans_fiat_git_sha1();

const char * ectrans_fiat_git_sha1_abbrev(unsigned int length);

#ifdef __cplusplus
}
#endif

#endif
