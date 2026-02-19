/*
 * (C) Copyright 2014- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#ifndef fiat_version_h
#define fiat_version_h

#ifdef __cplusplus
extern "C" {
#endif

const char * fiat_version();

unsigned int fiat_version_int();

const char * fiat_version_str();

const char * fiat_git_sha1();

const char * fiat_git_sha1_abbrev(unsigned int length);

#ifdef __cplusplus
}
#endif

#endif
