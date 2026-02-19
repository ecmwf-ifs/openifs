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

/* crc.h */
/* Calculates 32-bit  Cyclic Redundancy Check as in Unix cksum command */
/* Also calculates 64-bit  Cyclic Redundancy Check */

/* Sami Saarinen, 17-Feb-2005 : crc32       */
/*                24-Jun-2005 : Added crc64 */
/*                29-Dec-2005 : Added protos for direct C-calls */
/*                29-Dec-2005 : Length argument for crc64 now 64-bit int */

/* C callables */

extern unsigned int 
pp_cksum32(int nbuf, unsigned int nCRC);

extern unsigned int 
pp_cksum32but64len(long long int nbuf, unsigned int nCRC);

extern unsigned long long int 
pp_cksum64(long long int nbuf, unsigned long long int nCRC);

extern unsigned int 
cksum32(const char *buf, int nbuf, unsigned int nCRC);

extern unsigned long long int
cksum64(const char *buf, long long int nbuf, unsigned long long int nCRC);

/* Fortran callable */

extern void 
crc32_(const void *vbuf, const int *pnbuf, 
       unsigned int *pnCRC /* Note: An in & out -variable */);

extern void 
crc64_(const void *vbuf, const long long int *pnbuf, 
       unsigned long long int *pnCRC /* Note: An in & out -variable */);
