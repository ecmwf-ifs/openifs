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
#ifndef FORT2C_H
#define FORT2C_H

#ifdef CRAY
/* #include <fortran.h> */
typedef int _f_int6;
typedef long _f_int8;
typedef _f_int6 _f_int;
typedef _f_int8 _f_int4;
typedef _f_int8 _f_int2;
typedef _f_int8 _f_int1;
typedef long _f_log8;
typedef _f_log8 _f_log;
typedef _f_log8 _f_log4;
typedef _f_log8 _f_log2;
typedef _f_log8 _f_log1;
typedef double _f_real8;
typedef long double _f_real16;
typedef _f_real8 _f_real;
typedef _f_real16 _f_dble;
typedef _f_real8 _f_real4;
typedef _Complex double _f_comp8;

typedef _f_comp8 _f_comp;



typedef union _FCD {
    char *c_pointer;
    struct {
    unsigned bit_offset : 6,
        fcd_len : 26,
        word_addr : 32;
    } _F;
} _dcf;
typedef void *_fcd;




typedef void *_GPTR;
typedef const void *_GPTR2CONST;

extern _fcd _cptofcd (char *_Ccp, unsigned _Len);
extern char * _fcdtocp (_fcd _Fcd);
extern unsigned int _fcdlen (_fcd _Fcd);
extern _f_log _btol (long _BV);
extern long _lvtob (_f_log _LV);
extern long _ltob (_f_log *_LP);
extern char * _f2ccpy (_fcd f, ...);
extern char * _fc_copy (_fcd f, char *s, int slen);
extern char * _fc_acopy (_fcd f);
extern int _c2fcpy(char *c, _fcd f);
extern int _isfcd (long _P);

static _fcd
__cptofcd(char *c, unsigned int l);
#pragma _CRI inline __cptofcd
static _fcd
__cptofcd(char *c, unsigned int l)
{
    _dcf f;

    f.c_pointer = c;
    f._F.fcd_len = l << 3;

    return ((*(_fcd *) &f));
}
static char *
__fcdtocp(_fcd f);
#pragma _CRI inline __fcdtocp
static char *
__fcdtocp(_fcd f)
{
    char *c;
    _dcf d;
    d = (*(_dcf *) &f);
    d._F.fcd_len = 0;
    c = d.c_pointer;
    return (c);
}
#else
#define _fcd char *
#define _fcdtocp(a) a
#define _fcdlen(a) strlen(a)
#endif

#ifdef VAX
typedef struct {short length; short magic; char * address;}DESC;
#define _fcdtocp(a) fcdtocp(a)
#define _fcd DESC *
#endif

char *fcd2char(); /* fortran to c string convertion (alloc memory) */

#endif /* end of FORT2C_H */
