/*
 * (C) Copyright 2021- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#ifndef FIAT_PREPROCESSOR_H
#define FIAT_PREPROCESSOR_H


#if __GNUC__ >= 7
#define FIAT_PP_FALLTHROUGH __attribute__((fallthrough))
#else
#define FIAT_PP_FALLTHROUGH
#endif


#define FIAT_PP_CAT_(v1, v2) v1 ## v2
#define FIAT_PP_CAT(v1, v2) FIAT_PP_CAT_(v1, v2)

#define FIAT_PP_CAT5_(_0, _1, _2, _3, _4) _0 ## _1 ## _2 ## _3 ## _4

#define FIAT_PP_IDENTITY_(x) x
#define FIAT_PP_IDENTITY(x) FIAT_PP_IDENTITY_(x)

#define FIAT_PP_VA_ARGS_(...) __VA_ARGS__
#define FIAT_PP_VA_ARGS(...) FIAT_PP_VA_ARGS_(__VA_ARGS__)

#define FIAT_PP_IDENTITY_VA_ARGS_(x, ...) x, __VA_ARGS__
#define FIAT_PP_IDENTITY_VA_ARGS(x, ...) FIAT_PP_IDENTITY_VA_ARGS_(x, __VA_ARGS__)

#define FIAT_PP_IIF_0(x, ...) __VA_ARGS__
#define FIAT_PP_IIF_1(x, ...) x
#define FIAT_PP_IIF(c) FIAT_PP_CAT_(FIAT_PP_IIF_, c)

#define FIAT_PP_HAS_COMMA(...) FIAT_PP_IDENTITY(FIAT_PP_VA_ARGS_TAIL(__VA_ARGS__, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0))
#define FIAT_PP_IS_EMPTY_TRIGGER_PARENTHESIS_(...) ,
 
#define FIAT_PP_IS_EMPTY(...) FIAT_PP_IS_EMPTY_( \
    /* test if there is just one argument, eventually an empty one */ \
    FIAT_PP_HAS_COMMA(__VA_ARGS__),                                \
    /* test if _TRIGGER_PARENTHESIS_ together with the argument adds a comma */ \
    FIAT_PP_HAS_COMMA(FIAT_PP_IS_EMPTY_TRIGGER_PARENTHESIS_ __VA_ARGS__), \
    /* test if the argument together with a parenthesis adds a comma */ \
    FIAT_PP_HAS_COMMA(__VA_ARGS__ ()),                             \
    /* test if placing it between _TRIGGER_PARENTHESIS_ and the parenthesis adds a comma */ \
    FIAT_PP_HAS_COMMA(FIAT_PP_IS_EMPTY_TRIGGER_PARENTHESIS_ __VA_ARGS__ ()))
 
#define FIAT_PP_IS_EMPTY_(_0, _1, _2, _3) FIAT_PP_HAS_COMMA(FIAT_PP_CAT5_(FIAT_PP_IS_EMPTY_IS_EMPTY_CASE_, _0, _1, _2, _3))
#define FIAT_PP_IS_EMPTY_IS_EMPTY_CASE_0001 ,

#define FIAT_PP_VARIADIC_SIZE(...) FIAT_PP_IIF(FIAT_PP_IS_EMPTY(__VA_ARGS__))(0, FIAT_PP_VARIADIC_SIZE_(__VA_ARGS__, FIAT_PP_VA_ARGS_SEQ16()))
#define FIAT_PP_VARIADIC_SIZE_(...) FIAT_PP_IDENTITY(FIAT_PP_VA_ARGS_TAIL(__VA_ARGS__))

#define FIAT_PP_VA_ARGS_TAIL(_0,_1,_2,_3,_4,_5,_6,_7,_8,_9,_10,_11,_12,_13,_14, x, ...) x
#define FIAT_PP_VA_ARGS_SEQ16() 15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0

#endif
