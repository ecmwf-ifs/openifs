# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

ecbuild_info( "build type                   : [${CMAKE_BUILD_TYPE}]" )
set( Fortran_flags_str "Fortran flags" )
set( C_flags_str       "C flags      " )
set( CXX_flags_str     "C++ flags    " )
string( TOUPPER ${PROJECT_NAME} PNAME )
    foreach( lang Fortran C CXX )
        set( flags "${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${CMAKE_BUILD_TYPE_CAPS}} ${${PNAME}_${lang}_FLAGS} ${${PNAME}_${lang}_FLAGS_${CMAKE_BUILD_TYPE_CAPS}}" )
        string(REGEX REPLACE "[ ]+" " " flags ${flags})
        string(STRIP "${flags}" flags)
ecbuild_info( "${${lang}_flags_str}                : [${flags}]" )
    endforeach()
ecbuild_info( "OMP" )
    foreach( lang Fortran )
ecbuild_info( "    OpenMP_${lang}_FLAGS     : [${OpenMP_${lang}_FLAGS}]" )
    endforeach()
ecbuild_info( "BLAS/LAPACK" )
    if( HAVE_SINGLE_PRECISION AND HAVE_DOUBLE_PRECISION AND ECTRANS_CRAYHACK_DOUBLE_PRECISION_WITHOUT_MKL )
ecbuild_info( "    trans_dp                 : [${LAPACK_dp}]" )
ecbuild_info( "    trans_sp                 : [${LAPACK_sp}]" )
    else()
ecbuild_info( "    LAPACK_LIBRARIES         : [${LAPACK_LIBRARIES}]" )
    endif()
ecbuild_info( "FFTW" )
ecbuild_info( "    FFTW_LIBRARIES           : [${FFTW_LIBRARIES}]" )
ecbuild_info( "---------------------------------------------------------" )




