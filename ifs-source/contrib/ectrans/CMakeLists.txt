# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )
find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )

project( ectrans LANGUAGES C Fortran )
include( ectrans_macros )

ecbuild_enable_fortran( REQUIRED NO_MODULE_DIRECTORY )

### Find (optional) dependencies 

ecbuild_find_package( NAME fiat REQUIRED )

ecbuild_add_option( FEATURE MPI 
                    DESCRIPTION "Support for MPI distributed memory parallelism"
                    CONDITION fiat_HAVE_MPI )

ecbuild_add_option( FEATURE OMP
                    DEFAULT ON
                    DESCRIPTION "Support for OpenMP shared memory parallelism"
                    REQUIRED_PACKAGES "OpenMP COMPONENTS Fortran" )

ecbuild_add_option( FEATURE DOUBLE_PRECISION
                    DEFAULT ON
                    DESCRIPTION "Support for Double Precision" )

ecbuild_add_option( FEATURE SINGLE_PRECISION
                    DEFAULT ON
                    DESCRIPTION "Support for Single Precision" )

if( HAVE_SINGLE_PRECISION )
  set( single "single" )
endif()

ecbuild_add_option( FEATURE MKL
                    DESCRIPTION "Use MKL for BLAS and/or FFTW"
                    DEFAULT ON
                    REQUIRED_PACKAGES "MKL QUIET" )

if( NOT HAVE_MKL )
    option( FFTW_ENABLE_MKL OFF )
endif()

ecbuild_add_option( FEATURE FFTW
                    DEFAULT ON
                    DESCRIPTION "Support for fftw"
                    REQUIRED_PACKAGES "FFTW COMPONENTS double ${single}" )

ecbuild_add_option( FEATURE TRANSI
                    DEFAULT ON
                    DESCRIPTION "Compile TransI C-interface to trans"
                    CONDITION HAVE_DOUBLE_PRECISION )

ectrans_find_lapack()

### Add sources and tests
include( ectrans_compile_options )
add_subdirectory( src )
add_subdirectory( tests )

### Export
if( BUILD_SHARED_LIBS )
  set( PACKAGE_REQUIRES_PRIVATE_DEPENDENCIES 0 )
else()
  set( PACKAGE_REQUIRES_PRIVATE_DEPENDENCIES 1 )
endif()

ecbuild_install_project( NAME ${PROJECT_NAME} )

ecbuild_print_summary()

