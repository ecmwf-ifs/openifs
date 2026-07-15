# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )
find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )

project( fiat LANGUAGES Fortran C CXX )
include( fiat_macros )
ecbuild_enable_fortran( REQUIRED NO_MODULE_DIRECTORY )

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED TRUE)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

### Options

ecbuild_add_option( FEATURE OMP
                    DESCRIPTION "support for OpenMP shared memory parallelism"
                    REQUIRED_PACKAGES "OpenMP COMPONENTS Fortran" )

ecbuild_add_option( FEATURE MPI
                    DESCRIPTION "Support for MPI distributed parallelism"
                    REQUIRED_PACKAGES "MPI COMPONENTS Fortran" )

ecbuild_add_option( FEATURE FCKIT
                    DESCRIPTION "Support for fckit"
                    REQUIRED_PACKAGES "fckit QUIET" )

ecbuild_add_option( FEATURE DR_HOOK_MULTI_PRECISION_HANDLES
                    DESCRIPTION "[DEPRECATED] Support single precision handles for DR_HOOK"
                    DEFAULT OFF )

ecbuild_add_option( FEATURE DOUBLE_PRECISION
                    DEFAULT ON
                    DESCRIPTION "Support for Double Precision" )

ecbuild_add_option( FEATURE SINGLE_PRECISION
                    DEFAULT ON
                    DESCRIPTION "Support for Single Precision" )

ecbuild_add_option( FEATURE WARNINGS
                    DEFAULT ON
                    DESCRIPTION "Add warnings to compiler" )

ecbuild_find_package( NAME Realtime QUIET )

### Sources

include( fiat_compiler_warnings )
add_subdirectory(src)
add_subdirectory(share)

### Tests

add_subdirectory(tests)

### Export

ecbuild_install_project( NAME fiat )

ecbuild_print_summary()

