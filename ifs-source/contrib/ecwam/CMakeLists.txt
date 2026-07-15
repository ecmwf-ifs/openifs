# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

cmake_minimum_required( VERSION 3.24 FATAL_ERROR )
find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )

set(PARENT_PROJECT_NAME "${PROJECT_NAME}")
if( NOT ECWAM_PROJECT_NAME )
  set( ECWAM_PROJECT_NAME ecwam )
endif()

project( ${ECWAM_PROJECT_NAME} LANGUAGES Fortran C CXX )

include( ecwam_macros )
ecbuild_enable_fortran( REQUIRED NO_MODULE_DIRECTORY )

### Dependencies

ecbuild_find_package( fiat         REQUIRED )
ecbuild_find_package( eccodes      REQUIRED )

### Options

ecbuild_add_option( FEATURE OMP
                    DESCRIPTION "support for OpenMP shared memory parallelism"
                    REQUIRED_PACKAGES "OpenMP COMPONENTS Fortran" )

ecbuild_add_option( FEATURE MPI
                    DESCRIPTION "Support for MPI distributed parallelism"
                    REQUIRED_PACKAGES "MPI COMPONENTS Fortran" )

ecbuild_add_option( FEATURE MULTIO
                    DESCRIPTION "Support for IO to the FDB"
                    REQUIRED_PACKAGES "multio" )

ecbuild_add_option( FEATURE SINGLE_PRECISION
                    DEFAULT OFF
                    DESCRIPTION "Support for Single Precision" )

ecbuild_add_option( FEATURE UNWAM
                    DEFAULT OFF
                    DESCRIPTION "Support for UNWAM" )

ecbuild_add_option( FEATURE OCEAN_COUPLING 
                    DEFAULT ON
                    DESCRIPTION "Support ocean coupling"
                    CONDITION ${PNAME}_OCEANMODEL_LIBRARIES )

ecbuild_add_option( FEATURE ECWAM_LAUNCH
                    DEFAULT ON
                    DESCRIPTION "Use share/ecwam/scripts/ecwam-launch script to run ecwam tests" )

ecbuild_add_option( FEATURE ECFLOW
                    DESCRIPTION "ecflow meter updates"
                    REQUIRED_PACKAGES "ecflow_light")

if( HAVE_ECWAM_LAUNCH )
  # Overwrite meaning of mpiexec
  set( ECWAM_LAUNCH ${CMAKE_CURRENT_SOURCE_DIR}/share/ecwam/scripts/ecwam-launch )
  set( MPIEXEC_EXECUTABLE ${ECWAM_LAUNCH} )
  set( MPIEXEC_NUMPROC_FLAG   -n )
  set( MPIEXEC_NUMTHREAD_FLAG -c )
  set( CMAKE_CROSSCOMPILING_EMULATOR ${ECWAM_LAUNCH} ) # For serial runs
endif()

### Optional dependencies

ecwam_find_python_mods()
ecbuild_add_option( FEATURE GEN_DERIV_TYPES
                    DEFAULT ON
                    DESCRIPTION "Generate derived-type data structures using fypp pre-processor"
                    CONDITION FYPP_FOUND AND PYYAML_FOUND )

if( HAVE_SINGLE_PRECISION )
  set( prec sp )
else()
  set( prec dp )
endif()

ecbuild_find_package( field_api_${prec} REQUIRED )

### Set third party libraries

unset( OpenMP_Fortran_LIBRARIES )
if( HAVE_OMP )
  set( OpenMP_Fortran_LIBRARIES OpenMP::OpenMP_Fortran )
endif()

### Sources

include(ecwam_compile_flags)
add_subdirectory(src)
add_subdirectory(share)

### Tests

add_subdirectory(tests)

get_property(ecwam_test_configs GLOBAL PROPERTY ecwam_test_configs)
ecwam_cache_bathymetry( ${ecwam_test_configs} )

### Documentation

add_subdirectory(doc)

### Export

ecbuild_install_project( NAME ${PROJECT_NAME} )

ecbuild_print_summary()

