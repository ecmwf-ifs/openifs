# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )
project( ectrans_test_install VERSION 0.0.0 LANGUAGES Fortran )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if( COMPONENTS )
  find_package( ectrans REQUIRED COMPONENTS ${COMPONENTS} )
else()
  find_package( ectrans REQUIRED )
endif()

if( TARGET trans_dp )
  add_executable( main_dp main.F90 )
  target_link_libraries( main_dp trans_dp )
endif()

if( TARGET trans_sp )
  add_executable( main_sp main.F90 )
  target_link_libraries( main_sp trans_sp )
endif()

if( TARGET transi_dp )
  enable_language( C )
  add_executable( transi_sptogp transi_sptogp.c )
  target_link_libraries( transi_sptogp transi_dp )

  add_executable( transi_gptosp transi_gptosp.c )
  target_link_libraries( transi_gptosp transi_dp )

endif()
