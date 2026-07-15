# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set( HAVE_dp ${HAVE_DOUBLE_PRECISION} )
set( HAVE_sp ${HAVE_SINGLE_PRECISION} )

foreach( prec sp dp )
#  if( HAVE_${prec} )
 
  ecbuild_add_library(
    TARGET  parkind_${prec}
    SOURCES parkind1.F90 
            parkind2.F90
  )

  fiat_target_fortran_module_directory(
    TARGET            parkind_${prec}
    MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/parkind_${prec}
    INSTALL_DIRECTORY module/parkind_${prec}
  )

  if( ${prec} MATCHES "sp" )

    if( EC_OS_BITS EQUAL "64" )
      target_compile_definitions( parkind_${prec} PRIVATE ADDRESS64 )
    endif()

    target_compile_definitions( parkind_${prec} PRIVATE PARKIND1_SINGLE )

  endif()
  
#  endif()
endforeach()


