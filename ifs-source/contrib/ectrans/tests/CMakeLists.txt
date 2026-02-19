# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

### trans_test_install

if( HAVE_TESTS )

  find_package( MPI )
  set( LAUNCH )
  if( MPI_FOUND )
    set( LAUNCH "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1" )
  endif()
  configure_file( test-install.sh.in ${CMAKE_CURRENT_BINARY_DIR}/test-install.sh @ONLY )

  unset( _test_args )
  if( CMAKE_TOOLCHAIN_FILE )
    list( APPEND _test_args "-DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}" )
  endif()
  foreach( lang C CXX Fortran )
    if( CMAKE_${lang}_COMPILER )
      list( APPEND _test_args "-DCMAKE_${lang}_COMPILER=${CMAKE_${lang}_COMPILER}" )
    endif()
  endforeach()
  if( NOT HAVE_DOUBLE_PRECISION )
    list( APPEND _test_args "-DCOMPONENTS=single" )
  endif() 

  add_test( NAME ectrans_test_install
            COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test-install.sh ${_test_args} )
endif()


if( HAVE_DOUBLE_PRECISION )
  set( trans trans_dp )
else()
  set( trans trans_sp )
endif()

ecbuild_add_test(TARGET ectrans_test_adjoint
  SOURCES trans/test_adjoint.F90
  LIBS ${trans}
  LINKER_LANGUAGE Fortran
)
if( TEST ectrans_test_adjoint AND HAVE_OMP )
  target_link_libraries( ectrans_test_adjoint OpenMP::OpenMP_Fortran )
endif()


if( HAVE_TRANSI )

  check_include_files( malloc.h       EC_HAVE_MALLOC_H      )
  ecbuild_debug_var( EC_HAVE_MALLOC_H )


  if( EC_HAVE_MALLOC_H AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
    list( APPEND  TEST_DEFINITIONS
      TRANSI_HAVE_MEMORY
    )
  else()
    ecbuild_warn( "ectrans tests checking memory leaks are disabled as malloc.h was not found" )
  endif()

  ecbuild_add_library( TARGET ectrans_test
    SOURCES     transi/transi_test.h transi/transi_test.c
    PUBLIC_LIBS transi_dp
    NOINSTALL
  )
  target_compile_definitions( ectrans_test PUBLIC ${TEST_DEFINITIONS} )

  ecbuild_add_test( TARGET ectrans_test_transi_program
    SOURCES   transi/transi_test_program.c
    LIBS      ectrans_test
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_timings
    SOURCES   transi/transi_test_timings.c
    LIBS      ectrans_test
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_lonlat
    SOURCES   transi/transi_test_lonlat.c
    LIBS      ectrans_test
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_io
    SOURCES   transi/transi_test_io.c
    LIBS      ectrans_test
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_memory
    SOURCES   transi/transi_test_memory.c
    LIBS      ectrans_test
    CONDITION EC_HAVE_MALLOC_H
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_memory_lonlat
    SOURCES   transi/transi_test_memory_lonlat.c
    LIBS      ectrans_test
    CONDITION EC_HAVE_MALLOC_H
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_vordiv_to_UV
    SOURCES   transi/transi_test_vordiv_to_UV.c
    LIBS      ectrans_test
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_invtrans_adjoint
    SOURCES   transi/transi_test_invtrans_adjoint.c
    LIBS      ectrans_test
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_lonlat_diff_incr
    SOURCES   transi/transi_test_lonlat_diff_incr.c
    LIBS      ectrans_test
    ENVIRONMENT TRANS_USE_MPI=0 )

endif()
