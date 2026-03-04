# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if( HAVE_TESTS )

if( HAVE_MPI AND MPIEXEC )
  set( LAUNCH ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 )
else()
  set( LAUNCH ${CMAKE_CROSSCOMPILING_EMULATOR} )
endif()


if( HAVE_OMP )
  set( OMP OpenMP::OpenMP_Fortran )
endif()

if( HAVE_FCKIT )
  list( APPEND FCKIT_DEFINITIONS WITH_FCKIT )
  list( APPEND FCKIT_LIB fckit )
endif()

if( HAVE_FCKIT AND fckit_VERSION VERSION_GREATER_EQUAL 0.9.1 )
ecbuild_add_executable(
  TARGET  fiat-test-abort-exception-handler
  SOURCES test_abort_exception_handler.F90 test_abort_exception_handler.cc 
  LIBS    fiat fckit
  LINKER_LANGUAGE Fortran
  NOINSTALL)
add_test(NAME fiat_test_abort_exception_handler
  COMMAND ${CMAKE_COMMAND}
  "-DEXECUTABLE=$<TARGET_FILE:fiat-test-abort-exception-handler>"
  "-DLAUNCH=${LAUNCH}"
  "-DPASS_REGULAR_EXPRESSION=problem in function_2"
  -P ${CMAKE_CURRENT_SOURCE_DIR}/test_program_output.cmake )
endif()


add_subdirectory( drhook )




# ----------------------------------------------------------------------------------------
# Test: fiat_test_drhook_no_output
#
# Ensures that no output is generated with DR_HOOK_SILENT=1

ecbuild_add_executable(
    TARGET  fiat-test-drhook-no-output
    SOURCES test_drhook_no_output.c 
    LIBS    fiat
    NOINSTALL )
add_test( NAME fiat_test_drhook_no_output
          COMMAND ${CMAKE_COMMAND}
          "-DLAUNCH=${LAUNCH}"
          "-DEXECUTABLE=$<TARGET_FILE:fiat-test-drhook-no-output>"
          -P ${CMAKE_CURRENT_SOURCE_DIR}/test_drhook_no_output.cmake )
set_tests_properties(fiat_test_drhook_no_output PROPERTIES
    ENVIRONMENT "DR_HOOK_ASSERT_MPI_INITIALIZED=0"
)


# ----------------------------------------------------------------------------------------
# Tests: fiat_test_drhook_abort_my_abort, fiat_test_drhook_abort_abor1fl
#
# Test that we can use drhook within C and register a custom abort

ecbuild_add_executable(
    TARGET  fiat-test-drhook-abort
    SOURCES test_drhook_abort.c
    LIBS    fiat
    NOINSTALL )


add_test( NAME fiat_test_drhook_abort_abor1fl
          COMMAND fiat-test-drhook-abort )
set_tests_properties(fiat_test_drhook_abort_abor1fl PROPERTIES
    ENVIRONMENT "DR_HOOK_ASSERT_MPI_INITIALIZED=0"
    PASS_REGULAR_EXPRESSION "ABOR1.*EC_DRHOOK.*\[DrHookCallTree\]"
)

add_test( NAME fiat_test_drhook_abort_custom_abort
          COMMAND fiat-test-drhook-abort )
set_tests_properties(fiat_test_drhook_abort_custom_abort PROPERTIES
    ENVIRONMENT "CUSTOM_ABORT=1;DR_HOOK_ASSERT_MPI_INITIALIZED=0"
    PASS_REGULAR_EXPRESSION "custom_abort.*EC_DRHOOK.*\[DrHookCallTree\]"
)

# ----------------------------------------------------------------------------------------

ecbuild_add_test(
    TARGET  fiat_test_drhook_fortran
    SOURCES test_drhook_fortran.F90 
    LIBS    fiat ${OMP} ${FCKIT_LIB}
    LINKER_LANGUAGE Fortran
    DEFINITIONS ${FCKIT_DEFINITIONS})
set_tests_properties(fiat_test_drhook_fortran PROPERTIES
    ENVIRONMENT "MPL=0;DR_HOOK_ASSERT_MPI_INITIALIZED=0;DR_HOOK_OPT=NOPROPAGATE_SIGNALS"
    PASS_REGULAR_EXPRESSION "EC_DRHOOK.*\[DrHookCallTree\]" )

# ----------------------------------------------------------------------------------------
# Tests: fiat_test_ec_args_fortran

ecbuild_add_test(
    TARGET  fiat_test_ec_args_fortran 
    SOURCES test_ec_args_fortran.F90 
    ARGS arg1 arg2 arg3
    LIBS    fiat ${OMP}
    LINKER_LANGUAGE Fortran )

# ----------------------------------------------------------------------------------------
# Tests: fiat_test_bytes_io

ecbuild_add_test(
    TARGET  fiat_test_bytes_io 
    SOURCES test_bytes_io.F90 
    LIBS    fiat
    LINKER_LANGUAGE Fortran )

# ----------------------------------------------------------------------------------------
# Tests: fiat_test_mpl_no_output
#   ensures that no output is generated in this test upon mpl_init() and mpl_end()

ecbuild_add_executable(
    TARGET  fiat-test-mpl-no-output
    SOURCES test_mpl_no_output.F90 
    LIBS    fiat
    LINKER_LANGUAGE Fortran
    NOINSTALL )

add_test(NAME fiat_test_mpl_no_output
         COMMAND ${CMAKE_COMMAND}
         "-DEXECUTABLE=$<TARGET_FILE:fiat-test-mpl-no-output>"
         "-DLAUNCH=${LAUNCH}"
         -P ${CMAKE_CURRENT_SOURCE_DIR}/test_mpl_no_output.cmake )

# ----------------------------------------------------------------------------------------
# Tests: fiat_test_abor1
#
# Test abor1 function, and check the error message is printed

ecbuild_add_executable(
    TARGET  fiat-test-abor1
    SOURCES test_abor1.F90 
    LIBS    fiat ${OMP}
    LINKER_LANGUAGE Fortran
    NOINSTALL )
add_test(NAME fiat_test_abor1
      COMMAND ${CMAKE_COMMAND}
      "-DEXECUTABLE=$<TARGET_FILE:fiat-test-abor1>"
      "-DPASS_REGULAR_EXPRESSION=ABOR1.*aborting from OpenMP parallel region"
      -P ${CMAKE_CURRENT_SOURCE_DIR}/test_program_output.cmake )

# ----------------------------------------------------------------------------------------
# Test installation of fiat is working

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

add_test( NAME fiat_test_install
          COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test-install.sh ${_test_args} )

# ----------------------------------------------------------------------------------------

endif(HAVE_TESTS)
