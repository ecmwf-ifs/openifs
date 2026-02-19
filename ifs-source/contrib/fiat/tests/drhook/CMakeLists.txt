find_package( CMath )

ecbuild_add_executable( TARGET drhook_ex1
  SOURCES drhook_ex1.F90 
  LIBS fiat
  NOINSTALL )

ecbuild_add_executable( TARGET drhook_ex2
  SOURCES drhook_ex2.F90 
  LIBS fiat
  NOINSTALL )

ecbuild_add_executable( TARGET drhook_ex3
  SOURCES drhook_ex3.F90 
  LIBS fiat
  NOINSTALL )

ecbuild_add_executable( TARGET drhook_ex4
  SOURCES drhook_ex4.c
  LIBS fiat ${CMATH_LIBRARIES}
  NOINSTALL )

ecbuild_add_executable( TARGET drhook_ex5
  SOURCES drhook_ex5.F90
  LIBS fiat
  NOINSTALL )


set( FPE_TRAPPING_SUPPORTED 1 )
if( APPLE )
  set( FPE_TRAPPING_SUPPORTED 0 )
endif()

# Test floating point trapping
ecbuild_add_test( TARGET fiat_test_drhook_ex1
                  COMMAND drhook_ex1
                  ENVIRONMENT DR_HOOK=1 DR_HOOK_TRAPFPE=1 DR_HOOK_ASSERT_MPI_INITIALIZED=0 DR_HOOK_OPT=NOPROPAGATE_SIGNALS )

if( FPE_TRAPPING_SUPPORTED )
  set_tests_properties( fiat_test_drhook_ex1 PROPERTIES WILL_FAIL TRUE )
endif()

# Test watchpoint point trapping
ecbuild_add_test( TARGET fiat_test_drhook_ex2
                  COMMAND drhook_ex2
                  ENVIRONMENT DR_HOOK=1 DR_HOOK_ASSERT_MPI_INITIALIZED=0 DR_HOOK_OPT=NOPROPAGATE_SIGNALS )
set_tests_properties( fiat_test_drhook_ex2 PROPERTIES WILL_FAIL TRUE )

# Play with profiling options
ecbuild_add_test( TARGET fiat_test_drhook_ex3
                  COMMAND drhook_ex3
                  ENVIRONMENT DR_HOOK=1 DR_HOOK_ASSERT_MPI_INITIALIZED=0 )

# Test floating point trapping from C program
ecbuild_add_test( TARGET fiat_test_drhook_ex4
                  COMMAND drhook_ex4
                  ENVIRONMENT DR_HOOK=1 SIGFPE=1 DR_HOOK_ASSERT_MPI_INITIALIZED=0 DR_HOOK_OPT=NOPROPAGATE_SIGNALS )
if( FPE_TRAPPING_SUPPORTED )
  set_tests_properties( fiat_test_drhook_ex4 PROPERTIES WILL_FAIL TRUE )
endif()

# Play MPI and dr_hook together
ecbuild_add_test( TARGET fiat_test_drhook_ex5
                  COMMAND drhook_ex5
                  MPI 2
                  ENVIRONMENT DR_HOOK=1 DR_HOOK_SILENT=1 )

# TODO:
# Better parse output to see if it matches.
