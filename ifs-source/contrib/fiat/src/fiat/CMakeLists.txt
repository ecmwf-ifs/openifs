# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

### Assemble list of definitions

if( TARGET fckit )
  list( APPEND FIAT_DEFINITIONS WITH_FCKIT )
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  list( APPEND FIAT_DEFINITIONS LINUX )
endif()

if( HAVE_DR_HOOK_MULTI_PRECISION_HANDLES )
   list( APPEND FIAT_DEFINITIONS DR_HOOK_MULTI_PRECISION_HANDLES )
endif()

if( FIAT_ATTRIBUTE_CONSTRUCTOR_SUPPORTED )
  list( APPEND FIAT_DEFINITIONS FIAT_ATTRIBUTE_CONSTRUCTOR_SUPPORTED )
endif()

### Compilation flags

if( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  ## To disable checking of argument correctness of dummy mpi symbols
  ecbuild_add_fortran_flags( -nowarn nointerfaces )
endif()

if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  ecbuild_add_fortran_flags( -ffree-line-length-none )
  if( CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL 10 )
      ecbuild_add_fortran_flags( -fallow-argument-mismatch )
  endif()
endif()

if( NOT fiat_VERSION_PATCH )
    set( fiat_VERSION_PATCH 0 )
endif()
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/library/version.c.in ${CMAKE_CURRENT_BINARY_DIR}/version.c @ONLY )

ecbuild_list_add_pattern( LIST fiat_src GLOB *.c *.F* *.cc )
set( fiat_src ${fiat_src} PARENT_SCOPE )

ecbuild_add_library( TARGET fiat
    LINKER_LANGUAGE Fortran
    SOURCES ${fiat_src} ${CMAKE_CURRENT_BINARY_DIR}/version.c
    PRIVATE_LIBS 
        ${CMAKE_DL_LIBS}
        ${RT_LIB}
    PRIVATE_INCLUDES
        ${CMAKE_CURRENT_SOURCE_DIR}/drhook/internal
        ${CMAKE_CURRENT_SOURCE_DIR}/ecsort/internal
        ${CMAKE_CURRENT_SOURCE_DIR}/util/internal
        ${CMAKE_CURRENT_SOURCE_DIR}/system/internal
    PUBLIC_INCLUDES
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/fiat>
        $<INSTALL_INTERFACE:include>
        $<INSTALL_INTERFACE:include/fiat>
)


if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  # Following should not be necessary; 
  # Probably a bug in the M1 prerelease of gfortran 10.2.0.4
  target_compile_definitions( fiat PRIVATE $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:__APPLE__>)
endif()

target_compile_definitions( fiat PRIVATE ${FIAT_DEFINITIONS} )

if( TARGET fckit )
    target_link_libraries( fiat PRIVATE fckit )
endif()

if( HAVE_MPI )
    target_link_libraries( fiat PUBLIC ${MPI_Fortran_LIBRARIES} )
    target_include_directories( fiat PRIVATE ${MPI_Fortran_INCLUDE_DIRS} )
else()
    target_link_libraries( fiat PRIVATE  mpi_serial )
endif()

if( HAVE_OMP )
    target_link_libraries( fiat PRIVATE OpenMP::OpenMP_Fortran )
endif()

fiat_target_ignore_missing_symbols( TARGET fiat SYMBOLS
  _MallocExtension_GetNumericProperty
  _MPI_Get_version
  _MPI_Get_library_version
)

## Install and Export ##

execute_process( COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/include )
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/include/fiat ${PROJECT_BINARY_DIR}/include/fiat )
install(
  DIRECTORY   ${CMAKE_CURRENT_SOURCE_DIR}/include/fiat
  DESTINATION include
)
install(
  FILES       ${CMAKE_CURRENT_SOURCE_DIR}/version.h
  DESTINATION include/fiat
)

fiat_target_fortran_module_directory(
  TARGET fiat
  MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/fiat
  INSTALL_DIRECTORY module/fiat
)

set( FIAT_DEFINITIONS ${FIAT_DEFINITIONS} PARENT_SCOPE )
