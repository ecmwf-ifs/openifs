# (C) Copyright 2022- ECMWF.
# (C) Copyright 2022- Meteo-France.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

cmake_minimum_required (VERSION 3.24)

## find ecbuild
find_package( ecbuild 3.7 QUIET HINTS ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )
if( NOT ecbuild_FOUND )
  include(cmake/field_api_fetchcontent_ecbuild.cmake)
endif()

## find ecbuild toolchain
if( DEFINED ENV{ECBUILD_TOOLCHAIN} AND NOT CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE ${CMAKE_CURRENT_SOURCE_DIR}/$ENV{ECBUILD_TOOLCHAIN})
endif()

# Inherit precision setting from parent project
if( HAVE_SINGLE_PRECISION )
  set( PROJECT_NAME field_api_sp )
else()
  set( PROJECT_NAME field_api_dp )
endif()

## define project and field_api target name
project(${PROJECT_NAME} VERSION 0.2.1 LANGUAGES Fortran)
set(LIBNAME field_api)

## find required packages
ecbuild_find_package(OpenMP COMPONENTS Fortran REQUIRED)

## find fypp
find_program(FYPP fypp HINTS ${fypp_ROOT})
if( NOT FYPP)
  include(cmake/field_api_fetchcontent_fypp.cmake)
  set(FYPP ${fypp_SOURCE_DIR}/bin/fypp)
  set(fypp_ROOT ${fypp_SOURCE_DIR}/bin PARENT_SCOPE)
  ecbuild_info("fypp downloaded to: ${FYPP}")
endif()

## find OpenACC
### Workaround for CMake bug (fixed in CMake 3.25)
find_package(OpenACC COMPONENTS Fortran)
if( OpenACC_Fortran_FOUND )
   set( OpenACC_FOUND ON )
endif()
ecbuild_add_option( FEATURE ACC
    DESCRIPTION "OpenACC" DEFAULT ON
    CONDITION OpenACC_FOUND )

## set general compiler flags
include(cmake/field_api_compile_options.cmake)

## choose to build dp or sp or both
ecbuild_add_option( FEATURE SINGLE_PRECISION
	DESCRIPTION "Compile field_api in single precision"
    DEFAULT OFF)
ecbuild_add_option( FEATURE DOUBLE_PRECISION
	DESCRIPTION "Compile field_api in double precision"
    DEFAULT ON
    CONDITION NOT HAVE_SINGLE_PRECISION)

if(HAVE_SINGLE_PRECISION)
    set(prec sp)
elseif(HAVE_DOUBLE_PRECISION)
    set(prec dp)
endif()

## find fiat
ecbuild_find_package(NAME fiat)
if( NOT fiat_FOUND )
  if(NOT UTIL_MODULE_PATH)
    ecbuild_critical("If not building with fiat, then the path for utility modules must be specified")
  endif()

  ecbuild_info( "Checking for FIAT components in ${UTIL_MODULE_PATH}" )

  find_file( ABOR1_PATH abor1.F90 REQUIRED
      HINTS ${UTIL_MODULE_PATH} ${UTIL_MODULE_PATH}/src/fiat/util
  )
  ecbuild_info( "Found ABOR1: ${ABOR1_PATH}" )

  find_file( OML_PATH oml_mod.F90 REQUIRED
      HINTS ${UTIL_MODULE_PATH} ${UTIL_MODULE_PATH}/src/fiat/oml
  )
  ecbuild_info( "Found OML: ${OML_PATH}" )

  find_file( PARKIND1_PATH parkind1.F90 REQUIRED
      HINTS ${UTIL_MODULE_PATH} ${UTIL_MODULE_PATH}/src/parkind
  )
  ecbuild_info( "Found PARKIND1: ${PARKIND1_PATH}" )
  list(APPEND srcs ${ABOR1_PATH} ${OML_PATH} ${PARKIND1_PATH})
endif()

list(APPEND srcs field_basic_module.F90 field_init_debug_value_module.F90 dev_alloc.c)

## preprocess fypp files
foreach (SUFF IN ITEMS IM RM RB RD LM)
  string (TOLOWER ${SUFF} suff)
  foreach (RANK RANGE 2 5)
    foreach (FUNC IN ITEMS "" _gathscat _access _util _array_util)
      add_custom_command (OUTPUT field_${RANK}${suff}${FUNC}_module.F90
	      COMMAND ${FYPP} -n -DRANK=${RANK} -DSUFF='${SUFF}' -m os -M ${CMAKE_CURRENT_SOURCE_DIR} -m fieldType 
        ${CMAKE_CURRENT_SOURCE_DIR}/field_RANKSUFF${FUNC}_module.fypp > field_${RANK}${suff}${FUNC}_module.F90
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/field_RANKSUFF${FUNC}_module.fypp
        VERBATIM)
      list(APPEND srcs "field_${RANK}${suff}${FUNC}_module.F90")
    endforeach ()
  endforeach ()
endforeach ()

foreach (SRC IN ITEMS dev_alloc_module field_factory_module field_access_module 
         field_array_module field_module field_gathscat_module field_util_module field_array_util_module)
  add_custom_command (OUTPUT ${SRC}.F90
	  COMMAND ${FYPP} -n -m os -M ${CMAKE_CURRENT_SOURCE_DIR} -m fieldType ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.fypp > ${SRC}.F90
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.fypp
    VERBATIM)
  list(APPEND srcs "${SRC}.F90")
endforeach ()

## add field_api target
ecbuild_add_library(
    TYPE STATIC
    TARGET ${LIBNAME}_${prec}
    SOURCES ${srcs}
    DEFINITIONS $<$<NOT:${fiat_FOUND}>:${FIELD_API_DEFINITIONS}>
    PUBLIC_LIBS
       $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
    PRIVATE_LIBS
       $<${fiat_FOUND}:fiat>
       $<${fiat_FOUND}:parkind_${prec}>
       OpenMP::OpenMP_Fortran
    )
target_include_directories( ${LIBNAME}_${prec} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET ${LIBNAME}_${prec} PROPERTY C_STANDARD 99)
set_target_properties( ${LIBNAME}_${prec} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include/${LIBNAME}_${prec} )

# export target usage interface
target_include_directories( ${LIBNAME}_${prec}
                            INTERFACE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include/${LIBNAME}_${prec}> 
                            INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${LIBNAME}_${prec}> 
                          )

# set install location for .mod files
install(DIRECTORY ${CMAKE_BINARY_DIR}/include/${LIBNAME}_${prec} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

## build tests
ecbuild_add_option(
    FEATURE TESTS
    DEFAULT ON
    DESCRIPTION "Enable field_api testing"
    CONDITION fiat_FOUND
)
if( HAVE_TESTS )
  include(CTest)
  add_subdirectory(tests)
endif()

## export the project
ecbuild_install_project(NAME ${LIBNAME}_${prec})
