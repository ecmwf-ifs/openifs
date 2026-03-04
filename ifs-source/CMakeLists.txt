# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

cmake_minimum_required( VERSION 3.15.0 FATAL_ERROR )

find_package( ecbuild 3.1 REQUIRED )

# NOTE: Cannot use ecbuild_add_option() here. These options are used in
# compile_flags.cmake which seems to be included in ecbuild_system.cmake.
option(IFS_CHECK_BOUNDS "Enable bounds checking" OFF)
option(IFS_INIT_SNAN "Initialize REALs with signaling NaNs" OFF)
option(IFS_SHARED_LIBS "Build shared libraries (overrides BUILD_SHARED_LIBS)" ON)

set(BUILD_SHARED_LIBS ${IFS_SHARED_LIBS})

set( PROJECT_NAME ifs )

ecbuild_parse_version_file( VERSION.txt PREFIX ifs )

# If the directory name matches ifs_dp or ifs_sp, that will be chosen as PROJECT_NAME instead
get_filename_component( ifs-source_dir ${CMAKE_CURRENT_SOURCE_DIR} NAME )
if( ifs-source_dir MATCHES "ifs_(dp|sp)" )
  set( PROJECT_NAME ${ifs-source_dir} )
  if( ifs-source_dir MATCHES "ifs_dp" )
    set( PROJECT_PRECISION_SUFFIX _dp )
  else()
    set( PROJECT_PRECISION_SUFFIX _sp )
  endif()
endif()

# The custom compile_flags have to be set before project() is called
string( TOUPPER ${PROJECT_NAME} PNAME )
if(NOT ${PNAME}_ECBUILD_COMPILE_FLAGS)
  set( ${PNAME}_ECBUILD_COMPILE_FLAGS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compile_flags.cmake CACHE FILEPATH "")
endif()

project( ${PROJECT_NAME} VERSION ${ifs_VERSION} LANGUAGES Fortran C CXX )

include( add_symlink )
include( ifs_propagate_flags )

# ecbuild_enable_fortran: Set Fortran module directory. With ECBUILD_2_COMPAT=ON this automatically gets added as a global
# include as well, but we don't want to export it as include_directory as it may otherwise
# lead to sp builds including dp modules. Therefore temporarily set ECBUILD_2_COMPAT=OFF until 
# all of ifs-source works with ECBUILD_2_COMPAT=OFF
set( ECBUILD_2_COMPAT.backup ${ECBUILD_2_COMPAT} )
set( ECBUILD_2_COMPAT OFF )
ecbuild_enable_fortran ( REQUIRED MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/module )
set( ECBUILD_2_COMPAT ${ECBUILD_2_COMPAT.backup} )

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
  include_directories(${CMAKE_BINARY_DIR}/include) # for grib_api.mod
endif()

set( IFS_CYCLE ${${PROJECT_NAME}_VERSION_MAJOR} )
if( ${PROJECT_NAME}_VERSION_MINOR GREATER 0 )
  set( IFS_CYCLE ${IFS_CYCLE}R${${PROJECT_NAME}_VERSION_MINOR} )
endif()

ecbuild_add_option(FEATURE MPI        REQUIRED_PACKAGES "MPI COMPONENTS Fortran" )
ecbuild_add_option(FEATURE HDF5       REQUIRED_PACKAGES "HDF5 COMPONENTS Fortran Fortran_HL")
ecbuild_add_option(FEATURE IFS_NETCDF REQUIRED_PACKAGES "NetCDF COMPONENTS Fortran C") # C component for few prepdata tools

ecbuild_add_option(FEATURE NAG        REQUIRED_PACKAGES "NAG")

ecbuild_add_option(FEATURE AEOLUS     DESCRIPTION "Build Aeolus project"   DEFAULT ON)
ecbuild_add_option(FEATURE NEMO       DESCRIPTION "Build NEMO project"     DEFAULT ON)
ecbuild_add_option(FEATURE OBSTAT     DESCRIPTION "Build Obstat project"   DEFAULT ON)
ecbuild_add_option(FEATURE OOPS       REQUIRED_PACKAGES "oops"             DEFAULT ON)

ecbuild_add_option(FEATURE SINGLE_PRECISION DESCRIPTION "Build IFS in single precision" DEFAULT OFF)
ecbuild_add_option(FEATURE FORECAST_ONLY DESCRIPTION "Build IFS in forecast-only configuration" DEFAULT OFF)
ecbuild_add_option(FEATURE OPENIFS_ONLY
                   DESCRIPTION "Build IFS in OpenIFS-only configuration"
                   DEFAULT OFF
                   CONDITION HAVE_FORECAST_ONLY)


ecbuild_use_package(PROJECT eccodes  VERSION 2.0   REQUIRED)
ecbuild_use_package(PROJECT fdb5     VERSION 5.0   REQUIRED)



ecbuild_use_package(PROJECT eckit    VERSION 0.17  REQUIRED)
ecbuild_use_package(PROJECT fckit    VERSION 0.5   REQUIRED)
ecbuild_use_package(PROJECT multio   VERSION 0.3   REQUIRED)

ecbuild_add_option(FEATURE ATLAS     REQUIRED_PACKAGES "PROJECT atlas VERSION 0.17.1" DEFAULT ON)

ecbuild_add_option(FEATURE MGRIDS DESCRIPTION "Build mgrids" DEFAULT ON)

# Disable MKL for double precision on Cray for cycles less than 48R1 to stay bit-identical with operations.
# This code must be removed in CY48R1!
if( ENABLE_MKL AND NOT HAVE_SINGLE_PRECISION )
  if( CMAKE_Fortran_COMPILER_ID MATCHES "Cray"
      AND BUILD_SITE MATCHES "cc[abt].*"
      AND ${PROJECT_VERSION} VERSION_LESS 48.1 )
    ecbuild_warn("Explicitly disabled MKL for ${IFS_CYCLE} in double precision to be bit-identical with operations of CY47R1.\n"
                 "To FORCE MKL yet (not bit-identical), set ${PNAME}_ENABLE_MKL=ON")
    set( ECTRANS_ENABLE_MKL ON )
    set( ECTRANS_CRAYHACK_DOUBLE_PRECISION_WITHOUT_MKL ON )
    if( IFS_DP_ENABLE_MKL )
      set( ECTRANS_CRAYHACK_DOUBLE_PRECISION_WITHOUT_MKL OFF )
    endif()
    set( ENABLE_MKL OFF )
  endif()
endif()

# Intel 18.* has problems compiling ifs/oops/fields_io_mod, which is only used by OOPS.
# Deactivate OOPS-IFS compilation altogether for this compiler version.
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  if( CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 19)
    set( HAVE_OOPS 0 )
    set( ifs_HAVE_OOPS 0 )
    set( ${PROJECT_NAME}_HAVE_OOPS 0 )
  endif()
endif()

ecbuild_find_package( NAME MKL QUIET )
ecbuild_add_option( FEATURE MKL
                    DESCRIPTION "Use MKL for BLAS dgemm"
                    DEFAULT OFF
                    CONDITION MKL_FOUND )

if( HAVE_MKL )
  set( LAPACK_LIBRARIES ${MKL_LIBRARIES} )
  set( HAVE_LAPACK ON )
else()
  if( NOT ( BUILD_SITE MATCHES "cc[abt].*" ) )
    ecbuild_add_option( FEATURE LAPACK
                        DESCRIPTION "Linear Algebra PACKage"
                        REQUIRED_PACKAGES "LAPACK" )
    if( HAVE_LAPACK )
      if( TARGET lapack )
        ecbuild_debug( "LAPACK found as CMake target lapack" )
        set( LAPACK_LIBRARIES lapack )
      endif()
    endif()
    ecbuild_debug_var( LAPACK_LIBRARIES )
  else()
    set( _cray_libsci_loaded $ENV{CRAY_LIBSCI_DIR} )
    ecbuild_add_option( FEATURE LAPACK
                        CONDITION _cray_libsci_loaded
                        DESCRIPTION "Linear Algebra PACKage")
    if( HAVE_LAPACK )
      set( _CRAY_PRGENV $ENV{PE_ENV} )
      string( TOLOWER "${_CRAY_PRGENV}" _cray_prgenv  )
      set( LAPACK_LIBRARIES sci_${_cray_prgenv} )
      ecbuild_debug( "LAPACK found, already loaded as part of Cray's libsci" )
    endif()
  endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
  ecbuild_add_cxx_flags( "-L /usr/local/apps/pgi/18.4/lib -lpgf90 -lpgf90_prof -lpgf902 -lpgf90_rpm1 -lpgf90rtl -lpgftnrtl -lpgmath -lpgnuma -lpgmp -lpgf90_rpm1_p -lpgc -lpgf90_prof -lpgf90_rpm1_p ")
  ecbuild_add_fortran_flags( "-pgc++libs -L /usr/local/apps/pgi/18.4/lib -lpgf90_prof -lpgf90_rpm1_p")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel" AND NOT BUILD_SHARED_LIBS)
  # This prevents CMake from (wrongly) using icpc as a linker for Fortran programs
  set(CMAKE_Fortran_LINKER_PREFERENCE_PROPAGATES True)
  set(CMAKE_C_LINKER_PREFERENCE_PROPAGATES False)
  set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES False)
endif()

if(multio_FOUND)
  set(MULTIO_LIBRARIES)
  if(multio_HAVE_FDB5)
    # Force the linkage of multio-fdb5 for static builds
    if(NOT BUILD_SHARED_LIBS)
      list(APPEND MULTIO_LIBRARIES -Wl,--push-state,--no-as-needed multio-fdb5 -Wl,--pop-state)
    else()
      list(APPEND MULTIO_LIBRARIES multio-fdb5)
    endif()
  endif()
  if(multio_HAVE_FDB4)
    list(APPEND MULTIO_LIBRARIES multio-fdb4 ifsio)
  endif()
  list(APPEND MULTIO_LIBRARIES multio)
endif()

include(ProcessorCount)
ProcessorCount(N)
if(N EQUAL 0)
  set(N 8)
endif()
set(FCM_PARALLEL ${N} CACHE STRING "Generate interfaces in parallel")

list(APPEND IFS_DEFINITIONS
  BLAS
  ECMWF
  INTEGER_IS_INT
  PATHSCALE
  ADDRESS64
  _ABI64
  F90
  PARAL
  NONCRAYF
  _RTTOV_DO_DISTRIBCOEF
  _RTTOV_NETCDF
  WITH_WAVE
  WITH_FCKIT
  UseCDF_CMF
)

if(HAVE_NEMO)
  list(APPEND IFS_DEFINITIONS WITH_NEMO)
endif()

if(HAVE_ATLAS)
  list(APPEND IFS_DEFINITIONS WITH_ATLAS)
endif()

if( HAVE_SINGLE_PRECISION ) 
  list(APPEND IFS_DEFINITIONS PARKIND1_SINGLE )
  set(PREC SP)
else()
  set(PREC DP)
endif()
string( TOLOWER "${PREC}" prec )



if( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  if($ENV{CRAY_FTN_VERSION} VERSION_LESS 8.7.2) # cdt/18.07
    list(APPEND IFS_DEFINITIONS WITH_IEEE_HALT)
  endif()
endif()


if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  list(APPEND IFS_DEFINITIONS LINUX)
endif()

if(EC_LITTLE_ENDIAN)
 list(APPEND IFS_DEFINITIONS LITTLE)
endif()

if(NOT HAVE_MPI)
  set(MPI_Fortran_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/mpi_serial/source)
  set(MPI_Fortran_LIBRARIES mpi_serial)
endif()

iF (ENABLE_COMPO_DRHOOK)
  list(APPEND IFS_DEFINITIONS WITH_COMPO_DR_HOOK )
endif()

if( HAVE_MGRIDS )
  list(APPEND IFS_DEFINITIONS WITH_MGRIDS )
  add_subdirectory( atlas_apps )
endif()
include(dummy)
include(fiat)
include(intfb)
include(nemo)
include(ifsaux)
include(algor)
include(surf)
include(climfield)
include(trans)
include(wam)
include(prepdata)
include(arpifs)
include(python)

if ( ENABLE_SCMEC )
  include(scmec)
endif()

# Link mpi_serial library to targets with SERIAL property.
if(HAVE_MPI)
  foreach(target ${${PROJECT_NAME}_ALL_EXES})
    get_target_property(serial ${target} SERIAL)
    if(serial)
      target_link_libraries(${target} ${MPI_SERIAL_LIBRARIES})
    endif()
  endforeach()
endif()

# Tell linker to generate target-specfic link map files.
foreach(target ${${PROJECT_NAME}_ALL_EXES})
  if( NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_property(TARGET ${target} APPEND_STRING PROPERTY LINK_FLAGS " -Wl,-Map=load.map.${target}")
  endif()
endforeach()

if(NOT HAVE_SINGLE_PRECISION)
  include(autopromote)
endif()

# Print compile times using COMPILER_LAUNCHER properties.
# Prepare output file:
set(timing_file ${CMAKE_BINARY_DIR}/compile_times.log)
string(TIMESTAMP tstamp UTC)
file(APPEND ${timing_file} "===== ${tstamp} =====\n")
foreach(target ${EC_ALL_EXES} ${EC_ALL_LIBS})
  set( compiler_launcher /usr/bin/time -ao ${timing_file} -f "(%e s) %C" )
  set_property(TARGET ${target} PROPERTY Fortran_COMPILER_LAUNCHER ${compiler_launcher})
  set_property(TARGET ${target} PROPERTY C_COMPILER_LAUNCHER ${compiler_launcher})
  set_property(TARGET ${target} PROPERTY CXX_COMPILER_LAUNCHER ${compiler_launcher})
endforeach()

### finalize project

# Set any remaining variables for ecBuild to pass via ifs-import.cmake.in:
set( IFS_EXECUTABLE               ${CMAKE_BINARY_DIR}/bin/ifsMASTER.${PREC} )
set( SCM_EXECUTABLE               ${CMAKE_BINARY_DIR}/bin/MASTER_scm.${PREC} )
set( IFS_BLACKLIST_LIBRARY_PATH   ${PROJECT_BINARY_DIR} )

ecbuild_install_project( NAME ${PROJECT_NAME} )

ecbuild_print_summary()

# ifs-test will do find_package( ifs ), therefore, also export ifs_dp as ifs
if( NOT ifs_DIR )
  set( ifs_DIR        ${${PROJECT_NAME}_BINARY_DIR} CACHE INTERNAL "alias to ${PROJECT_NAME}_BINARY_DIR" )
  set( ifs_BINARY_DIR ${${PROJECT_NAME}_BINARY_DIR} CACHE INTERNAL "alias to ${PROJECT_NAME}_BINARY_DIR" )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake ${PROJECT_BINARY_DIR}/ifs-config.cmake )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake ${PROJECT_BINARY_DIR}/ifs-config-version.cmake )
endif()
