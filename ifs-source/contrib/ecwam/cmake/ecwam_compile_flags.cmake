# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  set(autopromote_flags   "-sreal64")
  set(checkbounds_flags   "-Rb")
  set(fpe_flags           "-Ktrap=fp")
  set(initsnan_flags      "-ei")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(autopromote_flags   "-fdefault-real-8 -fdefault-double-8")
  set(checkbounds_flags   "-fcheck=bounds")
  set(fpe_flags           "-ffpe-trap=invalid,zero,overflow")
  set(initsnan_flags      "-finit-real=snan")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(autopromote_flags   "-real-size 64")
  set(checkbounds_flags   "-check bounds")
  set(initsnan_flags      "-init=snan")
  set(fpe_flags           "-fpe0")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
  set(autopromote_flags   "-r8")
  set(fpe_flags           "-Ktrap=fp")
#  set(checkbounds_flags   "-Mbounds") # Added by default by CMake in NVHPC debug builds

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
  set(autopromote_flags   "-fdefault-real-8")
  set(fpe_flags           "-ffp-exception-behavior=strict")

endif()

if( NOT HAVE_SINGLE_PRECISION )
  ecbuild_add_fortran_flags( "${autopromote_flags}"   NAME autopromote )
endif()

## Debug flags for NVHPC are applied selectively to sourcefiles in src/ecwam/CMakeLists.txt
if( CMAKE_BUILD_TYPE MATCHES "Debug" AND NOT CMAKE_Fortran_COMPILER_ID MATCHES PGI|NVHPC )
  foreach( debug_flag    fpe initsnan checkbounds )
    if( ${debug_flag}_flags )
      ecbuild_add_fortran_flags( "${${debug_flag}_flags}" NAME ${debug_flag} )
    endif()
  endforeach()
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    # In case '-check all' has been added, we need to remove the '-check arg_temp_created' warnings
    ecbuild_add_fortran_flags( "-check noarg_temp_created" NAME check_noarg_temp_created BUILD DEBUG ) # the BUILD DEBUG argument makes sure it is appended after '-check all'
  endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  if( NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10 )
    ecbuild_add_fortran_flags( "-fallow-argument-mismatch" NAME argument_mismatch )
  endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
  # Linker complains of unknown arguments:
  #    warning: argument unused during compilation: '-fdefault-real-8' [-Wunused-command-line-argument]
  foreach( LINKER_FLAGS CMAKE_EXE_LINKER_FLAGS CMAKE_SHARED_LINKER_FLAGS CMAKE_STATIC_LINKER_FLAGS )
    set( ${LINKER_FLAGS} "${${LINKER_FLAGS}} -Wno-unused-command-line-argument")
  endforeach()
endif()
