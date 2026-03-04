# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


if( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )
  ecbuild_add_fortran_flags("-qextname -qnobindcextname")
endif()

# gfortran 10 has become stricter with argument matching
if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU"
    AND NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10 )
  ecbuild_add_fortran_flags("-fallow-argument-mismatch")
endif()


macro( ectrans_add_compile_options )
  set( options  )
  set( single_value_args FLAGS )
  set( multi_value_args SOURCES )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )
  if(_PAR_UNPARSED_ARGUMENTS)
    ecbuild_critical("Unknown keywords given to ectrans_add_compile_flags(): \"${_PAR_UNPARSED_ARGUMENTS}\"")
  endif()
  if(NOT _PAR_SOURCES)
    ecbuild_critical("SOURCES keyword missing to ectrans_add_compile_flags()")
  endif()
  if(NOT _PAR_FLAGS)
    ecbuild_critical("FLAGS keyword missing to ectrans_add_compile_flags()")
  endif()
  foreach( _file ${_PAR_SOURCES} )
    ecbuild_warn("Adding custom compile flags for file ${_file} : [${_PAR_FLAGS}]")
    if( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_file} )
        ecbuild_error("${_file} does not exist")
    endif()
    set_source_files_properties( ${_file} PROPERTIES COMPILE_FLAGS "${_PAR_FLAGS}" )
  endforeach()
endmacro()

