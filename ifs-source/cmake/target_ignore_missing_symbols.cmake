# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

function( target_ignore_missing_symbols )
  set( options  )
  set( single_value_args TARGET )
  set( multi_value_args SYMBOLS )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )
if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  unset( link_lazy )
  foreach(symbol ${_PAR_SYMBOLS} )
    set(link_lazy "${link_lazy} -Wl,-U,${symbol}")
    set(link_lazy "${link_lazy} -Wl,-U,_${symbol}")
  endforeach()
  set_property( TARGET ${_PAR_TARGET} APPEND_STRING PROPERTY LINK_FLAGS " ${link_lazy}" )
endif()
endfunction()
