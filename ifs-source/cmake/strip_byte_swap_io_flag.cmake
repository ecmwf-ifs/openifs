# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

# CMake script to strip byte-swap I/O flag from the link line

if(CMAKE_GENERATOR MATCHES "Makefiles")
  set(filename "${CMAKE_BINARY_DIR}/${PROJECT_NAME}/CMakeFiles/${TARGET}.dir/link.txt")
elseif(CMAKE_GENERATOR MATCHES "Ninja")
  set(filename "${CMAKE_BINARY_DIR}/${PROJECT_NAME}/CMakeFiles/${TARGET}.rsp")
else()
  message(FATAL_ERROR "Unsupported generator: ${CMAKE_GENERATOR}")
endif()

file(STRINGS "${filename}" link_line NEWLINE_CONSUME)

if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  string(REGEX REPLACE "-h[ ]*byteswapio" "" link_line "${link_line}")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  string(REGEX REPLACE "-fconvert=big-endian" "" link_line "${link_line}")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  string(REGEX REPLACE "-convert[ ]+big_endian" "" link_line "${link_line}")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  string(REGEX REPLACE "-[M]?byteswapio" "" link_line "${link_line}")
else()
  message(FATAL_ERROR "Unsupported Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

file(WRITE "${filename}" "${link_line}")
