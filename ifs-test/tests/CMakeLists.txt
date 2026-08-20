# (C) Copyright 2011- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

set( IFS_TEST_TIMEOUT -1 CACHE STRING "ifs-test timeout limit (time in seconds, negative value for defaults)" )
if( NOT ${IFS_TEST_TIMEOUT} MATCHES ^-?[0-9]+$ )
  ecbuild_error("IFS_TEST_TIMEOUT should be an integer value")
endif()

# Set timeout property for ifs-tests
if( IFS_TEST_TIMEOUT GREATER_EQUAL 0 )
  # Set local variable from cache variable
  set(IFS_TEST_TIMEOUT ${IFS_TEST_TIMEOUT})
else()
  # Use the provided defaults, based on build settings
  if( IFS_CHECK_BOUNDS )
    set(IFS_TEST_TIMEOUT $<IF:$<CONFIG:Bit>,300,600>)
  else()
    set(IFS_TEST_TIMEOUT $<IF:$<CONFIG:Bit>,120,360>)
  endif()
endif()

# Path to the ifsdata directory on disk.
set (STATIC_DATA_PATH "${DATA_CACHE_DIR}/ifsdata")

set (COMPO_DATA_PATH "${DATA_CACHE_DIR}/compo_data")

add_subdirectory(t21)
add_subdirectory(t42)

if( ENABLE_SCMEC )
  add_subdirectory(scm)
endif()

if( ifsbench_FOUND )
  add_subdirectory(ifsbench)
endif()

