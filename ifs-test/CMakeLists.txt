# (C) Copyright 1996-2015 ECMWF.

############################################################################################
# ifs_test project

cmake_minimum_required( VERSION 3.17 FATAL_ERROR )

find_package( ecbuild 3.1 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )

project( ifs_test )

find_package( ifs REQUIRED )
find_package( eccodes )

ecbuild_find_package( ifsbench )


ecbuild_info( "IFS_EXECUTABLE: ${IFS_EXECUTABLE}" )
ecbuild_info( "IFS_CYCLE: ${IFS_CYCLE}" )
ecbuild_info( "IFS_BLACKLIST_LIBRARY_PATH: ${IFS_BLACKLIST_LIBRARY_PATH}")
ecbuild_info( "IFS_PRECISION: ${IFS_PRECISION}")
ecbuild_info( "ECCODES_DEFINITION_PATH: ${ECCODES_DEFINITION_PATH}" )
ecbuild_info( "ECCODES_IFS_SAMPLES_PATH: ${ECCODES_IFS_SAMPLES_PATH}" )
ecbuild_info( "IFSBENCH_Python3_EXECUTABLE: ${IFSBENCH_Python3_EXECUTABLE}" )

#################################################################
# Configuration


set( IFSTEST_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( IFSTEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/bin )
file(COPY 
  ${CMAKE_CURRENT_SOURCE_DIR}/bin 
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

set( SHARE_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/share )
set( SHARE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/ifs_test )

execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${SHARE_BINARY_DIR} )
set( IFS_NAMELIST ${SHARE_SOURCE_DIR}/namelist )

if( IFS_CHECK_BOUNDS)
  set( IFS_NUNDEFLD "1")
else()
  set( IFS_NUNDEFLD "-99999")
endif()

if( CMAKE_SYSTEM_NAME MATCHES "Darwin")
  set( LD_LIBRARY_PATH_VAR "DYLD_LIBRARY_PATH" )
else()
  set( LD_LIBRARY_PATH_VAR "LD_LIBRARY_PATH" )
endif()

configure_file(${SHARE_SOURCE_DIR}/ifs_env.in ${SHARE_BINARY_DIR}/ifs_env )

# Setup the Python virtual env for fetching data.
include( setup_python )
python_venv()

set( CONFIG_FILE ${PROJECT_SOURCE_DIR}/storage.yaml )

# Put initial data into the build directory.
# Putting this in the source directory (or any other shared place) may cause
# race conditions when running two builds in parallel.
set( DATA_CACHE_DIR ${PROJECT_BINARY_DIR}/.cache )

# Allow user to set a custom path from which initial data gets pulled from.
set( IFS_TEST_CUSTOM_DATA_SOURCE "" CACHE PATH "Custom inidata path.")

# By default, force-fetch initial data when configuring.
set( IFS_TEST_FORCE_FETCH ON CACHE BOOL "Force update of initial data.")

if ( "${IFS_TEST_FORCE_FETCH}" )
  set( STORAGE_FETCH_FLAGS "--force" )
else()
  set( STORAGE_FETCH_FLAGS "" )
endif()

if ( NOT "${IFS_TEST_CUSTOM_DATA_SOURCE}" STREQUAL "" )
  set( STORAGE_FETCH_FLAGS ${STORAGE_FETCH_FLAGS} --custom-path ${IFS_TEST_CUSTOM_DATA_SOURCE} )
endif()

execute_process(
  COMMAND ${IFSTEST_PYTHON3} ${PROJECT_SOURCE_DIR}/bin/storage.py -vv fetch
    ${STORAGE_FETCH_FLAGS} ${CONFIG_FILE} ${DATA_CACHE_DIR}
  COMMAND_ECHO STDOUT
  RESULT_VARIABLE FETCH_EXIT_CODE
)

if( NOT (${FETCH_EXIT_CODE} EQUAL 0) )
  ecbuild_critical("Fetching the initial data failed with return code "
    "${FETCH_EXIT_CODE}.")
endif()

# Little bit hacky. Improve, if you find a better solution. We (ab-)use
# IFSBENCH_RUNNER_LIST as a global variable by putting it in the cache.
# Initially, it is empty but will in the end contain all the ifs-run
# subscripts.
set(IFSBENCH_RUNNER_LIST "" CACHE INTERNAL "")

# Actually update the value (set doesn't necessarily do this for cached
# variables).
mark_as_advanced(IFSBENCH_RUNNER_LIST)

# Copy the config file which is required by the runners.
file(COPY
  ${CMAKE_CURRENT_SOURCE_DIR}/tests/common/config.py
  DESTINATION ${IFSTEST_BINARY_DIR}/ifsbench_runner/
)

#################################################################
# Tests

add_subdirectory(tests)

if( ifsbench_FOUND )
  ecbuild_debug( "IFSBENCH_RUNNER_LIST: ${IFSBENCH_RUNNER_LIST}" )
endif()
configure_file( ${SHARE_SOURCE_DIR}/ifs-run.py ${IFSTEST_BINARY_DIR}/ifs-run.py )

