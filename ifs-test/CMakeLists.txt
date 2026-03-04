# (C) Copyright 1996-2015 ECMWF.

############################################################################################
# ifs_test project

cmake_minimum_required( VERSION 3.6 FATAL_ERROR )

find_package( ecbuild 3.1 REQUIRED )

project( ifs_test )

find_package( ifs REQUIRED )
find_package( eccodes )

ecbuild_info( "IFS_EXECUTABLE: ${IFS_EXECUTABLE}" )
ecbuild_info( "IFS_CYCLE: ${IFS_CYCLE}" )
ecbuild_info( "IFS_BLACKLIST_LIBRARY_PATH: ${IFS_BLACKLIST_LIBRARY_PATH}")
ecbuild_info( "IFS_PRECISION: ${IFS_PRECISION}")
ecbuild_info( "ECCODES_DEFINITION_PATH: ${ECCODES_DEFINITION_PATH}" )
ecbuild_info( "ECCODES_IFS_SAMPLES_PATH: ${ECCODES_IFS_SAMPLES_PATH}" )

#################################################################
# Configuration

set( SHARE_BINARY_DIR ${CMAKE_BINARY_DIR}/share/ifs_test )
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


# Detect only the system installed Python3 interpreter.
set( Python3_FIND_VIRTUALENV STANDARD )
find_package( Python3 3.8 COMPONENTS Interpreter REQUIRED )

set( VENV_NAME "test_venv" )
set( VENV_PATH "${CMAKE_CURRENT_BINARY_DIR}/venv/${VENV_NAME}" )


# Create a Python virtual environment in the current binary directory.
ecbuild_info( "Create Python virtual environment ${VENV_NAME}" )

execute_process( COMMAND ${CMAKE_COMMAND} -E make_directory ${VENV_PATH} )
execute_process( COMMAND ${Python3_EXECUTABLE} -m venv ${VENV_PATH} )

set( IFSTEST_PYTHON3 "${VENV_PATH}/bin/python3")

# Install the missing Python packages.
execute_process(
  COMMAND ${IFSTEST_PYTHON3} -m pip install click pyyaml requests
)

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

# Tests

add_subdirectory(tests)

