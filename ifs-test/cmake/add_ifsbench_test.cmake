# (C) Copyright 2011- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

##############################################################################
#.rst:
#
# add_ifsbench_test
# ===================
#
# Add an ifsbench test script to the test system.
#
#   add_ifsbench_test( FILENAME IS_LARGE LABELS )
#
# :param FILEPATH: The path to the experiment .yaml file. 
# :param IS_LARGE: Either "true" or "false" (defaults to "false"). 
#                  If "true", the test is designated as "large" and will not
#                  run unless the "IFSTEST_RUN_LARGE" environment variable is 
#                  set to 1.
# :param LABELS:   A list of labels that is added to the test in CTest.
#
# This has the following effects:
#
# * A directory is created for this test at ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}.
# * The experiment file is copied to this directory.
# * The runner script is configured and copied to this directory.
# * The test itself is registered in ecbuild/ctest, using the given labels.
#
##############################################################################

function( add_ifsbench_test )
  set( oneValueArgs FILEPATH )
  set( options IS_LARGE )
  set( multiValueArgs LABELS )

  cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if( _PAR_IS_LARGE)
    set( IS_LARGE ${_PAR_IS_LARGE} )
  else()
    set( IS_LARGE false)
  endif()

  if( NOT _PAR_FILEPATH )
    ecbuild_critical( "No FILEPATH specified for add_ifsbench_test()" )
  endif()

  set( FILEPATH ${_PAR_FILEPATH} )
  set( LABELS ${_PAR_LABELS} )

  get_filename_component(TESTNAME ${FILEPATH} NAME_WE)

  set(RUN_DIR "${CMAKE_CURRENT_BINARY_DIR}/${TESTNAME}")
  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${RUN_DIR})

  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink 
    ${IFSTEST_BINARY_DIR}/ifs-run.py ${RUN_DIR}/ifs-run.py)

  configure_file(
    "${SHARE_SOURCE_DIR}/ifsbench_runner.sh.in"
    "${RUN_DIR}/ifsbench_runner.sh"
    @ONLY
  )

  # At least for CTest we wrap the call to ifs-run.py in a wrapper script to
  # handle bit-identicality checks in the same way as in the other tests.
  set(RUN_CMD "${RUN_DIR}/ifsbench_runner.sh")

  list(APPEND RUN_OPTIONS 
    "--run-dir=${RUN_DIR}"
  )

  # Set NEXUS_CACHE_DIR (which was the old name for INIDATA_CACHE_DIR) for
  # backwards compatibility.
  set(NEXUS_CACHE_DIR "${INIDATA_CACHE_DIR}")

  # Configure the actual experiment file. This essentially sets the path to the
  # local cache directory.
  configure_file( 
    ${FILEPATH}
    ${RUN_DIR}
    @ONLY
  )

  # Use the get_exp_relpath script to extract the relative path of the inidata
  # files from the experiment file.
  execute_process(
    COMMAND ${IFSTEST_SITES_PYTHON3} ${IFSTEST_CMAKE_DIR}/get_exp_relpath.py ${FILEPATH}
    OUTPUT_VARIABLE REL_PATHS
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # For each required file, check if it exists and download it, if necessary.
  foreach(REL_PATH ${REL_PATHS})
    if(EXISTS "${INIDATA_CACHE_DIR}/${REL_PATH}")
      continue()
    endif()

   file(DOWNLOAD "${INIDATA_URL}/${REL_PATH}" "${INIDATA_CACHE_DIR}/${REL_PATH}")
  endforeach()

  ecbuild_add_test(
    TARGET "ifs_ifsbench_test_${TESTNAME}"
    COMMAND "${RUN_CMD}"
    ARGS "${RUN_DIR}/${TESTNAME}.yml" ${RUN_OPTIONS}
    LABELS ${LABELS}
    WORKING_DIRECTORY ${RUN_DIR}
  )


endfunction()
