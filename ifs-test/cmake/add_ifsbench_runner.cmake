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
# add_ifsbench_runner
# ===================
#
# Add an ifsbench runner script to ``${IFSBENCH_RUNNER_LIST}``::
#
#   add_ifsbench_runner( FILEPATH script_path DO_CONFIGURE do_configure )
#
# This adds the Python script to run a test case to the list of runner
# scripts, which is used to make them available in the common ifs-run.py
# executable.
# The script is copied to ${IFSTEST_BINARY_DIR}/ifsbench_runner. If DO_CONFIGURE
# is set to true, the file is not only copied but also configured (using the
# CMake configure_file functionality.
#
#
# Example
# -------
#
# code-block::
#   add_ifsbench_runner(FILEPATH ifs_run.py DO_CONFIGURE false )
#
##############################################################################

function( add_ifsbench_runner )
    set( oneValueArgs FILEPATH )
    set( options DO_CONFIGURE)

    cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( _PAR_DO_CONFIGURE)
      set( DO_CONFIGURE ${_PAR_DO_CONFIGURE} )
    else()
      set( DO_CONFIGURE false)
    endif()

    if( NOT _PAR_FILEPATH )
      ecbuild_critical( "No FILEPATH specified for add_ifsbench_test()" )
    endif()

    set( FILEPATH ${_PAR_FILEPATH} )

    # The path to which we will copy all run scripts.
    set( SCRIPT_DIR ${IFSTEST_BINARY_DIR}/ifsbench_runner)

    # Actually create the SCRIPT_DIR directory if it doesn't exist yet.
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${SCRIPT_DIR})

    # Get the name of the script (just the name, not the directory).
    get_filename_component(FILENAME ${FILEPATH} NAME)

    if( DO_CONFIGURE)
      configure_file(
        ${FILEPATH}
        ${SCRIPT_DIR}/${FILENAME}
        @ONLY
      )
    else()
      file(COPY
        ${FILEPATH}
        DESTINATION ${SCRIPT_DIR}/${FILENAME}
      )

    endif()

    # Get the name of the Python script without the .py suffix. We need this
    # here, as we will later do a "import FILENAME_NO_SUFFIX" in one of the
    # Python scripts.
    get_filename_component(FILENAME_NO_SUFFIX ${FILEPATH} NAME_WE)

    # Just a safety layer. Check that a file with the same name hasn't been
    # added to the runner list yet. We could probably circumvent this problem
    # by adding some kind of unique ID to each script, so that we could import
    # scripts from different places but with the same name. 
    list( FIND IFSBENCH_RUNNER_LIST ${FILENAME_NO_SUFFIX} ENTRY_FOUND)
    if(NOT ENTRY_FOUND EQUAL -1)
      message(FATAL_ERROR "Found two ifsbench runner scripts with the same name (${FILENAME}).")
    endif()


    list( APPEND IFSBENCH_RUNNER_LIST ${FILENAME_NO_SUFFIX} )
    set(IFSBENCH_RUNNER_LIST ${IFSBENCH_RUNNER_LIST} CACHE INTERNAL "")
    mark_as_advanced(IFSBENCH_RUNNER_LIST)
endfunction()
