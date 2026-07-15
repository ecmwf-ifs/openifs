# (C) Copyright 2011- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

# =================
#
# This file sets up a Python virtual environment that is needed for various
# data handling operations.
# The corresponding Python executable is forwarded to the parent scope as
# IFSTEST_PYTHON3.
#


function( python_venv )

    # Detect only the system installed Python3 interpreter.
    set( Python3_FIND_VIRTUALENV STANDARD )
    find_package( Python3 COMPONENTS Interpreter REQUIRED )

    set( VENV_NAME "ifstest_venv" )
    set( VENV_PATH "${CMAKE_CURRENT_BINARY_DIR}/venv/${VENV_NAME}" )


    # Create a Python virtual environment in the current binary directory.
    ecbuild_info( "Create Python virtual environment ${VENV_NAME}" )

    execute_process( COMMAND ${CMAKE_COMMAND} -E make_directory ${VENV_PATH} )
    execute_process( COMMAND ${Python3_EXECUTABLE} -m venv ${VENV_PATH} )

    set( IFSTEST_PYTHON3 "${VENV_PATH}/bin/python3")
    # Propagate the IFSTEST_PYTHON3 variable to the parent scope such
    # that is available even outside the function.
    set( IFSTEST_PYTHON3 ${IFSTEST_PYTHON3} PARENT_SCOPE)

    # Install the missing Python packages.
    execute_process(
      COMMAND ${IFSTEST_PYTHON3} -m pip install click pyyaml requests
    )
   
endfunction()
