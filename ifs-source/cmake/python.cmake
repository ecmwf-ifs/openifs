# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

ecbuild_find_python(VERSION 3.6 REQUIRED NO_LIBS)

# Loop over all python packages and install them
file( GLOB python_packages RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/python python/* )
foreach( python_package ${python_packages} )

  # Copy Python source code to the build directory.
  file( COPY python/${python_package} DESTINATION "./python" )

  # Install the python package using pip, targeting the lib/python directory
  # of the install.
  install( CODE "message( \"Installing Python package ${python_package}:
      ${PYTHON_EXECUTABLE} -m pip install --target=${CMAKE_INSTALL_PREFIX}/lib/python\" )
                 execute_process( COMMAND ${PYTHON_EXECUTABLE} -m pip install --target=${CMAKE_INSTALL_PREFIX}/lib/python ./${python_package}
                                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/python )" )

endforeach()
