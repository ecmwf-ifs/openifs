# (C) Copyright 1996-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

# - Try to find ecRegrid includes and library
# Once done this will define
#  ECREGRID_FOUND        - System has ecRegrid
#  ECREGRID_INCLUDE_DIRS - The ecRegrid include directories
#  ECREGRID_LIBRARIES    - The libraries needed to use ecRegrid
#
# also defined internally:
#  ECREGRID_LIBRARY, where to find the ecregrid library
#  ECREGRID_INCLUDE_DIR, where to find the ecregrid/ecregrid_api.h header

find_path(ECREGRID_INCLUDE_DIR ecregrid/ecregrid_api.h
  HINTS ${ECREGRID_PATH}/include $ENV{ECREGRID_DIR}/include $ENV{IFS_INSTALL_DIR}/include
)

find_library(ECREGRID_LIBRARY ecregrid
  HINTS ${ECREGRID_PATH}/lib $ENV{ECREGRID_DIR}/lib $ENV{IFS_INSTALL_DIR}/lib)

set(ECREGRID_INCLUDE_DIRS ${ECREGRID_INCLUDE_DIR})
set(ECREGRID_LIBRARIES ${ECREGRID_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ecRegrid DEFAULT_MSG ECREGRID_INCLUDE_DIR ECREGRID_LIBRARY)

mark_as_advanced(ECREGRID_INCLUDE_DIR ECREGRID_LIBRARY)
