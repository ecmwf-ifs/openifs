# (C) Copyright 1996-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

# - Try to find ParMETIS includes and library
# Once done this will define
#  PARMETIS_FOUND        - System has ParMETIS
#  PARMETIS_INCLUDE_DIRS - The ParMETIS include directories
#  PARMETIS_LIBRARIES    - The libraries needed to use ParMETIS
#
# also defined internally:
#  PARMETIS_LIBRARY, where to find the ParMETIS library
#  PARMETIS_INCLUDE_DIR, where to find the parmetis.h header
#  METIS_LIBRARY, where to find the METIS library

find_path(PARMETIS_INCLUDE_DIR parmetis.h
    HINTS ${PARMETIS_PATH}/include $ENV{PARMETIS_DIR}/include $ENV{PARMETIS_ROOT}/include)

find_library(PARMETIS_LIBRARY parmetis
    HINTS ${PARMETIS_PATH}/lib $ENV{PARMETIS_DIR}/lib $ENV{PARMETIS_ROOT}/lib)

find_library(METIS_LIBRARY metis
    HINTS ${PARMETIS_PATH}/lib $ENV{PARMETIS_DIR}/lib $ENV{PARMETIS_ROOT}/lib)

set(PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR})
set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ParMETIS DEFAULT_MSG PARMETIS_INCLUDE_DIR PARMETIS_LIBRARY METIS_LIBRARY)

mark_as_advanced(PARMETIS_INCLUDE_DIR PARMETIS_LIBRARY METIS_LIBRARY)
