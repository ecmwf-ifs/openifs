# (C) Copyright 2011- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

#Sets:
# RT_LIB      = the library to link against

if( DEFINED REALTIME_PATH )
    find_library(RT_LIB rt PATHS ${REALTIME_PATH}/lib NO_DEFAULT_PATH )
endif()

find_library( RT_LIB rt )

mark_as_advanced( RT_LIB )

include(FindPackageHandleStandardArgs)
# Handle the QUIET and REQUIRED arguments and set REALTIME_FOUND to TRUE
# if all listed variables are TRUE
# Note: capitalisation of the package name must be the same as in the file name
find_package_handle_standard_args(Realtime  DEFAULT_MSG RT_LIB )
