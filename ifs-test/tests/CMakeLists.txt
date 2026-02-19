# (C) Copyright 2011- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

# Path to the ifsdata directory on disk.
set (STATIC_DATA_PATH "${DATA_CACHE_DIR}/ifsdata")

set (COMPO_DATA_PATH "${DATA_CACHE_DIR}/compo_data")

add_subdirectory(t21)
add_subdirectory(t42)
if( ENABLE_SCMEC )
  add_subdirectory(scm)
endif()


