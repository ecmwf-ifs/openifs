# (C) Copyright 2011- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

# Path to the ifsbench inidata directories.
set( IFSBENCH_TL159_DATA "${DATA_CACHE_DIR}/tl159_ifsbench" )
set( IFSBENCH_TCO399_DATA "${DATA_CACHE_DIR}/tco399_ifsbench" )

include( add_ifsbench_test )

# NOTE: 49R1 tl159/tco399 ifsbench inidata has not been uploaded to ECPDS/Atos yet, so
# these large tests are not registered here (registering them would trigger a fetch
# for data that doesn't exist and fail configure). Once storage.yaml gains
# tl159_ifsbench/tco399_ifsbench entries for 49r1, uncomment this list to re-enable.
# list(APPEND EXP_FILES
#   tl159.yml
#   tco399.yml
# )

foreach(EXP_FILE ${EXP_FILES})
  add_ifsbench_test(
    FILEPATH "${CMAKE_CURRENT_SOURCE_DIR}/${EXP_FILE}"
    IS_LARGE true
    LABELS "IFSBENCH;LARGE"
  )
endforeach()

include( add_ifsbench_runner )
add_ifsbench_runner(
  FILEPATH from_yaml.py
  DO_CONFIGURE true
)
