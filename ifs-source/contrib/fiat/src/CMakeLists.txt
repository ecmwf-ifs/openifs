# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

add_subdirectory( mpi_serial )
add_subdirectory( fiat )
add_subdirectory( parkind )
add_subdirectory( programs )

foreach( filepath ${fiat_src} parkind1.F90 parkind2.F90 )
  get_filename_component( filename ${filepath} NAME )
  list( APPEND fiat_SOURCE_FILENAMES ${filename} )
endforeach()
set( fiat_SOURCE_FILENAMES ${fiat_SOURCE_FILENAMES} PARENT_SCOPE )

# needs to be available for tests
set( FIAT_DEFINITIONS ${FIAT_DEFINITIONS} PARENT_SCOPE )
