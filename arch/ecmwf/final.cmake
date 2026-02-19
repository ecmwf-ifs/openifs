# (C) Copyright 2011- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

# Try to guess the module install prefix and the module version from CMAKE_INSTALL_PREFIX.

set( MODULE_NAME ifs CACHE STRING "Module name (as in `module load <MODULE_NAME>`)." )

set( prefix ${CMAKE_INSTALL_PREFIX} )
get_filename_component( name ${prefix} NAME )
get_filename_component( prefix ${prefix} DIRECTORY )

while( TRUE )
  set( version ${name} )
  get_filename_component( name ${prefix} NAME )
  get_filename_component( prefix ${prefix} DIRECTORY )
  if ( NOT name )
    set( prefix )
    set( version )
    break()
  elseif( name STREQUAL MODULE_NAME )
    break()
  endif()
endwhile()

set( MODULE_VERSION ${version} CACHE STRING "Module version (as in `module load <MODULE_NAME>/<MODULE_VERSION>`)." )
set( MODULE_PREFIX ${prefix} CACHE FILEPATH "Path where module packages are installed (e.g. /usr/local/apps)." )

if( MODULE_PREFIX AND MODULE_NAME AND MODULE_VERSION )
  configure_file( ${CMAKE_CURRENT_BINARY_DIR}/modulefile.in ${MODULE_VERSION} @ONLY )
  install( FILES ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_VERSION} DESTINATION share/modulefiles/${MODULE_NAME} )
endif()
