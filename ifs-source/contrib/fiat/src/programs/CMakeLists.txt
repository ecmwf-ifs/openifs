# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

### fiat information tool

set( tool fiat )

get_property( langs GLOBAL PROPERTY ENABLED_LANGUAGES )

foreach( lang ${langs} )
  set( EC_${lang}_FLAGS "${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${CMAKE_BUILD_TYPE_CAPS}}" )
endforeach()

configure_file( ${tool}.in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${tool} @ONLY )

file(COPY ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${tool}
  DESTINATION ${CMAKE_BINARY_DIR}/bin
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(
  FILES  ${CMAKE_BINARY_DIR}/bin/${tool}
  DESTINATION   ${INSTALL_BIN_DIR}
  PERMISSIONS   OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

### drhook-merge tool

list( APPEND files
  drhook-merge
)

install(
  FILES       ${files}
  DESTINATION ${INSTALL_BIN_DIR}
  PERMISSIONS   OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)
foreach( file ${files} )
  execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
      "${CMAKE_CURRENT_SOURCE_DIR}/${file}"
      "${CMAKE_BINARY_DIR}/${INSTALL_BIN_DIR}/${file}")
endforeach()


### Program fiat-printbinding
# Note: it only depends on MPI and OpenMP
ecbuild_add_executable(TARGET fiat-printbinding
  SOURCES fiat-printbinding.c )
find_package( OpenMP COMPONENTS C )
if( TARGET OpenMP::OpenMP_C )
  target_link_libraries( fiat-printbinding OpenMP::OpenMP_C )
endif()
find_package( MPI COMPONENTS C )
if( TARGET MPI::MPI_C )
  target_link_libraries( fiat-printbinding MPI::MPI_C )
else()
  target_compile_definitions( fiat-printbinding PRIVATE NOMPI )
endif()
