# (C) Copyright 2022 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

find_package(FORD QUIET)
ecbuild_add_option(
  FEATURE DOCS
  DESCRIPTION "Generate reference documentation"
  CONDITION FORD_FOUND )

if( HAVE_DOCS )

  set( ECWAM_DOC ecwam_doc )
  if( PROJECT_NAME STREQUAL CMAKE_PROJECT_NAME )
    set( ECWAM_DOC doc )
  endif()

  if( NOT TARGET ${ECWAM_DOC} )

  set( FORDFILE ford.md CACHE INTERNAL "FORD project filename" )
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${FORDFILE}.in "${CMAKE_CURRENT_BINARY_DIR}/${FORDFILE}" @ONLY )

  file( GLOB_RECURSE FORD_DEPENDS ${PROJECT_SOURCE_DIR}/src/ecwam/* )
  list( APPEND FORD_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${FORDFILE}.in )

  add_custom_command(OUTPUT "ford/index.html"
    COMMAND ${FORD_EXECUTABLE} --debug "${CMAKE_CURRENT_BINARY_DIR}/${FORDFILE}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    DEPENDS ${FORD_DEPENDS}
    COMMENT "Building HTML documentation for ${CMAKE_PROJECT_NAME} using FORD (${CMAKE_CURRENT_BINARY_DIR}/ford/index.html)" )
  add_custom_target( ${ECWAM_DOC} DEPENDS "ford/index.html" ${FORD_DEPENDS} )

  endif()
endif()


