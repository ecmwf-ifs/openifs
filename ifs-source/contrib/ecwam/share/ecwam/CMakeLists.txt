# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

set( destination share/${PROJECT_NAME} )

###############################################################################################
# Files 
###############################################################################################

list( APPEND files parameters.yml )

install(
    FILES       ${files}
    DESTINATION ${destination}
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_WRITE )

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${destination})

foreach( file ${files} )
  execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
      "${CMAKE_CURRENT_SOURCE_DIR}/${file}"
      "${CMAKE_BINARY_DIR}/${destination}/${file}")
endforeach()

###############################################################################################
# scripts
###############################################################################################

file (GLOB files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} scripts/*.sh scripts/*.py scripts/ecwam-*)

install(
    FILES       ${files}
    DESTINATION ${destination}/scripts
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE WORLD_EXECUTE OWNER_WRITE )

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${destination}/scripts)

foreach( file ${files} )
  execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
      "${CMAKE_CURRENT_SOURCE_DIR}/${file}"
      "${CMAKE_BINARY_DIR}/${destination}/${file}")
endforeach()

###############################################################################################
# run scripts symlinks
###############################################################################################

foreach( program preproc preset model )
  execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
    "../share/${PROJECT_NAME}/scripts/ecwam_run_${program}.sh"
    "${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}-run-${program}"
  )
  install(FILES ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}-run-${program} DESTINATION bin)
endforeach()

