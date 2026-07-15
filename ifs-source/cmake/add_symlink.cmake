# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

# Creates a symlink to target's output file.
# If the symlink is a relative path, it's considered relative to target's output directory.

function(add_symlink target symlink)

  function(_get_target_location target_path target)
    # This function wraps the CMake function "get_target_property( <VAR> <TARGET> LOCATION )"
    # and disabled deprecation warnings of the policy change to OLD
    cmake_policy( PUSH )
    # Remove deprecation warnings in this scope
    if( NOT DEFINED CMAKE_WARN_DEPRECATED OR CMAKE_WARN_DEPRECATED )
      set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)
      set( restore ON )
    endif()
    cmake_policy( SET CMP0026 OLD ) # allow property LOCATION on target
    if( restore )
      set(CMAKE_WARN_DEPRECATED ON CACHE BOOL "" FORCE)
    endif()
    get_target_property(location ${target} LOCATION)
    set( ${target_path} ${location} PARENT_SCOPE )
    cmake_policy( POP )
  endfunction()

  cmake_parse_arguments(_ARG "NOINSTALL" "" "" ${ARGN})

  if(_ARG_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Unknown keywords given to add_symlink: \"${_ARG_UNPARSED_ARGUMENTS}\"")
  endif()

  _get_target_location(target_path ${target})
  get_filename_component(target_dir ${target_path} DIRECTORY)

  if(IS_ABSOLUTE ${symlink})
    set(symlink_path ${symlink})
  else()
    set(symlink_path ${target_dir}/${symlink})
  endif()

  get_filename_component(symlink_dir ${symlink_path} DIRECTORY)
  file(RELATIVE_PATH target_relative_path ${symlink_dir} ${target_path})

  add_custom_command(OUTPUT ${symlink_path}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${symlink_dir}
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${target_relative_path} ${symlink_path}
    COMMENT "Creating symlink ${symlink} to target ${target}")

  get_filename_component(symlink_name ${symlink_path} NAME)
  add_custom_target(_${symlink_name} ALL DEPENDS ${symlink_path})

  if(NOT ${_ARG_NOINSTALL})
    file(RELATIVE_PATH install_relative_dir ${CMAKE_BINARY_DIR}/bin ${symlink_dir})
    install(PROGRAMS ${symlink_path} DESTINATION ${INSTALL_BIN_DIR}/${install_relative_dir})
  endif()

endfunction()

