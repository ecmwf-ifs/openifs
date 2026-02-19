# (C) Copyright 2020- ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

### Workaround to extract GIT_SHA1 from parent directory
if( NOT ${PROJECT_NAME}_GIT_SHA1 )
    get_filename_component( PARENT_DIR ${PROJECT_SOURCE_DIR} DIRECTORY )
    if( EXISTS ${PARENT_DIR}/.git )
        get_filename_component( PARENT_REPOSITORY_NAME ${PARENT_DIR} NAME_WE )
        get_git_head_revision( GIT_REFSPEC ${PROJECT_NAME}_GIT_SHA1 )
        string( SUBSTRING "${${PROJECT_NAME}_GIT_SHA1}" 0 7 ${PROJECT_NAME}_GIT_SHA1_SHORT )
        set( ${PROJECT_NAME}_GIT_SHA1_SHORT "${PARENT_REPOSITORY_NAME}/${${PROJECT_NAME}_GIT_SHA1_SHORT}" )
        set( ${PROJECT_NAME}_GIT_SHA1       "${PARENT_REPOSITORY_NAME}/${${PROJECT_NAME}_GIT_SHA1}" )
    endif()
endif()


try_run( execute_result compile_result
         ${CMAKE_CURRENT_BINARY_DIR}
         ${PROJECT_SOURCE_DIR}/cmake/test_attribute_constructor.c
         COMPILE_OUTPUT_VARIABLE compile_output
         RUN_OUTPUT_VARIABLE execute_output )

ecbuild_debug("Compiling and running ${PROJECT_SOURCE_DIR}/cmake/test_attribute_constructor.c")
ecbuild_debug_var( compile_result )
ecbuild_debug_var( compile_output )
ecbuild_debug_var( execute_result )
ecbuild_debug_var( execute_output )


if( NOT DEFINED FIAT_ATTRIBUTE_CONSTRUCTOR_SUPPORTED )
  set( FIAT_ATTRIBUTE_CONSTRUCTOR_SUPPORTED 0 )
  if( compile_result )
    if( execute_result MATCHES 0 )
      set( FIAT_ATTRIBUTE_CONSTRUCTOR_SUPPORTED 1 )
    else()
      ecbuild_info("Compiler failed to correctly run program with \"__attribute__((constructor))\".")
    endif()
  endif()
endif()

include( fiat_target_fortran_module_directory )
include( fiat_target_ignore_missing_symbols )

