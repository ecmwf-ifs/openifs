# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

add_custom_target( ${PROJECT_NAME}_cache_bathymetry )
set_property(GLOBAL PROPERTY JOB_POOLS may_download=1)

function( ecwam_cache_bathymetry )
    set( configs ${ARGV} )

    unset( LAUNCH )
    if( CMAKE_CROSSCOMPILING_EMULATOR )
      set( LAUNCH ${CMAKE_CROSSCOMPILING_EMULATOR} )
    endif()

    foreach(config ${configs} )

      if( NOT IS_ABSOLUTE "${config}")
        set(config "${CMAKE_CURRENT_SOURCE_DIR}/${config}")
      endif()

      file(RELATIVE_PATH config_relpath "${PROJECT_SOURCE_DIR}" "${config}")
      STRING( REPLACE ".yml" "" config_relpath ${config_relpath} )

      set( binary_dir ${PROJECT_BINARY_DIR}/ecwam_cache_bathymetry/${config_relpath} )
      set( run_dir ${binary_dir} )
      if( ECWAM_WORKSPACE )
        set( run_dir ${ECWAM_WORKSPACE}/ecwam_cache_bathymetry/${config_relpath} )
      endif()
      if( ECWAM_WORKSPACE )
        get_filename_component(binary_dir_parent ${binary_dir} DIRECTORY)
        file( MAKE_DIRECTORY  ${binary_dir_parent} )
        execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${run_dir}" "${binary_dir}")
      endif()
      file( MAKE_DIRECTORY  ${run_dir} )
      execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${config}" "${run_dir}/config.yml")

      # Set environment
      file( WRITE  ${run_dir}/run.sh "#!/usr/bin/env bash\n" )
      file( APPEND ${run_dir}/run.sh "export PATH=${CMAKE_BINARY_DIR}/bin:\$PATH\n" )
      if( CMAKE_CROSSCOMPILING_EMULATOR )
      file( APPEND ${run_dir}/run.sh "export LAUNCH=\"\${LAUNCH:-${CMAKE_CROSSCOMPILING_EMULATOR}}\"\n" )
      endif()
      file( APPEND ${run_dir}/run.sh "export OMP_NUM_THREADS=1\n" )
      file( APPEND ${run_dir}/run.sh "${CMAKE_BINARY_DIR}/share/${ECWAM_PROJECT_NAME}/scripts/ecwam_run_create_bathymetry.sh --run-dir=${run_dir} \"\${@}\" \n" )
      file( APPEND ${run_dir}/run.sh "touch ${run_dir}/done" )
      file( CHMOD  ${run_dir}/run.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )

      add_custom_command( OUTPUT ${run_dir}/done 
        COMMAND ${run_dir}/run.sh
        WORKING_DIRECTORY ${run_dir}
        JOB_POOL may_download
        DEPENDS ${PROJECT_NAME}-create_wam_bathymetry ${PROJECT_NAME}-create_wam_bathymetry_ETOPO1 ${config}
      )

      STRING( REPLACE "/" "." config_relpath ${config_relpath} )
      add_custom_target( ${PROJECT_NAME}_cache_bathymetry.${config_relpath}
        DEPENDS ${run_dir}/done
      )
      add_dependencies( ${PROJECT_NAME}_cache_bathymetry ${PROJECT_NAME}_cache_bathymetry.${config_relpath} )

    endforeach()

endfunction()
