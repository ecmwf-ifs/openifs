# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set_property( GLOBAL PROPERTY ecwam_test_configs )
function( register_ecwam_test_config )
  get_property( ecwam_test_configs GLOBAL PROPERTY ecwam_test_configs )
  foreach( arg ${ARGV} )
    list( APPEND ecwam_test_configs "${arg}")
  endforeach()
  list(REMOVE_DUPLICATES ecwam_test_configs)
  set_property( GLOBAL PROPERTY ecwam_test_configs "${ecwam_test_configs}")
endfunction()

function( ecwam_add_test test_name )

    set( one_value_keywords MPI OMP CONFIG )
    set( multi_value_keywords ENVIRONMENT )

    include(CMakeParseArguments)
    cmake_parse_arguments( PARSE_ARGV 1 ecwam_add_test "${options}" "${one_value_keywords}" "${multi_value_keywords}" )

    set( test_name_fixture ${test_name} )
    set( test_name_preproc ${test_name}_preproc )
    set( test_name_preset  ${test_name}_preset  )

    if( ecwam_add_test_MPI )
      set( test_name "${test_name}_mpi${ecwam_add_test_MPI}" )
    endif()
    if( NOT ecwam_add_test_OMP )
      set( ecwam_add_test_OMP 1 )
    else()
      set( test_name "${test_name}_omp${ecwam_add_test_OMP}" )
    endif()

    if( ecwam_add_test_OMP GREATER 1 AND NOT HAVE_OMP )
      ecbuild_info("Skipping test ${test_name} because OMP is disabled")
      return()
    endif()

    unset( LAUNCH )
    if( ecwam_add_test_MPI )
      if( MPIEXEC_EXECUTABLE )
        set( MPIEXEC ${MPIEXEC_EXECUTABLE} )
      endif()

      if( NOT MPIEXEC )
        find_program( MPIEXEC NAMES mpiexec mpirun lamexec srun
                      DOC "Executable for running MPI programs." )
      endif()

      if( MPIEXEC )
        if ( NOT MPIEXEC_NUMPROC_FLAG )
            set( MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "Flag used by MPI to specify the number of processes for MPIEXEC" )
        endif()

        # MPI_ARGS is left for users to define @ configure time e.g. -DMPI_ARGS="--oversubscribe"
        if( MPI_ARGS )
          set( LAUNCH "${MPIEXEC} ${MPI_ARGS} ${MPIEXEC_NUMPROC_FLAG} ${ecwam_add_test_MPI}" )
        else()
          set( LAUNCH "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${ecwam_add_test_MPI}" )
        endif()
        if( MPIEXEC_NUMTHREAD_FLAG )
          set( LAUNCH "${LAUNCH} ${MPIEXEC_NUMTHREAD_FLAG} ${ecwam_add_test_OMP}" )
        endif()
      endif()
    elseif( CMAKE_CROSSCOMPILING_EMULATOR )
      set( LAUNCH ${CMAKE_CROSSCOMPILING_EMULATOR} )
    endif()

    set( test_preproc_binary_dir ${CMAKE_CURRENT_BINARY_DIR}/${test_name_preproc} )
    set( run_dir_preproc ${test_preproc_binary_dir} )
    if( ECWAM_WORKSPACE )
      set( run_dir_preproc ${ECWAM_WORKSPACE}/tests/${test_name_preproc} )
    endif()

    if( NOT TEST ${test_name_preproc} )
      if( ECWAM_WORKSPACE )
        execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
          "${run_dir_preproc}"
          "${test_preproc_binary_dir}"
        )
      endif()
      file( MAKE_DIRECTORY  ${run_dir_preproc} )
      execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
          "${CMAKE_CURRENT_SOURCE_DIR}/${ecwam_add_test_CONFIG}"
          "${run_dir_preproc}/config.yml"
      )

      file( WRITE  ${run_dir_preproc}/run.sh "#!/usr/bin/env bash\n" )
      file( APPEND ${run_dir_preproc}/run.sh "export PATH=${CMAKE_BINARY_DIR}/bin:\$PATH\n" )
      file( APPEND ${run_dir_preproc}/run.sh "if $(${ECWAM_PROJECT_NAME}-run-preproc --completed); then\n" )
      file( APPEND ${run_dir_preproc}/run.sh "  echo \"${test_name_preproc} already completed, SKIPPING...\"; echo \n" )
      file( APPEND ${run_dir_preproc}/run.sh "  echo \"    To force rerun, execute\"; echo \n" )
      file( APPEND ${run_dir_preproc}/run.sh "  echo \"        ${run_dir_preproc}/clean.sh \"; echo \n" )
      file( APPEND ${run_dir_preproc}/run.sh "  exit 255\n" )
      file( APPEND ${run_dir_preproc}/run.sh "fi\n\n" )
      if( CMAKE_CROSSCOMPILING_EMULATOR )
      file( APPEND ${run_dir_preproc}/run.sh "export LAUNCH=\"\${LAUNCH:-${CMAKE_CROSSCOMPILING_EMULATOR}}\"\n" )
      endif()
      file( APPEND ${run_dir_preproc}/run.sh "export OMP_NUM_THREADS=1\n" )
      file( APPEND ${run_dir_preproc}/run.sh "${ECWAM_PROJECT_NAME}-run-preproc --run-dir=${run_dir_preproc} \"\${@}\" \n" )

      file( WRITE  ${run_dir_preproc}/clean.sh "#!/usr/bin/env bash\n" )
      file( APPEND ${run_dir_preproc}/clean.sh "find ${run_dir_preproc} -maxdepth 1 -mindepth 1 \\\n" )
      file( APPEND ${run_dir_preproc}/clean.sh "     -not -name 'run.sh' \\\n" )
      file( APPEND ${run_dir_preproc}/clean.sh "  -a -not -name 'clean.sh' \\\n" )
      file( APPEND ${run_dir_preproc}/clean.sh "  -a -not -name 'config.yml' \\\n" )
      file( APPEND ${run_dir_preproc}/clean.sh "  -print0 | xargs -0 rm -rf -- \n" )
  
      file( CHMOD  ${run_dir_preproc}/run.sh ${run_dir_preproc}/clean.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )

      add_test(NAME ${test_name_preproc}
        WORKING_DIRECTORY ${run_dir_preproc}
        COMMAND ${run_dir_preproc}/run.sh
      )
      set_tests_properties( ${test_name_preproc} PROPERTIES FIXTURES_SETUP ${test_name_preset} )
      set_tests_properties( ${test_name_preproc} PROPERTIES LABELS "setup" )
      set_tests_properties( ${test_name_preproc} PROPERTIES SKIP_RETURN_CODE 255 )

    endif()

    set( test_preset_binary_dir ${CMAKE_CURRENT_BINARY_DIR}/${test_name_preset} )
    set( run_dir_preset ${test_preset_binary_dir} )
    if( ECWAM_WORKSPACE )
      set( run_dir_preset ${ECWAM_WORKSPACE}/tests/${test_name_preset} )
    endif()

    if( NOT TEST ${test_name_preset} )
      if( ECWAM_WORKSPACE )
        execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
          "${run_dir_preset}"
          "${test_preset_binary_dir}"
        )
      endif()
      file( MAKE_DIRECTORY  ${run_dir_preset} )
      execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
          "${CMAKE_CURRENT_SOURCE_DIR}/${ecwam_add_test_CONFIG}"
          "${run_dir_preset}/config.yml"
      )

      file( WRITE  ${run_dir_preset}/run.sh "#!/usr/bin/env bash\n" )
      file( APPEND ${run_dir_preset}/run.sh "export PATH=${CMAKE_BINARY_DIR}/bin:\$PATH\n" )
      file( APPEND ${run_dir_preset}/run.sh "if $(${ECWAM_PROJECT_NAME}-run-preset --completed); then\n" )
      file( APPEND ${run_dir_preset}/run.sh "  echo \"${test_name_preset} already completed, SKIPPING...\"; echo \n" )
      file( APPEND ${run_dir_preset}/run.sh "  echo \"    To force rerun, execute\"; echo \n" )
      file( APPEND ${run_dir_preset}/run.sh "  echo \"        ${run_dir_preset}/clean.sh \"; echo \n" )
      file( APPEND ${run_dir_preset}/run.sh "  exit 255\n" )
      file( APPEND ${run_dir_preset}/run.sh "fi\n\n" )
      if( CMAKE_CROSSCOMPILING_EMULATOR )
      file( APPEND ${run_dir_preset}/run.sh "export LAUNCH=\"\${LAUNCH:-${CMAKE_CROSSCOMPILING_EMULATOR}}\"\n" )
      endif()
      file( APPEND ${run_dir_preset}/run.sh "export OMP_NUM_THREADS=1\n" )
      file( APPEND ${run_dir_preset}/run.sh "export ECWAM_PREPROC_RUN_DIR=${run_dir_preproc}\n")
      file( APPEND ${run_dir_preset}/run.sh "${ECWAM_PROJECT_NAME}-run-preset --run-dir=${run_dir_preset} \"\${@}\" \n" )

      file( WRITE  ${run_dir_preset}/clean.sh "#!/usr/bin/env bash\n" )
      file( APPEND ${run_dir_preset}/clean.sh "find ${run_dir_preset} -maxdepth 1 -mindepth 1 \\\n" )
      file( APPEND ${run_dir_preset}/clean.sh "     -not -name 'run.sh' \\\n" )
      file( APPEND ${run_dir_preset}/clean.sh "  -a -not -name 'clean.sh' \\\n" )
      file( APPEND ${run_dir_preset}/clean.sh "  -a -not -name 'config.yml' \\\n" )
      file( APPEND ${run_dir_preset}/clean.sh "  -print0 | xargs -0 rm -rf -- \n" )

      file( CHMOD  ${run_dir_preset}/run.sh ${run_dir_preset}/clean.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )

      add_test(NAME ${test_name_preset}
        WORKING_DIRECTORY ${run_dir_preset}
        COMMAND ${run_dir_preset}/run.sh
      )
      set_tests_properties( ${test_name_preset} PROPERTIES FIXTURES_REQUIRED ${test_name_preset} )
      set_tests_properties( ${test_name_preset} PROPERTIES FIXTURES_SETUP ${test_name_fixture} )
      set_tests_properties( ${test_name_preset} PROPERTIES LABELS "setup" )
      set_tests_properties( ${test_name_preset} PROPERTIES SKIP_RETURN_CODE 255 )
    endif()

    set( test_binary_dir ${CMAKE_CURRENT_BINARY_DIR}/${test_name} )
    set( run_dir ${test_binary_dir} )
    if( ECWAM_WORKSPACE )
      set( run_dir ${ECWAM_WORKSPACE}/tests/${test_name} )
      execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
        "${run_dir}"
        "${test_binary_dir}"
      )
    endif()
    file( MAKE_DIRECTORY  ${run_dir} )
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
      "${CMAKE_CURRENT_SOURCE_DIR}/${ecwam_add_test_CONFIG}"
      "${run_dir}/config.yml"
    )

    file( WRITE  ${run_dir}/run.sh "#!/usr/bin/env bash\n" )
    file( APPEND ${run_dir}/run.sh "export PATH=${CMAKE_BINARY_DIR}/bin:\$PATH\n" )
    if( CMAKE_CROSSCOMPILING_EMULATOR )
    file( APPEND ${run_dir}/run.sh "export LAUNCH=\"\${LAUNCH:-${LAUNCH}}\"\n" )
    endif()
    file( APPEND ${run_dir}/run.sh "export OMP_NUM_THREADS=${ecwam_add_test_OMP}\n" )
    file( APPEND ${run_dir}/run.sh "export ECWAM_PREPROC_RUN_DIR=${run_dir_preproc}\n")
    file( APPEND ${run_dir}/run.sh "export ECWAM_PRESET_RUN_DIR=${run_dir_preset}\n")
    file( APPEND ${run_dir}/run.sh "${ECWAM_PROJECT_NAME}-run-model --run-dir=${run_dir} \"\${@}\" \n" )

    file( WRITE  ${run_dir}/clean.sh "#!/usr/bin/env bash\n" )
    file( APPEND ${run_dir}/clean.sh "find ${run_dir} -maxdepth 1 -mindepth 1 \\\n" )
    file( APPEND ${run_dir}/clean.sh "     -not -name 'run.sh' \\\n" )
    file( APPEND ${run_dir}/clean.sh "  -a -not -name 'clean.sh' \\\n" )
    file( APPEND ${run_dir}/clean.sh "  -a -not -name 'config.yml' \\\n" )
    file( APPEND ${run_dir}/clean.sh "  -print0 | xargs -0 rm -rf -- \n" )

    file( CHMOD  ${run_dir}/run.sh ${run_dir}/clean.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )

    add_test(NAME ${test_name}
      WORKING_DIRECTORY ${run_dir}
      COMMAND ${run_dir}/run.sh
    )
    set_tests_properties( ${test_name} PROPERTIES FIXTURES_REQUIRED ${test_name_fixture} )

    if( ecwam_add_test_MPI )
      set_tests_properties( ${test_name} PROPERTIES LABELS "mpi" )
    endif()

    # Set test environment
    file( WRITE  ${run_dir}/ctest_env.sh "export LAUNCH_PARALLEL=\"${LAUNCH}\"\n" )
    file( APPEND ${run_dir}/ctest_env.sh "export LAUNCH_SERIAL=\"${CMAKE_CROSSCOMPILING_EMULATOR}\"\n" )
    file( APPEND ${run_dir}/ctest_env.sh "export OMP_NUM_THREADS=${ecwam_add_test_OMP}\n" )
    file( APPEND ${run_dir}/ctest_env.sh "export ECWAM_PREPROC_RUN_DIR=${run_dir_preproc}\n")
    file( APPEND ${run_dir}/ctest_env.sh "export ECWAM_PRESET_RUN_DIR=${run_dir_preset}\n")
    file( APPEND ${run_dir}/ctest_env.sh "export DR_HOOK=1\n")
    file( APPEND ${run_dir}/ctest_env.sh "export DR_HOOK_OPT=PROF\n")
    file( APPEND ${run_dir}/ctest_env.sh "export GSTATS=1\n")
    file( APPEND ${run_dir}/ctest_env.sh "export PATH=${CMAKE_BINARY_DIR}/share/ecwam/scripts:\$PATH\n")


    set( ecwam_environment OMP_NUM_THREADS=${ecwam_add_test_OMP} )
    list( APPEND ecwam_environment ECWAM_ENVIRONMENT=${run_dir}/ctest_env.sh )
    set_tests_properties( ${test_name} PROPERTIES ENVIRONMENT "${ecwam_environment}" )

    register_ecwam_test_config( "${CMAKE_CURRENT_SOURCE_DIR}/${ecwam_add_test_CONFIG}" )

endfunction()
