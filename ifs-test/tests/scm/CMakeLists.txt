# (C) Copyright 1996-2019 ECMWF.

get_filename_component(SCM_EXECUTABLE_NAME ${SCM_EXECUTABLE} NAME)

if( TARGET ${SCM_EXECUTABLE_NAME} )
  list( APPEND TEST_DEPENDS ${SCM_EXECUTABLE_NAME} )
endif()

if( SCM_EXECUTABLE_NAME STREQUAL MASTER_SCM.DP )
  set(DOUBLE_PRECISION ON)
endif()

add_subdirectory( common )

#################################################################
# Helper function
function( add_ifs_test TARGET )

    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/${TARGET} )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/ifs-run ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs-run )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/ifs-debug ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs-debug )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/ifs-clean ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs-clean )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/ifs-check-bitid ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs-check-bitid )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/ifs-check-tracers ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs-check-tracers )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/ifs-grep-norms.pl ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs-grep-norms.pl )

  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/set_launcher.bash ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/set_launcher.bash )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${IFSTEST_BINARY_DIR}/ifs-run.py
      ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs-run.py
    )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/set_launcher.bash ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/set_launcher.bash )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${SHARE_BINARY_DIR}/ifs_env ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs_env )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink 
      ${STATIC_DATA_PATH} ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifsdata )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink 
      ${DATA_CACHE_DIR}/scm/scm_in.nc ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/scm_in.nc )

    if( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/params" )
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink 
        ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/params ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/params )
    endif()

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink 
      ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/setup ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/setup )

    if( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/postprocessing" )
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink 
        ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/postprocessing ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/postprocessing )
    endif()

    ecbuild_add_test( TARGET ifs_scm_${TARGET}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${TARGET} 
              COMMAND "./ifs-run"
              DEPENDS ${TEST_DEPENDS}
            )

    set_property(TEST ifs_scm_${TARGET} APPEND PROPERTY LABELS scm)

endfunction()

#################################################################
# add tests 

add_ifs_test( test_fc_TWPICE )


if( ENABLE_TESTS )
  set_property(TEST ifs_scm_test_fc_TWPICE                      APPEND PROPERTY LABELS forecast serial)


endif()
