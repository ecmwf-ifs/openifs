# (C) Copyright 1996-2019 ECMWF.

get_filename_component(IFS_EXECUTABLE_NAME ${IFS_EXECUTABLE} NAME)

if( TARGET ${IFS_EXECUTABLE_NAME} )
  list( APPEND TEST_DEPENDS ${IFS_EXECUTABLE_NAME} )
endif()

if( IFS_EXECUTABLE_NAME STREQUAL ifsMASTER.DP )
  set(DOUBLE_PRECISION ON)
endif()

# Path to the t42 inidata directory on disk (fetched from ECPDS).
set( INIDATA_PATH "${DATA_CACHE_DIR}/t42/inidata")

# Path to the t21 inidata directory on disk (fetched from ECPDS).
set( INIDATA_T21_PATH "${DATA_CACHE_DIR}/t21/inidata")

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
      ${SHARE_BINARY_DIR}/ifs_env ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs_env )

    ### for test_fc
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${INIDATA_PATH}/ICMGGhafbINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGhafbINIT )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${INIDATA_PATH}/ICMGGhafbINIUA ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGhafbINIUA )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${INIDATA_PATH}/ICMSHhafbINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMSHhafbINIT )



    ### for test_compo_fc
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${INIDATA_T21_PATH}/ICMGGhmecINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGhmecINIT )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${INIDATA_T21_PATH}/ICMGGhmecINIUA ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGhmecINIUA )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${INIDATA_T21_PATH}/ICMSHhmecINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMSHhmecINIT )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${INIDATA_T21_PATH}/ICMCLhmecINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMCLhmecINIT )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${COMPO_DATA_PATH}/OMI.data.extraterrest ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/OMI.data.extraterrest )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${COMPO_DATA_PATH}/aerosol_reduce.dat ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/aerosol_reduce.dat )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${COMPO_DATA_PATH}/tropo_look_up_cbmhybrid.dat ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/tropo_look_up_cbmhybrid.dat )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${COMPO_DATA_PATH}/uars_ratio.txt ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/uars_ratio.txt )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${COMPO_DATA_PATH}/haloe_ch4clim.dat ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/haloe_ch4clim.dat )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${COMPO_DATA_PATH}/21_full ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/21_full )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${COMPO_DATA_PATH}/42_full ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/42_full )


    ### for test compo_fc_dualconf
    if(${TARGET} MATCHES compo_fc_dualconf)
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/namelist_geom_part ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/namelist_geom_part )
    endif()


    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${STATIC_DATA_PATH} ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifsdata )

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

    ecbuild_add_test( TARGET ifs_t42_${TARGET}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}
              COMMAND "./ifs-run"
              DEPENDS ${TEST_DEPENDS}
            )

    set_property(TEST ifs_t42_${TARGET} APPEND PROPERTY LABELS t42)

endfunction()

#################################################################
# add tests

# Test which can only be run if OOPS-IFS has been compiled:
if( ifs_HAVE_OOPS )
  add_ifs_test( test_compo_fc_dualconf )
endif()


if( ENABLE_TESTS )

  if( ifs_HAVE_OOPS )
    set_property(TEST ifs_t42_test_compo_fc_dualconf APPEND PROPERTY LABELS forecast mpi threading composition dualconf)
  endif()

endif()
