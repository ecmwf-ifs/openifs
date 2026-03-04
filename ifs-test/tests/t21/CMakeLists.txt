# (C) Copyright 1996-2019 ECMWF.

get_filename_component(IFS_EXECUTABLE_NAME ${IFS_EXECUTABLE} NAME)

if( TARGET ${IFS_EXECUTABLE_NAME} )
  list( APPEND TEST_DEPENDS ${IFS_EXECUTABLE_NAME} )
endif()

if( IFS_EXECUTABLE_NAME STREQUAL ifsMASTER.DP )
  set(DOUBLE_PRECISION ON)
endif()

# Path to the t21 inidata directory on disk (fetched from ECPDS).
set( INIDATA_PATH "${DATA_CACHE_DIR}/t21/inidata")

# Path to the 4dvar data directory on disk (fetched from ECPDS).
set( 4DVAR_PATH "${DATA_CACHE_DIR}/t21/4dvardata")

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
      ${INIDATA_PATH}/ICMGGepc8INIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGepc8INIT )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${INIDATA_PATH}/ICMGGepc8INIUA ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGepc8INIUA )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${INIDATA_PATH}/ICMSHepc8INIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMSHepc8INIT )

    ### for test_adj
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${INIDATA_PATH}/ICMGGepc8INIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGepc8IMIN )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${INIDATA_PATH}/ICMGGepc8INIUA ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGepc8IMIUA )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${INIDATA_PATH}/ICMSHepc8INIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMSHepc8IMIN )

    ### for test_4dvar
    foreach(test_name IN ITEMS airep_t ecv_airep_t ecv_iasi ecv_orog ecv_bfas compo)
      if(${TARGET} MATCHES 4dvar_${test_name})
        #
        foreach(file ifs-check-cvg ifsobs_convert_ascii_CCMA.sh ifsobs_reformat_ascii.py mklinks)
          execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${PROJECT_SOURCE_DIR}/bin/${file} ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/${file} )
        endforeach()
        #
        foreach(file grib_filter create_ioassign convert_dbase)
          execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${CMAKE_BINARY_DIR}/bin/${file} ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/${file} )
        endforeach()
        #
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${DATA_CACHE_DIR}/sat_data/ ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/xdata )
        #
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${INIDATA_PATH}/ICMSHepc8INIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMRFepc80000 )
        #
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${4DVAR_PATH}/ifs ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs )
        #
        foreach(err_name IN ITEMS btmp gh lnsp o3 q r t u ucdv ucln uctp v vo)
          execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${4DVAR_PATH}/errgrib.${err_name} ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/errgrib.${err_name} )
        endforeach()
        #
        foreach(file single_obs.txt dummy_obs.txt)
          execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/${file} ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/${file} )
        endforeach()
        #
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${4DVAR_PATH}/wavelet_T21_L19.cv  ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/wavelet.cv )
      endif()
    endforeach()
    #
    foreach(test_name IN ITEMS ecv_airep_t ecv_iasi ecv_bfas)
      if(${TARGET} MATCHES 4dvar_${test_name})
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink 
          ${4DVAR_PATH}/errgrib.ecv3d001 ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/errgrib.ecv3d001 )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink 
          ${4DVAR_PATH}/errgrib.ecv2d001 ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/errgrib.ecv2d001 )
      endif()
    endforeach()
    #
    foreach(test_name IN ITEMS ecv_airep_t ecv_iasi ecv_orog ecv_bfas)
      if(${TARGET} MATCHES 4dvar_${test_name})
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/pre_ecv ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/pre_ecv )
      endif()
    endforeach()
    #
    if(${TARGET} MATCHES 4dvar_ecv_airep_t )
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${4DVAR_PATH}/acv_blh ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/acv_blh )
      foreach(pert RANGE 1 5)
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${4DVAR_PATH}/enspert_T21_00${pert} ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/enspert_T21_00${pert} )
      endforeach()
    endif()
    #
    if(${TARGET} MATCHES 4dvar_ecv_orog )
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${4DVAR_PATH}/errgrib.ecv3d001_sdfor ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/errgrib.ecv3d001 )
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${4DVAR_PATH}/errgrib.ecv2d001_sdfor ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/errgrib.ecv2d001 )
    endif()
    #
    if(${TARGET} MATCHES 4dvar_ecv_iasi )
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/add_o3.ksh ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/add_o3.ksh )
    endif()

    ### for test_compo_fc etc.
    if(${TARGET} MATCHES compo OR ${TARGET} MATCHES glomap OR ${TARGET} MATCHES ecv_bfas)
      if(${TARGET} MATCHES mozart)
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${INIDATA_PATH}/ICMGGhmerINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGhmerINIT )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${INIDATA_PATH}/ICMGGhmerINIUA ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGhmerINIUA )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${INIDATA_PATH}/ICMSHhmerINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMSHhmerINIT )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${INIDATA_PATH}/ICMCLhmerINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMCLhmerINIT )
      elseif(${TARGET} MATCHES mocage)
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${INIDATA_PATH}/ICMGGhmeqINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGhmeqINIT )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${INIDATA_PATH}/ICMGGhmeqINIUA ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGhmeqINIUA )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${INIDATA_PATH}/ICMSHhmeqINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMSHhmeqINIT )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${INIDATA_PATH}/ICMCLhmeqINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMCLhmeqINIT )
      else()
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${INIDATA_PATH}/ICMGGhmecINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGhmecINIT )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${INIDATA_PATH}/ICMGGhmecINIUA ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGhmecINIUA )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${INIDATA_PATH}/ICMSHhmecINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMSHhmecINIT )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${INIDATA_PATH}/ICMCLhmecINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMCLhmecINIT )
        if(${TARGET} MATCHES 4dvar )
          execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${INIDATA_PATH}/ICMSHhmecINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMRFhmec0000 )
          execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${INIDATA_PATH}/ICMSHhmecINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMSHhmecIMIN )
          execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${INIDATA_PATH}/ICMGGhmecINIUA ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGhmecIMIUA )
          execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${INIDATA_PATH}/ICMGGhmecINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGhmecIMIN )
        endif()
      endif()
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
    endif()

    ### for test_glomap_fc/test_glomap_edge_fc
    if(${TARGET} MATCHES glomap)
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${COMPO_DATA_PATH}/glomap_prep ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/glomap_prep )
    endif()

    ### for test_updclie/test_updclie2
    if(${TARGET} MATCHES updclie2)
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${INIDATA_PATH}/ICMCLgufnINIT.LMCCIECtrue ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMCLepc8INIT )
    else()
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${INIDATA_PATH}/ICMCLgu6bINIT.LMCCIECfalse ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMCLepc8INIT )
    endif()

    ### for test_compo_cb05bascoe_fc
    if(${TARGET} MATCHES bascoe)
      foreach(file info_j.dat j_ccl4.dat j_ch2br2.dat j_ch4.dat j_cloo.dat j_hcfc22.dat j_n2o5.dat j_o3-o1d.dat j_cclfo.dat j_ch2o-co.dat j_chbr3.dat j_co2.dat j_hcl.dat j_n2o.dat j_o3-o.dat j_br2.dat j_cf2o.dat j_ch2o-hco.dat j_cl2.dat j_h2o2-h2o.dat j_hno3.dat j_no2.dat j_oclo.dat j_brcl.dat j_cfc113.dat j_ch3br.dat j_cl2o2.dat j_h2o2-ho2.dat j_ho2.dat j_no3-o2.dat j_pan.dat j_bro.dat j_cfc114.dat j_ch3ccl3.dat j_clno2.dat j_h2o2-oh.dat j_ho2no2-ho2.dat j_no3-o.dat j_brono2-br.dat j_cfc115.dat j_ch3cl.dat j_clo.dat j_h2o.dat j_ho2no2-oh.dat j_no.dat j_brono2-bro.dat j_cfc11.dat j_ch3co3.dat j_clono2-cl.dat j_ha1211.dat j_hobr.dat j_o2-o1d.dat j_ccl2o.dat j_cfc12.dat j_ch3ooh.dat j_clono2-clo.dat j_ha1301.dat j_hocl.dat j_o2-o.dat CO2_vertprof.dat NO_vertprof.dat NO2_vertprof.dat O3_vertprof.dat solflux_socrates.dat solflux_socrates-lean97.dat srb_kockarts94.dat crs_o2-o.dat crs_o3-o.dat crs_o3-o1d.dat crs_ho2.dat crs_h2o2-oh.dat crs_no2.dat crs_no3-o.dat crs_no3-o2.dat crs_n2o5.dat crs_hno3.dat crs_ho2no2-ho2.dat crs_cl2.dat crs_br2.dat crs_oclo.dat crs_cl2o2.dat crs_hocl.dat crs_clono2-cl.dat crs_brcl.dat crs_brono2-br.dat crs_cclfo.dat crs_ccl4.dat crs_ccl2o.dat crs_brono2-bro.dat crs_bro.dat crs_cfc114.dat crs_cfc113.dat crs_cfc11.dat crs_cf2o.dat crs_ch2o-co.dat crs_ch2br2.dat crs_cfc12.dat crs_cfc115.dat crs_ch3co3.dat crs_ch3cl.dat crs_ch3ccl3.dat crs_ch3br.dat crs_ch2o-hco.dat crs_clno2.dat crs_chbr3.dat crs_ch4.dat crs_ch3ooh.dat crs_h2o.dat crs_co2.dat crs_cloo.dat crs_clono2-clo.dat crs_clo.dat crs_hcfc22.dat crs_ha1301.dat crs_ha1211.dat crs_h2so4.dat crs_h2o2-ho2.dat crs_no.dat crs_n2o.dat crs_hobr.dat crs_ho2no2-oh.dat crs_hcl.dat crs_so3.dat crs_pan.dat crs_ocs.dat crs_o2-o1d.dat)
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${COMPO_DATA_PATH}/bascoe_input/${file} ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/${file} )
      endforeach()
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${COMPO_DATA_PATH}/bascoe_input/lbc_CMIP6-1750-2100_15deg_CIFS_v1.dat ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/BASCOE_LBC.dat )
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${COMPO_DATA_PATH}/bascoe_input/ICBG_SAD_clim.txt ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICBG_SAD_clim.txt )
    endif()

    ### for test_compo_mocage_fc
    if(${TARGET} MATCHES mocage)
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${COMPO_DATA_PATH}/mocage_input/Jdata08_CIFS_v2.bin ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/mocage_jdata_v2 )
    endif()

    ### for test_compo_mozart_fc
    if(${TARGET} MATCHES mozart)
      foreach(file LBC_CH4_CO2.dat LBC_1850-2100_REF2_za_c080114.dat wasolar_1948-2140_daily_c100805.dat jacet_long.dat jcf2cl2_long.dat jch2o_a_long.dat jch4_b_long.dat jclono2_b_long.dat jh2o2_long.dat jhcl_long.dat jhocl_long.dat jn2o_long.dat jo3_b_long.dat jbr2_long.dat jcf2clbr_long.dat jch2o_b_long.dat jchbr2cl_long.dat jclono_long.dat jh2o_a_long.dat jhf_long.dat jhono_long.dat jno2_long.dat joclo_long.dat jbrcl_long.dat jcf3br_long.dat jch3br_long.dat jchbr3_long.dat jco2_long.dat jh2o_b_long.dat jhno3_long.dat jhyac_long.dat jno3_a_long.dat jocs_long.dat jbrno2_long.dat jcfc113_long.dat jch3ccl3_long.dat jchbrcl2_long.dat jcof2_long.dat jh2o_c_long.dat jhno3nat_long.dat jmacr_a_long.dat jno3_b_long.dat jpan_long.dat jbro_long.dat jcfc114_long.dat jch3cho_long.dat jcl2_long.dat jcofcl_long.dat jh2so4_long.dat jhno3sts_long.dat jmacr_b_long.dat jno_i_long.dat jsf6_long.dat jbrono2_a_long.dat jcfc115_long.dat jch3cl_long.dat jcl2o2_long.dat jcs2_long.dat jhbr_long.dat jho2_long.dat jmgly_long.dat jno_long.dat jso2_long.dat jbrono2_b_long.dat jcfcl3_long.dat jch3o2no2_long.dat jclno2_long.dat jglyald_long.dat jhcfc141b_long.dat jho2no2_a_long.dat jmvk_long.dat jo2_a_long.dat jso3_long.dat jbrono_long.dat jch2br2_long.dat jch3ooh_long.dat jclo_long.dat jh1202_long.dat jhcfc142b_long.dat jho2no2_b_long.dat jn2o5_a_long.dat jo2_b_long.dat jso_long.dat jccl4_long.dat jch2brcl_long.dat jch4_a_long.dat jclono2_a_long.dat jh2402_long.dat jhcfc22_long.dat jhobr_long.dat jn2o5_b_long.dat jo3_a_long.dat)
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${COMPO_DATA_PATH}/mozart_input/${file} ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/${file} )
      endforeach()
    endif()

    ### for test_tl_taylor
    if(${TARGET} MATCHES tl_taylor)
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/ifs-test_tl.py ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs-test_tl.py )
    endif()

    ### for test compo_fc_dualconf
    if(${TARGET} MATCHES compo_fc_dualconf)
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/namelist_geom_part ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/namelist_geom_part )
    endif()

    ### for test_compo_updclie
    if(${TARGET} MATCHES compo_updclie)
    #execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
    #  ${INIDATA_PATH}/ICMGGhilqINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGhilqINIT )
    #execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
    #  ${INIDATA_PATH}/ICMGGhilqINIUA ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGhilqINIUA )
    #execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
    #  ${INIDATA_PATH}/ICMSHhilqINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMSHhilqINIT )
    #execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
    #  ${INIDATA_PATH}/ICMCLhilqINIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMCLhilqINIT )
    #execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
    #  ${INIDATA_PATH}/ICMCLhilqINIT_COMPO ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMCLhilqINIT_COMPO )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${INIDATA_PATH}/ICMCLhmecINIT_COMPO ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMCLhmecINIT_COMPO )
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

    ecbuild_add_test( TARGET ifs_t21_${TARGET}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}
              COMMAND "./ifs-run"
              DEPENDS ${TEST_DEPENDS}
            )

    set_property(TEST ifs_t21_${TARGET} APPEND PROPERTY LABELS t21)

endfunction()

#################################################################
# add tests

add_ifs_test( test_fc )
add_ifs_test( test_np2 )
add_ifs_test( test_nt2 )
add_ifs_test( test_sppt )
add_ifs_test( test_spp )
add_ifs_test( test_hybrid )
add_ifs_test( test_fc_tegen )
add_ifs_test( test_updclie )
add_ifs_test( test_updclie2 )
add_ifs_test( test_compo_fc )
# add_ifs_test( test_compo_fc_tegen ) # This one is not currently expected to work
add_ifs_test( test_compo_fc_climrad )
add_ifs_test( test_compo_edge_fc )
add_ifs_test( test_compo_aersuonly_fc )
add_ifs_test( test_compo_chemonly_fc )
add_ifs_test( test_compo_ghgonly_fc )
add_ifs_test( test_compo_cb05bascoe_fc )
add_ifs_test( test_compo_chemonly_RnPb_fc )
add_ifs_test( test_compo_chemonly_nwpo3_fc )
add_ifs_test( test_compo_updclie )


# Tests which will not succeed in OpenIFS-only builds.
if( NOT ifs_HAVE_OPENIFS_ONLY )

  add_ifs_test( test_fc_nh )
  add_ifs_test( test_compo_mocage_fc )
  add_ifs_test( test_compo_mozart_fc )
  add_ifs_test( test_glomap_fc )
  add_ifs_test( test_glomap_edge_fc )

  # Tests which will not succeed in forecast-only builds:
  if( NOT ifs_HAVE_FORECAST_ONLY )
    # Tests which will only succeed at double precision:
    if( DOUBLE_PRECISION )
      add_ifs_test( test_adj )
      add_ifs_test( test_tl_taylor )
      add_ifs_test( test_sv )
      add_ifs_test( test_4dvar_airep_t)
      add_ifs_test( test_4dvar_ecv_airep_t)
      add_ifs_test( test_4dvar_ecv_iasi)
      add_ifs_test( test_4dvar_ecv_orog)
      add_ifs_test( test_4dvar_ecv_bfas)
      add_ifs_test( test_4dvar_compo)
    endif()
  endif()

  # Test which can only be run if OOPS-IFS has been compiled:
  if( ifs_HAVE_OOPS )
    add_ifs_test( test_compo_fc_dualconf )
  endif()

endif()

if( ifs_HAVE_MGRIDS )
  if( DOUBLE_PRECISION )
    add_ifs_test( test_mgrids_advection_sladv )
  endif()
  add_ifs_test( test_mgrids_advection_mpdata )
endif()

if( ENABLE_TESTS )
  set_property(TEST ifs_t21_test_fc                      APPEND PROPERTY LABELS forecast serial)
  set_property(TEST ifs_t21_test_np2                     APPEND PROPERTY LABELS forecast mpi)
  set_property(TEST ifs_t21_test_nt2                     APPEND PROPERTY LABELS forecast threading)
  set_property(TEST ifs_t21_test_sppt                    APPEND PROPERTY LABELS forecast mpi threading)
  set_property(TEST ifs_t21_test_spp                     APPEND PROPERTY LABELS forecast mpi threading)
  set_property(TEST ifs_t21_test_hybrid                  APPEND PROPERTY LABELS forecast mpi threading)
  set_property(TEST ifs_t21_test_fc_tegen                APPEND PROPERTY LABELS forecast mpi threading)
  set_property(TEST ifs_t21_test_updclie                 APPEND PROPERTY LABELS forecast mpi threading climate)
  set_property(TEST ifs_t21_test_updclie2                APPEND PROPERTY LABELS forecast mpi threading climate)
  set_property(TEST ifs_t21_test_compo_fc                APPEND PROPERTY LABELS forecast mpi threading composition)
  # set_property(TEST ifs_t21_test_compo_fc_tegen          APPEND PROPERTY LABELS forecast mpi threading composition)
  set_property(TEST ifs_t21_test_compo_fc_climrad        APPEND PROPERTY LABELS forecast mpi threading composition)
  set_property(TEST ifs_t21_test_compo_edge_fc           APPEND PROPERTY LABELS forecast mpi threading composition)
  set_property(TEST ifs_t21_test_compo_aersuonly_fc      APPEND PROPERTY LABELS forecast mpi threading composition)
  set_property(TEST ifs_t21_test_compo_chemonly_fc       APPEND PROPERTY LABELS forecast mpi threading composition)
  set_property(TEST ifs_t21_test_compo_ghgonly_fc        APPEND PROPERTY LABELS forecast mpi threading composition)
  set_property(TEST ifs_t21_test_compo_cb05bascoe_fc     APPEND PROPERTY LABELS forecast mpi threading composition cb05bascoe)
  set_property(TEST ifs_t21_test_compo_chemonly_RnPb_fc  APPEND PROPERTY LABELS forecast mpi threading composition RnPb)
  set_property(TEST ifs_t21_test_compo_chemonly_nwpo3_fc APPEND PROPERTY LABELS forecast mpi threading composition nwpo3)
  set_property(TEST ifs_t21_test_compo_updclie           APPEND PROPERTY LABELS forecast mpi threading composition climate mcc_compo)

  # Tests to include if this is not an OpenIFS-only build:
  if( NOT ifs_HAVE_OPENIFS_ONLY )
    set_property(TEST ifs_t21_test_fc_nh                   APPEND PROPERTY LABELS forecast serial)
    set_property(TEST ifs_t21_test_compo_mocage_fc         APPEND PROPERTY LABELS forecast mpi threading composition mocage)
    set_property(TEST ifs_t21_test_compo_mozart_fc         APPEND PROPERTY LABELS forecast mpi threading composition mozart)
    set_property(TEST ifs_t21_test_glomap_fc               APPEND PROPERTY LABELS forecast mpi threading composition glomap)
    set_property(TEST ifs_t21_test_glomap_edge_fc          APPEND PROPERTY LABELS forecast mpi threading composition glomap)

    # Tests to include if this is not a forecast-only build:
    if( NOT ifs_HAVE_FORECAST_ONLY )
      # Tests which will only succeed at double precision:
      if( DOUBLE_PRECISION )
        set_property(TEST ifs_t21_test_adj                   APPEND PROPERTY LABELS adjoint  mpi threading)
        set_property(TEST ifs_t21_test_tl_taylor             APPEND PROPERTY LABELS tangent  mpi threading)
        set_property(TEST ifs_t21_test_sv                    APPEND PROPERTY LABELS adjoint  mpi threading)
        set_property(TEST ifs_t21_test_4dvar_airep_t         APPEND PROPERTY LABELS 4dvar    mpi threading)
        set_property(TEST ifs_t21_test_4dvar_ecv_airep_t     APPEND PROPERTY LABELS 4dvar    mpi threading ecv)
        set_property(TEST ifs_t21_test_4dvar_ecv_iasi        APPEND PROPERTY LABELS 4dvar    mpi threading ecv)
        set_property(TEST ifs_t21_test_4dvar_ecv_orog        APPEND PROPERTY LABELS 4dvar    mpi threading ecv)
        set_property(TEST ifs_t21_test_4dvar_ecv_bfas        APPEND PROPERTY LABELS 4dvar    mpi threading ecv)
        set_property(TEST ifs_t21_test_4dvar_compo           APPEND PROPERTY LABELS 4dvar    mpi threading composition)
      endif()
    endif()

    if( ifs_HAVE_OOPS )
      set_property(TEST ifs_t21_test_compo_fc_dualconf     APPEND PROPERTY LABELS forecast mpi threading composition dualconf)
    endif()

  endif()

  if( ifs_HAVE_MGRIDS )
    if( DOUBLE_PRECISION )
      set_property(TEST ifs_t21_test_mgrids_advection_sladv  APPEND PROPERTY LABELS forecast serial)
    endif()
    set_property(TEST ifs_t21_test_mgrids_advection_mpdata APPEND PROPERTY LABELS forecast serial)
  endif()

endif()
