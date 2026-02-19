# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

ecbuild_info("[wam]")

# find_package(ParMETIS QUIET)   -->  ParMETIS was disabled in source code
# FIXME[IFS-GGG]: This is a workaround to link against static library libparmetis.a
# which was compiled *without* -fPIC option.

# if(IFS_SHARED_LIBS AND PARMETIS_FOUND)
#   set(PARMETIS_NON_PIC_LIBRARIES ${PARMETIS_LIBRARIES})
#   set(PARMETIS_LIBRARIES)
# endif()

set( WAM_DEFINITIONS ${IFS_DEFINITIONS} )
if( multio_FOUND )
  list(APPEND WAM_DEFINITIONS WAM_HAVE_MULTIO)
endif()

# FIXME: cyclic dependency between ifs & wam
list(APPEND wam_ifs_srcs  # Callback functions for handlers in coupled mode
  ifstowam.F90         # via yommp0, yomtag
  outwspec_io_serv.F90 # via wv_io_serv_suiosctmpl (interface still used by wam in yowcoup )
  outint_io_serv.F90   # via wv_io_serv_suiosctmpl (interface still used by wam in yowcoup )
)

list( APPEND wam_srcs
  wam/Wam_oper/abort1.F90
  wam/Wam_oper/adjust.F90
  wam/Wam_oper/airsea.F90
  wam/Wam_oper/aki.F90
  wam/Wam_oper/aki_ice.F90
  wam/Wam_oper/alphap_tail.F90
  wam/Wam_oper/bouinpt.F90
  wam/Wam_oper/buildstress.F90
  wam/Wam_oper/cal_second_order_spec.F90
  wam/Wam_oper/cdustarz0.F90
  wam/Wam_oper/check.F90
  wam/Wam_oper/checkcfl.F90
  wam/Wam_oper/chnkmin.F90
  wam/Wam_oper/cigetdeac.F90
  wam/Wam_oper/cimsstrn.F90
  wam/Wam_oper/cireduce.F90
  wam/Wam_oper/ciwabr.F90
  wam/Wam_oper/ciwaf.F90
  wam/Wam_oper/confile.F90
  wam/Wam_oper/ctuw.F90
  wam/Wam_oper/ctuwdrv.F90
  wam/Wam_oper/ctuwini.F90
  wam/Wam_oper/ctuwupdt.F90
  wam/Wam_oper/current2wam.F90
  wam/Wam_oper/difdate.F90
  wam/Wam_oper/dominant_period.F90
  wam/Wam_oper/expand_string.F90
  wam/Wam_oper/femean.F90
  wam/Wam_oper/femeanws.F90
  wam/Wam_oper/file_transfer.F90
  wam/Wam_oper/findb.F90
  wam/Wam_oper/fkmean.F90
  wam/Wam_oper/fldinter.F90
  wam/Wam_oper/fndprt.F90
  wam/Wam_oper/frcutindex.F90
  wam/function/gc_dispersion.h
  wam/Wam_oper/get_preset_wgrib_template.F90
  wam/Wam_oper/getcurr.F90
  wam/Wam_oper/getfrstwnd.F90
  wam/Wam_oper/getspec.F90
  wam/Wam_oper/getstress.F90
  wam/Wam_oper/getwnd.F90
  wam/Wam_oper/gradi.F90
  wam/Wam_oper/grib2wgrid.F90
  wam/Wam_oper/grstname.F90
  wam/Wam_oper/gsfile_new.F90
  wam/Wam_oper/h_max.F90
  wam/Wam_oper/halphap.F90
  wam/Wam_oper/headbc.F90
  wam/Wam_oper/imphftail.F90
  wam/Wam_oper/implsch.F90
  wam/Wam_oper/incdate.F90
  wam/Wam_oper/inisnonlin.F90
  wam/Wam_oper/init_fieldg.F90
  wam/Wam_oper/init_sdiss_ardh.F90
  wam/Wam_oper/init_x0tauhf.F90
  wam/Wam_oper/initdpthflds.F90
  wam/Wam_oper/initgc.F90
  wam/Wam_oper/initialint.F90
  wam/Wam_oper/initmdl.F90
  wam/Wam_oper/initnemocpl.F90
  wam/Wam_oper/iniwcst.F90
  wam/Wam_oper/intpol.F90
  wam/Wam_oper/intspec.F90
  wam/Wam_oper/inwgrib.F90
  wam/Wam_oper/iwam_get_unit.F90
  wam/Wam_oper/jafu.F90
  wam/Wam_oper/jonswap.F90
  wam/Wam_oper/kerkei.F90
  wam/Wam_oper/kgribsize.F90
  wam/Wam_oper/kurtosis.F90
  wam/Wam_oper/kzeone.F90
  wam/Wam_oper/makegrid.F90
  wam/Wam_oper/mblock.F90
  wam/Wam_oper/mbounc.F90
  wam/Wam_oper/mbounf.F90
  wam/Wam_oper/mboxb.F90
  wam/Wam_oper/mchunk.F90
  wam/Wam_oper/mcout.F90
  wam/Wam_oper/means.F90
  wam/Wam_oper/meansqs.F90
  wam/Wam_oper/meansqs_gc.F90
  wam/Wam_oper/meansqs_lf.F90
  wam/Wam_oper/mfredir.F90
  wam/Wam_oper/mgrid.F90
  wam/Wam_oper/micep.F90
  wam/Wam_oper/mintf.F90
  wam/Wam_oper/mnintw.F90
  wam/Wam_oper/mpabort.F90
  wam/Wam_oper/mpbcastgrid.F90
  wam/Wam_oper/mpbcastintfld.F90
  wam/Wam_oper/mpbcastsarin.F90
  wam/Wam_oper/mpclose_unit.F90
  wam/Wam_oper/mpcrtbl.F90
  wam/Wam_oper/mpdecomp.F90
  wam/Wam_oper/mpdistribfl.F90
  wam/Wam_oper/mpdistribscfld.F90
  wam/Wam_oper/mpexchng.F90
  wam/Wam_oper/mpexchngsarin.F90
  wam/Wam_oper/mpfldtoifs.F90
  wam/Wam_oper/mpgatherbc.F90
  wam/Wam_oper/mpgatherfl.F90
  wam/Wam_oper/mpgatherscfld.F90
  wam/Wam_oper/mpminmaxavg.F90
  wam/Wam_oper/mpuserin.F90
  wam/Wam_oper/mstart.F90
  wam/Wam_oper/mswell.F90
  wam/Wam_oper/mtabs.F90
  wam/Wam_oper/mubuf.F90
  wam/Wam_oper/mwp1.F90
  wam/Wam_oper/mwp2.F90
  wam/Wam_oper/newwind.F90
  wam/Wam_oper/nlweigt.F90
  wam/Wam_oper/notim.F90
  wam/Wam_oper/ns_gc.F90
  wam/Wam_oper/omegagc.F90
  wam/Wam_oper/out_onegrdpt.F90
  wam/Wam_oper/out_onegrdpt_sp.F90
  wam/Wam_oper/outbc.F90
  wam/Wam_oper/outbeta.F90
  wam/Wam_oper/outblock.F90
  wam/Wam_oper/outbs.F90
  wam/Wam_oper/outcom.F90
  wam/Wam_oper/outgrid.F90
  wam/Wam_oper/outint.F90
  wam/Wam_oper/outmdldcp.F90
  wam/Wam_oper/outnam.F90
  wam/Wam_oper/outpp.F90
  wam/module/output_struct.F90
  wam/Wam_oper/outsetwmask.F90
  wam/Wam_oper/outspec.F90
  wam/Wam_oper/outstep0.F90
  wam/Wam_oper/outwint.F90
  wam/Wam_oper/outwnorm.F90
  wam/Wam_oper/outwpsp.F90
  wam/Wam_oper/outwspec.F90
  wam/Wam_oper/packi.F90
  wam/Wam_oper/packr.F90
  wam/Wam_oper/parmean.F90
  wam/Wam_oper/peak.F90
  wam/Wam_oper/peak_ang.F90
  wam/Wam_oper/peak_freq.F90
  wam/Wam_oper/peakfri.F90
  wam/Wam_oper/preset_wgrib_template.F90
  wam/Wam_oper/prewind.F90
  wam/Wam_oper/proenvhalo.F90
  wam/Wam_oper/propag_wam.F90
  wam/Wam_oper/propags.F90
  wam/Wam_oper/propags1.F90
  wam/Wam_oper/propags2.F90
  wam/Wam_oper/propdot.F90
  wam/Wam_oper/readbou.F90
  wam/Wam_oper/readfl.F90
  wam/Wam_oper/readpre.F90
  wam/Wam_oper/readsta.F90
  wam/Wam_oper/readstress.F90
  wam/Wam_oper/readwgrib.F90
  wam/Wam_oper/readwind.F90
  wam/Wam_oper/recvnemofields.F90
  wam/Wam_oper/rotspec.F90
  wam/Wam_oper/runwam.F90
  wam/Wam_oper/savspec.F90
  wam/Wam_oper/savstress.F90
  wam/Wam_oper/sbottom.F90
  wam/Wam_oper/scosfl.F90
  wam/Wam_oper/sdepthlim.F90
  wam/Wam_oper/sdissip.F90
  wam/Wam_oper/sdissip_ard.F90
  wam/Wam_oper/sdissip_jan.F90
  wam/Wam_oper/sdiwbk.F90
  wam/Wam_oper/se10mean.F90
  wam/Wam_oper/sebtmean.F90
  wam/Wam_oper/second_order_lib.F90
  wam/Wam_oper/secondhh.F90
  wam/Wam_oper/secondhh_gen.F90
  wam/Wam_oper/secspom.F90
  wam/Wam_oper/semean.F90
  wam/Wam_oper/sep3tr.F90
  wam/Wam_oper/sepwisw.F90
  wam/Wam_oper/set_wflags.F90
  wam/Wam_oper/setice.F90
  wam/Wam_oper/setmarstype.F90
  wam/Wam_oper/setwavphys.F90
  wam/Wam_oper/sinflx.F90
  wam/Wam_oper/sinput.F90
  wam/Wam_oper/sinput_ard.F90
  wam/Wam_oper/sinput_jan.F90
  wam/Wam_oper/skewness.F90
  wam/Wam_oper/smoothsarspec.F90
  wam/Wam_oper/snonlin.F90
  wam/Wam_oper/spectra.F90
  wam/Wam_oper/spr.F90
  wam/Wam_oper/stat_nl.F90
  wam/Wam_oper/sthq.F90
  wam/Wam_oper/stokesdrift.F90
  wam/Wam_oper/stokestrn.F90
  wam/Wam_oper/stress_gc.F90
  wam/Wam_oper/stresso.F90
  wam/Wam_oper/strspec.F90
  wam/Wam_oper/tables_2nd.F90
  wam/Wam_oper/tabu_swellft.F90
  wam/Wam_oper/tau_phi_hf.F90
  wam/Wam_oper/taut_z0.F90
  wam/Wam_oper/topoar.F90
  wam/Wam_oper/transf.F90
  wam/Wam_oper/transf_bfi.F90
  wam/Wam_oper/transf_r.F90
  wam/Wam_oper/transf_snl.F90
  wam/Wam_oper/uibou.F90
  wam/Wam_oper/uiprep.F90
  wam/Wam_oper/unblkrord.F90
  wam/Wam_oper/unsetice.F90
  wam/module/unstruct_bound.F90
  wam/module/unstruct_curr.F90
  wam/module/unwam.F90
  wam/Wam_oper/updnemofields.F90
  wam/Wam_oper/updnemostress.F90
  wam/Wam_oper/userin.F90
  wam/Wam_oper/vmin.F90
  wam/Wam_oper/vmin_d.F90
  wam/Wam_oper/vplus.F90
  wam/Wam_oper/vplus_d.F90
  wam/Wam_oper/w_maxh.F90
  wam/Wam_oper/w_mode_st.F90
  wam/Wam_oper/wam_nproma.F90
  wam/Wam_oper/wam_sorti.F90
  wam/Wam_oper/wam_sortini.F90
  wam/Wam_oper/wam_u2l1cr.F90
  wam/Wam_oper/wam_user_clock.F90
  wam/Wam_oper/wamadswstar.F90
  wam/Wam_oper/wamcur.F90
  wam/Wam_oper/wamintgr.F90
  wam/Wam_oper/wamodel.F90
  wam/Wam_oper/wamwnd.F90
  wam/Wam_oper/wavemdl.F90
  wam/Wam_oper/wdfluxes.F90
  wam/Wam_oper/wdirspread.F90
  wam/Wam_oper/weflux.F90
  wam/Wam_oper/wgrib2fdb.F90
  wam/Wam_oper/wgribencode.F90
  wam/Wam_oper/wgribencode_model.F90
  wam/Wam_oper/wgribenout.F90
  wam/Wam_oper/wgribout.F90
  wam/Wam_oper/wnfluxes.F90
  wam/Wam_oper/wposnam.F90
  wam/Wam_oper/writefl.F90
  wam/Wam_oper/writestress.F90
  wam/Wam_oper/writsta.F90
  wam/Wam_oper/wsigstar.F90
  wam/Wam_oper/wsmfen.F90
  wam/Wam_oper/wstream_strg.F90
  wam/Wam_oper/wvalloc.F90
  wam/Wam_oper/wvdealloc.F90
  wam/Wam_oper/wvfricvelo.F90
  wam/Wam_oper/wvwamdecomp.F90
  wam/Wam_oper/wvwaminit.F90
  wam/Wam_oper/wvwaminit1.F90
  wam/module/yow_rank_gloloc.F90
  wam/module/yowaltas.F90
  wam/module/yowcard.F90
  wam/module/yowchecksmodule.F90
  wam/module/yowcinp.F90
  wam/module/yowcoer.F90
  wam/module/yowconst_2nd.F90
  wam/module/yowcoup.F90
  wam/module/yowcout.F90
  wam/module/yowcpbo.F90
  wam/module/yowcurg.F90
  wam/module/yowcurr.F90
  wam/module/yowdatapool.F90
  wam/module/yowdes.F90
  wam/module/yowdrvtype.F90
  wam/module/yowelementpool.F90
  wam/module/yowerror.F90
  wam/module/yowexchangeModule.F90
  wam/module/yowfpbo.F90
  wam/module/yowfred.F90
  wam/module/yowgrib_handles.F90
  wam/module/yowgribhd.F90
  wam/module/yowgrid.F90
  wam/module/yowice.F90
  wam/module/yowincludes.h
  wam/module/yowindn.F90
  wam/module/yowintp.F90
  wam/module/yowjons.F90
  wam/module/yowmap.F90
  wam/module/yowmean.F90
  wam/module/yowmespas.F90
  wam/module/yowmpiModule.F90
  wam/module/yowmpp.F90
  wam/module/yownemoflds.F90
  wam/module/yownodepool.F90
  wam/module/yowparam.F90
  wam/module/yowpcons.F90
  wam/module/yowpd.F90
  wam/module/yowpdlibmain.F90
  wam/module/yowphys.F90
  wam/module/yowprproc.F90
  wam/module/yowrankModule.F90
  wam/module/yowrefd.F90
  wam/module/yowsaras.F90
  wam/module/yowshal.F90
  wam/module/yowsidepool.F90
  wam/module/yowspec.F90
  wam/module/yowstat.F90
  wam/module/yowtabl.F90
  wam/module/yowtemp.F90
  wam/module/yowtest.F90
  wam/module/yowtext.F90
  wam/module/yowtrains.F90
  wam/module/yowubuf.F90
  wam/module/yowunit.F90
  wam/module/yowunpool.F90
  wam/module/yowwami.F90
  wam/module/yowwind.F90
  wam/module/yowwndg.F90
  wam/Wam_oper/z0wave.F90
  wam/Wam_oper/chkoops.F90
  wam/module/wam_multio_mod.F90
  wam/module/yowgrib.F90
  wam/module/yowabort.F90
  wam/module/yowassi.F90
  wam/module/yownemoio.F90
  wam/module/parkind_wave.F90

  wam/Alt/getclo.F # only for certain command-line tools
)
ecbuild_add_library(TARGET wam.${PREC}
  PRIVATE_DEFINITIONS ${WAM_DEFINITIONS}
  PUBLIC_INCLUDES
    #$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/wam/Alt/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/wam/function>
    $<BUILD_INTERFACE:${wam_intfb_includes}>

  PUBLIC_LIBS
    fiat
    parkind_${prec}
  PRIVATE_LIBS
    OpenMP::OpenMP_Fortran
    eccodes_f90
    wam_intfb
    ${MULTIO_LIBRARIES}
    ${NEMOVAR_LIBRARIES}
    # ${PARMETIS_LIBRARIES}
    SOURCES
      ${wam_srcs}
)

# ----------------------------------------------------------------------------
# Tools that don't need data assimilation
# ----------------------------------------------------------------------------

foreach(name IN ITEMS

    bouint
    preproc
    preset
    intwaminput
    write_mpdecomp
    create_wam_bathymetry
    create_wam_bathymetry_ETOPO1)

  if( NOT TARGET ${name} )
    ecbuild_add_executable(TARGET ${name}
      SOURCES wam/Wam_oper/${name}.F90
      LIBS
        wam.${PREC}
        wam_intfb
        MPI::MPI_Fortran
        OpenMP::OpenMP_Fortran
        # ${PARMETIS_NON_PIC_LIBRARIES}
      LINKER_LANGUAGE Fortran)
  endif()

endforeach()

# ----------------------------------------------------------------------------

# Standalone tools available with suffix "-standalone" that can be used
# in forecast-only or uncoupled mode, and have no further dependencies.
# These tools are compiled again, full-featured without suffix in "wamassitools.cmake"

foreach(name IN ITEMS
    chief)
  if( NOT TARGET ${name}-standalone )
    ecbuild_add_executable(TARGET ${name}-standalone
      SOURCES wam/Wam_oper/${name}.F90
      LIBS
        wam.${PREC}
        wam_intfb
        MPI::MPI_Fortran
        OpenMP::OpenMP_Fortran
        # ${PARMETIS_NON_PIC_LIBRARIES}
      LINKER_LANGUAGE Fortran)
  endif()
endforeach()

# ----------------------------------------------------------------------------
