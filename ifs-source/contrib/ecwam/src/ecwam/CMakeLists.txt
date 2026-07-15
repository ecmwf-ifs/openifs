# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set( ecwam ${PROJECT_NAME} )

list( APPEND unwam_srcs
  unwam/output_struct.F90
  unwam/unstruct_bound.F90
  unwam/unstruct_curr.F90
  unwam/unwam.F90
  unwam/yow_rank_gloloc.F90
  unwam/yowchecksmodule.F90
  unwam/yowdatapool.F90
  unwam/yowelementpool.F90
  unwam/yowerror.F90
  unwam/yowexchangeModule.F90
  unwam/yowincludes.h
  unwam/yowmpiModule.F90
  unwam/yownodepool.F90
  unwam/yowpd.F90
  unwam/yowpdlibmain.F90
  unwam/yowrankModule.F90
  unwam/yowsidepool.F90
  unwam/yowunblkrord.F90
  unwam/yowunpool.F90
)

ecbuild_generate_fortran_interfaces(
    TARGET ${ecwam}_intfb
    DIRECTORIES ecwam
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..
    DESTINATION intfb
    INCLUDE_DIRS ecwam_intfb_includes
    PARALLEL ${FCM_PARALLEL}
)
install( DIRECTORY ${ecwam_intfb_includes}/
         DESTINATION include/${ecwam}
         COMPONENT interfaces )

list( APPEND ecwam_srcs
  abort1.F90
  adjust.F90
  airsea.F90
  aki.F90
  aki_ice.F90
  alphap_tail.F90
  bouinpt.F90
  buildstress.F90
  cal_second_order_spec.F90
  cdustarz0.F90
  check.F90
  checkcfl.F90
  chkoops.F90
  chnkmin.F90
  cigetdeac.F90
  cimsstrn.F90
  cireduce.F90
  ciwabr.F90
  ciwaf.F90
  ctcor.F90
  ctuw.F90
  ctuwdrv.F90
  ctuwini.F90
  ctuwupdt.F90
  current2wam.F90
  depthprpt.F90
  difdate.F90
  dominant_period.F90
  expand_string.F90
  femean.F90
  femeanws.F90
  file_transfer.F90
  findb.F90
  fkmean.F90
  fldinter.F90
  fndprt.F90
  frcutindex.F90
  gc_dispersion.h
  get_preset_wgrib_template.F90
  getcurr.F90
  getfrstwnd.F90
  getspec.F90
  getstress.F90
  getwnd.F90
  gradi.F90
  grib2wgrid.F90
  grstname.F90
  gsfile_new.F90
  h_max.F90
  halphap.F90
  headbc.F90
  imphftail.F90
  implsch.F90
  incdate.F90
  inisnonlin.F90
  init_fieldg.F90
  init_sdiss_ardh.F90
  init_x0tauhf.F90
  initdpthflds.F90
  initgc.F90
  initialint.F90
  initmdl.F90
  initnemocpl.F90
  iniwcst.F90
  intpol.F90
  intrpolchk.F90
  intspec.F90
  inwgrib.F90
  iwam_get_unit.F90
  jafu.F90
  jonswap.F90
  kerkei.F90
  kgribsize.F90
  kurtosis.F90
  kzeone.F90
  makegrid.F90
  mblock.F90
  mbounc.F90
  mbounf.F90
  mboxb.F90
  mchunk.F90
  mcout.F90
  means.F90
  meansqs.F90
  meansqs_gc.F90
  meansqs_lf.F90
  mfr.F90
  mfredir.F90
  mgrid.F90
  micep.F90
  mintf.F90
  mnintw.F90
  mpabort.F90
  mpbcastgrid.F90
  mpbcastintfld.F90
  mpclose_unit.F90
  mpcrtbl.F90
  mpdecomp.F90
  mpdistribfl.F90
  mpdistribscfld.F90
  mpexchng.F90
  mpfldtoifs.F90
  mpgatherbc.F90
  mpgatherfl.F90
  mpgatherscfld.F90
  mpminmaxavg.F90
  mpuserin.F90
  mstart.F90
  mswell.F90
  mubuf.F90
  mwp1.F90
  mwp2.F90
  newwind.F90
  nlweigt.F90
  notim.F90
  ns_gc.F90
  omegagc.F90
  out_onegrdpt.F90
  out_onegrdpt_sp.F90
  outbc.F90
  outbeta.F90
  outblock.F90
  outbs.F90
  outcom.F90
  outgrid.F90
  outint.F90
  outmdldcp.F90
  outnam.F90
  outpp.F90
  outsetwmask.F90
  outspec.F90
  outstep0.F90
  outwint.F90
  outwnorm.F90
  outwpsp.F90
  outwspec.F90
  packi.F90
  packr.F90
  parkind_wave.F90
  parmean.F90
  peak.F90
  peak_ang.F90
  peak_freq.F90
  peakfri.F90
  preset_wgrib_template.F90
  prewind.F90
  proenvhalo.F90
  propag_wam.F90
  propags.F90
  propags1.F90
  propags2.F90
  propconnect.F90
  propdot.F90
  readbou.F90
  readfl.F90
  readmdlconf.F90
  readpre.F90
  readsta.F90
  readstress.F90
  readwgrib.F90
  readwind.F90
  recvnemofields.F90
  rotspec.F90
  runwam.F90
  savspec.F90
  savstress.F90
  sbottom.F90
  scosfl.F90
  sdepthlim.F90
  sdissip.F90
  sdissip_ard.F90
  sdissip_jan.F90
  sdiwbk.F90
  se10mean.F90
  sebtmean.F90
  second_order_lib.F90
  secondhh.F90
  secondhh_gen.F90
  secspom.F90
  semean.F90
  sep3tr.F90
  sepwisw.F90
  set_wflags.F90
  setice.F90
  setmarstype.F90
  setwavphys.F90
  sinflx.F90
  sinput.F90
  sinput_ard.F90
  sinput_jan.F90
  skewness.F90
  snonlin.F90
  spectra.F90
  spr.F90
  stat_nl.F90
  sthq.F90
  stokesdrift.F90
  stokestrn.F90
  stress_gc.F90
  stresso.F90
  strspec.F90
  tables_2nd.F90
  tabu_swellft.F90
  tau_phi_hf.F90
  taut_z0.F90
  topoar.F90
  transf.F90
  transf_bfi.F90
  transf_r.F90
  transf_snl.F90
  uibou.F90
  uiprep.F90
  unsetice.F90
  updnemofields.F90
  updnemostress.F90
  userin.F90
  vmin.F90
  vmin_d.F90
  vplus.F90
  vplus_d.F90
  w_cdf.F90
  w_cdfmax.F90
  w_maxh.F90
  w_mode_st.F90
  w_pdf.F90
  w_pmax.F90
  wam_multio_mod.F90
  wam_nproma.F90
  wam_sorti.F90
  wam_sortini.F90
  wam_u2l1cr.F90
  wam_user_clock.F90
  wamadswstar.F90
  wamcur.F90
  wamintgr.F90
  wamodel.F90
  wamwnd.F90
  wavemdl.F90
  wdfluxes.F90
  wdirspread.F90
  weflux.F90
  wgrib2fdb.F90
  wgrib_edition.F90
  wgribencode.F90
  wgribencode_model.F90
  wgribenout.F90
  wgribout.F90
  wnfluxes.F90
  wposnam.F90
  writefl.F90
  writestress.F90
  writsta.F90
  wsigstar.F90
  wsmfen.F90
  wstream_strg.F90
  wvalloc.F90
  wvdealloc.F90
  wvfricvelo.F90
  wvwamdecomp.F90
  wvwaminit.F90
  wvwaminit1.F90
  yowabort.F90
  yowaltas.F90
  yowassi.F90
  yowcard.F90
  yowcinp.F90
  yowcoer.F90
  yowconst_2nd.F90
  yowcoup.F90
  yowcout.F90
  yowcpbo.F90
  yowcurg.F90
  yowcurr.F90
  yowdes.F90
  yowdrvtype.F90
  yowfpbo.F90
  yowfred.F90
  yowgrib.F90
  yowgrib_handles.F90
  yowgribhd.F90
  yowgrid.F90
  yowgstats.F90
  yowice.F90
  yowindn.F90
  yowintp.F90
  yowjons.F90
  yowmap.F90
  yowmean.F90
  yowmespas.F90
  yowmpp.F90
  yownemoflds.F90
  yownemoio.F90
  yowparam.F90
  yowpcons.F90
  yowphys.F90
  yowprproc.F90
  yowrefd.F90
  yowshal.F90
  yowspec.F90
  yowsphere.F90
  yowstat.F90
  yowtabl.F90
  yowtemp.F90
  yowtest.F90
  yowtext.F90
  yowtrains.F90
  yowubuf.F90
  yowunit.F90
  yowwami.F90
  yowwind.F90
  yowwndg.F90
  z0wave.F90
)

# check for mixed precision build
set(MIXED_PRECISION 0)
if(${HAVE_SINGLE_PRECISION} AND NOT ${PNAME}_OCEANMODEL_HAVE_SINGLE_PRECISION)
   set(MIXED_PRECISION 1)
endif()

if( HAVE_GEN_DERIV_TYPES )
   list( APPEND ecwam_srcs ${CMAKE_CURRENT_BINARY_DIR}/yowfield_mod.F90)
   add_custom_command(
       OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/yowfield_mod.F90
       COMMAND ${FYPP} -m io -m yaml -m os -DMIXED_PRECISION=${MIXED_PRECISION} ${CMAKE_CURRENT_SOURCE_DIR}/yowfield_mod.fypp > yowfield_mod.F90
       DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/yowfield_mod.fypp
       VERBATIM)
else()
   list( APPEND ecwam_srcs yowfield_mod.F90)
endif()

unset( MPI_Fortran_LIBRARIES )
if( HAVE_UNWAM )
  list(APPEND ecwam_srcs ${unwam_srcs})
  list(APPEND ECWAM_DEFINITIONS WAM_HAVE_UNWAM)
  if( HAVE_MPI )
      if( NOT DEFINED HAVE_MPI_F08 )
         set( HAVE_MPI_F08 1 )
      endif()
      if( DEFINED ENV{OPENMPI_VERSION} )
        set( OPENMPI_VERSION $ENV{OPENMPI_VERSION} )
        if( OPENMPI_VERSION VERSION_LESS 4.0 )
           ecbuild_warn("There are compilation problems with F08 MPI API with OpenMPI < 4.0. Defaulting to F77 API.")
           set( HAVE_MPI_F08 0 )
        endif()
      endif()
      if( HAVE_MPI_F08 )
        list(APPEND ECWAM_DEFINITIONS WAM_HAVE_MPI_F08 )
      endif()
      set( MPI_Fortran_LIBRARIES MPI::MPI_Fortran )
  else()
      set( MPI_Fortran_LIBRARIES ${MPI_SERIAL_LIBRARIES} )
  endif()
endif()

set(MULTIO_LIBRARIES)
if( HAVE_MULTIO )
  if(multio_FOUND)
    if(multio_HAVE_FDB5)
      # Force the linkage of multio-fdb5 for static builds
      if(NOT BUILD_SHARED_LIBS)
        list(APPEND MULTIO_LIBRARIES -Wl,--push-state,--no-as-needed multio-fapi -Wl,--pop-state)
      else()
        list(APPEND MULTIO_LIBRARIES multio-fapi)
      endif()
    endif()
    list(APPEND MULTIO_LIBRARIES multio-fapi)
    list(APPEND ECWAM_DEFINITIONS WAM_HAVE_MULTIO)
  endif()
endif()

if( HAVE_OCEAN_COUPLING )
  list(APPEND ECWAM_DEFINITIONS WITH_NEMO )
  if( ${PNAME}_OCEANMODEL_HAVE_SINGLE_PRECISION )
    list(APPEND ECWAM_DEFINITIONS PARKIND1_SINGLE_NEMO )
  endif()
else()
  unset( ${PNAME}_OCEANMODEL_LIBRARIES )
  unset( ${PNAME}_OCEANMODEL_INCLUDE_DIRS )
endif()

if( HAVE_ECFLOW )
  list(APPEND ECWAM_PRIVATE_DEFINITIONS WAM_HAVE_ECFLOW)
endif()

ecbuild_add_library(
    TARGET           ${ecwam}
    SOURCES          ${ecwam_srcs}
    PUBLIC_LIBS      fiat parkind_${prec} ${ecwam}_intfb
                     ${MPI_Fortran_LIBRARIES}
                     ${${PNAME}_OCEANMODEL_LIBRARIES}
    PRIVATE_LIBS     eccodes_f90
                     ${MULTIO_LIBRARIES}
                     ${OpenMP_Fortran_LIBRARIES}
                     field_api_${prec}
                     $<${HAVE_ECFLOW}:ecflow_lightf>
    PUBLIC_INCLUDES  $<INSTALL_INTERFACE:include>
    PRIVATE_INCLUDES ${${PNAME}_OCEANMODEL_INCLUDE_DIRS}
    PRIVATE_DEFINITIONS ${ECWAM_PRIVATE_DEFINITIONS}
    PUBLIC_DEFINITIONS ${ECWAM_DEFINITIONS}
)

ecwam_target_fortran_module_directory(
    TARGET ${ecwam}
    MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/${ecwam}
    INSTALL_DIRECTORY module/${ecwam}
)

ecwam_target_compile_definitions_FILENAME( ${ecwam} )

### The file mubuf.F90, which is only used for "preproc" is sensitive to optimisations
#   possibly leading to different wam_grid_<1,2,3> files.
#   This in turn leads to non-neglibible differences
#   of average 'swh' when running "chief".

if( CMAKE_Fortran_COMPILER_ID MATCHES Intel )
  set_source_files_properties( mubuf.F90 PROPERTIES COMPILE_OPTIONS "-fp-model;strict" )
  set_source_files_properties( propconnect.F90 PROPERTIES COMPILE_OPTIONS "-fp-model;strict" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES GNU )
  set_source_files_properties( mubuf.F90 PROPERTIES COMPILE_OPTIONS "-ffp-contract=off" )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
  set_source_files_properties(
      sbottom.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros -Mstack_arrays "
  )
endif()

### The file grib2wgrid.F90 is sensitive to optimizations in single precision builds.
#   This leads to non-neglibible differences
#   of average 'swh' when running "chief".

if( CMAKE_Fortran_COMPILER_ID MATCHES Intel AND HAVE_SINGLE_PRECISION )
  set_source_files_properties( grib2wgrid.F90 PROPERTIES COMPILE_OPTIONS "-fp-model;strict" )
endif()

### NVHPC debug flags applied selectively to sourcefiles to avoid single-precision overflow
if( CMAKE_Fortran_COMPILER_ID MATCHES PGI|NVHPC AND CMAKE_BUILD_TYPE MATCHES "Debug" )
  foreach( debug_flag    initsnan checkbounds )
    if( ${debug_flag}_flags )
      ecbuild_add_fortran_flags( "${${debug_flag}_flags}" NAME ${debug_flag} )
    endif()
  endforeach()
  
  list(REMOVE_ITEM ecwam_srcs outbeta.F90 secondhh.F90)
  if( fpe_flags )
     set_source_files_properties( ${ecwam_srcs} PROPERTIES COMPILE_OPTIONS ${fpe_flags} )
     set_source_files_properties( outbeta.F90 PROPERTIES COMPILE_OPTIONS "-Ktrap=divz" )
     set_source_files_properties( secondhh.F90 PROPERTIES COMPILE_OPTIONS "-Ktrap=inv,ovf" )
  endif()
  list(APPEND ecwam_srcs outbeta.F90 secondhh.F90)
endif()
