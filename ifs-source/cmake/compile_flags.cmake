# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction


# noopt:        -DCMAKE_BUILD_TYPE=Debug
# NANS_C:       -DCMAKE_BUILD_TYPE=Bit   -DIFS_CHECK_BOUNDS=ON -DIFS_INIT_SNAN=ON
# noopt_NANS_C: -DCMAKE_BUILD_TYPE=Debug -DIFS_CHECK_BOUNDS=ON -DIFS_INIT_SNAN=ON

if(IFS_CHECK_BOUNDS)
  # Files that fail bounds checking across multiple compilers.
  # Only files generating false positives should be added to this list.
  list(APPEND no_bounds_checking
    arpifs/adiab/cpg5_gp.F90
    arpifs/adiab/cpg_gp.F90
    arpifs/adiab/cpg_gp_ad.F90
    arpifs/adiab/cpg_gp_hyd.F90
    arpifs/adiab/cpg_gp_tl.F90
    arpifs/adiab/lapinea.F90
    arpifs/adiab/lapinea5.F90
    arpifs/adiab/lapineaad.F90
    arpifs/adiab/lapineatl.F90
    arpifs/adiab/larcinb.F90
    arpifs/adiab/larmes.F90
    arpifs/adiab/larmes5.F90
    arpifs/adiab/larmesad.F90
    arpifs/adiab/larmestl.F90
    arpifs/adiab/lattex.F90
    arpifs/adiab/postphy.F90
    arpifs/chem/chem_bascoetm5.F90
    arpifs/chem/chem_main.F90
    arpifs/chem/chem_massdia.F90
    arpifs/chem/tm5_rbud.F90
    arpifs/phys_radi/radintg.F90
  )
endif()


if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  set(autopromote_flags -sreal64)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(autopromote_flags -fdefault-real-8 -fdefault-double-8)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(autopromote_flags -real-size 64)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
  set(autopromote_flags -r8)
endif()

set( IFS_LINK_FLAGS ${CMAKE_SHARED_LINKER_FLAGS_INIT} )

if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")

  find_package( OpenMP COMPONENTS Fortran REQUIRED )
  if( TARGET OpenMP::OpenMP_Fortran )
    set( IFS_OMP_Fortran_LIBRARIES OpenMP::OpenMP_Fortran )
  endif()

  set(IFS_Fortran_FLAGS "${IFS_Fortran_FLAGS} -ram -emf") # common flags for all build types

  if(IFS_CHECK_BOUNDS)
    set(IFS_Fortran_FLAGS "${IFS_Fortran_FLAGS} -Rb")
  endif()

  if(IFS_INIT_SNAN)
    set(IFS_Fortran_FLAGS "${IFS_Fortran_FLAGS} -ei")
  endif()

  if(IFS_SHARED_LIBS)
    set(PIC_FLAGS "-fPIC -hPIC")
  endif()

  set(IFS_Fortran_FLAGS_BIT "-hflex_mp=conservative -Othread1 -hfp1 -hadd_paren -homp") # no debug symbols
  set(IFS_Fortran_FLAGS_DEBUG "-G0 -O0 -hflex_mp=conservative -hfp0 -hadd_paren -homp")

  if($ENV{CRAY_FTN_VERSION} MATCHES "8.7.0")
    # FMA with 8.7.0 breaks the adjoint test
    set(IFS_Fortran_FLAGS_BIT "${IFS_Fortran_FLAGS_BIT} -hnofma") # no fused muiltply-add instructions
  endif()

  if($ENV{CRAY_FTN_VERSION} MATCHES "8.5.6|8.5.8")
    set(IFS_Fortran_FLAGS_BIT "${IFS_Fortran_FLAGS_BIT} -hipa0") # disable interprocedural analysis
  endif()

  set(IFS_C_FLAGS_BIT "-O2 -fopenmp")
  set(IFS_C_FLAGS_DEBUG "-g -O0 -fopenmp")

  set(IFS_CXX_FLAGS_BIT "-g -O1 -fopenmp")
  set(IFS_CXX_FLAGS_DEBUG "-g -O0 -fopenmp")

  #set(ECBUILD_SHARED_LINKER_FLAGS "-Wl,--eh-frame-hdr -Wl,--disable-new-dtags -Ktrap=fp -homp")
  set(ECBUILD_SHARED_LINKER_FLAGS "-Wl,--eh-frame-hdr -Wl,--disable-new-dtags")
  set(ECBUILD_EXE_LINKER_FLAGS "${ECBUILD_SHARED_LINKER_FLAGS} -Wl,--as-needed -hbyteswapio")

  # Make sure we link with hugetlbfs for hugepages on Cray/HPE systems
  if($ENV{CRAY_FTN_VERSION} VERSION_GREATER_EQUAL 15.0.1)
    set(ECBUILD_EXE_LINKER_FLAGS "${ECBUILD_SHARED_LINKER_FLAGS} -lhugetlbfs")
  endif()

  # Some macro expansions generate very long lines
  file(GLOB_RECURSE long_line_srcs RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
       "${CMAKE_CURRENT_SOURCE_DIR}/satrad/rttov/*/*.F*" "${CMAKE_CURRENT_SOURCE_DIR}/odb/*/*.F*")
  set_source_files_properties(${long_line_srcs}
    PROPERTIES COMPILE_FLAGS "-N1023")

  set_source_files_properties(arpifs/mwave/mwave_emis.F90
    arpifs/mwave/mwave_emis_cmem.F90
    arpifs/phys_radi/raddiag.F90 arpifs/setup/suvert.F90
    arpifs/phys_ec/vdfexcus.F90 arpifs/phys_ec/vdfexcusad.F90 arpifs/phys_ec/vdfexcustl.F90
    surf/external/surftstp.F90 surf/module/surftstp_ctl_mod.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -O0,fp1,omp ${PIC_FLAGS}")

  set_source_files_properties(arpifs/module/spectral_arp_mod.F90 arpifs/module/varbc_allsky.F90
    arpifs/utility/deallo.F90 arpifs/cma2odb/create_averaged_values.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "-g ${IFS_Fortran_FLAGS} -hflex_mp=conservative -hfp1 -hadd_paren ${PIC_FLAGS}")

  set_source_files_properties(arpifs/adiab/larmes.F90
    arpifs/adiab/larmes_xyz.F90
    arpifs/phys_ec/phys_arrays_ini.F90 arpifs/phys_ec/cloudsc.F90
    PROPERTIES COMPILE_FLAGS "-hcontiguous")

  set_source_files_properties(arpifs/phys_ec/vdfouter.F90
    PROPERTIES COMPILE_FLAGS "-Oloop_trips=small")

  set_source_files_properties(arpifs/phys_ec/local_state_ini.F90
    PROPERTIES COMPILE_FLAGS "-hcontiguous -hnopattern")

  set_source_files_properties(arpifs/adiab/larche.F90
    PROPERTIES COMPILE_FLAGS "-Oshortcircuit1 -Oloop_trips=small -hcontiguous")

  set_source_files_properties(wam/Wam_oper/propags2.F
    PROPERTIES COMPILE_FLAGS "-Oloop_trips=small")

  set_source_files_properties(arpifs/adiab/cpg.F90
    arpifs/phys_ec/state_update.F90 arpifs/phys_ec/state_increment.F90 arpifs/phys_ec/state_copy.F90
    PROPERTIES COMPILE_FLAGS "-Onopattern -hcontiguous")

  set_source_files_properties(arpifs/phys_radi/radlswr.F90
    PROPERTIES COMPILE_FLAGS "-Ovector0")

  set_source_files_properties(arpifs/phys_radi/lwad.F90
    arpifs/phys_radi/rrtm_ecrt_140gp_mcica.F90 arpifs/phys_radi/srtm_srtm_224gp_mcica.F90
    PROPERTIES COMPILE_FLAGS "-Onopattern")

  set_source_files_properties(wam/Wam_oper/secspom.F
    PROPERTIES COMPILE_FLAGS "-hloop_trips=small")

  set_source_files_properties(arpifs/interpol/lascaw.F90
    PROPERTIES COMPILE_FLAGS "-Othread2 -Oshortcircuit1 -hcontiguous")

  set_source_files_properties(arpifs/parallel/brptob.F90
    PROPERTIES COMPILE_FLAGS "-hnopattern -hcontiguous")

  set_source_files_properties(arpifs/phys_radi/radlswad.F90
    PROPERTIES COMPILE_FLAGS "-Onopattern -hcontiguous")

  # Pattern matching injects libsci routines, against which these tools are not linked
  set_source_files_properties(prepdata/mc_tools/find_target_area.F90
    PROPERTIES COMPILE_FLAGS "-hnopattern")

  # Optimised code takes wrong branch leading to division by zero in odb/lib/aggr.c
  set_source_files_properties(lib/aggr.c
    PROPERTIES COMPILE_FLAGS "-O0")

  # IFS-864 Chemistry solver is not bit-reproducible between NPES settings without
  # either -fp0 or -hcpu=ivybridge (as at 46r1 and crayftn 8.5.8)
  set_source_files_properties(arpifs/chem/tm5_do_ebi.F90 arpifs/chem/tm5_do_ebi_tc02b.F90
    PROPERTIES COMPILE_FLAGS "-hcpu=ivybridge")

  # Fix internal compiler error
  set_source_files_properties(satrad/programs/bufr_grid_screen.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -hflex_mp=intolerant -Othread1 -hfp1 -hadd_paren"
               OVERRIDE_COMPILE_FLAGS_DEBUG "${IFS_Fortran_FLAGS} -g -hflex_mp=intolerant -hfp0 -hadd_paren")

  # Compile MPL sources with -g for better traceback
  file(GLOB mpl_srcs RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/ifsaux/module/mpl_*.F90")
  set_source_files_properties(${mpl_srcs} PROPERTIES COMPILE_FLAGS "-g")

  # Fix for "CONGRAD: SPTSV/DPTSV returned non-zero info with crayftn 8.7.7 (cdt/18.12)
  if($ENV{CRAY_FTN_VERSION} VERSION_GREATER 8.7.6)
    set_source_files_properties(trans/module/ftinv_ctlad_mod.F90
    PROPERTIES COMPILE_FLAGS "-O0,fp1,omp")
  endif()

  # Fix for FP Overflow in CY47R2
  if($ENV{CRAY_FTN_VERSION} VERSION_GREATER 8.7.6)
    set_source_files_properties(arpifs/dia/cpdyddh.F90
    PROPERTIES COMPILE_FLAGS "-O0,fp1,omp")
  endif()

  # Fix for FPE in CY47R2 ifs-test/tests/t21/test_4dvar_airep_t
  if($ENV{CRAY_FTN_VERSION} VERSION_GREATER 8.7.6)
    set_source_files_properties(arpifs/dia/spnorm.F90
    PROPERTIES COMPILE_FLAGS "-O0,fp1,omp")
  endif()

  # Fix for internal compiler errors with crayftn 8.7.7 (cdt/18.12)
  if($ENV{CRAY_FTN_VERSION} VERSION_GREATER 8.7.6)
    set_source_files_properties(arpifs/module/surface_views_module.fypp
    PROPERTIES COMPILE_FLAGS "-O0,fp0,noomp")

    set_source_files_properties(
    satrad/rttov/main/rttov_alloc_sunglint.F90 satrad/rttov/main/rttov_alloc_transmission_aux.F90
    PROPERTIES COMPILE_FLAGS "-O0,fp1,omp")

    string(TOUPPER ${CMAKE_BUILD_TYPE} btype)
    string(REGEX REPLACE "-Rb" "" flags "${IFS_Fortran_FLAGS} ${IFS_Fortran_FLAGS_${btype}}")
    set_source_files_properties(arpifs/module/traj_semilag_mod.F90 arpifs/module/traj_physics_mod.F90 enkf/module/obs_base_mod.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${flags} ${PIC_FLAGS}")
 endif()

  # IFS-1933 Fix for adjoint test
  set_source_files_properties(arpifs/phys_ec/cuascn.F90
    PROPERTIES COMPILE_FLAGS "-hcpu=ivybridge")

  # Fix for adjoint problem with crayftn 8.7.7 (cdt/18.12)
  if($ENV{CRAY_FTN_VERSION} VERSION_GREATER 8.7.6)
    set_source_files_properties(
    surf/module/surfexcdrivers_ctl_mod.F90 surf/module/surfexcdriverstl_ctl_mod.F90
    surf/module/surfexcdriversad_ctl_mod.F90 surf/module/surfexcdriver_ctl_mod.F90
    PROPERTIES COMPILE_FLAGS "-O0,fp1,omp")
  endif()

  # Fix for an exception with ifs-test/tests/t21/test_glomap_edge_fc with crayftn 8.7.7 (cdt/18.23)
  if($ENV{CRAY_FTN_VERSION} VERSION_GREATER 8.7.6)
    set_source_files_properties(arpifs/glomap_mode/ukca_so2so4.F90
    PROPERTIES COMPILE_FLAGS "-O0,fp1")
  endif()

  # Routines with very long compile times with crayftn 8.5.8 (cdt/17.03)
  if($ENV{CRAY_FTN_VERSION} VERSION_GREATER 8.5.7)
    file(GLOB fa_srcs RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/ifsaux/fa/*.F90")
    set_source_files_properties(${fa_srcs} arpifs/op_obs/bgobs.F90
    arpifs/control/cnt0.F90 arpifs/control/cprep3.F90
    arpifs/module/field_registry_mod.fypp
    arpifs/phys_dmn/aplpar.F90 arpifs/phys_dmn/apl_arome.F90
    arpifs/op_obs/departure_jo.F90 arpifs/op_obs/departure_joad.F90 arpifs/op_obs/departure_jotl.F90
    arpifs/op_obs/hop.F90 arpifs/op_obs/hretr_aeolus.F90 arpifs/obs_preproc/mw_clearsky_screen_wrapper.F90
    arpifs/setup/suafn.F90 arpifs/setup/suafn1.F90 arpifs/setup/su0yomb.F90 arpifs/var/sujbtest.F90
    arpifs/oops/fields_io_mod.F90 arpifs/control/forecast_error.F90 arpifs/setup/su0yomb.F90
    arpifs/obs_preproc/gen_corr_pert.F90 ifsobs/src/misc/sql_expression_mod.F90
    ifsobs/src/misc/sql_query_mod.F90 ifsobs/src/misc/sql_dbase_mod.F90
    arpifs/obs_preproc/screen.F90 arpifs/op_obs/hretr_aeolus.F90
    arpifs/obs_preproc/mw_clearsky_screen_wrapper.F90
    arpifs/obs_preproc/decis.F90
    nemo/NEMOGCM/NEMO/OPA_SRC/LBC/mppini.F90 nemo/NEMOGCM_V40/src/OCE/LBC/mppini.F90
    PROPERTIES COMPILE_FLAGS "-O0,fp1")
  endif()

  if($ENV{CRAY_FTN_VERSION} VERSION_GREATER 8.5.7)
    set_source_files_properties(arpifs/control/cprep4.F90
    PROPERTIES COMPILE_FLAGS "-O0,fp1 -hnoomp")
  endif()

  # obs_preproc suffers from very long compile times with Cray 15.0.0/15.0.1 on LUMI
  if($ENV{CRAY_FTN_VERSION} VERSION_GREATER_EQUAL 15.0.0)
    file(GLOB obs_preproc_srcs RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/arpifs/obs_preproc/*.F90")
    set_source_files_properties(${obs_preproc_srcs}
      PROPERTIES COMPILE_FLAGS "-O0,fp1")
  endif()

  # Fix for LREPRO4DVAR in crayftn 8.6.2 (cdt/17.09)
  if($ENV{CRAY_FTN_VERSION} VERSION_GREATER 8.6.1)
    set_source_files_properties(arpifs/adiab/lattex.F90
      PROPERTIES COMPILE_FLAGS "-O0,fp1,omp")
  endif()

  # Fix for LREPRO4DVAR in crayftn 8.5.8 (cdt/17.03)
  if($ENV{CRAY_FTN_VERSION} VERSION_GREATER 8.5.7)
    set_source_files_properties(arpifs/phys_ec/cumastrn.F90 arpifs/phys_ec/cuadjtqs.F90
    arpifs/phys_ec/cuadjtqsad.F90 arpifs/phys_ec/cuadjtqstl.F90
    arpifs/phys_radi/swniad.F90 arpifs/phys_radi/swnitl.F90 arpifs/var/pregprh.F90
    PROPERTIES COMPILE_FLAGS "-O0,fp1,omp,scalar2")
  endif()

  # Fix for bugs in crayftn 8.4.5 (cdt/16.03) and 8.5.8 (cdt/17.03)
  if($ENV{CRAY_FTN_VERSION} MATCHES "8.4.5|8.5.8")
    file(GLOB_RECURSE prepdata_srcs RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/prepdata/*.F*")
    set_source_files_properties(${prepdata_srcs}
      ssa/plot/coordinates.F90 ssa/plot/fdb_output.F90 ssa/plot/getfields.F90 ssa/plot/print_nml.F90
      ssa/sub/inisnw.F90 ssa/sub/inisst.F90 ssa/sub/init2m.F90  ssa/sub/reg_to_gg.F90 ssa/util/setcomssa.F90
      satrad/programs/bufr_screen_nexrad.F90 satrad/programs/bufr_screen_nexrad.F90 satrad/programs/bufr_screen_opera.F90
      satrad/programs/bufr_screen_synop_rain_gauges.F90 satrad/programs/calc_radiance_fields.F90
      satrad/programs/eda_rad_scale.F90 satrad/programs/gensatim.F90
      PROPERTIES COMPILE_FLAGS "-hipa0")
    set_source_files_properties(arpifs/interpol/slcomm2a.F90 arpifs/pp_obs/pos.F90
      PROPERTIES COMPILE_FLAGS "-g")
  endif()

  if($ENV{CRAY_FTN_VERSION} VERSION_GREATER 8.4.0) # cdt/15.11

    # Fix for adjoint test
    set_source_files_properties(surf/module/vsurf_mod.F90 surf/module/vsurfs_mod.F90
      surf/module/vsurfstl_mod.F90 surf/module/vsurfsad_mod.F90 surf/module/vupdz0_mod.F90
      surf/module/vupdz0s_mod.F90 surf/module/vupdz0stl_mod.F90 surf/module/vupdz0sad_mod.F90
      surf/module/surfsebs_ctl_mod.F90 surf/module/surfsebstl_ctl_mod.F90 surf/module/surfsebsad_ctl_mod.F90
      surf/module/surfrad_ctl_mod.F90
      PROPERTIES COMPILE_FLAGS "-hcpu=ivybridge")

    # Fix for hang at MPI_finialize in 43R3 with crayftn 8.4.1 (in cdt/15.11)
    set_source_files_properties(arpifs/io_serv/io_serv_flush.F90
      PROPERTIES COMPILE_FLAGS "-g")

    # Fix for double free error
    string(TOUPPER ${CMAKE_BUILD_TYPE} btype)
    string(REGEX REPLACE "-g|-G[0-2]|-Gfast|-Rb" "" flags "${IFS_Fortran_FLAGS} ${IFS_Fortran_FLAGS_${btype}}")
    set_source_files_properties(arpifs/adiab/gp_derivatives.F90
      arpifs/setup/sump.F90 arpifs/dia/wrspeca.F90 arpifs/module/iogrida_mod.F90
      PROPERTIES OVERRIDE_COMPILE_FLAGS "${flags} ${PIC_FLAGS}")

  endif()

  if($ENV{CRAY_FTN_VERSION} VERSION_GREATER 8.4.5) # cdt/16.04

    # Fix for crash in  __pgas_runtime_error_checking with crayftn 8.4.6 (in cdt/16.04)
    set_source_files_properties(ifsaux/module/oml_mod.F90 arpifs/module/varbc_pred.F90
      PROPERTIES COMPILE_FLAGS "-hnocaf")

    set_source_files_properties(radiation/module/easy_netcdf_read_mpi.F90 radiation/module/easy_netcdf.F90
      PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} ${IFS_Fortran_FLAGS_DEBUG} ${PIC_FLAGS}")

  endif()

  set_source_files_properties(
    climfield/ifs_tools/grib_set_vtable.F90
    climfield/src/checkgg.F90
    climfield/ifs_tools/spint_special_filter.F90
    climfield/ifs_tools/cheminterpol.F90
    climfield/ifs_tools/cheminterpol_vms.F90
    PROPERTIES
      OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -hflex_mp=intolerant -hfp0"
      OVERRIDE_COMPILE_FLAGS_DEBUG "${IFS_Fortran_FLAGS} -g -hflex_mp=intolerant -hfp0")

  set_source_files_properties(
    climfield/ifs_tools/depth_mode_filter.F90
    PROPERTIES COMPILE_FLAGS "-N 255")

  # Following compiled without "-hnopattern" flag leads to undefined reference to libsci_cray symbol _idmax_ e.g. on LUMI with cce/14
  # A compiler optimisation tries to replace a pattern with a blas routine but we are not linking these routines with BLAS.
  set_source_files_properties (
    aeolus/OpticalProperties/opticalproperties.F90
    aeolus/HLOS_retrieval/hlos_retrieval.F90
    prepdata/kpp/nearneighbor_639l_2.F90
    prepdata/kpp/nearneighbor_399l_2.F90
    prepdata/kpp/nearneighbor_255l_2.F90
    PROPERTIES COMPILE_FLAGS "-hnopattern" )

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

  if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "aarch64")
    set(IFS_GNU64_FLAG "")
  else()
    set(IFS_GNU64_FLAG "-m64")
  endif()

  string(CONCAT IFS_Fortran_FLAGS "${IFS_Fortran_FLAGS} ${IFS_GNU64_FLAG} "
      "-fconvert=big-endian -fPIC -fopenmp "
      "-fno-range-check -ffree-line-length-none -fbacktrace -fno-second-underscore "
      "-ffpe-trap=invalid,zero,overflow -fconvert=swap")

  # gfortran 10 has become stricter with argument matching
  if( NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10 )
    string(CONCAT IFS_Fortran_FLAGS "${IFS_Fortran_FLAGS} -fallow-argument-mismatch -fallow-invalid-boz")
  endif()

  if(IFS_CHECK_BOUNDS)
    set(IFS_Fortran_FLAGS "${IFS_Fortran_FLAGS} -fcheck=bounds")

    # Disable bounds checking on files that would fail the check.
    # Only files generating false positives should be added to this list.
    set_source_files_properties(
      ${no_bounds_checking}
      arpifs/adiab/call_sl.F90
      arpifs/adiab/cpg_drv.F90
      arpifs/adiab/cpg_drv_ad.F90
      arpifs/adiab/cpg_drv_tl.F90
      arpifs/glomap_mode/aer_dust.F90
      arpifs/module/ecphys_perturb_type_mod.F90
      arpifs/phys_ec/ec_phys_ad.F90
      arpifs/phys_ec/ec_phys_drv.F90
      arpifs/phys_ec/ec_phys_drv_tl.F90
      arpifs/phys_ec/ec_phys_tl.F90
      PROPERTIES COMPILE_FLAGS "-fcheck=no-bounds")
  endif()

  if(IFS_INIT_SNAN)
    set(IFS_Fortran_FLAGS "${IFS_Fortran_FLAGS} -finit-real=snan")
  endif()

  set(IFS_Fortran_FLAGS_BIT "-g -O2")
  set(IFS_Fortran_FLAGS_DEBUG "-g -O0")

  set(IFS_C_FLAGS_BIT "-g ${IFS_GNU64_FLAG}")
  set(IFS_C_FLAGS_DEBUG "-g -O0 ${IFS_GNU64_FLAG}")

  if( CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
    set(IFS_CXX_FLAGS_BIT "-g -O1 ${IFS_GNU64_FLAG}")
    set(IFS_CXX_FLAGS_DEBUG "-g -O0 ${IFS_GNU64_FLAG}")
  endif()

  find_package( OpenMP COMPONENTS Fortran REQUIRED )
  if( TARGET OpenMP::OpenMP_Fortran )
    set( IFS_OMP_Fortran_LIBRARIES OpenMP::OpenMP_Fortran )
  endif()

  if( APPLE )

    if( NOT IFS_OMP_Fortran_LIBRARIES )
      ecbuild_error( "OpenMP libraries not found" )
    endif()

    # Remove linker warnings of sort "ld: warning: could not create compact unwind for ..."
    set( IFS_LINK_FLAGS "${IFS_LINK_FLAGS} -Wl,-keep_dwarf_unwind -Wl,-no_compact_unwind")

  else()

    set( IFS_LINK_FLAGS "${IFS_LINK_FLAGS} -Wl,--eh-frame-hdr -Wl,--disable-new-dtags -fopenmp" )

  endif()

  set( ECBUILD_SHARED_LINKER_FLAGS "${IFS_LINK_FLAGS}" )
  set( ECBUILD_EXE_LINKER_FLAGS    "${IFS_LINK_FLAGS}" )
  if( NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin") # Darwin linker does not support --as-needed
    set( ECBUILD_EXE_LINKER_FLAGS "${ECBUILD_EXE_LINKER_FLAGS} -Wl,--as-needed" )
    set( no_as_needed "-Wl,--no-as-needed" )  # For individual execs as required
  endif()

  if( $ENV{HOSTNAME} MATCHES "lxc.*|lxg.*")
    set( ECBUILD_EXE_LINKER_FLAGS "${ECBUILD_EXE_LINKER_FLAGS} -Wl,--allow-shlib-undefined" )
    set( ECBUILD_MODULE_LINKER_FLAGS "${ECBUILD_MODULE_LINKER_FLAGS} -Wl,--allow-shlib-undefined" )
    set( ECBUILD_SHARED_LINKER_FLAGS "${ECBUILD_SHARED_LINKER_FLAGS} -Wl,--allow-shlib-undefined" )
  endif()

  # Fixes for Ubuntu Docker container running OpenMPI 4.1
  # Without these changes, ifstest fails on the ubuntu dockers
  # None of these routines use OpenMP directly, so apply generally.

  set_source_files_properties(
    # Fix most ifs-test configurations
    arpifs/control/cnt0.F90 arpifs/setup/su0yomb.F90 arpifs/fullpos/sufpdyn.F90
    # Fix adjoint & tangent tests
    arpifs/var/suscal.F90
    # Fix 4dvar tests
    arpifs/var/suecges.F90 arpifs/var/suinfce.F90
    # Fix SCM test
    scmec/source/cnt1c.F90
    PROPERTIES COMPILE_FLAGS "-fno-openmp")

  # -O1 flag to prevent runtime hang
  # -fno-openmp for Docker issues (see below)
  set_source_files_properties(arpifs/control/cnt0.F90
    PROPERTIES COMPILE_FLAGS "-O1 -fno-openmp")

  # Fixes for Ubuntu Docker container running OpenMPI 4.1
  # None of these routines use OpenMP directly, so apply generally.
  set_source_files_properties(
    # Fix most ifs-test configurations
    arpifs/setup/su0yomb.F90 arpifs/fullpos/sufpdyn.F90
    # Fix adjoint & tangent tests
    arpifs/var/suscal.F90
    # Fix 4dvar tests
    arpifs/var/suecges.F90 arpifs/var/suinfce.F90
    # Fix SCM
    scmec/source/cnt1c.F90
    PROPERTIES COMPILE_FLAGS "-fno-openmp")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")

  # "classic" ifort compiler, else LLVM-based ifx compiler
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
    string(CONCAT IFS_Fortran_FLAGS "${IFS_Fortran_FLAGS} -fpe0 -convert big_endian -align array64byte "
      "-assume byterecl,noold_maxminloc,realloc_lhs "
      "-traceback -fpic -qopenmp -qopenmp-threadprivate compat -fp-model precise "
      "-fp-speculation=safe -qopt-report=2 -qopt-report-phase=all -fast-transcendentals -ftz "
      "-finline-functions -finline-limit=1500 -Winline "
      "-diag-disable=7713 "   # disable unused statement functions warning
      "-diag-disable=11021 "  # disable unresolved libraries in ipo warning
      "-diag-disable=10397" ) # disable message reporting location of opt report files
  else()
    string(CONCAT IFS_Fortran_FLAGS "${IFS_Fortran_FLAGS} -fpe0 -convert big_endian -align array64byte "
      "-assume byterecl,noold_maxminloc,realloc_lhs "
      "-traceback -fpic -qopenmp -qopenmp-threadprivate compat -fp-model precise "
      "-fp-speculation=safe -qopt-report=2 -ftz "
      "-finline-functions  "
      "-diag-disable=7713 "   # disable unused statement functions warning
      "-diag-disable=11021 "  # disable unresolved libraries in ipo warning
      "-diag-disable=10397" ) # disable message reporting location of opt report files
  endif()

  if(IFS_CHECK_BOUNDS)
    set(IFS_Fortran_FLAGS "${IFS_Fortran_FLAGS} -check bounds")

    # Disable bounds checking on files that would fail the check.
    # Only files generating false positives should be added to this list.
    set_source_files_properties(
      ${no_bounds_checking}
      PROPERTIES COMPILE_FLAGS "-check nobounds")
  endif()

  if(IFS_INIT_SNAN)
    set(IFS_Fortran_FLAGS "${IFS_Fortran_FLAGS} -init=snan")
  endif()

  set(IFS_Fortran_FLAGS_BIT "-g -O2 -march=core-avx2 -no-fma ")
  set(IFS_Fortran_FLAGS_DEBUG "-g -O0")

  set(IFS_C_FLAGS_BIT "-g -O2 -march=core-avx2 -no-fma ")
  set(IFS_C_FLAGS_DEBUG "-g -O0")

  set(IFS_CXX_FLAGS_BIT "-g -O2 -march=core-avx2 -no-fma ")
  set(IFS_CXX_FLAGS_DEBUG "-g -O0")

  set(ECBUILD_SHARED_LINKER_FLAGS "${IFS_LINK_FLAGS} -Wl,--eh-frame-hdr -Wl,--disable-new-dtags -qopenmp -O2 -L$ENV{TBBROOT}/lib/intel64/gcc4.8 -ltbbmalloc_proxy")
  if( $ENV{HOSTNAME} MATCHES "cc.*")
    # On Cray using the Intel compiler cmake links with the Intel ifcore run-time library without
    # multi-threaded support. This adds the version with multi-threaded support. See JIRA issue ECBUILD-464.
    set(ECBUILD_SHARED_LINKER_FLAGS "${ECBUILD_SHARED_LINKER_FLAGS} -lifcoremt" )
  endif()
  set(ECBUILD_EXE_LINKER_FLAGS "${ECBUILD_SHARED_LINKER_FLAGS}")

  # Use heap-arrays on leap42 desktops which have small default stack size
  if( $ENV{ECPLATFORM} MATCHES "desktop-leap42")
    set(IFS_Fortran_FLAGS "${IFS_Fortran_FLAGS} -heap-arrays 64 ")
  endif()

  # Workaround for Intel 18.0.1 bug
  string(TOUPPER ${CMAKE_BUILD_TYPE} btype)
  string(REPLACE "-qopenmp " "" flags "${IFS_Fortran_FLAGS} ${IFS_Fortran_FLAGS_${btype}}")
  set_source_files_properties(arpifs/fullpos/sufpdyn.F90 PROPERTIES OVERRIDE_COMPILE_FLAGS "${flags}")

  # Source file specific optimisations
  set_source_files_properties(arpifs/phys_ec/cloudsc.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O3 -march=core-avx2 -no-fma ")

  set_source_files_properties(arpifs/adiab/laitli.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O3 -march=core-avx2 -no-fma -qopt-prefetch=0 ")

  set_source_files_properties(arpifs/phys_ec/cloudvar.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O3 -march=core-avx2 -no-fma ")

  set_source_files_properties(radiation/module/radiation_mcica_sw.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O2 -march=core-avx2 -no-fma -vecabi=cmdtarget ")

  set_source_files_properties(radiation/module/radiation_cloud_generator.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O2 -march=core-avx2 -no-fma -vecabi=cmdtarget ")

  set_source_files_properties(arpifs/phys_radi/radintg.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O3 -march=core-avx2 -no-fma ")

  set_source_files_properties(arpifs/phys_ec/radiation_aerosol_optics.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O3 -march=core-avx2 -no-fma ")

  set_source_files_properties(arpifs/adiab/larmes.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O3 -march=core-avx2 -no-fma ")

  set_source_files_properties(ifs/adiab/larmes_xyz.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O3 -march=core-avx2 -no-fma ")


  # Workaround for Intel 18.0.0 bug IFS-1425
  set_source_files_properties(satrad/programs/calc_radiance_fields.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -O1 -g ")

  # Workaround for Intel 2021.4.0 bug in 49R1 v1
  set_source_files_properties(arpifs/mwave/mwave_tovscv.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -O0 -g ")

  # Workaround for ifx 2024.0.0 ICE
  set_source_files_properties(arpifs/control/cprep4.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -O0 -g ")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")

#  # Global flags
  set(IFS_Fortran_FLAGS "${IFS_Fortran_FLAGS} -Mbyteswapio -Kieee -mp")

#  set(IFS_Fortran_FLAGS_BIT "-g -O2 -Mvect=simd:256 -Mstack_arrays -tp=zen2")
  set(IFS_Fortran_FLAGS_BIT "-g -O3 -fast -Mstack_arrays")
  set(IFS_Fortran_FLAGS_DEBUG "-O0 -g  -Mchkstk -Ktrap=fp -Mcoff -Mdwarf1 -Mdwarf2 -Mdwarf3 -Melf  -traceback ")

  if(IFS_CHECK_BOUNDS)
    set(IFS_Fortran_FLAGS "${IFS_Fortran_FLAGS} -Mbounds")

    # Disable bounds checking on files that would fail the check.
    # Only files generating false positives should be added to this list.
    set_source_files_properties(
      ${no_bounds_checking}
      arpifs/adiab/larmes_xyz.F90
      arpifs/phys_radi/radintg.F90
      arpifs/phys_ec/turbulence_layer.F90
      arpifs/phys_ec/surftstp_layer.F90
      arpifs/phys_ec/postphy_layer.F90
      arpifs/control/scan2m.F90
      arpifs/adiab/cpg_drv.F90
      PROPERTIES COMPILE_FLAGS "-Mnobounds")
  endif()


#  # Single file flags
  set_source_files_properties(arpifs/phys_radi/surrtab.F90 PROPERTIES COMPILE_FLAGS "-Kieee -Mbyteswapio")

  set_source_files_properties(arpifs/module/gom_plus.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O0 ")
  set_source_files_properties(arpifs/control/cfcsens2obs.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O0 ")

  set_source_files_properties(arpifs/setup/sunhbmat.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros")

   set_source_files_properties(arpifs/setup/suvfe_basis.F90
	   PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros")


# next ones triggered by compo forecast
   set_source_files_properties(arpifs/module/spectral_arp_mod.F90 arpifs/chem/tm5_chem_ini.F90
                               arpifs/fullpos/hpos.F90 arpifs/phys_radi/uvradi.F90
			       arpifs/chem/tm5_wetchem.F90 arpifs/chem/chem_scav.F90
			       arpifs/phys_ec/aerc_scav.F90 arpifs/phys_radi/radvis.F90
	   PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1  ")
   set_source_files_properties(arpifs/phys_ec/aer_wind.F90 arpifs/phys_ec/aer_src.F90
	                       arpifs/phys_ec/aer_rad.F90 arpifs/phys_ec/aer_phy3.F90
			       arpifs/phys_ec/aer_scavin.F90 arpifs/phys_radi/radiation_scheme.F90
			       arpifs/phys_ec/aerc_scavl19.F90 arpifs/phys_ec/aer_phy2.F90
			       arpifs/phys_ec/aer_diagglomap.F90
			       arpifs/chem/bascoetm5_chem_ini.F90 arpifs/chem/bascoe_gliq.F90
			       arpifs/chem/mocage_nebj.F90 arpifs/chem/mocage_hana.F90
			       arpifs/chem/mocage_hsed.F90 arpifs/chem/mocage_soufreaq.F90
			       arpifs/chem/mocage_phot.F90 arpifs/chem/mo_sad_strat_calc.F90
			       arpifs/chem/mo_photo.F90 arpifs/chem/mo_imp_lu_fac.F90
			       arpifs/chem/mo_strat_aero_settling.F90
			       arpifs/glomap_mode/ukca_volume_mode.F90 arpifs/glomap_mode/aer_dust.F90
			       arpifs/glomap_mode/ukca_aero_step.F90
			       arpifs/glomap_mode/ukca_aero_step_ifs.F90
			       arpifs/glomap_mode/ukca_aero_diag_ifs.F90
			       arpifs/glomap_mode/ukca_chem_ifs.F90
			       arpifs/glomap_mode/ukca_rainout.F90
			       arpifs/glomap_mode/ukca_vgrav.F90
			       arpifs/glomap_mode/ukca_sedim.F90
			       arpifs/glomap_mode/ukca_coag_coff_v.F90
			       arpifs/glomap_mode/ukca_cond_coff_v.F90
			       arpifs/glomap_mode/ukca_conden.F90
			       arpifs/glomap_mode/ukca_solvecoagnucl_v.F90
			       arpifs/glomap_mode/ukca_coagwithnucl.F90
			       arpifs/glomap_mode/ukca_water_content_v.F90
			       surf/external/surfexcdriverstl.F90
			       surf/module/surfexcdriverstl_ctl_mod.F90
	   PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1  ")


   set_source_files_properties(arpifs/pp_obs/ppsta.F90
	   PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros")

  set_source_files_properties(arpifs/pp_obs/ppt_old.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1  -Mstack_arrays ")

  set_source_files_properties(arpifs/phys_radi/ice_effective_radius.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros -Mstack_arrays ")

  set_source_files_properties(radiation/module/radiation_cloud_optics.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros  ")

  set_source_files_properties(radiation/module/radiation_mcica_lw.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros -Mstack_arrays ")

   set_source_files_properties(radiation/module/radiation_two_stream.F90
	   PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros -Mstack_arrays ")

  set_source_files_properties(arpifs/phys_ec/vdfhghtn.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros -Mstack_arrays ")

  set_source_files_properties(arpifs/phys_ec/cubasen.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros -Mstack_arrays ")

  set_source_files_properties(arpifs/phys_ec/cudtdqn.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros ")

  set_source_files_properties(arpifs/phys_ec/cududv.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mstack_arrays  ")

  set_source_files_properties(arpifs/phys_ec/cumastrn.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mstack_arrays ")

  set_source_files_properties(arpifs/phys_ec/cuancape2.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect ")

  set_source_files_properties(arpifs/phys_ec/cloud_satadj.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros -Mstack_arrays ")

  set_source_files_properties(arpifs/phys_ec/cloudsc.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros ")

  set_source_files_properties(arpifs/phys_ec/cuddrafn.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros ")

  set_source_files_properties(arpifs/phys_ec/gwsetup.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros ")

  set_source_files_properties(arpifs/pp_obs/ppitpq.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mflushz -Mno-signed-zeros -Mstack_arrays ")

  set_source_files_properties(
          arpifs/pp_obs/ppgeop.F90 arpifs/pp_obs/pppmer.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mflushz -Mno-signed-zeros -Mstack_arrays ")

  set_source_files_properties(surf/external/surfexcdriver.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O2 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros ")

  set_source_files_properties(surf/module/surfexcdriver_ctl_mod.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros ")

  set_source_files_properties(surf/module/vexcs_mod.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros -Mstack_arrays ")

  set_source_files_properties(surf/module/srfwdif_mod.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mflushz -Mno-signed-zeros -Mstack_arrays ")

  set_source_files_properties(surf/module/sppcfl_mod.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros -Mstack_arrays ")

  set_source_files_properties(surf/module/vupdz0_mod.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -Mvect=gather,idiom,tile -Mflushz -Mno-signed-zeros -Mstack_arrays ")

  set_source_files_properties(surf/module/srfcotwo_mod.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1  -Mstack_arrays ")

  set_source_files_properties(ifsobs/src/dbase/ifsobs_dbase_mod.F90
	  PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O0 ")

  set_source_files_properties(arpifs/op_obs/hop.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O0 -Wlopt,'-passes=default<O0>' -Wllc,-O0 ")
  set_source_files_properties(arpifs/op_obs/taskobtl.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O0 -Wlopt,'-passes=default<O0>' -Wllc,-O0 ")
  set_source_files_properties(arpifs/op_obs/taskobtl_thread.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O0 -Wlopt,'-passes=default<O0>' -Wllc,-O0 ")
#
# (24.1) FPEs
  set_source_files_properties(arpifs/phys_ec/gwsetup.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -tp=znver1 -Wlopt,'-passes=default<O2>' -Wllc,-O1 ")
  set_source_files_properties(surf/module/surfexcdriver_ctl_mod.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -tp=znver1 -Wlopt,'-passes=default<O2>' -Wllc,-O1 ")
  set_source_files_properties(surf/module/vsurf_mod.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -tp=znver1 -Wlopt,'-passes=default<O2>' -Wllc,-O1 ")
  set_source_files_properties(arpifs/phys_ec/vdfhghtn.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -tp=znver1 -Wlopt,'-passes=default<O2>' -Wllc,-O1 ")
  set_source_files_properties(surf/module/voskin_mod.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -tp=znver1 -Wlopt,'-passes=default<O2>' -Wllc,-O1 ")
  set_source_files_properties(surf/module/srft_mod.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -tp=znver1 -Wlopt,'-passes=default<O2>' -Wllc,-O1 ")
  set_source_files_properties(surf/module/vexcs_mod.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -tp=znver1 -Wlopt,'-passes=default<O2>' -Wllc,-O1 ")
  set_source_files_properties(arpifs/phys_ec/cloud_satadj.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -tp=znver1 -Wlopt,'-passes=default<O2>' -Wllc,-O1 ")
  set_source_files_properties(arpifs/phys_ec/cloudsc.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -tp=znver1 -Mstack_arrays -Wlopt,'-passes=default<O2>' -Wllc,-O1 ")
  set_source_files_properties(arpifs/phys_ec/updecaec.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 -tp=znver1 -Wlopt,'-passes=default<O2>' -Wllc,-O1 ")

# (24.1) get the adjoint test to pass
  set_source_files_properties( arpifs/phys_ec/cuascn2ad.F90
                               arpifs/phys_ec/cuascn2tl.F90
                               arpifs/phys_ec/cuascn2.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 ")
  set_source_files_properties( arpifs/phys_ec/cuaflx2ad.F90
                               arpifs/phys_ec/cuaflx2tl.F90
                               arpifs/phys_ec/cuaflx2.F90
    PROPERTIES OVERRIDE_COMPILE_FLAGS "${IFS_Fortran_FLAGS} -g -O1 ")


  set(IFS_C_FLAGS_BIT "-g -O1")
  set(IFS_C_FLAGS_DEBUG "-g -O0")
#
#  set(ECBUILD_SHARED_LINKER_FLAGS "-Wl,--eh-frame-hdr -Wl,--disable-new-dtags -openmp")
#  set(ECBUILD_EXE_LINKER_FLAGS "${ECBUILD_SHARED_LINKER_FLAGS}")


endif()

if( CMAKE_C_COMPILER_ID MATCHES "Clang")
  set(CMAKE_C_LINK_FLAGS "") # no openmp possible
endif()
if( CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_LINK_FLAGS "") # no openmp possible
endif()

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(ECBUILD_EXE_LINKER_FLAGS "${ECBUILD_EXE_LINKER_FLAGS} -Wl,-stack_size,12000000") # 1000000 = 16MB
endif()

macro( add_precision_suffix _var )
  if( DEFINED ${_var} )
    string( REPLACE "IFS" "${PNAME}" _var_with_suffix "${_var}" )
    set( ${_var_with_suffix} ${${_var}} )
  endif()
endmacro()

add_precision_suffix( IFS_Fortran_FLAGS )
add_precision_suffix( IFS_Fortran_FLAGS_DEBUG )
add_precision_suffix( IFS_Fortran_FLAGS_BIT )

add_precision_suffix( IFS_C_FLAGS )
add_precision_suffix( IFS_C_FLAGS_DEBUG )
add_precision_suffix( IFS_C_FLAGS_BIT )

add_precision_suffix( IFS_CXX_FLAGS )
add_precision_suffix( IFS_CXX_FLAGS_DEBUG )
add_precision_suffix( IFS_CXX_FLAGS_BIT )
