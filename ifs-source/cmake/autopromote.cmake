# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

# Auto-promote REALs without explicit kind declaration to 64 bit
#
# A file "autopromoted_files" will be written in the build-dir with files that have been autopromoted.
# Two options are available, toggled with ENABLE_MANUAL_AUTOPROMOTE=ON|OFF (ON=default)
#  1) OFF : All files are parsed for detection of REALS without explicit kind
#  2) ON  : List has been prepared with files that need to be fixed
#
# The content of "autopromoted_files" should be identical for both options, after sorting and removing duplicates:
#    sort -u autopromoted_files
#
# TODO:
#   - The listed files under "ENABLE_MANUAL_AUTOPROMOTE" should be adapted with explicit REAL kinds
#   - Apply remaining autopromote flags explicitly to wam, scat, ifsaux/ddh, satrad/emos_interpolation sources
#   - Investigate intentionally ignored files for autopromotion, whichg parser triggered, but we refuse to autopromote
#       # Otherwise ifstraj/ifsmin *not* bit-reproducible. Please fully qualify REALs even if it breaks bit-reproducibility.
#          satrad/module/yoamsu.F90
#          arpifs/module/yomclmicst.F90
#          radiation/module/radiation_thermodynamics.F90
#       # Intentionally ingored otherwise EDA *not* bit-reproducible. Please fully qualify REALs even if it breaks bit-reproducibility
#          Spectral_Filter.F90
#          Spread_Skill_Time_Avg.F90
#          Wavelet_Decomp.F90
#       # Still needs porting to fully qualified REALs and should not affect forecast.
#          ^obstat/.*

option( ENABLE_MANUAL_AUTOPROMOTE "Use manual list of files to be autopromoted." ON )

### Collect all files that possibly have to be autopromoted.
set(all_srcs)
foreach(lib ${${PROJECT_NAME}_ALL_LIBS} ${${PROJECT_NAME}_ALL_EXES})

    get_target_property(sources ${lib} SOURCES)

    foreach(src ${sources})

      # Ignore target object generator expressions
      if(src MATCHES "\\\$<TARGET_OBJECTS:")
        continue()
      endif()

      # Ignore non-Fortran sources
      get_source_file_property(lang ${src} LANGUAGE)
      if(NOT lang MATCHES "Fortran")
        continue()
      endif()

      # Skip generated files which may not exist yet
      get_source_file_property(gen ${src} GENERATED)
      if (gen)
        continue()
      endif()

      list(APPEND all_srcs ${src})
  endforeach()
endforeach()

### Fill in list of autopromoted sources
set(autopromote_srcs)
if( ENABLE_MANUAL_AUTOPROMOTE )

  ecbuild_log("Manual created list used for autopromote")
  # 1) Following files were found by parsing all files for REAL without its kind specified.
  #       These should get priority in fixing.
  list(APPEND autopromote_srcs
      aeolus/support/lexer.F90
      aeolus/support/numerics.F90
      arpifs/chem/mo_jshort.F90
      arpifs/module/mo_mod_photo.F90
      ifsaux/module/xrd_unix_env.F90
      ifsaux/utilities/ssort.F
      odb/module/binterpol.F90
      odb/module/datstat.F90
      odb/module/merge_model_info.F90
      odb/module/odb2.F90
      odb/tools/Bufr2hex.F90
      odb/tools/Convert_varbcfile.F90
      odb/tools/Kind.F90
      odb/tools/Load_balancing.F90
      odb/tools/Merge_gmi_swaths.F90
      prepdata/kpp/conv2grib_geo_nemo.F90
      prepdata/ma_tools/ma_replacefields.F90
      prepdata/ma_tools/ma_rhtoq.F90
      prepdata/ma_tools/ma_vinterp.F90
      prepdata/module/asind.F90
      prepdata/module/create_regular_reduced_grid.F90
      prepdata/module/fill_latlon.F90
      prepdata/module/gettrk_modules.F90
      prepdata/module/grb_hdr.F90
      prepdata/module/grib2gridspecs.F90
      prepdata/module/griddef.F90
      prepdata/module/initgrid.F90
      prepdata/module/init_neighbours.F90
      prepdata/module/matutils.F90
      prepdata/module/print_gridspecs.F90
      prepdata/module/sind.F90
      prepdata/module/yomgrib_si.F90
      prepdata/programs/cint.F90
      prepdata/programs/ensms_veps.F90
      prepdata/programs/gen_pert.F90
      prepdata/programs/GH_RH.F90
      prepdata/programs/invlap.F90
      prepdata/programs/mmcalc.F90
      prepdata/programs/prob_perc.F90
      prepdata/programs/prob_thr.F90
      prepdata/programs/reord_veps.F90
      prepdata/programs/signi.F90
      prepdata/programs/smmacalc.F90
      prepdata/programs/wmem.F90
      prepdata/programs/wm.F90
      prepdata/tcyc/max.F90
      prepdata/tcyc/radii.F90
      prepdata/tcyc/storm.F90
      prepdata/tcyc/traj.F90
      satrad/module/mod_mwatlas_m2.F90
      satrad/module/rttov_unix_env.F90
      ssa/util/geo2ocem_ssa.F90
      ssa/util/netcdffield_ssa.F90
      ssa/util/orca_grid_ssa.F90
      surf/offline/util/conv_forcing.F90
      surf/offline/util/create_init_clim.F90
  )

  # 2) Add files that we know need to be promoted.
  #    These are sub-project-wide globs that need to be taken more care of.
  #    Adding autopromote flags can better be done in respective <subproject>.cmake files

  foreach(src ${all_srcs})

    # FIXME: WAM and SCAT are using implicit typing!
    if(src MATCHES "^wam/.*|^scat/.*")
      list(APPEND autopromote_srcs ${src})
      continue()
    endif()

    # FIXME: DDH decoding is using implicit typing!
    if(src MATCHES "^ifsaux/ddh/.*|ifsaux/programs/ddhrun.F")
      list(APPEND autopromote_srcs ${src})
      continue()
    endif()

    # FIXME: Until these routines can be replaced by atlas
    if(src MATCHES "^satrad/emoslib_interpolation/.*")
      list(APPEND autopromote_srcs ${src})
      continue()
    endif()

  endforeach()

else()
  ### Once the files under "1)" above are fixed, it should be safe to remove this else() branch

  ecbuild_log("Parsing all files for autopromote")

  foreach(src ${all_srcs})

    # Intentionally ignored, as it contains
    #    REAL, PRIVATE :: Z_DEFAULT_REAL      ! intentionally not REAL(KIND=JPRB)
    if(src MATCHES "arpifs/module/yommsc.F90")
      continue()
    endif()

    # FIXME: intentionally ingored otherwise ifstraj/ifsmin *not* bit-reproducible
    if(src MATCHES "satrad/module/yoamsu.F90|arpifs/module/yomclmicst.F90|radiation/module/radiation_thermodynamics.F90")
      continue()
    endif()

    # FIXME: intentionally ingored otherwise EDA *not* bit-reproducible
    if(src MATCHES "Spectral_Filter.F90|Spread_Skill_Time_Avg.F90|Wavelet_Decomp.F90")
      continue()
    endif()

    # FIXME: use fully qualified REALs
    if(src MATCHES "^obstat/.*")
      continue()
    endif()

    # FIXME: WAM and SCAT are using implicit typing!
    if(src MATCHES "^wam/.*|^scat/.*")
      list(APPEND autopromote_srcs ${src})
      continue()
    endif()

    # FIXME: DDH decoding is using implicit typing!
    if(src MATCHES "^ifsaux/ddh/.*|ifsaux/programs/ddhrun.F")
      list(APPEND autopromote_srcs ${src})
      continue()
    endif()

    # FIXME: Until these routines can be replaced by atlas
    if(src MATCHES "^satrad/emoslib_interpolation/.*")
      list(APPEND autopromote_srcs ${src})
      continue()
    endif()

    file(STRINGS ${src} real_no_kind REGEX "^[ ]*[Rr][Ee][Aa][Ll]([ ]+[:a-zA-Z]|[ ]*,)" LIMIT_COUNT 1)

    if(real_no_kind)
      list(APPEND autopromote_srcs ${src})
    endif()

  endforeach()
endif()

if( autopromote_srcs )
  file(WRITE "${PROJECT_BINARY_DIR}/autopromoted_files" "")
  list(LENGTH autopromote_srcs autopromote_srcs_size)
  ecbuild_info("WARNING!!! autopromote applied to ${autopromote_srcs_size} sources.")
  foreach( src ${autopromote_srcs} )
    ecbuild_debug(" - ${src}")
    file(APPEND "${PROJECT_BINARY_DIR}/autopromoted_files" "${src}\n")
  endforeach()
endif()

### Apply autopromote_flags to autopromote_srcs
# The autopromote_flags variable is defined in cmake/compile_flags.cmake
set_property(SOURCE ${autopromote_srcs} PROPERTY COMPILE_OPTIONS ${autopromote_flags})
