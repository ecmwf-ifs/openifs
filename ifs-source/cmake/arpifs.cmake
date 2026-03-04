# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

ecbuild_info("[arpifs]")

if( NOT HAVE_MGRIDS )
  list( APPEND arpifs_exclude arpifs/mgrids/* )
endif()

ecbuild_list_add_pattern(LIST arpifs.${PREC}_src GLOB

    arpifs/*

    # FIXME: circular dependency between arpifs & algor
    algor/internal/minim/*
    algor/external/minim/*

    # FIXME[IFS-DDD]: circular dependency between arpifs & satrad
    ${satrad_ifs_srcs}

    # FIXME[IFS-HHH]: circular dependency between arpifs & radiation
    radiation/module/*

    # Circular dependency sources from wam. Defined in wam.cmake
    ${wam_ifs_srcs}
    ${wamassi_ifs_srcs}
)

ecbuild_list_exclude_pattern(LIST arpifs.${PREC}_src REGEX

    arpifs/programs/*

    ${arpifs_exclude}

    # Not used, contain undefined references to {push,pop}{integer,real}array_
    arpifs/op_obs/vertdisc_ad.F90
    arpifs/op_obs/kernel_pbp_ad.F90
    arpifs/op_obs/pushreal8.F90
    arpifs/op_obs/popreal8.F90
    arpifs/op_obs/pushinteger4.F90
    arpifs/op_obs/popinteger4.F90
    arpifs/op_obs/popboolean.F90

    # FIXME: remove scat dependency on arpifs
    arpifs/module/yomersca.F90

    # FIXME[IFS-DDD]: circular dependency between arpifs & satrad
    arpifs/module/yommwave.F90

    # Included in grib_mean.x
    arpifs/utility/grib_mean.f90
    arpifs/utility/link.f90
)

list(APPEND arpifs.${PREC}_src
    openifs/dummy_ifsobs/dbase_mod.F90
)
ecbuild_list_add_pattern(LIST arpifs.${PREC}_src GLOB
# need to add the local copy of emos libs so that arpifs can build
# with both forecast-only and openifs-only
  openifs/emos/common/*
)

# Some #include dependencies need to be satisfied,
# even though the actual satrad routine is replaced with a dummy.
list(APPEND arpifs_private_includes satrad/interface openifs/emos)


list(APPEND arpifs_public_libs wam.${PREC})

# Add the openifs "smart" dummies, which are built in arpifs, rather than 
  # dummy to ensure consistent generation of fortran interface blocks, when 
  # compared to the full build. 
  ecbuild_list_add_pattern(LIST arpifs.${PREC}_src GLOB
    openifs/dummy/*
  )

  list(APPEND arpifs_public_libs openifs_intfb) 

# Intel 18.* has problems compiling arpifs/oops/fields_io_mod, which is only used by OOPS.
# OOPS not being tested with Intel 18, we exclude the file for this compiler major version
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  if( CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 19)
    ecbuild_list_exclude_pattern(LIST ifs.${PREC}_src REGEX
      arpifs/oops/fields_io_mod*
      arpifs/control/cprep4.F90
    )
  endif()
endif()

ecbuild_add_library(
  TARGET  arpifs.${PREC}
  SOURCES ${arpifs.${PREC}_src}

  DEFINITIONS ${IFS_DEFINITIONS}

  PUBLIC_INCLUDES
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/arpifs/common>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/arpifs/function>

  PRIVATE_INCLUDES
    arpifs/namelist
    arpifs/ald_inc/namelist
    arpifs/ald_inc/interface
    arpifs/ald_inc/function
    arpifs/var
    ${arpifs_private_includes}

  PUBLIC_LIBS arpifs_intfb surf.${PREC} trans.${PREC}
    ${arpifs_public_libs}
    algor.${PREC} ${IFSAUX_LIBRARIES} fckit
    ${ECCODES_LIBRARIES} ${ATLAS_LIBRARIES}
    ${MULTIO_LIBRARIES} ${FDB_LIBRARIES}
    ${NEMOVAR_LIBRARIES}
    NetCDF::NetCDF_Fortran # [IFS-HHH] for radiation/module/easy_netcdf.F90

  PRIVATE_LIBS
    ${arpifs_private_libs}
    wam_intfb
    ${LAPACK_LIBRARIES}
)

if( HAVE_MGRIDS )
  target_link_libraries( arpifs.${PREC} PUBLIC dwarf_mpdata.${PREC} dwarf_sladv.${PREC} )
endif()

fckit_target_preprocess_fypp( arpifs.${PREC}
  FYPP_ARGS -m os -m field_config -M ${CMAKE_CURRENT_SOURCE_DIR}/arpifs/scripts
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/arpifs/module/field_config.yaml
    ${CMAKE_CURRENT_SOURCE_DIR}/arpifs/module/surface_fields_config.yaml
)

if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray" AND CMAKE_GENERATOR STREQUAL "Ninja")
  add_custom_command(TARGET arpifs.${PREC} PRE_LINK
    COMMAND ${CMAKE_COMMAND} -D filename="${CMAKE_BINARY_DIR}/CMakeFiles/arpifs.${PREC}.rsp"
      -P ${PROJECT_SOURCE_DIR}/cmake/patch_arpifs_rsp.cmake
    COMMENT "Patching CMakeFiles/arpifs.${PREC}.rsp")
endif()

ecbuild_add_executable( TARGET ifsMASTER.${PREC}
  DEFINITIONS ${IFS_DEFINITIONS}
  SOURCES arpifs/programs/master.F90
  INCLUDES ${FCKIT_INCLUDE_DIRS}
  LIBS arpifs.${PREC} wam.${PREC}
  LINKER_LANGUAGE Fortran
  CONDITION HAVE_MPI
 )


if( NOT TARGET grib_mean.x )
  ecbuild_add_executable(TARGET grib_mean.x
    SOURCES arpifs/utility/grib_mean.f90 arpifs/utility/link.f90
    INCLUDES ${ECCODES_INCLUDE_DIRS}
    LIBS ${IFSAUX_LIBRARIES} ${ECCODES_LIBRARIES})
endif()
