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

    #Meteo-France routine, which could also be handled as a dummy
    aladin/nudging/nudglhprecip.F90
)

# Static arpifs source exclude list removed for OpenIFS minimisation.

# Selectively add back files that are required for forecast-only or
# OpenIFS-only but have been removed in the folder exclusions above

include(arpifs_fc_include)

list(APPEND arpifs_public_libs fc_only_intfb)

# Some #include dependencies need to be satisfied,
# even though the actual satrad routine is replaced with a dummy.
list(APPEND arpifs_private_includes satrad/interface openifs/emos)


include(arpifs_oifs_include)

# Needs to be a public libn so that master can access openifs_intfb 
# which permits the dummies inclusion in master.F90
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
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/etrans/interface>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mse/interface>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mse/externals>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/biper/interface>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpa/conv/interface>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpa/micro/interface>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpa/turb/interface>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpa/chem/interface>

  PRIVATE_INCLUDES
    arpifs/namelist
    arpifs/ald_inc/namelist
    arpifs/ald_inc/function
    arpifs/var
    blacklist/include
    ${arpifs_private_includes}

  PUBLIC_LIBS arpifs_intfb surf.${PREC} trans.${PREC}
    ${arpifs_public_libs}
    algor.${PREC} ${ECWAM_LIBRARIES} ${IFSAUX_LIBRARIES} fckit
    field_api_${prec}
    ${ECCODES_LIBRARIES} ${ATLAS_LIBRARIES}
    ${MULTIO_LIBRARIES} ${FDB_LIBRARIES}
    ${NEMOVAR_LIBRARIES}
    NetCDF::NetCDF_Fortran # [IFS-HHH] for radiation/module/easy_netcdf.F90
    ecflow_lightf

  PRIVATE_LIBS
    ${arpifs_private_libs}
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
  LIBS arpifs.${PREC} ${BLACKLIST_LIBRARIES}
  LINKER_LANGUAGE Fortran
  CONDITION HAVE_MPI
 )

