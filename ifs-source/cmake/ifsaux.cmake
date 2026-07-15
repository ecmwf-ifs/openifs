# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

ecbuild_info("[ifsaux]")

ecbuild_add_library( TARGET ifsaux.${PREC}
  DEFINITIONS ${IFS_DEFINITIONS}
  LINKER_LANGUAGE Fortran

  SOURCES_GLOB
    ifsaux/module/*
    ifsaux/fi_libc/*
    ifsaux/ddh/*
    ifsaux/support/*
    ifsaux/utilities/*

  PUBLIC_INCLUDES
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ifsaux/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ifsaux/fa>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ifsaux/ddh>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ifsaux/misc>
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    ${MPI_Fortran_INCLUDE_PATH}

  PRIVATE_INCLUDES
    ifsaux/fi_libc
    ifsaux/lfi
    ifsaux/lfi_alt

  PUBLIC_LIBS
    ${IFS_OMP_Fortran_LIBRARIES}
    ${ECCODES_LIBRARIES}
    fiat parkind_${prec}
    ifsdummy

  PRIVATE_LIBS
    ${MULTIO_LIBRARIES} ${FDB_LIBRARIES}
)

foreach(target ddhrun)
  if( NOT TARGET ${target} )
    ecbuild_add_executable(TARGET ${target} SOURCES ifsaux/programs/${target}.F LIBS ifsaux.${PREC})
    set_property(TARGET ${target} PROPERTY SERIAL TRUE)
  endif()
endforeach()

set( IFSAUX_LIBRARIES ifsaux.${PREC} )
