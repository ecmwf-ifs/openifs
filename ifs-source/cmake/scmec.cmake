# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

ecbuild_list_add_pattern(LIST scmec.${PREC}_src GLOB
  scmec/module/* scmec/source/* )

ecbuild_list_exclude_pattern(LIST scmec.${PREC}_src REGEX
  scmec/source/master1c* )

ecbuild_add_library(
  TARGET scmec.${PREC}
  SOURCES ${scmec.${PREC}_src}

  DEFINITIONS ${IFS_DEFINITIONS}

  PRIVATE_INCLUDES 
    ifsaux/include
    arpifs/namelist
    scmec/include

  PUBLIC_INCLUDES 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/arpifs/common>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/arpifs/function>
    $<BUILD_INTERFACE:${include_directories}>

  PUBLIC_LIBS 
      scmec_intfb
      arpifs.${PREC}
      ${IFSAUX_LIBRARIES}
      fiat
)

# Generate Single-Column-specific Fortran modules in a separate directory to avoid clashes with ifs modules.
ifs_target_fortran_module_directory(
    TARGET scmec.${PREC}
    MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/scmec${PROJECT_PRECISION_SUFFIX} )

# ----------------------------------------------------------------------------

ecbuild_add_executable(TARGET MASTER_scm.${PREC}
  DEFINITIONS ${IFS_DEFINITIONS}
  SOURCES scmec/source/master1c.F90
  INCLUDES
    scmec/namelist
    ${FCKIT_INCLUDE_DIRS}
  LIBS scmec.${PREC} NetCDF::NetCDF_Fortran)

# ----------------------------------------------------------------------------

