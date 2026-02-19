# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

ecbuild_info("[algor]")

ecbuild_add_library(TARGET algor.${PREC}
  LINKER_LANGUAGE Fortran
  DEFINITIONS ${IFS_DEFINITIONS}
  SOURCES_GLOB algor/*
  SOURCES_EXCLUDE_REGEX algor/external/minim/* algor/internal/minim/*
  PUBLIC_INCLUDES $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/algor/interface>
  PUBLIC_LIBS  ${IFSAUX_LIBRARIES}
  PRIVATE_LIBS ${LAPACK_LIBRARIES})
