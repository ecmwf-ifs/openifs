# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

ecbuild_info("[radiation]")


ecbuild_add_library(TARGET radiation.${PREC}
  SOURCES_GLOB radiation/module/*
  PRIVATE_INCLUDES ${NETCDF_INCLUDE_DIRS}
  LIBS arpifs.${PREC} ${IFSAUX_LIBRARIES} NetCDF::NetCDF_Fortran
  CONDITION FALSE) # FIXME[IFS-HHH]: circular dependency between ifs & radiation
