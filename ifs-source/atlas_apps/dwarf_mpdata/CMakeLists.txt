# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

ecbuild_add_library(
    TARGET  dwarf_mpdata.${PREC}
    SOURCES
        advection_MPDATA_auxiliary_module.F90
        advection_MPDATA_geometry_module.F90
        advection_MPDATA_limiter_module.F90
        advection_MPDATA_mappings_module.F90
        advection_MPDATA_MPDATA_module.F90
        advection_MPDATA_nabla_module.F90
        advection_MPDATA_topology_module.F90
        dwarf_mpdata_module.F90
    PUBLIC_LIBS
        fckit
        atlas_f
        dwarf_util
        ${IFS_OMP_Fortran_LIBRARIES}
)

if( HAVE_SINGLE_PRECISION )
  target_compile_definitions( dwarf_mpdata.${PREC} PRIVATE "-DHAVE_SINGLE_PRECISION" )
endif()
