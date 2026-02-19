# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

ecbuild_add_library(
    TARGET  dwarf_sladv.${PREC}
    SOURCES
        dwarf_sladv_module.F90
        dwarf_sladv_kind_module.F90
        dwarf_sladv_parameters_module.F90
        dwarf_sladv_use_module.F90
        dwarf_sladv_geometry_module.F90
        dwarf_sladv_compute_dp_noniter.F90
        dwarf_sladv_compute_dp_iter.F90
        dwarf_sladv_departurepoint_method_module.F90
        dwarf_sladv_dp_cart2latlon.F90
        dwarf_sladv_dp_rotmat_trans.F90
        dwarf_sladv.h
    PUBLIC_LIBS
        fckit
        atlas_f
        dwarf_util
        ${IFS_OMP_Fortran_LIBRARIES}
)

if( HAVE_SINGLE_PRECISION )
  target_compile_definitions( dwarf_sladv.${PREC} PRIVATE "-DHAVE_SINGLE_PRECISION" )
endif()

target_compile_definitions_FILENAME( dwarf_sladv.${PREC} )
