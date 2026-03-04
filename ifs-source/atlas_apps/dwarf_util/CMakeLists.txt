# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

if( NOT TARGET dwarf_util )
    ecbuild_add_library(
        TARGET dwarf_util
        SOURCES
            dwarf_module.F90
    )
    # Make sure that the module directory is precision indendent
    include(target_fortran_module_directory)
    target_fortran_module_directory( TARGET dwarf_util MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/atlas_apps )
endif()
