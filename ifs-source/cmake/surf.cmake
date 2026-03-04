# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

ecbuild_info("[surf]")

ecbuild_list_add_pattern(LIST surf.${PREC}_src GLOB
  surf/module/* 
  surf/external/*
)


ecbuild_add_library(TARGET surf.${PREC}
  DEFINITIONS ${IFS_DEFINITIONS}
  #SOURCES_GLOB surf/module/* surf/external/*
  SOURCES ${surf.${PREC}_src} 
  PUBLIC_INCLUDES $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/surf/interface>
  PRIVATE_INCLUDES surf/function 
  LIBS ${IFSAUX_LIBRARIES})

# ----------------------------------------------------------------------------

if( NOT TARGET surf_offline_intfb )
ecbuild_generate_fortran_interfaces(TARGET surf_offline_intfb
  DIRECTORIES driver
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/surf/offline
  DESTINATION surf_offline
  INCLUDE_DIRS surf_offline_intfb_includes
  PARALLEL ${FCM_PARALLEL})
endif()

ecbuild_add_library(TARGET cmflood.${PREC}
  DEFINITIONS ${IFS_DEFINITIONS}
  SOURCES_GLOB surf/cmflood/*
  PRIVATE_INCLUDES ${NETCDF_INCLUDE_DIRS}
  LIBS ${IFSAUX_LIBRARIES} ${NETCDF_LIBRARIES})

if( NOT TARGET cmfMASTER )
ecbuild_add_executable(TARGET cmfMASTER
  SOURCES surf/offline/cmfld1s.F90
  INCLUDES ${NETCDF_INCLUDE_DIRS}
  LIBS cmflood.${PREC} ${NETCDF_LIBRARIES}
  LINKER_LANGUAGE Fortran)
endif()


# Generate OSM's Fortran modules in a separate directory to avoid clashes with arpifs/ifsaux modules.
get_directory_property(include_directories INCLUDE_DIRECTORIES)
list(REMOVE_ITEM include_directories "${CMAKE_Fortran_MODULE_DIRECTORY}")
set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "${include_directories}")
set(CMAKE_Fortran_MODULE_DIRECTORY_tmp "${CMAKE_Fortran_MODULE_DIRECTORY}")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/module/osm${PROJECT_PRECISION_SUFFIX}")

if( NOT TARGET osmMASTER )
ecbuild_add_executable(TARGET osmMASTER
  DEFINITIONS ${IFS_DEFINITIONS}
  SOURCES surf/offline/master1s.F90
  SOURCES_GLOB surf/offline/driver/* surf/cmflood/*
  INCLUDES
    ${CMAKE_Fortran_MODULE_DIRECTORY}
    surf/offline/function
    surf/offline/namelist
    ${NETCDF_INCLUDE_DIRS}
  LIBS surf_offline_intfb surf.${PREC} ${IFSAUX_LIBRARIES} ${NETCDF_LIBRARIES})
endif()

# Restore the default Fortran module directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY_tmp})
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# ----------------------------------------------------------------------------

if( NOT TARGET surf_offline_util_objs )
ecbuild_add_library(TARGET surf_offline_util_objs
  SOURCES surf/offline/util/lib_dates.F90
          surf/offline/util/get_points_mod.F90
  PRIVATE_INCLUDES ${NETCDF_INCLUDE_DIRS}
  TYPE OBJECT)

foreach(program IN ITEMS

    create_grid_info
    convNetcdf2Grib
    create_init_clim
    adjust_forc
    caldtdz
    conv_forcing
    find_points)

  ecbuild_add_executable(TARGET ${program}
    SOURCES surf/offline/util/${program}.F90
    OBJECTS surf_offline_util_objs
    INCLUDES ${NETCDF_INCLUDE_DIRS} ${ECCODES_INCLUDE_DIRS}
    LIBS ${NETCDF_LIBRARIES} ${ECCODES_LIBRARIES})

endforeach()
endif()
