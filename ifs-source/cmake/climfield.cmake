# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

ecbuild_info("[climfield]")

find_package(ecRegrid)
find_package(TIFF)

if(NOT ECREGRID_FOUND)
  ecbuild_warn("Could not find ecRegrid library - skipping climfield project")
elseif(NOT TIFF_FOUND)
  ecbuild_warn("Could not find TIFF library - skipping climfield project")
else()

if( NOT TARGET climfield_filter2d )
  ecbuild_add_library(TARGET climfield_filter2d
    SOURCES climfield/src/filter2d.F90
    TYPE OBJECT)
endif()

ecbuild_add_library(TARGET climfield.${PREC}
  SOURCES_GLOB climfield/src/[A-Z]*.cpp
               climfield/src/[A-Z]*.h
  OBJECTS climfield_filter2d # for filter_latlon, filter_orog
  PUBLIC_INCLUDES $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/climfield/src>
  PRIVATE_INCLUDES ${ECCODES_INCLUDE_DIRS} ${TIFF_INCLUDE_DIRS} ${ECREGRID_INCLUDE_DIRS}
  LIBS ifsaux.${PREC} ${ECCODES_LIBRARIES} ${TIFF_LIBRARIES} ${ECREGRID_LIBRARIES})

# ----------------------------------------------------------------------------

file(GLOB program_srcs climfield/src/[a-z]*.cpp)

foreach(src IN ITEMS ${program_srcs})

  get_filename_component(program ${src} NAME_WE)

  # The following programs are not needed (also don't compile)
  if(program MATCHES "size_test|read_globe|diff_orog|filter_orog|read_lai_or_albedo|generate_subscale_param|generate_mask|read_modis|magics")
    continue()
  endif()

  if( NOT TARGET ${program} )
    ecbuild_add_executable(TARGET ${program} SOURCES ${src}
      INCLUDES ${ECCODES_INCLUDE_DIRS} ${NETCDF_INCLUDE_DIRS}
      LIBS climfield.${PREC} ${ECCODES_LIBRARIES} ${NETCDF_LIBRARIES})
  endif()

endforeach()

# ----------------------------------------------------------------------------

file(GLOB program_srcs climfield/src/*.F90)
list(FILTER program_srcs EXCLUDE REGEX climfield/src/filter2d.F90) # not a program

foreach(src IN ITEMS ${program_srcs})

  get_filename_component(program ${src} NAME_WE)

  if( NOT TARGET ${program} )
    ecbuild_add_executable(TARGET ${program} SOURCES ${src}
      OBJECTS climfield_filter2d # for filter_globe
      INCLUDES ${ECCODES_INCLUDE_DIRS}
      LIBS ifsaux.${PREC} ${ECCODES_LIBRARIES})
  endif()

endforeach()

# ----------------------------------------------------------------------------

file(GLOB program_srcs climfield/ifs_tools/*.F90)


foreach(src IN ITEMS ${program_srcs})

  get_filename_component(program ${src} NAME_WE)

  if( NOT TARGET ${program} )
    ecbuild_add_executable(TARGET ${program} SOURCES ${src}
      INCLUDES ${ECCODES_INCLUDE_DIRS}
      LIBS trans.${PREC} ifsaux.${PREC} ${ECCODES_LIBRARIES})


    if( program MATCHES "depth_mode_filter" )
      # Disable byte-swapping for depth_mode_filter program.
      # (Note: simply appending -hnobyteswapio via LINK_FLAGS property didn't work!)
      ecbuild_warn("depth_mode_filter")
      add_custom_command(TARGET depth_mode_filter PRE_LINK
        COMMAND ${CMAKE_COMMAND}
          -D TARGET=depth_mode_filter
          -D CMAKE_GENERATOR="${CMAKE_GENERATOR}"
          -D CMAKE_Fortran_COMPILER_ID="${CMAKE_Fortran_COMPILER_ID}"
          -P ${PROJECT_SOURCE_DIR}/cmake/strip_byte_swap_io_flag.cmake
        COMMENT "Stripping byte-swap I/O flag from depth_mode_filter")
    endif()

  endif()

endforeach()

endif()
