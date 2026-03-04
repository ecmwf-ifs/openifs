# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

if( NOT TARGET mpi_serial )
    ecbuild_add_library(TARGET mpi_serial
      SOURCES_GLOB mpi_serial/source/*
      TYPE STATIC)

    if(NOT HAVE_MPI)
      target_include_directories(mpi_serial
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpi_serial/source>)
    endif()

    if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      set(_whole_archive "-Wl,-force_load")
      set(_no_whole_archive "")
    else()
      set(_whole_archive "-Wl,--whole-archive")
      set(_no_whole_archive "-Wl,--no-whole-archive")
    endif()
    set(MPI_SERIAL_LIBRARIES ${_whole_archive} mpi_serial ${_no_whole_archive})
    set(MPI_SERIAL_LIBRARIES ${MPI_SERIAL_LIBRARIES} PARENT_SCOPE)
endif()