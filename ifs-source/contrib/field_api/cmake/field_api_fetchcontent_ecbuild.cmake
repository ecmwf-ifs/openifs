# (C) Copyright 2022- ECMWF.
# (C) Copyright 2022- Meteo-France.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

include(FetchContent)
FetchContent_Populate(
    ecbuild
    URL            https://github.com/ecmwf/ecbuild/archive/refs/tags/3.7.0.tar.gz
    SOURCE_DIR     ${CMAKE_BINARY_DIR}/ecbuild
    BINARY_DIR     ${CMAKE_BINARY_DIR}/_deps/ecbuild-build
    SUBBUILD_DIR   ${CMAKE_BINARY_DIR}/_deps/ecbuild-subbuild
  )
find_package( ecbuild 3.7 REQUIRED HINTS ${CMAKE_BINARY_DIR} )
