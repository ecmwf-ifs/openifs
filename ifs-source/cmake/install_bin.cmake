# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

ecbuild_info("[install_bin]")

foreach(name IN ITEMS
  fff
  mkabs)

  # Place in build directory but do not install
  file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/install_bin/${name}
                   ${CMAKE_BINARY_DIR}/${name} SYMBOLIC)
endforeach()
