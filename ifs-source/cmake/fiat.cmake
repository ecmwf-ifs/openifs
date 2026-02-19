# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction


### Add contributed fiat if not already added via bundle
if( NOT TARGET fiat )
  ifs_propagate_flags( fiat )
  add_subdirectory( contrib/fiat )
endif()

### Find fiat

ecbuild_find_package( fiat REQUIRED )
