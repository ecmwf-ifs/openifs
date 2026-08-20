# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction


### Add contributed field_api if not already added via bundle
if( NOT TARGET field_api_${prec} )
  ifs_propagate_flags( field_api )
  set(FCKIT_FYPP ${FYPP})
  unset(FYPP)
  set( FIELD_API_${PREC}_ENABLE_SINGLE_PRECISION ${HAVE_SINGLE_PRECISION} )
  add_subdirectory( contrib/field_api )
  set(FYPP ${FCKIT_FYPP})
endif()

### Find field_api

ecbuild_find_package( field_api_${prec} REQUIRED )
