# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

if(HAVE_NEMO)
  # Make sure NEMO's precision matches the IFS, unless overridden by specifying
  # NEMO_[SP/DP]_ENABLE_SINGLE_PRECISION
  if( NOT DEFINED NEMO_${PREC}_ENABLE_SINGLE_PRECISION )
    set( NEMO_${PREC}_ENABLE_SINGLE_PRECISION ${HAVE_SINGLE_PRECISION} )
  endif()

  add_subdirectory(nemo)

  set( NEMOVAR_LIBRARIES ${NEMO_${PREC}_LIBRARIES} )

  # Finally, if we are using single precision for NEMO, add the
  # PARKIND1_SINGLE_NEMO preprocessor definition to set the ocean working
  # precision within the IFS coupling routines (JPRO)
  if( NEMO_${PREC}_HAVE_SINGLE_PRECISION )
    list(APPEND IFS_DEFINITIONS PARKIND1_SINGLE_NEMO)
  endif()
endif()
