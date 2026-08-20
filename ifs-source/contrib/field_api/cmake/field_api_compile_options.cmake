# (C) Copyright 2022- ECMWF.
# (C) Copyright 2022- Meteo-France.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI|NVIDIA|NVHPC)
  ecbuild_add_fortran_flags("-Mlarge_arrays")
  ecbuild_add_fortran_flags("-gopt")

  ecbuild_add_fortran_flags("-Minfo=accel,all,ccff" BUILD DEBUG)

# These are architecture/compiler/offload-library specific options 
# that should really be coming from external input
# set(CMAKE_Fortran_FLAGS "-gpu=cc70")
endif ()
