# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction
# 
# (C) Copyright 1989- Meteo-France.
# 


ecbuild_info( "build type         : [${CMAKE_BUILD_TYPE}]" )
ecbuild_info( "precision          : [${PREC}]")
ecbuild_info( "LAPACK_LIBRARIES   : [${LAPACK_LIBRARIES}]")
foreach( lang Fortran C CXX )
  ecbuild_info( "${lang} -- ${CMAKE_${lang}_COMPILER_ID} ${CMAKE_${lang}_COMPILER_VERSION}"  )
  ecbuild_info( "    compiler     : ${CMAKE_${lang}_COMPILER}" )
  ecbuild_info( "    flags        : ${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${CMAKE_BUILD_TYPE_CAPS}} ${${PNAME}_${lang}_FLAGS} ${${PNAME}_${lang}_FLAGS_${CMAKE_BUILD_TYPE_CAPS}}" )
  ecbuild_info( "    OpenMP flags : ${OpenMP_${lang}_FLAGS}" )
  ecbuild_info( "    link flags   : ${CMAKE_${lang}_LINK_FLAGS}" )
endforeach()

ecbuild_info( "link flags" )
ecbuild_info( "    executable [${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${CMAKE_BUILD_TYPE_CAPS}}]" )
ecbuild_info( "    shared lib [${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_${CMAKE_BUILD_TYPE_CAPS}}]" )
ecbuild_info( "    static lib [${CMAKE_MODULE_LINKER_FLAGS} ${CMAKE_MODULE_LINKER_FLAGS_${CMAKE_BUILD_TYPE_CAPS}}]" )