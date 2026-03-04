# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction


if( NOT TARGET arpifs_intfb )
  ecbuild_generate_fortran_interfaces(

    TARGET arpifs_intfb

    DIRECTORIES
        
        ${aladin_include}

        arpifs/adiab
        arpifs/c9xx
        arpifs/canari
        arpifs/chem
        arpifs/climate
        arpifs/clradlid
        arpifs/cma2odb
        arpifs/control
        arpifs/dfi
        arpifs/dia
        arpifs/fullpos
        arpifs/gbrad
        arpifs/glomap_mode
        arpifs/interpol
        arpifs/kalman
        arpifs/mwave
        arpifs/nemo
        arpifs/obs_error
        arpifs/obs_preproc
        arpifs/ocean
        arpifs/oops
        arpifs/op_obs
        arpifs/parallel
        arpifs/phys_dmn
        arpifs/phys_ec
        arpifs/phys_radi
        arpifs/pp_obs
        arpifs/raingg
        arpifs/sekf
        arpifs/setup
        arpifs/sinvect
        arpifs/smos
        arpifs/transform
        arpifs/utility

    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DESTINATION arpifs
    INCLUDE_DIRS arpifs_intfb_includes
    PARALLEL ${FCM_PARALLEL}
  )
endif()

if( HAVE_OPENIFS_ONLY )
  if( NOT TARGET openifs_intfb ) 
    ecbuild_generate_fortran_interfaces(
      TARGET openifs_intfb
      DIRECTORIES dummy var
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/openifs
      DESTINATION openifs
      INCLUDE_DIRS openifs_intfb_includes
      PARALLEL ${FCM_PARALLEL}  
    )
  endif() 
endif()


if( NOT TARGET wam_intfb )
  ecbuild_generate_fortran_interfaces(
    TARGET wam_intfb
    DIRECTORIES Wam_oper
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/wam
    DESTINATION wam
    INCLUDE_DIRS wam_intfb_includes
    PARALLEL ${FCM_PARALLEL}
  )
endif()

if( NOT TARGET scmec_intfb )
  ecbuild_generate_fortran_interfaces(
    TARGET scmec_intfb
    DIRECTORIES source 
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/scmec
    DESTINATION scmec
    INCLUDE_DIRS scmec_intfb_includes
    PARALLEL ${FCM_PARALLEL}
  )
endif()
