# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

## Apply workarounds for some known compilers

if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  if( CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8.7 )

    # Fix for IFS "CONGRAD: SPTSV/DPTSV returned non-zero info with crayftn 8.7.7 (cdt/18.12)
    ectrans_add_compile_options(
        SOURCES internal/ftinv_ctlad_mod.F90
        FLAGS   "-O0,fp1,omp")

  endif()
endif()

## Assemble sources

ecbuild_list_add_pattern( LIST trans_src
                          GLOB
                                sharedmem/*
                                algor/*
                                internal/*
                                external/*
                          QUIET
                        )

if( NOT HAVE_FFTW )
  ecbuild_list_exclude_pattern( LIST trans_src REGEX tpm_fftw.F90 )
endif()

set( HAVE_dp ${HAVE_DOUBLE_PRECISION} )
set( HAVE_sp ${HAVE_SINGLE_PRECISION} )

foreach( prec sp dp )
  if( HAVE_${prec} )

    ecbuild_add_library(
      TARGET           trans_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          ${trans_src}
      PUBLIC_INCLUDES  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/ectrans>
                       $<INSTALL_INTERFACE:include/ectrans>
                       $<INSTALL_INTERFACE:include>
      PUBLIC_LIBS      fiat parkind_${prec}
    )
    ectrans_target_fortran_module_directory(
      TARGET            trans_${prec}
      MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/trans_${prec}
      INSTALL_DIRECTORY module/trans_${prec}
    )
    target_link_libraries( trans_${prec} PUBLIC fiat parkind_${prec} )
    if( HAVE_FFTW )
      set( FFTW_LINK PRIVATE )
      if( LAPACK_${prec} MATCHES "mkl" AND NOT FFTW_LIBRARIES MATCHES "mkl" )
          ecbuild_warn( "Danger: Both MKL and FFTW are linked in trans_${prec}. "
                        "No guarantees on link order can be made for the final executable.")
          set( FFTW_LINK PUBLIC ) # Attempt anyway to give FFTW precedence
      endif()
      target_link_libraries( trans_${prec} ${FFTW_LINK} ${FFTW_LIBRARIES} )
      target_include_directories( trans_${prec} PRIVATE ${FFTW_INCLUDE_DIRS} )
      target_compile_definitions( trans_${prec} PRIVATE WITH_FFTW )
    endif()
    target_link_libraries( trans_${prec} ${LAPACK_LINK} ${LAPACK_${prec}} )
    if( HAVE_OMP )
      target_link_libraries( trans_${prec} PRIVATE OpenMP::OpenMP_Fortran )
    endif()

  endif()
endforeach()

## Install trans interface

file( GLOB trans_interface include/ectrans/* )
install(
  FILES        ${trans_interface}
  DESTINATION  include/ectrans
)
