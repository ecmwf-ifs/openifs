# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

macro( ectrans_find_lapack )
  # This macro sets the LAPACK_LIBRARIES variable
  # IF MKL is preferred, unless ENABLE_MKL=OFF

  if( HAVE_MKL )
    set( LAPACK_LIBRARIES ${MKL_LIBRARIES} )
  else()
    # Following libsci code should disappear soon, with more recent cmake versions (needs more investigation)
    set( _cray_libsci_loaded $ENV{CRAY_LIBSCI_DIR} )
    if( _cray_libsci_loaded )
        set( _CRAY_PRGENV $ENV{PE_ENV} )
        string( TOLOWER "${_CRAY_PRGENV}" _cray_prgenv  )
        set( LAPACK_LIBRARIES sci_${_cray_prgenv} )
        ecbuild_debug( "LAPACK found, already loaded as part of Cray's libsci" )
    else()
      ecbuild_find_package( NAME LAPACK REQUIRED )
      if( TARGET lapack )
        ecbuild_debug( "LAPACK found as CMake target lapack" )
        set( LAPACK_LIBRARIES lapack )
      endif()
    endif()
  endif()
  ecbuild_debug_var( LAPACK_LIBRARIES )

  set( LAPACK_sp ${LAPACK_LIBRARIES} CACHE STRING "ectrans: Double precision LAPACK libraries" )
  set( LAPACK_dp ${LAPACK_LIBRARIES} CACHE STRING "ectrans: Single precision LAPACK libraries" )

  set( LAPACK_LINK PRIVATE )

  ### Following is a hack that should be removed when there is no more Cray computer at ECMWF
  #   It allows to use a different LAPACK library for single and double precision, to be able to
  #   stay bitreproducible for double precision in operations of CY47R1

  set( _cray_libsci_loaded $ENV{CRAY_LIBSCI_DIR} )

  if( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
    if( HAVE_MKL AND ECTRANS_CRAYHACK_DOUBLE_PRECISION_WITHOUT_MKL )
      # Following libsci code should disappear soon, with more recent cmake versions (needs more investigation)
      if( _cray_libsci_loaded )
        set( _CRAY_PRGENV $ENV{PE_ENV} )
        string( TOLOWER "${_CRAY_PRGENV}" _cray_prgenv  )
        set( LAPACK_dp sci_${_cray_prgenv} )
        ecbuild_debug( "LAPACK found, already loaded as part of Cray's libsci" )
      else()
        ecbuild_find_package( NAME LAPACK REQUIRED )
        set( LAPACK_dp ${LAPACK_LIBRARIES} )
        if( TARGET lapack )
          ecbuild_debug( "LAPACK found as CMake target lapack" )
          set( LAPACK_dp lapack )
        endif()
      endif()
    endif()
  endif()

  if( _cray_libsci_loaded )
    if( NOT LAPACK_sp MATCHES "sci" OR NOT LAPACK_dp MATCHES "sci" )
      ecbuild_warn( "Danger! Cray's libsci is loaded, which is different from selected LAPACK. "
                    "No guarantees on link order can be made for the final executable.")
      set( LAPACK_LINK PUBLIC )
    endif()
  endif()

endmacro()
