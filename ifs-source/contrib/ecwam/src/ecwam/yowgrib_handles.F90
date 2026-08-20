! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWGRIB_HANDLES

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE
!     GRIB HANDLES
      SAVE
      INTEGER(KIND=JWIM) :: NGRIB_HANDLE_WAM_I  = -99   ! DEFAULT GRIB TEMPLATE FOR WAM INTEGRATED PARAMETERS.
      INTEGER(KIND=JWIM) :: NGRIB_HANDLE_WAM_I2 = -99   ! GRIB2 TEMPLATE FOR WAM INTEGRATED PARAMETERS.
                                                        ! THAT CAN ONLY BE PRODUCED IN GRIB2 WHILE GRIB1 IS STILL USED
      INTEGER(KIND=JWIM) :: NGRIB_HANDLE_WAM_S  = -99   ! GRIB TEMPLATE FOR WAM SPECTRA.
      INTEGER(KIND=JWIM) :: NGRIB_HANDLE_IFS    = -99   ! FOR IFS FIELDS PASSED TO WAM (GRIB1 or GRIB2, depending on NGRIB_VERSION.
      INTEGER(KIND=JWIM) :: NGRIB_HANDLE_IFS2   = -99   ! FOR IFS FIELDS PASSED TO WAM GRIB2 ONLY.
      END MODULE YOWGRIB_HANDLES
