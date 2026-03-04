! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE YOMMSC

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    Basic parameters

! NINTLEN             - length of INTEGER_M variables in bytes
! NLOGLEN             - length of logicals in bytes
! NREALEN             - length of REAL_B variables in bytes
! NDBLLEN             - length of REAL_H variables in bytes
! N_DEFAULT_REAL_KIND - The KIND value of the default REAL kind
! N_DOUBLE_KIND       - The KIND value of double precision

INTEGER(KIND=JPIM) :: NINTLEN
INTEGER(KIND=JPIM) :: NREALEN
INTEGER(KIND=JPIM) :: NLOGLEN
INTEGER(KIND=JPIM) :: NDBLLEN

REAL, PRIVATE             :: Z_DEFAULT_REAL      ! intentionally not REAL(KIND=JPRB)
INTEGER(KIND=JPIM), PARAMETER :: N_DEFAULT_REAL_KIND = KIND(Z_DEFAULT_REAL)

DOUBLE PRECISION, PRIVATE :: DL_DOUBLE_PRECISION ! intentionally not REAL(KIND=JPRH)
INTEGER(KIND=JPIM), PARAMETER :: N_DOUBLE_KIND       = KIND(DL_DOUBLE_PRECISION)

!     ------------------------------------------------------------------
END MODULE YOMMSC
