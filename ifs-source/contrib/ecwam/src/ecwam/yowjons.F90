! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWJONS

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    **  *JONS* - JONSWAP PARAMETERS.

      REAL(KIND=JWRB), PARAMETER :: AJONS = 2.84_JWRB
      REAL(KIND=JWRB), PARAMETER :: BJONS = 0.033_JWRB
      REAL(KIND=JWRB), PARAMETER :: DJONS = -3.0_JWRB/10.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: EJONS = 2.0_JWRB/3.0_JWRB

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *AJONS*     REAL      CONSTANT USED IN FETCH LIMITED GROWTH
!                            SEE ROUTINE PEAK
!      *BJONS*     REAL      CONSTANT USED IN FETCH LIMITED GROWTH
!                            SEE ROUTINE PEAK
!      *DJONS*     REAL      CONSTANT USED IN FETCH LIMITED GROWTH
!                            SEE ROUTINE PEAK
!      *EJONS*     REAL      CONSTANT USED IN FETCH LIMITED GROWTH
!                            SEE ROUTINE PEAK
! ----------------------------------------------------------------------
      END MODULE YOWJONS
