! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUNDISTCORR(KMONTH, KDAY, PFACTOR)
!
! Purpose:
!   Computes a dimensionless scaling factor that represents the annual
!   variation of the Earth-Sun distance. This factor multiplies photolysis
!   rates for a given calendar date. 
!    
!   Called from IFS chemistry (e.g. chem_bascoe/tm5), and required for 
!   photolysis evaluation.
!
! Interface:
!   
!   Inputs:  Calendar month and day (no year)
!   Output:  Earth-Sun distance correction factor (dimensionless)
!
! Scientific basis:
!   Uses Fourier-series approximation from:
!   Spencer, J.W. (1971), "Fourier series representation of the position
!   of the sun", Search, 2, 172.
!
! Arguments:
!   KMONTH  (input)  : Calendar month (1-12)
!   KDAY    (input)  : Day of month (1-31)
!   PFACTOR (output) : Earth-Sun distance correction factor (dimensionless)
!
!------------------------------------------------------------------------------

USE PARKIND1, ONLY: JPIM, JPRB
USE YOMCST,   ONLY: RPI
USE YOMHOOK,  ONLY: LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

! Arguments
INTEGER(KIND=JPIM), INTENT(IN)  :: KMONTH
INTEGER(KIND=JPIM), INTENT(IN)  :: KDAY
REAL(KIND=JPRB),    INTENT(OUT) :: PFACTOR

! Local variables
INTEGER(KIND=JPIM) :: IDOY          ! Day of year
INTEGER(KIND=JPIM) :: IMONTH        ! Loop counter
INTEGER(KIND=JPIM), PARAMETER :: IDAYS(12) = &
  (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

REAL(KIND=JPRB) :: ZTHETA           ! Day angle (radians)
REAL(KIND=JPRB) :: ZSTHETA          ! sin(theta)
REAL(KIND=JPRB) :: ZCTHETA          ! cos(theta)
REAL(KIND=JPRB) :: ZS2THETA         ! sin(2*theta)
REAL(KIND=JPRB) :: ZC2THETA         ! cos(2*theta)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUNDISTCORR', 0, ZHOOK_HANDLE)

! Convert calendar date to day-of-year (1-based)
IDOY = KDAY
DO IMONTH = 1, KMONTH - 1
  IDOY = IDOY + IDAYS(IMONTH)
ENDDO

! Compute day angle: theta = 2*pi/365 * (doy - 1 + 0.5)
! The +0.5 offset places the angle at the middle of the day
ZTHETA = (2.0_JPRB * RPI / 365.0_JPRB) * (REAL(IDOY - 1, JPRB) + 0.5_JPRB)

! Compute trigonometric functions
ZSTHETA = SIN(ZTHETA)
ZCTHETA = COS(ZTHETA)

! Compute double-angle terms using identities to avoid extra intrinsic calls
! sin(2*theta) = 2 * sin(theta) * cos(theta)
! cos(2*theta) = cos^2(theta) - sin^2(theta)
ZS2THETA = 2.0_JPRB * ZSTHETA * ZCTHETA
ZC2THETA = ZCTHETA * ZCTHETA - ZSTHETA * ZSTHETA

! Evaluate Spencer Fourier expansion for Earth-Sun distance correction
PFACTOR = 1.000110_JPRB &
  & + 0.034221_JPRB * ZCTHETA + 0.001280_JPRB * ZSTHETA &
  & + 0.000719_JPRB * ZC2THETA + 0.000077_JPRB * ZS2THETA

IF (LHOOK) CALL DR_HOOK('SUNDISTCORR', 1, ZHOOK_HANDLE)

END SUBROUTINE SUNDISTCORR