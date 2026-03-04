! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_SUNDIS(KIDIA, KFDIA, KLON, KLEV, KMONTH, KDAY, PRJ)
!
! Purpose:
!   Applies a date-dependent, dimensionless Earth-Sun distance scaling 
!   factor to the TM5 photolysis-rate array used by IFS chemistry. The 
!   correction accounts for the seasonal variation in solar irradiance 
!   due to the changing Earth-Sun distance (commonly treated as an 
!   inverse-square distance effect).
!
! Interface:
!   Called from the TM5 chemistry interface (e.g., CHEM_tm5) during 
!   photolysis evaluation to scale photolysis rates for the current 
!   model date.
!
!   Inputs:  Horizontal index range, array dimensions, calendar month/day
!   In/Out:  PRJ array - photolysis rates modified in-place
!
! Scientific basis:
!   Uses Fourier-series approximation (Spencer, 1971) widely used in 
!   atmospheric photolysis schemes:
!
!   F = 1.000110 + 0.034221*cos(theta) + 0.001280*sin(theta)
!                + 0.000719*cos(2*theta) + 0.000077*sin(2*theta)
!
!   where theta = 2*pi/365 * (day_of_year - 1 + 0.5)
!
! Arguments:
!   KIDIA  (input)  : First active horizontal index in PRJ to update
!   KFDIA  (input)  : Last active horizontal index in PRJ to update
!   KLON   (input)  : Declared horizontal dimension of PRJ
!   KLEV   (input)  : Number of vertical levels in PRJ
!   KMONTH (input)  : Calendar month (1-12)
!   KDAY   (input)  : Day of month (1-31)
!   PRJ    (in/out) : Photolysis rates (KLON, KLEV, NPHOTO)
!
!------------------------------------------------------------------------------

USE PARKIND1,      ONLY: JPIM, JPRB
USE YOMCST,        ONLY: RPI
USE TM5_PHOTOLYSIS, ONLY: NPHOTO
USE YOMHOOK,       ONLY: LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

! Arguments
INTEGER(KIND=JPIM), INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)    :: KLON
INTEGER(KIND=JPIM), INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM), INTENT(IN)    :: KMONTH
INTEGER(KIND=JPIM), INTENT(IN)    :: KDAY
REAL(KIND=JPRB),    INTENT(INOUT) :: PRJ(KLON, KLEV, NPHOTO)

! Local variables
INTEGER(KIND=JPIM) :: IDOY          ! Day of year
INTEGER(KIND=JPIM) :: IMONTH        ! Loop counter
INTEGER(KIND=JPIM) :: JL            ! Horizontal loop index
INTEGER(KIND=JPIM) :: JK            ! Vertical loop index
INTEGER(KIND=JPIM) :: JN            ! Photolysis channel loop index
INTEGER(KIND=JPIM), PARAMETER :: IDAYS(12) = &
  (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

REAL(KIND=JPRB) :: ZTHETA           ! Day angle (radians)
REAL(KIND=JPRB) :: ZSTHETA          ! sin(theta)
REAL(KIND=JPRB) :: ZCTHETA          ! cos(theta)
REAL(KIND=JPRB) :: ZS2THETA         ! sin(2*theta)
REAL(KIND=JPRB) :: ZC2THETA         ! cos(2*theta)
REAL(KIND=JPRB) :: ZFACTOR          ! Earth-Sun distance correction factor

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TM5_SUNDIS', 0, ZHOOK_HANDLE)

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
ZFACTOR = 1.000110_JPRB &
  & + 0.034221_JPRB * ZCTHETA + 0.001280_JPRB * ZSTHETA &
  & + 0.000719_JPRB * ZC2THETA + 0.000077_JPRB * ZS2THETA

! Apply correction factor to photolysis rates
DO JN = 1, NPHOTO
  DO JK = 1, KLEV
    DO JL = KIDIA, KFDIA
      PRJ(JL, JK, JN) = PRJ(JL, JK, JN) * ZFACTOR
    ENDDO
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('TM5_SUNDIS', 1, ZHOOK_HANDLE)

END SUBROUTINE TM5_SUNDIS