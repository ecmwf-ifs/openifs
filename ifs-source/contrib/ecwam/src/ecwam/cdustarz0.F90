! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE CDUSTARZ0(KIJS, KIJL, U10, ZREF, CD, USTAR, Z0)

! ----------------------------------------------------------------------

!**** *CDUSTARZ0* - COMPUTATION OF FRICTION VELOCITY OVER SEA


!*    PURPOSE.
!     ---------

!       TO COMPUTE CD, USTAR AND Z0 USING CD(U10) RELATION

!*    METHOD.
!     -------

!     FIND CD USIND CD v U10 RELATION, THEN USTAR.
!     Z0 IS DERIVED FROM THE NEUTRAL LOG PROFILE: U10 = (USTAR/XKAPPA)*LOG((ZREF+Z0)/Z0)

!     NOTE. THIS IS ONLY VALID FOR ZREF=10m.

!     REFERENCE.
!     ----------

!     USE HERSBACH 2011 FOR CD(U10) (SEE ALSO EDSON et al. 2013)


! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWPCONS , ONLY : C1CD, C2CD, P1CD, P2CD, CDMAX, EPSUS
USE YOWPHYS  , ONLY : XKAPPA
USE YOWWIND  , ONLY : WSPMIN

USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
REAL(KIND=JWRB), INTENT(IN) :: ZREF
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: U10
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) ::  CD, USTAR, Z0


INTEGER(KIND=JWIM) :: IJ

REAL(KIND=JWRB), PARAMETER :: Z0MIN = 0.000001_JWRB
REAL(KIND=JWRB) :: U10P
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CDUSTARZ0',0,ZHOOK_HANDLE)

DO IJ=KIJS,KIJL
  U10P = MAX(U10(IJ),WSPMIN)
  CD(IJ) = MIN((C1CD + C2CD*U10P**P1CD)*U10P**P2CD,CDMAX)
  USTAR(IJ) = MAX(SQRT(CD(IJ))*U10P, EPSUS)
  Z0(IJ) = MAX(ZREF /(EXP(XKAPPA*MIN(U10P/USTAR(IJ),100.0_JWRB))-1.0_JWRB), Z0MIN)
ENDDO

IF (LHOOK) CALL DR_HOOK('CDUSTARZ0',1,ZHOOK_HANDLE)

END SUBROUTINE CDUSTARZ0
