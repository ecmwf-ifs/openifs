! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE DEPTHPRPT(KIJS, KIJL, DEPTH,          &
 &                   WAVNUM, CINV, CGROUP,       &
 &                   XK2CG, OMOSNH2KD, STOKFAC)

! ----------------------------------------------------------------------

!**** *DEPTHPRPT* -


!*    PURPOSE.
!     -------

! COMPUTES WAVE PROPERTIES WHICH DEPEND ON WATER DEPTH AND FREQUENCIES.


! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWFRED  , ONLY : FR, ZPIFR 
USE YOWPARAM , ONLY : NFRE
USE YOWPCONS , ONLY : G, PI

USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

IMPLICIT NONE

#include "aki.intfb.h"

INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH           ! WATER DEPTH

REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(OUT) :: WAVNUM    ! WAVE NUMBER 
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(OUT) :: CINV      ! RECIPROCAL OF THE PHASE VELOCITY (1/c)
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(OUT) :: CGROUP    ! GROUP SPEED
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(OUT) :: XK2CG     ! (WAVE NUMBER)**2 * GROUP SPEED
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(OUT) :: OMOSNH2KD ! OMEGA / SINH(2KD)
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(OUT) :: STOKFAC   ! FACTOR TO COMPUTE SURFACE STOKES DRIFT FROM SPECTRUM 2*G*K**2/(OMEGA*TANH(2KD))


INTEGER(KIND=JWIM) :: M, IJ

REAL(KIND=JWRB) :: GH, OM, AK, AKD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DEPTHPRPT',0,ZHOOK_HANDLE)

GH = G/(4.0_JWRB*PI)

DO M = 1, NFRE
  OM=ZPIFR(M)
  DO IJ = KIJS, KIJL
    AK = AKI(OM,DEPTH(IJ))
    WAVNUM(IJ,M) = AK
    AKD = AK*DEPTH(IJ)
    IF (AKD <= 10.0_JWRB) THEN
      CGROUP(IJ,M) = 0.5_JWRB*SQRT(G*TANH(AKD)/AK)*(1.0_JWRB+2.0_JWRB*AKD/SINH(2.0_JWRB*AKD))
      OMOSNH2KD(IJ,M) = OM/SINH(2.0_JWRB*AKD)
      STOKFAC(IJ,M) = 2.0_JWRB*G*AK**2/(OM*TANH(2.0_JWRB*AKD))
    ELSE
      CGROUP(IJ,M) = GH/FR(M)
      OMOSNH2KD(IJ,M) = 0.0_JWRB
      STOKFAC(IJ,M) = 2.0_JWRB/G*OM**3
    ENDIF

    CINV(IJ,M) = WAVNUM(IJ,M)/OM
    XK2CG(IJ,M) = WAVNUM(IJ,M)**2 * CGROUP(IJ,M)

  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('DEPTHPRPT',1,ZHOOK_HANDLE)

END SUBROUTINE DEPTHPRPT

