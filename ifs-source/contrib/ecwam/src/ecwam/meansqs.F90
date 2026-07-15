! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MEANSQS(XKMSS, KIJS, KIJL, F, WAVNUM, USTAR, COSWDIF, XMSS)

! ----------------------------------------------------------------------

!**** *MEANSQS* - COMPUTATION OF MEAN SQUARE SLOPE UP TO WAVE NUMBER XKMSS.

!*    PURPOSE.
!     --------

!       COMPUTE MEAN SQUARE SLOPE AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *MEANSQS (XKMSS, KIJS, KIJL, F, WAVNUM, USTAR, COSWDIF, XMSS)*
!              *XKMSS*   - WAVE NUMBER CUT OFF
!              *KIJS*    - INDEX OF FIRST GRIDPOINT
!              *KIJL*    - INDEX OF LAST GRIDPOINT
!              *F*       - SPECTRUM.
!              *WAVNUM*  - WAVE NUMBER.
!              *USTAR*   - NEW FRICTION VELOCITY IN M/S (INPUT).
!              *COSWDIF* - COSINE (WIND SPEED DIRECTION - WAVES DIRECTIONS)
!              *XMSS*    - MEAN SQUARE SLOPE (OUTPUT).

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : G        ,GM1      ,ZPI
      USE YOWFRED  , ONLY : FR       ,ZPIFR    ,FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "halphap.intfb.h"
#include "meansqs_gc.intfb.h"
#include "meansqs_lf.intfb.h"

      REAL(KIND=JWRB), INTENT(IN) :: XKMSS 
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: F
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: USTAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG), INTENT(IN) :: COSWDIF 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: XMSS 


      INTEGER(KIND=JWIM) :: IJ, NFRE_MSS, NFRE_EFF

      REAL(KIND=JWRB) :: XLOGFS, FCUT
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: FD
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP1, TEMP2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: XMSSLF, XMSS_TAIL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: HALP, FRGC

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MEANSQS',0,ZHOOK_HANDLE)

!*    1. COMPUTE THE GRAVITY-CAPILLARY CONTRIBUTION TO MSS 
!        -------------------------------------------------

!     COMPUTE THE PHILLIPS PARAMETER
      CALL HALPHAP(KIJS, KIJL, WAVNUM, COSWDIF, F, HALP)

!     GRAVITY-CAPILLARY CONTRIBUTION TO MSS
      CALL MEANSQS_GC(XKMSS, KIJS, KIJL, HALP, USTAR, XMSS, FRGC)

!     MODEL SPECTRUM PART
      FCUT = SQRT(G*XKMSS)/ZPI
      NFRE_MSS = INT(LOG(FCUT/FR(1))/LOG(FRATIO))+1 
      NFRE_EFF = MIN(NFRE,NFRE_MSS)

      CALL MEANSQS_LF (NFRE_EFF, KIJS, KIJL, F, WAVNUM, XMSSLF)
      XMSS(:) = XMSS(:) + XMSSLF(:)

!     ADD TAIL CORRECTION TO MEAN SQUARE SLOPE (between FR(NFRE_EFF) and FRGC).

      XLOGFS = LOG(FR(NFRE_EFF))
      DO IJ=KIJS,KIJL
        XMSS_TAIL(IJ) = 2.0_JWRB*HALP(IJ)*MAX((LOG(MIN(FRGC(IJ),FCUT)) - XLOGFS), 0.0_JWRB)
        XMSS(IJ) = XMSS(IJ) + XMSS_TAIL(IJ)
      ENDDO

      IF (LHOOK) CALL DR_HOOK('MEANSQS',1,ZHOOK_HANDLE)

      END SUBROUTINE MEANSQS
