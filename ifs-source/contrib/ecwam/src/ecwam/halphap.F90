! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE HALPHAP(KIJS, KIJL, WAVNUM, COSWDIF, FL1, HALP)

! ----------------------------------------------------------------------

!**** *HALPHAP* - COMPUTATION OF 1/2 PHILLIPS PARAMETER


!**   INTERFACE.
!     ----------

!       *CALL* *HALPHAP(KIJS, KIJL, WAVNUM, UDIR, FL1, HALP)
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *WAVNUM* - WAVE NUMBER
!          *COSWDIF*- COSINE ( WIND SPEED DIRECTION - WAVE DIRECTIONS)
!          *FL1*    - SPECTRA
!          *HALP*   - 1/2 PHILLIPS PARAMETER 

!     METHOD.
!     -------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       , TH       , FR5      ,DELTH
      USE YOWPARAM , ONLY : NANG     , NFRE
      USE YOWPCONS , ONLY : G        , ZPI      ,ZPI4GM2
      USE YOWPHYS  , ONLY : ALPHAPMAX

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "femean.intfb.h"
#include "meansqs_lf.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG), INTENT(IN) :: COSWDIF 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: HALP

      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: ZLNFRNFRE
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ALPHAP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: XMSS, EM, FM, F1D 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: WD 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: FLWD

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('HALPHAP',0,ZHOOK_HANDLE)

      ZLNFRNFRE = LOG(FR(NFRE))

      ! Find spectrum in wind direction
      DO K = 1, NANG
        DO IJ = KIJS, KIJL
           WD(IJ,K) = 0.5_JWRB + 0.5_JWRB * SIGN(1.0_JWRB, COSWDIF(IJ,K))
        ENDDO
      ENDDO

      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = KIJS, KIJL
             FLWD(IJ,K,M) = FL1(IJ,K,M) * WD(IJ,K) 
          ENDDO
        ENDDO
      ENDDO

      CALL MEANSQS_LF(NFRE, KIJS, KIJL, FLWD, WAVNUM, XMSS)

      CALL FEMEAN (KIJS, KIJL, FLWD, EM, FM)

      DO IJ = KIJS, KIJL
        IF (EM(IJ) > 0.0_JWRB .AND. FM(IJ) < FR(NFRE-2) ) THEN
          ALPHAP(IJ) = XMSS(IJ) /(ZLNFRNFRE - LOG(FM(IJ)))
          IF ( ALPHAP(IJ) > ALPHAPMAX ) THEN
            ! some odd cases, revert to tail value
            F1D(IJ) = 0.0_JWRB
            DO K = 1, NANG
              F1D(IJ) = F1D(IJ) + FLWD(IJ,K,NFRE)*DELTH
            ENDDO
            ALPHAP(IJ) = ZPI4GM2*FR5(NFRE)*F1D(IJ)
          ENDIF
        ELSE
          F1D(IJ) = 0.0_JWRB
          DO K = 1, NANG
            F1D(IJ) = F1D(IJ) + FLWD(IJ,K,NFRE)*DELTH
          ENDDO
          ALPHAP(IJ) = ZPI4GM2*FR5(NFRE)*F1D(IJ)
        ENDIF
      ENDDO

!     1/2 ALPHAP:
      HALP(KIJS:KIJL) = 0.5_JWRB*MIN(ALPHAP(KIJS:KIJL), ALPHAPMAX)

IF (LHOOK) CALL DR_HOOK('HALPHAP',1,ZHOOK_HANDLE)

END SUBROUTINE HALPHAP
