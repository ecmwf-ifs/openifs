! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE SPECTRA (KIJS, KIJL, ZGAMMA, SA, SB, FP, ALPHAJ, THES, FL1)

! ----------------------------------------------------------------------

!**** *SPECTRA*  - COMPUTATION OF 2-D SPECTRA FOR ONE BLOCK.

!     S. HASSELMANN  - JULY 1990
!     H. GUNTHER     - DECEMBER 1990   MODIFIED FOR CYCLE_4.

!*    PURPOSE.
!     --------

!       INITIALISATION OF A BLOCK BY 2-D SPECTRA.

!**   INTERFACE.
!     ----------

!       *CALL* *SPECTRA (KIJS, KIJL, ZGAMMA, SA, SB, FP, ALPHAJ, THES, FL1)
!          *KIJS*    INTEGER  FIRST POINT IN BLOCK.
!          *KIJL*    INTEGER  LAST  POINT IN BLOCK.
!          *ZGAMMA*  REAL      OVERSHOOT FACTOR.
!          *SA*      REAL      LEFT PEAK WIDTH.
!          *SB*      REAL      RIGHT PEAK WIDTH.
!          *FP*      REAL      PEAK FREQUENCY OF SPECTRA IN A BLOCK (HZ).
!          *ALPHJ*   REAL      ALPHA PARAMETER OF SPECTRA IN A BLOCK.
!          *THES*    REAL      MEAN DIRECTION OF SPECTRA IN A BLOCK (RAD).

!          *FL1*      REAL      2-D SPECTRUM FOR EACH GRID POINT

!     METHOD.
!     -------
!       1-D JONSWAP SPECTRA AND COSINE**2 SPREADING FUNCTIONS ARE
!       COMPUTED FROM GIVEN WINDS AND PARAMETERS AT EACH GRID POINT.
!       THE 1-D SPECTRA ARE SPREAD OVER THE DIRECTIONS BY MULTIPLICATION
!       WITH THE SPREADING FUNCTION.

!     EXTERNALS.
!     ----------

!       *SPT*       - COMPUTATION OF COS**2 SPREADING FUNCTION.


!     REFERENCES.
!     -----------

!       K.HASSELMAN,D.B.ROSS,P.MUELLER AND W.SWELL
!          A PARAMETRIC WAVE PREDICTION MODEL
!          JOURNAL OF PHYSICAL OCEANOGRAPHY, VOL. 6, NO. 2, MARCH 1976

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "jonswap.intfb.h"
#include "spr.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), INTENT(IN)  :: ZGAMMA, SA, SB
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN):: FP, ALPHAJ, THES
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT):: FL1


      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NANG) :: STH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: ST
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: ET

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPECTRA',0,ZHOOK_HANDLE)

!*    1. COMPUTE JONSWAP SPECTRUM.
!        -------------------------

      CALL JONSWAP (ALPHAJ, ZGAMMA, SA, SB, FP, KIJS, KIJL, ET) 

! ----------------------------------------------------------------------

!*    2. COMPUTATION OF SPREADING FUNCTION.
!        ----------------------------------

      DO IJ=KIJS,KIJL
        CALL SPR (NANG, THES(IJ), TH, STH)
        DO K=1,NANG
          ST(IJ,K) = STH(K)
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    3. COMPUTATION OF 2-D SPECTRUM.
!        ----------------------------

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            FL1(IJ,K,M) = ET(IJ,M)*ST(IJ,K)
          ENDDO
        ENDDO
      ENDDO

IF (LHOOK) CALL DR_HOOK('SPECTRA',1,ZHOOK_HANDLE)

END SUBROUTINE SPECTRA
