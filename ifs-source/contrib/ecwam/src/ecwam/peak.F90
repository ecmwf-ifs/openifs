! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE PEAK (KIJS, KIJL, FETCH, FPMAX, U10, FP, ALPHAJ)

! ----------------------------------------------------------------------

!**** *PEAK* - COMPUTE JONSWAP PARAMETERS FROM WINDSPEED.

!     S. HASSELMANN  - JULY 1990
!     H. GUNTHER     - DECEMBER 1990   MODIFIED FOR CYCLE_4.
!     J. BIDLOT      - FEBRUARY 1995   MESSAGE PASSING

!*    PURPOSE.
!     --------

!       COMPUTES FOR EACH GRID POINT OF A BLOCK THE PEAK FREQUENCY
!       FROM A FETCH LAW AND THE JONSWAP ALPHAJONS FROM THE ALPHAJONS NY
!       RELATION.

!**   INTERFACE.
!     ----------

!       *CALL* *PEAK (KIJL, KIJS, FETCH, FPMAX, U10, FP, ALPHAJ)*
!          *KIJS*    INTEGER  FIRST POINT IN BLOCK.
!          *KIJL*    INTEGER  LAST  POINT IN BLOCK.
!          *FETCH*   REAL     FETCH TO BE USED (METRES).
!          *FPMAX*   REAL     MAXIMUM PEAK FREQUENCY (HERTZ).
!          *U10*     REAL     WIND SPEED.
!          *FP*      REAL     PEAK FREQUENCY.
!          *ALPHAJ*  REAL     ALPHA PARAMETER.

!     METHOD.
!     -------

!       FP = A * (G*FETCH/U_10**2)**D    A = 2.84
!       FP = MAX [FP, 0.13]              D = -3./10.
!       FP = MIN [FP, FRMAX*U_10/G]      b = 0.033
!       ALPHAJ = B * FP-**2/3            B = 0.033
!       ALPHAJ = MAX [ALPHAJ, 0.0081]
!       FP = G/U_10*FP

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCES.
!     -----------

!       K.HASSELMAN,D.B.ROOS,P.MUELLER AND W.SWELL
!          A PARAMETRIC WAVE PREDICTION MODEL
!          JOURNAL OF PHSICAL OCEANOGRAPHY, VOL. 6, NO. 2, MARCH 1976.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWJONS  , ONLY : AJONS    ,BJONS    ,DJONS    ,EJONS
      USE YOWPCONS , ONLY : G

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), INTENT(IN) :: FETCH, FPMAX
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: U10
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(OUT) :: FP
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(INOUT) :: ALPHAJ 


      INTEGER(KIND=JWIM) :: IJ
      REAL(KIND=JWRB) :: GX, GXU, UG

! ----------------------------------------------------------------------

!*    1. COMPUTE VALUES FROM FETCH LAWS.
!        -------------------------------

      GX = G * FETCH
      DO IJ = KIJS, KIJL
        IF (U10(IJ) > 0.1E-08_JWRB) THEN
          GXU = GX/(U10(IJ)*U10(IJ))
          UG = G/U10(IJ)
          FP(IJ) = AJONS * GXU ** DJONS
          FP(IJ) = MAX(0.13_JWRB, FP(IJ))
          FP(IJ) = MIN(FP(IJ), FPMAX/UG)
          ALPHAJ(IJ) = BJONS * FP(IJ)**EJONS
          ALPHAJ(IJ) = MAX(ALPHAJ(IJ), 0.0081_JWRB)
          FP(IJ) = FP(IJ)*UG
        ELSE
          FP(IJ) = 0.0_JWRB
        ENDIF
      ENDDO

      END SUBROUTINE PEAK
