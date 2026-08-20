! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE PEAK_FREQ(KIJS, KIJL, FL1, FP)
 
!***  *PEAK*   DETERMINES THE PEAK FREQUENCY OF THE SPECTRUM
 
!     PETER JANSSEN
 
!     PURPOSE.
!     --------
 
!              DETERMINATION OF PEAK PARAMETERS
 
!     INTERFACE.
!     ----------
!              *CALL*  *PEAK_FREQ(FL1,KIJS,KIJL,FP)*
 
!                       INPUT:
!                            *FL1*    -  2D-SPECTRUM
!                            *KIJS*   - FIRST GRIDPOINT              
!                            *KIJL*   - LAST GRIDPOINT  
 
!                       OUTPUT: 
!                            *FP*     - PEAK FREQUENCY
 
!     METHOD.
!     -------
!              NONE
 
!     EXTERNALS.
!     ----------
!              NONE
 
!-----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : G        ,PI       ,ZPI      
      USE YOWFRED  , ONLY : FR       ,FRATIO   ,DELTH
      USE YOWICE   , ONLY : FLMIN
      USE YOWPARAM , ONLY : NANG     ,NFRE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
 
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG, NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: FP


      INTEGER(KIND=JWIM) :: IJ, M, K
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL) ::  MMAX

      REAL(KIND=JWRB), PARAMETER :: WL=0.1_JWRB
      REAL(KIND=JWRB), PARAMETER :: WC=0.8_JWRB
      REAL(KIND=JWRB), PARAMETER :: WR=0.1_JWRB

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: A,B,C,XP1,XM1,F1P1,F1M1,F10,TF,XMAX_MIN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: F1D, F1DSM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: XMAX

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('PEAK_FREQ',0,ZHOOK_HANDLE)


!***  1. DETERMINE ONE-DIMENSIONAL FREQUENCY SPECTRUM.
!     ------------------------------------------------
 
      DO M=1,NFRE
        K=1
        DO IJ=KIJS,KIJL
          F1D(IJ,M) = FL1(IJ,K,M)*DELTH
        ENDDO
        DO K=2,NANG
          DO IJ=KIJS,KIJL
            F1D(IJ,M) = F1D(IJ,M)+FL1(IJ,K,M)*DELTH
          ENDDO
        ENDDO
      ENDDO

!     SMOOTH F1D FOR A BETTER EXTIMATE OF THE PEAK
      TF=1._JWRB/(FRATIO**5)
      M=1
      DO IJ=KIJS,KIJL
        F1DSM(IJ,M)=WC*F1D(IJ,M)+WR*F1D(IJ,M+1)
      ENDDO
      DO M=2,NFRE-1
        DO IJ=KIJS,KIJL
          F1DSM(IJ,M)=WL*F1D(IJ,M-1)+WC*F1D(IJ,M)+WR*F1D(IJ,M+1)
        ENDDO
      ENDDO
      M=NFRE
      DO IJ=KIJS,KIJL
        F1DSM(IJ,M)=WL*F1D(IJ,M-1)+(WC+WR*TF)*F1D(IJ,M)
      ENDDO
 
!***  2. DETERMINE F_PEAK USING QUADRATIC FIT.
!     ---------------------------------------
 
!     MAX OF 1-D SPECTRUM
 
      XMAX_MIN = NANG*FLMIN*DELTH
      DO IJ=KIJS,KIJL
        XMAX(IJ) = XMAX_MIN 
        MMAX(IJ) = NFRE-1
      ENDDO

      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          IF (F1DSM(IJ,M) >= XMAX(IJ)) THEN
            MMAX(IJ) = M
            XMAX(IJ) = F1DSM(IJ,M)
          ENDIF
        ENDDO
      ENDDO

      DO IJ=KIJS,KIJL
        FP(IJ) = FR(MMAX(IJ))
      ENDDO

!     DETERMINE QUADRATIC FIT.
      DO IJ=KIJS,KIJL
        IF (XMAX(IJ) > XMAX_MIN .AND. MMAX(IJ) > 1 .AND. MMAX(IJ) < NFRE ) THEN
          XP1 = FR(MMAX(IJ)+1)-FR(MMAX(IJ))
          XM1 = FR(MMAX(IJ)-1)-FR(MMAX(IJ))
          F10  = F1DSM(IJ,MMAX(IJ))
          F1P1 = (F1DSM(IJ,MMAX(IJ)+1)-F10)/XP1
          F1M1 = (F1DSM(IJ,MMAX(IJ)-1)-F10)/XM1

!!!          A = F10
          B = (XM1*F1P1-XP1*F1M1)/(XM1-XP1)
          C = (F1M1-F1P1)/(XM1-XP1)
        
          IF (C < 0._JWRB) THEN
            FP(IJ) = FR(MMAX(IJ))-B/(2._JWRB*C)
!!!            XMAX(IJ) = A-B**2/(4._JWRB*C)
          ENDIF
        ENDIF
      ENDDO 
 
      IF (LHOOK) CALL DR_HOOK('PEAK_FREQ',1,ZHOOK_HANDLE)

      END SUBROUTINE PEAK_FREQ
