! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE PEAK_ANG(KIJS, KIJL, FL1, XNU, SIG_TH)

!***  *PEAK_ANG*   DETERMINES ANGULAR WIDTH NEAR PEAK OF SPECTRUM

!     PETER JANSSEN

!     PURPOSE.
!     --------

!              DETERMINATION OF PEAK PARAMETERS

!     INTERFACE.
!     ----------
!              *CALL*  *PEAK_ANG(KIJS,KIJL,FL1,XNU,SIG_TH)*

!               INPUT:
!                  *KIJS*   - FIRST GRIDPOINT              
!                  *KIJL*   - LAST GRIDPOINT              
!                  *FL1*    - SPECTRUM
!               OUTPUT:
!                  *XNU*    - RELATIVE SPECTRAL WIDTH
!                  *SIG_TH* - RELATIVE WIDTH IN DIRECTION

!     METHOD.
!     -------
!              NONE

!     EXTERNALS.
!     ----------
!              NONE

!-----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMFR  ,DFIMOFR ,      &
     &               DFIMFR2  ,DELTH ,TH       ,SINTH    ,COSTH  ,      &
     &               WETAIL   ,WP1TAIL         ,WP2TAIL  ,FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
 
! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: XNU, SIG_TH


      INTEGER(KIND=JWIM) :: NSH
      INTEGER(KIND=JWIM) :: IJ, M, K
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL)::MMAX, MMSTART, MMSTOP
      REAL(KIND=JWRB), PARAMETER :: CONST_SIG = 1.0_JWRB
      REAL(KIND=JWRB) :: R1
      REAL(KIND=JWRB) :: DELT25, COEF_FR, COEF_FR2
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: ZEPSILON
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL) :: SUM0, SUM1, SUM2, SUM4, XMAX
      REAL,DIMENSION(KIJS:KIJL) :: TEMP
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL,NFRE) :: THMEAN, SUM_S, SUM_C

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('PEAK_ANG',0,ZHOOK_HANDLE)

!***  1. DETERMINE L-H SPECTRAL WIDTH OF THE 2-D SPECTRUM.
!     ---------------------------------------------------

      ZEPSILON=10._JWRB*EPSILON(ZEPSILON)

      NSH = 1 + INT(LOG(1.5_JWRB)/LOG(FRATIO)) 

      DO IJ=KIJS,KIJL
        SUM0(IJ)= ZEPSILON
        SUM1(IJ)= 0._JWRB
        SUM2(IJ)= 0._JWRB          
      ENDDO

      DO M=1,NFRE
        K=1
        DO IJ=KIJS,KIJL
          TEMP(IJ) = FL1(IJ,K,M)
        ENDDO
        DO K=2,NANG
          DO IJ=KIJS,KIJL
            TEMP(IJ) = TEMP(IJ)+FL1(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          SUM0(IJ) = SUM0(IJ)+TEMP(IJ)*DFIM(M)
          SUM1(IJ) = SUM1(IJ)+TEMP(IJ)*DFIMFR(M)
          SUM2(IJ) = SUM2(IJ)+TEMP(IJ)*DFIMFR2(M)
         ENDDO
      ENDDO

!     ADD TAIL CORRECTIONS
      DELT25 = WETAIL*FR(NFRE)*DELTH
      COEF_FR = WP1TAIL*DELTH*FR(NFRE)**2
      COEF_FR2 = WP2TAIL*DELTH*FR(NFRE)**3
      DO IJ=KIJS,KIJL
        SUM0(IJ) = SUM0(IJ)+DELT25*TEMP(IJ)
        SUM1(IJ) = SUM1(IJ)+COEF_FR*TEMP(IJ)
        SUM2(IJ) = SUM2(IJ)+COEF_FR2*TEMP(IJ)
      ENDDO

      DO IJ=KIJS,KIJL
        IF (SUM0(IJ) > ZEPSILON) THEN
          XNU(IJ) = SQRT(MAX(ZEPSILON,SUM2(IJ)*SUM0(IJ)/SUM1(IJ)**2-1._JWRB))
        ELSE
          XNU(IJ) = ZEPSILON 
        ENDIF 
      ENDDO

!***  2. DETERMINE ANGULAR WIDTH OF THE 2-D SPECTRUM.
!     ----------------------------------------------
 
!     MAX OF 2-D SPECTRUM
      DO IJ=KIJS,KIJL
        XMAX(IJ) = 0._JWRB
        MMAX(IJ) = 2
      ENDDO

      DO M=2,NFRE-1
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            IF (FL1(IJ,K,M) > XMAX(IJ)) THEN
              MMAX(IJ) = M
              XMAX(IJ) = FL1(IJ,K,M)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO IJ=KIJS,KIJL
        SUM1(IJ) = ZEPSILON
        SUM2(IJ) = 0._JWRB
      ENDDO
 
      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          SUM_S(IJ,M) = 0._JWRB
          SUM_C(IJ,M) = ZEPSILON 
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            SUM_S(IJ,M) = SUM_S(IJ,M) +SINTH(K)*FL1(IJ,K,M)
            SUM_C(IJ,M) = SUM_C(IJ,M) +COSTH(K)*FL1(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          THMEAN(IJ,M) = ATAN2(SUM_S(IJ,M),SUM_C(IJ,M))
        ENDDO
      ENDDO

      DO IJ=KIJS,KIJL
        MMSTART(IJ) = MAX(1,MMAX(IJ)-NSH) 
        MMSTOP(IJ)  = MIN(NFRE,MMAX(IJ)+NSH) 
      ENDDO

      DO IJ=KIJS,KIJL
        DO M=MMSTART(IJ),MMSTOP(IJ)
          DO K=1,NANG
            SUM1(IJ) = SUM1(IJ) +FL1(IJ,K,M)*DFIM(M)
            SUM2(IJ) = SUM2(IJ) +COS(TH(K)-THMEAN(IJ,M))*FL1(IJ,K,M)*DFIM(M)
          ENDDO
        ENDDO
      ENDDO

      DO IJ=KIJS,KIJL
        IF (SUM1(IJ) > ZEPSILON) THEN
          R1 = SUM2(IJ)/SUM1(IJ)
          SIG_TH(IJ) = CONST_SIG*SQRT(2._JWRB*(1._JWRB-R1))
        ELSE
          SIG_TH(IJ) = 0._JWRB
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('PEAK_ANG',1,ZHOOK_HANDLE)

      END SUBROUTINE PEAK_ANG
