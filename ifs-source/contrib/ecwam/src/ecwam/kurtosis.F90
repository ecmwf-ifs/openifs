! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE KURTOSIS(KIJS, KIJL, FL1,                        &
     &                    DEPTH,                                  &
     &                    C3, C4, BF2, QP, HMAX, TMAX,            &
     &                    ETA_M, R, XNSLC, SIG_TH, EPS, XNU)
 
!***  *KURTOSIS*   DETERMINES SKEWNESS C3, KURTOSIS C4, THE SQUARE OF THE
!                  BENJAMIN-FEIR INDEX BFI2, GODA'S PEAKEDNESS PARAMETER QP, 
!                  THE PARAMETER R, MEASURING THE RATIO OF DIRECTIONAL WIDTH 
!                  AND FREQUENCY WIDTH, NORMALIZED MAXIMUM WAVE HEIGHT 
!                  AND PERIOD AT MAXIMUM HEIGHT.
 
!     PETER JANSSEN       JULY 2007.
!     PETER JANSSEN       APRIL 2014. INCLUDES SKEWNESS EFFECTS, WHILE 
!                         NUMBER OF WAVES IS DETERMINED BY EWING (1973)
!     F. VANA             MARCH 2015 SINGLE PRECISION SUPPORT
!     JEAN BIDLOT/
!     PETER JANSSEN       NOVEMBER 2017. INTRODUCES SHALLOW WATER
!                         EFFECTS AND IMPROVED DYNAMICS PART OF
!                         KURTOSIS

!     PURPOSE.
!     --------
 
!           DETERMINATION OF SKEWNESS, KURTOSIS, B-F INDEX , GODA QP, NW,
!                              HMAXN AND TMAX
!     INTERFACE.
!     ----------
!           *CALL*  *KURTOSIS(KIJS, KIJL, FL1,
!                             DEPTH,
!                             C3, C4, BF2, QP, HMAX, TMAX,
!                             ETA_M, R, XNSLC, SIG_YH, EPS, XNU)
!                      INPUT:
!                           *FL1*    - 2-DIMENSIONAL SPECTRUM
!                           *KIJS*   - STARTING INDEX
!                           *KIJL*   - LAST INDEX
!                           *DEPTH*  - WATER DEPTH
!                      OUTPUT: 
!                           *C3*    - SKEWNESS 
!                           *C4*    - KURTOSIS
!                           *BF*    - SQUARE OF BENJAMIN-FEIR INDEX
!                           *QP*    - GODA'S QUALITY FACTOR
!                           *HMAX*  - MAXIMUM WAVE HEIGHT
!                           *TMAX*  - MAXIMUM WAVE PERIOD
!                           *ETA_M* - WAVE-INDUCED MEAN SURFACE ELEVATION
!                           *R*     - SPECTRAL WIDTH INDEX
!                           *XNSLC* - NUMBER OF EVENTS
!                           *SIG_TH*- RELATIVE WIDTH IN DIRECTION
!                           *EPS*   - WAVE STEEPNESS
!                           *XNU*   - RELATIVE SPECTRAL WIDTH


 
!     METHOD.
!     -------
!            JULY 2007:
!            ---------
 
!            USED NUMERICAL SIMULATIONS OF THE NONLINEAR SCHROEDINGER
!            EQUATION TO OBTAIN THE DEPENDENCE OF KURTOSIS ON THE
!            SQUARE OF THE BENJAMIN-FEIR INDEX, BF2, AND THE WIDTH OF
!            THE DIRECTIONAL DISTRIBUTION, SIG_TH. FITTING THE NUMERICAL 
!            RESULTS ONE FINDS (ONORATO, MORI AND JANSSEN, 2007)
                  
!              C4 = A_MORI/SIG_TH * BF2 * PI/(3*SQRT(3))
 
!            WHERE A_MORI=0.031. NOTE THAT FOR NARROW SPECTRA, I.E.
!            SIG_TH => A_MORI, THE ONE-DIMENSIONAL RESULT, PREVIOUSLY
!            IMPLEMENTED, IS OBTAINED.
 
!            IN ADDITION, IT IS WELL-KNOWN THAT WHEN A WAVE GROUP
!            APPROACHES SHALLOW WATER THE NONLINEAR FOCUSSING 
!            EFFECT DISAPPEARS FOR K*D=1.363 AND WHEN K*D BECOMES
!            LESS THAN THIS CRITICAL VALUE EVEN DEFOCUSSING OCCURS.
!            AS A CONSEQUENCE, FOR K*D>1.363 C4 > 0 WHILE IN THE
!            OPPOSITE CASE C4 < 0. HERE, FOLLOWING THE WORK
!            OF JANSSEN & ONORATO (2007) WE HAVE PARAMETRIZED THIS SHALLOW
!            WATER EFFECT BY MEANS OF THE FUNCTION TRANSF_BFI.
 
!            APRIL 2014:
!            ----------
 
!            1) FORMULATE DYNAMIC KURTOSIS C4 IN TERMS OF THE DIMENSIONLESS
!               NUMBER R = 0.5*(SIG_TH/SIG_OM)**2 (SEE TECH MEMO 588). 
!               RESULTS IN GOOD AGREEMENT WITH  MORI FORMULATION FOR
!                  C4 = XJ*BFI**2
!               WITH
!                  XJ = C4_CONST/SQRT(1+7.*R),
!               WHERE
!                  C4_CONST = PI/(3.*SQRT(3.))
!            2) BOUND WAVES HAVE ALSO FINITE SKEWNESS WHICH AFFECTS ENVELOPE
!               WAVE HEIGHT DISTRIBUTION (SEE MEMO ON UPDATES)
!            3) FOR MAXIMUM WAVE HEIGHT THE NUMBER NW OF INDEPENDENT EVENTS
!               IS REQUIRED. THUS FAR NUMBER OF INDEPENDENT EVENTS WAS OBTAINED 
!               BY SAMPLING WITH THE PEAK FREQUENCY. THIS IS NOT CORRECT. IT
!               IS MORE APPROPRIATE TO USE AS ESTIMATE OF EVENTS THE NUMBER OF 
!               WAVE GROUPS IN THE TIME SERIES. THIS IS OBTAINED USING JOHN
!               EWINGS (1973) RESULT ON THE LENGTH OF A WAVE GROUP REFERRED TO
!               SIGNIFICANT WAVE HEIGHT.
 
!            NOVEMBER 2017(FOR DETAILS SEE TECH MEMO 813):
!            ---------------------------------------------
 
!            1) FORMULATE DYNAMIC KURTOSIS C4 IN TERMS OF THE DIMENSIONLESS
!               NUMBER R = 0.5*(SIG_TH/SIG_OM)**2*f(K*D) (SEE TECH MEMO 813). 
!               NOW, DYNAMIC KURTOSIS IS OBTAINED FROM A PARAMETRIZATION OF 
!               THE EXACT SOLUTION FOR MAXIMUM KURTOSIS OBTAINED FROM NLS 
!               EQUATION. FORMULA FOR KURTOSIS READS
!                  C4 = XJ*BFI**2
!               WITH
!                  XJ = C4_CONST*R_0*(1.1-R)/(R+R_0),
!               WHERE
!                  C4_CONST = PI/(3.*SQRT(3.))
!               AND 
!                  R_0 = 7.44*SQRT(3)/(4*PI**3).
 
!            2) SIMPLE PARAMETRIZATIONS FOR THE DEPTH-DEPENDENCE OF THE BOUND
!               WAVE PART OF SKEWNESS AND KURTOSIS HAVE BEEN INTRODUCED. 
!               EFFECT OF SPECTRAL WIDTH HAS BEEN INTRODUCED THROUGH AN EFFECTIVE 
!               WAVENUMBER OBTAINED FROM THE INVERSE OF THE DISPERSION
!               RELATION AT 0.9*THE MEAN ANGULAR FREQUENCY.
 

!            3) THE MAXIMUM WAVE HEIGHT DISTRIBUTION IS NOW BASED ON THE WORK 
!               OF NAESS (1982). THIS PDF INTRODUCES AS THE NUMBER OF INDEPENDENT 
!               EVENTS THE NUMBER OF WAVE GROUPS (EWING, 1973) AT THE 
!               SIGNIFICANT LEVEL, CALLED XNSLC, IN A NATURAL WAY.
 
!            4) DEVIATIONS FROM NORMALITY ARE GIVEN BY THE EDGEWORTH 
!               DISTRIBUTION. HOWEVER, THIS SMALL STEEPNESS EXPANSION IS 
!               NOT UNIFORMLY VALID. IN PRACTICE THE TAIL OF THE DISTRIBUTION IS 
!               EXPONENTIAL. ACCOMODATE FOR THIS BY USING THE STRETCHED 
!               EXPONENTIAL DISTRIBUTION 
 
!                  P(E)= EXP{-AA+SQRT(AA**2+BB*E)
 
!               WITH E = WAVE ENERGY 2*H**2, AND H IS ENVELOPE WAVE HEIGHT 
!               NORMALIZED WITH SIGNIFICANT WAVE HEIGHT. THE PARAMETERS AA
!               AND BB=2*(AA+1) DEPEND ON SKEWNESS AND KURTOSIS AS OBTAINED
!               FROM A FIT OF THE STRETCHED EXPONENTIAL DISTRIBUTION TO THE 
!               EDGEWORTH DISTRIBUTION. THE CALCULATION OF THE EXPECTED
!               VALUE OF MAXIMUM ENVELOPE WAVE HEIGHT IS NOW OBTAINED USING THE
!               STRETCHED EXPONENTIAL DISTRIBUTION.  
                
!     EXTERNALS.
!     ----------
!             AKI,TRANSF_BFI,H_MAX,PEAK_ANG
 
!     REFERENCES:
!     ---------
 
!     NONLINEAR FOUR WAVE INTERACTIONS AND FREAK WAVES
!     PETER A.E.M. JANSSEN, JPO, 863-884, APRIL 2003
 
!     ON KURTOSIS AND OCCURRENCE PROBABILITY OF FREAK WAVES
!     NOBUHITO MORI AND PETER A.E.M. JANSSEN, JPO, 1471-1483, JULY 2006 
 
!     THE INTERMEDIATE WATER DEPTH LIMIT OF THE ZAKHAROV EQUATION AND 
!     CONSEQUENCES FOR WAVE PREDICTION
!     PETER A.E.M. JANSSEN AND MIGUEL ONORATO, ACCEPTED FOR PUBLICATION
!     IN JPO, 2007
 
!     EFFECTS OF DIRECTIONALITY ON THE NONLINEAR FOCUSSING IN RANDOM SEAS. !
!     MIGUEL ONORATO, NOBUHITO MORI AND PETER A.E.M. JANSSEN
!     IN PREPARATION, 2007
 
!     ON THE EXTENSION OF THE FREAK WAVE WARNING SYSTEM AND ITS VERIFICATION.
!     PETER A.E.M JANSSEN AND J.-R. BIDLOT, ECMWF TECH MEMO 588.
 
!     FURTHER UPDATES TO THE FREAK WAVE WARNING SYSTEM. PETER A.E.M. 
!     JANSSEN AND J.-R. BIDLOT, TO BE PUBLISHED AS ECMWF TECH MEMO
 
!-----------------------------------------------------------------------
 
!     MODULES:
!     --------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : G        ,PI       ,ZPI     ,ZPISQRT
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH   ,DFIMOFR ,      &
     &             DFIMFR   ,DFIMFR2,                                   &
     &             WETAIL   ,WP1TAIL ,WP2TAIL  ,QPTAIL  ,FRTAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
 
!-----------------------------------------------------------------------

      IMPLICIT NONE
#include "aki.intfb.h"
#include "peak_ang.intfb.h"
#include "stat_nl.intfb.h"
#include "transf_bfi.intfb.h"
#include "h_max.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: C3, C4, BF2, QP, HMAX, TMAX
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: ETA_M, R, XNSLC, SIG_TH, EPS, XNU

      INTEGER(KIND=JWIM) :: IJ, M, K

      REAL(KIND=JWRB), PARAMETER :: QPMIN = 0.5_JWRB
      REAL(KIND=JWRB), PARAMETER :: QPMAX = 15.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: BF2MIN = -5._JWRB
      REAL(KIND=JWRB), PARAMETER :: BF2MAX = 5._JWRB
      REAL(KIND=JWRB), PARAMETER :: FLTHRS = 0.4_JWRB
      REAL(KIND=JWRB), PARAMETER :: SQRT2 = SQRT(2._JWRB)

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: DELT25, COEF_FR1, COEF_FR2, DELT2
      REAL(KIND=JWRB) :: CONST_SIG_SQRTPIM1, CONST_OM_ZPI, OM_MEAN
      REAL(KIND=JWRB) :: TRANS, DUR, TAU, ZFAC, ZEPS, HS
      REAL(KIND=JWRB) :: ZEPSILON, ZSQREPSILON, FRMAX, FRMIN
 
      REAL(KIND=JWRB), DIMENSION(NFRE) :: FAC4
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: HMAXN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: SUM0,SUM1,SUM2,SUM4,SUM40,SUM6
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: XKP,SIG_OM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: SIG_HM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: F_M,OM_UP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: AA,BB,C4_DYN,C4_B
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: FFMAX
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: FF 
     
!-----------------------------------------------------------------------
 
      IF (LHOOK) CALL DR_HOOK('KURTOSIS',0,ZHOOK_HANDLE)

      ZEPSILON=10._JWRB*EPSILON(ZEPSILON)
      ZSQREPSILON=SQRT(ZEPSILON)

!!!      CONST_SIG_SQRTPIM1 = 0.55_JWRB/SQRT(PI)
      CONST_SIG_SQRTPIM1 = 1.0_JWRB/SQRT(PI)
      CONST_OM_ZPI  = 0.89_JWRB*ZPI
      FRMAX = FR(NFRE)
      FRMIN = FR(1)

!***  2. DETERMINE ANGULAR WIDTH PARAMETER SIG_TH.
!     -------------------------------------------
 
      CALL PEAK_ANG(KIJS, KIJL, FL1, XNU, SIG_TH)

!     COMPUTES THE DIFFERENT MOMENTS 

      DO M=1,NFRE
        K=1
        DO IJ=KIJS,KIJL
          FF(IJ,M) = FL1(IJ,K,M)
        ENDDO   
        DO K=2,NANG
          DO IJ=KIJS,KIJL
            FF(IJ,M) = FF(IJ,M)+FL1(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO 

      DO IJ=KIJS,KIJL
        FFMAX(IJ)=FF(IJ,1)
      ENDDO
      DO M=2,NFRE
        DO IJ=KIJS,KIJL
          FFMAX(IJ)=MAX(FFMAX(IJ),FF(IJ,M))
        ENDDO
      ENDDO 

      DO IJ=KIJS,KIJL
        SUM0(IJ)= ZEPSILON 
        SUM1(IJ)= 0._JWRB
        SUM2(IJ)= 0._JWRB
        SUM6(IJ)= 0._JWRB          
      ENDDO

      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          SUM0(IJ) = SUM0(IJ)+FF(IJ,M)*DFIM(M)
          SUM1(IJ) = SUM1(IJ)+FF(IJ,M)*DFIMFR(M)
          SUM2(IJ) = SUM2(IJ)+FF(IJ,M)*DFIMFR2(M)
          SUM6(IJ) = SUM6(IJ)+FF(IJ,M)*DFIMOFR(M)
        ENDDO
      ENDDO
!     ADD HIGH FREQUENCY TAIL CONTRIBUTIONS
      DELT25   = WETAIL*FR(NFRE)*DELTH
      COEF_FR1 = WP1TAIL*DELTH*FR(NFRE)**2
      COEF_FR2 = WP2TAIL*DELTH*FR(NFRE)**3
      DELT2    = FRTAIL*DELTH
      DO IJ=KIJS,KIJL
        SUM0(IJ) = SUM0(IJ)+DELT25*FF(IJ,NFRE)
        SUM1(IJ) = SUM1(IJ)+COEF_FR1*FF(IJ,NFRE)
        SUM2(IJ) = SUM2(IJ)+COEF_FR2*FF(IJ,NFRE)
        SUM6(IJ) = SUM6(IJ)+DELT2*FF(IJ,NFRE)
      ENDDO


      DO M=1,NFRE
        FAC4(M) = 2._JWRB*DELTH*DFIMFR(M)
      ENDDO
      DO IJ=KIJS,KIJL
        SUM40(IJ)= ZSQREPSILON 
        SUM4(IJ)= 0._JWRB
        FFMAX(IJ)=FLTHRS*FFMAX(IJ)
      ENDDO
      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          IF (FF(IJ,M) > FFMAX(IJ)) THEN
            SUM40(IJ)= SUM40(IJ)+FF(IJ,M)*DFIM(M)
            SUM4(IJ) = SUM4(IJ)+FF(IJ,M)**2*FAC4(M)
          ENDIF
        ENDDO
      ENDDO

!***  3. DETERMINE EPS AND BFI^2.
!     --------------------------
      DO IJ=KIJS,KIJL
        IF (SUM1(IJ) > ZSQREPSILON .AND. SUM0(IJ) > ZEPSILON) THEN

          F_M(IJ) = MAX(MIN(SUM1(IJ)/SUM0(IJ),FRMAX),FRMIN)
          QP(IJ) = MAX(MIN(SUM4(IJ)/SUM40(IJ)**2,QPMAX),QPMIN)
          SIG_OM(IJ) = CONST_SIG_SQRTPIM1/QP(IJ)

          OM_MEAN = CONST_OM_ZPI*MAX(MIN(SUM0(IJ)/SUM6(IJ),FRMAX),FRMIN)
          XKP(IJ) = AKI(OM_MEAN,DEPTH(IJ))
          EPS(IJ) = XKP(IJ)*SQRT(SUM0(IJ))

          TRANS  = TRANSF_BFI(XKP(IJ),DEPTH(IJ),XNU(IJ),SIG_TH(IJ))
          BF2(IJ)= 2._JWRB*TRANS*(EPS(IJ)/MAX(SIG_OM(IJ),ZEPSILON))**2
          BF2(IJ) = MAX(MIN(BF2(IJ),BF2MAX),BF2MIN)
        ELSE
          F_M(IJ)    = 0._JWRB
          QP(IJ)     = 0._JWRB
          SIG_OM(IJ) = 0._JWRB
          OM_MEAN    = CONST_OM_ZPI*FRMAX
          XKP(IJ)    = OM_MEAN**2/G
          EPS(IJ)    = 0._JWRB
          BF2(IJ)    = 0._JWRB
        ENDIF 
      ENDDO
     
!***  4. DETERMINE C3 AND C4.
!     ----------------------
 
      CALL STAT_NL(KIJS, KIJL,                                         &
     &             SUM0, XKP, BF2, XNU, SIG_TH, DEPTH,                 &
     &             C3, C4, ETA_M, R, C4_B, C4_DYN)

  
!***  5. DETERMINE HMAXN AND TMAX.
!     ---------------------------
 
 
!     DETERMINE NUMBER OF DEGREES OF FREEDOM

!     DURATION IN SECONDS: 
      DUR = 1200._JWRB

      ZFAC = 2._JWRB*ZPI/ZPISQRT
      DO IJ=KIJS,KIJL
        IF (F_M(IJ) > 0._JWRB) THEN
          OM_UP(IJ) = ZFAC*XNU(IJ)*F_M(IJ)

          XNSLC(IJ) = REAL(NINT(DUR*OM_UP(IJ)),JWRB)
        ELSE
          XNSLC(IJ) = 0._JWRB 
        ENDIF
      ENDDO

      CALL H_MAX(C3,C4,XNSLC,KIJS,KIJL,AA,BB,HMAXN,SIG_HM)

      DO IJ=KIJS,KIJL
        IF (SUM1(IJ) > ZEPSILON .AND. HMAXN(IJ) > ZEPSILON) THEN
          ZEPS = XNU(IJ)/(SQRT2*HMAXN(IJ))
          TMAX(IJ) = 1._JWRB+0.5_JWRB*ZEPS**2+0.75_JWRB*ZEPS**4
          TAU = SUM0(IJ)/SUM1(IJ)
          TMAX(IJ) = TAU*TMAX(IJ)
        ELSE
          TMAX(IJ) = 0._JWRB
        ENDIF
      ENDDO

      DO IJ=KIJS,KIJL
        IF (SUM0(IJ) > 0._JWRB) THEN
          HS = 4._JWRB*SQRT(SUM0(IJ))
          HMAX(IJ) = HMAXN(IJ)*HS
        ELSE
          HMAX(IJ) = 0._JWRB
        ENDIF
      ENDDO 

      DO IJ=KIJS,KIJL
        IF (SUM0(IJ) <= 0._JWRB) THEN
          XNU(IJ) = 0._JWRB
        ENDIF
      ENDDO 

      IF (LHOOK) CALL DR_HOOK('KURTOSIS',1,ZHOOK_HANDLE)
     
      END SUBROUTINE KURTOSIS
