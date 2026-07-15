! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE W_MAXH (KIJS, KIJL, F, DEPTH, WAVNUM,           &
 &                 CMAX_F, HMAX_N, CMAX_ST, HMAX_ST, PHIST)

 ! ---------------------------------------------------------------------------- !
 !                                                                              !
 !   W_MAXH DETERMINE EXPECTED MAXIMUM CREST AND CREST-TO-TROUGH WAVE HEIGHTS.  !
 !                                                                              !
 !     FRANCESCO BARBARIOL  CNR-ISMAR   adapted from WW3 5.16 (09/2018)         !
 !     PAOLO PEZZUTTO       CNR-ISMAR   min of autocov golden search (03/2019)  !
 !     JEAN BIDLOT          ECMWF ADAPTATION (Feb. 2020)
 !                                                                              !
 !     PURPOSE:                                                                 !
 !     -------                                                                  !
 !                                                                              !
 !           DETERMINE EXPECTED MAXIMUM CREST AND CREST-TO-TROUGH WAVE HEIGHTS  !
 !           ACCORDING TO TWO DIFFERENT STATISTICAL APPROACHES:                 !
 !           1. TIME EXTREMES:                                                  !
 !              1a. CREST H., FORRISTALL (2000), NONLINEAR 2ND ORDER            !
 !              1b. WAVE H., NAESS (1985), LINEAR (ARBITRARY BANDWIDTH)         !
 !           2. SPACE-TIME EXTREMES:                                            !
 !              2a. CREST H., BENETAZZO ET AL. (2015), NONLINEAR 2ND ORDER      !
 !              2b. WAVE H., FEDELE (2012), BOCCOTTI (2000), LINEAR             !
 !           See also
 !           Barbariol et al., 2019:
 !           Maximum wave height from global model reanalysis
 !                                                                              !
! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWPCONS , ONLY : G        ,ZPI
USE YOWFRED  , ONLY : FR       ,TH       ,DELTH   ,DFIM    ,DFIMFR  , DFIMFR2, &
&                     COSTH    ,SINTH    ,WETAIL
USE YOWPARAM , ONLY : NANG     ,NFRE

USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

!     INTERFACE VARIABLES.
!     --------------------

IMPLICIT NONE
#include "aki.intfb.h"
#include "w_mode_st.intfb.h"

INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: F  !! BLOCK OF SPECTRA
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN)    :: DEPTH     !! DEPTH
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM  !! WAVE NUMBER
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT)   :: CMAX_F    !! MAXIMUM CREST H.- TIME (FORRISTALL)
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT)   :: HMAX_N    !! MAXIMUM WAVE H.- TIME (NAESS)
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT)   :: CMAX_ST   !! MAXIMUM CREST H.- SPACE-TIME (STQD)
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT)   :: HMAX_ST   !! MAXIMUM WAVE H.- SPACE-TIME (STQD)
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT)   :: PHIST     !! 1st minimum of the aotocovariance function

!    LOCAL VARIABLES.
!    ----------------

INTEGER(KIND=JWIM), PARAMETER :: MAXIT = 10

REAL(KIND=JWRB), PARAMETER :: GAMMA_E = 0.57721566_JWRB  !! EULER CONSTANT 
REAL(KIND=JWRB), PARAMETER :: GRRM1 = 2._JWRB/(1._JWRB+SQRT(5._JWRB))  !! INVERSE OF GOLDEN RATIO
REAL(KIND=JWRB), PARAMETER :: THREEHALF = 3._JWRB/2._JWRB
REAL(KIND=JWRB), PARAMETER :: XNWVP = 100._JWRB  !! NUMBER OF WAVE PERIODS TO SET WMDUR_ST (if constant value is not used)
REAL(KIND=JWRB), PARAMETER :: TOL = 0.01_JWRB ! GOLDEN SEARCH TOLERANCE

INTEGER(KIND=JWIM) :: IJ, IT, M, K
INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL) :: K_THMAX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JWRB) :: Z0, RNW, ALFA, BETA, URSN, STEEP, WNUM1
REAL(KIND=JWRB) :: XK, XK2, AXYT, RN3, RN2, RN1, SQRTEM, TMIN, TMAX, WVLMIN, XX
REAL(KIND=JWRB) :: ZEPSILON
REAL(KIND=JWRB) :: DELT25, EM, XK2_DFIM, XK_ZPI_DFIM
REAL(KIND=JWRB), DIMENSION(4) :: TLGS, ACFS
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: HS, RLX, RLY, AXT, AYT, AXY
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: WMDX, WMDY ! SPACE TIME EXTREME OVER WMDX x WMDY m**2 AREA
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: WMDUR_TD ! TIME EXTREME OVER WMDUR_TD for TIME DOMAIN in sec.
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: WMDUR_ST ! TIME EXTREME OVER WMDUR_ST for SPACE TIME in sec.
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: FMAX
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ACF, T1, T2, EMEAN 
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP, TEMP_X, TEMP_Y, TEMP_X2, TEMP_Y2, TEMP_XY
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: RNI, RMU 
REAL(KIND=JWRB), DIMENSION(NFRE) :: OMEGA
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG) :: CX, CY, CX2, CY2, CXCY
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NFRE) :: TEMPDFIM

! ----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('W_MAXH',0,ZHOOK_HANDLE)

! ----------------------------------------------------------------------------

!     0. COMPUTE SPECTRAL PARAMETER BY INTEGRATION OF SPECTRA.
!        -----------------------------------------------------

ZEPSILON = 10._JWRB*EPSILON(ZEPSILON)
OMEGA(:) = ZPI*FR(:)
TMIN = 1._JWRB/FR(NFRE)
TMAX = 1._JWRB/FR(1)
WVLMIN = G / (ZPI*FR(NFRE)**2)

!!! see below, we are now using wave state dependant values
!!!WMDX(KIJS:KIJL) = 100._JWRB
!!!WMDY(KIJS:KIJL) = 100._JWRB

!!! use a constant duration for time domain
WMDUR_TD(KIJS:KIJL) = 1200._JWRB
!!! but not for space time
!!!WMDUR_ST(KIJS:KIJL) = 1200._JWRB

! INITIALIZE VARIABLES FOR INTEGRATION

ACF(KIJS:KIJL) = 0._JWRB
T1(KIJS:KIJL) = 0._JWRB
T2(KIJS:KIJL) = 0._JWRB
RLX(KIJS:KIJL) = 0._JWRB
RLY(KIJS:KIJL) = 0._JWRB
AXT(KIJS:KIJL) = 0._JWRB
AXY(KIJS:KIJL) = 0._JWRB
AYT(KIJS:KIJL) = 0._JWRB
EMEAN(KIJS:KIJL) = 0._JWRB

! FIND THE DIRECTION WHICH CORRESPONDS TO THE SPECTRAL MAXIMUM
!! in radians!
K_THMAX(KIJS:KIJL)=1
FMAX(KIJS:KIJL)=0._JWRB
DO M = 1,NFRE
  DO K = 1,NANG
    DO IJ = KIJS, KIJL 
      IF (F(IJ,K,M) > FMAX(IJ) ) THEN
         FMAX(IJ) = F(IJ,K,M)
         K_THMAX(IJ)=K
      ENDIF
    ENDDO
  ENDDO
ENDDO

! COMPUTE MEAN PERIODS (WITHOUT TAIL CONTRIBUTION)
! COMPUTE MEAN WAVE AND CREST LENGTH (WITHOUT TAIL CONTRIBUTION, WRT PEAK DIR.)
! COMPUTE IRREGULARITY PARAMETERS (WITHOUT TAIL CONTRIBUTION, WRT PEAK DIR.)

DO K = 1,NANG
  DO IJ = KIJS, KIJL 
    CX(IJ,K) = COSTH(K)*COSTH(K_THMAX(IJ))+SINTH(K)*SINTH(K_THMAX(IJ))
    CY(IJ,K) = SINTH(K)*COSTH(K_THMAX(IJ))-COSTH(K)*SINTH(K_THMAX(IJ))
    CX2(IJ,K) = CX(IJ,K)**2
    CY2(IJ,K) = CY(IJ,K)**2
    CXCY(IJ,K) = CX(IJ,K)*CY(IJ,K)
  ENDDO
ENDDO


DO M = 1,NFRE
   TEMP(KIJS:KIJL) = 0._JWRB
   TEMP_X(KIJS:KIJL) = 0._JWRB
   TEMP_Y(KIJS:KIJL) = 0._JWRB
   TEMP_X2(KIJS:KIJL) = ZEPSILON 
   TEMP_Y2(KIJS:KIJL) = ZEPSILON 
   TEMP_XY(KIJS:KIJL) = 0._JWRB

   DO K = 1, NANG
     DO IJ = KIJS, KIJL 
       TEMP(IJ) = TEMP(IJ) + F(IJ,K,M)
       TEMP_X(IJ) = TEMP_X(IJ) + F(IJ,K,M)*CX(IJ,K)
       TEMP_Y(IJ) = TEMP_Y(IJ) + F(IJ,K,M)*CY(IJ,K)
       TEMP_X2(IJ) = TEMP_X2(IJ) + F(IJ,K,M)*CX2(IJ,K)
       TEMP_Y2(IJ) = TEMP_Y2(IJ) + F(IJ,K,M)*CY2(IJ,K)
       TEMP_XY(IJ) = TEMP_XY(IJ) + F(IJ,K,M)*CXCY(IJ,K)
     ENDDO
   ENDDO

   DO IJ = KIJS, KIJL
     XK = WAVNUM(IJ,M)
     XK2 = XK**2
     T1(IJ) = T1(IJ) + TEMP(IJ)*DFIMFR(M)
     T2(IJ) = T2(IJ) + TEMP(IJ)*DFIMFR2(M)
     EMEAN(IJ) = EMEAN(IJ) + TEMP(IJ)*DFIM(M)
     TEMPDFIM(IJ,M) = TEMP(IJ)*DFIM(M)
     XK2_DFIM = XK2*DFIM(M)
     RLX(IJ) = RLX(IJ) + TEMP_X2(IJ)*XK2_DFIM
     RLY(IJ) = RLY(IJ) + TEMP_Y2(IJ)*XK2_DFIM
     AXY(IJ) = AXY(IJ) + TEMP_XY(IJ)*XK2_DFIM
     XK_ZPI_DFIM = XK*ZPI*DFIMFR(M)
     AXT(IJ) = AXT(IJ) + TEMP_X(IJ)*XK_ZPI_DFIM
     AYT(IJ) = AYT(IJ) + TEMP_Y(IJ)*XK_ZPI_DFIM
   ENDDO
ENDDO

! ADD HIGH FREQUENCY TAIL CONTRIBUTION AND COMPUTE HS
DELT25 = WETAIL*FR(NFRE)*DELTH
DO IJ = KIJS, KIJL
  IF ( EMEAN(IJ) > ZEPSILON ) THEN
   EM = EMEAN(IJ) + DELT25*TEMP(IJ)
   HS(IJ) = 4.0_JWRB*SQRT(EM)
  ELSE
    HS(IJ) = 0.0_JWRB
  ENDIF
ENDDO

WHERE (EMEAN > ZEPSILON)
   AXY = MIN(AXY/SQRT(RLX*RLY), 1._JWRB)
   AXT = MIN(AXT/(ZPI*SQRT(RLX*T2)), 1._JWRB)
   AYT = MIN(AYT/(ZPI*SQRT(RLY*T2)), 1._JWRB)
   RLX = ZPI*SQRT(EMEAN/RLX)
   RLY = ZPI*SQRT(EMEAN/RLY)
   RNI = SQRT(MAX(EMEAN*T2/T1**2 - 1._JWRB,ZEPSILON))
   RMU = (ZPI*T1)**2*(1._JWRB-RNI+RNI**2)/(G*EMEAN**THREEHALF)
   T1 = MIN(MAX(EMEAN/T1,TMIN),TMAX)
   T2 = MIN(MAX(SQRT(EMEAN/T2),TMIN),TMAX)
!!!   WMDUR_TD = XNWVP*T2
   WMDX = MAX(RLX,WVLMIN)
   WMDY = MAX(RLY,WVLMIN)
   WMDUR_ST = XNWVP*T2
ELSEWHERE
   AXY = 0._JWRB
   AXT = 0._JWRB
   AYT = 0._JWRB
   RLX = 1._JWRB
   RLY = 1._JWRB
   RNI = 0._JWRB
   RMU = 0._JWRB
   T1 = TMIN
   T2 = TMIN
!!!   WMDUR_TD = XNWVP*T2
   WMDX = MAX(RLX,WVLMIN)
   WMDY = MAX(RLY,WVLMIN)
   WMDUR_ST = XNWVP*T2
END WHERE 

!  MIN OF AUTOCOVARIANCE FUNCTION  VIA GOLDEN RATIO SEARCH

DO IJ = KIJS, KIJL
  TLGS(1) = 0.3_JWRB*T2(IJ)
  TLGS(4) = 1.3_JWRB*T2(IJ)
  TLGS(2) = TLGS(4) - (TLGS(4) - TLGS(1))*GRRM1
  TLGS(3) = TLGS(1) + (TLGS(4) - TLGS(1))*GRRM1
  ACFS(1) = SUM( COS( OMEGA(:)*TLGS(1) ) * TEMPDFIM(IJ,:) )
  ACFS(4) = SUM( COS( OMEGA(:)*TLGS(4) ) * TEMPDFIM(IJ,:) )
  ACFS(2) = SUM( COS( OMEGA(:)*TLGS(2) ) * TEMPDFIM(IJ,:) )
  ACFS(3) = SUM( COS( OMEGA(:)*TLGS(3) ) * TEMPDFIM(IJ,:) )

  DO IT = 1,MAXIT
    IF (ACFS(2) < ACFS(3)) THEN
      ACF(IJ) = ACFS(2)
      TLGS(4) = TLGS(3)
      ACFS(4) = ACFS(3)
      TLGS(3) = TLGS(2)
      ACFS(3) = ACFS(2)
      TLGS(2) = TLGS(4) - (TLGS(4) - TLGS(1))*GRRM1
      ACFS(2) = SUM( COS( OMEGA(:)*TLGS(2) ) * TEMPDFIM(IJ,:) )
    ELSE
      ACF(IJ) = ACFS(3)
      TLGS(1) = TLGS(2)
      ACFS(1) = ACFS(2)
      TLGS(2) = TLGS(3)
      ACFS(2) = ACFS(3)
      TLGS(3) = TLGS(1) + (TLGS(4) - TLGS(1))*GRRM1
      ACFS(3) = SUM( COS( OMEGA(:)*TLGS(3) ) * TEMPDFIM(IJ,:) )
    ENDIF
    IF ( (ABS(TLGS(4)-TLGS(1))) < (TOL*(ABS(TLGS(2))+ABS(TLGS(3)))) ) EXIT
  ENDDO
ENDDO

!     1. COMPUTE OUTPUT VARIABLES

DO IJ = KIJS, KIJL 

   IF ( EMEAN(IJ) > ZEPSILON ) THEN
      SQRTEM = 0.25_JWRB*HS(IJ)
      WNUM1 = AKI(ZPI/T1(IJ),DEPTH(IJ))
      STEEP = ZPI*HS(IJ)/(G*T1(IJ)**2)
      URSN = HS(IJ)/(WNUM1**2*DEPTH(IJ)**3)
      ALFA = 0.3536_JWRB+0.2568_JWRB*STEEP+0.08_JWRB*URSN
      BETA = 2._JWRB-1.7912_JWRB*STEEP-0.5302_JWRB*URSN+0.284_JWRB*URSN**2
      RNW = WMDUR_TD(IJ)/T2(IJ)
      Z0 = LOG(RNW)

!     1a. DETERMINE EXPECTED MAXIMUM CREST HEIGHT - TIME (FORRISTALL).
!     ----------------------------------------------------------------

      CMAX_F(IJ) = ALFA*Z0**(1._JWRB/BETA)*(1._JWRB+GAMMA_E/(BETA*Z0))*HS(IJ)

!     1b. DETERMINE EXPECTED MAXIMUM WAVE HEIGHT - TIME (NAESS).
!     ----------------------------------------------------------
      PHIST(IJ) = MIN(ACF(IJ)/EMEAN(IJ), 1._JWRB)
      HMAX_N(IJ) = 0.5_JWRB*SQRT(1._JWRB-PHIST(IJ))*SQRT(Z0)*(1._JWRB+0.5_JWRB*GAMMA_E/Z0)*HS(IJ)



!     2a. DETERMINE EXPECTED MAXIMUM CREST HEIGHT - SPACE-TIME (STQD).
!     ----------------------------------------------------------------

      AXYT = SQRT(1._JWRB+2._JWRB*AXT(IJ)*AXY(IJ)*AYT(IJ)-AXT(IJ)**2-AXY(IJ)**2-AYT(IJ)**2)
      RN3 = ZPI*WMDX(IJ)*WMDY(IJ)*WMDUR_ST(IJ)*AXYT/(RLX(IJ)*RLY(IJ)*T2(IJ))
      RN2 = SQRT(ZPI)*(WMDX(IJ)*WMDUR_ST(IJ)/(RLX(IJ)*T2(IJ))*SQRT(1._JWRB-AXT(IJ)**2) &
&          +WMDX(IJ)*WMDY(IJ)/(RLX(IJ)*RLY(IJ))*SQRT(1._JWRB-AXY(IJ)**2) &
&          +WMDY(IJ)*WMDUR_ST(IJ)/(RLY(IJ)*T2(IJ))*SQRT(1._JWRB-AYT(IJ)**2))
      RN1 = WMDX(IJ)/RLX(IJ) + WMDY(IJ)/RLY(IJ) + WMDUR_ST(IJ)/T2(IJ)

!      IF (WMDX(IJ) /= 0._JWRB .AND. WMDY(IJ)/= 0._JWRB .AND. WMDUR_ST(IJ) /= 0._JWRB) THEN
         Z0 = W_MODE_ST(RN3, RN2, RN1)
!      ELSE IF (WMDUR_ST(IJ) == 0._JWRB) THEN
!         Z0 = SQRT(2.JWRB*LOG(REAL(RN2))+LOG(2._JWRB*LOG(REAL(RN2))+LOG(2._JWRB*LOG(REAL(RN2)))))
!      ELSE IF (WMDX(IJ) == 0._JWRB .AND. WMDY(IJ) == 0._JWRB) THEN
!         Z0 = SQRT(2._JWRB*LOG(REAL(RN1)))
!      END IF

      XX = 1._JWRB/(Z0-(2._JWRB*RN3*Z0+RN2)/(RN3*Z0**2+RN2*Z0+RN1))

      CMAX_ST(IJ) = ((Z0+0.5_JWRB*RMU(IJ)*Z0**2)+GAMMA_E*((1._JWRB+RMU(IJ)*Z0)*XX) ) *SQRTEM

!     2b. DETERMINE EXPECTED MAXIMUM WAVE HEIGHT - SPACE-TIME (STQD).
!     ---------------------------------------------------------------

      HMAX_ST(IJ) = (Z0+GAMMA_E*XX)*SQRT(2._JWRB*(1._JWRB-PHIST(IJ))) *SQRTEM

   ELSE
      CMAX_F(IJ) = 0._JWRB
      HMAX_N(IJ) = 0._JWRB
      CMAX_ST(IJ) = 0._JWRB
      HMAX_ST(IJ) = 0._JWRB
   ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('W_MAXH',1,ZHOOK_HANDLE)

END SUBROUTINE W_MAXH
