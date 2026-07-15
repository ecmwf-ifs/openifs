! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE STAT_NL(KIJS, KIJL,                                    &
     &                   XM0, XK0, BF2, XNU, SIG_TH, DPTH,              &
     &                   C3, C4, ETA_M, R, C4_B, C4_DYN)
 
!***  DETERMINE SKEWNESS AND ENVELOPE KURTOSIS FOR A NARROW-BAND WAVE 
!     TRAIN IN ARBITRARY DEPTH.
 
!     PURPOSE:
!     -------
     
!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
 
!     XM0            REAL         WAVE VARIANCE
!     XK0            REAL         PEAK WAVE NUMBER
!     DEPTH          REAL         DEPTH
 
!     OUTPUT:
!     ------
 
!     C3             REAL         SKEWNESS
!     C4             REAL         ENVELOPE KURTOSIS
!     ETA_M          REAL         WAVE-INDUCED MEAN SURFACE ELEVATION
!     R              REAL         SPECTRAL WIDTH INDEX
!     XNU            REAL         RELATIVE SPECTRAL WIDTH
!     SIG_TH         REAL         RELATIVE WIDTH in DIRECTION

 
!     AUTHOR:
!     ------
 
!     P.A.E.M. JANSSEN, NOVEMBER 2017
 
!----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : G        ,PI       ,DKMAX
      USE YOWSHAL  , ONLY : XKDMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

!----------------------------------------------------------------------

      IMPLICIT NONE 

#include "transf_r.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: XM0, XK0, BF2, XNU, SIG_TH, DPTH
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(OUT) :: C3, C4, ETA_M, R, C4_B, C4_DYN

      REAL(KIND=JWRB), PARAMETER :: EPS = 0.0001_JWRB
      REAL(KIND=JWRB), PARAMETER :: RMIN = 0._JWRB
      REAL(KIND=JWRB), PARAMETER :: RMAX = 16.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: C3MIN = 0._JWRB
      REAL(KIND=JWRB), PARAMETER :: C3MAX = 0.25_JWRB
      REAL(KIND=JWRB), PARAMETER :: C4MIN = -0.25_JWRB
      REAL(KIND=JWRB), PARAMETER :: C4MAX = 0.25_JWRB

      REAL(KIND=JWRB), PARAMETER :: CONST_C3 = 1.12_JWRB*2._JWRB
      REAL(KIND=JWRB), PARAMETER :: CONST_C4 = 0.93_JWRB*8._JWRB  
      REAL(KIND=JWRB), PARAMETER :: C3DELTA_ADJ = 0.9_JWRB
      REAL(KIND=JWRB), PARAMETER :: SQRT3 = SQRT(3._JWRB)

      INTEGER(KIND=JWIM) :: IJ

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: X,XK,D,T0,T0_SQ,OM,ALPH,GAM,DELTA,DELTA_1D
      REAL(KIND=JWRB) :: C4_CONST,ZC1,ZC2,ZC3,ZR
      REAL(KIND=JWRB) :: DELTA_2D,C_0,C_S_SQ,V_G,V_G_SQ,ZFAC,ZFAC1,ZFAC2
      REAL(KIND=JWRB) :: XKAPPA1,ALPHA,XJ
      REAL(KIND=JWRB) :: ZEPSILON, ZSQREPSILON
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TRANSF

!-----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STAT_NL',0,ZHOOK_HANDLE)

!     CONSTANTS.

      ZEPSILON=10._JWRB*EPSILON(ZEPSILON)
      ZSQREPSILON=SQRT(ZEPSILON)
 
      C4_CONST  = 0.9_JWRB*PI/(3._JWRB*SQRT3)
      ZC1 = 4._JWRB*SQRT3/PI
      ZC2 = (1._JWRB/3._JWRB+2._JWRB*SQRT3/PI)
      ZC3 = 2._JWRB*SQRT3/PI-4._JWRB/3._JWRB

!     RESULTS FOR A NARROW BAND WAVE TRAIN
  
      DO IJ = KIJS,KIJL
        TRANSF(IJ) = TRANSF_R(XK0(IJ), DPTH(IJ))
      ENDDO

      DO IJ = KIJS,KIJL
        D   = DPTH(IJ)
        IF (XM0(IJ) > ZEPSILON .AND. D > 0._JWRB .AND. XK0(IJ) > 0._JWRB) THEN
          XK  = MAX(XK0(IJ),XKDMIN/D)
          X   = XK*D
          T0  = TANH(X)
          OM  = SQRT(G*XK*T0) 
          T0_SQ = T0**2
          ALPH = XK/(4._JWRB*T0_SQ*T0)*(3._JWRB-T0_SQ)
          GAM = -0.5_JWRB*ALPH**2
          C_0      = OM/XK
          C_S_SQ   = G*D
          IF (X > DKMAX) THEN
            V_G = 0.5_JWRB*C_0
          ELSEIF (X < EPS) THEN
            V_G = C_0
          ELSE
            V_G = 0.5_JWRB*C_0*(1._JWRB+2._JWRB*X/SINH(2._JWRB*X))
          ENDIF
          V_G_SQ=V_G**2

          ZFAC     = -0.25_JWRB*XK*C_S_SQ/(C_S_SQ-V_G_SQ)
          DELTA_1D = ZFAC*(2._JWRB*(1._JWRB-T0_SQ)/T0+1._JWRB/X)

          ZFAC1    = 0.5_JWRB*C_0*C_S_SQ*V_G/T0

          XKAPPA1  = ZFAC1*(2._JWRB*C_0+V_G*(1._JWRB-T0_SQ))/(C_S_SQ-V_G_SQ)
          ALPHA    = (1._JWRB-V_G_SQ/C_S_SQ)*C_0**2/V_G_SQ
          ZFAC2    = SIG_TH(IJ)**2/(SIG_TH(IJ)**2+ALPHA*XNU(IJ)**2)
          DELTA_2D = 0.5_JWRB*XK**2*XKAPPA1/(OM*C_S_SQ)*ZFAC2
  
          DELTA    = DELTA_1D+DELTA_2D
 
!         WAVE-INDUCED MEAN SEA LEVEL
 
          ETA_M(IJ) = 2._JWRB*XM0(IJ)*DELTA
 
!         SKEWNESS
 

          C3(IJ) = CONST_C3*SQRT(XM0(IJ))*(ALPH+C3DELTA_ADJ*DELTA)
          C3(IJ)  = MAX(MIN(C3MAX,C3(IJ)),C3MIN)
 
!         BOUND KURTOSIS
 
          C4_B(IJ) = CONST_C4*XM0(IJ)*(GAM+ALPH**2+(ALPH+DELTA)**2)
 
!         DYNAMIC KURTOSIS
 
!         USE JANSSEN AND JANSSEN (2018) PARAMETRIZATION
 
          R(IJ) = MAX(MIN(TRANSF(IJ)*(SIG_TH(IJ)/XNU(IJ))**2,RMAX),RMIN)
          ZR    = R(IJ)
          IF (ZR > 1._JWRB) THEN
            XJ =-C4_CONST/ZR*(1._JWRB-ZC1/SQRT(ZR)+ZC2/ZR+ZC3/ZR**2)
          ELSE
            XJ = C4_CONST*(1._JWRB-ZC1*SQRT(ZR)+ZC2*ZR+ZC3*ZR**2)
          ENDIF

          C4_DYN(IJ)  = XJ*BF2(IJ)
 
!         DETERMINE TOTAL KURTOSIS AND PROTECT
 
          C4(IJ)  = C4_DYN(IJ)+C4_B(IJ)
          C4(IJ)  = MAX(MIN(C4MAX,C4(IJ)),C4MIN)
        ELSE
          ETA_M(IJ)  = 0._JWRB
          C3(IJ)     = 0._JWRB
          C4(IJ)     = 0._JWRB
          R(IJ)      = 0._JWRB
          C4_DYN(IJ) = 0._JWRB
          C4_B(IJ)   = 0._JWRB
        ENDIF
      ENDDO
 
      IF (LHOOK) CALL DR_HOOK('STAT_NL',1,ZHOOK_HANDLE)

      END SUBROUTINE STAT_NL
