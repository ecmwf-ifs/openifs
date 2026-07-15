! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE CAL_SECOND_ORDER_SPEC(KIJS, KIJL, F1, WAVNUM, DEPTH, SIG)
 
!***  *CAL_SEC_ORDER_SPEC*   DETERMINES SECOND_ORDER SPECTRUM
 
!     PETER JANSSEN JULY 2010
 
!     PURPOSE. 
!     --------
 
!              DETERMINATION OF SECOND-ORDER SPECTRUM
 
!     INTERFACE.
!     ----------
!              *CALL*  *CAL_SEC_ORDER_SPEC(KIJS,KIJL,F1,WAVNUM,DEPTH,SIG)*
 
!                       INPUT:
!                            *KIJS*   - FIRST GRIDPOINT 
!                            *KIJL*   - LAST GRIDPOINT
!                            *F1*     - 2-D FREE WAVE SPECTRUM (at input)
!                            *WAVNUM* - WAVE NUMBER
!                            *DEPTH*  - DEPTH ARRAY
!                            *SIG*    - DIRECTION OF MAPPING:
!                                       FORWARD: SIG = +1
!                                       INVERSE: SIG = -1
 
!                       OUTPUT: 
!                            *F1*    - 2-D SPECTRUM INCLUDING SECOND-ORDER
!                                      CORRECTION (at output)
 
!     METHOD.
!     -------
!              IS DESCRIBED IN JANSSEN (2009), JFM, 637, 1-44.
 
!     EXTERNALS.
!     ----------
!              SECSPOM
 
!-----------------------------------------------------------------------
 
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED , ONLY : FR, DFIM, DELTH,TH,FRATIO
      USE YOWPARAM, ONLY : NANG, NFRE
      USE YOWPCONS, ONLY : G, PI, ZPI
      USE YOWSHAL , ONLY : NDEPTH, DEPTHA, DEPTHD
      USE YOWTABL , ONLY : MR, XMR, MA, XMA, NFREH, NANGH, NMAX,         &
     &                     OMEGA, DFDTH, THH, DELTHH, IM_P, IM_M,        &
     &                     TA, TB, TC_QL, TT_4M, TT_4P
      USE YOWTEST , ONLY : IU06

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!-----------------------------------------------------------------------

      IMPLICIT NONE
#include "fkmean.intfb.h"
#include "secspom.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: F1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH
      REAL(KIND=JWRB), INTENT(IN) :: SIG

      INTEGER(KIND=JWIM) :: IJ,M,K,K0,M0,MP,KP,MM,KM,KL,KLL,ML

      REAL(KIND=JWRB) :: FRAC,CO1,DEL,DELF,D1,D2,D3,D4,C1
      REAL(KIND=JWRB) :: C2,XM,XK,OMSTART,AREA,SUM,SUM1,SUM3,GAM_B_J,ZFAC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMEAN, FMEAN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: F1MEAN, AKMEAN, XKMEAN, EMAXL  
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: F3
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANGH,NFREH) :: PF1, PF3

!-----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CAL_SECOND_ORDER_SPEC',0,ZHOOK_HANDLE)
 
!***  1. DETERMINE SECOND ORDER CORRECTION TO THE SPECTRUM
!     ----------------------------------------------------
 
      GAM_B_J = 0.6_JWRB
      ZFAC = GAM_B_J**2/16.0_JWRB
      FRAC = FRATIO-1.0_JWRB
      OMSTART = ZPI*FR(1)

      CALL FKMEAN(KIJS, KIJL, F1, WAVNUM,                     &
     &            EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN)
 
!***  1.1 INTERPOLATION OR NOT.
!     ------------------------
 
      IF (MR == 1 .AND. MA == 1) THEN
 
!***     1.11 NO INTERPOLATION.
!        ----------------------
  
         CALL SECSPOM(F1,F3,KIJS,KIJL,NFRE,NANG,NMAX,NDEPTH,DEPTHA,       &
     &                DEPTHD,OMSTART,FRAC,MR,DFDTH,OMEGA,DEPTH,         &
     &                AKMEAN,TA,TB,TC_QL,TT_4M,TT_4P,IM_P,IM_M)
         DO M=1,NFRE
            DO K=1,NANG
               DO IJ=KIJS,KIJL
                  DELF = F3(IJ,K,M)
                  F1(IJ,K,M)=MAX(MIN(0.000001_JWRB,F1(IJ,K,M)),F1(IJ,K,M)+SIG*DELF)
               ENDDO
            ENDDO
         ENDDO

      ELSE

 
!***     1.12 ENERGY CONSERVING INTERPOLATION SCHEME
!        -------------------------------------------
 
         PF1(:,:,:) = 0._JWRB
         DO M=1,NFREH
            M0 = MR*M
            DO K=1,NANGH
               K0 = MA*K+1
               IF (K0 > NANG) K0 = K0-NANG
               IF (K0 < 1) K0 = K0+NANG
               DO IJ=KIJS,KIJL
                  PF1(IJ,K,M)=PF1(IJ,K,M)+F1(IJ,K0,M0)
               ENDDO
            ENDDO
         ENDDO

!***     1.13 DETERMINE SECOND-ORDER SPEC
!        --------------------------------

         CALL SECSPOM(PF1,PF3,KIJS,KIJL,NFREH,NANGH,NMAX,NDEPTH,DEPTHA,   &
     &                DEPTHD,OMSTART,FRAC,MR,DFDTH,OMEGA,DEPTH,           &
     &                AKMEAN,TA,TB,TC_QL,TT_4M,TT_4P,IM_P,IM_M)
 
!***     2.24 INTERPOLATE TOWARDS HIGH-RES GRID
!        --------------------------------------
 
         DO IJ=KIJS,KIJL
           IF (EMEAN(IJ) <= ZFAC*DEPTH(IJ)**2) THEN 
             EMAXL(IJ)=1._JWRB
           ELSE
             EMAXL(IJ)=0._JWRB
           ENDIF
         ENDDO

         DO M=1,NFRE
            XM = REAL(M/MR)
!!!            D1 = REAL(M)/REAL(MR)-XM
            M0 = INT(XM)
            IF(M0 < 1) THEN
              M0 = 1
              MP = M0+1
              D1 = 1._JWRB
            ELSE IF(M0 < NFREH) THEN
              MP = M0+1
              D1 = (FR(M)-FR(MR*M0))/(FR(MR*MP)-FR(MR*M0)) 
            ELSE
              M0 = NFREH
              MP = M0
              D1 = 0._JWRB
            ENDIF
            D2 = 1.0_JWRB-D1

            DO K=1,NANG
               XK = REAL((K-1)/MA)
               K0 = INT(XK)
               D3 = REAL(K-1)/REAL(MA)-XK
               D4 = 1.-D3
               IF (K0 < 1) K0 = K0+NANGH
               KP = K0+1
               IF (KP > NANGH) KP = KP-NANGH

               DO IJ=KIJS,KIJL
                  C1 = PF3(IJ,K0,M0)*D4+PF3(IJ,KP,M0)*D3
                  C2 = PF3(IJ,KP,MP)*D3+PF3(IJ,K0,MP)*D4
                  DELF = C1*D2+C2*D1

                  F1(IJ,K,M)=MAX(MIN(0.000001_JWRB,F1(IJ,K,M)),F1(IJ,K,M)+EMAXL(IJ)*SIG*DELF)

               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF (LHOOK) CALL DR_HOOK('CAL_SECOND_ORDER_SPEC',1,ZHOOK_HANDLE)

END SUBROUTINE CAL_SECOND_ORDER_SPEC
