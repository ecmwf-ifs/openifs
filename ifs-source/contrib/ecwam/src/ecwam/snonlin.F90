! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SNONLIN (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)

! ----------------------------------------------------------------------

!**** *SNONLIN* - COMPUTATION OF NONLINEAR TRANSFER RATE AND ITS
!****             FUNCTIONAL DERIVATIVE (DIAGONAL TERMS ONLY) AND
!****             ADDITION TO CORRESPONDING NET EXPRESSIONS.

!     S.D. HASSELMANN.  MPI

!     G. KOMEN, P. JANSSEN   KNMI             MODIFIED TO SHALLOW WATER
!     H. GUENTHER, L. ZAMBRESKY               OPTIMIZED
!     H. GUENTHER       GKSS/ECMWF  JUNE 1991 INTERACTIONS BETWEEN DIAG-
!                                             AND PROGNOSTIC PART.
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
!     P. JANSSEN  ECMWF  JUNE 2005       IMPROVED SCALING IN SHALLOW
!                                        WATER
!     J. BIDLOT   ECMWF  AUGUST 2006     KEEP THE OLD FORMULATION
!                                        UNDER A SWITCH (ISNONLIN = 0 for OLD
!                                                                 = 1 for NEW
!                                        BE AWARE THAT THE OLD FORMULATION
!                                        REQUIRES THE MEAN WAVE NUMBER AKMEAN.
!     J. BIDLOT   ECMWF  JANUARY 2012    ADD EXTENSION TO LOW FREQUENCIES 
!                                        OPTIMISATION FOR IBM.

!*    PURPOSE.
!     --------

!       SEE ABOVE.

!**   INTERFACE.
!     ----------

!       *CALL* *SNONLIN (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)*
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FL1*    - SPECTRUM.
!          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*     - TOTAL SOURCE FUNCTION ARRAY.
!          *WAVNUM* - WAVE NUMBER.
!          *DEPTH*  - WATER DEPTH.
!          *AKMEAN* - MEAN WAVE NUMBER  BASED ON sqrt(1/k)*F INTGRATION

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

      USE YOWFRED  , ONLY : FR       ,FRATIO   ,ZPIFR
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWINDN  , ONLY : IKP      ,IKP1     ,IKM      ,IKM1     ,    &
     &            K1W      ,K2W      ,K11W     ,K21W     ,AF11     ,    &
     &            FKLAP    ,FKLAP1   ,FKLAM    ,FKLAM1   ,              &
     &            DAL1     ,DAL2     ,                                  &
     &            MFRSTLW  ,MLSTHG   ,                                  &
     &            KFRH     ,INLCOEF  ,RNLCOEF
      USE YOWSTAT  , ONLY : ISNONLIN
      USE YOWPCONS , ONLY : GM1

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "peak_ang.intfb.h"
#include "transf.intfb.h"
#include "transf_snl.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT):: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH, AKMEAN


      INTEGER(KIND=JWIM) :: IJ, K, M, MC, KH, K1, K2, K11, K21
      INTEGER(KIND=JWIM) :: MP, MP1, MM, MM1, IC, IP, IP1, IM, IM1
      INTEGER(KIND=JWIM) :: MFR1STFR, MFRLSTFR

      REAL(KIND=JWRB), PARAMETER :: ENH_MAX=10.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: ENH_MIN=0.1_JWRB   ! to prevent ENH to become too small
      REAL(KIND=JWRB) :: XK 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,MLSTHG) :: ENH

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: FTAIL, FKLAMP, GW1, GW2, GW3, GW4
      REAL(KIND=JWRB) :: FKLAMPA, FKLAMPB, FKLAMP2, FKLAMP1, FKLAPA2
      REAL(KIND=JWRB) :: FKLAPB2, FKLAP12, FKLAP22, FKLAMM, FKLAMM1
      REAL(KIND=JWRB) :: GW5, GW6, GW7, GW8, FKLAMMA, FKLAMMB, FKLAMM2
      REAL(KIND=JWRB) :: FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
      REAL(KIND=JWRB) :: SAP, SAM, FIJ, FAD1, FAD2, FCEN

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: XNU, SIG_TH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: FTEMP, AD, DELAD, DELAP, DELAM, ENHFR


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SNONLIN',0,ZHOOK_HANDLE)


!*    1. SHALLOW WATER SCALING 
!        ---------------------

      SELECT CASE (ISNONLIN)
      CASE(0)      
        DO IJ=KIJS,KIJL
          ENHFR(IJ) = MAX(0.75_JWRB*DEPTH(IJ)*AKMEAN(IJ), 0.5_JWRB)
          ENHFR(IJ) = 1.0_JWRB+(5.5_JWRB/ENHFR(IJ)) * (1.0_JWRB-.833_JWRB*ENHFR(IJ)) * EXP(-1.25_JWRB*ENHFR(IJ))
        ENDDO
        DO MC=1,MLSTHG
          DO IJ=KIJS,KIJL
            ENH(IJ,MC) = ENHFR(IJ)
          ENDDO
        ENDDO

      CASE(1)      
        DO MC=1,NFRE
          DO IJ = KIJS, KIJL 
            ENH(IJ,MC) = MAX(MIN(ENH_MAX,TRANSF(WAVNUM(IJ,MC),DEPTH(IJ))),ENH_MIN)
          ENDDO
        ENDDO
        DO MC=NFRE+1,MLSTHG
          XK = GM1*(ZPIFR(NFRE)*FRATIO**(MC-NFRE))**2
          DO IJ = KIJS, KIJL 
            ENH(IJ,MC) = MAX(MIN(ENH_MAX,TRANSF(XK,DEPTH(IJ))),ENH_MIN)
          ENDDO
        ENDDO

      CASE(2)      
        CALL PEAK_ANG(KIJS, KIJL, FL1, XNU, SIG_TH)
        DO MC=1,NFRE
          DO IJ = KIJS, KIJL 
            ENH(IJ,MC) = TRANSF_SNL(WAVNUM(IJ,MC), DEPTH(IJ), XNU(IJ), SIG_TH(IJ))
          ENDDO
        ENDDO
        DO MC=NFRE+1,MLSTHG
          XK = GM1*(ZPIFR(NFRE)*FRATIO**(MC-NFRE))**2
          DO IJ = KIJS, KIJL 
            ENH(IJ,MC) = TRANSF_SNL(XK, DEPTH(IJ), XNU(IJ), SIG_TH(IJ))
          ENDDO
        ENDDO
      END SELECT


!*    2. FREQUENCY LOOP.
!        ---------------

      MFR1STFR=-MFRSTLW+1
      MFRLSTFR=NFRE-KFRH+MFR1STFR

      DO MC=1,MLSTHG
        MP  = IKP (MC)
        MP1 = IKP1(MC)
        MM  = IKM (MC)
        MM1 = IKM1(MC)
        IC  = INLCOEF(1,MC)
        IP  = INLCOEF(2,MC) 
        IP1 = INLCOEF(3,MC)
        IM  = INLCOEF(4,MC) 
        IM1 = INLCOEF(5,MC)

        FTAIL  = RNLCOEF(1,MC)

        FKLAMP  = FKLAP(MC)
        FKLAMP1 = FKLAP1(MC)
        GW1 = RNLCOEF(2,MC) 
        GW2 = RNLCOEF(3,MC) 
        GW3 = RNLCOEF(4,MC)
        GW4 = RNLCOEF(5,MC) 
        FKLAMPA = RNLCOEF(6,MC)
        FKLAMPB = RNLCOEF(7,MC)
        FKLAMP2 = RNLCOEF(8,MC) 
        FKLAMP1 = RNLCOEF(9,MC) 
        FKLAPA2 = RNLCOEF(10,MC) 
        FKLAPB2 = RNLCOEF(11,MC) 
        FKLAP12 = RNLCOEF(12,MC)
        FKLAP22 = RNLCOEF(13,MC) 

        FKLAMM  = FKLAM(MC)
        FKLAMM1 = FKLAM1(MC)
        GW5 = RNLCOEF(14,MC)
        GW6 = RNLCOEF(15,MC)
        GW7 = RNLCOEF(16,MC) 
        GW8 = RNLCOEF(17,MC)
        FKLAMMA = RNLCOEF(18,MC) 
        FKLAMMB = RNLCOEF(19,MC) 
        FKLAMM2 = RNLCOEF(20,MC) 
        FKLAMM1 = RNLCOEF(21,MC) 
        FKLAMA2 = RNLCOEF(22,MC)
        FKLAMB2 = RNLCOEF(23,MC) 
        FKLAM12 = RNLCOEF(24,MC) 
        FKLAM22 = RNLCOEF(25,MC) 

        DO IJ=KIJS,KIJL
          FTEMP(IJ) = AF11(MC)*ENH(IJ,MC)
        ENDDO


        IF (MC > MFR1STFR .AND. MC < MFRLSTFR ) THEN
!       the interactions for MC are all within the fully resolved spectral domain

          DO KH=1,2
            DO K=1,NANG
              K1  = K1W (K,KH)
              K2  = K2W (K,KH)
              K11 = K11W(K,KH)
              K21 = K21W(K,KH)

!*    2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND
!*            DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!             ----------------------------------------------
              DO IJ=KIJS,KIJL
                SAP = GW1*FL1(IJ,K1 ,IP ) + GW2*FL1(IJ,K11,IP )             &
     &              + GW3*FL1(IJ,K1 ,IP1) + GW4*FL1(IJ,K11,IP1)
                SAM = GW5*FL1(IJ,K2 ,IM ) + GW6*FL1(IJ,K21,IM )             &
     &              + GW7*FL1(IJ,K2 ,IM1) + GW8*FL1(IJ,K21,IM1)
!!!! not needed ftail always=1.                FIJ = FL1(IJ,K  ,IC )*FTAIL
                FIJ = FL1(IJ,K  ,IC )
                FAD1 = FIJ*(SAP+SAM)
                FAD2 = FAD1-2.0_JWRB*SAP*SAM
                FAD1 = FAD1+FAD2
                FCEN = FTEMP(IJ)*FIJ
                AD(IJ) = FAD2*FCEN
                DELAD(IJ) = FAD1*FTEMP(IJ)
                DELAP(IJ) = (FIJ-2.0_JWRB*SAM)*DAL1*FCEN
                DELAM(IJ) = (FIJ-2.0_JWRB*SAP)*DAL2*FCEN
              ENDDO

              DO IJ=KIJS,KIJL
                SL(IJ,K  ,MC ) = SL(IJ,K  ,MC ) - 2.0_JWRB*AD(IJ)
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K  ,MC ) = FLD(IJ,K  ,MC ) - 2.0_JWRB*DELAD(IJ)
              ENDDO
              DO IJ=KIJS,KIJL
                SL(IJ,K2 ,MM ) = SL(IJ,K2 ,MM ) + AD(IJ)*FKLAMM1
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K2 ,MM ) = FLD(IJ,K2 ,MM ) + DELAM(IJ)*FKLAM12
              ENDDO
              DO IJ=KIJS,KIJL
                SL(IJ,K21,MM ) = SL(IJ,K21,MM ) + AD(IJ)*FKLAMM2
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K21,MM ) = FLD(IJ,K21,MM ) + DELAM(IJ)*FKLAM22
              ENDDO
              DO IJ=KIJS,KIJL
                SL(IJ,K2 ,MM1) = SL(IJ,K2 ,MM1) + AD(IJ)*FKLAMMA
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K2 ,MM1) = FLD(IJ,K2 ,MM1) + DELAM(IJ)*FKLAMA2
              ENDDO
              DO IJ=KIJS,KIJL
                SL(IJ,K21,MM1) = SL(IJ,K21,MM1) + AD(IJ)*FKLAMMB
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K21,MM1) = FLD(IJ,K21,MM1) + DELAM(IJ)*FKLAMB2
              ENDDO
              DO IJ=KIJS,KIJL
                SL(IJ,K1 ,MP ) = SL(IJ,K1 ,MP ) + AD(IJ)*FKLAMP1
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K1 ,MP ) = FLD(IJ,K1 ,MP ) + DELAP(IJ)*FKLAP12
              ENDDO
              DO IJ=KIJS,KIJL
                SL(IJ,K11,MP ) = SL(IJ,K11,MP ) + AD(IJ)*FKLAMP2
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K11,MP ) = FLD(IJ,K11,MP ) + DELAP(IJ)*FKLAP22
              ENDDO
              DO IJ=KIJS,KIJL
                SL(IJ,K1 ,MP1) = SL(IJ,K1 ,MP1) + AD(IJ)*FKLAMPA
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K1 ,MP1) = FLD(IJ,K1 ,MP1) + DELAP(IJ)*FKLAPA2
              ENDDO
              DO IJ=KIJS,KIJL
                SL(IJ,K11,MP1) = SL(IJ,K11,MP1) + AD(IJ)*FKLAMPB
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K11,MP1) = FLD(IJ,K11,MP1) + DELAP(IJ)*FKLAPB2
              ENDDO
            ENDDO
          ENDDO

        ELSEIF (MC >= MFRLSTFR ) THEN
          DO KH=1,2
            DO K=1,NANG
              K1  = K1W (K,KH)
              K2  = K2W (K,KH)
              K11 = K11W(K,KH)
              K21 = K21W(K,KH)

              DO IJ=KIJS,KIJL
                SAP = GW1*FL1(IJ,K1 ,IP ) + GW2*FL1(IJ,K11,IP )         &
     &              + GW3*FL1(IJ,K1 ,IP1) + GW4*FL1(IJ,K11,IP1)
                SAM = GW5*FL1(IJ,K2 ,IM ) + GW6*FL1(IJ,K21,IM )         &
     &              + GW7*FL1(IJ,K2 ,IM1) + GW8*FL1(IJ,K21,IM1)
                FIJ = FL1(IJ,K  ,IC )*FTAIL
                FAD1 = FIJ*(SAP+SAM)
                FAD2 = FAD1-2.0_JWRB*SAP*SAM
                FAD1 = FAD1+FAD2
                FCEN = FTEMP(IJ)*FIJ
                AD(IJ) = FAD2*FCEN
                DELAD(IJ) = FAD1*FTEMP(IJ)
                DELAP(IJ) = (FIJ-2.0_JWRB*SAM)*DAL1*FCEN
                DELAM(IJ) = (FIJ-2.0_JWRB*SAP)*DAL2*FCEN
              ENDDO

              DO IJ=KIJS,KIJL
                SL(IJ,K2 ,MM ) = SL(IJ,K2 ,MM ) + AD(IJ)*FKLAMM1
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K2 ,MM ) = FLD(IJ,K2 ,MM ) + DELAM(IJ)*FKLAM12
              ENDDO
              DO IJ=KIJS,KIJL
                SL(IJ,K21,MM ) = SL(IJ,K21,MM ) + AD(IJ)*FKLAMM2
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K21,MM ) = FLD(IJ,K21,MM ) + DELAM(IJ)*FKLAM22
              ENDDO

              IF (MM1 <= NFRE) THEN
                DO IJ=KIJS,KIJL
                  SL(IJ,K2 ,MM1) = SL(IJ,K2 ,MM1) + AD(IJ)*FKLAMMA
                ENDDO
                DO IJ=KIJS,KIJL
                  FLD(IJ,K2 ,MM1) = FLD(IJ,K2 ,MM1) + DELAM(IJ)*FKLAMA2
                ENDDO
                DO IJ=KIJS,KIJL
                  SL(IJ,K21,MM1) = SL(IJ,K21,MM1) + AD(IJ)*FKLAMMB
                ENDDO
                DO IJ=KIJS,KIJL
                  FLD(IJ,K21,MM1) = FLD(IJ,K21,MM1) + DELAM(IJ)*FKLAMB2
                ENDDO

                IF (MC <= NFRE) THEN
                  DO IJ=KIJS,KIJL
                    SL(IJ,K  ,MC ) = SL(IJ,K  ,MC ) - 2.0_JWRB*AD(IJ)
                  ENDDO
                  DO IJ=KIJS,KIJL
                    FLD(IJ,K  ,MC ) = FLD(IJ,K  ,MC ) - 2.0_JWRB*DELAD(IJ)
                  ENDDO

                  IF (MP <= NFRE) THEN
                    DO IJ=KIJS,KIJL
                      SL(IJ,K1 ,MP ) = SL(IJ,K1 ,MP ) + AD(IJ)*FKLAMP1
                    ENDDO
                    DO IJ=KIJS,KIJL
                      FLD(IJ,K1 ,MP ) = FLD(IJ,K1 ,MP )                 &
     &                               + DELAP(IJ)*FKLAP12
                    ENDDO
                    DO IJ=KIJS,KIJL
                      SL(IJ,K11,MP ) = SL(IJ,K11,MP ) + AD(IJ)*FKLAMP2
                    ENDDO
                    DO IJ=KIJS,KIJL
                      FLD(IJ,K11,MP ) = FLD(IJ,K11,MP )                 &
     &                               + DELAP(IJ)*FKLAP22
                    ENDDO

                    IF (MP1 <= NFRE) THEN
                      DO IJ=KIJS,KIJL
                        SL(IJ,K1 ,MP1) = SL(IJ,K1 ,MP1)                 &
     &                                 + AD(IJ)*FKLAMPA
                      ENDDO
                      DO IJ=KIJS,KIJL
                        FLD(IJ,K1 ,MP1) = FLD(IJ,K1 ,MP1)               &
     &                                 + DELAP(IJ)*FKLAPA2
                      ENDDO
                      DO IJ=KIJS,KIJL
                        SL(IJ,K11,MP1) = SL(IJ,K11,MP1)                 &
     &                                 + AD(IJ)*FKLAMPB
                      ENDDO
                      DO IJ=KIJS,KIJL
                        FLD(IJ,K11,MP1) = FLD(IJ,K11,MP1)               &
     &                                 + DELAP(IJ)*FKLAPB2
                      ENDDO
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDDO

        ELSE

          DO KH=1,2
            DO K=1,NANG
              K1  = K1W (K,KH)
              K2  = K2W (K,KH)
              K11 = K11W(K,KH)
              K21 = K21W(K,KH)

              DO IJ=KIJS,KIJL
                SAP = GW1*FL1(IJ,K1 ,IP ) + GW2*FL1(IJ,K11,IP )        &
     &              + GW3*FL1(IJ,K1 ,IP1) + GW4*FL1(IJ,K11,IP1)
                SAM = GW5*FL1(IJ,K2 ,IM ) + GW6*FL1(IJ,K21,IM )        &
     &              + GW7*FL1(IJ,K2 ,IM1) + GW8*FL1(IJ,K21,IM1)
                FIJ = FL1(IJ,K  ,IC )*FTAIL
                FAD1 = FIJ*(SAP+SAM)
                FAD2 = FAD1-2.0_JWRB*SAP*SAM
                FAD1 = FAD1+FAD2
                FCEN = FTEMP(IJ)*FIJ
                AD(IJ) = FAD2*FCEN
                DELAD(IJ) = FAD1*FTEMP(IJ)
                DELAP(IJ) = (FIJ-2.0_JWRB*SAM)*DAL1*FCEN
                DELAM(IJ) = (FIJ-2.0_JWRB*SAP)*DAL2*FCEN
              ENDDO

              IF (MM1 >= 1) THEN
                DO IJ=KIJS,KIJL
                  SL(IJ,K2 ,MM1) = SL(IJ,K2 ,MM1) + AD(IJ)*FKLAMMA
                ENDDO
                DO IJ=KIJS,KIJL
                  FLD(IJ,K2 ,MM1) = FLD(IJ,K2 ,MM1) + DELAM(IJ)*FKLAMA2
                ENDDO
                DO IJ=KIJS,KIJL
                  SL(IJ,K21,MM1) = SL(IJ,K21,MM1) + AD(IJ)*FKLAMMB
                ENDDO
                DO IJ=KIJS,KIJL
                  FLD(IJ,K21,MM1) = FLD(IJ,K21,MM1) + DELAM(IJ)*FKLAMB2
                ENDDO
              ENDIF

              DO IJ=KIJS,KIJL
                SL(IJ,K  ,MC ) = SL(IJ,K  ,MC ) - 2.0_JWRB*AD(IJ)
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K  ,MC ) = FLD(IJ,K  ,MC ) - 2.0_JWRB*DELAD(IJ)
              ENDDO
              DO IJ=KIJS,KIJL
                SL(IJ,K1 ,MP ) = SL(IJ,K1 ,MP ) + AD(IJ)*FKLAMP1
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K1 ,MP ) = FLD(IJ,K1 ,MP ) + DELAP(IJ)*FKLAP12
              ENDDO
              DO IJ=KIJS,KIJL
                SL(IJ,K11,MP ) = SL(IJ,K11,MP ) + AD(IJ)*FKLAMP2
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K11,MP ) = FLD(IJ,K11,MP ) + DELAP(IJ)*FKLAP22
              ENDDO
              DO IJ=KIJS,KIJL
                SL(IJ,K1 ,MP1) = SL(IJ,K1 ,MP1) + AD(IJ)*FKLAMPA
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K1 ,MP1) = FLD(IJ,K1 ,MP1) + DELAP(IJ)*FKLAPA2
              ENDDO
              DO IJ=KIJS,KIJL
                SL(IJ,K11,MP1) = SL(IJ,K11,MP1) + AD(IJ)*FKLAMPB
              ENDDO
              DO IJ=KIJS,KIJL
                FLD(IJ,K11,MP1) = FLD(IJ,K11,MP1) + DELAP(IJ)*FKLAPB2
              ENDDO
            ENDDO
          ENDDO

        ENDIF

!*    BRANCH BACK TO 2. FOR NEXT FREQUENCY.

      ENDDO

      IF (LHOOK) CALL DR_HOOK('SNONLIN',1,ZHOOK_HANDLE)

      END SUBROUTINE SNONLIN
