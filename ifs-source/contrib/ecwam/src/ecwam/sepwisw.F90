! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SEPWISW (KIJS, KIJL, MIJ, FL1, XLLWS, CINV,      &
     &                    UFRIC    ,WSWAVE   ,WDWAVE, COSWDIF,    &
     &                    ESWELL   ,FSWELL   ,THSWELL  ,          &
     &                    P1SWELL  ,P2SWELL  ,SPRDSWELL,          &
     &                    ESEA     ,FSEA     ,THWISEA  ,          &
     &                    P1SEA    ,P2SEA    ,SPRDSEA  ,          &
     &                    EMTRAIN  ,THTRAIN  ,PMTRAIN)

! ----------------------------------------------------------------------

!**** *SEPWISW* - COMPUTES THE SWELL ENERGY, THE MEAN SWELL DIRECTION,
!****             THE MEAN SWELL FREQUENCY AND THE MEAN WINDSEA DIR.

!     P.LIONELLO     FEBRUARY 87

!     L.ZAMBRESKY    NOVEMBER 87   GKSS/ECMWF   OPTIMIZED SUB.
!     J. BIDLOT      FEBRARY 1996       ECMWF   MESSAGE PASSING
!     J. BIDLOT      MAY     1999       ECMWF   ADD HIGH FREQUENCY TAIL
!     J. BIDLOT      APRIL   2000       ECMWF   ADD  EXTRA PARAMETERS
!     J. BIDLOT      DECEMBER2003       ECMWF   MOVE ALL ALLOCATION TO
!                                               *OUTBS*
!                                               TO MAKE IT CALLABLE IN AN
!                                               OPENMP LOOP.
!     J. BIDLOT     JUNE 2013           ECMWF    ADD PARTITIONING SCHEME

!*    PURPOSE.
!     --------

!       TO SEPARATE THE SWELL FROM THE WIND INTERACTING SEA

!**   INTERFACE.
!     ----------

!       *CALL* SEPWISW (KIJS, KIJL, MIJ, FL1, XLLWS,
!    &                  UFRIC, WSWAVE, WDWAVE, COSWDIF,
!    &                  ESWELL   ,FSWELL   ,THSWELL  ,
!    &                  P1SWELL  ,P2SWELL  ,SPRDSWELL,
!    &                  ESEA     ,FSEA     ,THWISEA  ,
!    &                  P1SEA    ,P2SEA    ,SPRDSEA  ,
!    &                  EMTRAIN  ,THTRAIN  ,PMTRAIN)
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *FL1*    - BLOCK OF SPECTRA
!          *XLLWS*  - WINDSEA MASK FROM INPUT SOURCE TERM
!          *UFRIC*  - NEW FRICTION VELOCITY IN M/S.
!          *WSWAVE* - LATESt WIND SPEED.
!          *WDWAVE* - LATEST WIND DIRECTION.
!          *COSWDIF*- COSINE (WDWAVE - WAVES DIRECTIONS)
!          *        - SWELL and WINDSEA PARAMETERS.

!     METHOD.
!     -------

!       THE WAVES WHICH DO NOT INTERACT WITH THE WIND ARE
!       CONSIDERED SWELL.
!       IF LLPARTITION THEN A SPECTRAl PARTITIONING SCHEME WILL BE USED
!       TO SPLIT THE SWELL SPECTRUM INTO ITS TWO MAIN COMPONENTS.

!     EXTERNALS.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : NTRAIN   ,LLPARTITION
      USE YOWFRED  , ONLY : FR       ,TH       ,FRIC     ,OLDWSFC, ZPIFR
      USE YOWPCONS , ONLY : G        ,EPSMIN
      USE YOWMAP   , ONLY : CLDOMAIN
      USE YOWPARAM , ONLY : NANG     ,NFRE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! -----------------------------------------------------------------------

      IMPLICIT NONE

#include "femean.intfb.h"
#include "mwp1.intfb.h"
#include "mwp2.intfb.h"
#include "sep3tr.intfb.h"
#include "sthq.intfb.h"
#include "wdirspread.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: MIJ
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: XLLWS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: CINV 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UFRIC, WSWAVE, WDWAVE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG), INTENT(IN) :: COSWDIF
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: ESWELL ,FSWELL ,THSWELL, P1SWELL, P2SWELL, SPRDSWELL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: ESEA   ,FSEA   ,THWISEA, P1SEA  , P2SEA  , SPRDSEA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NTRAIN), INTENT(OUT) :: EMTRAIN, THTRAIN, PMTRAIN


      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: COEF
      REAL(KIND=JWRB) :: CHECKTA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: R
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: XINVWVAGE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: DIRCOEF
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: SWM, F1

      LOGICAL :: LLPEAKF

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SEPWISH',0,ZHOOK_HANDLE)

!*    1. THE SWELL DISTRIBUTION IS COMPUTED.
!        -----------------------------------

      COEF = OLDWSFC*FRIC

      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          XINVWVAGE(IJ,M)=UFRIC(IJ)*CINV(IJ,M)
        ENDDO
      ENDDO

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          DIRCOEF(IJ,K)=COEF*COSWDIF(IJ,K)
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            IF (XLLWS(IJ,K,M) /= 0.0_JWRB) THEN
              ! this is windsea 
              SWM(IJ,K,M)=0.0_JWRB
            ELSE
              CHECKTA=XINVWVAGE(IJ,M)*DIRCOEF(IJ,K)
              IF (CHECKTA >= 1.0_JWRB) THEN
                ! this is extra windsea 
                SWM(IJ,K,M)=0.0_JWRB
              ELSE
                ! this is swell
                SWM(IJ,K,M)=1.0_JWRB
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      IF (.NOT.(CLDOMAIN == 's')) THEN
!     CHECK THAT TOTAL SWELL MEAN FREQUENCY IS LOWER THAN WINDSEA ONE
!     OTHERWISE RESET WIND SECTOR TO WINDSEA

!     SWELL:
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            F1(IJ,K,M)=FL1(IJ,K,M)*SWM(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO
      CALL FEMEAN(KIJS, KIJL, F1, ESWELL, FSWELL)

!     WINDSEA:
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            F1(IJ,K,M)=MAX(FL1(IJ,K,M)-F1(IJ,K,M),0.0_JWRB)
          ENDDO
        ENDDO
      ENDDO
      CALL FEMEAN(KIJS, KIJL, F1, ESEA, FSEA)

!     CHECK:
      DO IJ=KIJS,KIJL
         IF (FSWELL(IJ) > 0.96_JWRB*FSEA(IJ)) THEN
           R(IJ)=1.0_JWRB
         ELSE
           R(IJ)=0.0_JWRB
         ENDIF
      ENDDO
      DO K=1,NANG
        DO IJ=KIJS,KIJL
          ! add factor to extend windsea area
          DIRCOEF(IJ,K)=R(IJ)*COEF*SIGN(1.0_JWRB,0.4_JWRB+COSWDIF(IJ,K))
        ENDDO
      ENDDO
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            CHECKTA = XINVWVAGE(IJ,M)*DIRCOEF(IJ,K)
            IF (CHECKTA >= 1.0_JWRB) THEN
              ! this is additional windsea 
              SWM(IJ,K,M)=0.0_JWRB
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!     ADJUST THE LOW FREQUENCY BOUNDARY OF THE WINDSEAS AREA 
!     TO TOPOLOGICALLY CONNECT IT TO THE WINDSEA
      DO IJ=KIJS,KIJL
        DO K=1,NANG
          DO M=NFRE,2,-1
            IF (SWM(IJ,K,M) == 1.0_JWRB .AND. SWM(IJ,K,M-1) == 1.0_JWRB) THEN
              EXIT
            ELSEIF (SWM(IJ,K,M) == 0.0_JWRB .AND. SWM(IJ,K,M-1) == 1.0_JWRB) THEN
               IF (FL1(IJ,K,M) >= FL1(IJ,K,M-1)) SWM(IJ,K,M-1)=0.0_JWRB
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ENDIF ! CLDOMAIN

!     SWELL SPECTRUM
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            F1(IJ,K,M)=MAX(FL1(IJ,K,M),EPSMIN)*SWM(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

!*    2.1 COMPUTATION OF THE PARTITIONED SWELL OUTPUT PARAMETERS
!         ------------------------------------------------------

      IF (LLPARTITION) THEN
        CALL FEMEAN(KIJS, KIJL, F1, ESWELL, FSWELL)
        CALL STHQ(KIJS, KIJL, F1, THSWELL)
        CALL SEP3TR (KIJS, KIJL, FL1, MIJ, WSWAVE, WDWAVE, COSWDIF, &
     &               ESWELL, FSWELL, THSWELL, FSEA,                 &
     &               F1, SWM,                                       &
     &               EMTRAIN  ,THTRAIN  ,PMTRAIN)
      ELSE
        EMTRAIN(:,:) = 0.0_JWRB
        THTRAIN(:,:) = 0.0_JWRB
        PMTRAIN(:,:) = 0.0_JWRB
      ENDIF

!*    2.2 COMPUTATION OF TOTAL SWELL OUTPUT PARAMETERS
!         --------------------------------------------

      CALL FEMEAN(KIJS, KIJL, F1, ESWELL, FSWELL)
      
      CALL STHQ(KIJS, KIJL, F1, THSWELL)

      CALL MWP1(KIJS, KIJL, F1, P1SWELL)

      CALL MWP2(KIJS, KIJL, F1, P2SWELL)

      LLPEAKF = .TRUE.
      CALL WDIRSPREAD (KIJS, KIJL, F1, ESWELL, LLPEAKF, SPRDSWELL)


!*    3. COMPUTATION OF WIND SEA OUTPUT PARAMETERS
!        -----------------------------------------

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            F1(IJ,K,M)=MAX(FL1(IJ,K,M)-F1(IJ,K,M)+EPSMIN,0.0_JWRB)
          ENDDO
        ENDDO
      ENDDO

      CALL FEMEAN(KIJS, KIJL, F1, ESEA, FSEA)

      CALL STHQ(KIJS, KIJL, F1, THWISEA)
!     if there isn't any windsea energy, set the windsea mean direction
!     to the wind direction.
      DO IJ=KIJS,KIJL
        IF (ESEA(IJ) <= 1.0E-9_JWRB) THEN
          THWISEA(IJ)=WDWAVE(IJ)
        ENDIF
      ENDDO

      CALL MWP1(KIJS, KIJL, F1, P1SEA)

      CALL MWP2(KIJS, KIJL, F1, P2SEA)

      LLPEAKF = .TRUE.
      CALL WDIRSPREAD(KIJS, KIJL, F1, ESEA, LLPEAKF, SPRDSEA)

      IF (LHOOK) CALL DR_HOOK('SEPWISH',1,ZHOOK_HANDLE)

      END SUBROUTINE SEPWISW
