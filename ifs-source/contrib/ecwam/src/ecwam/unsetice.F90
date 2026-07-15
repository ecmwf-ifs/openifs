! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE UNSETICE(KIJS, KIJL, DEPTH, EMAXDPT, WDWAVE, WSWAVE, CICOVER, FL1)
! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------
!     TRY TO CREATE A JONSWAP SPECTRUM AT POINTS WHICH WERE ICE IN THE
!     PREVIOUS RUN AND WHICH ARE NOW OPEN SEA.
!     ALSO REINSTATE THE MIMIMUM ENERGY LEVEL.

!**   INTERFACE.
!     ----------
!     *CALL* *UNSETICE(KIJS, KIJL, WVENVI, FF_NOW, FL1)
!     *KIJS*    - LOCAL INDEX OF FIRST GRIDPOINT
!     *KIJL*    - LOCAL  INDEX OF LAST GRIDPOINT
!     *WVENVI*  - WAVE ENVIRONMENT
!     *FF_NOW*  - FORCING FIELDS AT CURRENT TIME
!     *FL1*     - SPECTRA


!     METHOD.
!     -------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FORCING_FIELDS

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRID  , ONLY : DELPHI
      USE YOWJONS  , ONLY : AJONS    ,BJONS    ,DJONS    ,EJONS
      USE YOWICE   , ONLY : FLMIN    ,LICERUN  ,LMASKICE ,CITHRSH
      USE YOWMAP   , ONLY : CLDOMAIN
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : PI       ,G        ,R        ,EPSMIN
      USE YOWSTAT  , ONLY : LBIWBK
      USE YOWTEXT  , ONLY : LRESTARTED

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "sdepthlim.intfb.h"
#include "jonswap.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH, EMAXDPT
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WDWAVE, WSWAVE, CICOVER
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FL1


      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: CILIMIT
      REAL(KIND=JWRB) :: ZGAMMA, SA, SB, FETCH, GX, FPMAX, ZDP
      REAL(KIND=JWRB) :: U10, GXU, UG, FLLOWEST 
      REAL(KIND=JWRB) :: STK(NANG)
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: FPK, ALPHAV
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: ET
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: SPRD, COS2NOISE

      LOGICAL, DIMENSION(KIJS:KIJL) :: LICE2SEA

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('UNSETICE',0,ZHOOK_HANDLE)

      IF (.NOT.LRESTARTED .AND. CLDOMAIN /= 's') THEN

        ZDP=2.0_JWRB/PI
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            COS2NOISE(IJ,K)=MAX(0.0_JWRB,COS(TH(K)-WDWAVE(IJ)))**2
            SPRD(IJ,K)=ZDP*COS2NOISE(IJ,K)
          ENDDO
        ENDDO

        IF (LICERUN .AND. LMASKICE) THEN
          CILIMIT=CITHRSH
        ELSE
!         WHEN WE ARE NOT IMPOSING A FIXED ICE BOUNDARY THEN
!         ONLY POINTS WITH NO ICE WILL BE RESET
          CILIMIT=0.0_JWRB
        ENDIF

!       TRY TO CREATE A JONSWAP SPECTRUM AT POINTS WHICH WERE ICE IN THE
!       PREVIOUS RUN AND WHICH ARE NOW OPEN SEA.
!       THEY ARE CHARCTERISED BY HAVING NO WAVE ENERGY BESIDE NOISE
!       AND HAVE AN CICOVER > CITHRSH .
        DO IJ=KIJS,KIJL
          LICE2SEA(IJ)=.TRUE.
        ENDDO
        DO K=1,NANG
          DO M=1,NFRE
            DO IJ=KIJS,KIJL
               IF (FL1(IJ,K,M) > EPSMIN) LICE2SEA(IJ) = .FALSE. 
            ENDDO
          ENDDO
        ENDDO

        ZGAMMA=3.0_JWRB
        SA=7.0E-02_JWRB
        SB=9.0E-02_JWRB
        FPMAX=FR(NFRE-1)

        DO IJ=KIJS,KIJL
          IF (LICE2SEA(IJ) .AND. CICOVER(IJ) <= CILIMIT) THEN
            IF (DEPTH(IJ) <= 10.0_JWRB) THEN
              FETCH=MIN(0.5_JWRB*DELPHI,10000.0_JWRB)
            ELSEIF (DEPTH(IJ) <= 50.0_JWRB) THEN
              FETCH=MIN(DELPHI,50000.0_JWRB)
            ELSEIF (DEPTH(IJ) <= 150.0_JWRB) THEN
              FETCH=MIN(2*DELPHI,100000.0_JWRB)
            ELSEIF (DEPTH(IJ) <= 250.0_JWRB) THEN
              FETCH=MIN(3*DELPHI,200000.0_JWRB)
            ELSE
              FETCH=MIN(5*DELPHI,250000.0_JWRB)
            ENDIF

!           FIND PEAK PERIOD AND ENERGY LEVEL FROM FETCH LAW
!           THE SAME FORMULATION AS IN SUBROUTINE PEAK IS USED.
            IF (WSWAVE(IJ) > 0.1E-08_JWRB ) THEN
              GX = G * FETCH
              U10 = WSWAVE(IJ)
              GXU = GX/(U10*U10)
              UG = G/U10
              FPK(IJ) = AJONS * GXU ** DJONS
              FPK(IJ) = MAX(0.13_JWRB, FPK(IJ))
              FPK(IJ) = MIN(FPK(IJ), FPMAX/UG)
              ALPHAV(IJ) = BJONS * FPK(IJ)** EJONS
              ALPHAV(IJ) = MAX(ALPHAV(IJ), 0.0081_JWRB)
              FPK(IJ) = FPK(IJ)*UG
            ELSE
              FPK(IJ)=0.0_JWRB
              ALPHAV(IJ)=0.0_JWRB
            ENDIF
          ELSE
            FPK(IJ)=FPMAX 
            ALPHAV(IJ)=0.0_JWRB
          ENDIF
        ENDDO

        CALL JONSWAP(ALPHAV, ZGAMMA, SA, SB, FPK, KIJS, KIJL, ET)

        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              FL1(IJ,K,M) = MAX(FL1(IJ,K,M),ET(IJ,M)*SPRD(IJ,K))
              FLLOWEST = FLMIN*COS2NOISE(IJ,K)
              FL1(IJ,K,M) = MAX(FL1(IJ,K,M),FLLOWEST)
            ENDDO
          ENDDO
        ENDDO

        IF (LBIWBK) THEN
          CALL SDEPTHLIM(KIJS, KIJL, EMAXDPT, FL1)
        ENDIF

      ENDIF


IF (LHOOK) CALL DR_HOOK('UNSETICE',1,ZHOOK_HANDLE)

END SUBROUTINE UNSETICE
