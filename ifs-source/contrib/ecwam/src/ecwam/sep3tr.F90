! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SEP3TR (KIJS, KIJL, FL1, MIJ, WSWAVE, WDWAVE , COSWDIF, &
     &                   ESWELL, FSWELL, THSWELL, FSEA,                  &
     &                   FLSW, SWM,                                      &
     &                   EMTRAIN  ,THTRAIN  ,PMTRAIN)

! ----------------------------------------------------------------------

!**** *SEP3TR* -  COMPUTE ENERGY, DIRECTION AND PERIOD FOR THE 3 WAVE
!****             TRAINS: SWELL1, SWELL2, SWELL3

!     D.PETTENUZZO     MAY 2012
!     JEAN BIDLOT      JUNE 2013  SIMPLIFIED BY REMOVING THE WIND SEA PART
!                                 AND ONLY PARTITIONED THE SWELL SPECTRUM
!                                 AND ONLY LOOKING IN THE PROGNOSTIC RANGE


!*    PURPOSE.
!     --------

!     CREATED TO TEST WAVE SWELL SPECTRA PARTITIONING INTO 3 WAVE TRAINS 

!**   INTERFACE.
!     ----------

!       *CALL* *SEP3TR (KIJS, KIJL, FL1, MIJ, WSWAVE, WDWAVE, COSWDIF,
!                       ESWELL, FSWELL, THSWELL, FSEA,
!                       FLSW, SWM,
!                       EMTRAIN  ,THTRAIN  ,PMTRAIN)
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FL1*    - BLOCK OF FULL SPECTRA
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *WSWAVE* - LATESt WIND SPEED.
!          *WDWAVE* - LATEST WIND DIRECTION.
!          *COSWDIF*- COSINE (WDWAVE - WAVES DIRECTIONS)
!          *ESWELL* - TOTAL SWELL ENERGY
!          *FSWELL* - TOTAL SWELL MEAN FREQUENCY
!          *THSWELL*- TOTAL SWELL MEAN DIRECTION
!          *FSEA*   - WINDSEA MEAN FREQUENCY
!          *FLSW*   - SWELL SPECTRA
!          *SWM*    - ORIGINAL SWELL MASK
!                     IT MIGHT NEED TO BE ADJUSTED
!                     THIS IS POSSIBLE BECAUSE WE SUPPLY SWELL SPECTRA
!                     IN WHICH THE WIND SEA HAS BEEN REMOVED. SOME OF
!                     THE ENERGY LEAKING OUT FROM THE WINDSEA SPECTRUM
!                     MIGHT NOT BE PARTITIONED.
!          *EMTRAIN*- PARTITIONED ENERGY
!          *THTRAIN*- PARTITIONED MEAN DIRECTION
!          *PMTRAIN*- PARTITIONED MEAN PERIOD

!     METHOD.
!     -------

!       HANSON AND PHILLIPS 2001

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : NTRAIN
      USE YOWFRED  , ONLY : FR       ,DFIM     ,C        ,              &
     &            TH       ,COSTH    ,SINTH    ,FRIC
      USE YOWICE   , ONLY : FLMIN
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZPI      ,G        ,EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "fndprt.intfb.h"
#include "semean.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: MIJ
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WSWAVE, WDWAVE 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG), INTENT(IN) :: COSWDIF
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: ESWELL,  FSWELL, THSWELL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: FSEA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FLSW, SWM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NTRAIN), INTENT(OUT) :: EMTRAIN, THTRAIN, PMTRAIN


      INTEGER(KIND=JWIM), PARAMETER :: NPMAX=20
      INTEGER(KIND=JWIM) :: NPMAX_LOC
      INTEGER(KIND=JWIM) :: IJ, M, K, IP
      INTEGER(KIND=JWIM) :: ISORT, I, IPLOC
      INTEGER(KIND=JWIM) :: IFL, IFH, ITHL, ITHH
      INTEGER(KIND=JWIM) :: KM, KP
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL) :: NPEAK, NPK
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL) :: FRINVMIJ
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL) :: MMIN, MMAX 
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL) :: IPNOW 
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL,NTRAIN) :: IENERGY
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL,NPMAX) :: NFRP, NTHP

      ! relative value above max swell value that is considered above noise level
      REAL(KIND=JWRB), PARAMETER :: XNOISELEVEL=0.005_JWRB
      ! minimum Hs (m) value to keep a partition HSMIN=HSMIN_INTER+HSMIN_SLOPE*PERIOD
      REAL(KIND=JWRB), PARAMETER :: HSMIN_INTER=0.05_JWRB
      REAL(KIND=JWRB), PARAMETER :: HSMIN_SLOPE=-0.0017_JWRB
      REAL(KIND=JWRB) :: THRS
      REAL(KIND=JWRB) :: HSMIN

      REAL(KIND=JWRB) :: DELDW
      REAL(KIND=JWRB) :: COSDIR, FRLIMIT 
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: FLLOWEST
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ENEX, SUMETRAIN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ETT
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ENMAX, FLNOISE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: SPRD
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,0:NPMAX) :: DIR, PER, ENE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NTRAIN) :: TEMPDIR, TEMPPER, TEMPENE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: FL, FLLOW

      LOGICAL :: LLEPSMIN
      LOGICAL :: LLADDPART
      LOGICAL, DIMENSION(KIJS:KIJL,NTRAIN) :: LPWSECTOR
      LOGICAL, DIMENSION(KIJS:KIJL,NANG) :: LLCOSDIFF

! ----------------------------------------------------------------------

!*    LOOP THROUGH THE GRID POINTS
!     ----------------------------

      IF (LHOOK) CALL DR_HOOK('SEP3TR',0,ZHOOK_HANDLE)

      DO IJ=KIJS,KIJL
        FRINVMIJ(IJ)=1.0_JWRB/FR(MIJ(IJ))
      ENDDO

      DO IJ=KIJS,KIJL
        ENMAX(IJ)=0.0_JWRB
      ENDDO
!     SMOOTH INPUT SPECTRA (in direction only)
      DO M=1,NFRE
        DO K=1,NANG
          KM=K-1
          IF (KM < 1) KM=NANG
          KP=K+1
          IF (KP > NANG) KP=1
          DO IJ=KIJS,KIJL
!           RE-IMPOSE THE WINDSEA MASK
            IF (FLSW(IJ,K,M) <= 0.0_JWRB) THEN
              FL(IJ,K,M) = 0.0_JWRB
            ELSE
              FL(IJ,K,M) = 0.10_JWRB*(FLSW(IJ,KM,M)+FLSW(IJ,KP,M)) + 0.80_JWRB*FLSW(IJ,K,M) 
              ENMAX(IJ)=MAX(ENMAX(IJ),FL(IJ,K,M))
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO IJ=KIJS,KIJL
        FLNOISE(IJ) = XNOISELEVEL*ENMAX(IJ)
      ENDDO

      DO ISORT=1,NTRAIN
        DO IJ=KIJS,KIJL
          EMTRAIN(IJ,ISORT) = 0.0_JWRB
          PMTRAIN(IJ,ISORT) = 0.0_JWRB
          IF (ESWELL(IJ) > 0.0_JWRB) THEN
            THTRAIN(IJ,ISORT) = THSWELL(IJ)
          ELSE
            THTRAIN(IJ,ISORT) = WDWAVE(IJ)
          ENDIF
        ENDDO
      ENDDO

      DO IJ=KIJS,KIJL
        NPEAK(IJ)=0
      ENDDO

!*    1. DETERMINATES MAXIMA POSITIONS AND NUMBER
!        ----------------------------------------

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          LLCOSDIFF(IJ,K) = (COSWDIF(IJ,K) < -0.4_JWRB)
          SPRD(IJ,K) = MAX(0.0_JWRB, COSWDIF(IJ,K))**2
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            FLLOW(IJ,K,M)=FLMIN
          ENDDO
        ENDDO
      ENDDO


      DO IJ=KIJS,KIJL
        OUT: DO M= 2,MIJ(IJ)-1
          IFL    = MAX ( 1 , M-1 )
          IFH    = MIN ( NFRE , M+1 )
          DO K=1,NANG

            FLLOWEST = MAX(FLLOW(IJ,K,M),FLNOISE(IJ))

            IF (FL(IJ,K,M) > FLLOWEST) THEN
              ITHL   = 1 + MOD(NANG+K-2,NANG)
              ITHH   = 1 + MOD(K,NANG)
              IF ( FL(IJ,ITHL,M   ) > 0.0_JWRB .AND.                 &
     &             FL(IJ,ITHH,M   ) > 0.0_JWRB .AND.                 &
     &             FL(IJ,K   ,IFL ) > 0.0_JWRB .AND.                 &
     &             FL(IJ,K   ,IFH ) > 0.0_JWRB .AND.                 &
     &             FL(IJ,ITHL,IFL ) > 0.0_JWRB .AND.                 &
     &             FL(IJ,ITHL,IFH ) > 0.0_JWRB .AND.                 &
     &             FL(IJ,ITHH,IFL ) > 0.0_JWRB .AND.                 &
     &             FL(IJ,ITHH,IFH ) > 0.0_JWRB ) THEN


                IF ( FL(IJ,K,M) >= FL(IJ,K   ,IFL ) .AND.             &
     &               FL(IJ,K,M) >= FL(IJ,K   ,IFH ) .AND.             &
     &               FL(IJ,K,M) >= FL(IJ,ITHL,IFL ) .AND.             &
     &               FL(IJ,K,M) >= FL(IJ,ITHL,M   ) .AND.             &
     &               FL(IJ,K,M) >= FL(IJ,ITHL,IFH ) .AND.             &
     &               FL(IJ,K,M) >= FL(IJ,ITHH,IFL ) .AND.             &
     &               FL(IJ,K,M) >= FL(IJ,ITHH,M   ) .AND.             &
     &               FL(IJ,K,M) >= FL(IJ,ITHH,IFH ) ) THEN
                    NPEAK(IJ) = NPEAK(IJ) + 1
                    IF (NPEAK(IJ) > NPMAX) EXIT OUT
                    NFRP(IJ,NPEAK(IJ)) = M
                    NTHP(IJ,NPEAK(IJ)) = K
                ENDIF
              ENDIF
            ENDIF

          ENDDO

        ENDDO OUT
      ENDDO

      DO IJ=KIJS,KIJL
        NPEAK(IJ)=MIN(NPEAK(IJ),NPMAX)
      ENDDO

!*    2. GENERATE MASK FOR EACH PARTITION AND COMPUTE STATISTICS
!        -------------------------------------------------------

      CALL FNDPRT(KIJS, KIJL, NPMAX,                  &
     &            NPEAK, MIJ, NTHP, NFRP,             &
     &            FLLOW, LLCOSDIFF, FLNOISE,          &
     &            FL, SWM,                            &
     &            ENE, DIR, PER)

      ! UPDATE SWELL SPECTRUM ??????
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            FLSW(IJ,K,M)=MAX(FL1(IJ,K,M),EPSMIN)*SWM(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

!     TOTAL ENERGY IN THE UPDATED SWELL SPECTRA
      LLEPSMIN=.FALSE.
      CALL SEMEAN (FLSW, KIJS, KIJL, ETT, LLEPSMIN)

!     REMOVE PARTITION WITH Hs < HSMIN 
!     MEAN PERIOD > 1/FR(MIJ)
      DO IJ=KIJS,KIJL
        NPK(IJ)=NPEAK(IJ)
        DO IP=1,NPEAK(IJ)
          HSMIN=HSMIN_INTER+HSMIN_SLOPE*PER(IJ,IP)
          THRS=0.0625_JWRB*HSMIN**2
          IF (ENE(IJ,IP) < THRS .OR. PER(IJ,IP) < FRINVMIJ(IJ)) THEN
            ENE(IJ,IP)=0.
            DIR(IJ,IP)=0.
            PER(IJ,IP)=0.
            NPK(IJ)=NPK(IJ)-1
          ENDIF
        ENDDO
      ENDDO

!     IF NO SWELL PARTITION BUT TOTAL SWELL EXIST THEN
!     ASSIGN FIRST PARTITION TO TOTAL SWELL
!     IF IT HAS SWELL CHARACTERSITICS
      DO IJ=KIJS,KIJL
        IF (NPK(IJ) <= 0 .AND. ESWELL(IJ) > 0.0_JWRB) THEN
          IF (FSWELL(IJ) < FSEA(IJ)) THEN
            NPEAK(IJ)=1
            ENE(IJ,1)=ESWELL(IJ)
            DIR(IJ,1)=THSWELL(IJ)
            PER(IJ,1)=1.0_JWRB/FSWELL(IJ)
          ENDIF 
        ENDIF 
      ENDDO

!*    5. SORT PARTITIONS ACCORDING TO ENERGY AND ASSIGN THE FIRST NTRAIN TRAINS
!        WE ONLY NEED THE FIRST NTRAIN TO BE SORTED
!        -----------------------------------------------------------------

      NPMAX_LOC=MAXVAL(NPEAK(:))

      DO ISORT=1,NTRAIN
        DO IJ=KIJS,KIJL
          IPNOW(IJ)=0
          ENMAX(IJ)=0.0_JWRB
        ENDDO

        DO IP=1,NPMAX_LOC
          DO IJ=KIJS,KIJL
            IF (ENE(IJ,IP) > ENMAX(IJ)) THEN
              IPNOW(IJ) = IP
              ENMAX(IJ) = ENE(IJ,IP)
            ENDIF
          ENDDO
        ENDDO

        DO IJ=KIJS,KIJL
          EMTRAIN(IJ,ISORT)=ENE(IJ,IPNOW(IJ))
          THTRAIN(IJ,ISORT)=DIR(IJ,IPNOW(IJ))
          PMTRAIN(IJ,ISORT)=PER(IJ,IPNOW(IJ))
          ENE(IJ,IPNOW(IJ))=0.0_JWRB
          IENERGY(IJ,ISORT)=MIN(IPNOW(IJ),1)
        ENDDO
      ENDDO


!*    6. PRESERVE TOTAL ENERGY
!        ---------------------

!     6.1 Distribute extra energy proportionally to swell trains

      DO IJ=KIJS,KIJL
        SUMETRAIN(IJ)=MAX(EMTRAIN(IJ,1),EPSMIN)
      ENDDO
      DO ISORT=2,NTRAIN
        DO IJ=KIJS,KIJL
          SUMETRAIN(IJ)=SUMETRAIN(IJ)+EMTRAIN(IJ,ISORT)
        ENDDO
      ENDDO

      DO IJ=KIJS,KIJL
        ENEX(IJ)=MAX((ETT(IJ)-SUMETRAIN(IJ)),0.0_JWRB)/SUMETRAIN(IJ)
      ENDDO

      DO ISORT=1,NTRAIN
        DO IJ=KIJS,KIJL
          EMTRAIN(IJ,ISORT)=EMTRAIN(IJ,ISORT)+ENEX(IJ)*EMTRAIN(IJ,ISORT)
        ENDDO
      ENDDO

!     RE-ORDER SWELL SYSTEM ACCORDING TO TOTAL ENERGY

      DO ISORT=1,NTRAIN
        DO IJ=KIJS,KIJL
          TEMPENE(IJ,ISORT)=EMTRAIN(IJ,ISORT)
          TEMPDIR(IJ,ISORT)=THTRAIN(IJ,ISORT)
          TEMPPER(IJ,ISORT)=PMTRAIN(IJ,ISORT)
        ENDDO
      ENDDO

      DO ISORT=1,NTRAIN
        DO IJ=KIJS,KIJL
          IPNOW(IJ)=0
          ENMAX(IJ)=0.0_JWRB
        ENDDO

        DO IP=1,NTRAIN
          DO IJ=KIJS,KIJL
            IF (TEMPENE(IJ,IP) > ENMAX(IJ)) THEN
              IPNOW(IJ) = IP
              ENMAX(IJ) = TEMPENE(IJ,IP)
            ENDIF
          ENDDO
        ENDDO

        DO IJ=KIJS,KIJL
          IPLOC=MAX(IPNOW(IJ),1)
          EMTRAIN(IJ,ISORT)=TEMPENE(IJ,IPLOC)
          THTRAIN(IJ,ISORT)=TEMPDIR(IJ,IPLOC)
          PMTRAIN(IJ,ISORT)=TEMPPER(IJ,IPLOC)
          TEMPENE(IJ,IPLOC)=0.0_JWRB
        ENDDO
      ENDDO

!     PREPARE OUTPUT

      DO ISORT=1,NTRAIN
        DO IJ=KIJS,KIJL
          IF (IENERGY(IJ,ISORT) == 0) THEN
            EMTRAIN(IJ,ISORT) = 0.0_JWRB
            THTRAIN(IJ,ISORT) = WDWAVE(IJ)
            PMTRAIN(IJ,ISORT) = 0.0_JWRB
          ENDIF
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SEP3TR',1,ZHOOK_HANDLE)

      END SUBROUTINE SEP3TR
