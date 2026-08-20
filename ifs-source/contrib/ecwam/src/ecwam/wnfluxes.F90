! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE WNFLUXES (KIJS, KIJL,                       &
 &                   MIJ, RHOWGDFTH,                   &
 &                   CINV,                             &
 &                   SSURF, CICOVER,                   &
 &                   PHIWA,                            &
 &                   EM, F1, WSWAVE, WDWAVE,           &
 &                   USTRA, VSTRA,                     &
 &                   UFRIC, AIRD,                      &
 &                   NPHIEPS, NTAUOC, NSWH, NMWP,      &
 &                   NEMOTAUX, NEMOTAUY,               &
 &                   NEMOWSWAVE, NEMOPHIF,             &
 &                   TAUXD, TAUYD, TAUOCXD, TAUOCYD,   &
 &                   TAUOC, PHIOCD, PHIEPS, PHIAW,     &
 &                   LNUPD)

! ----------------------------------------------------------------------

!**** *WNFLUXES* - WAVE FLUXES CALCULATION

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *WNFLUXES* (KIJS, KIJL,
!    &                     MIJ, RHOWGDFTH,
!    &                     CINV,
!    &                     SSURF, CICOVER,
!    &                     PHIWA,
!    &                     EM, F1, WSWAVE, WDWAVE,
!    &                     UFRIC, AIRD, INTFLDS, WAM2NEMO,
!    &                     LNUPD)
!          *KIJS*    - INDEX OF FIRST GRIDPOINT.
!          *KIJL*    - INDEX OF LAST GRIDPOINT.
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE
!         *RHOWGDFTH    - WATER DENSITY * G * DF * DTHETA
!                         FOR TRAPEZOIDAL INTEGRATION BETWEEN FR(1) and FR(MIJ)
!                         !!!!!!!!  RHOWGDFTH=0 FOR FR > FR(MIJ)
!          *CINV*   - INVERSE PHASE SPEED.
!          *SSURF*  - CONTRIBUTION OF ALL SOURCE TERMS ACTING ON 
!                     THE SURFACE MOMENTUM AND ENERGY FLUXES.
!          *CICOVER*- SEA ICE COVER.
!          *PHIWA*  - ENERGY FLUX FROM WIND INTO WAVES INTEGRATED
!                     OVER THE FULL FREQUENCY RANGE.
!          *EM*     - MEAN WAVE VARIANCE.
!          *F1*     - MEAN WAVE FREQUENCY BASED ON f*F INTEGRATION.
!          *WSWAVE* - WIND SPEED IN M/S.
!          *WDWAVE* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC CONVENTION
!          *UFRIC*  - FRICTION VELOCITY IN M/S.
!          *AIRD*   - AIR DENSITY IN KG/M3.
!          *INTFLDS*-  INTEGRATED/DERIVED PARAMETERS
!          WAM2NEMO*- WAVE FIELDS PASSED TO NEMO
!          *LNUPD*  - UPDATE NEMO FIELDS.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS, INTGT_PARAM_FIELDS, WAVE2OCEAN

      USE YOWALTAS , ONLY : EGRCRV   ,AFCRV       ,BFCRV
      USE YOWCOUP  , ONLY : LWCOUAST, LWNEMOCOU, LWNEMOTAUOC
      USE YOWFRED  , ONLY : FR       ,COSTH       ,SINTH
      USE YOWICE   , ONLY : LICERUN  ,LWAMRSETCI, CITHRSH, CIBLOCK

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : TAUOCMIN ,TAUOCMAX ,PHIEPSMIN,PHIEPSMAX,    &
     &               EPSUS ,EPSU10   ,G        ,ZPI

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: MIJ

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: RHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: CINV 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: SSURF
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: CICOVER 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: PHIWA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: EM, F1, WSWAVE, WDWAVE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: USTRA, VSTRA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UFRIC, AIRD
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: PHIOCD, PHIEPS, PHIAW
      REAL(KIND=JWRO), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: NPHIEPS, NTAUOC, NSWH, NMWP, NEMOTAUX
      REAL(KIND=JWRO), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: NEMOTAUY, NEMOWSWAVE, NEMOPHIF
      LOGICAL, INTENT(IN) :: LNUPD


      INTEGER(KIND=JWIM) :: IJ, K, M

!     FICTITIOUS VALUE OF THE NORMALISED WAVE ENERGY FLUX UNDER THE SEA ICE 
!     (negative because it is defined as leaving the waves)
      REAL(KIND=JWRB), PARAMETER :: PHIOC_ICE=-3.75_JWRB
      REAL(KIND=JWRB), PARAMETER :: PHIAW_ICE=3.75_JWRB

!     USE HERSBACH 2011 FOR CD(U10) (SEE ALSO EDSON et al. 2013)
      REAL(KIND=JWRB), PARAMETER :: C1 = 1.03E-3_JWRB
      REAL(KIND=JWRB), PARAMETER :: C2 = 0.04E-3_JWRB
      REAL(KIND=JWRB), PARAMETER :: P1 = 1.48_JWRB
      REAL(KIND=JWRB), PARAMETER :: P2 = -0.21_JWRB
      REAL(KIND=JWRB), PARAMETER :: CDMAX = 0.003_JWRB

      REAL(KIND=JWRB), PARAMETER :: EFD_MIN = 0.0625_JWRB  ! corresponds to min Hs=1m under sea ice
      REAL(KIND=JWRB), PARAMETER :: EFD_MAX = 6.25_JWRB    ! corresponds to max Hs=10m under sea ice

      REAL(KIND=JWRB) :: TAU, XN, TAUO
      REAL(KIND=JWRB) :: U10P, CD_BULK, CD_WAVE, CD_ICE 
      REAL(KIND=JWRB) :: CNST
      REAL(KIND=JWRB) :: EPSUS3 
      REAL(KIND=JWRB) :: CITHRSH_INV
      REAL(KIND=JWRB) :: EFD, FFD, EFD_FAC, FFD_FAC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: XSTRESS, YSTRESS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: USTAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: PHILF
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: OOVAL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EM_OC, F1_OC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CMRHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: SUMT, SUMX, SUMY

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WNFLUXES',0,ZHOOK_HANDLE)


      EPSUS3 =  EPSUS*SQRT(EPSUS)

      CITHRSH_INV=1._JWRB/MAX(CITHRSH,0.01_JWRB)

      EFD_FAC = 4.0_JWRB*EGRCRV/G**2 
      FFD_FAC = (EGRCRV/AFCRV)**(1.0_JWRB/BFCRV) * G

!*    DETERMINE NORMALIZED FLUXES FROM AIR TO WAVE AND FROM WAVE TO OCEAN.
!     -------------------------------------------------------------------

!     ENERGY FLUX from SSURF
!     MOMENTUM FLUX FROM SSURF
      DO IJ=KIJS,KIJL
        PHILF(IJ) = 0.0_JWRB
        XSTRESS(IJ) = 0.0_JWRB
        YSTRESS(IJ) = 0.0_JWRB
      ENDDO

!     THE INTEGRATION ONLY UP TO FR=MIJ
      DO M=1,MAXVAL(MIJ(:))
        K=1
        DO IJ=KIJS,KIJL
          SUMT(IJ) = SSURF(IJ,K,M)
          SUMX(IJ) = SINTH(K)*SSURF(IJ,K,M)
          SUMY(IJ) = COSTH(K)*SSURF(IJ,K,M)
        ENDDO
        DO K=2,NANG
          DO IJ=KIJS,KIJL
            SUMT(IJ) = SUMT(IJ) + SSURF(IJ,K,M)
            SUMX(IJ) = SUMX(IJ) + SINTH(K)*SSURF(IJ,K,M)
            SUMY(IJ) = SUMY(IJ) + COSTH(K)*SSURF(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          PHILF(IJ)   = PHILF(IJ)   + SUMT(IJ)*RHOWGDFTH(IJ,M)
          CMRHOWGDFTH(IJ) = CINV(IJ,M)*RHOWGDFTH(IJ,M)
          XSTRESS(IJ) = XSTRESS(IJ) + SUMX(IJ)*CMRHOWGDFTH(IJ)
          YSTRESS(IJ) = YSTRESS(IJ) + SUMY(IJ)*CMRHOWGDFTH(IJ)
        ENDDO
      ENDDO

      IF (LICERUN .AND. LWAMRSETCI) THEN
        DO IJ=KIJS,KIJL
          IF(CICOVER(IJ) > CIBLOCK) THEN
            OOVAL(IJ)=EXP(-MIN((CICOVER(IJ)*CITHRSH_INV)**4,10._JWRB))
!           ADJUST USTAR FOR THE PRESENCE OF SEA ICE
            U10P = MAX(WSWAVE(IJ),EPSU10)
            CD_BULK = MIN((C1 + C2*U10P**P1)*U10P**P2, CDMAX)
            CD_WAVE = (UFRIC(IJ)/U10P)**2
            CD_ICE = OOVAL(IJ)*CD_WAVE + (1.0_JWRB-OOVAL(IJ))*CD_BULK
            USTAR(IJ) = MAX(SQRT(CD_ICE)*U10P,EPSUS)

            ! EM_OC and F1_OC with fully developed model ENERGY 
            ! The significant wave height derived from EM_OC will be used
            ! by NEMO as a scaling factor as if it was open ocean
            EFD = MIN(EFD_FAC*USTAR(IJ)**4, EFD_MAX)
            EM_OC(IJ) = MAX(OOVAL(IJ)*EM(IJ)+(1.0_JWRB-OOVAL(IJ))*EFD, EFD_MIN)
            FFD = FFD_FAC/USTAR(IJ)
            F1_OC(IJ) = OOVAL(IJ)*F1(IJ) + (1.0_JWRB-OOVAL(IJ))*FFD
            F1_OC(IJ) = MIN(MAX(F1_OC(IJ), FR(2)),FR(NFRE))
          ELSE
            OOVAL(IJ) = 1.0_JWRB
            USTAR(IJ) = UFRIC(IJ)
            EM_OC(IJ) = EM(IJ)
            F1_OC(IJ) = F1(IJ)
          ENDIF
        ENDDO
      ELSE
        OOVAL(KIJS:KIJL) = 1.0_JWRB
        USTAR(KIJS:KIJL) = UFRIC(KIJS:KIJL)
        EM_OC(KIJS:KIJL) = EM(KIJS:KIJL)
        F1_OC(KIJS:KIJL) = F1(KIJS:KIJL)
      ENDIF

      DO IJ=KIJS,KIJL
        TAU        = AIRD(IJ)*MAX(USTAR(IJ)**2,EPSUS)
        TAUXD(IJ)  = TAU*SIN(WDWAVE(IJ))
        TAUYD(IJ)  = TAU*COS(WDWAVE(IJ))

        TAUOCXD(IJ)= TAUXD(IJ)-OOVAL(IJ)*XSTRESS(IJ)
        TAUOCYD(IJ)= TAUYD(IJ)-OOVAL(IJ)*YSTRESS(IJ)
        TAUO       = SQRT(TAUOCXD(IJ)**2+TAUOCYD(IJ)**2)
        TAUOC(IJ)  = MIN(MAX(TAUO/TAU,TAUOCMIN),TAUOCMAX)
      ENDDO

      IF (LWCOUAST) THEN
!       USE THE SURFACE STRESSES AS PROVIDED BY THE ATMOSPHERE AND 
!       ONLY RESCALE THEM WITH TAUOC
        WHERE ( USTRA(KIJS:KIJL) /= 0.0_JWRB .OR. VSTRA(KIJS:KIJL) /= 0.0_JWRB )
          TAUXD(KIJS:KIJL)  = USTRA(KIJS:KIJL)
          TAUOCXD(KIJS:KIJL)= USTRA(KIJS:KIJL) * TAUOC(KIJS:KIJL)
          TAUYD(KIJS:KIJL)  = VSTRA(KIJS:KIJL)
          TAUOCYD(KIJS:KIJL)= VSTRA(KIJS:KIJL) * TAUOC(KIJS:KIJL)
        ENDWHERE
      ENDIF

      DO IJ=KIJS,KIJL

        XN        = AIRD(IJ)*MAX(USTAR(IJ)**3,EPSUS3)
        PHIOCD(IJ)= OOVAL(IJ)*(PHILF(IJ)-PHIWA(IJ))+(1.0_JWRB-OOVAL(IJ))*PHIOC_ICE*XN

        PHIEPS(IJ)= PHIOCD(IJ)/XN 
        PHIEPS(IJ)= MIN(MAX(PHIEPS(IJ),PHIEPSMIN),PHIEPSMAX)

        PHIOCD(IJ)= PHIEPS(IJ)*XN

        PHIAW(IJ) = PHIWA(IJ)/XN
        PHIAW(IJ) = OOVAL(IJ)*PHIWA(IJ)/XN+(1.0_JWRB-OOVAL(IJ))*PHIAW_ICE

        IF (LWNEMOCOU .AND. LNUPD) THEN
          NPHIEPS(IJ) = PHIEPS(IJ)
          NTAUOC(IJ)  = TAUOC(IJ)
          IF (EM_OC(IJ) /= 0.0_JWRB) THEN
             NSWH(IJ) = 4.0_JWRO*SQRT(EM_OC(IJ))
          ELSE
             NSWH(IJ) = 0.0_JWRO
          ENDIF
          IF (F1_OC(IJ) /= 0.0_JWRB) THEN
             NMWP(IJ) = 1.0_JWRO/F1_OC(IJ)
          ELSE
             NMWP(IJ) = 0.0_JWRO
          ENDIF

          IF (LWNEMOTAUOC) THEN
            NEMOTAUX(IJ) = NEMOTAUX(IJ) + TAUOCXD(IJ)
            NEMOTAUY(IJ) = NEMOTAUY(IJ) + TAUOCYD(IJ)
          ELSE
            NEMOTAUX(IJ) = NEMOTAUX(IJ) + TAUXD(IJ)
            NEMOTAUY(IJ) = NEMOTAUY(IJ) + TAUYD(IJ)
          ENDIF 
          NEMOWSWAVE(IJ) = NEMOWSWAVE(IJ) + WSWAVE(IJ)
          NEMOPHIF(IJ)   = NEMOPHIF(IJ) + PHIOCD(IJ)
        ENDIF
      ENDDO

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WNFLUXES',1,ZHOOK_HANDLE)

END SUBROUTINE WNFLUXES
