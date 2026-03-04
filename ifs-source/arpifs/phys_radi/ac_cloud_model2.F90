! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AC_CLOUD_MODEL2(                          &
 & YDPHY,YDPHY3,KIDIA    ,KFDIA    ,KLON     ,KTDIA    ,KLEV     ,&
 & KJN      ,KIIDIA   ,KIFDIA   ,KAUCR    ,PDELP    ,&
 & PNEB     ,PQI      ,PQL      ,PR       ,PAPRSF   ,&
 & PT       ,PDEOSA   ,PBSFSI   ,PBSFSL   ,PBSFTI   ,&
 & PBSFTL   ,PEOASI   ,PEOASL   ,PEOATI   ,PEOATL   ,&
 & PEODSI   ,PEODSL   ,PEODTI   ,PEODTL   ,PEOSIDIR ,&
 & PEOSLDIR ,PUSAI    ,PUSAL    ,PUSBI    ,PUSBL     &
 & )

! Purpose:
! --------
!   AC_CLOUD_MODEL2 - Computes cloud optical coefficients taking into
!   account liquid/ice water content and in case of solar absorption
!   also non-local saturation effect caused by cloud layers above.
!   This version is part of ACRANEB2 scheme, original ACRANEB uses
!   frozen version of AC_CLOUD_MODEL.

! Interface:
! ----------
! INPUT:
!   KIDIA   - initial index for horizontal loops
!   KFDIA   - final index for horizontal loops
!   KLON    - horizontal dimension of arrays
!   KTDIA   - initial index for vertical loops (usually 1)
!   KLEV    - vertical dimension of full level arrays
!   KJN     - dimension of arrays containing "daylight" intervals
!   KIIDIA  - array of indices marking start of "daylight" intervals
!   KIFDIA  - array of indices marking end of "daylight" intervals
!   KAUCR   - number of "daylight" intervals
!   PDELP   - pressure difference across layers
!   PNEB    - cloud fraction (0-1)
!   PQI     - specific mass of ice INSIDE CLOUD
!   PQL     - specific mass of liquid water INSIDE CLOUD
!   PR      - gas constant of air
!   PAPRSF  - full level pressure
!   PT      - temperature
!   PDEOSA  - total gaseous optical depths (solar descending)

! OUTPUT:
!   PBSF[X][Y] - back scatter fraction beta
!   PEOA[X][Y] - absorption coefficient k_abs
!   PEOD[X][Y] - scattering coefficient k_scat (D = diffuse)
!   PEOS[Y]DIR - SW extinction coefficient k_ext, not delta-scaled 
!   PUSA[Y]    - coefficient in numerator of upscatter fraction   (solar band)
!   PUSB[Y]    - coefficient in denominator of upscatter fraction (solar band)

!     [X] - spectral band: S - solar
!                          T - thermal
!     [Y] - water phase:   I - ice
!                          L - liquid

! Externals:
! ----------

! Method:
! -------


! Reference:
! ----------


! Author:
! -------
!   2006-02, J.-F. Geleyn, J. Masek (original AC_CLOUD_MODEL)

! Modifications:
! --------------
!   2006-06, R. Brozkova
!   Code optimization.
!
!   2009-07, K. Yessad
!   Remove CDLOCK + some cleanings.
!
!   2013-11, J. Masek
!   Ice clouds based on Edwards et al. 2007 data, unified saturation
!   formula for cloud solar absorption, revised geometry factors for
!   effective cloud optical depth. Upgraded to AC_CLOUD_MODEL2.
!
!   2014-11, J. Masek
!   New ice cloud optical properties, explicit conversion of cloud
!   water content to particle size, separated saturation for ice/liquid clouds,
!   SW gas-cloud overlap parameterisation.
!
!   2016-04, J. Masek
!   Unscaled SW extinction coefficients for true direct solar flux.
!
! End Modifications
!-------------------------------------------------------------------------------

USE PARKIND1 ,ONLY: JPIM     ,JPRB

USE YOMCST   ,ONLY: RG

USE YOMHOOK  ,ONLY: LHOOK    ,DR_HOOK, JPHOOK

USE YOMPHY   ,ONLY : TPHY

USE YOMPHY3  ,ONLY : TPHY3

!-------------------------------------------------------------------------------
                   
IMPLICIT NONE

TYPE(TPHY)        ,INTENT(INOUT):: YDPHY
TYPE(TPHY3)       ,INTENT(INOUT):: YDPHY3
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KTDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KJN
INTEGER(KIND=JPIM),INTENT(IN) :: KIIDIA(KJN)
INTEGER(KIND=JPIM),INTENT(IN) :: KIFDIA(KJN)
INTEGER(KIND=JPIM),INTENT(IN) :: KAUCR

REAL(KIND=JPRB),INTENT(IN)  :: PDELP   (KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)  :: PNEB    (KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)  :: PQI     (KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)  :: PQL     (KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)  :: PR      (KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)  :: PAPRSF  (KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)  :: PT      (KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)  :: PDEOSA  (KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PBSFSI  (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PBSFSL  (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PBSFTI  (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PBSFTL  (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PEOASI  (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PEOASL  (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PEOATI  (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PEOATL  (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PEODSI  (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PEODSL  (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PEODTI  (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PEODTL  (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PEOSIDIR(KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PEOSLDIR(KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PUSAI   (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PUSAL   (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PUSBI   (KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PUSBL   (KLON,KLEV)

!-------------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IAUCR
INTEGER(KIND=JPIM) :: JB,JLEV,JLEV2,JLON,JN

INTEGER(KIND=JPIM) :: IIDIA(KJN),IFDIA(KJN)

REAL(KIND=JPRB) :: ZARGI,ZARGL,ZA,ZB,ZCI,ZCL,ZGI0,ZGL0
REAL(KIND=JPRB) :: ZAI2L,ZDI2L,ZDE,ZRE,ZRHO,ZRRG,ZTRLI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZIWC      (KLON,KLEV)
REAL(KIND=JPRB) :: ZLWC      (KLON,KLEV)
REAL(KIND=JPRB) :: ZDEL0_EFFA(KLON,KLEV)
REAL(KIND=JPRB) :: ZDEL0_EFFD(KLON,KLEV)
REAL(KIND=JPRB) :: ZDE1      (KLON,KLEV,YDPHY3%N_SPBAND)
REAL(KIND=JPRB) :: ZRE1      (KLON,KLEV,YDPHY3%N_SPBAND)
REAL(KIND=JPRB) :: ZDEL0     (KLON,KLEV,YDPHY3%N_SPBAND)
REAL(KIND=JPRB) :: ZEOAI     (KLON,KLEV,YDPHY3%N_SPBAND)
REAL(KIND=JPRB) :: ZEOAL     (KLON,KLEV,YDPHY3%N_SPBAND)
REAL(KIND=JPRB) :: ZEODI     (KLON,KLEV,YDPHY3%N_SPBAND)
REAL(KIND=JPRB) :: ZEODL     (KLON,KLEV,YDPHY3%N_SPBAND)
REAL(KIND=JPRB) :: ZGI       (KLON,KLEV,YDPHY3%N_SPBAND)
REAL(KIND=JPRB) :: ZGL       (KLON,KLEV,YDPHY3%N_SPBAND)

LOGICAL :: LLQI(KLON,KLEV),LLQL(KLON,KLEV)

!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AC_CLOUD_MODEL2',0,ZHOOK_HANDLE)
ASSOCIATE(FCM_Q_DI=>YDPHY3%FCM_Q_DI, FCM_Q_DL=>YDPHY3%FCM_Q_DL, &
 & FCM_P_DL=>YDPHY3%FCM_P_DL, FCM_P_DI=>YDPHY3%FCM_P_DI, FCM_AI=>YDPHY3%FCM_AI, &
 & FCM_Q_GI=>YDPHY3%FCM_Q_GI, FCM_B_BI=>YDPHY3%FCM_B_BI, FCM_AL=>YDPHY3%FCM_AL, &
 & FCM_Q_GL=>YDPHY3%FCM_Q_GL, FCM_IWC2DE=>YDPHY3%FCM_IWC2DE, &
 & EODTN=>YDPHY3%EODTN, EODTI=>YDPHY3%EODTI, EOATI=>YDPHY3%EOATI, &
 & FCM_Q_AI=>YDPHY3%FCM_Q_AI, EOATN=>YDPHY3%EOATN, USBN=>YDPHY3%USBN, &
 & FCM_B_BL=>YDPHY3%FCM_B_BL, USBI=>YDPHY3%USBI, FCM_P_GI=>YDPHY3%FCM_P_GI, &
 & BSFSN=>YDPHY3%BSFSN, BSFSI=>YDPHY3%BSFSI, FCM_P_GL=>YDPHY3%FCM_P_GL, &
 & EOASN=>YDPHY3%EOASN, EOASI=>YDPHY3%EOASI, USAI=>YDPHY3%USAI, &
 & USAN=>YDPHY3%USAN, FCM_DEL_DI=>YDPHY3%FCM_DEL_DI, &
 & FCM_DEL_DL=>YDPHY3%FCM_DEL_DL, N_SPBAND=>YDPHY3%N_SPBAND, &
 & EODSN=>YDPHY3%EODSN, FCM_NU_DI=>YDPHY3%FCM_NU_DI, FCM_B_AI=>YDPHY3%FCM_B_AI, &
 & FCM_NU_DL=>YDPHY3%FCM_NU_DL, FCM_B_AL=>YDPHY3%FCM_B_AL, BSFTI=>YDPHY3%BSFTI, &
 & BSFTN=>YDPHY3%BSFTN, FCM_MU_DI=>YDPHY3%FCM_MU_DI, &
 & FCM_LWC2RE=>YDPHY3%FCM_LWC2RE, FCM_MU_DL=>YDPHY3%FCM_MU_DL, &
 & EODSI=>YDPHY3%EODSI, FCM_DEL_AL=>YDPHY3%FCM_DEL_AL, &
 & FCM_Q_AL=>YDPHY3%FCM_Q_AL, FCM_DEL_AI=>YDPHY3%FCM_DEL_AI, &
 & FCM_NU_AI=>YDPHY3%FCM_NU_AI, FCM_NU_AL=>YDPHY3%FCM_NU_AL, &
 & FCM_P_AL=>YDPHY3%FCM_P_AL, FCM_P_AI=>YDPHY3%FCM_P_AI, &
 & FCM_MU_AI=>YDPHY3%FCM_MU_AI, FCM_MU_AL=>YDPHY3%FCM_MU_AL, &
 & LCLSATUR=>YDPHY%LCLSATUR)
!-------------------------------------------------------------------------------

! 1. Computation of cloud optical properties
! ------------------------------------------

IF ( .NOT.LCLSATUR ) THEN

  ! 1.1 Old scheme
  ! --------------
  ! Old scheme ignores dependency of cloud optical coefficients on
  ! liquid/ice water content. Only mean saturation effect is taken 
  ! into account.

  ! fill arrays with namelist values
  DO JLEV=KTDIA,KLEV
    DO JLON=KIDIA,KFDIA

      ! absorption coefficient
      PEOASI(JLON,JLEV)=EOASI
      PEOASL(JLON,JLEV)=EOASN
      PEOATI(JLON,JLEV)=EOATI
      PEOATL(JLON,JLEV)=EOATN

      ! scattering coefficient
      PEODSI(JLON,JLEV)=EODSI
      PEODSL(JLON,JLEV)=EODSN
      PEODTI(JLON,JLEV)=EODTI
      PEODTL(JLON,JLEV)=EODTN

      ! back scatter fraction
      PBSFSI(JLON,JLEV)=BSFSI
      PBSFSL(JLON,JLEV)=BSFSN
      PBSFTI(JLON,JLEV)=BSFTI
      PBSFTL(JLON,JLEV)=BSFTN

      ! coefficients for computation of upscatter fraction
      PUSAI (JLON,JLEV)=USAI
      PUSAL (JLON,JLEV)=USAN
      PUSBI (JLON,JLEV)=USBI
      PUSBL (JLON,JLEV)=USBN

      PEOSIDIR(JLON,JLEV)=0._JPRB
      PEOSLDIR(JLON,JLEV)=0._JPRB

    ENDDO
  ENDDO

ELSE

  ! 1.2 New scheme
  ! --------------
  ! New scheme assumes dependency of cloud optical coefficients on
  ! liquid/ice water content. It introduces solar cloud saturation
  ! depending on cloud layers above and below.

  ! auxiliary parameters
  ZTRLI=EXP(-250._JPRB)
  ZRRG =1._JPRB/RG
  ZAI2L=LOG(FCM_DEL_AI/FCM_DEL_AL)
  ZDI2L=LOG(FCM_DEL_DI/FCM_DEL_DL)

  ! logical keys for skipping points without cloud liquid/ice
  LLQI(:,:)=PQI(:,:) > 0._JPRB
  LLQL(:,:)=PQL(:,:) > 0._JPRB

  ! 1.2.1 Convert IWC/LWC to D_e/R_e
  ! --------------------------------

  ! determine IWC/LWC [g/m^3]
  DO JLEV=KTDIA,KLEV
    DO JLON=KIDIA,KFDIA
      ZRHO=PAPRSF(JLON,JLEV)/(PR(JLON,JLEV)*PT(JLON,JLEV))
      ZIWC(JLON,JLEV)=ZRHO*PQI(JLON,JLEV)*1000._JPRB
      ZLWC(JLON,JLEV)=ZRHO*PQL(JLON,JLEV)*1000._JPRB
    ENDDO
  ENDDO

  ! determine D_e/R_e [micron] and scale it for fitting
  DO JB=1,N_SPBAND
    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA

        ! effective dimension D_e of ice particles [micron]
        IF ( LLQI(JLON,JLEV) ) THEN
          ZDE=FCM_IWC2DE(JB,0)+FCM_IWC2DE(JB,1)*&
           & (ZIWC(JLON,JLEV)+FCM_IWC2DE(JB,2))**FCM_IWC2DE(JB,3)
          ZDE=MAX(FCM_IWC2DE(JB,-2),MIN(FCM_IWC2DE(JB,-1),ZDE))
          ZDE1(JLON,JLEV,JB)=LOG(MAX(ZDE,ZTRLI))
        ELSE
          ZDE1(JLON,JLEV,JB)=0._JPRB
        ENDIF

        ! effective radius R_e of water droplets [micron]
        IF ( LLQL(JLON,JLEV) ) THEN
          ZRE=FCM_LWC2RE(JB,0)+FCM_LWC2RE(JB,1)*&
           & (ZLWC(JLON,JLEV)+FCM_LWC2RE(JB,2))**FCM_LWC2RE(JB,3)
          ZRE=MAX(FCM_LWC2RE(JB,-2),MIN(FCM_LWC2RE(JB,-1),ZRE))
          ZRE1(JLON,JLEV,JB)=LOG(MAX(ZRE,ZTRLI))
        ELSE
          ZRE1(JLON,JLEV,JB)=0._JPRB
        ENDIF

      ENDDO
    ENDDO
  ENDDO

  ! 1.2.2 Unsaturated cloud optical properties
  ! ------------------------------------------

  ! initialize arrays for the case of no clouds
  ZEOAI(:,:,:)=0._JPRB
  ZEOAL(:,:,:)=0._JPRB
  ZEODI(:,:,:)=0._JPRB
  ZEODL(:,:,:)=0._JPRB
  ZGI  (:,:,:)=0._JPRB
  ZGL  (:,:,:)=0._JPRB

  ! loop through spectral bands
  DO JB=1,N_SPBAND

    ! differentiate between solar/thermal bounds
    IF ( JB == 1 ) THEN
      IAUCR=KAUCR
      IIDIA(1:KAUCR)=KIIDIA(1:KAUCR)
      IFDIA(1:KAUCR)=KIFDIA(1:KAUCR)
    ELSE
      IAUCR=1
      IIDIA(1)=KIDIA
      IFDIA(1)=KFDIA
    ENDIF

    ! Pade approximants for scaled k_abs, k_scat and g
    CALL FIT1(ZDE1(1,1,JB),FCM_P_AI(JB,:),FCM_Q_AI(JB,:),ZEOAI(1,1,JB))
    CALL FIT1(ZRE1(1,1,JB),FCM_P_AL(JB,:),FCM_Q_AL(JB,:),ZEOAL(1,1,JB))
    CALL FIT1(ZDE1(1,1,JB),FCM_P_DI(JB,:),FCM_Q_DI(JB,:),ZEODI(1,1,JB))
    CALL FIT1(ZRE1(1,1,JB),FCM_P_DL(JB,:),FCM_Q_DL(JB,:),ZEODL(1,1,JB))
    CALL FIT1(ZDE1(1,1,JB),FCM_P_GI(JB,:),FCM_Q_GI(JB,:),ZGI  (1,1,JB))
    CALL FIT1(ZRE1(1,1,JB),FCM_P_GL(JB,:),FCM_Q_GL(JB,:),ZGL  (1,1,JB))

    DO JLEV=KTDIA,KLEV
      DO JN=1,IAUCR
        DO JLON=IIDIA(JN),IFDIA(JN)

          ! unscale k_abs, k_scat, convert units to [1/Pa],
          ! unscale asymmetry factor
          IF ( LLQI(JLON,JLEV) ) THEN
            ZEOAI(JLON,JLEV,JB)=EXP(ZEOAI(JLON,JLEV,JB))*ZRRG
            ZEODI(JLON,JLEV,JB)=EXP(ZEODI(JLON,JLEV,JB))*ZRRG
            ZGI  (JLON,JLEV,JB)=0.5_JPRB/&
             & (1.0_JPRB+EXP(-2.0_JPRB*ZGI(JLON,JLEV,JB)))
          ENDIF
          IF ( LLQL(JLON,JLEV) ) THEN
            ZEOAL(JLON,JLEV,JB)=EXP(ZEOAL(JLON,JLEV,JB))*ZRRG
            ZEODL(JLON,JLEV,JB)=EXP(ZEODL(JLON,JLEV,JB))*ZRRG
            ZGL  (JLON,JLEV,JB)=0.5_JPRB/&
             & (1.0_JPRB+EXP(-2.0_JPRB*ZGL(JLON,JLEV,JB)))
          ENDIF

          ! unsaturated optical depth delta0
          ZDEL0(JLON,JLEV,JB)=PDELP(JLON,JLEV)*(&
           & PQI(JLON,JLEV)*(ZEOAI(JLON,JLEV,JB)+ZEODI(JLON,JLEV,JB))+&
           & PQL(JLON,JLEV)*(ZEOAL(JLON,JLEV,JB)+ZEODL(JLON,JLEV,JB)))

        ENDDO
      ENDDO
    ENDDO

  ! end of loop through spectral bands
  ENDDO

  ! 1.2.3 Solar saturation
  ! ----------------------

  ! solar band
  JB=1

  ! loop through levels
  DO JLEV=KTDIA,KLEV

    ! initialize effective cloud optical depth with local one
    DO JN=1,KAUCR
      DO JLON=KIIDIA(JN),KIFDIA(JN)
        ZDEL0_EFFA(JLON,JLEV)=ZDEL0(JLON,JLEV,JB)
        ZDEL0_EFFD(JLON,JLEV)=ZDEL0(JLON,JLEV,JB)
      ENDDO
    ENDDO

    ! sum effective cloud optical depth: 1, ..., JLEV-1
!cdir outerunroll=8
    DO JLEV2=KTDIA,JLEV-1
      DO JN=1,KAUCR
        DO JLON=KIIDIA(JN),KIFDIA(JN)
          ZB=(FCM_B_AI*ZIWC(JLON,JLEV)+FCM_B_AL*ZLWC(JLON,JLEV))/&
           & (ZIWC(JLON,JLEV)+ZLWC(JLON,JLEV)+ZTRLI)
          ZDEL0_EFFA(JLON,JLEV)=ZDEL0_EFFA(JLON,JLEV)+&
           & PNEB(JLON,JLEV2)*ZDEL0(JLON,JLEV2,JB)*ZB
          ZDEL0_EFFD(JLON,JLEV)=ZDEL0_EFFD(JLON,JLEV)+&
           & PNEB(JLON,JLEV2)*ZDEL0(JLON,JLEV2,JB)
        ENDDO
      ENDDO
    ENDDO

    ! sum effective cloud optical depth: JLEV+1, ..., KLEV
!cdir outerunroll=8
    DO JLEV2=JLEV+1,KLEV
      DO JN=1,KAUCR
        DO JLON=KIIDIA(JN),KIFDIA(JN)
          ZB=(FCM_B_BI*ZIWC(JLON,JLEV)+FCM_B_BL*ZLWC(JLON,JLEV))/&
           & (ZIWC(JLON,JLEV)+ZLWC(JLON,JLEV)+ZTRLI)
          ZDEL0_EFFA(JLON,JLEV)=ZDEL0_EFFA(JLON,JLEV)+&
           & PNEB(JLON,JLEV2)*ZDEL0(JLON,JLEV2,JB)*ZB
          ZDEL0_EFFD(JLON,JLEV)=ZDEL0_EFFD(JLON,JLEV)+&
           & PNEB(JLON,JLEV2)*ZDEL0(JLON,JLEV2,JB)
        ENDDO
      ENDDO
    ENDDO

    ! saturation of k_abs
    DO JN=1,KAUCR
      DO JLON=KIIDIA(JN),KIFDIA(JN)

        ! saturation factors
        ZA   =(FCM_AI*ZIWC(JLON,JLEV)+FCM_AL*ZLWC(JLON,JLEV))/&
         &    (ZIWC(JLON,JLEV)+ZLWC(JLON,JLEV)+ZTRLI)
        ZARGI=LOG(MAX((ZDEL0_EFFA(JLON,JLEV)+&
         &    ZA*PDEOSA(JLON,JLEV))/FCM_DEL_AI,ZTRLI))
        ZARGL=ZARGI+ZAI2L
        ZCI  =EXP(-FCM_NU_AI*LOG(1.0_JPRB+EXP(FCM_MU_AI*ZARGI)))
        ZCL  =EXP(-FCM_NU_AL*LOG(1.0_JPRB+EXP(FCM_MU_AL*ZARGL)))

        ! saturated k_abs
        ZEOAI(JLON,JLEV,JB)=ZCI*ZEOAI(JLON,JLEV,JB)
        ZEOAL(JLON,JLEV,JB)=ZCL*ZEOAL(JLON,JLEV,JB)

      ENDDO
    ENDDO

    IF ( FCM_NU_DI /= 0._JPRB .OR. FCM_NU_DL /= 0._JPRB ) THEN

      ! saturation of k_scat
      DO JN=1,KAUCR
        DO JLON=KIIDIA(JN),KIFDIA(JN)

          ! saturation factors
          ZARGI=LOG(MAX(ZDEL0_EFFD(JLON,JLEV)/FCM_DEL_DI,ZTRLI))
          ZARGL=ZARGI+ZDI2L
          ZCI  =EXP(-FCM_NU_DI*LOG(1.0_JPRB+EXP(FCM_MU_DI*ZARGI)))
          ZCL  =EXP(-FCM_NU_DL*LOG(1.0_JPRB+EXP(FCM_MU_DL*ZARGL)))

          ! saturated k_scat
          ZEODI(JLON,JLEV,JB)=ZCI*ZEODI(JLON,JLEV,JB)
          ZEODL(JLON,JLEV,JB)=ZCL*ZEODL(JLON,JLEV,JB)

        ENDDO
      ENDDO

    ENDIF

  ! end of loop through levels
  ENDDO

  ! 1.2.4 Fill output arrays
  ! ------------------------

  ! solar band
  JB=1
  DO JLEV=KTDIA,KLEV
    DO JN=1,KAUCR
      DO JLON=KIIDIA(JN),KIFDIA(JN)
        PEOASI(JLON,JLEV)=ZEOAI(JLON,JLEV,JB)
        PEOASL(JLON,JLEV)=ZEOAL(JLON,JLEV,JB)
        PEODSI(JLON,JLEV)=ZEODI(JLON,JLEV,JB)
        PEODSL(JLON,JLEV)=ZEODL(JLON,JLEV,JB)
        PBSFSI(JLON,JLEV)=0.5_JPRB-0.375_JPRB*ZGI(JLON,JLEV,JB)
        PBSFSL(JLON,JLEV)=0.5_JPRB-0.375_JPRB*ZGL(JLON,JLEV,JB)

        ! extinction coefficients, not detla-scaled
        ZGI0=ZGI(JLON,JLEV,JB)/(1._JPRB-ZGI(JLON,JLEV,JB))
        ZGL0=ZGL(JLON,JLEV,JB)/(1._JPRB-ZGL(JLON,JLEV,JB))
        PEOSIDIR(JLON,JLEV)=PEOASI(JLON,JLEV)+&
         &                  PEODSI(JLON,JLEV)/(1._JPRB-ZGI0*ZGI0)
        PEOSLDIR(JLON,JLEV)=PEOASL(JLON,JLEV)+&
         &                  PEODSL(JLON,JLEV)/(1._JPRB-ZGL0*ZGL0)

        ! coefficients for computation of upscatter fraction
        PUSAI(JLON,JLEV)=2.0_JPRB*PBSFSI(JLON,JLEV)-1.0_JPRB
        PUSAL(JLON,JLEV)=2.0_JPRB*PBSFSL(JLON,JLEV)-1.0_JPRB
        PUSBI(JLON,JLEV)=0.0_JPRB
        PUSBL(JLON,JLEV)=0.0_JPRB
      ENDDO
    ENDDO
  ENDDO

  ! thermal band
  JB=2
  DO JLEV=KTDIA,KLEV
    DO JLON=KIDIA,KFDIA
      PEOATI(JLON,JLEV)=ZEOAI(JLON,JLEV,JB)
      PEOATL(JLON,JLEV)=ZEOAL(JLON,JLEV,JB)
      PEODTI(JLON,JLEV)=ZEODI(JLON,JLEV,JB)
      PEODTL(JLON,JLEV)=ZEODL(JLON,JLEV,JB)
      PBSFTI(JLON,JLEV)=0.5_JPRB-0.375_JPRB*ZGI(JLON,JLEV,JB)
      PBSFTL(JLON,JLEV)=0.5_JPRB-0.375_JPRB*ZGL(JLON,JLEV,JB)
    ENDDO
  ENDDO

ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AC_CLOUD_MODEL2',1,ZHOOK_HANDLE)

! ------------------------------------------------------------------------------

CONTAINS

! Private subroutines/functions
! -----------------------------

! P.1 Subroutine for evaluating Pade approximants
! -----------------------------------------------

SUBROUTINE FIT1(PSIZE,PP,PQ,POUT)

! Interface:
! ----------
! INPUT:
!   PSIZE - scaled D_e or R_e
!   PP    - Pade coefficients in numerator 
!   PQ    - Pade coefficients in denominator

! OUTPUT:
!   POUT - scaled fitted quantity

REAL(KIND=JPRB),INTENT(IN)  :: PSIZE(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)  :: PP(0:3)
REAL(KIND=JPRB),INTENT(IN)  :: PQ(1:3)
REAL(KIND=JPRB),INTENT(OUT) :: POUT(KLON,KLEV)

INTEGER(KIND=JPIM) :: JLON,JLEV,JN

REAL(KIND=JPRB) :: ZSIZE,ZP,ZQ

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('AC_CLOUD_MODEL2:FIT1',0,ZHOOK_HANDLE)
ASSOCIATE(FCM_Q_DI=>YDPHY3%FCM_Q_DI, FCM_Q_DL=>YDPHY3%FCM_Q_DL, &
 & FCM_P_DL=>YDPHY3%FCM_P_DL, FCM_P_DI=>YDPHY3%FCM_P_DI, FCM_AI=>YDPHY3%FCM_AI, &
 & FCM_Q_GI=>YDPHY3%FCM_Q_GI, FCM_B_BI=>YDPHY3%FCM_B_BI, FCM_AL=>YDPHY3%FCM_AL, &
 & FCM_Q_GL=>YDPHY3%FCM_Q_GL, FCM_IWC2DE=>YDPHY3%FCM_IWC2DE, &
 & EODTN=>YDPHY3%EODTN, EODTI=>YDPHY3%EODTI, EOATI=>YDPHY3%EOATI, &
 & FCM_Q_AI=>YDPHY3%FCM_Q_AI, EOATN=>YDPHY3%EOATN, USBN=>YDPHY3%USBN, &
 & FCM_B_BL=>YDPHY3%FCM_B_BL, USBI=>YDPHY3%USBI, FCM_P_GI=>YDPHY3%FCM_P_GI, &
 & BSFSN=>YDPHY3%BSFSN, BSFSI=>YDPHY3%BSFSI, FCM_P_GL=>YDPHY3%FCM_P_GL, &
 & EOASN=>YDPHY3%EOASN, EOASI=>YDPHY3%EOASI, USAI=>YDPHY3%USAI, &
 & USAN=>YDPHY3%USAN, FCM_DEL_DI=>YDPHY3%FCM_DEL_DI, &
 & FCM_DEL_DL=>YDPHY3%FCM_DEL_DL, N_SPBAND=>YDPHY3%N_SPBAND, &
 & EODSN=>YDPHY3%EODSN, FCM_NU_DI=>YDPHY3%FCM_NU_DI, FCM_B_AI=>YDPHY3%FCM_B_AI, &
 & FCM_NU_DL=>YDPHY3%FCM_NU_DL, FCM_B_AL=>YDPHY3%FCM_B_AL, BSFTI=>YDPHY3%BSFTI, &
 & BSFTN=>YDPHY3%BSFTN, FCM_MU_DI=>YDPHY3%FCM_MU_DI, &
 & FCM_LWC2RE=>YDPHY3%FCM_LWC2RE, FCM_MU_DL=>YDPHY3%FCM_MU_DL, &
 & EODSI=>YDPHY3%EODSI, FCM_DEL_AL=>YDPHY3%FCM_DEL_AL, &
 & FCM_Q_AL=>YDPHY3%FCM_Q_AL, FCM_DEL_AI=>YDPHY3%FCM_DEL_AI, &
 & FCM_NU_AI=>YDPHY3%FCM_NU_AI, FCM_NU_AL=>YDPHY3%FCM_NU_AL, &
 & FCM_P_AL=>YDPHY3%FCM_P_AL, FCM_P_AI=>YDPHY3%FCM_P_AI, &
 & FCM_MU_AI=>YDPHY3%FCM_MU_AI, FCM_MU_AL=>YDPHY3%FCM_MU_AL, &
 & LCLSATUR=>YDPHY%LCLSATUR)

DO JLEV=KTDIA,KLEV
  DO JN=1,IAUCR
    DO JLON=IIDIA(JN),IFDIA(JN)
      ZSIZE=PSIZE(JLON,JLEV)
      ZP=PP(0)   +ZSIZE*(PP(1)+ZSIZE*(PP(2)+ZSIZE*PP(3)))
      ZQ=1.0_JPRB+ZSIZE*(PQ(1)+ZSIZE*(PQ(2)+ZSIZE*PQ(3)))
      POUT(JLON,JLEV)=ZP/ZQ
    ENDDO
  ENDDO
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AC_CLOUD_MODEL2:FIT1',1,ZHOOK_HANDLE)

END SUBROUTINE FIT1

END SUBROUTINE AC_CLOUD_MODEL2
