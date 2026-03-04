! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AC_CLOUD_MODEL(                            &
 & YDPHY,YDPHY3,KIDIA    ,KFDIA    ,KLON     ,KTDIA    ,KLEV     , &
 & PEPSNEB  ,PDELP    ,PNEB     ,PQI      ,PQL      , &
 & PR       ,PAPRSF   ,PT       ,PBSFSI   ,PBSFSL   , &
 & PBSFTI   ,PBSFTL   ,PEOASI   ,PEOASL   ,PEOATI   , &
 & PEOATL   ,PEODSI   ,PEODSL   ,PEODTI   ,PEODTL   , &
 & PUSAI    ,PUSAL    ,PUSBI    ,PUSBL                &
 & )

! Purpose:
! --------
!   AC_CLOUD_MODEL - Computes cloud optical coefficients taking into account
!   liquid/ice water content and non-local saturation effect caused by cloud
!   layers above and below.

! Interface:
! ----------
! INPUT:
!   KIDIA   - initial index for horizontal loops
!   KFDIA   - final index for horizontal loops
!   KLON    - horizontal dimension of arrays
!   KTDIA   - initial index for vertical loops (usually 1)
!   KLEV    - vertical dimension of full level arrays
!   PDELP   - pressure difference across layers
!   PNEB    - cloud fraction (0-1)
!   PEPSNEB - treshold for cloud fraction in divisions
!   PQI     - specific mass of ice INSIDE CLOUD
!   PQL     - specific mass of liquid water INSIDE CLOUD
!   PR      - gas constant of air
!   PAPRSF  - full level pressure
!   PT      - temperature

! OUTPUT:
!   PBSF[X][Y] - back scatter fraction beta
!   PEOA[X][Y] - absorption coefficient k_abs
!   PEOD[X][Y] - scattering coefficient k_scat (D = diffuse)
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
!   2006-02, J.-F. Geleyn, J. Masek

! Modifications:
! --------------
!   2006-06, R. Brozkova
!   Code optimization.
!
!   2009-07, K. Yessad
!   Remove CDLOCK + some cleanings.
!   R. El Khatib 05-Feb-2018 fix uninitialized variables
! End Modifications
!-------------------------------------------------------------------------------

USE PARKIND1 , ONLY: JPIM     ,JPRB

USE YOMCST   , ONLY: RG

USE YOMHOOK  , ONLY: LHOOK    ,DR_HOOK, JPHOOK

USE YOMPHY   , ONLY : TPHY

USE YOMPHY3  , ONLY : TPHY3

!-------------------------------------------------------------------------------
                   
IMPLICIT NONE

TYPE(TPHY)           ,INTENT(INOUT):: YDPHY
TYPE(TPHY3)          ,INTENT(INOUT):: YDPHY3
INTEGER(KIND = JPIM), INTENT(IN)  :: KIDIA
INTEGER(KIND = JPIM), INTENT(IN)  :: KFDIA
INTEGER(KIND = JPIM), INTENT(IN)  :: KLON
INTEGER(KIND = JPIM), INTENT(IN)  :: KTDIA
INTEGER(KIND = JPIM), INTENT(IN)  :: KLEV

REAL   (KIND = JPRB), INTENT(IN)  :: PEPSNEB
REAL   (KIND = JPRB), INTENT(IN)  :: PDELP  (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(IN)  :: PNEB   (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(IN)  :: PQI    (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(IN)  :: PQL    (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(IN)  :: PR     (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(IN)  :: PAPRSF (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(IN)  :: PT     (KLON, KLEV)

REAL   (KIND = JPRB), INTENT(OUT) :: PBSFSI (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PBSFSL (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PBSFTI (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PBSFTL (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PEOASI (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PEOASL (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PEOATI (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PEOATL (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PEODSI (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PEODSL (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PEODTI (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PEODTL (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PUSAI  (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PUSAL  (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PUSBI  (KLON, KLEV)
REAL   (KIND = JPRB), INTENT(OUT) :: PUSBL  (KLON, KLEV)

!-------------------------------------------------------------------------------

INTEGER(KIND = JPIM) :: JB, JLEV, JLEV2, JLON

REAL   (KIND = JPRB) :: ZARG, ZARG2
REAL   (KIND = JPRB) :: ZCA, ZCD, ZDELA, ZDELD
REAL   (KIND = JPRB) :: ZDRYEOAI, ZDRYEOAL
REAL   (KIND = JPRB) :: ZDRYEODI, ZDRYEODL
REAL   (KIND = JPRB) :: ZDRYGI, ZDRYGL
REAL   (KIND = JPRB) :: ZEPS, ZEPSNEB1
REAL   (KIND = JPHOOK) :: ZHOOK_HANDLE
REAL   (KIND = JPRB) :: ZMUA, ZMUD
REAL   (KIND = JPRB) :: ZRHO, ZRRG

REAL   (KIND = JPRB) :: ZRNEB1    (KLON)
REAL   (KIND = JPRB) :: ZDEL0     (KLON, KLEV)
REAL   (KIND = JPRB) :: ZDEL0_EFF (KLON, KLEV)
REAL   (KIND = JPRB) :: ZDEL0_EFF2(KLON, KLEV)
REAL   (KIND = JPRB) :: ZEOAI     (KLON, KLEV)
REAL   (KIND = JPRB) :: ZEOAL     (KLON, KLEV)
REAL   (KIND = JPRB) :: ZEODI     (KLON, KLEV)
REAL   (KIND = JPRB) :: ZEODL     (KLON, KLEV)
REAL   (KIND = JPRB) :: ZGI       (KLON, KLEV)
REAL   (KIND = JPRB) :: ZGL       (KLON, KLEV)
REAL   (KIND = JPRB) :: ZIWC1     (KLON, KLEV)
REAL   (KIND = JPRB) :: ZLWC1     (KLON, KLEV)
REAL   (KIND = JPRB) :: ZNEB1     (KLON, KLEV)

REAL   (KIND = JPRB) :: ZGEOM    (KLON, KLEV, KLEV)
 
LOGICAL :: LLQI(KLON, KLEV), LLQL(KLON, KLEV)

!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AC_CLOUD_MODEL',0,ZHOOK_HANDLE)
ASSOCIATE(FCM_Q_DI=>YDPHY3%FCM_Q_DI, USAI=>YDPHY3%USAI, &
 & FCM_Q_AL=>YDPHY3%FCM_Q_AL, FCM_Q_DL=>YDPHY3%FCM_Q_DL, EODSI=>YDPHY3%EODSI, &
 & EODSN=>YDPHY3%EODSN, USAN=>YDPHY3%USAN, FCM_P_DL=>YDPHY3%FCM_P_DL, &
 & FCM_P_DI=>YDPHY3%FCM_P_DI, BSFTI=>YDPHY3%BSFTI, BSFTN=>YDPHY3%BSFTN, &
 & FCM_Q_GI=>YDPHY3%FCM_Q_GI, FCM_N_I=>YDPHY3%FCM_N_I, &
 & FCM_Q_GL=>YDPHY3%FCM_Q_GL, FCM_N_L=>YDPHY3%FCM_N_L, EODTN=>YDPHY3%EODTN, &
 & EODTI=>YDPHY3%EODTI, EOATI=>YDPHY3%EOATI, FCM_Q_AI=>YDPHY3%FCM_Q_AI, &
 & EOATN=>YDPHY3%EOATN, FCM_DEL_D=>YDPHY3%FCM_DEL_D, &
 & FCM_DEL_A=>YDPHY3%FCM_DEL_A, USBN=>YDPHY3%USBN, REXP_NEB=>YDPHY3%REXP_NEB, &
 & FCM_MU_A=>YDPHY3%FCM_MU_A, FCM_MU_D=>YDPHY3%FCM_MU_D, USBI=>YDPHY3%USBI, &
 & FCM_P_GI=>YDPHY3%FCM_P_GI, BSFSN=>YDPHY3%BSFSN, BSFSI=>YDPHY3%BSFSI, &
 & FCM_P_GL=>YDPHY3%FCM_P_GL, FCM_P_AL=>YDPHY3%FCM_P_AL, &
 & N_SPBAND=>YDPHY3%N_SPBAND,&
 & FCM_P_AI=>YDPHY3%FCM_P_AI, EOASN=>YDPHY3%EOASN, EOASI=>YDPHY3%EOASI, &
 & LCLSATUR=>YDPHY%LCLSATUR, LRNUMX=>YDPHY%LRNUMX)
!-------------------------------------------------------------------------------

! 1. Preliminary calculations
! ---------------------------

! logical keys for skipping points without cloud liquid/ice
LLQI(KIDIA:KFDIA, :) = PQI(KIDIA:KFDIA, :) > 0._JPRB
LLQL(KIDIA:KFDIA, :) = PQL(KIDIA:KFDIA, :) > 0._JPRB

! 2. Computation of cloud optical properties
! ------------------------------------------

IF ( .NOT.LCLSATUR ) THEN

  ! 2.1 Old scheme
  ! --------------
  ! Old scheme ignores dependency of cloud optical coefficients on
  ! liquid/ice water content. Only mean saturation effect is taken 
  ! into account.

  ! fill arrays with namelist values
  DO JLEV = KTDIA, KLEV
!DEC$ IVDEP
    DO JLON = KIDIA, KFDIA

      ! absorption coefficient
      PEOASI(JLON, JLEV) = EOASI
      PEOASL(JLON, JLEV) = EOASN
      PEOATI(JLON, JLEV) = EOATI
      PEOATL(JLON, JLEV) = EOATN

      ! scattering coefficient
      PEODSI(JLON, JLEV) = EODSI
      PEODSL(JLON, JLEV) = EODSN
      PEODTI(JLON, JLEV) = EODTI
      PEODTL(JLON, JLEV) = EODTN

      ! back scatter fraction
      PBSFSI(JLON, JLEV) = BSFSI
      PBSFSL(JLON, JLEV) = BSFSN
      PBSFTI(JLON, JLEV) = BSFTI
      PBSFTL(JLON, JLEV) = BSFTN

      ! coefficients for computation of upscatter fraction
      PUSAI (JLON, JLEV) = USAI
      PUSAL (JLON, JLEV) = USAN
      PUSBI (JLON, JLEV) = USBI
      PUSBL (JLON, JLEV) = USBN

    ENDDO
  ENDDO

ELSE

  ! 2.2 New scheme
  ! --------------
  ! New scheme assumes dependency of cloud optical coefficients on
  ! liquid/ice water content. Saturation effect depends on cloud
  ! layers above and below.
  ! It ovewrites the results only in presence of cloud water or ice

  ! auxiliary parameters
  ZEPS     = 1.0E-12_JPRB
  ZEPSNEB1 = PEPSNEB**REXP_NEB
  ZRRG     = 1.0_JPRB/RG

  ! 2.2.1 Geometry factors
  ! ----------------------

  DO JLEV = KTDIA, KLEV
    DO JLON = KIDIA, KFDIA
      ZNEB1(JLON, JLEV) = PNEB(JLON, JLEV)**REXP_NEB
    ENDDO
  ENDDO

  IF ( LRNUMX ) THEN

    DO JLEV = KTDIA, KLEV
      DO JLON = KIDIA, KFDIA
        ZRNEB1(JLON) = 1.0_JPRB/MAX(ZNEB1(JLON, JLEV), ZEPSNEB1)
      ENDDO
      DO JLEV2 = KTDIA, KLEV
        DO JLON = KIDIA, KFDIA
          ZGEOM(JLON, JLEV, JLEV2) = MIN( 1.0_JPRB,&
           & ZNEB1(JLON, JLEV2)*ZRNEB1(JLON) )
        ENDDO
      ENDDO
    ENDDO

  ELSE

    DO JLEV = KTDIA, KLEV
      DO JLEV2 = KTDIA, KLEV
        DO JLON = KIDIA, KFDIA
          ZGEOM(JLON, JLEV, JLEV2) = ZNEB1(JLON, JLEV2)
        ENDDO
      ENDDO
    ENDDO

  ENDIF

  ! reset diagonal elements to one
  DO JLEV = KTDIA, KLEV
    DO JLON = KIDIA, KFDIA
      ZGEOM(JLON, JLEV, JLEV) = 1.0_JPRB
    ENDDO
  ENDDO

  ! 2.2.2 Liquid/ice water content
  ! ------------------------------

  DO JLEV = KTDIA, KLEV
    DO JLON = KIDIA, KFDIA
      ZRHO = PAPRSF(JLON, JLEV)/(PR(JLON, JLEV)*PT(JLON, JLEV))

      ! cloud ice water content
      IF ( LLQI(JLON, JLEV) ) THEN
        ZIWC1(JLON, JLEV) = (ZRHO*PQI(JLON, JLEV))**FCM_N_I
      ELSE
        ZIWC1(JLON, JLEV) = 0._JPRB
      ENDIF

      ! cloud liquid water content
      IF ( LLQL(JLON, JLEV) ) THEN
        ZLWC1(JLON, JLEV) = (ZRHO*PQL(JLON, JLEV))**FCM_N_L
      ELSE
        ZLWC1(JLON, JLEV) = 0._JPRB
      ENDIF

    ENDDO
  ENDDO

  ! loop through spectral bands
  DO JB = 1, N_SPBAND

    ZMUA  = FCM_MU_A(JB)
    ZMUD  = FCM_MU_D(JB)
    ZDELA = ZMUA*LOG(FCM_DEL_A(JB))
    ZDELD = ZMUD*LOG(FCM_DEL_D(JB))

    ! when no ice and no liquid
    ZDRYEOAI = EXP(FCM_P_AI(JB, 0))*ZRRG
    ZDRYEOAL = EXP(FCM_P_AL(JB, 0))*ZRRG
    ZDRYEODI = EXP(FCM_P_DI(JB, 0))*ZRRG
    ZDRYEODL = EXP(FCM_P_DL(JB, 0))*ZRRG
    ZDRYGI =&
      & 1.0_JPRB/( 1.0_JPRB + EXP(-2.0_JPRB*FCM_P_GI(JB, 0)) )
    ZDRYGL =&
      & 1.0_JPRB/( 1.0_JPRB + EXP(-2.0_JPRB*FCM_P_GL(JB, 0)) )


    ! 2.2.3 Pade approximants for k_abs, k_scat and g
    ! -----------------------------------------------

    CALL FIT1(ZIWC1, FCM_P_AI(JB, :), FCM_Q_AI(JB, :), ZEOAI)
    CALL FIT1(ZLWC1, FCM_P_AL(JB, :), FCM_Q_AL(JB, :), ZEOAL)
    CALL FIT1(ZIWC1, FCM_P_DI(JB, :), FCM_Q_DI(JB, :), ZEODI)
    CALL FIT1(ZLWC1, FCM_P_DL(JB, :), FCM_Q_DL(JB, :), ZEODL)
    CALL FIT1(ZIWC1, FCM_P_GI(JB, :), FCM_Q_GI(JB, :), ZGI  )
    CALL FIT1(ZLWC1, FCM_P_GL(JB, :), FCM_Q_GL(JB, :), ZGL  )

    DO JLEV = KTDIA, KLEV
      DO JLON = KIDIA, KFDIA

        ! 2.2.4 Unscaling of k_abs, k_scat and g
        ! --------------------------------------

        ! unscale k_abs, k_scat, convert units to [1/Pa],
        ! unscale asymmetry factor
        IF ( LLQI(JLON, JLEV) ) THEN
          ZEOAI(JLON, JLEV) = EXP(ZEOAI(JLON, JLEV))*ZRRG
          ZEODI(JLON, JLEV) = EXP(ZEODI(JLON, JLEV))*ZRRG
          ZGI  (JLON, JLEV) =&
            & 1.0_JPRB/( 1.0_JPRB + EXP(-2.0_JPRB*ZGI(JLON, JLEV)) )
        ELSE
          ZEOAI(JLON, JLEV) = ZDRYEOAI
          ZEODI(JLON, JLEV) = ZDRYEODI
          ZGI  (JLON, JLEV) = ZDRYGI
        ENDIF
        IF ( LLQL(JLON, JLEV) ) THEN
          ZEOAL(JLON, JLEV) = EXP(ZEOAL(JLON, JLEV))*ZRRG
          ZEODL(JLON, JLEV) = EXP(ZEODL(JLON, JLEV))*ZRRG
          ZGL  (JLON, JLEV) =&
            & 1.0_JPRB/( 1.0_JPRB + EXP(-2.0_JPRB*ZGL(JLON, JLEV)) )
        ELSE
          ZEOAL(JLON, JLEV) = ZDRYEOAL
          ZEODL(JLON, JLEV) = ZDRYEODL
          ZGL  (JLON, JLEV) = ZDRYGL
        ENDIF


        ! 2.2.5 Unsaturated optical depth delta0
        ! --------------------------------------

        ZDEL0(JLON, JLEV) = PDELP(JLON, JLEV)*(&
         & PQI(JLON, JLEV)*(ZEOAI(JLON, JLEV) + ZEODI(JLON, JLEV)) +&
         & PQL(JLON, JLEV)*(ZEOAL(JLON, JLEV) + ZEODL(JLON, JLEV)) )

      ENDDO
    ENDDO

    ! Loop through levels have to be interrupted here since the
    ! computation of delta0_eff requires information from all levels.

    ! loop through levels
    DO JLEV = KTDIA, KLEV

      ! 2.2.6 Effective optical depth delta0_eff
      ! ----------------------------------------

      ! sum effective optical depth: 1, ..., JLEV
      DO JLON = KIDIA, KFDIA
        ZDEL0_EFF(JLON, JLEV) = 0.0_JPRB
      ENDDO
!cdir outerunroll=8
      DO JLEV2 = KTDIA, JLEV
        DO JLON = KIDIA, KFDIA
          ZDEL0_EFF(JLON, JLEV) = ZDEL0_EFF(JLON, JLEV) +&
           & ZGEOM(JLON, JLEV, JLEV2)*ZDEL0(JLON, JLEV2)
        ENDDO
      ENDDO

      ! sum effective optical depth: 1, ..., KLEV
      DO JLON = KIDIA, KFDIA
        ZDEL0_EFF2(JLON, JLEV) = ZDEL0_EFF(JLON, JLEV)
      ENDDO
!cdir outerunroll=8
      DO JLEV2 = JLEV + 1, KLEV
        DO JLON = KIDIA, KFDIA
          ZDEL0_EFF2(JLON, JLEV) = ZDEL0_EFF2(JLON, JLEV) +&
           & ZGEOM(JLON, JLEV, JLEV2)*ZDEL0(JLON, JLEV2)
        ENDDO
      ENDDO

      ! 2.2.7 Saturated k_abs, k_scat
      ! -----------------------------

      DO JLON = KIDIA, KFDIA

        ! saturation factors c_abs, c_scat
        ZARG2 = LOG(MAX(ZDEL0_EFF2(JLON, JLEV), ZEPS))
        IF ( JB == 1 ) THEN

          ! solar band
          ZARG = LOG(MAX(ZDEL0_EFF(JLON, JLEV), ZEPS))
          ZCA = 1.0_JPRB/(1.0_JPRB + EXP(ZMUA*ZARG  - ZDELA))
          ZCA = ZCA*(ZMUA*ZCA + 1.0_JPRB - ZMUA)
          ZCD = 1.0_JPRB/(1.0_JPRB + EXP(ZMUD*ZARG2 - ZDELD))

        ELSE

          ! thermal band
          ZCA = 1.0_JPRB/(1.0_JPRB + EXP(ZMUA*ZARG2 - ZDELA))
          ZCD = 1.0_JPRB/(1.0_JPRB + EXP(ZMUD*ZARG2 - ZDELD))

        ENDIF

        ! saturated absorption coefficient k_abs
        ZEOAI(JLON, JLEV) = ZCA*ZEOAI(JLON, JLEV)
        ZEOAL(JLON, JLEV) = ZCA*ZEOAL(JLON, JLEV)

        ! saturated scattering coefficient k_scat
        ZEODI(JLON, JLEV) = ZCD*ZEODI(JLON, JLEV)
        ZEODL(JLON, JLEV) = ZCD*ZEODL(JLON, JLEV)

      ENDDO
    ! end of loop through levels
    ENDDO

    ! 2.2.8 Fill output arrays
    ! ------------------------

    IF ( JB == 1 ) THEN

      ! solar band
      DO JLEV = KTDIA, KLEV
        DO JLON = KIDIA, KFDIA
          PEOASI(JLON, JLEV) = ZEOAI(JLON, JLEV)
          PEOASL(JLON, JLEV) = ZEOAL(JLON, JLEV)
          PEODSI(JLON, JLEV) = ZEODI(JLON, JLEV)
          PEODSL(JLON, JLEV) = ZEODL(JLON, JLEV)
          PBSFSI(JLON, JLEV) = (4.0_JPRB + ZGI(JLON, JLEV))/&
           & (8.0_JPRB + 8.0_JPRB*ZGI(JLON, JLEV))
          PBSFSL(JLON, JLEV) = (4.0_JPRB + ZGL(JLON, JLEV))/&
           & (8.0_JPRB + 8.0_JPRB*ZGL(JLON, JLEV))

          ! coefficients for computation of upscatter fraction
          PUSAI(JLON, JLEV) = 2.0_JPRB*PBSFSI(JLON, JLEV) - 1.0_JPRB
          PUSAL(JLON, JLEV) = 2.0_JPRB*PBSFSL(JLON, JLEV) - 1.0_JPRB
          PUSBI(JLON, JLEV) = 0.0_JPRB
          PUSBL(JLON, JLEV) = 0.0_JPRB
        ENDDO
      ENDDO

    ELSEIF ( JB == 2 ) THEN

      ! thermal band
      DO JLEV = KTDIA, KLEV
        DO JLON = KIDIA, KFDIA
          PEOATI(JLON, JLEV) = ZEOAI(JLON, JLEV)
          PEOATL(JLON, JLEV) = ZEOAL(JLON, JLEV)
          PEODTI(JLON, JLEV) = ZEODI(JLON, JLEV)
          PEODTL(JLON, JLEV) = ZEODL(JLON, JLEV)
          PBSFTI(JLON, JLEV) = (4.0_JPRB + ZGI(JLON, JLEV))/&
           & (8.0_JPRB + 8.0_JPRB*ZGI(JLON, JLEV))
          PBSFTL(JLON, JLEV) = (4.0_JPRB + ZGL(JLON, JLEV))/&
           & (8.0_JPRB + 8.0_JPRB*ZGL(JLON, JLEV))
        ENDDO
      ENDDO

    ENDIF

  ! end of loop through spectral bands
  ENDDO

ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AC_CLOUD_MODEL',1,ZHOOK_HANDLE)

! ------------------------------------------------------------------------------

CONTAINS

! Private subroutines/functions
! -----------------------------

! P.1 Subroutine for evaluating Pade approximants
! -----------------------------------------------

SUBROUTINE FIT1(PCWC1, PP, PQ, POUT)

! Interface:
! ----------
! INPUT:
!   PCWC1 - scaled liquid or ice water content
!   PP    - Pade coefficients in numerator
!   PQ    - Pade coefficients in denominator

! OUTPUT:
!   POUT - scaled fitted quantity

REAL(KIND = JPRB), INTENT(IN)  :: PCWC1(KLON, KLEV)
REAL(KIND = JPRB), INTENT(IN)  :: PP(0:3)
REAL(KIND = JPRB), INTENT(IN)  :: PQ(1:3)

REAL(KIND = JPRB), INTENT(OUT) :: POUT(KLON, KLEV)

INTEGER(KIND = JPIM) :: JLON, JLEV

REAL(KIND = JPRB) :: ZCWC1, ZP, ZQ

REAL(KIND = JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('AC_CLOUD_MODEL:FIT1',0,ZHOOK_HANDLE)

DO JLEV = KTDIA, KLEV
  DO JLON = KIDIA, KFDIA

    ZCWC1 = PCWC1(JLON, JLEV)

    ZP = PP(0)    + ZCWC1*( PP(1) + ZCWC1*( PP(2) + ZCWC1*PP(3) ) )
    ZQ = 1.0_JPRB + ZCWC1*( PQ(1) + ZCWC1*( PQ(2) + ZCWC1*PQ(3) ) )

    POUT(JLON, JLEV) = ZP/ZQ

  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('AC_CLOUD_MODEL:FIT1',1,ZHOOK_HANDLE)

END SUBROUTINE FIT1

END SUBROUTINE AC_CLOUD_MODEL
