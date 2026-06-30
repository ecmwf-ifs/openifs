! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------
MODULE FLAKE_DRIVER_MOD
CONTAINS
SUBROUTINE FLAKE_DRIVER( KIDIA, KFDIA, KLON,         &
                  & LDLAKEPOINT,                     &
                  & PSSRFLTI,PSLRFL,PAHFSTI,PEVAPTI, &
                  & PUSTRTI,PVSTRTI,PFRTI,           &
                  & PLDEPT,PGEMU,PTSTP,              &
                  & YDCST,YDFLAKE,                   & 

                  & PTLICEM1M,PTLMNWM1M,PTLWMLM1M,PTLBOTM1M,PTLSFM1M, & 
                  & PHLICEM1M,PHLMLM1M,                               &

                  & PTLICE,PTLMNW,PTLWML,                             &
                  & PTLBOT,PTLSF,PHLICE,PHLML )
                  
                  

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!------------------------------------------------------------------------------
!
! Description:
!
!  The FLAKE_DRIVER is a communication routine between HTESSEL
!  and a FLake routines.
!  It assigns the FLake variables at the previous time step 
!  to their input values given by the driving model,
!  calls a routine FLAKE_RADFLUX to compute the radiation fluxes,
!  calls FLAKE_ENE,
!  and returns the updated FLake variables to the driving model.
!  The FLAKE_DRIVER does not contain any Flake physics. 
!  It only serves as a convenient means to organize calls of Flake routines.

! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release
! 1.01T      2008/03/10 Victor Stepanenko
!  The code is accomodated to be used in HTESSEL model
!      F. Vana  05-Mar-2015  Support for single precision
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST  , ONLY : TCST
USE YOS_FLAKE, ONLY : TFLAKE, ROPTICPAR_MEDIUM

USE FLAKERAD_MOD
USE FLAKEENE_MOD

 
!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Input (procedure arguments)

INTEGER (KIND = JPIM), INTENT(IN):: KIDIA                            
INTEGER (KIND = JPIM), INTENT(IN):: KFDIA                            
INTEGER (KIND = JPIM), INTENT(IN):: KLON


LOGICAL, INTENT(IN):: LDLAKEPOINT  (:)  ! Indicates, if the point is a lake point (.TRUE.) or not (.FALSE.)

REAL (KIND = JPRB), INTENT(IN) :: PSSRFLTI (:) ! Net solar radiation flux at the surface [W m^{-2}]
REAL (KIND = JPRB), INTENT(IN) :: PSLRFL   (:) ! Net long-wave radiation flux at the surface [W m^{-2}]
REAL (KIND = JPRB), INTENT(IN) :: PAHFSTI  (:) ! Sensible heat flux downwards [W m^{-2}]
REAL (KIND = JPRB), INTENT(IN) :: PEVAPTI  (:) ! Latent heat flux downwards [W m^{-2}]
REAL (KIND = JPRB), INTENT(IN) :: PUSTRTI  (:) ! X-component of momentum flux [N m^{-2}]
REAL (KIND = JPRB), INTENT(IN) :: PVSTRTI  (:) ! Y-component of momentum flux [N m^{-2}]
REAL (KIND = JPRB), INTENT(IN) :: PLDEPT   (:) ! The lake depth [m]
REAL (KIND = JPRB), INTENT(IN) :: PGEMU    (:) ! The SIN of latitude
REAL (KIND = JPRB), INTENT(IN) :: PFRTI    (:,:)  ! TILE FRACTIONS 
REAL (KIND = JPRB), INTENT(IN) :: PTSTP        ! The model time step [s]
TYPE(TCST),         INTENT(IN) :: YDCST
TYPE(TFLAKE),       INTENT(IN) :: YDFLAKE

! The prognostic variables at previous time step  

REAL (KIND = JPRB), INTENT(IN) :: PTLICEM1M (:) ! Temperature at the air-ice or snow-ice interface [K]
REAL (KIND = JPRB), INTENT(IN) :: PTLMNWM1M (:) ! Mean temperature of the water column [K]
REAL (KIND = JPRB), INTENT(IN) :: PTLWMLM1M (:) ! Mixed-layer temperature [K]
REAL (KIND = JPRB), INTENT(IN) :: PTLBOTM1M (:) ! Temperature at the water-bottom sediment interface [K]
REAL (KIND = JPRB), INTENT(IN) :: PTLSFM1M  (:) ! Shape factor (thermocline)
REAL (KIND = JPRB), INTENT(IN) :: PHLICEM1M (:) ! Ice thickness [m]
REAL (KIND = JPRB), INTENT(IN) :: PHLMLM1M  (:) ! Thickness of the mixed-layer [m]

! The prognostic variables at next time step  

REAL (KIND = JPRD), INTENT(OUT) :: PTLICE   (:) ! Temperature at the air-ice or snow-ice  interface [K]
REAL (KIND = JPRD), INTENT(OUT) :: PTLMNW   (:) ! Mean temperature of the water column [K]
REAL (KIND = JPRD), INTENT(OUT) :: PTLWML   (:) ! Mixed-layer temperature [K]
REAL (KIND = JPRD), INTENT(OUT) :: PTLBOT   (:) ! Temperature at the water-bottom sediment interface [K]
REAL (KIND = JPRD), INTENT(OUT) :: PTLSF    (:) ! Shape factor (thermocline)
REAL (KIND = JPRD), INTENT(OUT) :: PHLICE   (:) ! Ice thickness [m]
REAL (KIND = JPRD), INTENT(OUT) :: PHLML    (:) ! Thickness of the mixed-layer [m]

! Local variables 

TYPE (ROPTICPAR_MEDIUM) :: ZOPTICPAR_WATER ! Optical characteristics of water
TYPE (ROPTICPAR_MEDIUM) :: ZOPTICPAR_ICE   ! Optical characteristics of ice

REAL (KIND = JPRD) :: ZQ_MOMENTUM   (KLON) ! Momentum flux downwards [N m^{-2}]
REAL (KIND = JPRD) :: ZPAR_CORIOLIS (KLON) ! Coriolis parameter [s^{-1}]
REAL (KIND = JPRD) :: ZT_SFC_P      (KLON) ! Surface temperature [K] at the previous time step
REAL (KIND = JPRD) :: ZT_SFC_N      (KLON) ! Surface temperature [K] at the next time step

REAL (KIND = JPRD) :: ZH_ML_P_FLK   (KLON) ! Thickness of the mixed-layer [m] at the previous time step 
REAL (KIND = JPRD) :: ZH_ML_N_FLK   (KLON) ! Thickness of the mixed-layer [m] at the next time step
REAL (KIND = JPRD) :: ZH_ICE_P_FLK  (KLON) ! Ice thickness [m] at the previous time step
REAL (KIND = JPRD) :: ZH_ICE_N_FLK  (KLON) ! Ice thickness [m] at the next time step

REAL (KIND = JPRD) :: ZT_ICE_P_FLK  (KLON) ! Temperature at the snow-ice or air-ice interface [K] at the previous time step
REAL (KIND = JPRD) :: ZT_ICE_N_FLK  (KLON) ! Temperature at the snow-ice or air-ice interface [K] at the next time step
REAL (KIND = JPRD) :: ZT_MNW_P_FLK  (KLON) ! Mean temperature of the water column [K] at the previous time step
REAL (KIND = JPRD) :: ZT_MNW_N_FLK  (KLON) ! Mean temperature of the water column [K] at the next time step
REAL (KIND = JPRD) :: ZT_WML_P_FLK  (KLON) ! Mixed-layer temperature [K] at the previous time step
REAL (KIND = JPRD) :: ZT_WML_N_FLK  (KLON) ! Mixed-layer temperature [K] at the next time step
REAL (KIND = JPRD) :: ZT_BOT_P_FLK  (KLON) ! Temperature at the water-bottom sediment interface [K] at the previous time step
REAL (KIND = JPRD) :: ZT_BOT_N_FLK  (KLON) ! Temperature at the water-bottom sediment interface [K] at the next time step
REAL (KIND = JPRD) :: ZC_T_P_FLK    (KLON) ! Shape factor (thermocline) at the previous time step
REAL (KIND = JPRD) :: ZC_T_N_FLK    (KLON) ! Shape factor (thermocline) at the next time step
                    
REAL (KIND = JPRD) :: ZQ_ICE_FLK    (KLON) ! Heat flux through the snow-ice or air-ice interface [W m^{-2}]
REAL (KIND = JPRD) :: ZQ_W_FLK      (KLON) ! Heat flux through the ice-water or air-water interface [W m^{-2}]
REAL (KIND = JPRD) :: ZQ_BOT_FLK    (KLON) ! Heat flux through the water-bottom sediment interface [W m^{-2}]
REAL (KIND = JPRD) :: ZI_ATM_FLK    (KLON) ! Radiation flux at the lower boundary of the atmosphere [W m^{-2}]
                                           ! i.e. the incident radiation flux with no regard for the surface albedo
REAL (KIND = JPRD) :: ZI_ICE_FLK    (KLON) ! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
REAL (KIND = JPRD) :: ZI_W_FLK      (KLON) ! Radiation flux through the ice-water or air-water interface [W m^{-2}]
REAL (KIND = JPRD) :: ZI_H_FLK      (KLON) ! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
REAL (KIND = JPRD) :: ZI_BOT_FLK    (KLON) ! Radiation flux through the water-bottom sediment interface [W m^{-2}]

REAL (KIND = JPRD) :: ZI_INTM_0_H_FLK (KLON) ! Mean radiation flux over the mixed layer [W m^{-1}]
REAL (KIND = JPRD) :: ZI_INTM_H_D_FLK (KLON) ! Mean radiation flux over the thermocline [W m^{-1}]
REAL (KIND = JPRD) :: ZQ_STAR_FLK     (KLON) ! A generalized heat flux scale [W m^{-2}]
REAL (KIND = JPRD) :: ZU_STAR_W_FLK   (KLON) ! Friction velocity in the surface layer of lake water [m s^{-1}]
REAL (KIND = JPRD) :: ZW_STAR_SFC_FLK (KLON) ! Convective velocity scale, using a generalized heat flux scale [m s^{-1}]

REAL (KIND = JPRD) :: ZPLDEPT (KLON),ZPTSTP ! double precision lake depth 
REAL (KIND = JPRD) :: ZOPTIC1 ! double precision REXTINCOEF_OPTIC
INTEGER (KIND = JPIM) :: JL                  ! Loop index 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FLAKE_DRIVER_MOD:FLAKE_DRIVER',0,ZHOOK_HANDLE)
ASSOCIATE(RDAY=>YDCST%RDAY, RPI=>YDCST%RPI, &
 & LEFLAKE=>YDFLAKE%LEFLAKE, RH_ICE_MIN_FLK=>YDFLAKE%RH_ICE_MIN_FLK, &
 & ROPTICPAR_ICE_OPAQUE=>YDFLAKE%ROPTICPAR_ICE_OPAQUE, &
 & ROPTICPAR_WATER_REF=>YDFLAKE%ROPTICPAR_WATER_REF, RTPL_L_F=>YDFLAKE%RTPL_L_F, &
 & RTPL_RHO_W_R=>YDFLAKE%RTPL_RHO_W_R, RTPL_T_F=>YDFLAKE%RTPL_T_F, &
 & RTPSF_L_EVAP=>YDFLAKE%RTPSF_L_EVAP)

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Set optical characteristics of the lake water, lake ice and snow
!------------------------------------------------------------------------------
! Use default values
  ZOPTICPAR_WATER = ROPTICPAR_WATER_REF
  ZOPTICPAR_ICE   = ROPTICPAR_ICE_OPAQUE   ! Opaque ice

!------------------------------------------------------------------------------
!  Set initial values
!------------------------------------------------------------------------------
  ZPTSTP = PTSTP
  DO JL = KIDIA, KFDIA
    ZPLDEPT(JL) = PLDEPT(JL)
   
    ZT_ICE_P_FLK  (JL) = PTLICEM1M    (JL)
    ZT_MNW_P_FLK  (JL) = PTLMNWM1M    (JL)
    ZT_WML_P_FLK  (JL) = PTLWMLM1M    (JL)
    ZT_BOT_P_FLK  (JL) = PTLBOTM1M    (JL)
    ZC_T_P_FLK    (JL) = PTLSFM1M     (JL)
    ZH_ICE_P_FLK  (JL) = PHLICEM1M    (JL)
    ZH_ML_P_FLK   (JL) = PHLMLM1M     (JL)
    IF (LDLAKEPOINT(JL)) THEN

      ZI_ATM_FLK    (JL) = PSSRFLTI     (JL)
      ZQ_W_FLK      (JL) = PSLRFL       (JL)

      ZQ_MOMENTUM   (JL) = SQRT(PUSTRTI(JL)**2+PVSTRTI(JL)**2)
      ZU_STAR_W_FLK (JL) = SQRT(ZQ_MOMENTUM(JL)/RTPL_RHO_W_R)
  
      ZPAR_CORIOLIS (JL) = 1.E-4_JPRD !2_JPRD*(2_JPRD*RPI/RDAY)*PGEMU(JL)

!------------------------------------------------------------------------------
!  Compute heat fluxes Q_snow_flk, Q_ice_flk, Q_w_flk
!------------------------------------------------------------------------------

      ZQ_W_FLK(JL) = ZQ_W_FLK(JL)  &
     & + PAHFSTI(JL) + RTPSF_L_EVAP*PEVAPTI(JL)    ! Add sensible and latent heat fluxes (notice the signs)
      IF(ZH_ICE_P_FLK(JL).GE.RH_ICE_MIN_FLK) THEN  ! Ice exists
        ZQ_ICE_FLK (JL) = ZQ_W_FLK(JL) + RTPL_L_F*PEVAPTI(JL)  
        ZQ_W_FLK   (JL) = 0._JPRD
        ZT_SFC_P   (JL) = PTLICEM1M (JL)
      ELSE                                         ! No ice cover
        ZQ_ICE_FLK (JL) = 0._JPRD
        ZT_SFC_P   (JL) = PTLWMLM1M (JL)
      END IF
    ENDIF
  ENDDO

! Debugging test
!  ZQ_W_FLK = 0._JPRD
!  IF (ZH_ICE_P_FLK == 0._JPRD) THEN
!   ZQ_W_FLK = ZQ_W_FLK - 2.E+3_JPRD
!  ELSE
!   ZQ_W_FLK = ZQ_W_FLK - 2.E+2_JPRD
!  ENDIF
!  ZI_ATM_FLK = 50._JPRD
!  ZQ_W_FLK = 0._JPRD
!  ZI_ATM_FLK = 0._JPRD
!------------------------------------------------------------------------------
!  Compute solar radiation fluxes (positive downward)
!------------------------------------------------------------------------------

  CALL FLAKERAD                                            &
    & (KIDIA        , KFDIA           , LDLAKEPOINT     ,  &
    &  ZPLDEPT       , ZOPTICPAR_WATER , ZOPTICPAR_ICE   ,  &
    &  YDFLAKE      , &
    &  ZI_ATM_FLK   , ZH_ICE_P_FLK    , ZH_ML_P_FLK     ,  &

    &  ZI_ICE_FLK   , ZI_BOT_FLK      , ZI_W_FLK        ,  &
    &  ZI_H_FLK     , ZI_INTM_0_H_FLK , ZI_INTM_H_D_FLK    )
    
    
!------------------------------------------------------------------------------
!  Advance FLake variables
!------------------------------------------------------------------------------
  ZOPTIC1 = ZOPTICPAR_WATER%REXTINCOEF_OPTIC(1)
  CALL FLAKEENE                                                                      &
    &  (KIDIA              , KFDIA                  , KLON                    ,      &
    &   ZPLDEPT             , ZPAR_CORIOLIS          , LDLAKEPOINT             ,      &
    &   ZOPTIC1         , YDFLAKE                 ,      &
    &   ZPTSTP              , ZT_SFC_P               , ZT_ICE_P_FLK            ,      &
    &   ZT_WML_P_FLK       , ZT_MNW_P_FLK           , ZT_BOT_P_FLK            ,      &
    &   ZH_ICE_P_FLK       , ZH_ML_P_FLK            , ZC_T_P_FLK              ,      &
    &   ZQ_W_FLK           , ZQ_ICE_FLK             , ZU_STAR_W_FLK           ,      &
    &   ZI_ICE_FLK         , ZI_BOT_FLK             , ZI_W_FLK                ,      &
    &   ZI_H_FLK           , ZI_INTM_0_H_FLK        , ZI_INTM_H_D_FLK         ,      &
  
    &   ZT_ICE_N_FLK       , ZT_WML_N_FLK           , ZT_MNW_N_FLK            ,      &
    &   ZT_BOT_N_FLK       , ZH_ICE_N_FLK           , ZH_ML_N_FLK             ,      &
    &   ZC_T_N_FLK         , ZT_SFC_N                                                )

!------------------------------------------------------------------------------
!  Set output values
!------------------------------------------------------------------------------

  PTLICE  = ZT_ICE_N_FLK 
  PTLMNW  = ZT_MNW_N_FLK 
  PTLWML  = ZT_WML_N_FLK 
  PTLBOT  = ZT_BOT_N_FLK 
  PTLSF   = ZC_T_N_FLK   
  PHLICE  = ZH_ICE_N_FLK 
  PHLML   = ZH_ML_N_FLK  


!------------------------------------------------------------------------------
!  End calculations
!==============================================================================
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FLAKE_DRIVER_MOD:FLAKE_DRIVER',1,ZHOOK_HANDLE)

END SUBROUTINE FLAKE_DRIVER
END MODULE FLAKE_DRIVER_MOD
