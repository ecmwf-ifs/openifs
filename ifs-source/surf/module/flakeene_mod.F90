MODULE FLAKEENE_MOD
CONTAINS
SUBROUTINE FLAKEENE                                                      &
  &  (KIDIA           , KFDIA             , KLON                  ,      &
  &   PDEPTH_W        , PAR_CORIOLIS      , LDLAKEPOINT           ,      &
  &   PEXTINCOEF_WATER_TYP                , YDFLAKE               ,     &
  &   PDEL_TIME       , PT_SFC_P          , PT_ICE_P_FLK          ,      &
  &   PT_WML_P_FLK    , PT_MNW_P_FLK      , PT_BOT_P_FLK          ,      &
  &   PH_ICE_P_FLK    , PH_ML_P_FLK       , PC_T_P_FLK            ,      &
  &   PQ_W_FLK        , PQ_ICE_FLK        , PU_STAR_W_FLK         ,      &
  &   PI_ICE_FLK      , PI_BOT_FLK        , PI_W_FLK              ,      &
  &   PI_H_FLK        , PI_INTM_0_H_FLK   , PI_INTM_H_D_FLK       ,      &
  &   PT_ICE_N_FLK    , PT_WML_N_FLK      , PT_MNW_N_FLK          ,      &
  &   PT_BOT_N_FLK    , PH_ICE_N_FLK      , PH_ML_N_FLK           ,      &
  &   PC_T_N_FLK      , PT_SFC_N                                         )         

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
!  The main driving routine of the lake model FLake 
!  where computations are performed.
!  Advances the surface temperature
!  and other FLake variables one time step.
!  At the moment, the Euler explicit scheme is used.
!
!  Lines embraced/marked with "!_dm" are DM's comments
!  that may be helpful to a user.
!
!
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
! 1.01T      10.03.2008     V. M. Stepanenko
!  The code and variables, relevant to snow and bottom sediments are omitted
! 1.02       11.06.2010     R. Salgado and G. Balsamo 
!  Set lake shape factor to constant 0.65 value for safety (instability issues)
! 1.03       11.11.2010     G. Balsamo 
!  Fixes for coupled atmospheric runs
!            17.12.2015     F. Vana
!  Support for single precision

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
USE YOS_FLAKE, ONLY : TFLAKE

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Input (procedure arguments)

INTEGER (KIND = JPIM), INTENT(IN):: KIDIA  
INTEGER (KIND = JPIM), INTENT(IN):: KFDIA  
INTEGER (KIND = JPIM), INTENT(IN):: KLON   

LOGICAL, INTENT(IN):: LDLAKEPOINT (:)                     ! Indicates, if it is lake point (.TRUE.) 
  
REAL (KIND = JPRD), INTENT(IN):: PDEPTH_W             (:) ! The lake depth [m]
REAL (KIND = JPRD), INTENT(IN):: PAR_CORIOLIS         (:) ! The Coriolis parameter [s^{-1}]
REAL (KIND = JPRD), INTENT(IN):: PT_SFC_P             (:) ! Surface temperature at the previous time step [K]  
                                                          ! (equal to either T_ice, T_snow or to T_wML)

REAL (KIND = JPRD), INTENT(IN):: PEXTINCOEF_WATER_TYP     ! "Typical" extinction coeff.of the lake water [m^{-1}]
                                                          ! used to compute the equilibrium CBL depth
TYPE(TFLAKE),       INTENT(IN):: YDFLAKE

REAL (KIND = JPRD), INTENT(IN):: PDEL_TIME                ! The model time step [s]


REAL (KIND = JPRD), INTENT(IN):: PT_ICE_P_FLK    (:) 
REAL (KIND = JPRD), INTENT(IN):: PT_WML_P_FLK    (:)
REAL (KIND = JPRD), INTENT(IN):: PT_MNW_P_FLK    (:)
REAL (KIND = JPRD), INTENT(IN):: PT_BOT_P_FLK    (:)
REAL (KIND = JPRD), INTENT(IN):: PH_ICE_P_FLK    (:)
REAL (KIND = JPRD), INTENT(IN):: PH_ML_P_FLK     (:)
REAL (KIND = JPRD), INTENT(IN):: PC_T_P_FLK      (:)
REAL (KIND = JPRD), INTENT(IN):: PQ_ICE_FLK      (:)
REAL (KIND = JPRD), INTENT(IN):: PU_STAR_W_FLK   (:)
REAL (KIND = JPRD), INTENT(IN):: PI_ICE_FLK      (:)
REAL (KIND = JPRD), INTENT(IN):: PI_BOT_FLK      (:)
REAL (KIND = JPRD), INTENT(IN):: PI_W_FLK        (:)
REAL (KIND = JPRD), INTENT(IN):: PI_H_FLK        (:)
REAL (KIND = JPRD), INTENT(IN):: PI_INTM_0_H_FLK (:)
REAL (KIND = JPRD), INTENT(IN):: PI_INTM_H_D_FLK (:)


REAL (KIND = JPRD), INTENT(INOUT):: PQ_W_FLK (:)

!  Output (procedure arguments)
REAL (KIND = JPRD), INTENT(OUT):: PT_SFC_N  (:)         ! Updated surface temperature [K] 
                                                        ! (equal to the updated value of either T_ice, T_snow or T_wML)

REAL (KIND = JPRD), INTENT(OUT) :: PT_ICE_N_FLK (:)
REAL (KIND = JPRD), INTENT(OUT) :: PT_WML_N_FLK (:)
REAL (KIND = JPRD), INTENT(OUT) :: PT_MNW_N_FLK (:)
REAL (KIND = JPRD), INTENT(OUT) :: PT_BOT_N_FLK (:)
REAL (KIND = JPRD), INTENT(OUT) :: PH_ICE_N_FLK (:)
REAL (KIND = JPRD), INTENT(OUT) :: PH_ML_N_FLK  (:)
REAL (KIND = JPRD), INTENT(OUT) :: PC_T_N_FLK   (:)


!  Local variables of type LOGICAL
LOGICAL:: LL_ICE_CREATE     ! Switch, .TRUE. = ice does not exist but should be created
LOGICAL:: LL_ICE_MELTABOVE  ! Switch, .TRUE. = snow/ice melting from above takes place

!  Local variables of type INTEGER
INTEGER (KIND = JPIM):: I, JL  ! Loop indexes

!  Local variables of type REAL
REAL (KIND = JPRD):: ZD_T_MNW_DT      ! Time derivative of T_mnw [K s^{-1}] 
REAL (KIND = JPRD):: ZD_T_ICE_DT      ! Time derivative of T_ice [K s^{-1}] 
REAL (KIND = JPRD):: ZD_T_BOT_DT      ! Time derivative of T_bot [K s^{-1}] 
REAL (KIND = JPRD):: ZD_H_ICE_DT      ! Time derivative of h_ice [m s^{-1}]
REAL (KIND = JPRD):: ZD_H_ML_DT       ! Time derivative of h_ML [m s^{-1}]
REAL (KIND = JPRD):: ZD_C_T_DT        ! Time derivative of C_T [s^{-1}]

!  Local variables of type REAL
REAL (KIND = JPRD):: ZN_T_MEAN            ! The mean buoyancy frequency in the thermocline [s^{-1}] 
REAL (KIND = JPRD):: ZM_H_SCALE           ! The ZM96 equilibrium SBL depth scale [m] 
REAL (KIND = JPRD):: ZCONV_EQUIL_H_SCALE  ! The equilibrium CBL depth scale [m]

!  Local variables of type REAL

REAL (KIND = JPRD):: ZH_ICE_THRESHOLD  ! If h_ice<ZH_ICE_THRESHOLD, use quasi-equilibrium ice model 
REAL (KIND = JPRD):: ZFLK_STR_1        ! Help storage variable
REAL (KIND = JPRD):: ZFLK_STR_2        ! Help storage variable
REAL (KIND = JPRD):: ZR_H_ICESNOW      ! Dimensionless ratio, used to store intermediate results
REAL (KIND = JPRD):: ZR_RHO_C_ICESNOW  ! Dimensionless ratio, used to store intermediate results
REAL (KIND = JPRD):: ZR_TI_ICESNOW     ! Dimensionless ratio, used to store intermediate results
REAL (KIND = JPRD):: ZR_TSTAR_ICESNOW  ! Dimensionless ratio, used to store intermediate results

!  The shape factor(s) at the previous time step ("p") and the updated value(s) ("n") 
REAL (KIND = JPRD):: ZC_TT_FLK       ! Dimensionless parameter (thermocline)
REAL (KIND = JPRD):: ZC_Q_FLK        ! Shape factor with respect to the heat flux (thermocline)
REAL (KIND = JPRD):: ZC_I_FLK        ! Shape factor (ice)

!  Derivatives of the shape functions
REAL (KIND = JPRD):: ZPHI_T_PR0_FLK  ! d\Phi_T(0)/d\zeta   (thermocline)
REAL (KIND = JPRD):: ZPHI_I_PR0_FLK  ! d\Phi_I(0)/d\zeta_I (ice)
REAL (KIND = JPRD):: ZPHI_I_PR1_FLK  ! d\Phi_I(1)/d\zeta_I (ice)

!  Heat and radiation fluxes
REAL (KIND = JPRD):: ZQ_BOT_FLK      ! Heat flux through the water-bottom sediment interface [W m^{-2}]
REAL (KIND = JPRD):: ZQ_STAR_FLK     ! A generalized heat flux scale [W m^{-2}]

!  Velocity scales
REAL (KIND = JPRD):: ZW_STAR_SFC_FLK ! Convective velocity scale, 
                                     ! using a generalized heat flux scale [m s^{-1}]

REAL (KIND = JPRD):: ZFLAKE_BUOYPAR  ! Buoyancy parameter
REAL (KIND = JPRD):: ZT_WATER        ! Temperature of water
REAL (KIND = JPRD):: ZDEPTH_W(KLON)  ! Depth really used in the calculation
REAL (KIND = JPRD):: ZDEPTH_W_MAX    ! Maximum lake depth allowed
REAL (KIND = JPRD):: ZDEPTH_W_MIN    ! Minimum lake depth allowed
REAL (KIND = JPRD):: ZC_T_CON        ! Mean lake shape factor
REAL (KIND = JPRD):: ZD_C_T_DT_MAX   ! Maximum abs lake shape factor tendency
REAL (KIND = JPRB):: ZEPS            ! Small "epsilon" in jprb precision

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! Declaration of FLAKE_BOUYPAR function (calculates bouyancy parameter)
ZFLAKE_BUOYPAR(ZT_WATER) = YDFLAKE%RTPL_GRAV*YDFLAKE%RTPL_A_T*(ZT_WATER-YDFLAKE%RTPL_T_R)

IF (LHOOK) CALL DR_HOOK('FLAKEENE_MOD:FLAKEENE',0,ZHOOK_HANDLE)
ASSOCIATE(RC_CBL_1=>YDFLAKE%RC_CBL_1, RC_CBL_2=>YDFLAKE%RC_CBL_2, &
 & RC_I_LIN=>YDFLAKE%RC_I_LIN, RC_I_MR=>YDFLAKE%RC_I_MR, &
 & RC_RELAX_C=>YDFLAKE%RC_RELAX_C, RC_RELAX_H=>YDFLAKE%RC_RELAX_H, &
 & RC_SBL_ZM_I=>YDFLAKE%RC_SBL_ZM_I, RC_SBL_ZM_N=>YDFLAKE%RC_SBL_ZM_N, &
 & RC_SBL_ZM_S=>YDFLAKE%RC_SBL_ZM_S, RC_SMALL_FLK=>YDFLAKE%RC_SMALL_FLK, &
 & RC_TT_1=>YDFLAKE%RC_TT_1, RC_TT_2=>YDFLAKE%RC_TT_2, &
 & RC_T_MAX=>YDFLAKE%RC_T_MAX, RC_T_MIN=>YDFLAKE%RC_T_MIN, &
 & RH_ICE_MAX=>YDFLAKE%RH_ICE_MAX, RH_ICE_MIN_FLK=>YDFLAKE%RH_ICE_MIN_FLK, &
 & RH_ML_MAX_FLK=>YDFLAKE%RH_ML_MAX_FLK, RH_ML_MIN_FLK=>YDFLAKE%RH_ML_MIN_FLK, &
 & RPHI_I_AST_MR=>YDFLAKE%RPHI_I_AST_MR, RPHI_I_PR0_LIN=>YDFLAKE%RPHI_I_PR0_LIN, &
 & RPHI_I_PR1_LIN=>YDFLAKE%RPHI_I_PR1_LIN, RPHI_T_PR0_1=>YDFLAKE%RPHI_T_PR0_1, &
 & RPHI_T_PR0_2=>YDFLAKE%RPHI_T_PR0_2, RTPL_A_T=>YDFLAKE%RTPL_A_T, &
 & RTPL_C_I=>YDFLAKE%RTPL_C_I, RTPL_C_W=>YDFLAKE%RTPL_C_W, &
 & RTPL_GRAV=>YDFLAKE%RTPL_GRAV, RTPL_KAPPA_I=>YDFLAKE%RTPL_KAPPA_I, &
 & RTPL_KAPPA_W=>YDFLAKE%RTPL_KAPPA_W, RTPL_L_F=>YDFLAKE%RTPL_L_F, &
 & RTPL_RHO_I=>YDFLAKE%RTPL_RHO_I, RTPL_RHO_W_R=>YDFLAKE%RTPL_RHO_W_R, &
 & RTPL_T_F=>YDFLAKE%RTPL_T_F, RTPL_T_R=>YDFLAKE%RTPL_T_R, &
 & RU_STAR_MIN_FLK=>YDFLAKE%RU_STAR_MIN_FLK)

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------
ZEPS=10._JPRB*EPSILON(ZEPS)           ! Small epsilon
ZDEPTH_W_MAX = 50.0_JPRD              ! Maximum lake depth simulated by FLAKE
ZDEPTH_W_MIN =  2.0_JPRD              ! Minimum lake depth simulated by FLAKE
ZC_T_CON=(RC_T_MAX+RC_T_MIN)/2._JPRD  ! Average lake shape factor
ZD_C_T_DT_MAX=(RC_T_MAX-RC_T_MIN)/(86400._JPRD*30._JPRD) ! Maximum abs C_T tendency

!_dm 
! Security. Set time-rate-of-change of prognostic variables to zero.
! Set prognostic variables to their values at the previous time step.
! (This is to avoid spurious changes of prognostic variables 
! when FLake is used within a 3D model, e.g. to avoid spurious generation of ice 
! at the neighbouring lake points as noticed by Burkhardt Rockel.)
!_dm 

PT_ICE_N_FLK  = PT_ICE_P_FLK    
PT_WML_N_FLK  = PT_WML_P_FLK   
PT_MNW_N_FLK  = PT_MNW_P_FLK     
PT_BOT_N_FLK  = PT_BOT_P_FLK  
PH_ICE_N_FLK  = PH_ICE_P_FLK   
PH_ML_N_FLK   = PH_ML_P_FLK    
PC_T_N_FLK    = PC_T_P_FLK    

MAIN_LOOP: DO JL = KIDIA, KFDIA

LAKEPOINT: IF (LDLAKEPOINT(JL)) THEN ! Prognostic variables are updated only
                                     ! for points, that are indicated as lakes

! Limit lake depth to MAX MIN values (note: FLAKE is a shallow-lake model)
ZDEPTH_W(JL)=MIN(ZDEPTH_W_MAX ,MAX(ZDEPTH_W_MIN,PDEPTH_W(JL)))
! Security for thickness (prevents ICs imbalance)
PH_ML_N_FLK(JL)=MIN(PH_ML_N_FLK(JL),ZDEPTH_W(JL))
PH_ICE_N_FLK(JL)=MIN(PH_ICE_N_FLK(JL),ZDEPTH_W(JL))
! Security for temperature (prevents ICs imbalance)
PT_WML_N_FLK(JL) = MAX(PT_WML_N_FLK(JL),RTPL_T_F)  
PT_MNW_N_FLK(JL) = MAX(PT_MNW_N_FLK(JL),RTPL_T_F)  
PT_BOT_N_FLK(JL) = MAX(PT_BOT_N_FLK(JL),RTPL_T_F)  
PT_ICE_N_FLK(JL) = MIN(PT_ICE_N_FLK(JL),RTPL_T_F) 

ZD_T_MNW_DT   = 0._JPRD 
ZD_T_ICE_DT   = 0._JPRD 
ZD_T_BOT_DT   = 0._JPRD 
ZD_H_ICE_DT   = 0._JPRD 
ZD_H_ML_DT    = 0._JPRD 
ZD_C_T_DT     = 0._JPRD 

!------------------------------------------------------------------------------
!  Compute fluxes, using variables from the previous time step.
!------------------------------------------------------------------------------

!_dm
! At this point, the heat and radiation fluxes, namely,
! PQ_ICE_FLK, PQ_W_FLK,PI_ICE_FLK, PI_W_FLK, PI_H_FLK, PI_BOT_FLK,     
! the mean radiation flux over the mixed layer, PI_INTM_0_H_FLK, 
! and the mean radiation flux over the thermocline, PI_INTM_H_D_FLK, 
! should be known.
! They are computed within routine FLAKE_RADFLUX, called by FLAKE_DRIVER.
! In case a lake is ice-covered, PQ_W_FLK is re-computed below.
!_dm

! Heat flux through the ice-water interface
IF(PH_ICE_P_FLK(JL).GE.RH_ICE_MIN_FLK) THEN    ! Ice exists 
  IF(PH_ML_P_FLK(JL).LE.RH_ML_MIN_FLK) THEN    ! Mixed-layer depth is zero, compute flux 
    PQ_W_FLK(JL) = -RTPL_KAPPA_W*(PT_BOT_P_FLK(JL)-PT_WML_P_FLK(JL))/ &
   & ZDEPTH_W(JL)                             ! Flux with linear T(z) 
    ZPHI_T_PR0_FLK = RPHI_T_PR0_1*PC_T_P_FLK(JL)-RPHI_T_PR0_2         ! d\Phi(0)/d\zeta (thermocline)
    PQ_W_FLK(JL)   = PQ_W_FLK(JL)*MAX(ZPHI_T_PR0_FLK, 1._JPRD)      ! Account for an increased d\Phi(0)/d\zeta 
  ELSE                    
    PQ_W_FLK(JL)   = 0._JPRD                  ! Mixed-layer depth is greater than zero, set flux to zero
  END IF   
END IF   

! A generalized heat flux scale 
ZQ_STAR_FLK = PQ_W_FLK(JL) + PI_W_FLK(JL) + PI_H_FLK(JL) - 2._JPRD*PI_INTM_0_H_FLK(JL)

! Heat flux through the water-bottom sediment interface
ZQ_BOT_FLK = 0._JPRD   ! The bottom-sediment scheme is not used



!------------------------------------------------------------------------------
!  Check if ice exists or should be created.
!  If so, compute the thickness and the temperature of ice.
!------------------------------------------------------------------------------

!_dm
! Notice that a quasi-equilibrium ice model is used 
! to avoid numerical instability when the ice is thin.
! This is always the case when new ice is created.
!_dm


! Default values
LL_ICE_CREATE    = .FALSE.  
LL_ICE_MELTABOVE = .FALSE.  

ICE_EXIST: IF(PH_ICE_P_FLK(JL).LT.RH_ICE_MIN_FLK) THEN   ! Ice does not exist 

  LL_ICE_CREATE = PT_WML_P_FLK(JL).LE.(RTPL_T_F+RC_SMALL_FLK).AND. &
 & PQ_W_FLK(JL).LT.0._JPRD
  IF(LL_ICE_CREATE) THEN                            ! Ice does not exist but should be created
    ZD_H_ICE_DT = -PQ_W_FLK(JL)/RTPL_RHO_I/RTPL_L_F
    PH_ICE_N_FLK(JL) = PH_ICE_P_FLK(JL) + ZD_H_ICE_DT*PDEL_TIME   ! Advance h_ice 
    PT_ICE_N_FLK(JL) = RTPL_T_F + PH_ICE_N_FLK(JL)*PQ_W_FLK(JL)/ &
   & RTPL_KAPPA_I/RPHI_I_PR0_LIN                                    ! Ice temperature
    ZPHI_I_PR1_FLK = RPHI_I_PR1_LIN                              & 
                   + RPHI_I_AST_MR*MIN(1._JPRD, PH_ICE_N_FLK(JL)/RH_ICE_MAX) ! d\Phi_I(1)/d\zeta_I (ice)
  END IF

ELSE ICE_EXIST                                     ! Ice exists

  MELTING: IF(PT_ICE_P_FLK(JL).GE.(RTPL_T_F-RC_SMALL_FLK)) THEN    ! T_sfc = T_f, check for melting from above
      ZFLK_STR_1 = PQ_ICE_FLK(JL) + PI_ICE_FLK(JL) - PI_W_FLK(JL) - PQ_W_FLK(JL)  ! Atmospheric forcing + heating from the water
    IF(ZFLK_STR_1.GE.0._JPRD) THEN  ! Melting of ice from above, snow accumulation may occur
      LL_ICE_MELTABOVE = .TRUE.
      ZD_H_ICE_DT  = -ZFLK_STR_1/RTPL_L_F/RTPL_RHO_I 
    END IF 
    IF(LL_ICE_MELTABOVE) THEN  ! Melting from above takes place
      PH_ICE_N_FLK(JL)  = PH_ICE_P_FLK(JL)  + ZD_H_ICE_DT *PDEL_TIME  ! Advance h_ice
      PT_ICE_N_FLK(JL)  = RTPL_T_F                                     ! Set T_ice to the freezing point
    END IF

  END IF MELTING

  NO_MELTING: IF(.NOT.LL_ICE_MELTABOVE) THEN                 ! No melting from above
    
    ZPHI_I_PR0_FLK = PH_ICE_P_FLK(JL)/RH_ICE_MAX                          ! h_ice relative to its maximum value
    ZC_I_FLK = RC_I_LIN - RC_I_MR*(1._JPRD+RPHI_I_AST_MR)*ZPHI_I_PR0_FLK    ! Shape factor (ice)
    ZPHI_I_PR1_FLK = RPHI_I_PR1_LIN + RPHI_I_AST_MR*ZPHI_I_PR0_FLK         ! d\Phi_I(1)/d\zeta_I (ice)
    ZPHI_I_PR0_FLK = RPHI_I_PR0_LIN - ZPHI_I_PR0_FLK                      ! d\Phi_I(0)/d\zeta_I (ice)

    ZH_ICE_THRESHOLD = MAX(1._JPRD, 2._JPRD*ZC_I_FLK*RTPL_C_I*(RTPL_T_F-PT_ICE_P_FLK(JL))/ &
   & RTPL_L_F)
    ZH_ICE_THRESHOLD = ZPHI_I_PR0_FLK/ZC_I_FLK*RTPL_KAPPA_I/RTPL_RHO_I/RTPL_C_I* &
   & ZH_ICE_THRESHOLD
    ZH_ICE_THRESHOLD = SQRT(MAX(0.0_JPRD,ZH_ICE_THRESHOLD*PDEL_TIME))     ! Threshold value of h_ice
    ZH_ICE_THRESHOLD = MIN(0.9_JPRD*RH_ICE_MAX, MAX(ZH_ICE_THRESHOLD, RH_ICE_MIN_FLK))
                                                                          ! h_ice(threshold) < 0.9*RH_Ice_max

    IF(PH_ICE_P_FLK(JL).LT.ZH_ICE_THRESHOLD) THEN  ! Use a quasi-equilibrium ice model

      ZFLK_STR_1 = PQ_ICE_FLK(JL) + PI_ICE_FLK(JL) - PI_W_FLK(JL)
      ZD_H_ICE_DT = -(ZFLK_STR_1-PQ_W_FLK(JL))/RTPL_L_F/RTPL_RHO_I
      PH_ICE_N_FLK(JL) = PH_ICE_P_FLK(JL) + ZD_H_ICE_DT *PDEL_TIME       ! Advance h_ice
      PT_ICE_N_FLK(JL) = RTPL_T_F + PH_ICE_N_FLK(JL)*ZFLK_STR_1/RTPL_KAPPA_I/ &
     & ZPHI_I_PR0_FLK  ! Ice temperature

    ELSE                                     ! Use a complete ice model

      ZD_H_ICE_DT  = RTPL_KAPPA_I*(RTPL_T_F-PT_ICE_P_FLK(JL))/PH_ICE_P_FLK(JL)* &
     & ZPHI_I_PR0_FLK
      ZD_H_ICE_DT  = (PQ_W_FLK(JL)+ZD_H_ICE_DT)/RTPL_L_F/RTPL_RHO_I
      PH_ICE_N_FLK(JL) = PH_ICE_P_FLK(JL)  + ZD_H_ICE_DT*PDEL_TIME  ! Advance h_ice

      ZR_TI_ICESNOW = RTPL_C_I*(RTPL_T_F-PT_ICE_P_FLK(JL))/RTPL_L_F   ! Dimensionless parameter
      ZR_TSTAR_ICESNOW = 1._JPRD - ZC_I_FLK                        ! Dimensionless parameter
      ZR_TSTAR_ICESNOW = ZR_TSTAR_ICESNOW*ZR_TI_ICESNOW            ! Dimensionless parameter

      ZFLK_STR_2 = PQ_ICE_FLK(JL)+PI_ICE_FLK(JL)-PI_W_FLK(JL)      ! Atmospheric fluxes
      ZFLK_STR_1  = ZC_I_FLK*PH_ICE_P_FLK(JL)
      ZD_T_ICE_DT = 0._JPRD
      ZD_T_ICE_DT = ZD_T_ICE_DT + RTPL_KAPPA_I*(RTPL_T_F-PT_ICE_P_FLK(JL))/ &
     & PH_ICE_P_FLK(JL)*ZPHI_I_PR0_FLK * (1._JPRD-ZR_TSTAR_ICESNOW)  ! Add flux due to heat conduction
      ZD_T_ICE_DT = ZD_T_ICE_DT - ZR_TSTAR_ICESNOW*PQ_W_FLK(JL)      ! Add flux from water to ice
      ZD_T_ICE_DT = ZD_T_ICE_DT + ZFLK_STR_2                         ! Add atmospheric fluxes
      ZD_T_ICE_DT = ZD_T_ICE_DT/RTPL_RHO_I/RTPL_C_I                    ! Total forcing
      ZD_T_ICE_DT = ZD_T_ICE_DT/ZFLK_STR_1                           ! dT_ice/dt 
      PT_ICE_N_FLK(JL) = PT_ICE_P_FLK(JL) + ZD_T_ICE_DT*PDEL_TIME    ! Advance T_ice
    END IF

    ZPHI_I_PR1_FLK = MIN(1._JPRD, PH_ICE_N_FLK(JL)/RH_ICE_MAX)        ! h_ice relative to its maximum value
    ZPHI_I_PR1_FLK = RPHI_I_PR1_LIN + RPHI_I_AST_MR*ZPHI_I_PR1_FLK     ! d\Phi_I(1)/d\zeta_I (ice)

  END IF NO_MELTING

END IF ICE_EXIST   

! Security, limit h_ice by its maximum value
PH_ICE_N_FLK(JL) = MIN(PH_ICE_N_FLK(JL), RH_ICE_MAX)      

! Security, limit the ice and snow temperatures by the freezing point 
PT_ICE_N_FLK(JL) = MIN(PT_ICE_N_FLK(JL),  RTPL_T_F)    

! Security, avoid too low values (these constraints are used for debugging purposes)
PT_ICE_N_FLK(JL) = MAX(PT_ICE_N_FLK(JL),  73.15_JPRD)    

! Remove too thin ice and/or snow
IF(PH_ICE_N_FLK(JL).LT.RH_ICE_MIN_FLK)  THEN        ! Check ice
  PH_ICE_N_FLK(JL) = 0._JPRD       ! Ice is too thin, remove it, and
  PT_ICE_N_FLK(JL) = RTPL_T_F       ! set T_ice to the freezing point.
  LL_ICE_CREATE    = .FALSE.       ! "Exotic" case, ice has been created but proved to be too thin
END IF


!------------------------------------------------------------------------------
!  Compute the mean temperature of the water column.
!------------------------------------------------------------------------------

IF(LL_ICE_CREATE) PQ_W_FLK(JL) = 0._JPRD     ! Ice has just been created, set Q_w to zero
ZD_T_MNW_DT = (PQ_W_FLK(JL) - ZQ_BOT_FLK + PI_W_FLK(JL) - PI_BOT_FLK(JL))/ &
& RTPL_RHO_W_R/RTPL_C_W/ZDEPTH_W(JL)
PT_MNW_N_FLK(JL) = PT_MNW_P_FLK(JL) + ZD_T_MNW_DT*PDEL_TIME   ! Advance T_mnw
PT_MNW_N_FLK(JL) = MAX(PT_MNW_N_FLK(JL), RTPL_T_F)             ! Limit T_mnw by the freezing point 


!------------------------------------------------------------------------------
!  Compute the mixed-layer depth, the mixed-layer temperature, 
!  the bottom temperature and the shape factor
!  with respect to the temperature profile in the thermocline. 
!  Different formulations are used, depending on the regime of mixing. 
!------------------------------------------------------------------------------

HTC_WATER: IF(PH_ICE_N_FLK(JL).GE.RH_ICE_MIN_FLK) THEN    ! Ice exists

  PT_MNW_N_FLK(JL) = MIN(PT_MNW_N_FLK(JL), RTPL_T_R) ! Limit the mean temperature under the ice by T_r 
  PT_WML_N_FLK(JL) = RTPL_T_F                        ! The mixed-layer temperature is equal to the freezing point 

  IF(LL_ICE_CREATE) THEN                       ! Ice has just been created 
    IF(PH_ML_P_FLK(JL).GE.ZDEPTH_W(JL)-RH_ML_MIN_FLK) THEN    ! h_ML=D when ice is created 
      PH_ML_N_FLK(JL) = 0._JPRD                ! Set h_ML to zero 
      PC_T_N_FLK(JL)  = ZC_T_CON               ! Set C_T to constant value
    ELSE                                       ! h_ML<D when ice is created 
      PH_ML_N_FLK(JL) = PH_ML_P_FLK(JL)        ! h_ML remains unchanged 
      PC_T_N_FLK(JL)  = PC_T_P_FLK(JL)         ! C_T (thermocline) remains unchanged 
    END IF 
    !*GA: change condition .NE. to a more strict one
    !*IF(PH_ML_N_FLK(JL) .NE. ZDEPTH_W(JL)) THEN !avoid 1 singularity
    IF(ABS(PH_ML_N_FLK(JL)-ZDEPTH_W(JL)) > ZEPS) THEN
      PT_BOT_N_FLK(JL) = PT_WML_N_FLK(JL) - &
   & (PT_WML_N_FLK(JL)-PT_MNW_N_FLK(JL))/PC_T_N_FLK(JL)/ &
   & (1._JPRD-PH_ML_N_FLK(JL)/ZDEPTH_W(JL))    ! Update the bottom temperature 
    ENDIF
  ELSE IF(PT_BOT_P_FLK(JL).LT.RTPL_T_R) THEN   ! Ice exists and T_bot < T_r, molecular heat transfer 
    PH_ML_N_FLK(JL) = PH_ML_P_FLK(JL)          ! h_ML remains unchanged 
    PC_T_N_FLK(JL) = PC_T_P_FLK(JL)            ! C_T (thermocline) remains unchanged 
    !*GA: change condition .NE. to a more strict one
    !*IF(PH_ML_N_FLK(JL) .NE. ZDEPTH_W(JL)) THEN !avoid 1 singularity
    IF(ABS(PH_ML_N_FLK(JL)-ZDEPTH_W(JL)) > ZEPS) THEN
      PT_BOT_N_FLK(JL) = PT_WML_N_FLK(JL) - &
   & (PT_WML_N_FLK(JL)-PT_MNW_N_FLK(JL))/PC_T_N_FLK(JL)/ &
   & (1._JPRD-PH_ML_N_FLK(JL)/ZDEPTH_W(JL))    ! Update the bottom temperature 
    ENDIF
  ELSE                                         ! Ice exists and T_bot = T_r, convection due to bottom heating 
    PT_BOT_N_FLK(JL) = RTPL_T_R                ! T_bot is equal to the temperature of maximum density 
    IF(PH_ML_P_FLK(JL).GE.RC_SMALL_FLK) THEN   ! h_ML > 0 
      PC_T_N_FLK(JL) = PC_T_P_FLK(JL)          ! C_T (thermocline) remains unchanged 
      PH_ML_N_FLK(JL) = ZDEPTH_W(JL)*(1._JPRD-(PT_WML_N_FLK(JL)-PT_MNW_N_FLK(JL))/ &
     & (PT_WML_N_FLK(JL)-PT_BOT_N_FLK(JL))/PC_T_N_FLK(JL))
      PH_ML_N_FLK(JL) = MAX(PH_ML_N_FLK(JL), 0._JPRD) ! Update the mixed-layer depth  
      ! 2019 fix from Mironov
      ! Security, PH_ML_N_FLK  is equal or very close to the lake depth
      IF(PH_ML_N_FLK(JL) >= (ZDEPTH_W(JL)-RH_ML_MIN_FLK)) THEN
          PH_ML_N_FLK(JL) = 0._JPRD                 ! Set mixed-layer depth to zero
          PT_BOT_N_FLK(JL) = PT_WML_N_FLK(JL) - &
                           & (PT_WML_N_FLK(JL)-PT_MNW_N_FLK(JL))/PC_T_N_FLK(JL) ! Adjust the bottom temperature
      END IF
      ! 2019 fix from Mironov
    ELSE                                       ! h_ML = 0 
      PH_ML_N_FLK(JL) = PH_ML_P_FLK(JL)        ! h_ML remains unchanged 
      PC_T_N_FLK(JL)  = (PT_WML_N_FLK(JL)-PT_MNW_N_FLK(JL))/ &
     & (PT_WML_N_FLK(JL)-PT_BOT_N_FLK(JL)) 
      PC_T_N_FLK(JL)  = MIN(RC_T_MAX, MAX(PC_T_N_FLK(JL), RC_T_MIN)) ! Update the shape factor (thermocline)  
    END IF 
! Limit lake shape factor tendency to avoid instabilities
    ZD_C_T_DT=(PC_T_N_FLK(JL)-PC_T_P_FLK(JL))/PDEL_TIME
    IF ( ABS (ZD_C_T_DT) .GT. ZD_C_T_DT_MAX ) THEN
      PC_T_N_FLK(JL)=PC_T_P_FLK(JL)+SIGN(ZD_C_T_DT_MAX,ZD_C_T_DT)*PDEL_TIME
    END IF
  END IF 

  PT_BOT_N_FLK(JL) = MIN(PT_BOT_N_FLK(JL), RTPL_T_R) ! Security, limit the bottom temperature by T_r 

ELSE HTC_WATER                                      ! Open water

! Generalised buoyancy flux scale and convective velocity scale
  ZFLK_STR_1 = ZFLAKE_BUOYPAR(PT_WML_P_FLK(JL))*ZQ_STAR_FLK/RTPL_RHO_W_R/RTPL_C_W
! Fix for PH_ML_P_FLK(JL) negative
  IF(ZFLK_STR_1.LT.0._JPRD .AND. PH_ML_P_FLK(JL)>0._JPRD) THEN       
    ZW_STAR_SFC_FLK = (-ZFLK_STR_1*PH_ML_P_FLK(JL))**(1._JPRD/3._JPRD)  ! Convection     
  ELSE 
    ZW_STAR_SFC_FLK = 0._JPRD                                           ! Neutral or stable stratification
  END IF 

!_dm
! The equilibrium depth of the CBL due to surface cooling with the volumetric heating
! is not computed as a solution to the transcendental equation.
! Instead, an algebraic formula is used
! that interpolates between the two asymptotic limits.
!_dm

  ZCONV_EQUIL_H_SCALE = -PQ_W_FLK(JL)/MAX(PI_W_FLK(JL), RC_SMALL_FLK)
  IF(ZCONV_EQUIL_H_SCALE.GT.0._JPRD .AND. ZCONV_EQUIL_H_SCALE.LT.1._JPRD  &
    .AND. PT_WML_P_FLK(JL).GT.RTPL_T_R) THEN       ! The equilibrium CBL depth scale is only used above T_r
    ZCONV_EQUIL_H_SCALE = SQRT(6._JPRD*ZCONV_EQUIL_H_SCALE)                 &
                       + 2._JPRD*ZCONV_EQUIL_H_SCALE/(1._JPRD-ZCONV_EQUIL_H_SCALE)
    ZCONV_EQUIL_H_SCALE = MIN(ZDEPTH_W(JL), ZCONV_EQUIL_H_SCALE/PEXTINCOEF_WATER_TYP)
  ELSE
    ZCONV_EQUIL_H_SCALE = 0._JPRD       ! Set the equilibrium CBL depth to zero
  END IF

! Mean buoyancy frequency in the thermocline
  ZN_T_MEAN = ZFLAKE_BUOYPAR(0.5_JPRD*(PT_WML_P_FLK(JL)+PT_BOT_P_FLK(JL)))* &
 & (PT_WML_P_FLK(JL)-PT_BOT_P_FLK(JL))
  IF(PH_ML_P_FLK(JL).LE.ZDEPTH_W(JL)-RH_ML_MIN_FLK) THEN
    ZN_T_MEAN = SQRT(MAX(ZN_T_MEAN/(ZDEPTH_W(JL)-PH_ML_P_FLK(JL)),0._JPRD))  ! Compute N
  ELSE 
    ZN_T_MEAN = 0._JPRD                            ! h_ML=D, set N to zero
  END IF 

! The rate of change of C_T
  ZD_C_T_DT = MAX(ZW_STAR_SFC_FLK, PU_STAR_W_FLK(JL), RU_STAR_MIN_FLK)**2_JPIM
  ZD_C_T_DT = ZN_T_MEAN*(ZDEPTH_W(JL)-PH_ML_P_FLK(JL))**2_JPIM       &
           /RC_RELAX_C/ZD_C_T_DT                                 ! Relaxation time scale for C_T
  ZD_C_T_DT = (RC_T_MAX-RC_T_MIN)/MAX(ZD_C_T_DT, RC_SMALL_FLK)     ! Rate-of-change of C_T 

! Compute the shape factor and the mixed-layer depth, 
! using different formulations for convection and wind mixing

  ZC_TT_FLK = RC_TT_1*PC_T_P_FLK(JL)-RC_TT_2       ! C_TT, using C_T at the previous time step
  ZC_Q_FLK = 2._JPRD*ZC_TT_FLK/PC_T_P_FLK(JL)      ! C_Q using C_T at the previous time step

  MIXING_REGIME: IF(ZFLK_STR_1.LT.0._JPRD) THEN    ! Convective mixing 
  
    IF ( ABS (ZD_C_T_DT) .GT. ZD_C_T_DT_MAX ) ZD_C_T_DT=SIGN(ZD_C_T_DT_MAX,ZD_C_T_DT)
    PC_T_N_FLK(JL) = PC_T_P_FLK(JL) + ZD_C_T_DT*PDEL_TIME                ! Update C_T, assuming dh_ML/dt>0
    PC_T_N_FLK(JL) = MIN(RC_T_MAX, MAX(PC_T_N_FLK(JL), RC_T_MIN))        ! Limit C_T 
    ZD_C_T_DT = (PC_T_N_FLK(JL)-PC_T_P_FLK(JL))/PDEL_TIME                ! Re-compute dC_T/dt


    IF(PH_ML_P_FLK(JL).LE.ZDEPTH_W(JL)-RH_ML_MIN_FLK) THEN       ! Compute dh_ML/dt
      IF(PH_ML_P_FLK(JL).LE.RH_ML_MIN_FLK) THEN                  ! Use a reduced entrainment equation (spin-up)
        ZD_H_ML_DT = RC_CBL_1/RC_CBL_2*MAX(ZW_STAR_SFC_FLK, RC_SMALL_FLK)
      ELSE                                   ! Use a complete entrainment equation 
        ZR_H_ICESNOW     = ZDEPTH_W(JL)/PH_ML_P_FLK(JL)
        ZR_RHO_C_ICESNOW = ZR_H_ICESNOW-1._JPRD
        ZR_TI_ICESNOW    = PC_T_P_FLK(JL)/ZC_TT_FLK
        ZR_TSTAR_ICESNOW = (ZR_TI_ICESNOW/2._JPRD-1._JPRD)*ZR_RHO_C_ICESNOW + 1._JPRD
        ZD_H_ML_DT = -ZQ_STAR_FLK*(ZR_TSTAR_ICESNOW*(1._JPRD+RC_CBL_1)-1._JPRD) - ZQ_BOT_FLK
        ZD_H_ML_DT = ZD_H_ML_DT/RTPL_RHO_W_R/RTPL_C_W           ! Q_* and Q_b flux terms
        ZFLK_STR_2 = (ZDEPTH_W(JL)-PH_ML_P_FLK(JL))*(PT_WML_P_FLK(JL)-PT_BOT_P_FLK(JL))* &
       & RC_TT_2/ZC_TT_FLK*ZD_C_T_DT 
        ZD_H_ML_DT = ZD_H_ML_DT + ZFLK_STR_2                    ! Add dC_T/dt term
        ZFLK_STR_2 = PI_BOT_FLK(JL) + (ZR_TI_ICESNOW-1._JPRD)*PI_H_FLK(JL) - &
       & ZR_TI_ICESNOW*PI_INTM_H_D_FLK(JL)
        ZFLK_STR_2 = ZFLK_STR_2 + (ZR_TI_ICESNOW-2._JPRD)*ZR_RHO_C_ICESNOW* &
       & (PI_H_FLK(JL)-PI_INTM_0_H_FLK(JL))
        ZFLK_STR_2 = ZFLK_STR_2/RTPL_RHO_W_R/RTPL_C_W
        ZD_H_ML_DT = ZD_H_ML_DT + ZFLK_STR_2                    ! Add radiation terms
        ZFLK_STR_2 = -RC_CBL_2*ZR_TSTAR_ICESNOW*ZQ_STAR_FLK/RTPL_RHO_W_R/RTPL_C_W/ &
       & MAX(ZW_STAR_SFC_FLK, RC_SMALL_FLK)
        ZFLK_STR_2 = ZFLK_STR_2 + PC_T_P_FLK(JL)*(PT_WML_P_FLK(JL)-PT_BOT_P_FLK(JL))
        ZD_H_ML_DT = ZD_H_ML_DT/SIGN(MAX(ABS(ZFLK_STR_2),RC_SMALL_FLK),ZFLK_STR_2) ! dh_ML/dt = r.h.s.
      END IF 
!_dm
! Notice that dh_ML/dt may appear to be negative  
! (e.g. due to buoyancy loss to bottom sediments and/or
! the effect of volumetric radiation heating),
! although a negative generalized buoyancy flux scale indicates 
! that the equilibrium CBL depth has not yet been reached
! and convective deepening of the mixed layer should take place.
! Physically, this situation reflects an approximate character of the lake model.
! Using the self-similar temperature profile in the thermocline, 
! there is always communication between the mixed layer, the thermocline 
! and the lake bottom. As a result, the rate of change of the CBL depth
! is always dependent on the bottom heat flux and the radiation heating of the thermocline.
! In reality, convective mixed-layer deepening may be completely decoupled
! from the processes underneath. In order to account for this fact,
! the rate of CBL deepening is set to a small value
! if dh_ML/dt proves to be negative.
! This is "double insurance" however, 
! as a negative dh_ML/dt is encountered very rarely.
!_dm

      ZD_H_ML_DT = MAX(ZD_H_ML_DT, RC_SMALL_FLK)    
      PH_ML_N_FLK(JL) = PH_ML_P_FLK(JL) + ZD_H_ML_DT*PDEL_TIME     ! Update h_ML 
      PH_ML_N_FLK(JL) = MAX(RH_ML_MIN_FLK, MIN(PH_ML_N_FLK(JL), ZDEPTH_W(JL)))  ! Security, limit h_ML
    ELSE                                              ! Mixing down to the lake bottom
      PH_ML_N_FLK(JL) = ZDEPTH_W(JL)
    END IF

  ELSE MIXING_REGIME                                  ! Wind mixing

    ZD_H_ML_DT = MAX(PU_STAR_W_FLK(JL), RU_STAR_MIN_FLK)     ! The surface friction velocity
    ZM_H_SCALE = (ABS(PAR_CORIOLIS(JL))/RC_SBL_ZM_N + &
   & ZN_T_MEAN/RC_SBL_ZM_I)*ZD_H_ML_DT**2_JPIM
    ZM_H_SCALE = ZM_H_SCALE + ZFLK_STR_1/RC_SBL_ZM_S
    ZM_H_SCALE = MAX(ZM_H_SCALE, RC_SMALL_FLK)
    ZM_H_SCALE = ZD_H_ML_DT**3_JPIM/ZM_H_SCALE 
    ZM_H_SCALE = MAX(RH_ML_MIN_FLK, MIN(ZM_H_SCALE, RH_ML_MAX_FLK)) ! The ZM96 SBL depth scale 
    ZM_H_SCALE = MAX(ZM_H_SCALE, ZCONV_EQUIL_H_SCALE)               ! Equilibrium mixed-layer depth 

!_dm 
! In order to avoid numerical discretization problems,
! an analytical solution to the evolution equation 
! for the wind-mixed layer depth is used.
! That is, an exponential relaxation formula is applied
! over the time interval equal to the model time step.
!_dm 

    ZD_H_ML_DT = RC_RELAX_H*ZD_H_ML_DT/ZM_H_SCALE*PDEL_TIME
    PH_ML_N_FLK(JL) = ZM_h_scale - (ZM_h_scale-PH_ML_P_FLK(JL))* &
   & EXP(-ZD_H_ML_DT)    ! Update h_ML 
    PH_ML_N_FLK(JL) = MAX(RH_ML_MIN_FLK, MIN(PH_ML_N_FLK(JL), ZDEPTH_W(JL)))  ! Limit h_ML 
    ZD_H_ML_DT = (PH_ML_N_FLK(JL)-PH_ML_P_FLK(JL))/PDEL_TIME                  ! Re-compute dh_ML/dt

    IF(PH_ML_N_FLK(JL).LE.PH_ML_P_FLK(JL))           &
      ZD_C_T_DT = -ZD_C_T_DT            ! Mixed-layer retreat or stationary state, dC_T/dt<0
    IF ( ABS (ZD_C_T_DT) .GT. ZD_C_T_DT_MAX ) ZD_C_T_DT=SIGN(ZD_C_T_DT_MAX,ZD_C_T_DT)
    PC_T_N_FLK(JL) = PC_T_P_FLK(JL) + ZD_C_T_DT*PDEL_TIME          ! Update C_T
    PC_T_N_FLK(JL) = MIN(RC_T_MAX, MAX(PC_T_N_FLK(JL), RC_T_MIN))  ! Limit C_T 
    ZD_C_T_DT = (PC_T_N_FLK(JL)-PC_T_P_FLK(JL))/PDEL_TIME          ! Re-compute dC_T/dt

  END IF MIXING_REGIME

! Compute the time-rate-of-change of the the bottom temperature, 
! depending on the sign of dh_ML/dt 
! Update the bottom temperature and the mixed-layer temperature

  IF(PH_ML_N_FLK(JL).LE.ZDEPTH_W(JL)-RH_ML_MIN_FLK) THEN  ! Mixing did not reach the bottom 

    IF(PH_ML_N_FLK(JL).GT.PH_ML_P_FLK(JL)) THEN   ! Mixed-layer deepening 
      ZR_H_ICESNOW     = PH_ML_P_FLK(JL)/ZDEPTH_W(JL)
      ZR_RHO_C_ICESNOW = 1._JPRD-ZR_H_ICESNOW 
      ZR_TI_ICESNOW    = 0.5_JPRD*PC_T_P_FLK(JL)*ZR_RHO_C_ICESNOW+ZC_TT_FLK* &
     & (2._JPRD*ZR_H_ICESNOW-1._JPRD)
      ZR_TSTAR_ICESNOW = (0.5_JPRD+ZC_TT_FLK-ZC_Q_FLK)/ZR_TI_ICESNOW
      ZR_TI_ICESNOW    = (1._JPRD-PC_T_P_FLK(JL)*ZR_RHO_C_ICESNOW)/ZR_TI_ICESNOW
     
      ZD_T_BOT_DT = (PQ_W_FLK(JL)-ZQ_BOT_FLK+PI_W_FLK(JL)-PI_BOT_FLK(JL))/ &
     & RTPL_RHO_W_R/RTPL_C_W
      ZD_T_BOT_DT = ZD_T_BOT_DT - PC_T_P_FLK(JL)*(PT_WML_P_FLK(JL)-PT_BOT_P_FLK(JL))* &
     & ZD_H_ML_DT
      ZD_T_BOT_DT = ZD_T_BOT_DT*ZR_TSTAR_ICESNOW/ZDEPTH_W(JL) ! Q+I fluxes and dh_ML/dt term

      ZFLK_STR_2 = PI_INTM_H_D_FLK(JL) - (1._JPRD-ZC_Q_FLK)*PI_H_FLK(JL) - &
     & ZC_Q_FLK*PI_BOT_FLK(JL)
      ZFLK_STR_2 = ZFLK_STR_2*ZR_TI_ICESNOW/(ZDEPTH_W(JL)-PH_ML_P_FLK(JL))/ &
     & RTPL_RHO_W_R/RTPL_C_W
      ZD_T_BOT_DT = ZD_T_BOT_DT + ZFLK_STR_2       ! Add radiation-flux term

      ZFLK_STR_2 = (1._JPRD-RC_TT_2*ZR_TI_ICESNOW)/PC_T_P_FLK(JL)
      ZFLK_STR_2 = ZFLK_STR_2*(PT_WML_P_FLK(JL)-PT_BOT_P_FLK(JL))*ZD_C_T_DT
      ZD_T_BOT_DT = ZD_T_BOT_DT + ZFLK_STR_2       ! Add dC_T/dt term
      
    ELSE                                ! Mixed-layer retreat or stationary state
      ZD_T_BOT_DT = 0._JPRD             ! dT_bot/dt=0
    END IF

    PT_BOT_N_FLK(JL) = PT_BOT_P_FLK(JL) + ZD_T_BOT_DT*PDEL_TIME  ! Update T_bot  
    PT_BOT_N_FLK(JL) = MAX(PT_BOT_N_FLK(JL), RTPL_T_F)           ! Security, limit T_bot by the freezing point
    ZFLK_STR_2 = (PT_BOT_N_FLK(JL)-RTPL_T_R)*ZFLAKE_BUOYPAR(PT_MNW_N_FLK(JL))
    IF(ZFLK_STR_2.LT.0._JPRD) PT_BOT_N_FLK(JL) = RTPL_T_R        ! Security, avoid T_r crossover 
    !*! 2019 fix from Mironov Not active
    !*ZFLK_STR_2 = (PT_BOT_N_FLK(JL)-RTPL_T_R)*(PT_MNW_N_FLK(JL)-RTPL_T_R)
    !*IF(ZFLK_STR_2 <= 0._JPRD) PT_BOT_N_FLK(JL) =  RTPL_T_R ! Security, avoid T_r crossover 
    !*! 2019 fix from Mironov

    PT_WML_N_FLK(JL) = PC_T_N_FLK(JL)*(1._JPRD-PH_ML_N_FLK(JL)/ZDEPTH_W(JL))
    PT_WML_N_FLK(JL) = (PT_MNW_N_FLK(JL)-PT_BOT_N_FLK(JL)*PT_WML_N_FLK(JL))/ &
   & (1._JPRD-PT_WML_N_FLK(JL))
    PT_WML_N_FLK(JL) = MAX(PT_WML_N_FLK(JL), RTPL_T_F)           ! Security, limit T_wML by the freezing point

  ELSE                                                           ! Mixing down to the lake bottom 

    PH_ML_N_FLK(JL)  = ZDEPTH_W(JL)
    PT_WML_N_FLK(JL) = PT_MNW_N_FLK(JL)
    PT_BOT_N_FLK(JL) = PT_MNW_N_FLK(JL)
    PC_T_N_FLK(JL)   = ZC_T_CON
  END IF

!*  ! 2019 fix from Mironov - not active, also a repetition of what below ln 696
!*  ! In case of unstable stratification, force mixing down to the bottom
!*  ZFLK_STR_2 = (PT_WML_N_FLK(JL)-PT_BOT_N_FLK(JL))*(PT_MNW_N_FLK(JL)-RTPL_T_R)
!*  IF(ZFLK_STR_2 < 0._JPRD) THEN
!*    PH_ML_N_FLK(JL)  = ZDEPTH_W(JL)
!*    PT_WML_N_FLK(JL) = PT_MNW_N_FLK(JL)
!*    PT_BOT_N_FLK(JL) = PT_MNW_N_FLK(JL)
!*    PC_T_N_FLK(JL)   = ZC_T_CON
!*  END IF
!*  ! 2019 fix from Mironov
END IF HTC_WATER


!------------------------------------------------------------------------------
!  Impose additional constraints.
!------------------------------------------------------------------------------

! In case of unstable stratification, force mixing down to the bottom
ZFLK_STR_2 = (PT_WML_N_FLK(JL)-PT_BOT_N_FLK(JL))*ZFLAKE_BUOYPAR(PT_MNW_N_FLK(JL))
IF(ZFLK_STR_2.LT.0._JPRD) THEN 
  PH_ML_N_FLK(JL)  = ZDEPTH_W(JL)
  PT_WML_N_FLK(JL) = PT_MNW_N_FLK(JL)
  PT_BOT_N_FLK(JL) = PT_MNW_N_FLK(JL)
  PC_T_N_FLK(JL)   = ZC_T_CON
END IF

!------------------------------------------------------------------------------
!  Update the surface temperature.
!------------------------------------------------------------------------------

IF(PH_ICE_N_FLK(JL).GE.RH_ICE_MIN_FLK) THEN
  PT_SFC_N(JL) = PT_ICE_N_FLK(JL) ! Ice exists but there is no snow, use the ice temperature
ELSE 
  PT_SFC_N(JL) = PT_WML_N_FLK(JL) ! No ice-snow cover, use the mixed-layer temperature
END IF

ENDIF LAKEPOINT

ENDDO MAIN_LOOP

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FLAKEENE_MOD:FLAKEENE',1,ZHOOK_HANDLE)

END SUBROUTINE FLAKEENE
END MODULE FLAKEENE_MOD
