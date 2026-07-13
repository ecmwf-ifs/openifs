! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE SUSFLAKE_MOD

CONTAINS
SUBROUTINE SUSFLAKE(LD_LEFLAKE,YDFLAKE)

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
!  Initialization of YOS_FLAKE 
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
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!            17-Dec-2015 F. Vana
!    Single precision support
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
USE YOS_FLAKE, ONLY : TFLAKE, NBAND_OPTIC_MAX, ROPTICPAR_MEDIUM

!==============================================================================

IMPLICIT NONE

LOGICAL,      INTENT(IN)    :: LD_LEFLAKE 
TYPE(TFLAKE), INTENT(INOUT) :: YDFLAKE

!==============================================================================
!
! Declarations

INTEGER (KIND = JPIM) :: J ! DO loop index
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUSFLAKE_MOD:SUSFLAKE',0,ZHOOK_HANDLE)
ASSOCIATE(LEFLAKE=>YDFLAKE%LEFLAKE, RC_B1=>YDFLAKE%RC_B1, RC_B2=>YDFLAKE%RC_B2, &
 & RC_CBL_1=>YDFLAKE%RC_CBL_1, RC_CBL_2=>YDFLAKE%RC_CBL_2, &
 & RC_I_LIN=>YDFLAKE%RC_I_LIN, RC_I_MR=>YDFLAKE%RC_I_MR, &
 & RC_RELAX_C=>YDFLAKE%RC_RELAX_C, RC_RELAX_H=>YDFLAKE%RC_RELAX_H, &
 & RC_SBL_ZM_I=>YDFLAKE%RC_SBL_ZM_I, RC_SBL_ZM_N=>YDFLAKE%RC_SBL_ZM_N, &
 & RC_SBL_ZM_S=>YDFLAKE%RC_SBL_ZM_S, RC_SMALL_FLK=>YDFLAKE%RC_SMALL_FLK, &
 & RC_S_LIN=>YDFLAKE%RC_S_LIN, RC_TT_1=>YDFLAKE%RC_TT_1, &
 & RC_TT_2=>YDFLAKE%RC_TT_2, RC_T_MAX=>YDFLAKE%RC_T_MAX, &
 & RC_T_MIN=>YDFLAKE%RC_T_MIN, RH_ICE_MAX=>YDFLAKE%RH_ICE_MAX, &
 & RH_ICE_MIN_FLK=>YDFLAKE%RH_ICE_MIN_FLK, RH_ML_MAX_FLK=>YDFLAKE%RH_ML_MAX_FLK, &
 & RH_ML_MIN_FLK=>YDFLAKE%RH_ML_MIN_FLK, &
 & RH_SNOW_MIN_FLK=>YDFLAKE%RH_SNOW_MIN_FLK, &
 & ROPTICPAR_BLUEICE_REF=>YDFLAKE%ROPTICPAR_BLUEICE_REF, &
 & ROPTICPAR_ICE_OPAQUE=>YDFLAKE%ROPTICPAR_ICE_OPAQUE, &
 & ROPTICPAR_WATER_REF=>YDFLAKE%ROPTICPAR_WATER_REF, &
 & ROPTICPAR_WATER_TRANS=>YDFLAKE%ROPTICPAR_WATER_TRANS, &
 & ROPTICPAR_WHITEICE_REF=>YDFLAKE%ROPTICPAR_WHITEICE_REF, &
 & RPHI_I_AST_MR=>YDFLAKE%RPHI_I_AST_MR, RPHI_I_PR0_LIN=>YDFLAKE%RPHI_I_PR0_LIN, &
 & RPHI_I_PR1_LIN=>YDFLAKE%RPHI_I_PR1_LIN, &
 & RPHI_S_PR0_LIN=>YDFLAKE%RPHI_S_PR0_LIN, RPHI_T_PR0_1=>YDFLAKE%RPHI_T_PR0_1, &
 & RPHI_T_PR0_2=>YDFLAKE%RPHI_T_PR0_2, RTPL_A_T=>YDFLAKE%RTPL_A_T, &
 & RTPL_C_I=>YDFLAKE%RTPL_C_I, RTPL_C_W=>YDFLAKE%RTPL_C_W, &
 & RTPL_GRAV=>YDFLAKE%RTPL_GRAV, RTPL_KAPPA_I=>YDFLAKE%RTPL_KAPPA_I, &
 & RTPL_KAPPA_W=>YDFLAKE%RTPL_KAPPA_W, RTPL_L_F=>YDFLAKE%RTPL_L_F, &
 & RTPL_RHO_I=>YDFLAKE%RTPL_RHO_I, RTPL_RHO_W_R=>YDFLAKE%RTPL_RHO_W_R, &
 & RTPL_T_F=>YDFLAKE%RTPL_T_F, RTPL_T_R=>YDFLAKE%RTPL_T_R, &
 & RTPSF_L_EVAP=>YDFLAKE%RTPSF_L_EVAP, &
 & RU_STAR_MIN_FLK=>YDFLAKE%RU_STAR_MIN_FLK)

!  Dimensionless constants 
!  in the equations for the mixed-layer depth 
!  and for the shape factor with respect to the temperature profile in the thermocline
RC_CBL_1       = 0.17_JPRD             ! Constant in the CBL entrainment equation
RC_CBL_2       = 1._JPRD               ! Constant in the CBL entrainment equation
RC_SBL_ZM_N    = 0.5_JPRD              ! Constant in the ZM1996 equation for the equilibrium SBL depth
RC_SBL_ZM_S    = 10._JPRD              ! Constant in the ZM1996 equation for the equilibrium SBL depth
RC_SBL_ZM_I    = 20._JPRD              ! Constant in the ZM1996 equation for the equilibrium SBL depth
RC_RELAX_H     = 0.030_JPRD            ! Constant in the relaxation equation for the SBL depth
RC_RELAX_C     = 0.030_JPRD           ! Constant with respect to the temperature profile in the thermocline

!  Parameters of the shape functions 
!  Indices refer to T - thermocline, S - snow, I - ice,
!  B1 - upper layer of the bottom sediments, B2 - lower layer of the bottom sediments.
!  "pr0" and "pr1" denote zeta derivatives of the corresponding shape function 
!  at "zeta=0" ad "zeta=1", respectively.

RC_T_MIN       = 0.5_JPRD            ! Minimum value of the shape factor C_T (thermocline)
RC_T_MAX       = 0.8_JPRD            ! Maximum value of the shape factor C_T (thermocline)
RPHI_T_PR0_1   = 40._JPRD/3._JPRD    ! Constant in the expression for the T shape-function derivative 
RPHI_T_PR0_2   = 20._JPRD/3._JPRD    ! Constant in the expression for the T shape-function derivative 
RC_TT_1        = 11._JPRD/18._JPRD   ! Constant in the expression for C_TT (thermocline)
RC_TT_2        = 7._JPRD/45._JPRD    ! Constant in the expression for C_TT (thermocline)
RC_B1          = 2._JPRD/3._JPRD     ! Shape factor (upper layer of bottom sediments)
RC_B2          = 3._JPRD/5._JPRD     ! Shape factor (lower layer of bottom sediments)
RC_S_LIN       = 0.5_JPRD            ! Shape factor (linear temperature profile in the snow layer)
RPHI_S_PR0_LIN = 1._JPRD             ! S shape-function derivative (linear profile) 
RC_I_LIN       = 0.5_JPRD            ! Shape factor (linear temperature profile in the ice layer)
RPHI_I_PR0_LIN = 1._JPRD             ! I shape-function derivative (linear profile) 
RPHI_I_PR1_LIN = 1._JPRD             ! I shape-function derivative (linear profile) 
RPHI_I_AST_MR  = 2._JPRD             ! Constant in the MR2004 expression for I shape factor
RC_I_MR        = 1._JPRD/12._JPRD    ! Constant in the MR2004 expression for I shape factor
RH_ICE_MAX     = 3._JPRD             ! Maximum ice tickness in 
                                    ! the Mironov and Ritter (2004, MR2004) ice model [m] 

!  Security constants

RH_SNOW_MIN_FLK = 1.0E-5_JPRD          ! Minimum snow thickness [m]
RH_ICE_MIN_FLK  = 1.0E-9_JPRD          ! Minimum ice thickness [m]
RH_ML_MIN_FLK   = 1.0E-2_JPRD          ! Minimum mixed-layer depth [m]
RH_ML_MAX_FLK   = 1.0E+3_JPRD          ! Maximum mixed-layer depth [m]
RU_STAR_MIN_FLK = 1.0E-6_JPRD          ! Minimum value of the surface friction velocity [m s^{-1}]

!  Security constant(s)

RC_SMALL_FLK    = 100._JPRD*EPSILON(RC_SMALL_FLK)  ! A small number

!  Thermodynamic parameters

RTPL_GRAV          = 9.81_JPRD        ! Acceleration due to gravity [m s^{-2}]
RTPL_T_R           = 277.13_JPRD      ! Temperature of maximum density of fresh water [K]
RTPL_T_F           = 273.15_JPRD      ! Fresh water freezing point [K]
RTPL_A_T           = 1.6509E-05_JPRD  ! Constant in the fresh-water equation of state [K^{-2}]
RTPL_RHO_W_R       = 1.0E+03_JPRD     ! Maximum density of fresh water [kg m^{-3}]
RTPL_RHO_I         = 9.1E+02_JPRD     ! Density of ice [kg m^{-3}]
RTPL_L_F           = 3.3E+05_JPRD     ! Latent heat of fusion [J kg^{-1}]
RTPL_C_W           = 4.2E+03_JPRD     ! Specific heat of water [J kg^{-1} K^{-1}]
RTPL_C_I           = 2.1E+03_JPRD     ! Specific heat of ice [J kg^{-1} K^{-1}]
RTPL_KAPPA_W       = 5.46E-01_JPRD    ! Molecular heat conductivity of water [J m^{-1} s^{-1} K^{-1}]
RTPL_KAPPA_I       = 2.29_JPRD        ! Molecular heat conductivity of ice [J m^{-1} s^{-1} K^{-1}]
RTPSF_L_EVAP       = 2.501E+06_JPRD   ! Specific heat of evaporation [J kg^{-1}]


! Optical characteristics for water, ice and snow.
! The simplest one-band approximation is used as a reference.

ROPTICPAR_WATER_REF = ROPTICPAR_MEDIUM(1,                          & ! Water (reference)
& (/1._JPRD, (0._JPRD    , J=2,NBAND_OPTIC_MAX)/),               &
& (/3._JPRD, (1.E+10_JPRD, J=2,NBAND_OPTIC_MAX)/))                  ! 0.5 m-1 in Lemoigne, 2015, (pers. comm.) 1.0 m-1 in Mironov 2015 (pers.comm.)
!& (/1._JPRD, (1.E+10_JPRD, J=2,NBAND_OPTIC_MAX)/))
ROPTICPAR_WATER_TRANS = ROPTICPAR_MEDIUM(2,                        & ! Transparent Water (two-band)
& (/0.10_JPRD, 0.90_JPRD, (0._JPRD    , J=3,NBAND_OPTIC_MAX)/),  &
& (/2.0_JPRD , 0.20_JPRD, (1.E+10_JPRD, J=3,NBAND_OPTIC_MAX)/))
ROPTICPAR_WHITEICE_REF = ROPTICPAR_MEDIUM(1,                       & ! White ice
& (/1._JPRD  , (0._JPRD    , J=2,NBAND_OPTIC_MAX)/),             &
& (/17.1_JPRD, (1.E+10_JPRD, J=2,NBAND_OPTIC_MAX)/))
ROPTICPAR_BLUEICE_REF = ROPTICPAR_MEDIUM(1,                        & ! Blue ice
& (/1._JPRD , (0._JPRD    , J=2,NBAND_OPTIC_MAX)/),              &
& (/8.4_JPRD, (1.E+10_JPRD, J=2,NBAND_OPTIC_MAX)/))
ROPTICPAR_ICE_OPAQUE = ROPTICPAR_MEDIUM(1,                         & ! Opaque ice
& (/1._JPRD     , (0._JPRD    , J=2,NBAND_OPTIC_MAX)/),          &
& (/1.0E+07_JPRD, (1.E+10_JPRD, J=2,NBAND_OPTIC_MAX)/))

LEFLAKE=LD_LEFLAKE
 
!==============================================================================

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSFLAKE_MOD:SUSFLAKE',1,ZHOOK_HANDLE)

END SUBROUTINE SUSFLAKE
END MODULE SUSFLAKE_MOD

