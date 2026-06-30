! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE YOS_FLAKE

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
!  Values of empirical constants of the lake model FLake 
!  and of several thermodynamic parameters are set.
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

USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD

!==============================================================================

IMPLICIT NONE

SAVE

!==============================================================================
!

INTEGER (KIND = JPIM), PARAMETER :: NBAND_OPTIC_MAX = 6_JPIM

! Define TYPE "ropticpar_medium"
TYPE ROPTICPAR_MEDIUM
  ! Number of wave-length bands
  INTEGER (KIND = JPIM):: NBAND_OPTIC  
  ! Fractions of total shortwave radiation flux 
  REAL (KIND = JPRB), DIMENSION (NBAND_OPTIC_MAX):: RFRAC_OPTIC
  ! Extinction coefficients 
  REAL (KIND = JPRB), DIMENSION (NBAND_OPTIC_MAX):: REXTINCOEF_OPTIC 
END TYPE ROPTICPAR_MEDIUM


TYPE :: TFLAKE

! Declarations

!  Dimensionless constants 
!  in the equations for the mixed-layer depth 
!  and for the shape factor with respect to the temperature profile in the thermocline

REAL (KIND = JPRB):: RC_CBL_1         ! Constant in the CBL entrainment equation
REAL (KIND = JPRB):: RC_CBL_2         ! Constant in the CBL entrainment equation
REAL (KIND = JPRB):: RC_SBL_ZM_N      ! Constant in the ZM1996 equation for the equilibrium SBL depth
REAL (KIND = JPRB):: RC_SBL_ZM_S      ! Constant in the ZM1996 equation for the equilibrium SBL depth
REAL (KIND = JPRB):: RC_SBL_ZM_I      ! Constant in the ZM1996 equation for the equilibrium SBL depth
REAL (KIND = JPRB):: RC_RELAX_H       ! Constant in the relaxation equation for the SBL depth
REAL (KIND = JPRB):: RC_RELAX_C       ! Constant in the relaxation equation for the shape factor
                                     ! with respect to the temperature profile in the thermocline

!  Parameters of the shape functions 
!  Indices refer to T - thermocline, S - snow, I - ice,
!  B1 - upper layer of the bottom sediments, B2 - lower layer of the bottom sediments.
!  "pr0" and "pr1" denote zeta derivatives of the corresponding shape function 
!  at "zeta=0" ad "zeta=1", respectively.

REAL (KIND = JPRB):: RC_T_MIN         ! Minimum value of the shape factor C_T (thermocline)
REAL (KIND = JPRB):: RC_T_MAX         ! Maximum value of the shape factor C_T (thermocline)
REAL (KIND = JPRB):: RPHI_T_PR0_1     ! Constant in the expression for the T shape-function derivative 
REAL (KIND = JPRB):: RPHI_T_PR0_2     ! Constant in the expression for the T shape-function derivative 
REAL (KIND = JPRB):: RC_TT_1          ! Constant in the expression for C_TT (thermocline)
REAL (KIND = JPRB):: RC_TT_2          ! Constant in the expression for C_TT (thermocline)
REAL (KIND = JPRB):: RC_B1            ! Shape factor (upper layer of bottom sediments)
REAL (KIND = JPRB):: RC_B2            ! Shape factor (lower layer of bottom sediments)
REAL (KIND = JPRB):: RC_S_LIN         ! Shape factor (linear temperature profile in the snow layer)
REAL (KIND = JPRB):: RPHI_S_PR0_LIN   ! S shape-function derivative (linear profile) 
REAL (KIND = JPRB):: RC_I_LIN         ! Shape factor (linear temperature profile in the ice layer)
REAL (KIND = JPRB):: RPHI_I_PR0_LIN   ! I shape-function derivative (linear profile) 
REAL (KIND = JPRB):: RPHI_I_PR1_LIN   ! I shape-function derivative (linear profile) 
REAL (KIND = JPRB):: RPHI_I_AST_MR    ! Constant in the MR2004 expression for I shape factor
REAL (KIND = JPRB):: RC_I_MR          ! Constant in the MR2004 expression for I shape factor
REAL (KIND = JPRB):: RH_ICE_MAX       ! Maximum ice tickness in 
                                     ! the Mironov and Ritter (2004, MR2004) ice model [m] 

!  Security constants

REAL (KIND = JPRB):: RH_SNOW_MIN_FLK  ! Minimum snow thickness [m]
REAL (KIND = JPRB):: RH_ICE_MIN_FLK   ! Minimum ice thickness [m]
REAL (KIND = JPRB):: RH_ML_MIN_FLK    ! Minimum mixed-layer depth [m]
REAL (KIND = JPRB):: RH_ML_MAX_FLK    ! Maximum mixed-layer depth [m]
REAL (KIND = JPRB):: RU_STAR_MIN_FLK  ! Minimum value of the surface friction velocity [m s^{-1}]

!  Security constant(s)

REAL (KIND = JPRD):: RC_SMALL_FLK  ! A small number

!  Thermodynamic parameters

REAL (KIND = JPRB):: RTPL_GRAV     ! Acceleration due to gravity [m s^{-2}]
REAL (KIND = JPRB):: RTPL_T_R      ! Temperature of maximum density of fresh water [K]
REAL (KIND = JPRB):: RTPL_T_F      ! Fresh water freezing point [K]
REAL (KIND = JPRB):: RTPL_A_T      ! Constant in the fresh-water equation of state [K^{-2}]
REAL (KIND = JPRB):: RTPL_RHO_W_R  ! Maximum density of fresh water [kg m^{-3}]
REAL (KIND = JPRB):: RTPL_RHO_I    ! Density of ice [kg m^{-3}]
REAL (KIND = JPRB):: RTPL_L_F      ! Latent heat of fusion [J kg^{-1}]
REAL (KIND = JPRB):: RTPL_C_W      ! Specific heat of water [J kg^{-1} K^{-1}]
REAL (KIND = JPRB):: RTPL_C_I      ! Specific heat of ice [J kg^{-1} K^{-1}]
REAL (KIND = JPRB):: RTPL_KAPPA_W  ! Molecular heat conductivity of water [J m^{-1} s^{-1} K^{-1}]
REAL (KIND = JPRB):: RTPL_KAPPA_I  ! Molecular heat conductivity of ice [J m^{-1} s^{-1} K^{-1}]
REAL (KIND = JPRB):: RTPSF_L_EVAP  ! Specific heat of evaporation [J kg^{-1}]

! Optical characteristics for water and ice. 
! The simplest one-band approximation is used as a reference.
TYPE (ROPTICPAR_MEDIUM):: ROPTICPAR_WATER_REF      ! Water (reference)
TYPE (ROPTICPAR_MEDIUM):: ROPTICPAR_WATER_TRANS    ! Transparent Water (two-band)
TYPE (ROPTICPAR_MEDIUM):: ROPTICPAR_WHITEICE_REF   ! White ice
TYPE (ROPTICPAR_MEDIUM):: ROPTICPAR_BLUEICE_REF    ! Blue ice
TYPE (ROPTICPAR_MEDIUM):: ROPTICPAR_ICE_OPAQUE     ! Opaque ice

LOGICAL :: LEFLAKE
    
END TYPE TFLAKE

!==============================================================================

END MODULE YOS_FLAKE

