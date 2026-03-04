! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE YOMSEKF

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!     N_OBSTIME_MAX : Maximum number of observation times in one sekf
!                     analysis window
!     N_ANALTIME_MAX: Maximum number of analysis times in one sekf
!                     analysis window
!     N_SOILM_MAX:    Maximum number of analysed soil layers (moisture)

!     ------------------------------------------------------------------
!*    Variables for the simplfied EKF soil moisture analysis.
!*    ------------------------------------------------------------------

!     Control variables for the EKF
!     -----------------------------
!     N_SEKF_CV  :    Number of control variables in the SEKF
!     N_SEKF_PT  :    Number of perturbation run
!     N_SEKF_OUPT:    Number of unperturbed runs
!     NSLAY      :    Number of perturbed (= analysed soil layers)
!     NCOM       :    Number of ocean mix-layer model levels

!     Debugg variables for output
!     ---------------------------
!     NREF_BLOCK
!     NREF_GP

!     Perturbation runs
!     -----------------
!     VSM_PERT_INC:   Increments for soil moisture perturbation
!     PERT_CHESS:     Use chessboard type perturbations

!     Save tendencies and selected surface forcing from the reference fc
!     ------------------------------------------------------------------
!     FKF_TENT    :    Temperature tendencies from the moist physics parts
!                     (CALLPAR.F90 sections 5 and 6)
!     FKF_TENQ    :    Moisture tendencies from the moist physics parts
!     FKF_TENU    :    Wind tendencies from the moist physics parts
!     FKF_TENV    :    Wind tendencies from the moist physics parts
!                      KF_TEN* are global 4D (NPROMA, NLEV, NBLOCK, NSTEP)
!     FKF_SURF_SO :    Net surface radiation (solar) at the surface (W/m2)
!     FKF_SURF_TH :    Net surface radiation (thermal) at the surface (W/m2)
!     FKF_SURF_CR :    Convective rain at the surface
!     FKF_SURF_LR :    Large scale rain at the surface
!                      KF_SURF_* are global 3D (NPROMA, NBLOCK, NSTEP)

!     Save variable for off-line Jacobians in callpar
!     -----------------------------------------------

! FOR VDFOUTER:
! ~~~~~~~~~~~~
! SKF_ZWSA       : MULTI-LAYER SOIL WETNESS
! SKF_ZTSA       : MULTI-LAYER SOIL TEMPERATURE

! SKF_PVDIS      : TURBULENT DISSIPATION
! SKF_PVDISG     : GRAVITY WAVE DRAG DISSIPATION
! SKF_PFTLHEV    : SURFACE LATENT HEAT FLUX (SNOW FREE FRACTION)
! SKF_PFTLHSB    : SURFACE LATENT HEAT FLUX (SNOW COVERED FRACTION)
! SKF_PFWSB      : SURFACE WATER VAPOUR FLUX (SUBLIMATION)
! SKF_PUCFL      : U-COMPONENT OF WIND AT 10 METERS (DIAGNOSTIC).
! SKF_PVCFL      : V-COMPONENT OF WIND AT 10 METERS (DIAGNOSTIC).
! SKF_PTCFL      : TEMPERATURE AT 2 METERS (DIAGNOSTIC).
! SKF_PDCFL      : DEW POINT TEMPERATURE AT 2 METERS (DIAGNOSTIC).
! SKF_PTENT_CUM  : TEMPERATURE PROFILE cumulated (sum)
! SKF_PTENQ_CUM  : Humidity profile cumulated (sum)
! SKF_PQCFL      : SPECIFIC HUMIDITY AT 2 METERS (DIAGNOSTIC).
! SKF_PBLH       : BOUNDARY LAYER HEIGHT
! SKF_PU10N      : U-COMPONENT OF NEUTRAL WIND AT 10 METERS (DIAGNOSTIC).
! SKF_PV10N      : V-COMPONENT OF NEUTRAL WIND AT 10 METERS (DIAGNOSTIC).
! SKF_PUST       : FRICTION VELOCITY (DIAGNOSTIC).
! SKF_PTENT      : TENDENCY OF TEMPERATURE.
! SKF_PTENQ      : TENDENCY OF HUMIDITY
! SKF_PTENU      : TENDENCY OF U-COMP. OF WIND.
! SKF_PTENV      : TENDENCY OF V-COMP. OF WIND.
! SKF_PTLE1      : OF SKIN TEMPERATURE
! SKF_PTL        : SKIN TEMPERATURE
! SKF_PUSTRTI    : E-W (INSTANTANEOUS) SURFACE STRESS FOR EACH TILE
! SKF_PVSTRTI    : N-S (INSTANTANEOUS) SURFACE STRESS FOR EACH TILE
! SKF_PAHFSTI    : (INSTANTANEOUS) SURFACE SENSIBLE HEAT FLUX FOR EACH TILE
! SKF_PEVAPTI       : (INSTANTANEOUS) EVAPORATION FOR EACH TILE
! SKF_PTSKTI     : SKIN TEMPERATURE FOR EACH TILE
! SKF_PDIFTS     : TURBULENT FLUX OF ENTHALPY (OR DRY STATIC ENERGY).
! SKF_PDIFTQ     : TURBULENT FLUX (INC. Q NEGATIVE) OF SPECIFIC HUMIDITY.
! SKF_PDIFTL     : TURBULENT FLUX (INC. Q NEGATIVE) OF LIQUID WATER.
! SKF_PDIFTI     : TURBULENT FLUX (INC. Q NEGATIVE) OF SOLID WATER.
! SKF_PSTRTU     : TURBULENT FLUX OF MOMENTUM "U".
! SKF_PSTRTV     : TURBULENT FLUX OF MOMENTUM "V".
! SKF_PSTRDU     : GRAVITY WAVE DRAG FLUX "U".
! SKF_PSTRDV     : GRAVITY WAVE DRAG FLUX "V".
! SKF_ZKH_VDF    : TURBULENT DIFFUSION COEFFICIENT FOR HEAT
!! SKF_PDHTLS     : Diagnostic array for tiles (see module yomcdh)
!! SKF_PDHTSS     : Diagnostic array for snow T (see module yomcdh)
!! SKF_PDHTTS     : Diagnostic array for soil T (see module yomcdh)
!! SKF_PDHTIS     : Diagnostic array for ice T (see module yomcdh)

! FOR SURFTSTP
!~~~~~~~~~~~~~
! SKF_PZO,SKF_PHO,SKF_PHO_INV,SKF_PDO,SKF_PADVT,SKF_PADVS, : Ocean model
! SKF_PTRI0,SKF_PTRI1,SKF_PSWDK_SAVE,SKF_PUO0,SKF_PVO0,SKF_PTO0,SKF_PSO0 : Ocean model
! SKF_PFTG12     UPWARD FLUX BETWEEN SURFACE AND DEEP LAYER   W/M**2
! SKF_PFWROD     DEEP LAYER RUN-OFF                          kg/m**2/s
! SKF_PFWRO1     SURFACE RUN-OFF                             kg/m**2/s
! SKF_PFWG12     WATER FLUX BETWEEN LAYER 1 AND 2            kg/m**2/s
! SKF_PFWMLT     WATER FLUX CORRESPONDING TO SNOW MELT       kg/m**2/s
! SKF_PFWEV      EVAPORATION OVER LAND SURFACE               kg/m**2/s
! SKF_PENES      SOIL ENERGY per unit area                   J/M**2
! KPP OCEAN MODEL AND FLAKE
! SKF_PDIFM      KPP Ocean model
! SKF_PDIFT
! SKF_PDIFS
! SKF_PTLICEE1,SKF_PTLMNWE1,SKF_PTLWMLE1,SKF_PTLBOTE1,SKF_PTLSFE1&
! SKF_PHLICEE1,SKF_PHLMLE1,SKF_PUOE1,SKF_PVOE1,SKF_PTOE1,SKF_PSOE1
! SKF_PUOC,SKF_PVOC,SKF_PUSTRC,SKF_PVSTRC


! SKF_PDHTSS     Diagnostic array for snow T (see module yomcdh)
! SKF_PDHTTS     Diagnostic array for soil T (see module yomcdh)
! SKF_PDHTIS     Diagnostic array for ice T (see module yomcdh)
! SKF_PDHSSS     Diagnostic array for snow mass (see module yomcdh)
! SKF_PDHIIS     Diagnostic array for interception layer (see module yomcdh)
! SKF_PDHWLS     Diagnostic array for soil water (see module yomcdh)

!     Store bg fields from the reference and the perturbed runs
!     ---------------------------------------------------------
!     VOL_SM_FG  :     Volumetric soil moisture from the fg & pert. runs
!     T2M_FG     :     Temperature at 2 m from the fg & pert. runs
!     TD2M_FG    :     Dew point temperature at 2 m from the fg & pert. runs
!                               *FG are global 5D (NPROMA, NSLAY, NBLOCK, NSTEP, N_SEKF_CV+1)
!     TTB_FG     :     CMEM radiances from the fg & pert. runs
!     STYPE      :     Soil type for 3D-error structure in the SEKF background error matrix 


!     Observation screening and spatial / temporal matching
!     -----------------------------------------------------
!     NSSA_TSLOT     : Number of screen level observations time slots within time window
!     NSMOS_SM_FIL   : Number of SMOS neural network soil moisture observations files within time window
!     LSMOS_SM_OBS   : presence of SMOS SM observation for a given point and time
!     NOBS_SMOS_SM   : Number of SMOS neural network soil moisture observations on each grid point
!     NOBS_SMOS_TB   : Number of SMOS TB observations on each grid point
!     NOBS_ASCAT     : Number of ASCAT soil moisture observations on each grid point
!     N_EDAMEMBERS   : Number of EDA members used to compute the jacobians
!     NOBS_SEKF_MAX  : Max. number of observations on each grid point
!     NOBS_SIM_EDA   : Number of observation model equivalent fields needed from the EDA to compute the Jacobians 
!     NOBS_ASCAT_MAX : Max. number of ASCAT observations
!     NGRIB_FIELDS_SCREEN : number of distinct grib fields input screen level observations in the SEKF
!     NGRIB_FIELDS_SMOS : number of distinct grib fields input SMOS SM observations in the SEKF
!     NANAL          : Number of analyses
!     NSCREENTIME    : 'Observation times' for the screen level parameters (HH)
!     NANATIME       : Analyses times within the time window (HH)
!     NSCREENDATE    : 'Observation times' for the screen level parameters (YYYYMMDD)
!     NANADATE       : Analyses times within the time window (YYYYMMDD)
!     NOBS_TSTEP     : Model time step for observations
!     NANA_TSTEP     : Model time step for analyses
!     T2M_OBS        : 2m temperature observations (analysis)
!     TD2M_OBS       : dew point temperature observations (analysis)
!                               global 3D (NPROMA, NBLOCK, NOBS_SEKF_T/H_MAX)
!     SMOS_SM_OBS    : SMOS neural network soil moisture observation value
!     SMOS_SM_ERR    : SMOS neural network soil moisture observation error
!     SMOS_SM_TST   : SMOS neural network soil moisture observation time step
!     COVAR_RH2M_CV_SEKF: covariance of EDA soil moisture and RH2m (at each obs time and for each layer)
!     COVAR_T2M_CV_SEKF : covariance of EDA soil moisture and T2m (at each obs time and for each layer)
!     COVAR_SSM_CV_SEKF : covariance of EDA soil moisture and surface soil moisture (for each time slot and for each layer)
!     VAR_CV_SEKF    : variance of EDA soil moisture 

!     Output from the EKF soil moisture analysis
!     ------------------------------------------
!     NSEKFQ         : Quality flag
!     VSM_ANAL_INC   : Soil moisture analysis increment
!                      (local) 1D (N_SEKF_CV)
!     SEKF_INC       : Soil moisture analysis increments at t0
!                      (global) 3D (NPROMA,NSLAY,1)
!     SEKF_F         : Soil moisture analyses at analyses times
!                      (global) 3D (NPROMA,NSLAY*NANAL,1)
!     VSM_INC_L1     : Increment layer 1 (NPROMA, NBLOCK)
!     VSM_INC_L2     : Increment layer 2 (NPROMA, NBLOCK)
!     VSM_INC_L3     : Increment layer 3 (NPROMA, NBLOCK)
!     SEKF_Q         : SEKF quality flag (NGPTOT,1,1)
!     SEKF_DEBUG     : Debug fields from the local SEKF analysis (32)
!     SEKF_DEBUG_3D  : Debug fields global (NGPTOT,32,1)
!     SEKF_A         : Analysis error (NGPTOT, N_SEKF_CV*N_SEKF_CV,1)
!     SEKF_B         : Background error (NGPTOT, N_SEKF_CV*N_SEKF_CV,1)
!     SEKF_G         : Kalman Gain (NGPTOT, N_SEKF_CV*NOBS_SEKF,1)
!     SEKF_J         : Jacobians (NGPTOT, NOBS_SEKF*N_SEKF_CV,1)

!     Logicals
!     --------
!     LUSEKF_REF     : True if the SEKF soil moisture analysis is used
!                      (nconf = 302)
!     LUSE_T2M       : True if analysed 2 m temperature data are used for the SEKF
!     LUSE_RH2M      : True if analysed 2 m relativ humidity data are used for the SEKF
!     LUSE_ASCAT     : True if ASCAT soil moisture is used for the SEKF
!     LUSE_SMOS_TB   : True if SMOS brightness temperature is used for the SEKF
!     LUSE_SMOS_SM   : True if SMOS neural network soil moisture is used for the SEKF
!     LUSE_JATM      : True is J computed from full 3D perturbed run
!     LUSE_EDA_JACOB : True is H computed from EDA spread
!     LREAD_EDA_COV_SSA : True to use variance and covariance fields from EDA postprocessing

!     Author:
!     M. Drusch *ECMWF* 17/07/07
!     2008-06-12 - Klaus Scipal - extended namelist (error estimates) 
!     2009-02-03 - Patricia de Rosnay - Offline Jacobian computation
!     2009-10-12 - Patricia de Rosnay - Offline Jacobian completed
!     2011-01-10 - Joaquin Munoz Sabater - Introduce SMOS variables
!     2014-12-22 - Joaquin Munoz Sabater - 3D structure B error matrix
!     2017-12-19 - Patricia de Rosnay - EDA Jacobians
!     2018-04-14 - Patricia de Rosnay - EDA Jacobians, add the EDA post processing to archive var and covar needed for the EDA Jacobians
!     2018-05-   - Patricia de Rosnay - introduce SMOS soil moisture variables

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: N_OBSTIME_MAX=10
INTEGER(KIND=JPIM), PARAMETER :: N_ANALTIME_MAX=10
INTEGER(KIND=JPIM), PARAMETER :: N_SOILM_MAX=4

INTEGER(KIND=JPIM) :: N_SEKF_OUPT
INTEGER(KIND=JPIM) :: N_SEKF_CV
INTEGER(KIND=JPIM) :: N_SEKF_PT
INTEGER(KIND=JPIM) :: NSLAY

INTEGER(KIND=JPIM) :: NSSA_TSLOT
INTEGER(KIND=JPIM) :: NSMOS_SM_FIL
INTEGER(KIND=JPIM) :: N_EDAMEMBERS
INTEGER(KIND=JPIM), ALLOCATABLE :: NOBS_ASCAT(:,:), NOBS_SMOS_TB(:,:), NOBS_SMOS_SM(:,:)
INTEGER(KIND=JPIM) :: NOBS_ASCAT_TOT,NOBS_ASCAT_MAX
INTEGER(KIND=JPIM) :: NOBS_SEKF_MAX
INTEGER(KIND=JPIM) :: NOBS_SIM_EDA
INTEGER(KIND=JPIM) :: NGRIB_FIELDS_SCREEN
INTEGER(KIND=JPIM) :: NGRIB_FIELDS_SMOS
INTEGER(KIND=JPIM) :: NANAL

INTEGER(KIND=JPIM) :: NREF_BLOCK
INTEGER(KIND=JPIM) :: NREF_GP
INTEGER(KIND=JPIM) :: NSEKFQ

INTEGER(KIND=JPIM) :: NSCREENTIME(N_OBSTIME_MAX)
INTEGER(KIND=JPIM) :: NANATIME(N_ANALTIME_MAX)
INTEGER(KIND=JPIM) :: NSCREENDATE(N_OBSTIME_MAX)
INTEGER(KIND=JPIM) :: NANADATE(N_ANALTIME_MAX)
INTEGER(KIND=JPIM), ALLOCATABLE :: NOBS_TSTEP(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: NANA_TSTEP(:)

INTEGER(KIND=JPIM), ALLOCATABLE :: STYPE(:,:)

REAL(KIND=JPRB)    :: VSM_PERT_INC(N_SOILM_MAX)

REAL(KIND=JPRB)    :: MOD_ERR
REAL(KIND=JPRB)    :: BACK_ERR
REAL(KIND=JPRB)    :: BACK_ERR_L1
REAL(KIND=JPRB)    :: BACK_ERR_L2
REAL(KIND=JPRB)    :: BACK_ERR_L3
REAL(KIND=JPRB)    :: ANA_ERR
REAL(KIND=JPRB)    :: T2M_ERR
REAL(KIND=JPRB)    :: RH2M_ERR
REAL(KIND=JPRB)    :: ASCAT_ERR
REAL(KIND=JPRB)    :: SMOS_TB_COEF_ERR
REAL(KIND=JPRB)    :: SMOS_SM_COEF_ERR
REAL(KIND=JPRB)    :: SMOS_SM_MAX_RFI
REAL(KIND=JPRB)    :: EDAH_TAPER
REAL(KIND=JPRB)    :: SMOS_SM_MIN_VOL_ERR

! For coupled Jacobians: save trajectories
REAL(KIND=JPRB), ALLOCATABLE :: FKF_TENT(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: FKF_TENQ(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: FKF_TENU(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: FKF_TENV(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: FKF_SURF_SO(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: FKF_SURF_TH(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: FKF_SURF_CR(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: FKF_SURF_LR(:,:,:)

! Offline Jacobians: 1 time step, last dim on pert dim
! for vdfouter
REAL(KIND=JPRB), ALLOCATABLE :: SKF_ZWSA(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_ZTSA(:,:,:,:)

INTEGER(KIND=JPIM), ALLOCATABLE :: ISKF_IPBLTYPE(:,:,:)

REAL(KIND=JPRB), ALLOCATABLE :: SKF_ZAZ0M (:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_ZAZ0H(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PVDIS(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PVDISG(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_ZDISSGW(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PFTLHEV(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PFTLHSB(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PFWSB(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PUCFL(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PVCFL(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTCFL(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PDCFL(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PQCFL(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTENT_CUM(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTENQ_CUM(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_ZZINV(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PBLH(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PU10N(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PV10N(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PUST(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_ZFRSOTI(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_ZEVAPSNW(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_ZGUST(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTENT(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTENQ(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_ZTENA(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_ZTENI(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_ZTENL(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTENU(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTENV(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTLE1(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTL(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PUSTRTI(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PVSTRTI(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PAHFSTI(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PEVAPTI(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTSKTI(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PDIFTS(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PDIFTQ(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PDIFTL(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PDIFTI(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PSTRTU(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PSTRTV(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTOFDU(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTOFDV(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PSTRDU(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PSTRDV(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_ZKH_VDF(:,:,:,:)
!REAL(KIND=JPRB), ALLOCATABLE :: SKF_PDHTLS(:,:,:,:,:)
!REAL(KIND=JPRB), ALLOCATABLE :: SKF_PDHTSS(:,:,:,:,:)
!REAL(KIND=JPRB), ALLOCATABLE :: SKF_PDHTTS(:,:,:,:,:)
!REAL(KIND=JPRB), ALLOCATABLE :: SKF_PDHTIS(:,:,:,:,:)

! For surftstp
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PZO(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PHO(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PHO_INV(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PDO(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PADVT(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PADVS(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTRI0(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTRI1(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PSWDK_SAVE(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PUO0(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PVO0(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTO0(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PSO0(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PUOC(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PVOC(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PUSTRC(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PVSTRC(:,:,:)

REAL(KIND=JPRB), ALLOCATABLE :: SKF_PFTG12(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PFWROD(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PFWRO1(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PFWG12(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PFWMLT(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PFWEV(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PENES(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PDIFM(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PDIFT(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PDIFS(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PSNSE1(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTSNE1(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PASNE1(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PRSNE1(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTSAE1(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTIAE1(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PWLE1(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PWSAE1(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTLICEE1(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTLMNWE1(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTLWMLE1(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTLBOTE1(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTLSFE1(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PHLICEE1(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PHLMLE1(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PUOE1(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PVOE1(:,:,:,:)     
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PTOE1(:,:,:,:)    
REAL(KIND=JPRB), ALLOCATABLE :: SKF_PSOE1(:,:,:,:)   


! SEKF variables
REAL(KIND=JPRB), ALLOCATABLE :: VOL_SM_FG(:,:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: T2M_FG(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: TD2M_FG(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: TTB_FG(:,:,:,:,:,:)

REAL(KIND=JPRB), ALLOCATABLE :: T2M_OBS(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: TD2M_OBS(:,:,:)

REAL(KIND=JPRB), ALLOCATABLE :: SMOS_SM_OBS(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SMOS_SM_ERR(:,:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: SMOS_SM_TST(:,:,:)

REAL(KIND=JPRB), ALLOCATABLE :: COVAR_RH2M_CV_SEKF(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: COVAR_T2M_CV_SEKF(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: COVAR_SSM_CV_SEKF(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: VAR_CV_SEKF(:,:,:,:)

REAL(KIND=JPRB), ALLOCATABLE :: VSM_ANAL_INC(:)
REAL(KIND=JPRB), ALLOCATABLE :: SEKF_INC(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SEKF_F(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: VSM_INC_L1(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: VSM_INC_L2(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: VSM_INC_L3(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SEKF_Q(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SEKF_DEBUG_3D(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SEKF_DEBUG(:)

REAL(KIND=JPRB), ALLOCATABLE :: SEKF_A(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SEKF_B(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SEKF_G(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SEKF_J(:,:,:)

LOGICAL :: LUSEKF_REF
LOGICAL :: LUSE_T2M, LUSE_RH2M, LUSE_ASCAT, LUSE_SMOS_TB, LUSE_SMOS_SM
LOGICAL :: LUSE_JATM  ! switch for offline jacobians
LOGICAL :: PERT_CHESS
LOGICAL :: LUSE_EDA_JACOB ! use the EDA to compute the SEKF Jacobians
LOGICAL :: LREAD_EDA_COV_SSA ! use the EDA post processed fields instead of raw fields
LOGICAL, ALLOCATABLE  :: LSMOS_SM_OBS(:,:,:)


TYPE TYPE_ASCAT_OBS
  REAL(KIND=JPRB)     :: ZOBS
  INTEGER(KIND=JPIM)  :: ISTEP, IPROMA, IGPBLCK, IBODY, IHDR, IGP
END TYPE TYPE_ASCAT_OBS

TYPE(TYPE_ASCAT_OBS), ALLOCATABLE :: YLASCATOBS(:)

TYPE TYPE_SMOS_OBS
  REAL(KIND=JPRB)     :: ZOBS
  REAL(KIND=JPRB)     :: ZOBSBC
  INTEGER(KIND=JPIM)  :: ISTEP
  INTEGER(KIND=JPIM)  :: IPOL
  REAL(KIND=JPRB)     :: ZANG
  REAL(KIND=JPRB)     :: ZOERR
END TYPE TYPE_SMOS_OBS

TYPE(TYPE_SMOS_OBS), ALLOCATABLE :: YLSMOSOBS(:,:,:)


!     ------------------------------------------------------------------
END MODULE YOMSEKF
