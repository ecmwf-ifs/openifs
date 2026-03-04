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

MODULE YOMPARAR

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE TPARAR
!*
!     ------------------------------------------------------------------

!     VARIABLES pour utiliser la PHYSIQUE de meso_NH :
!     VARIABLES to use the MESO-NH physics:

CHARACTER(LEN=4)   :: CMICRO      ! Microphysics scheme ('OLD3', 'OLD4', 'ICE3',
                                  ! 'ICE4' or 'LIMA')
CHARACTER(LEN=4)   :: CSEDIM      ! Microphysics sedimentation scheme
                                  ! (SPLI or STAT)
INTEGER(KIND=JPIM) :: NSPLITR     ! Time splitting for Eulerian sedimentation
INTEGER(KIND=JPIM) :: NSPLITG     ! Time splitting for Eulerian sedimentation
INTEGER(KIND=JPIM) :: NRR, NRRL, NRRI   !number of microphysical species
INTEGER(KIND=JPIM) :: NSV         !number of passiv variables in MesoNH,
                                  ! always 0 in AROME
INTEGER(KIND=JPIM) :: NSWB_MNH    !number of SW bands for surface
                                  ! (must be equal to NSW !!)
INTEGER(KIND=JPIM) :: NGPAR       !number of fields in the buffer containing
                                  ! the 2D pseudo-historical variables.
INTEGER(KIND=JPIM) :: MINPRR      !pointer on INPRR 
INTEGER(KIND=JPIM) :: MACPRR      !pointer on ACPRR 
INTEGER(KIND=JPIM) :: MINPRS      !pointer on INPRS
INTEGER(KIND=JPIM) :: MACPRS      !pointer on ACPRS
INTEGER(KIND=JPIM) :: MINPRG      !pointer on INPRG
INTEGER(KIND=JPIM) :: MACPRG      !pointer on ACPRG
INTEGER(KIND=JPIM) :: MALBDIR     !pointer on ALBDIR
INTEGER(KIND=JPIM) :: MALBSCA     !pointer on ALBSCA
INTEGER(KIND=JPIM) :: MRAIN       !pointer on surface rain
INTEGER(KIND=JPIM) :: MSNOW       !pointer on surface snow
INTEGER(KIND=JPIM) :: MGZ0        !pointer on GZ0 
INTEGER(KIND=JPIM) :: MGZ0H       !pointer on GZ0H 
INTEGER(KIND=JPIM) :: MVQS        !pointer on surface moisture
INTEGER(KIND=JPIM) :: MVTS        !pointer on surface temperature 
INTEGER(KIND=JPIM) :: MVEMIS      !pointer on surface emissivity
INTEGER(KIND=JPIM) :: MSWDIR      !pointer on SW direct surface flux
INTEGER(KIND=JPIM) :: MSWDIF      !pointer on SW surface diffuse flux
INTEGER(KIND=JPIM) :: MFRTHDS     !pointer on IR downward surface flux
INTEGER(KIND=JPIM) :: MLSM        !pointer on land-sea mask
INTEGER(KIND=JPIM) :: MCD         !pointer on drag coefficient

REAL(KIND=JPRB), DIMENSION(:), ALLOCATABLE  :: XSW_BANDS  !SW spectral bands
! for ext. surface scheme
LOGICAL :: LOSUBG_COND ! switch to activate subgrid condensation
LOGICAL :: LOSEDIC     ! activate cloud sedimentation
LOGICAL :: LOWARM      ! see OWARM in mesoNH
LOGICAL :: LOSIGMAS    ! activate calculation of variance of departure to 
                       ! saturation in turb scheme (to be used in subgrid condensation)
LOGICAL :: LOLSMC      ! Land/sea mask for cloud droplet number conc. 
LOGICAL :: LOTOWNC     ! Town mask for cloud droplet number conc.

LOGICAL :: LOCND2      ! Separate solid and liquid phase
LOGICAL :: LGRSN       ! Turn graupel to snow for high supersaturation wrt ice 
REAL(KIND=JPRB) :: RADGR  ! Tuning of ice for radiation, TO BE REMOVED
REAL(KIND=JPRB) :: RADSN  ! Tuning of ice for radiation, TO BE REMOVED

REAL(KIND=JPRB) :: VSIGQSAT ! coeff applied to qsat variance contribution
                            ! for subgrid condensation

! switches for MF scheme (Pergaud et al)
CHARACTER (LEN=4)  :: CMF_UPDRAFT  ! Type of Mass Flux Scheme 
                                     ! 'NONE','DUAL', 'EDKF', 'RHCJ' or 'RAHA' 
CHARACTER (LEN=4)  :: CMF_CLOUD  ! type of cloud scheme associated with MF Scheme
                                     ! 'NONE', 'DIRE' or 'STAT'
LOGICAL            :: LMIXUV    ! True if mixing of momentum

LOGICAL            :: LLCRIT    ! True if temperature dependent 
                                ! critical condensation in EDMFm 
LOGICAL            :: LTOTPREC  ! True if precipitation tendencies
                                ! from the sub-grid scheme are
                                ! added to the total precip tendencies.

! Tuning variables for MF scheme 

REAL(KIND=JPRB) :: XALP_PERT    ! coefficient for the perturbation of
                                ! theta_l and r_t at the first level of 
                                ! the updraft
REAL(KIND=JPRB) :: XABUO        ! coefficient of the buoyancy term in the w_up equation  
REAL(KIND=JPRB) :: XBENTR       ! coefficient of the entrainment term in the w_up equation
REAL(KIND=JPRB) :: XBDETR       ! coefficient of the detrainment term in the w_up equation
REAL(KIND=JPRB) :: XCMF         ! coefficient for the mass flux at the first level 
                                ! of the updraft (closure) 
REAL(KIND=JPRB) :: XENTR_MF     ! entrainment constant (m/Pa) = 0.2 (m)
REAL(KIND=JPRB) :: XCRAD_MF     ! cloud radius in cloudy part
REAL(KIND=JPRB) :: XENTR_DRY    ! coefficient for entrainment in dry part 
REAL(KIND=JPRB) :: XDETR_DRY    ! coefficient for detrainment in dry part
REAL(KIND=JPRB) :: XDETR_LUP    ! coefficient for detrainment in dry part
REAL(KIND=JPRB) :: XKCF_MF      ! coefficient for cloud fraction
REAL(KIND=JPRB) :: XKRC_MF      ! coefficient for convective rc
REAL(KIND=JPRB) :: XTAUSIGMF    ! typical eddy turnover time for the STAT cloud scheme of EDMF
REAL(KIND=JPRB) :: XPRES_UV     ! coefficient for pressure term in wind mixing
REAL(KIND=JPRB) :: XFRAC_UP_MAX ! maximum Updraft fraction
REAL(KIND=JPRB) :: XALPHA_MF    ! coefficient for updraft fraction in STA2 cloud scheme
REAL(KIND=JPRB) :: XSIGMA_MF    ! coefficient for sigma in STA2 cloud scheme

! Tuning variables for RHCJ10 updraft :

REAL(KIND=JPRB) :: XA1  ! Tuning variables for RHCJ10 updraft 
REAL(KIND=JPRB) :: XB   ! Tuning variables for RHCJ10 updraft 
REAL(KIND=JPRB) :: XC   ! Tuning variables for RHCJ10 updraft 
REAL(KIND=JPRB) :: XBETA1 ! Tuning variables for RHCJ10 updraft 

!  Tuning parameter for Hourdin et al closure

REAL(KIND=JPRB) :: XR !  Tuning parameter for Hourdin et al closure

! Thermodynamic constant to compute thetas from thetal

REAL(KIND=JPRB) :: XLAMBDA  ! Thermodynamic constant to compute thetas from thetal
LOGICAL :: LTHETAS  ! TRUE to use Thetas, FALSE to use Thetal

! * for the squall line case:
LOGICAL :: LSQUALL ! use for the squall line case
INTEGER(KIND=JPIM) :: NREFROI1 !starting point for cooling (lsquall case) 
INTEGER(KIND=JPIM) :: NREFROI2 !end point for cooling (lsquall case)
REAL(KIND=JPRB) :: VSQUALL ! mean velocity displacement of the squall line (lsquall case)

! * for the MESO-NH physics printings:
INTEGER(KIND=JPIM) :: NPTP ! index in NPROMA paquet where the print will be done
INTEGER(KIND=JPIM) :: NPRINTFR !frequency of physical prints in apl_arome

!* for other diagnostics
!  wmax per vertical level 
LOGICAL :: LDIAGWMAX !activate print of WMAX in apl_arome
INTEGER(KIND=JPIM) :: NDIAGWMAX ! frequency of preceding prints (in time step)

!* for chemical scheme
! time step factor
INTEGER(KIND=JPIM) :: NDTCHEM ! time step factor for chemical scheme
!* for MNH budget anlysis
LOGICAL :: LAROBU_ENABLE ! for MNH budget anlysis
!* for turbulence scheme
REAL(KIND=JPRB) :: XLINI ! minimum bl89 mixing length
LOGICAL :: LHARATU !  if true RACMO turbulence is used
                   !  (yet only to be used in combination with EDMFm convection (DUAL))
!* Subgrid precipitation scheme
CHARACTER (LEN=4) :: CSUBG_AUCV_RC ! type of rc->rr autoconversion scheme
                                   ! 'CLFR', 'PDF' or 'NONE'
CHARACTER(LEN=80) :: CSUBG_RC_RR_ACCR ! type of rc rc autoconversion
                                   ! 'PRFR' or 'NONE'
CHARACTER(LEN=80) :: CSUBG_RR_EVAP ! type of evaporation scheme
                                   ! 'PRFR, 'CLFR' or 'NONE'
CHARACTER(LEN=80) :: CSUBG_PR_PDF  ! PDF chosen for precipitation production
                                   ! (NONE, SIGM, HLCRECTPD, HLCTRIANGPDF, HLCQUADRAPDF or HLCISOTRIPDF)

!* for autoconversion qi,qc
REAL(KIND=JPRB) :: RCRIAUTI ! ice autoconversion threshold
REAL(KIND=JPRB) :: RCRIAUTC ! cloud water autoconversion threshold
REAL(KIND=JPRB) :: RT0CRIAUTI ! temperature threshold for ice autoconversion
LOGICAL :: LCRIAUTI ! to activate use of namelists parameters for ice autoconversion
REAL(KIND=JPRB) :: XTSTEP_TS ! Approximative time step for microphysics time-splitting (0 for no time-splitting)
REAL(KIND=JPRB) :: XMRSTEP ! Maximum mixing ratio change before computing again microphysical tendencies
INTEGER(KIND=JPIM) :: NMAXITER_MICRO ! Maximum number of iterations for mixing ratio or time splitting
!* for the snow riming
CHARACTER (LEN=4) :: CSNOWRIMING ! choice of the snow riming parameterisation
                                 ! 'OLD' or 'M90' (Murakami 1990)
REAL(KIND=JPRB) :: XFRACM90 ! Fraction used for the Murakami 1990 formulation
!* for the graupel and the hail growth
LOGICAL :: LNULLWETG  ! When .TRUE. graupel wet growth is activated with null rate (to allow water shedding)
LOGICAL :: LWETGPOST  ! When .TRUE. graupel wet growth is activated with positive temperature (to allow water shedding)
LOGICAL :: LNULLWETH  ! Same as LNULLWETG but for hail
LOGICAL :: LWETHPOST  ! Same as LWETGPOST but for hail
!
LOGICAL :: LFEEDBACKT ! When .TRUE. feed back on temperature is taken into account
LOGICAL :: LEVLIMIT   ! When .TRUE. water vapour pressure is limited by saturation
LOGICAL :: LCONVHG ! TRUE to allow the conversion from hail to graupel
!
LOGICAL :: LCRFLIMIT !True to limit rain contact freezing to possible heat exchange
CHARACTER(LEN=1) :: CFRAC_ICE_ADJUST ! Rule to compute ice/liquid fraction in adjustments
CHARACTER(LEN=1) :: CFRAC_ICE_SHALLOW_MF ! Rule to compute ice/liquid fraction in shallow_mf
LOGICAL :: LSEDIM_AFTER !Sedimentation done after microphysics (.TRUE.) or before (.FALSE.)
!
REAL(KIND=JPRB) :: XSPLIT_MAXCFL ! Maximum CFL number allowed for SPLIT sedimentation scheme
!
!* For total cumulative 3D prec flux for MOCAGE
LOGICAL :: LFPREC3D ! Switch on total cumulative 3D prec flux output (for MOCAGE use)
!* For radiation :
REAL(KIND=JPRB) :: XCQVR ! reduction factor of water vapour used for radiation computation.
REAL(KIND=JPRB) :: GQVPLIM ! pressure value over which qv is damped towards 0 for radiation.
REAL(KIND=JPRB) :: GQVTOP ! qv value at the top of the atmopshere.
LOGICAL :: LQVTOP ! to activate modification of qv in input to radiation.

END TYPE TPARAR

!!TYPE(TPARAR), POINTER :: YRPARAR => NULL()


!     ------------------------------------------------------------------
END MODULE YOMPARAR
