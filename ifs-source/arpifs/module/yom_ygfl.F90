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

MODULE YOM_YGFL

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOE_AERODIAG, ONLY : NPAERAOT, NPAERLISI_VAR, NPAERLISI_WVL, NPAERLISI, &
  & TYPE_AERO_WVL_DIAG, NPAERO_WVL_DIAG

USE PAR_GFL

IMPLICIT NONE
SAVE

!-------------------------------------------------------------------------
! Contains the descriptors of GFL arrays
!-------------------------------------------------------------------------

TYPE TYPE_GFL_COMP ! Individual field descriptor

SEQUENCE ! daand: necessary to avoid memory corruption with gfortran 4.3.3

CHARACTER(LEN=16)  :: CNAME     = ''        ! ARPEGE field name 
INTEGER(KIND=JPIM) :: IGRBCODE  = -999      ! GRIB code
LOGICAL            :: LADV      = .FALSE.   ! Field advected or not
LOGICAL            :: LADV5     = .FALSE.   ! Field advected without wind increments
LOGICAL            :: LTDIABLIN = .FALSE.   ! Diabatic tendency is interpolated by lin. int.
LOGICAL            :: LHORTURB  = .FALSE.   ! Horizontal part affected by 3D turbulence
INTEGER(KIND=JPIM) :: NREQIN    = 0         ! 1 if field requiered in input, 0 if not, -1 if initialised
                                            ! with a reference value REFVALI
LOGICAL            :: LREQOUT   = .FALSE.   ! T if field requiered in output
LOGICAL            :: LGPINGP   = .TRUE.    ! GP field input as GP
LOGICAL            :: LMGRID    = .FALSE.   ! GP field to run on a different resolution grid
LOGICAL            :: LGP       = .FALSE.   ! Field exists and of grid-point type
LOGICAL            :: LSP       = .FALSE.   ! Field exists and of spectral type
LOGICAL            :: LCDERS    = .FALSE.   ! Derivatives required (spectral only)
LOGICAL            :: LACTIVE   = .FALSE.   ! Field in use
LOGICAL            :: LTHERMACT = .FALSE.   ! Field thermodynamically active
REAL(KIND=JPRB)    :: R         = 0.0_JPRB
REAL(KIND=JPRB)    :: RCP       = 0.0_JPRB
LOGICAL            :: LT9       = .FALSE.   ! Field in t-dt GFL
LOGICAL            :: LT1       = .FALSE.   ! Field in t+dt GFL
LOGICAL            :: LT5       = .FALSE.   ! Field in trajectory GFL
LOGICAL            :: LPHY      = .FALSE.   ! Field in physics GFL
LOGICAL            :: LPT       = .FALSE.   ! Field in PC phy. tend. GFL (GFLPT)
LOGICAL            :: LTRAJIO   = .FALSE.   ! Field written to and from trajectory structure
LOGICAL            :: LDIAG     = .FALSE.   ! Field is "diagnostic" at t; e.g. cloud fraction 
LOGICAL            :: LPC       = .FALSE.   ! Field in predictor/corrector time stepping (GFLPC)
REAL(KIND=JPRB)    :: REFVALI   = 0.0_JPRB  ! Reference value for init, used in case NREQIN==-1
! LAM specific attributes (Arome/Aladin)
LOGICAL            :: LADJUST0  = .FALSE.   ! True if field is thermodynamically adjusted at t
                                            ! (immediatly after inverse spectral transforms)
LOGICAL            :: LADJUST1  = .FALSE.   ! True if field is thermodynamically adjusted at t+dt
                                            ! (after SL interpolations and NL residuals)
INTEGER(KIND=JPIM) :: NCOUPLING = 0         ! 1 if field is coupled by Davies relaxation, 0 if not,
                                            ! -1 if coupled with reference value for coupling REFVALC
REAL(KIND=JPRB)    :: REFVALC   = 0.0_JPRB  ! Reference value for coupling, used in case NCOUPLING==-1
LOGICAL            :: LBIPER    = .FALSE.   ! True if field must be biperiodised inside the transforms
! End LAM specific attributes (Arome/Aladin)
CHARACTER(LEN=12)  :: CSLINT    = ''        ! S.L interpolaion "type"
INTEGER(KIND=JPIM) :: MP        = -99999999 ! Basic field "pointer"
INTEGER(KIND=JPIM) :: MPL       = -99999999 ! zonal derivative "pointer"
INTEGER(KIND=JPIM) :: MPM       = -99999999 ! Meridional derivative "pointer"
INTEGER(KIND=JPIM) :: MP9       = -99999999 ! Basic field "pointer" t-dt
INTEGER(KIND=JPIM) :: MP9_PH    = -99999999 ! Basic field "pointer" for Physics
INTEGER(KIND=JPIM) :: MP1       = -99999999 ! Basic field "pointer" t+dt
INTEGER(KIND=JPIM) :: MP5       = -99999999 ! Basic field "pointer" trajectory
INTEGER(KIND=JPIM) :: MP5L      = -99999999 ! zonal derivative "pointer" trajectory
INTEGER(KIND=JPIM) :: MP5M      = -99999999 ! Meridional derivative "pointer" trajectory
INTEGER(KIND=JPIM) :: MPSLP     = -99999999 ! Basic field "pointer" physics
INTEGER(KIND=JPIM) :: MPSP      = -99999999 ! Basic field "pointer" spectral space
INTEGER(KIND=JPIM) :: MP_SPL    = -99999999 ! Basic field "pointer" spline interpolation
INTEGER(KIND=JPIM) :: MP_SL1    = -99999999 ! Basic field "pointer" in SLBUF1
INTEGER(KIND=JPIM) :: MP_SLX    = -99999999 ! Basic field "pointer" in SLBUF1 for CPG_PT
INTEGER(KIND=JPIM) :: MPPT      = -99999999 ! Physics tendency "pointer"
INTEGER(KIND=JPIM) :: MPPC      = -99999999 ! Predictor/corrector auxiliary array "pointer"

! daand: INTFLEX attributes
LOGICAL            :: LWATER                ! TRUE for water species
LOGICAL            :: LPRECIP               ! TRUE for precipitating water species
REAL(KIND=JPRB)    :: RLZER                 ! Latent heat change at 0K

! gems nl ext
! general
INTEGER(KIND=JPIM) :: NCOUPLO4              ! Coupled to CTM by OASIS4 interface (obsolete)
! chemistry-only
LOGICAL            :: LASSIM                ! use as Control Variable (either monitored or assimilated)
INTEGER(KIND=JPIM) :: IGRIBDV               ! GRIB code of deposition velocity
INTEGER(KIND=JPIM) :: IGRIBTC               ! GRIB code of Total Column
INTEGER(KIND=JPIM) :: IGRBFLXO              ! GRIB code of Surface Flux (output)
INTEGER(KIND=JPIM) :: IGRBWDFLX             ! GRIB code of wet deposition flux (output)
INTEGER(KIND=JPIM) :: IGRBDDFLX             ! GRIB code of dry deposition flux (output)
LOGICAL            :: LDIFF                 ! Diffusion  on
LOGICAL            :: LCONV                 ! Convection on
REAL(KIND=JPRB)    :: RMOLMASS              ! Molar Mass
REAL(KIND=JPRB)    :: REFOLD                ! Efolding decay time
REAL(KIND=JPRB)    :: HENRYA                ! Henry constant a
REAL(KIND=JPRB)    :: HENRYB                ! Henry constant b
! end gems nl ext

LOGICAL            :: LNEGFIX               ! Cut off negative values in sugridug an
LOGICAL            :: LCOMAD                ! COMAD interpolation for GFL component
REAL(KIND=JPRB)    :: WENO_ALPHA            ! Alpha param (known as p) controlling WENO interpolation
LOGICAL            :: LMASSFIX              ! Correct mass error of sl advection in gpmodel (if LTRCMFIX)
REAL(KIND=JPRB)    :: BETAMFBC              ! beta parameters for BC mass fixer (if LMASSFIX)
INTEGER(KIND=JPIM) :: IFID = -1             ! FIELD identifier (all, not just COMPO)
TYPE(TYPE_GFL_COMP),POINTER :: PREVIOUS=>NULL()     ! Pointer to previously def. field

END TYPE TYPE_GFL_COMP

TYPE TYPE_GFL_NAML ! Individual field descriptor for namelist input

SEQUENCE ! daand: necessary to avoid memory corruption with gfortran 4.3.3

CHARACTER(LEN=16)  :: CNAME     ! ARPEGE field name 
INTEGER(KIND=JPIM) :: IGRBCODE  ! GRIB code
INTEGER(KIND=JPIM) :: NREQIN    ! 1 if field required in input, 0 if not, -1 if initialised
                                ! with a reference value REFVALI
REAL(KIND=JPRB) :: REFVALI      ! Reference value for initialisation, used in case NREQIN==-1
LOGICAL :: LREQOUT   = .FALSE.  ! T if field requiered in output
LOGICAL :: LGPINGP   = .FALSE.  ! GP field input as GP
LOGICAL :: LMGRID    = .FALSE.  ! GP field to run on a different resolution grid
LOGICAL :: LGP       = .FALSE.  ! Field exists and of grid-point type
LOGICAL :: LSP       = .FALSE.  ! Field exists and of spectral type
LOGICAL :: LCDERS    = .FALSE.  ! Derivatives required (spectral only)
LOGICAL :: LT9       = .FALSE.  ! Field in t-dt GFL
LOGICAL :: LT1       = .FALSE.  ! Field in t+dt GFL
LOGICAL :: LT5       = .FALSE.  ! Field in trajectory GFL
LOGICAL :: LPHY      = .FALSE.  ! Field with physics tendencies GFL
LOGICAL :: LPT       = .FALSE.  ! Field in PC physics tendency GFLPT
LOGICAL :: LTRAJIO   = .FALSE.  ! Field written to and from trajectory structure
LOGICAL :: LDIAG     = .FALSE.  ! Field is "diagnostic" at t; e.g. cloud fraction 
LOGICAL :: LPC       = .FALSE.  ! Field in predictor/corrector time stepping GFLPC
LOGICAL :: LADV      = .FALSE.  ! Field advected or not
LOGICAL :: LADV5     = .FALSE.  ! Field advected without wind increments
LOGICAL :: LINTLIN   = .FALSE.  ! Linear interpolation for field
LOGICAL :: LTDIABLIN = .FALSE.  ! Diabatic tendency is interpolated by linear int.
LOGICAL :: LHORTURB  = .FALSE.  ! Horizontal part affected by 3D turbulence
LOGICAL :: LQM       = .FALSE.  ! quasi-monotonous interpolation for field
LOGICAL :: LQMH      = .FALSE.  ! quasi-monotonous interpolation in horizontal for field
LOGICAL :: LQM3D     = .FALSE.  ! quasi-monotone interpolation applied directly in 3 dimensions
LOGICAL :: LQML3D    = .FALSE.  ! quasi-monotone cubi-linear hybrid interpolation
LOGICAL :: LVWENO    = .FALSE.  ! vertical quintic interpolation of WENO family
REAL(KIND=JPRB) :: WENO_ALPHA = 0.0_JPRB   ! alpha parameter defining the WENO interpolation
LOGICAL :: LSLHD     = .FALSE.  ! Semi-lagrangian horizontal diffusion used for field
LOGICAL :: LCOMAD    = .FALSE.  ! COMAD weights used for SL interpolation of field
LOGICAL :: LHV       = .FALSE.  ! Hermite vertical interpolation used for field (only ozone sofar)
LOGICAL :: LVSPLIP   = .FALSE.  ! vertical spline interpolation used for field (only ozone sofar)
INTEGER(KIND=JPIM) :: NCOUPLING ! 1 if field is coupled by Davies relaxation, 0 if not,
                                ! -1 if coupled with reference value for coupling REFVALC
REAL(KIND=JPRB) :: REFVALC      ! Reference value for coupling, used in case 
                                ! NCOUPLING==-1
! gems nl ext
! general
INTEGER(KIND=JPIM)  :: NCOUPLO4              ! Coupled to CTM by OASIS4 intefrace =1 input,=2 in&output,=-1 none (obsolete)
! chemistry-only
LOGICAL             :: LASSIM  = .FALSE.     ! use as Control Variable (either monitored or assimilated)
INTEGER(KIND=JPIM)  :: IGRIBDV               ! GRIB code of deposition velocity 
INTEGER(KIND=JPIM)  :: IGRIBTC               ! GRIB code of Total Column
INTEGER(KIND=JPIM)  :: IGRBFLXO              ! GRIB code of Surface Flux (output)
INTEGER(KIND=JPIM)  :: IGRBWDFLX             ! GRIB code of wet deposition flux (output)
INTEGER(KIND=JPIM)  :: IGRBDDFLX             ! GRIB code of dry deposition flux (output)
LOGICAL             :: LDIFF   = .FALSE.     ! Diffusion  on
LOGICAL             :: LCONV   = .FALSE.     ! Convection on
LOGICAL             :: LNEGFIX = .FALSE.     ! Cut off negative values in sugridug and callpar
LOGICAL             :: LMASSFIX= .FALSE.     ! Correct mass error of sl advection in gpmodel (if LTRCMFIX)
REAL(KIND=JPRB)     :: BETAMFBC              ! beta parameters for BC mass fixer (if LMASSFIX)
REAL(KIND=JPRB)     :: RMOLMASS              ! Molar Mass 
REAL(KIND=JPRB)     :: REFOLD                ! Efolding  decay time 
REAL(KIND=JPRB)     :: HENRYA                ! Henry constant a 
REAL(KIND=JPRB)     :: HENRYB                ! Henry constant b 
! end gems nl ext

END TYPE TYPE_GFL_NAML

!-------------------------------------------------------------------------
! Derived types for describing the GFL structure.
!-------------------------------------------------------------------------
! Modifications:
! 03/07/09 C. Fischer - add Arome/Aladin attributes
! 03/10/01 C. Moussy  - add Arome/Aladin attributes coupling
! 03/10/31 M. Tudor   - add physics tendencies for predictor-corrector
! 05/10/10 J. Haseler - switch for I/O to trajectory structure
! 2004-Nov F. Vana    - update of CSLINT attribute
! 20-Feb-2005 Vivoda  - 3TL Eul PC scheme (GFLPC)
! 07/06/27 E. Holm    - TL/AD advection without wind increments LADV5
! 12/04/08 J. Flemming - GFL attribute extention for GEMS 
! 22-Feb-11 F. Vana   - LTDIABLIN and LHORTURB
! spring 2011 ECMWF   - LINTLIN
! Nov. 2013           - LCOMAD
! 2013-11, D. Degrauwe - INTFLEX attributes
! 06-Feb-2015 F. Vana & M. Kharoutdinov - super-parametrization fields.
! 02-Feb-2017 S.-J. Lock - LSPPTGFL, switch for SPPT perturbations to chemistry/aerosol tendencies
! 2017-11-11 M Ahlgrimm - add cloud heterogeneity FSD
! 14-Feb-2019 P Bechtold- add turbulence diagnostics EDR Parameter
! 20-Feb-2019 F. Vana  - WENO quintic interpolation
! 11-Sep-2020 F. Vana  - SLAVEPP in TL/AD

TYPE TYPE_GFLD

SEQUENCE ! daand: necessary to avoid memory corruption with gfortran 4.3.3

! Overall descriptor,dimensioning etc.
INTEGER(KIND=JPIM) :: NUMFLDS     = 0  ! Number of GFL fields
INTEGER(KIND=JPIM) :: NDERS       = 0  ! Number of horizontal derivatives fields
INTEGER(KIND=JPIM) :: NUMSPFLDS   = 0  ! Number of spectrally represented GFL fields
INTEGER(KIND=JPIM) :: NUMGPFLDS   = 0  ! Number of grid-point GFL fields
INTEGER(KIND=JPIM) :: NUMFLDS9    = 0  ! Number of GFL fields in (t-dt) part
INTEGER(KIND=JPIM) :: NUMFLDS1    = 0  ! Number of GFL fields in (t+dt) array
INTEGER(KIND=JPIM) :: NUMSPFLDS1  = 0  ! Number of spectrally represented GFL fields (t+dt)
INTEGER(KIND=JPIM) :: NUMFLDS5    = 0  ! Number of GFL fields (trajectory)
INTEGER(KIND=JPIM) :: NUMFLDSPHY  = 0  ! Number of GFL fields (phys.)
INTEGER(KIND=JPIM) :: NUMFLDS_SPL = 0  ! Number of GFL fields (S.L. spline interpolation)
INTEGER(KIND=JPIM) :: NUMFLDS_SL1 = 0  ! Number of GFL fields in S.L. buffer 1
INTEGER(KIND=JPIM) :: NUMFLDSPC   = 0  ! Number of GFL fields (predictor/corrector)
INTEGER(KIND=JPIM) :: NDIM        = 0  ! Dimension of main array holding GFL fields(GFL)
INTEGER(KIND=JPIM) :: NUMFLDSPT   = 0  ! Number of GFL fields (phy. tend.)
INTEGER(KIND=JPIM) :: NDIM0       = 0  ! Dimension of t0 part of GFL
INTEGER(KIND=JPIM) :: NDIM9       = 0  ! Dimension of t-dt part of GFL
INTEGER(KIND=JPIM) :: NDIM1       = 0  ! Dimension of t+dt array (GFLT1)
INTEGER(KIND=JPIM) :: NDIM5       = 0  ! Dimension of traj. GFL array (GFL5)
INTEGER(KIND=JPIM) :: NDIMSLP     = 0  ! Diminsion of S.L. phys. GFL/GFL5 array (GFLSLP)
INTEGER(KIND=JPIM) :: NDIM_SPL    = 0  ! Dim. of arrays holding GFL fields (S.L.spline int.)
INTEGER(KIND=JPIM) :: NDIMPT      = 0  ! Dimension of phy. tend. GFL array (GFLPT)
INTEGER(KIND=JPIM) :: NDIMPC      = 0  ! Dimension of iterative scheme auxiliary array (GFLPC)

INTEGER(KIND=JPIM) :: NGFL_EXT
INTEGER(KIND=JPIM) :: NGFL_FORC
INTEGER(KIND=JPIM) :: NGFL_EZDIAG
INTEGER(KIND=JPIM) :: NGHG
INTEGER(KIND=JPIM) :: NGHG_ASSIM
INTEGER(KIND=JPIM) :: NAERO
INTEGER(KIND=JPIM) :: NACTAERO
INTEGER(KIND=JPIM) :: NDDHAERO
INTEGER(KIND=JPIM) :: NERA40
INTEGER(KIND=JPIM) :: NNOGW
INTEGER(KIND=JPIM) :: NEDRP
INTEGER(KIND=JPIM) :: NAEROUT
INTEGER(KIND=JPIM) :: NAEROCLIM
INTEGER(KIND=JPIM) :: NUVP
INTEGER(KIND=JPIM) :: NSLDIA
INTEGER(KIND=JPIM) :: NSLDIAGP
INTEGER(KIND=JPIM) :: NFSD

LOGICAL :: LAERCHEM    = .FALSE.   ! Couple sulphur cycle between chemistry and aerosol
LOGICAL :: LTRCMFIX    = .FALSE.   ! Activates tracer mass fixer
LOGICAL :: LTRCMFIX_PS = .FALSE.   ! Adjust pressure to conserve dry mass in mass fixer calculations
LOGICAL :: LTRCMFQM    = .FALSE.   ! Ensure global conservation is fully restored after limiter is used
                                   ! For NWP tracers limiter effect is negligible and option not essential 
LOGICAL :: LAERAOT     = .FALSE.
LOGICAL :: LAERLISI    = .FALSE.
LOGICAL :: LAEROUT     = .FALSE.
LOGICAL :: LAEROCLIM   = .FALSE.
LOGICAL :: LUVPOUT     = .FALSE.
LOGICAL :: LSPPTGFL    = .FALSE.   ! Activates SPPT perturbations to CHEM tendencies (C-IFS ensemble)

INTEGER(KIND=JPIM) :: NGEMS   ! The total number of "GEMS" fields.
INTEGER(KIND=JPIM) :: NCHEM
INTEGER(KIND=JPIM) :: NCHEM_ASSIM
INTEGER(KIND=JPIM) :: NCHEM_FLXO 
INTEGER(KIND=JPIM) :: NCHEM_WDFLX 
INTEGER(KIND=JPIM) :: NCHEM_DDFLX 
INTEGER(KIND=JPIM) :: NCHEM_DV
INTEGER(KIND=JPIM) :: NCHEM_TC
INTEGER(KIND=JPIM) :: NCHEM_SCV
INTEGER(KIND=JPIM) :: NEMIS2D ! Number of 2D emission fields for composition
INTEGER(KIND=JPIM) :: NEMIS2DAUX ! Number of 2D emission auxiliary fields for composition
INTEGER(KIND=JPIM) :: NEMIS3D ! Number of 3D emission fields for composition
INTEGER(KIND=JPIM) :: NLIMA  ! Total number of LIMA fields
INTEGER(KIND=JPIM) :: NAERO_WVL_DIAG
INTEGER(KIND=JPIM) :: NAERO_WVL_DIAG_TYPES

!     ------------------------------------------------------------------
!      Mass fixers
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: NNEGAFIX     ! Num of fields to apply -ve fixer
INTEGER(KIND=JPIM) :: NOPTNEGFIX   ! 1: simple negative fixer (reset to 0)
                                   ! 2: reset to local minimum

LOGICAL :: LQM3DCONS    = .FALSE.  ! Bermejo & Staniforth quasi-monotone limiter with improved
                                   ! conservation option. When true, applied to all GFL s.t. LQM3D=true
LOGICAL :: LADVNEGFIX   = .FALSE.  ! Activates negative fixer for advection
LOGICAL :: LTRCMFBC     = .FALSE.  ! Activate Bermejo & Conde mass fixer if true
LOGICAL :: LTRCMFPR     = .FALSE.  ! Activate Priestley mass fixer if true
LOGICAL :: LTRCMFMG     = .FALSE.  ! Activate Mac Gregor's mass fixer if true
LOGICAL :: LTRCMFP      = .FALSE.  ! Activate proportional mass fixer  
LOGICAL :: LTRCMFA_DIF  = .FALSE.  ! Activate additive mass fixer based on t - t+1 diff 
LOGICAL :: LTRCMFA_LAP  = .FALSE.  ! Activate additive mass fixer based on horizontal laplace 
LOGICAL :: LTRCMFA_VER  = .FALSE.  ! Activate additive mass fixer based on vertical laplace

LOGICAL :: LEXTRADF     = .FALSE.  ! Extra diagnostics 


INTEGER(KIND=JPIM) :: NFLDSFIX     ! Number of fields to be fixed
INTEGER(KIND=JPIM) :: NOPTMFBC     ! Bermejo & Conde fixer version
INTEGER(KIND=JPIM) :: NOPTMFPR     ! Priestley fixer version
INTEGER(KIND=JPIM) :: NOPTVFE      ! Use Vertical FE in calculation of column mass total
INTEGER(KIND=JPIM) :: NMFDIAGLEV   ! Determines global diagnostic output level for fixer:
                                   ! 0 - nothing, 1 - norms printed, 2 - norms + monotonicity
INTEGER(KIND=JPIM) :: NMFIXFLDS(JPNAMED_GFL+JPGHG+JPCHEM+JPAERO) 
                                   ! Index of fields to be corrected by mass fixers
INTEGER(KIND=JPIM) :: NNEGFLDS(JPNAMED_GFL+JPGHG+JPCHEM+JPAERO)  
                                   ! Index of fields to be corrected by SL -ve fixer
REAL(KIND=JPRB)    :: ZMFIXEPS     ! Threshold for mass fixing scheme
INTEGER(KIND=JPIM) :: NFMG         ! Number of fields for mgrid simulations
INTEGER(KIND=JPIM) :: NMGFLDS(JPNAMED_GFL+JPGHG+JPCHEM+JPAERO)
                                   ! Index of fields for mgrid simulations


TYPE(TYPE_GFL_COMP) :: YCOMP(JPGFL)    ! General descriptor of all components

TYPE(TYPE_GFL_COMP),POINTER  :: YLASTCOMP   => NULL() ! General descriptor of last component

TYPE(TYPE_GFL_COMP),POINTER  :: YQ          => NULL() ! Specific humidity
TYPE(TYPE_GFL_COMP),POINTER  :: YI          => NULL() ! Ice water
TYPE(TYPE_GFL_COMP),POINTER  :: YL          => NULL() ! Liquid water
TYPE(TYPE_GFL_COMP),POINTER  :: YLCONV      => NULL() ! Liquid water (CONV. PART)
TYPE(TYPE_GFL_COMP),POINTER  :: YICONV      => NULL() ! Ice    water (CONV. PART)
TYPE(TYPE_GFL_COMP),POINTER  :: YRCONV      => NULL() ! Rain         (CONV. PART)
TYPE(TYPE_GFL_COMP),POINTER  :: YSCONV      => NULL() ! Snow         (CONV. PART)
TYPE(TYPE_GFL_COMP),POINTER  :: YIRAD       => NULL() ! Radiative cloud Ice water
TYPE(TYPE_GFL_COMP),POINTER  :: YLRAD       => NULL() ! Radiative cloud Liquid water
TYPE(TYPE_GFL_COMP),POINTER  :: YS          => NULL() ! Snow
TYPE(TYPE_GFL_COMP),POINTER  :: YR          => NULL() ! Rain
TYPE(TYPE_GFL_COMP),POINTER  :: YG          => NULL() ! Graupel
TYPE(TYPE_GFL_COMP),POINTER  :: YH          => NULL() ! Hail
TYPE(TYPE_GFL_COMP),POINTER  :: YTKE        => NULL() ! Turbulent Kinetic Energy
TYPE(TYPE_GFL_COMP),POINTER  :: YTTE        => NULL() ! Turbulent Total Energy
TYPE(TYPE_GFL_COMP),POINTER  :: YEFB1       => NULL() ! First variable EFB scheme
TYPE(TYPE_GFL_COMP),POINTER  :: YEFB2       => NULL() ! Second variable EFB scheme
TYPE(TYPE_GFL_COMP),POINTER  :: YEFB3       => NULL() ! Third variable EFB scheme
TYPE(TYPE_GFL_COMP),POINTER  :: YA          => NULL() ! Cloud fraction
TYPE(TYPE_GFL_COMP),POINTER  :: YO3         => NULL() ! Ozone
TYPE(TYPE_GFL_COMP),POINTER  :: YSRC        => NULL() ! Second-order flux for AROME s'rc'/2Sigma_s2 multiplied by Lambda_3
TYPE(TYPE_GFL_COMP),POINTER  :: YMXL        => NULL() ! Prognostic mixing length
TYPE(TYPE_GFL_COMP),POINTER  :: YSHTUR      => NULL() ! Shear source term for turbulence.
TYPE(TYPE_GFL_COMP),POINTER  :: YFQTUR      => NULL() ! Flux form source term for turbulence -moisture.
TYPE(TYPE_GFL_COMP),POINTER  :: YFSTUR      => NULL() ! Flux form source term for turbulence -enthalpy.
TYPE(TYPE_GFL_COMP),POINTER  :: YCPF        => NULL() ! Convective precipitation flux
TYPE(TYPE_GFL_COMP),POINTER  :: YSPF        => NULL() ! Stratiform precipitation flux
TYPE(TYPE_GFL_COMP),POINTER  :: YCVGQ       => NULL() ! Moisture Convergence for french physics
TYPE(TYPE_GFL_COMP),POINTER  :: YQVA        => NULL() ! total humidity variation
TYPE(TYPE_GFL_COMP),POINTER  :: YGHG(:)     => NULL() ! Greenhouse Gases
TYPE(TYPE_GFL_COMP),POINTER  :: YCHEM(:)    => NULL() ! Chemistry
TYPE(TYPE_GFL_COMP),POINTER  :: YAERO(:)    => NULL() ! Aerosols
TYPE(TYPE_GFL_COMP),POINTER  :: YLRCH4      => NULL() ! CH4 loss rate (instantaneous field)
TYPE(TYPE_GFL_COMP),POINTER  :: YEMIS3D(:)  => NULL() ! 3D emissions for atmospheric composition
TYPE(TYPE_GFL_COMP),POINTER  :: YFORC(:)    => NULL() ! large scale forcing
TYPE(TYPE_GFL_COMP),POINTER  :: YEZDIAG(:)  => NULL() ! easy diagnostics
TYPE(TYPE_GFL_COMP),POINTER  :: YERA40(:)   => NULL() ! ERA40 diagnostic fields
TYPE(TYPE_GFL_COMP),POINTER  :: YNOGW(:)    => NULL() ! NORO GWD SCHEME
TYPE(TYPE_GFL_COMP),POINTER  :: YEDRP(:)    => NULL() ! turbulence diagnostics EDR Parameter
TYPE(TYPE_GFL_COMP),POINTER  :: YSLDIA(:)   => NULL() ! SL dynamics diagnostics
TYPE(TYPE_GFL_COMP),POINTER  :: YAERAOT(:)  => NULL() ! Aerosol optical thicknesses
TYPE(TYPE_GFL_COMP),POINTER  :: YAERLISI(:) => NULL() ! Aerosol lidar simulator
TYPE(TYPE_GFL_COMP),POINTER  :: YAEROUT(:)  => NULL() ! Aerosol outputs
TYPE(TYPE_GFL_COMP),POINTER  :: YAEROCLIM(:) => NULL() ! Aerosol climatology
TYPE(TYPE_GFL_COMP),POINTER  :: YUVP(:)     => NULL() ! UV-processor output
TYPE(TYPE_GFL_COMP),POINTER  :: YPHYS(:)    => NULL() ! PHYS output
TYPE(TYPE_GFL_COMP),POINTER  :: YPHYCTY     => NULL() ! PHYS input for MassCTY

TYPE(TYPE_GFL_COMP),POINTER  :: YRSPEC      => NULL() ! Specific gas constant


TYPE(TYPE_GFL_COMP),POINTER  :: YSDSAT      => NULL() ! Standard Deviation of the
                                                      ! SATuration Depression (Sigma_s) 
TYPE(TYPE_GFL_COMP),POINTER  :: YCVV        => NULL() ! Convective Vertical Velocity
TYPE(TYPE_GFL_COMP),POINTER  :: YRKTH       => NULL() ! Rasch-Kristjansson H tendency
TYPE(TYPE_GFL_COMP),POINTER  :: YRKTQV      => NULL() ! Rasch-Kristjansson Qv tendency
TYPE(TYPE_GFL_COMP),POINTER  :: YRKTQC      => NULL() ! Rasch-Kristjansson Qc tendency

! Prognostic convection variables: add 6 named components
TYPE(TYPE_GFL_COMP),POINTER  :: YUOM        => NULL() ! Updraught vert velocity
TYPE(TYPE_GFL_COMP),POINTER  :: YUAL        => NULL() ! Updraught mesh fraction
TYPE(TYPE_GFL_COMP),POINTER  :: YDOM        => NULL() ! Downdraught vert velocity
TYPE(TYPE_GFL_COMP),POINTER  :: YDAL        => NULL() ! Downdraught mesh fraction
TYPE(TYPE_GFL_COMP),POINTER  :: YUEN        => NULL() ! Updraught entrainment
TYPE(TYPE_GFL_COMP),POINTER  :: YUNEBH      => NULL() ! pseudo-historic convective

TYPE(TYPE_GFL_COMP),POINTER  :: YLIMA(:)    => NULL() ! LIMA prognostic fields

! cloud heterogeneity FSD
TYPE(TYPE_GFL_COMP),POINTER  :: YFSD        => NULL() ! PHYS output

! Extra fields

TYPE(TYPE_GFL_COMP),POINTER  :: YEXT(:)     => NULL() ! Extra fields

TYPE(TYPE_GFL_NAML)  :: YQ_NL                 ! Specific humidity
TYPE(TYPE_GFL_NAML)  :: YI_NL                 ! Ice water
TYPE(TYPE_GFL_NAML)  :: YL_NL                 ! Liquid water
TYPE(TYPE_GFL_NAML)  :: YLCONV_NL             ! Liquid water (CONV. PART)
TYPE(TYPE_GFL_NAML)  :: YICONV_NL             ! Ice    water (CONV. PART)
TYPE(TYPE_GFL_NAML)  :: YRCONV_NL             ! Rain         (CONV. PART)
TYPE(TYPE_GFL_NAML)  :: YSCONV_NL             ! Snow         (CONV. PART)
TYPE(TYPE_GFL_NAML)  :: YIRAD_NL              ! Radiative cloud Ice water
TYPE(TYPE_GFL_NAML)  :: YLRAD_NL              ! Radiative cloud Liquid water
TYPE(TYPE_GFL_NAML)  :: YS_NL                 ! Snow
TYPE(TYPE_GFL_NAML)  :: YR_NL                 ! Rain
TYPE(TYPE_GFL_NAML)  :: YG_NL                 ! Graupels
TYPE(TYPE_GFL_NAML)  :: YH_NL                 ! Hail
TYPE(TYPE_GFL_NAML)  :: YTKE_NL               ! Turbulent Kinetic Energy
TYPE(TYPE_GFL_NAML)  :: YTTE_NL               ! Turbulent Total Energy
TYPE(TYPE_GFL_NAML)  :: YEFB1_NL              ! First variable EFB scheme
TYPE(TYPE_GFL_NAML)  :: YEFB2_NL              ! Second variable EFB scheme
TYPE(TYPE_GFL_NAML)  :: YEFB3_NL              ! Third variable EFB scheme
TYPE(TYPE_GFL_NAML)  :: YA_NL                 ! Cloud fraction
TYPE(TYPE_GFL_NAML)  :: YO3_NL                ! Ozone
TYPE(TYPE_GFL_NAML)  :: YSRC_NL               ! Second-order flux for AROME
                                              ! s'rc'/2Sigma_s2
                                              ! multiplied by Lambda_3
TYPE(TYPE_GFL_NAML)  :: YMXL_NL               ! Prognostic mixing length
TYPE(TYPE_GFL_NAML)  :: YSHTUR_NL             ! Shear source term for turbulence.
TYPE(TYPE_GFL_NAML)  :: YFQTUR_NL             ! Flux form source term for turbulence -moisture.
TYPE(TYPE_GFL_NAML)  :: YFSTUR_NL             ! Flux form source term for turbulence -enthalpy.
TYPE(TYPE_GFL_NAML)  :: YCPF_NL               ! Convective precipitation flux
TYPE(TYPE_GFL_NAML)  :: YSPF_NL               ! Stratiform precipitation flux
TYPE(TYPE_GFL_NAML)  :: YCVGQ_NL              ! Moisture Convergence for french physics
TYPE(TYPE_GFL_NAML)  :: YQVA_NL               ! Total humidity variation

TYPE(TYPE_GFL_NAML)  :: YGHG_NL(JPGHG)        ! Greenhouse Gases
TYPE(TYPE_GFL_NAML)  :: YCHEM_NL(JPCHEM)      ! Chemical species
TYPE(TYPE_GFL_NAML)  :: YAERO_NL(JPAERO)      ! Aerosol fields
TYPE(TYPE_GFL_NAML)  :: YERA40_NL(JPERA40)    ! ERA40 diagnostic fields
TYPE(TYPE_GFL_NAML)  :: YNOGW_NL(JPNOGW)      ! NORO GWD SCHEME
TYPE(TYPE_GFL_NAML)  :: YEDRP_NL(JPEDRP)      ! EDR Parameter turbulence diagnostics
TYPE(TYPE_GFL_NAML)  :: YSLDIA_NL(JPSLDIA)    ! SL dynamics diagnostics
TYPE(TYPE_GFL_NAML)  :: YLRCH4_NL             ! CH4 loss rate
TYPE(TYPE_GFL_NAML)  :: YAERAOT_NL(NPAERAOT)  ! Aerosol optical thicknesses
TYPE(TYPE_GFL_NAML)  :: YAERLISI_NL(NPAERLISI_WVL,NPAERLISI_VAR) ! Aerosol lidar simulator
TYPE(TYPE_GFL_NAML)  :: YAEROUT_NL(JPAEROUT)  ! Aerosol outputs
TYPE(TYPE_GFL_NAML)  :: YAEROCLIM_NL(JPAEROCLIM)  ! Aerosol climatology
TYPE(TYPE_GFL_NAML)  :: YUVP_NL(JPUVP)        ! UV-processor outputs
TYPE(TYPE_GFL_NAML)  :: YRKTH_NL              ! Rasch-Kristjansson H tendency
TYPE(TYPE_GFL_NAML)  :: YRKTQV_NL             ! Rasch-Kristjansson Qv tendency
TYPE(TYPE_GFL_NAML)  :: YRKTQC_NL             ! Rasch-Kristjansson Qc tendency
TYPE(TYPE_GFL_NAML)  :: YPHYCTY_NL            ! PHYS input for MassCTY
TYPE(TYPE_GFL_NAML)  :: YFSD_NL               ! cloud heterogeneity
TYPE(TYPE_GFL_NAML)  :: YRSPEC_NL             ! Specific humidity

! Extra fields
TYPE(TYPE_GFL_NAML)  :: YSDSAT_NL             ! Standard Deviation of the
                                              ! SATuration Depression (Sigma_s) 
TYPE(TYPE_GFL_NAML)  :: YCVV_NL               ! Convective Vertical Velocity
TYPE(TYPE_GFL_NAML)  :: YFORC_NL(JPFORC)      ! Forcing precursor
TYPE(TYPE_GFL_NAML)  :: YEZDIAG_NL(JPEZDIAG)  ! Easy diagnostics
TYPE(TYPE_GFL_NAML)  :: YEXT_NL(JPGFL-JPNAMED_GFL-JPGHG-JPFORC-JPEZDIAG-JPAERO-JPERA40-&
 &                              JPNOGW-JPSLDIA-JPCH4S-JPAEROUT-JPAEROCLIM-NPAERAOT-NPAERLISI-&
 &                              JPUVP-JPCHEM-JPLIMA-JPPHYCTY-JPFSD-JPEDRP) ! Extra fields

! Prognostic convection variables: 6 more namelist components
TYPE(TYPE_GFL_NAML)  :: YUOM_NL               ! Updraught vert velocity
TYPE(TYPE_GFL_NAML)  :: YUAL_NL               ! Updraught mesh fraction
TYPE(TYPE_GFL_NAML)  :: YDOM_NL               ! Downdraught vert velocity
TYPE(TYPE_GFL_NAML)  :: YDAL_NL               ! Downdraught mesh fraction
TYPE(TYPE_GFL_NAML)  :: YUEN_NL               ! Updraught entrainment
TYPE(TYPE_GFL_NAML)  :: YUNEBH_NL             ! Pseudi Hist Conv cloud fraction

TYPE(TYPE_AERO_WVL_DIAG) :: YAERO_WVL_DIAG_NL(NPAERO_WVL_DIAG) ! Per-wavelength aerosol optical diagnostics
TYPE(TYPE_GFL_NAML)  :: YLIMA_NL(JPLIMA)      ! LIMA prognostic fields

END TYPE TYPE_GFLD

! GFL general descriptor
TYPE(TYPE_GFLD), POINTER :: YGFL => NULL()

END MODULE YOM_YGFL
