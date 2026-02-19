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

MODULE YOMJG

USE PARKIND1                    , ONLY : JPIM     ,JPRB
USE YOMWAVELET                  , ONLY : TYPE_WAVELETJB_CONFIG, TYPE_WAVELETJB_DATA,     &
                                       & TYPE_WAVELETJB_VCOR_STRUCT, TYPE_WAVELETJB_GRID_STRUCT
USE YEMWAVELET                  , ONLY : TYPE_LAMWAVELETJB_CONFIG, TYPE_LAMWAVELETJB_INFO,   &
                                       & TYPE_LAMWAVELET_VCOR_STRUCT, &
                                       & TYPE_LAMWAVELET_GRID_STRUCT
USE YOMJBCHVAR                  , ONLY : TYPE_JBCHVAR
USE SPECTRAL_FIELDS_MOD         , ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)
USE SPECTRAL_VARIABLES_MOD      , ONLY : SPECTRAL_VARIABLES

IMPLICIT NONE

SAVE

PRIVATE
PUBLIC :: TYPE_JB_STRUCT, TYPE_JB_CONFIG, TYPE_JB_CONFIG2, &
        & TYPE_JBCHVAR, JB_STRUCT, JPAEROCV, &
        & TYPE_JB_DATA, TYPE_FCEMN_STRUCT, TYPE_VCORS_STRUCT, &
        & TYPE_SPJB_VARS_INFO_STRUCT, TYPE_WAVELETJB_CONFIG, &
        & TYPE_WAVELETJB_DATA, TYPE_WAVELETJB_VCOR_STRUCT, &
        & TYPE_WAVELETJB_GRID_STRUCT, TYPE_LAMWAVELETJB_CONFIG, &
        & TYPE_LAMWAVELETJB_INFO, TYPE_LAMWAVELET_VCOR_STRUCT, &
        & TYPE_LAMWAVELET_GRID_STRUCT

!*
!      /YOMJG/ - Jb RELEVANT PARAMETERS

!     BILL HECKLEY     *ECMWF*

!     --------- definition of Jb through namelist NAMJG --------------------

!     LCFCE  : .T. ===> Use horizontally constant background errors   T= YES
!     LSPFCE : .T. ===> background error normalization in spectral space
!     L3DBGERR:.T. ===> Use 3d (nonseperable) variances of background error
!     REDNMC   Rate of standard deviation reduction from the NMC covariances
!     REDQGLOLAM: Factor between sigmab(q) from GLObal model Arpege and the LAM Aladin  
!     LFACHRO: Increase climate variance in tropics using ad hoc factor
!     LCLMSFCE: Use climatalogical variances in error growth model
!     LRDQERR: T=Read background error for Q from file. F=Calculate it.
!     LSBLATLONG: T=errgrib on lat/lon grid, F=errgrib on Gaussian grid
!     NGEOCOR:   1 ===> (default) homogeneous correlations on work sphere
!                2 ===> homogeneous correls on geographical sphere using SPDICO
!                 OPTION REMOVED
!     NJBVERB  level of diagnostics printouts in Jb setups : (default is 1)
!              0=terse  1=normal  2=verbose
!     LJBENER  make Jb equal to a total energy scalar product
!     RAMENER  amount of energy vs. standard Jb (LJBENER option)
!     RSCENER  scaling of the total energy Jb (LJBENER option)
!     LCORCOSU enforce compactly supported correlations
!     LCORCOSU_VCOR compact support on vertical correlations (instead of covariances)
!     LJGNRSGP normalize forecast gridpoint fields (q and o3) in spectral space
!     L_JBVCOORD:  .T. ===> correlations in "Jb vertical coordinates"
!     LJB_UNIVARIATE:    .T. => Jb is univariate. No balance operators.
!     LJB_OMEGA_BALANCE: .T. => balance operator includes an omega equation
!     LJB_OMEGA_BALANCE_RELAX : .T. => Apply optional scale-dependent relaxation of omega balance
!     LJB_NONLINEAR_BALANCE: .T. => use the linearized nonlinear balance eqn.
!     LJB_NONLINEAR_BALANCE_RELAX : .T. => Apply optional scale-dependent relaxation of non linear balance
!     LJB_NONLINEAR_BALANCE_THR : .T. => Apply optional threshold on the background vorticity used to linearize non linear balance
!     LJB_NONLINEAR_CVHUM: .T. => use the nonlinearly normalized relative
!                                 humidity control variable
!                                 (linearized in the inner loops)
!     LJB_SUPSAT: .T. => humidity cv with supersaturation and q-T balance
!     LJB_EDAQERR: .T. => relative errors read in and used in humidity control variable 
!                         instead of normalization look-up table. Option under LJB_NONLINEAR_CVHUM=.T.
!     LJB_SEPQERR: .T. => standard modification of errors to prescribed global mean profile
!     LJBBALBETA: .T. => use the beta-plane balance for Aladin Jb
!     NPSMAX: truncation of the bi-Fourier f for the beta-plane in Aladin Jb
!     TMEANUVER: Vertical profiles of std dev for Aladin mean wind components (U,V)
!     ELC is a length-scale at which the omega and non-linear balances
!     are (optionally) relaxed to a zero value
!     THRVORTNL is the threshold value applied (optionnaly) on the background vorticity of non linear balance
!     LJB_TLBALSTAT: .T. => Use TL instead of full balance operators when 
!                           generating BG statistics
!     RJB_CSLEN: Compact support length scale when generating BG statistics
!     LBALSTRATO: .T. => use the balance operators in stratosphere
!     NLEVBAL0, NLEVBAL1: levels for switching off balance in the stratoshere
!     COEQTERMJB: Weight given to humidity term in total-energy Jb
!     --------------------  Jb Internal arrays ------------------------

!     FGSCNM3D: INVERSE HORIZONTAL PREDICTION ERROR CORRELATION MATRIX
!               IN SPECTRAL SPACE - FOR VARIABLES IN SPA3JB
!     FGSCNM2D: INVERSE HORIZONTAL PREDICTION ERROR CORRELATION MATRIX
!               IN SPECTRAL SPACE - FOR VARIABLES IN SPA2JB

!     FP2TPS  : P TO T,PS OPERATOR - NEEDED FOR CVAR2IN/AD

!     FCEMN3D : global mean profiles of bg errors for fields in SPA3JB
!     FCEMN2D : global mean profiles of bg errors for fields in SPA2JB
!     FCEMN   : global mean bg errors for fields not in JB
!                     (Used in the forecast-error calculation, only.)
!         %PS : Global rms values of background errors (lnsp)
!         %U  : Global rms values of background errors (zonal wind component)
!         %Z  : Global rms values of background errors (height)
!         %T  : Global rms values of background errors (temperature)

!     FCEIMN  : Global rms values of input background errors

!     VCORS   : Structure containing vertical correlation matrices
!     VCORS%EVALS : MATRIX OF EIGENVALUES OF VERTICAL CORRELATIONS
!     VCORS%EVECS : MATRIX OF EIGENVECTORS OF VERTICAL CORRELATIONS

!     FGMWNE : MASS-WEIGHTED NORMALIZED FC ERRORS FROM OI

!     STPS(:,:,:),SDIV(:,:,:),SO3(:,:,:),BFACT1(:),BFACT2(:) : coefficients
!        of the generalized statistical balance operator (J. Derber)
!     SHUM(:,:,:),BFACT(:) - extension of new JB for ALADIN
!     SCORI_R(:),SCORI_I(:) : bi-Fourier f for JB in ALADIN

!     N_SPJB_VARS     : number of model variables that have a
!                       spectral representation in SPA3JB/SPA2JB.
!     M_GRIBCODES     : Used for NAMELIST input of the GRIB codes that
!                       identify each variable in Jb. These codes are
!                       used to match Jb variable with model variables.
!     M_GRIBCODES_EP  : Used for NAMELIST input of the GRIB codes that
!                       identify the ensemble perturbation fields.
!     M_GRIBCODES_FCE : Used for NAMELIST input of the GRIB codes that
!                       identify the forecast error fields.
!     C_COR_STRINGS   : Used for NAMELIST input of they strings that
!                       identify each variable in Jb. These strings are
!                       used to match Jb variables with data in the ".cv"
!                       statistics file.
!     SPJB_VARS_INFO  : information about the fields in SPA3JB/SPA2JB.

!     FCE%BUF         : buffer for forecast error standard deviations
!     FCEREN%BUF      : buffer for optional fce normalisation factors
!     EP%BUF          : buffer for ensemble perturbations
!     JBVCOORD%BUF    : buffer to hold Jb vertical levels


!     --------------- Jb internal variables ----------------------------

!     LCTLFG   Keep FG in control space (i.e. do not use LSUBBG)
!     LSUBBG   .T. ===> CHAVAR subtracts background
!              .F. ===> CHAVAR assumes background already subtracted
!              This is an internal switch used to control the
!              operation of routines CVAR2, CVAR3, CVARBL, FLTRG and
!              their inverses and adjoints.

!     NOFEP :  Number of fields in ensemble perturbation file
!     NOFCEF : Number of fields on forecast error file
!     NOMSPEC : Number of variables spectral in model space

!     NSMAX_JB max meridian wave number for B stats.
!     NMSMAX_JB max zonal wave number for B stats.
!-------------- GEMS specific variables ------------------------
!     NAEROCV: Number of aerosol-related control variables
!     JPAEROCV: Max number of aerosol-related control variables

!-------------- Ensemble bg statistics calculation ----------------

!     LFLTBGVARENS  .T. Filter the ensemble variances
!     LFLTBGVARENS_STORE  .T. Store members or calculate variances on the fly
!     LFLTBGCALC_CRT      .T. Filter the variances on the fly with criterion
!     LFLTBGCALC_GAUSS    .T. Use Gaussian version of the criterion (general otherwise)
!     LFLTBGCALC    .T. Calculate objectives filters for ens variances
!     LOBJTRUNC     .T. Calculate objective truncation of the filter
!     NSMAXT        Calculations are performed on a grid NSMAX, while
!                   results are truncated to NSMAXT (to avoid aliasing)
!     NFLEVTR       Vertical geometry of the input truncation file

!-------------- Inflation of ensemble members ---------------------
!     LINFLATCALC   .T. Calculate inflation factors
!     LENSMEAN_CALC .T. Calculate ensemble mean

!-------------- Specific Humidity variances inflation --------------
!     LREDNMCQ
!     REDNMC_Q
!     LSBQ_STRATOLOW .T. Set stratospheric values to 1.E-8
!
! Modifications
! -------------
!     M. Fisher   7-March-2012 Moved DEALLOCATES out of DEALGES into module
!-------------------------------------------------------------------------------


INTEGER(KIND=JPIM), PARAMETER :: JPAEROCV=1 ! Change to 2 if running with dual control variable

TYPE TYPE_JB_CONFIG   ! Namelist variables
  LOGICAL :: LRDQERR
  LOGICAL :: LREDNMCQ
  LOGICAL :: L3DBGERR
  LOGICAL :: LCFCE
  LOGICAL :: LCTLFG
  LOGICAL :: LSBLATLONG
  LOGICAL :: LSPFCE
  LOGICAL :: LJBENER
  LOGICAL :: LFACHRO
  LOGICAL :: LCLMSFCE
  LOGICAL :: LCORCOSU
  LOGICAL :: LCORCOSU_VCOR
  LOGICAL :: L_JBVCOORD
  LOGICAL :: LJB_UNIVARIATE
  LOGICAL :: LJB_OMEGA_BALANCE
  LOGICAL :: LJB_OMEGA_BALANCE_RELAX
  LOGICAL :: LJB_NONLINEAR_BALANCE
  LOGICAL :: LJB_NONLINEAR_BALANCE_RELAX
  LOGICAL :: LJB_NONLINEAR_BALANCE_THR
  LOGICAL :: LJB_NONLINEAR_CVHUM
  LOGICAL :: LJB_SUPSAT
  LOGICAL :: LJB_EDAQERR
  LOGICAL :: LJB_SEPQERR
  LOGICAL :: LJBBALBETA
  LOGICAL :: LFLTBGVARENS
  LOGICAL :: LFLTBGVARENS_STORE
  LOGICAL :: LFLTBGCALC
  LOGICAL :: LFLTBGCALC_CRT
  LOGICAL :: LFLTBGCALC_GAUSS
  LOGICAL :: LJBSIGB_UV
  LOGICAL :: LGP_FLTBGCALC 
  LOGICAL :: LOBJTRUNC
  LOGICAL :: LJB_TLBALSTAT
  LOGICAL :: LBALSTRATO
  LOGICAL :: LSBQ_STRATOLOW
  LOGICAL :: LINFLATCALC
  LOGICAL :: LENSMEAN_CALC
  LOGICAL :: LJB_NETCDF_IO
  LOGICAL :: LCV_NETCDF_IO
  
  INTEGER(KIND=JPIM) :: NGEOCOR
  INTEGER(KIND=JPIM) :: NJBVERB
  INTEGER(KIND=JPIM) :: NPSMAX
  INTEGER(KIND=JPIM) :: NSMAXT
  INTEGER(KIND=JPIM) :: NFLEVTR
  INTEGER(KIND=JPIM) :: NLEVBAL0
  INTEGER(KIND=JPIM) :: NLEVBAL1
  INTEGER(KIND=JPIM) :: NLEV_REDINC
  INTEGER(KIND=JPIM) :: N_SPJB_VARS
  INTEGER(KIND=JPIM) :: NSMAX_JB
  INTEGER(KIND=JPIM) :: NMSMAX_JB
  
  REAL(KIND=JPRB) :: REDNMC
  REAL(KIND=JPRB) :: REDNMC_Q
  REAL(KIND=JPRB) :: REDQGLOLAM
  REAL(KIND=JPRB) :: RAMENER
  REAL(KIND=JPRB) :: RSCENER
  REAL(KIND=JPRB) :: RJB_CSLEN
  REAL(KIND=JPRB) :: ELC
  REAL(KIND=JPRB) :: THRVORTNL
  REAL(KIND=JPRB) :: COEQTERMJB

  TYPE(SPECTRAL_VARIABLES) :: SPVARS
END TYPE TYPE_JB_CONFIG

TYPE TYPE_JB_CONFIG2   ! Namelist variables that must be set up after the
                       ! variables in TYPE_JB_CONFIG
  INTEGER(KIND=JPIM), DIMENSION(100) :: M_GRIBCODES
  INTEGER(KIND=JPIM), DIMENSION(100) :: M_GRIBCODES_EP
  INTEGER(KIND=JPIM), DIMENSION(100) :: M_GRIBCODES_FCE
  CHARACTER(LEN=40) , DIMENSION(100) :: C_COR_STRINGS
END TYPE TYPE_JB_CONFIG2

TYPE TYPE_JB_BUF_STRUCT
  REAL(KIND=JPRB), ALLOCATABLE :: BUF(:,:,:)
  LOGICAL :: LFILLED = .FALSE.
END TYPE TYPE_JB_BUF_STRUCT
 
TYPE TYPE_JB_DATA
  LOGICAL :: LJGNRSGP
  LOGICAL :: LSUBBG
  INTEGER(KIND=JPIM) :: NOFEP
  INTEGER(KIND=JPIM) :: NOFCEF
  INTEGER(KIND=JPIM) :: NAEROCV
  INTEGER(KIND=JPIM) :: NOMSPEC
  INTEGER(KIND=JPIM), ALLOCATABLE :: ISPECFIDS(:)

  REAL(KIND=JPRB), ALLOCATABLE :: REDINC(:)
  REAL(KIND=JPRB), ALLOCATABLE :: FP2TPS(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: FGSCNM3D(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: FGSCNM2D(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: TMEANUVER(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: FCEMN3D(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: FCEMN2D(:)
  REAL(KIND=JPRB), ALLOCATABLE :: FCEIMN(:)
  REAL(KIND=JPRB), ALLOCATABLE :: FGMWNE(:)
  REAL(KIND=JPRB), ALLOCATABLE :: STPS(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: SDIV(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: SO3(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: BFACT1(:)
  REAL(KIND=JPRB), ALLOCATABLE :: BFACT2(:)
  REAL(KIND=JPRB), ALLOCATABLE :: SHUM(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: BFACT(:)
  REAL(KIND=JPRB), ALLOCATABLE :: SCORI_R(:)
  REAL(KIND=JPRB), ALLOCATABLE :: SCORI_I(:)

  TYPE (TYPE_JB_BUF_STRUCT), ALLOCATABLE :: EP(:) 
  TYPE (TYPE_JB_BUF_STRUCT) :: FCE
  TYPE (TYPE_JB_BUF_STRUCT) :: FCEREN
  TYPE (TYPE_JB_BUF_STRUCT) :: JBVCOORD
  TYPE(SPECTRAL_FIELD)      :: SPJB
END TYPE TYPE_JB_DATA

TYPE TYPE_FCEMN_STRUCT
  REAL(KIND=JPRB) :: PS
  REAL(KIND=JPRB) :: PSU
  REAL(KIND=JPRB), ALLOCATABLE :: U (:)
  REAL(KIND=JPRB), ALLOCATABLE :: Z (:)
  REAL(KIND=JPRB), ALLOCATABLE :: T (:)
  REAL(KIND=JPRB), ALLOCATABLE :: O3(:)
END TYPE TYPE_FCEMN_STRUCT

TYPE TYPE_VCORS_STRUCT
  REAL(KIND=JPRB), ALLOCATABLE :: EVALS(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: EVECS(:,:,:)
END TYPE TYPE_VCORS_STRUCT

TYPE TYPE_SPJB_VARS_INFO_STRUCT
  INTEGER(KIND=JPIM) :: IFID          ! Field identifier
  INTEGER(KIND=JPIM) :: IGRIBCODE     ! GRIB code for model field
  INTEGER(KIND=JPIM) :: IGRIBCODE_FCE ! GRIB code for background error
  INTEGER(KIND=JPIM) :: IPT           ! index into model (SPA3/SPA2/SPGFL/GLF)
  INTEGER(KIND=JPIM) :: IPTJB         ! index into SPA3JB/SPA2JB
  INTEGER(KIND=JPIM) :: IPTECV        ! index into RSPECV2D/RSPECV3D/RGPECV3D
  INTEGER(KIND=JPIM) :: IPTFCE        ! index into FCE%BUF
  INTEGER(KIND=JPIM) :: IPTBG         ! index into SP7A3, GP7A3 or SP7A2
  LOGICAL            :: L_IN_SPA3     ! .T. => maps to model SPA3 array
  LOGICAL            :: L_IN_SPA2     ! .T. => maps to model SAP2 array
  LOGICAL            :: L_IN_SPGFL    ! .T. => maps to model SPGFL array
  LOGICAL            :: L_IN_GPGFL    ! .T. => maps to model GPGFL array
  LOGICAL            :: L_SPEC_MODEL  ! True if spectral field in model
  LOGICAL            :: L_IN_ECVSP_2D ! .T. => maps to model SPA2 array
  LOGICAL            :: L_IN_ECVGP_2D ! .T. => maps to model GPGFL2 array
  LOGICAL            :: L_IN_ECVSP_3D ! .T. => maps to model SPA3 array
  LOGICAL            :: L_IN_ECVGP_3D ! .T. => maps to model GPGFL array
  LOGICAL            :: L_IN_ECV      ! .T. => any of the above ECV .T.
  CHARACTER(LEN=40)  :: COR_STRING    ! Identifying string in ".cv" file
END TYPE TYPE_SPJB_VARS_INFO_STRUCT

TYPE TYPE_HYJB_VARS_INFO_STRUCT
  INTEGER(KIND=JPIM) :: IPTEP         ! index into EP(:)%BUF
  INTEGER(KIND=JPIM) :: IGRIBCODE_EP  ! GRIB code for ensemble perturbation
END TYPE TYPE_HYJB_VARS_INFO_STRUCT


TYPE TYPE_JB_STRUCT
  TYPE (TYPE_JB_CONFIG)                          :: CONFIG
  TYPE (TYPE_JB_CONFIG2)                         :: CONFIG2
  TYPE (TYPE_JB_DATA)                            :: JB_DATA
  TYPE (TYPE_FCEMN_STRUCT)                       :: FCEMN
  TYPE (TYPE_VCORS_STRUCT),          ALLOCATABLE :: VCORS(:)
  TYPE (TYPE_SPJB_VARS_INFO_STRUCT), ALLOCATABLE :: SPJB_VARS_INFO(:)
  TYPE (TYPE_HYJB_VARS_INFO_STRUCT), ALLOCATABLE :: HYJB_VARS_INFO(:)
  TYPE (TYPE_WAVELETJB_CONFIG)                   :: WJBCONF
  TYPE (TYPE_WAVELETJB_DATA)                     :: WJBDATA
  TYPE (TYPE_WAVELETJB_VCOR_STRUCT), ALLOCATABLE :: WAVELET_VCORS(:,:)
  TYPE (TYPE_WAVELETJB_GRID_STRUCT), ALLOCATABLE :: GRID_DEFINITION(:)
  TYPE (TYPE_LAMWAVELETJB_CONFIG)                :: LAMWAVELET_CONFIG
  TYPE (TYPE_LAMWAVELETJB_INFO)                  :: LAMWAVELET_INFO
  TYPE (TYPE_LAMWAVELET_VCOR_STRUCT)             :: LAMWAVELET_VCOR
  TYPE (TYPE_LAMWAVELET_GRID_STRUCT)             :: LAMWAVELET_GRID
  TYPE (TYPE_JBCHVAR)                            :: JBCHVAR
END TYPE TYPE_JB_STRUCT

TYPE(TYPE_JB_STRUCT), POINTER :: JB_STRUCT => NULL()

END MODULE YOMJG

