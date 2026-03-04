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

MODULE YOMSPSDT

! Purpose :
! -------
!    Define types for the SPPT scheme, which represents model uncertainties
!    with Stochastically Perturbed Parameterisation Tendencies.

! Interface :
! ---------
!    Empty.

! External :
! --------
!    None.

! Method :
! ------
!    See Documentation.

! Reference :
! ---------
!    Tech Memo 598: Stochastic parameterisations, Palmer et al. (2009)

! Author :
! ------
!
!    Original :

! Modifications :
! -------------
!    SJ. Lock : Jan-2016 : General tidying & enabling multiple patterns ("iSPPT")
!    SJ. Lock : Oct-2016 : Option to exclude clear-skies radiation from SPPT
!     S. Lang : Feb-2017 : Modifications to cycle pattern in data assimilation
!    SJ  Lock : Oct 2017 : Options to reduce frequency of pattern updates
!   M. Chrust : Aug 2021 : Remove global SPPT objects -> part of the model object
!
!-----------------------------------------------------------------------------


USE PARKIND1,ONLY : JPIM, JPRB
USE SPECTRAL_ARP_MOD, ONLY : SPECTRAL_ARP
USE GRIDPOINT_FIELDS_MIX, ONLY : GRIDPOINT_FIELD

IMPLICIT NONE
SAVE

!     ------------------------------------------------------

!*    control parameters for Stochastically Perturbed Parametrization Tendency
!*    with independent patterns option: "iSPPT"

INTEGER(KIND=JPIM), PARAMETER :: NMAXSCALES_SDT=5  !> max number of pattern-scales
INTEGER(KIND=JPIM), PARAMETER :: NMAXINDPAT_SDT=6  !> max number of independent patterns
INTEGER(KIND=JPIM), PARAMETER :: NWRMAX_SDT=20     !> max number of output times
!=============================================================================

!      1.    TYPE DEFINITIONS
!            ----------------

!      1.1   TYPE holding the configuration of the SPPT scheme
TYPE TSPPT_CONFIG
  LOGICAL  :: LSPSDT=.FALSE.                         !> activates SPPT scheme
  LOGICAL  :: LRADCLR_SDT                            !> .F.: exclude clear-skies radiation tendencies from SPPT
  LOGICAL  :: LSATADJ_SDT                            !> .F.: exclude saturation adjustment tendencies from SPPT
  LOGICAL  :: LSPPT1                                 !> .T.: "Standard SPPT": MPSDT=(1,1,1,1,1,1)
  LOGICAL  :: LSPPTGFIX                              !> .T. constrain global integrals of perturbed tendencies
                                                     !     to those from the unperturbed tendencies
  LOGICAL  :: LCLIP_SPEC_SDT                         !> clip spherical harmonics coefficients
  LOGICAL  :: LCLIP_GRID_SDT                         !> clip grid-point pattern field
  LOGICAL  :: LQPERTLIMIT2                           !> quasi-symmetric limiting of humidity perturbations
                                                     !>   between 0 and saturation
  LOGICAL  :: LTAPER_BL0                             !> compute taper functions for the boundary layer (BL)
  LOGICAL, DIMENSION(NMAXINDPAT_SDT)  :: LTAPER_BLI  !> apply BL tapering for each physics tendency
  LOGICAL  :: LTAPER_ST0                             !> compute taper functions for the stratosphere (ST)
  LOGICAL, DIMENSION(NMAXINDPAT_SDT)  :: LTAPER_STI  !> apply ST tapering for each physics tendency

  LOGICAL  :: LUSESETRAN_SDT=.TRUE.      !> use setran to modify the seed value
  LOGICAL  :: LRESETSEED_SDT=.FALSE.     !> reset the random number seed to get reproducible results
                                         !     when restart is performed
  LOGICAL  :: LABSTIMSEED_SDT            !> use absolute time for random number seed (initial time + step)

  !I/O:
  LOGICAL  :: LWRITE_ARP                 !> write ARP to GRIB file
  LOGICAL  :: LRDPATINIT_SDT=.FALSE.     !> read initial state of pattern from GRIB file
  LOGICAL  :: LWRPATTRUN_SDT=.FALSE.     !> write state of pattern at truncation time to GRIB file

  !KEYMASK switches for use in CALLPAR
  LOGICAL, POINTER  :: LACTIVEQ_SDT => NULL()   !> .T.: Q arrays are active
  LOGICAL, POINTER  :: LACTIVEO_SDT => NULL()   !> .F.: all others (O3, A, Qx) are inactive


  INTEGER(KIND=JPIM) :: NSCALES_SDT=1                      !> number of pattern scales

  INTEGER(KIND=JPIM) :: NPATFR=1                           !> frequency of pattern updates:
                                                           !>  >0: every NPATFR timesteps
                                                           !>  <0: every -NPATFR hours

  INTEGER(KIND=JPIM) :: RESETSEEDFRQ_SDT                   !> LRESETSEED_SDT: frequency (timesteps) for
                                                           !     resetting the random number seed

  INTEGER(KIND=JPIM) :: NWRPATTRUN_SDT                     !> number of output times
  INTEGER(KIND=JPIM), DIMENSION(NWRMAX_SDT) :: NHOUR_PATTRUN  !> hour that pattern is written to grib file
  INTEGER(KIND=JPIM), DIMENSION(NWRMAX_SDT) :: NSTEP_PATTRUN  !> step that pattern is written to grib file

  INTEGER(KIND=JPIM) :: NCLIP_KIND_SDT                     !> kind of threshold used to define clipping:
                                                           !>   0: use XCLIP_RATIO_SDT
                                                           !>   1: compute XCLIP_RATIO_SDT from XCLIPM_SDT

  INTEGER(KIND=JPIM) :: NQSAT_SDT=3                        !> treatment of supersaturation:
                                                           !>   0: original treatment
                                                           !>   3: modified for vertical consistency

  INTEGER(KIND=JPIM), DIMENSION(4000) :: NWN_TESTSPEC_SDT  !> wavenumbers with nonzero power in test spectrum
  INTEGER(KIND=JPIM) :: NWAV_TESTSPEC_SDT                  !> count of nonzero power wavenumbers


  REAL(KIND=JPRB) :: SDEV_SDT(NMAXSCALES_SDT,NMAXINDPAT_SDT)      !> standard deviation of patterns
  REAL(KIND=JPRB), DIMENSION(NMAXSCALES_SDT)    :: TAU_SDT        !> correlation time scale (s)
  REAL(KIND=JPRB), DIMENSION(NMAXSCALES_SDT)    :: XLCOR_SDT      !> correlation length scale (m)
  REAL(KIND=JPRB), DIMENSION(NMAXSCALES_SDT)    :: XWC_KAPPA_T    !> controls Weaver&Courtier (2001) spectrum

  REAL(KIND=JPRB) :: XCLIPM_SDT(NMAXINDPAT_SDT)                   !> max modulus of perturbation in gridpoint space
  REAL(KIND=JPRB) :: XCLIP_RATIO_SDT(NMAXINDPAT_SDT)              !> relative amplitude of clipping
  REAL(KIND=JPRB) :: XMEANRED_THICKNESS_SDT                       !> layer thickness (Pa) for reduction factor
  REAL(KIND=JPRB) :: SDEVTOT_SDT(NMAXINDPAT_SDT)                  !> net std dev of multiscale pattern

  REAL(KIND=JPRB) :: XSIGMATOP                   !> normalized pressure of top of BL tapering layer
  REAL(KIND=JPRB) :: XSIGMABOT                   !> normalized pressure of bottom of BL tapering layer
  REAL(KIND=JPRB) :: XPRESSTOP_ST0               !> pressure at top of STratosphere tapering layer
  REAL(KIND=JPRB) :: XPRESSBOT_ST0               !> pressure at bottom of ST tapering layer
  REAL(KIND=JPRB) :: XTAPER3, XTAPER2, XTAPER1, XTAPER0                   !> for the BL tapering function
  REAL(KIND=JPRB) :: XTAPER3_ST0, XTAPER2_ST0, XTAPER1_ST0, XTAPER0_ST0   !> for the ST tapering function

  REAL(KIND=JPRB) :: XSPPTGFIXRLX                !> relaxation factor

  CHARACTER(LEN=128) :: CSPEC_SHAPE_SDT='WeaverCourtier'     !> shape of variance spectrum
  CHARACTER(LEN=256), DIMENSION(:), POINTER :: CIPATINIT_SDT  => NULL()   !> filename for initial state of pattern
  CHARACTER(LEN=256), DIMENSION(:), POINTER :: COPATTRUN_SDT  => NULL()   !> filename for pattern at VAREPS truncation timestep
  CHARACTER(LEN=256), DIMENSION(:), POINTER :: COPATSP_SDT    => NULL()   !> filename for spectral representation of AR(1)
  CHARACTER(LEN=256), DIMENSION(:), POINTER :: COPATGP_SDT    => NULL()   !> filename for gridpoint representation of AR(1),
END TYPE TSPPT_CONFIG

TYPE TSPPT_DATA
  INTEGER(KIND=JPIM), DIMENSION(NMAXINDPAT_SDT)       :: MPSDT       !> indices referring to "SPPT" patterns

  INTEGER(KIND=JPIM) :: N2D                                          !> number of (independent) 2D patterns
  INTEGER(KIND=JPIM), DIMENSION(:), POINTER           :: NSEED_SDT => NULL()  !> seeds for (independent) 2D patterns

  TYPE (SPECTRAL_ARP),     DIMENSION(:), POINTER :: YSPSDT_AR1 => NULL()      !> array with patterns
  TYPE (GRIDPOINT_FIELD),  DIMENSION(:), POINTER :: YGPSDT     => NULL()      !> array with grid-point fields
  TYPE (GRIDPOINT_FIELD),  DIMENSION(:), POINTER :: YGPSDT0    => NULL()      !> array with first  time-level
  TYPE (GRIDPOINT_FIELD),  DIMENSION(:), POINTER :: YGPSDT1    => NULL()      !> array with second time-level
END TYPE TSPPT_DATA

!=============================================================================

END MODULE YOMSPSDT
