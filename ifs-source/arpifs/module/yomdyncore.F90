! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE YOMDYNCORE

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

! MAIN SWITCHES

! LAQUA=.T.       switches to initialize the surface fields 
!                 to run an aqua-planet (all planet covered by sea)
!                 with "forced" SST.
! The SST setup is selected by the value of MSSTSCHEME (see below)
! LAPE=.T.  Extra modifications for radiation/geoide
!           needed to reproduce the AquaPlanet Experiment
!           from Neale and Hoskins, 2001 (needs both LAPE and LAQUA=T)
! LHELDSUAREZ=.T. switches on idealized Held-Suarez

! LPPSTEPS=.T.    specify postprocessing/fdb output in steps (emulating hours)

LOGICAL         :: LAQUA, LAPE
LOGICAL         :: LHELDSUAREZ

LOGICAL         :: LPPSTEPS

! small planet earth radius factor
REAL(KIND=JPRB)     :: RPLRADI, RUSPLRADI
! Coriolis acceleration factor
REAL(KIND=JPRB)     :: RCORIOI, RUSCORIOI
! small planet gravity factor
REAL(KIND=JPRB)     :: RPLRG, RUSPLRG
! small planet DARE factor
REAL(KIND=JPRB)     :: RPLDARE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Test case selection
!!!!!!!!!!!!!!!!!!!!!
! Each test case is defined by a family (CTEST_FAMILY)
! and a case number (NTESTCASE)
CHARACTER(LEN=12) :: CTEST_FAMILY
INTEGER(KIND=JPIM) :: NTESTCASE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Amplitude (can be used for perturbation, mountain height etc
REAL(KIND=JPRB)     :: RAMP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOISEVOR : if 1 add initial noise in vorticity
INTEGER(KIND=JPIM) :: NOISEVOR

!     PARAMETER USED FOR AQUAPLANET

! defined in GP_SSTAQUA()
! MSSTSCHEME: defines the prescribed SST distribution as described in
! Atmos. Sc. Letters (2001) Vol.1
! == 1 Control
! == 2 Peaked
! == 3 Flat
! == 4 Qobs
! == 5 Control5N
! == 6 1KEQ, chi=1.0_JPRB
! == 7 3KEQ, chi=3.0_JPRB
! == 8 3KW1, chi=3.0_JPRB
! == 9 constant SST everywhere

! RLAMBDA0       : Longitude of maximum SST anomaly
! RLAMBDAD, RPHID: Half width in long/lat of SST anomaly
! RCHI           : Maximum magnitude of SST anomaly 

INTEGER(KIND=JPIM)  :: MSSTSCHEME
REAL(KIND=JPRB)     :: RCHI
REAL(KIND=JPRB)     :: RLAMBDA0
REAL(KIND=JPRB)     :: RLAMBDAD
REAL(KIND=JPRB)     :: RPHID
!-----------------------------------------------------------------------------

!    INITIAL PROFILES

! RT00_DYN  -  Initial idealized temperature
! RP00_DYN  -  Initial idealized pressure
! RU00_DYN  -  Initial idealized zonal wind
! RST0_DYN  -  Initial idealized wind shear (du/dz)
REAL(KIND=JPRB)  :: RT00_DYN
REAL(KIND=JPRB)  :: RP00_DYN
REAL(KIND=JPRB)  :: RU00_DYN
REAL(KIND=JPRB)  :: RST0_DYN

!     PARAMETER USED FOR THE ORIGINAL HELD/SUAREZ FORCING

! RSIGMAB          - HS test: vertical extent of friction layer expressed in
!                    sigma (between 0..1)
! REFERENCE_LEVEL  - HS test: 1.0_JPRB-RSIGMAB
! RDELTA_T         - HS test: pole - equator temperature difference
! RDELTA_THETA     - HS test: tropical heating differential
! RPRESSURE_SCALE  - HS test: pressure scaling
! RHS_KF           - HS test: frictional relaxation time scale
! RHS_KA           - HS test: heating time scale (atmosphere) 0.025_JPRB*RHS_KF
! RHS_KS           - HS test: heating time scale (surface) 0.25_JPRB*RHS_KF
! RHS_KS_KA        - HS test: RHS_KS-RHS_KA
! RHS_KAPPA        - HS test: Kappa = 2/7 
REAL(KIND=JPRB) :: RSIGMAB
REAL(KIND=JPRB) :: REFERENCE_LEVEL
REAL(KIND=JPRB) :: RDELTA_T
REAL(KIND=JPRB) :: RDELTA_THETA
REAL(KIND=JPRB) :: RPRESSURE_SCALE
REAL(KIND=JPRB) :: RHS_KF
REAL(KIND=JPRB) :: RHS_KA
REAL(KIND=JPRB) :: RHS_KS
REAL(KIND=JPRB) :: RHS_KS_KA
REAL(KIND=JPRB) :: RHS_KAPPA

!     PARAMETER USED FOR THE MODIFIED FORCING IN UPPER LEVELS
!     (AFTER WILLIAMSON)

! L_HS_WILLIAMSON=.FALSE. - .T. use the Williamson extension of the Held-Suarez
!  test in the stratosphere (MWR 126:1001-1012)
LOGICAL         :: L_HS_WILLIAMSON=.FALSE.

! RGAMMA_D         - HSW extension: lapse rate active when pressure <= p_d,
!                    decreases temperature in stratospheric polar regions
! RGAMMA_I         - HSW extension: lapse rate active when pressure <= p_i,
!                    increases temperature in stratospheric equatorial regions
! RTHSW_0          - HSW extension: stratospheric reference temperature
! RPHSW_D          - HSW extension: pressure threshold p_d
! RPHSW_EQ         - HSW extension: equatorial stratospheric reference pressure
! RPHSW_PL         - HSW extension: vertical transition threshold in hPa
! RPHSW_EQ_PL_HALF - HSW extension: 0.5_JPRB*(RPHSW_EQ-RPHSW_PL)
! RPHSW_D_REV      - HSW extension: 1.0_JPRB/RPHSW_D
REAL(KIND=JPRB), PARAMETER :: RGAMMA_D = 2.E-3_JPRB
REAL(KIND=JPRB), PARAMETER :: RGAMMA_I = -3.345E-3_JPRB
REAL(KIND=JPRB) :: RTHSW_0
REAL(KIND=JPRB) :: RPHSW_D
REAL(KIND=JPRB) :: RPHSW_EQ
REAL(KIND=JPRB) :: RPHSW_PL
REAL(KIND=JPRB) :: RPHSW_EQ_PL_HALF
REAL(KIND=JPRB) :: RPHSW_D_REV

!     PARAMETER USED FOR SAMPLING STATISTICS

! RSAMPLING_START     -  start of time mean calculations
! RSAMPLING_INTERVAL  -  interval used in time mean calculations
REAL(KIND=JPRB), PARAMETER :: RSAMPLING_START=0.0_JPRB
REAL(KIND=JPRB), PARAMETER :: RSAMPLING_INTERVAL=12._JPRB*3600._JPRB

! PARAMETER USED FOR SEMI-LAGRANGIAN DIAGNOSTICS

! RJACSLDIA == THRESHOLD FOR OUTPUT OF JACOBIAN
! RLIPSLDIA == THRESHOLD FOR OUTPUT OF LIPSCHITZ NUMBER

REAL(KIND=JPRB) :: RJACSLDIA
REAL(KIND=JPRB) :: RLIPSLDIA

END MODULE YOMDYNCORE
