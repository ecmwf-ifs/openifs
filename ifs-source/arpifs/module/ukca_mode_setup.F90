! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module to store MODE setup arrays after call to UKCA_MODE_SETUP_ALL.
!    Contains public subroutines:
!      UKCA_MODE_IMSCAVCOFF
!      UKCA_MODE_ALLCP_4MODE
!      UKCA_MODE_SUSS_4MODE
!      UKCA_MODE_SUSSBCOC_4MODE
!      UKCA_MODE_SUSSBCOC_5MODE
!      UKCA_MODE_SUSSBCOCSO_5MODE
!      UKCA_MODE_SUSSBCOCSO_4MODE
!      UKCA_MODE_DUonly_2MODE
!      UKCA_MODE_DUonly_3MODE (needs to be added at some point)
!      UKCA_MODE_SUSSBCOCDU_7MODE
!      UKCA_MODE_SUSSBCOCDU_4MODE
!    which define modes and components for different components/modes setup.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
! Subroutine Interface:
      MODULE UKCA_MODE_SETUP
! ---------------------------------------------------------------------|
!+ Module to contain modes and components
!
! Description:
! To allow use throughout UM, module stores MODE setup arrays after
! call to UKCA_MODE_SETUP_ALL.
!
! Note: currently code is hard-coded so that ordering of modes must
! 1) nucln, 2)   soluble Aitken, 3)   soluble accum, 4)   soluble coarse
!           5) insoluble Aitken, 6) insoluble accum, 7) insoluble coarse
!
! History:
! Version    Date      Author         Description
! -----      --------  ---------      ------------
! vn6.1      1/XII/06  Graham Mann    Original code
! vn6.1      23/IV/07  Colin Johnson  Converted to MODULE procedure
! vn6.1      25/VII/07 Graham Mann    Modified to include extra fields
!                                     and optional setups by module
!                                     procedures (also includes coeffs
!                                     for impaction scavenging)
!
! ---------------------------------------------------------------------|

      USE YOMHOOK,  ONLY: LHOOK, DR_HOOK, JPHOOK
      USE PARKIND1, ONLY: JPRB, JPIM
      IMPLICIT NONE
      SAVE

      INTEGER(KIND=JPIM), PARAMETER :: NMODES=7  ! No of modes
      INTEGER(KIND=JPIM), PARAMETER :: NCP=6     ! No of components
      INTEGER(KIND=JPIM), PARAMETER :: NCATION=3 ! No possible cation species
      INTEGER(KIND=JPIM), PARAMETER :: NANION =4 ! No possible anion species
!
      INTEGER(KIND=JPIM), PARAMETER :: CP_SU=1  ! Index to store SO4    cpt
      INTEGER(KIND=JPIM), PARAMETER :: CP_BC=2  ! Index to store BC     cpt
      INTEGER(KIND=JPIM), PARAMETER :: CP_OC=3  ! Index to store 1st OC cpt
      INTEGER(KIND=JPIM), PARAMETER :: CP_CL=4  ! Index to store NaCl   cpt
      INTEGER(KIND=JPIM), PARAMETER :: CP_DU=5  ! Index to store dust   cpt
      INTEGER(KIND=JPIM), PARAMETER :: CP_SO=6  ! Index to store 2nd OC cpt
      INTEGER(KIND=JPIM), PARAMETER :: CP_NI=7  ! Index to store NI cpt
      INTEGER(KIND=JPIM), PARAMETER :: CP_AM=8  ! Index to store AM cpt

! Mode switches (1=on, 0=0ff)
      INTEGER(KIND=JPIM), DIMENSION(NMODES) :: MODE_CHOICE
! Component switches (1=on, 0=off)
      INTEGER(KIND=JPIM), DIMENSION(NCP)    :: COMPONENT_CHOICE
! Components that are soluble
      INTEGER(KIND=JPIM), DIMENSION(NCP)    :: SOLUBLE_CHOICE
! Components allowed in each mode (must be consistent with coag_mode)
      INTEGER(KIND=JPIM), DIMENSION(NMODES,NCP) :: COMPONENT_MODE
! Modes resulting when two modes coagulate
      INTEGER(KIND=JPIM), DIMENSION(NMODES,NMODES) :: COAG_MODE
! Specify which modes are soluble
      INTEGER(KIND=JPIM), DIMENSION(NMODES) :: MODESOL
!
! Variables for impaction scavenging (as in Pringle, 2006 PhD thesis)
      INTEGER(KIND=JPIM), PARAMETER  :: NCOLL=20 ! # of columns in LUT (aer. bins)
      INTEGER(KIND=JPIM), PARAMETER  :: NROW =19 ! # of rows in LUT (raindrop bins)

! Tracer indices set in ukca_aero_ctl:
      INTEGER(KIND=JPIM) :: II_ND(NMODES)        ! indices in mode_tracers of NCONC
      INTEGER(KIND=JPIM) :: II_MD(NMODES,NCP)    ! indices in mode_tracers of CPTMMR

! Molar masses of components (kg mol-1)
      REAL(KIND=JPRB), DIMENSION(NCP)    :: MM
! Mass density of components (kg m^-3)
      REAL(KIND=JPRB), DIMENSION(NCP)    :: RHOCOMP
! Number of dissociating ions in soluble components
      REAL(KIND=JPRB), DIMENSION(NCP)    :: NO_IONS
! Lower size limits of geometric mean radius for each mode
      REAL(KIND=JPRB), DIMENSION(NMODES) :: FRACBCEM
! Fraction of bc ems to go into each mode
      REAL(KIND=JPRB), DIMENSION(NMODES) :: FRACOCEM
! Fraction of oc ems to go into each mode
      REAL(KIND=JPRB), DIMENSION(NMODES) :: DDPLIM0
! Upper size limits of geometric mean radius for each mode
      REAL(KIND=JPRB), DIMENSION(NMODES) :: DDPLIM1
! Mid-point of size mode (m)
      REAL(KIND=JPRB), DIMENSION(NMODES) :: DDPMID
! Mid-point masses for initial radius grid
      REAL(KIND=JPRB), DIMENSION(NMODES) :: MMID
! Lo-interf masses for initial radius grid
      REAL(KIND=JPRB), DIMENSION(NMODES) :: MLO
! Hi-interf masses for initial radius grid
      REAL(KIND=JPRB), DIMENSION(NMODES) :: MHI
! Fixed geometric standard deviation for each mode
      REAL(KIND=JPRB), DIMENSION(NMODES) :: SIGMAG
! EXP((9/2)*LOG^2(SIGMA_G))
      REAL(KIND=JPRB), DIMENSION(NMODES) :: X
! Threshold for number in mode to carry out calculations
      REAL(KIND=JPRB), DIMENSION(NMODES) :: NUM_EPS
! Initial fractions of mass in each mode among components
      REAL(KIND=JPRB), DIMENSION(NMODES,NCP) :: MFRAC_0

      REAL(KIND=JPRB), DIMENSION(NROW)       :: RADDROP  ! raindrop bins
      REAL(KIND=JPRB), DIMENSION(NCOLL,NROW) :: COLLEFF4 ! collision efficiency

! Mode names
      CHARACTER(LEN=7),DIMENSION(NMODES) :: MODE_NAMES
! Component names
      CHARACTER(LEN=7),DIMENSION(NCP)    :: COMPONENT_NAMES

! Modes (T/F)
      LOGICAL, DIMENSION(NMODES)     :: MODE
! Components set in each mode (T/F)
      LOGICAL, DIMENSION(NMODES,NCP) :: COMPONENT
! Components which are soluble (T/F)
      LOGICAL, DIMENSION(NCP)        :: SOLUBLE

      CONTAINS

      SUBROUTINE UKCA_MODE_IMSCAVCOFF

      IMPLICIT NONE
      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_IMSCAVCOFF',0,ZHOOK_HANDLE)

      RADDROP =(/   1.0_JPRB, 1.587_JPRB,  2.52_JPRB,   4.0_JPRB,  6.35_JPRB, 10.08_JPRB,             &
       &           16.0_JPRB,  25.4_JPRB, 40.32_JPRB,  64.0_JPRB, 101.6_JPRB, 161.3_JPRB,             &
       &          256.0_JPRB, 406.4_JPRB, 645.1_JPRB,1024.0_JPRB,1625.0_JPRB,2580.0_JPRB,4096.0_JPRB/)

      COLLEFF4( 1,1:NROW)=(/ 0.522E+05_JPRB,0.139E+05_JPRB,0.328E+04_JPRB,0.775E+03_JPRB,   &
       &                     0.183E+03_JPRB,0.432E+02_JPRB,0.102E+02_JPRB,0.291E+01_JPRB,   &
       &                     0.108E+01_JPRB,0.439E+00_JPRB,0.201E+00_JPRB,0.110E+00_JPRB,   &
       &                     0.633E-01_JPRB,0.366E-01_JPRB,0.815E-02_JPRB,0.168E-02_JPRB,   &
       &                     0.394E-03_JPRB,0.591E-04_JPRB,0.132E-04_JPRB/)

      COLLEFF4( 2,1:NROW)=(/ 0.126E+05_JPRB,0.373E+04_JPRB,0.985E+03_JPRB,0.260E+03_JPRB,   &
       &                     0.687E+02_JPRB,0.182E+02_JPRB,0.480E+01_JPRB,0.150E+01_JPRB,   &
       &                     0.536E+00_JPRB,0.224E+00_JPRB,0.107E+00_JPRB,0.608E-01_JPRB,   &
       &                     0.364E-01_JPRB,0.236E-01_JPRB,0.731E-02_JPRB,0.167E-02_JPRB,   &
       &                     0.395E-03_JPRB,0.592E-04_JPRB,0.133E-04_JPRB/)

      COLLEFF4( 3,1:NROW)=(/ 0.445E+04_JPRB,0.139E+04_JPRB,0.390E+03_JPRB,0.110E+03_JPRB,   &
       &                     0.308E+02_JPRB,0.864E+01_JPRB,0.243E+01_JPRB,0.783E+00_JPRB,   &
       &                     0.270E+00_JPRB,0.117E+00_JPRB,0.564E-01_JPRB,0.325E-01_JPRB,   &
       &                     0.203E-01_JPRB,0.136E-01_JPRB,0.595E-02_JPRB,0.166E-02_JPRB,   &
       &                     0.396E-03_JPRB,0.594E-04_JPRB,0.133E-04_JPRB/)

      COLLEFF4( 4,1:NROW)=(/ 0.259E+04_JPRB,0.810E+03_JPRB,0.227E+03_JPRB,0.639E+02_JPRB,   &
       &                     0.179E+02_JPRB,0.503E+01_JPRB,0.141E+01_JPRB,0.440E+00_JPRB,   &
       &                     0.146E+00_JPRB,0.662E-01_JPRB,0.309E-01_JPRB,0.177E-01_JPRB,   &
       &                     0.115E-01_JPRB,0.728E-02_JPRB,0.446E-02_JPRB,0.164E-02_JPRB,   &
       &                     0.399E-03_JPRB,0.597E-04_JPRB,0.134E-04_JPRB/)

      COLLEFF4( 5,1:NROW)=(/ 0.196E+04_JPRB,0.602E+03_JPRB,0.166E+03_JPRB,0.457E+02_JPRB,   &
       &                     0.126E+02_JPRB,0.346E+01_JPRB,0.953E+00_JPRB,0.280E+00_JPRB,   &
       &                     0.927E-01_JPRB,0.410E-01_JPRB,0.183E-01_JPRB,0.990E-02_JPRB,   &
       &                     0.655E-02_JPRB,0.444E-02_JPRB,0.389E-02_JPRB,0.164E-02_JPRB,   &
       &                     0.402E-03_JPRB,0.604E-04_JPRB,0.135E-04_JPRB/)

      COLLEFF4( 6,1:NROW)=(/ 0.192E+04_JPRB,0.569E+03_JPRB,0.151E+03_JPRB,0.401E+02_JPRB,   &
       &                     0.106E+02_JPRB,0.282E+01_JPRB,0.749E+00_JPRB,0.206E+00_JPRB,   &
       &                     0.689E-01_JPRB,0.280E-01_JPRB,0.119E-01_JPRB,0.580E-02_JPRB,   &
       &                     0.384E-02_JPRB,0.304E-02_JPRB,0.393E-02_JPRB,0.165E-02_JPRB,   &
       &                     0.406E-03_JPRB,0.612E-04_JPRB,0.138E-04_JPRB/)

      COLLEFF4( 7,1:NROW)=(/ 0.208E+04_JPRB,0.604E+03_JPRB,0.156E+03_JPRB,0.405E+02_JPRB,   &
       &                     0.105E+02_JPRB,0.271E+01_JPRB,0.703E+00_JPRB,0.185E+00_JPRB,   &
       &                     0.572E-01_JPRB,0.217E-01_JPRB,0.903E-02_JPRB,0.462E-02_JPRB,   &
       &                     0.280E-02_JPRB,0.208E-02_JPRB,0.398E-02_JPRB,0.168E-02_JPRB,   &
       &                     0.415E-03_JPRB,0.628E-04_JPRB,0.142E-04_JPRB/)

      COLLEFF4( 8,1:NROW)=(/ 0.233E+04_JPRB,0.636E+03_JPRB,0.162E+03_JPRB,0.410E+02_JPRB,   &
       &                     0.104E+02_JPRB,0.264E+01_JPRB,0.671E+00_JPRB,0.174E+00_JPRB,   &
       &                     0.513E-01_JPRB,0.184E-01_JPRB,0.747E-02_JPRB,0.384E-02_JPRB,   &
       &                     0.222E-02_JPRB,0.161E-02_JPRB,0.408E-02_JPRB,0.173E-02_JPRB,   &
       &                     0.430E-03_JPRB,0.657E-04_JPRB,0.149E-04_JPRB/)

      COLLEFF4( 9,1:NROW)=(/ 0.235E+04_JPRB,0.659E+03_JPRB,0.165E+03_JPRB,0.412E+02_JPRB,   &
       &                     0.103E+02_JPRB,0.257E+01_JPRB,0.643E+00_JPRB,0.168E+00_JPRB,   &
       &                     0.490E-01_JPRB,0.168E-01_JPRB,0.661E-02_JPRB,0.326E-02_JPRB,   &
       &                     0.188E-02_JPRB,0.140E-02_JPRB,0.422E-02_JPRB,0.180E-02_JPRB,   &
       &                     0.452E-03_JPRB,0.698E-04_JPRB,0.160E-04_JPRB/)

      COLLEFF4(10,1:NROW)=(/ 0.165E+04_JPRB,0.457E+03_JPRB,0.112E+03_JPRB,0.277E+02_JPRB,   &
       &                     0.680E+01_JPRB,0.167E+01_JPRB,0.412E+00_JPRB,0.106E+00_JPRB,   &
       &                     0.304E-01_JPRB,0.999E-02_JPRB,0.386E-02_JPRB,0.186E-02_JPRB,   &
       &                     0.124E-02_JPRB,0.140E-02_JPRB,0.447E-02_JPRB,0.193E-02_JPRB,   &
       &                     0.491E-03_JPRB,0.771E-04_JPRB,0.179E-04_JPRB/)

      COLLEFF4(11,1:NROW)=(/ 0.899E+03_JPRB,0.246E+03_JPRB,0.597E+02_JPRB,0.145E+02_JPRB,   &
       &                     0.352E+01_JPRB,0.856E+00_JPRB,0.208E+00_JPRB,0.524E-01_JPRB,   &
       &                     0.145E-01_JPRB,0.466E-02_JPRB,0.179E-02_JPRB,0.860E-03_JPRB,   &
       &                     0.719E-03_JPRB,0.165E-02_JPRB,0.486E-02_JPRB,0.213E-02_JPRB,   &
       &                     0.554E-03_JPRB,0.891E-04_JPRB,0.211E-04_JPRB/)

      COLLEFF4(12,1:NROW)=(/ 0.117E+04_JPRB,0.326E+03_JPRB,0.807E+02_JPRB,0.200E+02_JPRB,   &
       &                     0.496E+01_JPRB,0.123E+01_JPRB,0.305E+00_JPRB,0.777E-01_JPRB,   &
       &                     0.219E-01_JPRB,0.720E-02_JPRB,0.281E-02_JPRB,0.137E-02_JPRB,   &
       &                     0.941E-03_JPRB,0.330E-02_JPRB,0.563E-02_JPRB,0.255E-02_JPRB,   &
       &                     0.686E-03_JPRB,0.116E-03_JPRB,0.283E-04_JPRB/)

      COLLEFF4(13,1:NROW)=(/ 0.130E+04_JPRB,0.371E+03_JPRB,0.938E+02_JPRB,0.237E+02_JPRB,   &
       &                     0.601E+01_JPRB,0.152E+01_JPRB,0.385E+00_JPRB,0.979E-01_JPRB,   &
       &                     0.276E-01_JPRB,0.926E-02_JPRB,0.364E-02_JPRB,0.173E-02_JPRB,   &
       &                     0.101E-02_JPRB,0.406E-02_JPRB,0.694E-02_JPRB,0.327E-02_JPRB,   &
       &                     0.930E-03_JPRB,0.167E-03_JPRB,0.429E-04_JPRB/)

      COLLEFF4(14,1:NROW)=(/ 0.118E+04_JPRB,0.333E+03_JPRB,0.842E+02_JPRB,0.213E+02_JPRB,   &
       &                     0.537E+01_JPRB,0.136E+01_JPRB,0.342E+00_JPRB,0.876E-01_JPRB,   &
       &                     0.250E-01_JPRB,0.841E-02_JPRB,0.330E-02_JPRB,0.153E-02_JPRB,   &
       &                     0.801E-03_JPRB,0.260E-02_JPRB,0.973E-02_JPRB,0.490E-02_JPRB,   &
       &                     0.152E-02_JPRB,0.303E-03_JPRB,0.842E-04_JPRB/)

      COLLEFF4(15,1:NROW)=(/ 0.774E+03_JPRB,0.223E+03_JPRB,0.572E+02_JPRB,0.147E+02_JPRB,   &
       &                     0.378E+01_JPRB,0.970E+00_JPRB,0.249E+00_JPRB,0.636E-01_JPRB,   &
       &                     0.180E-01_JPRB,0.606E-02_JPRB,0.238E-02_JPRB,0.110E-02_JPRB,   &
       &                     0.658E-03_JPRB,0.107E-02_JPRB,0.167E-01_JPRB,0.940E-02_JPRB,   &
       &                     0.335E-02_JPRB,0.791E-03_JPRB,0.249E-03_JPRB/)

      COLLEFF4(16,1:NROW)=(/ 0.372E+01_JPRB,0.177E+01_JPRB,0.781E+00_JPRB,0.345E+00_JPRB,   &
       &                     0.153E+00_JPRB,0.675E-01_JPRB,0.299E-01_JPRB,0.130E-01_JPRB,   &
       &                     0.624E-02_JPRB,0.346E-02_JPRB,0.179E-02_JPRB,0.142E-02_JPRB,   &
       &                     0.164E-02_JPRB,0.401E-02_JPRB,0.413E-01_JPRB,0.277E-01_JPRB,   &
       &                     0.124E-01_JPRB,0.389E-02_JPRB,0.152E-02_JPRB/)

      COLLEFF4(17,1:NROW)=(/ 0.234E-18_JPRB,0.108E-16_JPRB,0.705E-15_JPRB,0.462E-13_JPRB,   &
       &                     0.302E-11_JPRB,0.198E-09_JPRB,0.129E-07_JPRB,0.844E-06_JPRB,   &
       &                     0.568E-04_JPRB,0.386E-02_JPRB,0.286E-01_JPRB,0.448E-01_JPRB,   &
       &                     0.569E-01_JPRB,0.859E-01_JPRB,0.183E+00_JPRB,0.165E+00_JPRB,   &
       &                     0.108E+00_JPRB,0.536E-01_JPRB,0.295E-01_JPRB/)

      COLLEFF4(18,1:NROW)=(/ 0.902E-37_JPRB,0.794E-33_JPRB,0.160E-28_JPRB,0.324E-24_JPRB,   &
       &                     0.655E-20_JPRB,0.132E-15_JPRB,0.267E-11_JPRB,0.540E-07_JPRB,   &
       &                     0.816E-03_JPRB,0.482E-01_JPRB,0.245E+00_JPRB,0.372E+00_JPRB,   &
       &                     0.436E+00_JPRB,0.473E+00_JPRB,0.493E+00_JPRB,0.497E+00_JPRB,   &
       &                     0.414E+00_JPRB,0.299E+00_JPRB,0.225E+00_JPRB/)

      COLLEFF4(19,1:NROW)=(/ 0.275E-30_JPRB,0.136E-26_JPRB,0.146E-22_JPRB,0.156E-18_JPRB,   &
       &                     0.167E-14_JPRB,0.178E-10_JPRB,0.191E-06_JPRB,0.202E-02_JPRB,   &
       &                     0.203E+00_JPRB,0.427E+00_JPRB,0.586E+00_JPRB,0.669E+00_JPRB,   &
       &                     0.708E+00_JPRB,0.730E+00_JPRB,0.746E+00_JPRB,0.738E+00_JPRB,   &
       &                     0.679E+00_JPRB,0.588E+00_JPRB,0.520E+00_JPRB/)

      COLLEFF4(20,1:NROW)=(/ 0.136E-33_JPRB,0.238E-29_JPRB,0.102E-24_JPRB,0.436E-20_JPRB,   &
       &                     0.186E-15_JPRB,0.797E-11_JPRB,0.341E-06_JPRB,0.143E-01_JPRB,   &
       &                     0.722E+00_JPRB,0.805E+00_JPRB,0.869E+00_JPRB,0.902E+00_JPRB,   &
       &                     0.915E+00_JPRB,0.932E+00_JPRB,0.104E+01_JPRB,0.927E+00_JPRB,   &
       &                     0.904E+00_JPRB,0.871E+00_JPRB,0.842E+00_JPRB/)

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_IMSCAVCOFF',1,ZHOOK_HANDLE)

      END SUBROUTINE UKCA_MODE_IMSCAVCOFF
      SUBROUTINE UKCA_MODE_ALLCP_4MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components with all components
!+ switched on but only 4 modes used.
!
! History:
! Version    Date      Author         Description
! -----      --------  ---------      ------------
! vn6.1      1/XII/06  Graham Mann    Original code
!
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: AVC, RHOSUL, MMSUL, PPI
      IMPLICIT NONE

      INTEGER(KIND=JPIM) :: JMODE
      INTEGER(KIND=JPIM) :: JCP

      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_ALLCP_4MODE',0,ZHOOK_HANDLE)

! Mode names
      MODE_NAMES=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
       &                     'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      MODE_CHOICE=(/1,1,1,1,0,0,0/)
! Specify which modes are soluble
      MODESOL=(/1,1,1,1,0,0,0/)
! Component names
      COMPONENT_NAMES=                                                  &
     & (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      COMPONENT_CHOICE=(/1,1,1,1,1,0/)
! ***n.b. in above have kept all cpts on (not SO) for UM test***
! Components that are soluble
      SOLUBLE_CHOICE=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      COMPONENT_MODE(1,1:NCP)=(/1,0,0,0,0,1/)       !allowed in nuc_sol
      COMPONENT_MODE(2,1:NCP)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      COMPONENT_MODE(3,1:NCP)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      COMPONENT_MODE(4,1:NCP)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      COMPONENT_MODE(5,1:NCP)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      COMPONENT_MODE(6,1:NCP)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      COMPONENT_MODE(7,1:NCP)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!! set dplim34 here to be 500nm to match value found from bin-mode comparison
!!      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,1.0e-6,1.0e-8,1.0e-7,1.0e-6/)
!!      ddplim1=(/1.0e-8,1.0e-7,1.0e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)
      DDPLIM0=(/1.0E-9_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,1.0E-6_JPRB/)
      DDPLIM1=(/1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-5_JPRB,1.0E-7_JPRB,1.0E-6_JPRB,1.0E-5_JPRB/)

! Specify fixed geometric standard deviation for each mode
      SIGMAG=(/1.59_JPRB,1.59_JPRB,1.40_JPRB,2.0_JPRB,1.59_JPRB,1.59_JPRB,2.0_JPRB/) ! to match M7 but sigacc=1.4
!!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!
SIGMAG=(/1.59_JPRB,1.45_JPRB,1.4_JPRB,2.0_JPRB,1.45_JPRB,1.4_JPRB,2.0_JPRB/) ! to match UM scheme

      DO JMODE=1,NMODES
        X(JMODE)=EXP(4.5_JPRB*LOG(SIGMAG(JMODE))*LOG(SIGMAG(JMODE)))
      ENDDO
!
! Specify threshold for ND (per cc) below which don't do calculations
!      NUM_EPS=(/1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB/)
      NUM_EPS=(/1.0E-1,1.0E-2,1.0E-3,1.0E-4,1.0E-2,1.0E-3,1.0E-4/)
!       NUM_EPS=(/1.0e-1,1.0e-1,1.0e-1,1.0e-1,1.0e-1,1.0e-1,1.0e-1/)
!
      DO JMODE=1,NMODES
       DDPMID(JMODE)=                                                   &
        & EXP(0.5_JPRB*(LOG(DDPLIM0(JMODE))+LOG(DDPLIM1(JMODE))))
       MMID(JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPMID(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MLO (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM0(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MHI (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM1(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
      ENDDO

! Initial fractions of mass in each mode among components
      MFRAC_0(1,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !nucln. soluble
      MFRAC_0(2,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken soluble
      MFRAC_0(3,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !accum. soluble
      MFRAC_0(4,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/) !coarse soluble
      MFRAC_0(5,1:NCP)=(/0.0_JPRB,0.5_JPRB,0.5_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken insoluble
      MFRAC_0(6,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !accum. insoluble
      MFRAC_0(7,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !coarse insoluble

! Modes resulting when two modes coagulate
      COAG_MODE(1,1:NMODES)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      COAG_MODE(2,1:NMODES)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(3,1:NMODES)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(4,1:NMODES)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(5,1:NMODES)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(6,1:NMODES)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(7,1:NMODES)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      MM=(/0.098_JPRB,0.012_JPRB,0.0168_JPRB,0.05844_JPRB,0.100_JPRB,0.0168_JPRB/)
!          h2so4  bc     oc    nacl   dust   so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 C-H ratio)
! Mass density of components (kg m^-3)
      RHOCOMP=(/1769.0_JPRB,1500.0_JPRB,1500.0_JPRB,1600.0_JPRB,2650.0_JPRB,1500.0_JPRB/)
! Assume other components have same mass density as H2SO4
! number of dissociating ions in soluble components
      NO_IONS=(/3.0_JPRB,0.0_JPRB,0.0_JPRB,2.0_JPRB,0.0_JPRB,0.0_JPRB/)
!
      FRACBCEM=(/0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/)
      FRACOCEM=(/0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).
!
! Set logical variables
      MODE=(MODE_CHOICE > 0)
      COMPONENT=.FALSE.
      SOLUBLE=.FALSE.
      DO JMODE=1,NMODES
        DO JCP=1,NCP
          IF(((COMPONENT_MODE(JMODE,JCP) == 1).AND.                     &
           &    (COMPONENT_CHOICE(JCP) == 1)).AND.                      &
           &    (MODE_CHOICE(JMODE) == 1)) THEN
             COMPONENT(JMODE,JCP)=.TRUE.
          ENDIF
          IF(SOLUBLE_CHOICE(JCP) == 1) SOLUBLE(JCP)=.TRUE.
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_ALLCP_4MODE',1,ZHOOK_HANDLE)

      END SUBROUTINE UKCA_MODE_ALLCP_4MODE
      SUBROUTINE UKCA_MODE_SUSS_4MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ sulfate and sea-salt only in 4 modes.
!+ Uses 10 aerosol tracers
!
! History:
! Version    Date      Author         Description
! -----      --------  ---------      ------------
! vn6.1      1/XII/06  Graham Mann    Original code
!
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: AVC, RHOSUL, MMSUL, PPI
      IMPLICIT NONE

      INTEGER(KIND=JPIM) :: JMODE
      INTEGER(KIND=JPIM) :: JCP

      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_SUSS_4MODE',0,ZHOOK_HANDLE)

! Mode names
      MODE_NAMES=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
       &                     'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      MODE_CHOICE=(/1,1,1,1,0,0,0/)
! Specify which modes are soluble
      MODESOL=(/1,1,1,1,0,0,0/)
! Component names
      COMPONENT_NAMES=                                                  &
     & (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      COMPONENT_CHOICE=(/1,0,0,1,0,0/)
! *** n.b. only have h2so4 and nacl cpts on for this setup ***
! Components that are soluble
      SOLUBLE_CHOICE=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      COMPONENT_MODE(1,1:NCP)=(/1,0,0,0,0,1/)       !allowed in nuc_sol
      COMPONENT_MODE(2,1:NCP)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      COMPONENT_MODE(3,1:NCP)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      COMPONENT_MODE(4,1:NCP)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      COMPONENT_MODE(5,1:NCP)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      COMPONENT_MODE(6,1:NCP)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      COMPONENT_MODE(7,1:NCP)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!! set dplim34 here to be 500nm to match value found from bin-mode comparison
!!      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,1.0e-6,1.0e-8,1.0e-7,1.0e-6/)
!!      ddplim1=(/1.0e-8,1.0e-7,1.0e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)
      DDPLIM0=(/1.0E-9_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,1.0E-6_JPRB/)
      DDPLIM1=(/1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-5_JPRB,1.0E-7_JPRB,1.0E-6_JPRB,1.0E-5_JPRB/)

! Specify fixed geometric standard deviation for each mode
      SIGMAG=(/1.59_JPRB,1.59_JPRB,1.40_JPRB,2.0_JPRB,1.59_JPRB,1.59_JPRB,2.0_JPRB/) ! to match M7 but sigacc=1.4
!!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO JMODE=1,NMODES
        X(JMODE)=EXP(4.5_JPRB*LOG(SIGMAG(JMODE))*LOG(SIGMAG(JMODE)))
      ENDDO
!
! Specify threshold for ND (per cc) below which don't do calculations
!      NUM_EPS=(/1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB/)
      !NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-4,1.0e-3,1.0e-3,1.0e-4/)
      NUM_EPS=(/1.0E-1,1.0E-2,1.0E-3,1.0E-4,1.0E-2,1.0E-3,1.0E-4/)
!!      NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3/)
!       NUM_EPS=(/1.0e-1,1.0e-1,1.0e-1,1.0e-1,1.0e-1,1.0e-1,1.0e-1/)
!
      DO JMODE=1,NMODES
       DDPMID(JMODE)=                                                   &
        & EXP(0.5_JPRB*(LOG(DDPLIM0(JMODE))+LOG(DDPLIM1(JMODE))))
       MMID(JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPMID(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MLO (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM0(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MHI (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM1(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
      ENDDO

! Initial fractions of mass in each mode among components
      MFRAC_0(1,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !nucln. soluble
      MFRAC_0(2,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken soluble
      MFRAC_0(3,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !accum. soluble
      MFRAC_0(4,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/) !coarse soluble
      MFRAC_0(5,1:NCP)=(/0.0_JPRB,0.5_JPRB,0.5_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken insoluble
      MFRAC_0(6,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !accum. insoluble
      MFRAC_0(7,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !coarse insoluble

! Modes resulting when two modes coagulate
      COAG_MODE(1,1:NMODES)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      COAG_MODE(2,1:NMODES)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(3,1:NMODES)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(4,1:NMODES)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(5,1:NMODES)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(6,1:NMODES)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(7,1:NMODES)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      MM=(/0.098_JPRB,0.012_JPRB,0.0168_JPRB,0.05844_JPRB,0.100_JPRB,0.0168_JPRB/)
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 C-H ratio)
! Mass density of components (kg m^-3)
      RHOCOMP=(/1769.0_JPRB,1500.0_JPRB,1500.0_JPRB,1600.0_JPRB,2650.0_JPRB,1500.0_JPRB/)
! Assume other components have same mass density as H2SO4
! number of dissociating ions in soluble components
      NO_IONS=(/3.0_JPRB,0.0_JPRB,0.0_JPRB,2.0_JPRB,0.0_JPRB,0.0_JPRB/)
!
      FRACBCEM=(/0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/)
      FRACOCEM=(/0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).
!
! Set logical variables
      MODE=(MODE_CHOICE > 0)
      COMPONENT=.FALSE.
      SOLUBLE=.FALSE.
      DO JMODE=1,NMODES
        DO JCP=1,NCP
          IF(((COMPONENT_MODE(JMODE,JCP) == 1).AND.                     &
           &    (COMPONENT_CHOICE(JCP) == 1)).AND.                      &
           &    (MODE_CHOICE(JMODE) == 1)) THEN
             COMPONENT(JMODE,JCP)=.TRUE.
          ENDIF
          IF(SOLUBLE_CHOICE(JCP) == 1) SOLUBLE(JCP)=.TRUE.
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_SUSS_4MODE',1,ZHOOK_HANDLE)

      END SUBROUTINE UKCA_MODE_SUSS_4MODE
      SUBROUTINE UKCA_MODE_SUSSBCOCDU_4MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ SO4, sea-salt, bc, oc (secondary & primary combined) & du in 4 modes.
!+ Uses 19 aerosol tracers
!
! History:
! Version    Date      Author         Description
! -----      --------  ---------      ------------
! vn6.1      1/XII/06  Graham Mann    Original code
!
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: AVC, RHOSUL, MMSUL, PPI
      IMPLICIT NONE

      INTEGER(KIND=JPIM) :: JMODE
      INTEGER(KIND=JPIM) :: JCP

      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCDU_4MODE',0,ZHOOK_HANDLE)

! Mode names
      MODE_NAMES=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
       &                     'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      MODE_CHOICE=(/1,1,1,1,0,0,0/)
! Specify which modes are soluble
      MODESOL=(/1,1,1,1,0,0,0/)
! Component names
      COMPONENT_NAMES=                                                  &
     & (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      COMPONENT_CHOICE=(/1,1,1,1,1,0/)
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
      SOLUBLE_CHOICE=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      COMPONENT_MODE(1,1:NCP)=(/1,0,1,0,0,1/)       !allowed in nuc_sol
      COMPONENT_MODE(2,1:NCP)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      COMPONENT_MODE(3,1:NCP)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      COMPONENT_MODE(4,1:NCP)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      COMPONENT_MODE(5,1:NCP)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      COMPONENT_MODE(6,1:NCP)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      COMPONENT_MODE(7,1:NCP)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!! set dplim34 here to be 500nm to match value found from bin-mode comparison
!!      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,1.0e-6,1.0e-8,1.0e-7,1.0e-6/)
!!      ddplim1=(/1.0e-8,1.0e-7,1.0e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)
      DDPLIM0=(/1.0E-9_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,1.0E-6_JPRB/)
      DDPLIM1=(/1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-5_JPRB,1.0E-7_JPRB,1.0E-6_JPRB,1.0E-5_JPRB/)

! Specify fixed geometric standard deviation for each mode
      SIGMAG=(/1.59_JPRB,1.59_JPRB,1.40_JPRB,2.0_JPRB,1.59_JPRB,1.59_JPRB,2.0_JPRB/) ! to match M7 but sigacc=1.4
!!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO JMODE=1,NMODES
        X(JMODE)=EXP(4.5_JPRB*LOG(SIGMAG(JMODE))*LOG(SIGMAG(JMODE)))
      ENDDO
!
! Specify threshold for ND (per cc) below which don't do calculations
!      NUM_EPS=(/1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB/)
      NUM_EPS=(/1.0E-1,1.0E-2,1.0E-3,1.0E-4,1.0E-2,1.0E-3,1.0E-4/)
     ! NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-4,1.0e-3,1.0e-3,1.0e-4/)
!!      NUM_EPS=(/1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8/)
!
      DO JMODE=1,NMODES
       DDPMID(JMODE)=                                                   &
        & EXP(0.5_JPRB*(LOG(DDPLIM0(JMODE))+LOG(DDPLIM1(JMODE))))
       MMID(JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPMID(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MLO (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM0(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MHI (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM1(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
      ENDDO

! Initial fractions of mass in each mode among components
      MFRAC_0(1,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !nucln. soluble
      MFRAC_0(2,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken soluble
      MFRAC_0(3,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !accum. soluble
      MFRAC_0(4,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/) !coarse soluble
      MFRAC_0(5,1:NCP)=(/0.0_JPRB,0.5_JPRB,0.5_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken insoluble
      MFRAC_0(6,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !accum. insoluble
      MFRAC_0(7,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !coarse insoluble

! Modes resulting when two modes coagulate
      COAG_MODE(1,1:NMODES)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      COAG_MODE(2,1:NMODES)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(3,1:NMODES)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(4,1:NMODES)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(5,1:NMODES)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(6,1:NMODES)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(7,1:NMODES)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      MM=(/0.098_JPRB,0.012_JPRB,0.0168_JPRB,0.05844_JPRB,0.100_JPRB,0.0168_JPRB/)
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 C-H ratio)
! Mass density of components (kg m^-3)
      RHOCOMP=(/1769.0_JPRB,1500.0_JPRB,1500.0_JPRB,1600.0_JPRB,2650.0_JPRB,1500.0_JPRB/)
! Assume other components have same mass density as H2SO4
! number of dissociating ions in soluble components
      NO_IONS=(/3.0_JPRB,0.0_JPRB,0.0_JPRB,2.0_JPRB,0.0_JPRB,0.0_JPRB/)
!
      FRACBCEM=(/0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/)
      FRACOCEM=(/0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).
!
! Set logical variables
      MODE=(MODE_CHOICE > 0)
      COMPONENT=.FALSE.
      SOLUBLE=.FALSE.
      DO JMODE=1,NMODES
        DO JCP=1,NCP
          IF(((COMPONENT_MODE(JMODE,JCP) == 1).AND.                     &
           &    (COMPONENT_CHOICE(JCP) == 1)).AND.                      &
           &    (MODE_CHOICE(JMODE) == 1)) THEN
             COMPONENT(JMODE,JCP)=.TRUE.
          ENDIF
          IF(SOLUBLE_CHOICE(JCP) == 1) SOLUBLE(JCP)=.TRUE.
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCDU_4MODE',1,ZHOOK_HANDLE)

      END SUBROUTINE UKCA_MODE_SUSSBCOCDU_4MODE
      SUBROUTINE UKCA_MODE_SUSSBCOCDU_7MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ SO4, sea-salt, bc, oc (secondary & primary combined) & du in 7 modes.
!+ Uses 26 aerosol tracers
!
! History:
! Version    Date      Author         Description
! -----      --------  ---------      ------------
! vn6.1      1/XII/06  Graham Mann    Original code
!
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: AVC, RHOSUL, MMSUL, PPI
      IMPLICIT NONE

      INTEGER(KIND=JPIM) :: JMODE
      INTEGER(KIND=JPIM) :: JCP

      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCDU_7MODE',0,ZHOOK_HANDLE)

! Mode names
      MODE_NAMES=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
       &                     'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      MODE_CHOICE=(/1,1,1,1,1,1,1/)
! Specify which modes are soluble
      MODESOL=(/1,1,1,1,0,0,0/)
! Component names
      COMPONENT_NAMES=                                                  &
     & (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      COMPONENT_CHOICE=(/1,1,1,1,1,0/)
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
      SOLUBLE_CHOICE=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      COMPONENT_MODE(1,1:NCP)=(/1,0,1,0,0,1/)       !allowed in nuc_sol
      COMPONENT_MODE(2,1:NCP)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      COMPONENT_MODE(3,1:NCP)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      COMPONENT_MODE(4,1:NCP)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      COMPONENT_MODE(5,1:NCP)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      COMPONENT_MODE(6,1:NCP)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      COMPONENT_MODE(7,1:NCP)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!! set dplim34 here to be 500nm to match value found from bin-mode comparison
!!      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,1.0e-6,1.0e-8,1.0e-7,1.0e-6/)
!!      ddplim1=(/1.0e-8,1.0e-7,1.0e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)
      DDPLIM0=(/1.0E-9_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,1.0E-6_JPRB/)
      DDPLIM1=(/1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-5_JPRB,1.0E-7_JPRB,1.0E-6_JPRB,1.0E-5_JPRB/)

! Specify fixed geometric standard deviation for each mode
      SIGMAG=(/1.59_JPRB,1.59_JPRB,1.40_JPRB,2.0_JPRB,1.59_JPRB,1.59_JPRB,2.0_JPRB/) ! to match M7 but sigacc=1.4
!!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO JMODE=1,NMODES
        X(JMODE)=EXP(4.5_JPRB*LOG(SIGMAG(JMODE))*LOG(SIGMAG(JMODE)))
      ENDDO
!
! Specify threshold for ND (per cc) below which don't do calculations
!      NUM_EPS=(/1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB/)
      !NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-4,1.0e-3,1.0e-3,1.0e-4/)
      NUM_EPS=(/1.0E-1,1.0E-2,1.0E-3,1.0E-4,1.0E-2,1.0E-3,1.0E-4/)
!!      NUM_EPS=(/1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8/)
!
      DO JMODE=1,NMODES
       DDPMID(JMODE)=                                                   &
        & EXP(0.5_JPRB*(LOG(DDPLIM0(JMODE))+LOG(DDPLIM1(JMODE))))
       MMID(JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPMID(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MLO (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM0(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MHI (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM1(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
      ENDDO

! Initial fractions of mass in each mode among components
      MFRAC_0(1,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !nucln. soluble
      MFRAC_0(2,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken soluble
      MFRAC_0(3,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !accum. soluble
      MFRAC_0(4,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/) !coarse soluble
      MFRAC_0(5,1:NCP)=(/0.0_JPRB,0.5_JPRB,0.5_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken insoluble
      MFRAC_0(6,1:NCP)=(/0.0_JPRB,0.5_JPRB,0.5_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !accum. insoluble
      MFRAC_0(7,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !coarse insoluble

! Modes resulting when two modes coagulate
      COAG_MODE(1,1:NMODES)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      COAG_MODE(2,1:NMODES)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(3,1:NMODES)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(4,1:NMODES)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(5,1:NMODES)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(6,1:NMODES)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(7,1:NMODES)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      MM=(/0.098_JPRB,0.012_JPRB,0.0168_JPRB,0.05844_JPRB,0.100_JPRB,0.0168_JPRB/)
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 C-H ratio)
! Mass density of components (kg m^-3)
      RHOCOMP=(/1769.0_JPRB,1500.0_JPRB,1500.0_JPRB,1600.0_JPRB,2650.0_JPRB,1500.0_JPRB/)
! Assume other components have same mass density as H2SO4
! number of dissociating ions in soluble components
      NO_IONS=(/3.0_JPRB,0.0_JPRB,0.0_JPRB,2.0_JPRB,0.0_JPRB,0.0_JPRB/)
!
      FRACBCEM=(/0.0_JPRB,0.2_JPRB,0.0_JPRB,0.0_JPRB,0.8_JPRB,0.0_JPRB,0.0_JPRB/)
      FRACOCEM=(/0.0_JPRB,0.5_JPRB,0.0_JPRB,0.0_JPRB,0.5_JPRB,0.0_JPRB,0.0_JPRB/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
!
! Set logical variables
      MODE=(MODE_CHOICE > 0)
      COMPONENT=.FALSE.
      SOLUBLE=.FALSE.
      DO JMODE=1,NMODES
        DO JCP=1,NCP
          IF(((COMPONENT_MODE(JMODE,JCP) == 1).AND.                     &
           &    (COMPONENT_CHOICE(JCP) == 1)).AND.                      &
           &    (MODE_CHOICE(JMODE) == 1)) THEN
             COMPONENT(JMODE,JCP)=.TRUE.
          ENDIF
          IF(SOLUBLE_CHOICE(JCP) == 1) SOLUBLE(JCP)=.TRUE.
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCDU_7MODE',1,ZHOOK_HANDLE)

      END SUBROUTINE UKCA_MODE_SUSSBCOCDU_7MODE
      SUBROUTINE UKCA_MODE_SUSSBCOC_4MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ sulfate, sea-salt, bc & oc (secondary & primary combined) in 4 modes.
!+ Uses 17 aerosol tracers
!
! History:
! Version    Date      Author         Description
! -----      --------  ---------      ------------
! vn6.1      1/XII/06  Graham Mann    Original code
!
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: AVC, RHOSUL, MMSUL, PPI
      IMPLICIT NONE

      INTEGER(KIND=JPIM) :: JMODE
      INTEGER(KIND=JPIM) :: JCP

      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOC_4MODE',0,ZHOOK_HANDLE)

! Mode names
      MODE_NAMES=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
       &                     'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      MODE_CHOICE=(/1,1,1,1,0,0,0/)
! Specify which modes are soluble
      MODESOL=(/1,1,1,1,0,0,0/)
! Component names
      COMPONENT_NAMES=                                                  &
     & (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      COMPONENT_CHOICE=(/1,1,1,1,0,0/)
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
      SOLUBLE_CHOICE=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      COMPONENT_MODE(1,1:NCP)=(/1,0,1,0,0,1/)       !allowed in nuc_sol
      COMPONENT_MODE(2,1:NCP)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      COMPONENT_MODE(3,1:NCP)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      COMPONENT_MODE(4,1:NCP)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      COMPONENT_MODE(5,1:NCP)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      COMPONENT_MODE(6,1:NCP)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      COMPONENT_MODE(7,1:NCP)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!! set dplim34 here to be 500nm to match value found from bin-mode comparison
!!      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,1.0e-6,1.0e-8,1.0e-7,1.0e-6/)
!!      ddplim1=(/1.0e-8,1.0e-7,1.0e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)
      DDPLIM0=(/1.0E-9_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,1.0E-6_JPRB/)
      DDPLIM1=(/1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-5_JPRB,1.0E-7_JPRB,1.0E-6_JPRB,1.0E-5_JPRB/)

! Specify fixed geometric standard deviation for each mode
      SIGMAG=(/1.59_JPRB,1.59_JPRB,1.40_JPRB,2.0_JPRB,1.59_JPRB,1.59_JPRB,2.0_JPRB/) ! to match M7 but sigacc=1.4
!!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO JMODE=1,NMODES
        X(JMODE)=EXP(4.5*LOG(SIGMAG(JMODE))*LOG(SIGMAG(JMODE)))
      ENDDO
!
! Specify threshold for ND (per cc) below which don't do calculations
!      NUM_EPS=(/1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB/)
      !NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-4,1.0e-3,1.0e-3,1.0e-4/)
      NUM_EPS=(/1.0E-1,1.0E-2,1.0E-3,1.0E-4,1.0E-2,1.0E-3,1.0E-4/)
!!      NUM_EPS=(/1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8/)
!!      NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3/)
!
      DO JMODE=1,NMODES
       DDPMID(JMODE)=                                                   &
        & EXP(0.5_JPRB*(LOG(DDPLIM0(JMODE))+LOG(DDPLIM1(JMODE))))
       MMID(JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPMID(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MLO (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM0(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MHI (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM1(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
      ENDDO

! Initial fractions of mass in each mode among components
      MFRAC_0(1,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !nucln. soluble
      MFRAC_0(2,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken soluble
      MFRAC_0(3,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !accum. soluble
      MFRAC_0(4,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/) !coarse soluble
      MFRAC_0(5,1:NCP)=(/0.0_JPRB,0.5_JPRB,0.5_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken insoluble
      MFRAC_0(6,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !accum. insoluble
      MFRAC_0(7,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !coarse insoluble

! Modes resulting when two modes coagulate
      COAG_MODE(1,1:NMODES)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      COAG_MODE(2,1:NMODES)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(3,1:NMODES)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(4,1:NMODES)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(5,1:NMODES)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(6,1:NMODES)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(7,1:NMODES)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      MM=(/0.098_JPRB,0.012_JPRB,0.0168_JPRB,0.05844_JPRB,0.100_JPRB,0.0168_JPRB/)
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 C-H ratio)
! Mass density of components (kg m^-3)
      RHOCOMP=(/1769.0_JPRB,1500.0_JPRB,1500.0_JPRB,1600.0_JPRB,2650.0_JPRB,1500.0_JPRB/)
! Assume other components have same mass density as H2SO4
! number of dissociating ions in soluble components
      NO_IONS=(/3.0_JPRB,0.0_JPRB,0.0_JPRB,2.0_JPRB,0.0_JPRB,0.0_JPRB/)
!
      FRACBCEM=(/0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/)
      FRACOCEM=(/0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).
!
! Set logical variables
      MODE=(MODE_CHOICE > 0)
      COMPONENT=.FALSE.
      SOLUBLE=.FALSE.
      DO JMODE=1,NMODES
        DO JCP=1,NCP
          IF(((COMPONENT_MODE(JMODE,JCP) == 1).AND.                     &
           &    (COMPONENT_CHOICE(JCP) == 1)).AND.                      &
           &    (MODE_CHOICE(JMODE) == 1)) THEN
             COMPONENT(JMODE,JCP)=.TRUE.
          ENDIF
          IF(SOLUBLE_CHOICE(JCP) == 1) SOLUBLE(JCP)=.TRUE.
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOC_4MODE',1,ZHOOK_HANDLE)

      END SUBROUTINE UKCA_MODE_SUSSBCOC_4MODE
      SUBROUTINE UKCA_MODE_SUSSBCOC_5MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ sulfate, sea-salt, bc & oc (secondary & primary combined) in 5 modes.
!+ Uses 20 aerosol tracers
!
! History:
! Version    Date      Author         Description
! -----      --------  ---------      ------------
! vn6.1      1/XII/06  Graham Mann    Original code
!
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: AVC, RHOSUL, MMSUL, PPI
      IMPLICIT NONE

      INTEGER(KIND=JPIM) :: JMODE
      INTEGER(KIND=JPIM) :: JCP

      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOC_5MODE',0,ZHOOK_HANDLE)

! Mode names
      MODE_NAMES=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
       &                     'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      MODE_CHOICE=(/1,1,1,1,1,0,0/)
! Specify which modes are soluble
      MODESOL=(/1,1,1,1,0,0,0/)
! Component names
      COMPONENT_NAMES=                                                  &
     & (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      COMPONENT_CHOICE=(/1,1,1,1,0,0/)
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
      SOLUBLE_CHOICE=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      COMPONENT_MODE(1,1:NCP)=(/1,0,1,0,0,1/)       !allowed in nuc_sol
      COMPONENT_MODE(2,1:NCP)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      COMPONENT_MODE(3,1:NCP)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      COMPONENT_MODE(4,1:NCP)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      COMPONENT_MODE(5,1:NCP)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      COMPONENT_MODE(6,1:NCP)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      COMPONENT_MODE(7,1:NCP)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!! set dplim34 here to be 500nm to match value found from bin-mode comparison
!!      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,1.0e-6,1.0e-8,1.0e-7,1.0e-6/)
!!      ddplim1=(/1.0e-8,1.0e-7,1.0e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)
      DDPLIM0=(/1.0E-9_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,1.0E-6_JPRB/)
      DDPLIM1=(/1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-5_JPRB,1.0E-7_JPRB,1.0E-6_JPRB,1.0E-5_JPRB/)

! Specify fixed geometric standard deviation for each mode
      SIGMAG=(/1.59_JPRB,1.59_JPRB,1.40_JPRB,2.0_JPRB,1.59_JPRB,1.59_JPRB,2.0_JPRB/) ! to match M7 but sigacc=1.4
!!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO JMODE=1,NMODES
        X(JMODE)=EXP(4.5_JPRB*LOG(SIGMAG(JMODE))*LOG(SIGMAG(JMODE)))
      ENDDO
!
! Specify threshold for ND (per cc) below which don't do calculations
!      NUM_EPS=(/1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB/)
      NUM_EPS=(/1.0E-1,1.0E-2,1.0E-3,1.0E-4,1.0E-2,1.0E-3,1.0E-4/)
     ! NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-4,1.0e-3,1.0e-3,1.0e-4/)
!!      NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3/)
!
      DO JMODE=1,NMODES
       DDPMID(JMODE)=                                                   &
        & EXP(0.5_JPRB*(LOG(DDPLIM0(JMODE))+LOG(DDPLIM1(JMODE))))
       MMID(JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPMID(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MLO (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM0(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MHI (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM1(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
      ENDDO

! Initial fractions of mass in each mode among components
      MFRAC_0(1,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !nucln. soluble
      MFRAC_0(2,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken soluble
      MFRAC_0(3,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !accum. soluble
      MFRAC_0(4,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/) !coarse soluble
      MFRAC_0(5,1:NCP)=(/0.0_JPRB,0.5_JPRB,0.5_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken insoluble
      MFRAC_0(6,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !accum. insoluble
      MFRAC_0(7,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !coarse insoluble

! Modes resulting when two modes coagulate
      COAG_MODE(1,1:NMODES)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      COAG_MODE(2,1:NMODES)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(3,1:NMODES)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(4,1:NMODES)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(5,1:NMODES)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(6,1:NMODES)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(7,1:NMODES)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      MM=(/0.098_JPRB,0.012_JPRB,0.0168_JPRB,0.05844_JPRB,0.100_JPRB,0.0168_JPRB/)
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 C-H ratio)
! Mass density of components (kg m^-3)
      RHOCOMP=(/1769.0_JPRB,1500.0_JPRB,1500.0_JPRB,1600.0_JPRB,2650.0_JPRB,1500.0_JPRB/)
! Assume other components have same mass density as H2SO4
! number of dissociating ions in soluble components
      NO_IONS=(/3.0_JPRB,0.0_JPRB,0.0_JPRB,2.0_JPRB,0.0_JPRB,0.0_JPRB/)
!
      FRACBCEM=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/)
      FRACOCEM=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
!
! Set logical variables
      MODE=(MODE_CHOICE > 0)
      COMPONENT=.FALSE.
      SOLUBLE=.FALSE.
      DO JMODE=1,NMODES
        DO JCP=1,NCP
          IF(((COMPONENT_MODE(JMODE,JCP) == 1).AND.                     &
           &    (COMPONENT_CHOICE(JCP) == 1)).AND.                      &
           &    (MODE_CHOICE(JMODE) == 1)) THEN
             COMPONENT(JMODE,JCP)=.TRUE.
          ENDIF
          IF(SOLUBLE_CHOICE(JCP) == 1) SOLUBLE(JCP)=.TRUE.
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOC_5MODE',1,ZHOOK_HANDLE)

      END SUBROUTINE UKCA_MODE_SUSSBCOC_5MODE
      SUBROUTINE UKCA_MODE_SUSSBCOCSO_4MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ sulfate, sea-salt, bc, primary oc & secondary oc cpts in 5 modes.
!+ Uses 20 aerosol tracers
!
! History:
! Version    Date      Author         Description
! -----      --------  ---------      ------------
! vn6.1      1/XII/06  Graham Mann    Original code
!
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: AVC, RHOSUL, MMSUL, PPI
      IMPLICIT NONE

      INTEGER(KIND=JPIM) :: JMODE
      INTEGER(KIND=JPIM) :: JCP

      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCSO_4MODE',0,ZHOOK_HANDLE)

! Mode names
      MODE_NAMES=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
       &                     'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      MODE_CHOICE=(/1,1,1,1,0,0,0/)
! Specify which modes are soluble
      MODESOL=(/1,1,1,1,0,0,0/)
! Component names
      COMPONENT_NAMES=                                                  &
     & (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      COMPONENT_CHOICE=(/1,1,1,1,0,1/) ! ***all cpts on except dust***
! Components that are soluble
      SOLUBLE_CHOICE=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      COMPONENT_MODE(1,1:NCP)=(/1,0,0,0,0,1/)       !allowed in nuc_sol
      COMPONENT_MODE(2,1:NCP)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      COMPONENT_MODE(3,1:NCP)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      COMPONENT_MODE(4,1:NCP)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      COMPONENT_MODE(5,1:NCP)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      COMPONENT_MODE(6,1:NCP)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      COMPONENT_MODE(7,1:NCP)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!! set dplim34 here to be 500nm to match value found from bin-mode comparison
!!      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,1.0e-6,1.0e-8,1.0e-7,1.0e-6/)
!!      ddplim1=(/1.0e-8_JPRB,1.0e-7_JPRB,1.0e-6_JPRB,1.0e-5_JPRB,1.0e-7_JPRB,1.0e-6_JPRB,1.0e-5/)
      DDPLIM0=(/1.0E-9_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,1.0E-6_JPRB/)
      DDPLIM1=(/1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-5_JPRB,1.0E-7_JPRB,1.0E-6_JPRB,1.0E-5_JPRB/)

! Specify fixed geometric standard deviation for each mode
      SIGMAG=(/1.59_JPRB,1.59_JPRB,1.40_JPRB,2.0_JPRB,1.59_JPRB,1.59_JPRB,2.0_JPRB/) ! to match M7 but sigacc=1.4
!!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO JMODE=1,NMODES
        X(JMODE)=EXP(4.5_JPRB*LOG(SIGMAG(JMODE))*LOG(SIGMAG(JMODE)))
      ENDDO
!
! Specify threshold for ND (per cc) below which don't do calculations
!!      NUM_EPS=(/1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB/)
      !NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-4,1.0e-3,1.0e-3,1.0e-4/)
      NUM_EPS=(/1.0E-1,1.0E-2,1.0E-3,1.0E-4,1.0E-2,1.0E-3,1.0E-4/)
!!      NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3/)
!
      DO JMODE=1,NMODES
       DDPMID(JMODE)=                                                   &
        & EXP(0.5_JPRB*(LOG(DDPLIM0(JMODE))+LOG(DDPLIM1(JMODE))))
       MMID(JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPMID(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MLO (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM0(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MHI (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM1(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
      ENDDO

! Initial fractions of mass in each mode among components
      MFRAC_0(1,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !nucln. soluble
      MFRAC_0(2,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken soluble
      MFRAC_0(3,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !accum. soluble
      MFRAC_0(4,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/) !coarse soluble
      MFRAC_0(5,1:NCP)=(/0.0_JPRB,0.5_JPRB,0.5_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken insoluble
      MFRAC_0(6,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !accum. insoluble
      MFRAC_0(7,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !coarse insoluble

! Modes resulting when two modes coagulate
      COAG_MODE(1,1:NMODES)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      COAG_MODE(2,1:NMODES)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(3,1:NMODES)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(4,1:NMODES)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(5,1:NMODES)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(6,1:NMODES)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(7,1:NMODES)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      MM=(/0.098_JPRB,0.012_JPRB,0.0168_JPRB,0.05844_JPRB,0.100_JPRB,0.0168_JPRB/)
!          h2so4  bc     oc    nacl   dust  sec_org
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 C-H ratio)
! assume secondary organic species has mm=mm_oc
! Mass density of components (kg m^-3)
      RHOCOMP=(/1769.0_JPRB,1500.0_JPRB,1500.0_JPRB,1600.0_JPRB,2650.0_JPRB,1500.0_JPRB/)
! number of dissociating ions in soluble components
      NO_IONS=(/3.0_JPRB,0.0_JPRB,0.0_JPRB,2.0_JPRB,0.0_JPRB,0.0_JPRB/)
!
      FRACBCEM=(/0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/)
      FRACOCEM=(/0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into   soluble Aitken for this setup).
!
! Set logical variables
      MODE=(MODE_CHOICE > 0)
      COMPONENT=.FALSE.
      SOLUBLE=.FALSE.
      DO JMODE=1,NMODES
        DO JCP=1,NCP
          IF(((COMPONENT_MODE(JMODE,JCP) == 1).AND.                     &
           &    (COMPONENT_CHOICE(JCP) == 1)).AND.                      &
           &    (MODE_CHOICE(JMODE) == 1)) THEN
             COMPONENT(JMODE,JCP)=.TRUE.
          ENDIF
          IF(SOLUBLE_CHOICE(JCP) == 1) SOLUBLE(JCP)=.TRUE.
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCSO_4MODE',1,ZHOOK_HANDLE)

      END SUBROUTINE UKCA_MODE_SUSSBCOCSO_4MODE
      SUBROUTINE UKCA_MODE_SUSSBCOCSO_5MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ sulfate, sea-salt, bc, primary oc & secondary oc cpts in 5 modes.
!+ Uses 23 aerosol tracers
!
! History:
! Version    Date      Author         Description
! -----      --------  ---------      ------------
! vn6.1      1/XII/06  Graham Mann    Original code
!
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: AVC, RHOSUL, MMSUL, PPI
      IMPLICIT NONE

      INTEGER(KIND=JPIM) :: JMODE
      INTEGER(KIND=JPIM) :: JCP

      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCSO_5MODE',0,ZHOOK_HANDLE)

! Mode names
      MODE_NAMES=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
       &                     'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      MODE_CHOICE=(/1,1,1,1,1,0,0/)
! Specify which modes are soluble
      MODESOL=(/1,1,1,1,0,0,0/)
! Component names
      COMPONENT_NAMES=                                                  &
     & (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      COMPONENT_CHOICE=(/1,1,1,1,0,1/) ! ***all cpts on except dust***
! Components that are soluble
      SOLUBLE_CHOICE=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      COMPONENT_MODE(1,1:NCP)=(/1,0,0,0,0,1/)       !allowed in nuc_sol
      COMPONENT_MODE(2,1:NCP)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      COMPONENT_MODE(3,1:NCP)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      COMPONENT_MODE(4,1:NCP)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      COMPONENT_MODE(5,1:NCP)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      COMPONENT_MODE(6,1:NCP)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      COMPONENT_MODE(7,1:NCP)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!! set dplim34 here to be 500nm to match value found from bin-mode comparison
!!      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,1.0e-6,1.0e-8,1.0e-7,1.0e-6/)
!!      ddplim1=(/1.0e-8,1.0e-7,1.0e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)
      DDPLIM0=(/1.0E-9_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,1.0E-6_JPRB/)
      DDPLIM1=(/1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-5_JPRB,1.0E-7_JPRB,1.0E-6_JPRB,1.0E-5_JPRB/)

! Specify fixed geometric standard deviation for each mode
      SIGMAG=(/1.59_JPRB,1.59_JPRB,1.40_JPRB,2.0_JPRB,1.59_JPRB,1.59_JPRB,2.0_JPRB/) ! to match M7 but sigacc=1.4
!!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO JMODE=1,NMODES
        X(JMODE)=EXP(4.5_JPRB*LOG(SIGMAG(JMODE))*LOG(SIGMAG(JMODE)))
      ENDDO
!
! Specify threshold for ND (per cc) below which don't do calculations
!!      NUM_EPS=(/1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB/)
      NUM_EPS=(/1.0E-1,1.0E-2,1.0E-3,1.0E-4,1.0E-2,1.0E-3,1.0E-4/)
      !NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-4,1.0e-3,1.0e-3,1.0e-4/)
!!      NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3/)
!
      DO JMODE=1,NMODES
       DDPMID(JMODE)=                                                   &
        & EXP(0.5_JPRB*(LOG(DDPLIM0(JMODE))+LOG(DDPLIM1(JMODE))))
       MMID(JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPMID(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MLO (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM0(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MHI (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM1(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
      ENDDO

! Initial fractions of mass in each mode among components
      MFRAC_0(1,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !nucln. soluble
      MFRAC_0(2,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken soluble
      MFRAC_0(3,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !accum. soluble
      MFRAC_0(4,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/) !coarse soluble
      MFRAC_0(5,1:NCP)=(/0.0_JPRB,0.5_JPRB,0.5_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken insoluble
      MFRAC_0(6,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !accum. insoluble
      MFRAC_0(7,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !coarse insoluble

! Modes resulting when two modes coagulate
      COAG_MODE(1,1:NMODES)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      COAG_MODE(2,1:NMODES)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(3,1:NMODES)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(4,1:NMODES)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(5,1:NMODES)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(6,1:NMODES)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(7,1:NMODES)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      MM=(/0.098_JPRB,0.012_JPRB,0.0168_JPRB,0.05844_JPRB,0.100_JPRB,0.0168_JPRB/)
!          h2so4  bc     oc    nacl   dust  sec_org
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 C-H ratio)
! assume secondary organic species has mm=mm_oc
! Mass density of components (kg m^-3)
      RHOCOMP=(/1769.0_JPRB,1500.0_JPRB,1500.0_JPRB,1600.0_JPRB,2650.0_JPRB,1500.0_JPRB/)
! number of dissociating ions in soluble components
      NO_IONS=(/3.0_JPRB,0.0_JPRB,0.0_JPRB,2.0_JPRB,0.0_JPRB,0.0_JPRB/)
!
      FRACBCEM=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/)
      FRACOCEM=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
!
! Set logical variables
      MODE=(MODE_CHOICE > 0)
      COMPONENT=.FALSE.
      SOLUBLE=.FALSE.
      DO JMODE=1,NMODES
        DO JCP=1,NCP
          IF(((COMPONENT_MODE(JMODE,JCP) == 1).AND.                     &
           &    (COMPONENT_CHOICE(JCP) == 1)).AND.                      &
           &    (MODE_CHOICE(JMODE) == 1)) THEN
             COMPONENT(JMODE,JCP)=.TRUE.
          ENDIF
          IF(SOLUBLE_CHOICE(JCP) == 1) SOLUBLE(JCP)=.TRUE.
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCSO_5MODE',1,ZHOOK_HANDLE)

      END SUBROUTINE UKCA_MODE_SUSSBCOCSO_5MODE
      SUBROUTINE UKCA_MODE_DUONLY_2MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ only du cpt in 2 (insoluble) modes.
!+ Uses  4 aerosol tracers
!
! History:
! Version    Date      Author         Description
! -----      --------  ---------      ------------
! vn6.?      ?/???/??  Graham Mann    Original code
!
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: AVC, RHOSUL, MMSUL, PPI
      IMPLICIT NONE

      INTEGER(KIND=JPIM) :: JMODE
      INTEGER(KIND=JPIM) :: JCP

      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_DUONLY_2MODE',0,ZHOOK_HANDLE)

! Mode names
      MODE_NAMES=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
       &                     'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      MODE_CHOICE=(/0,0,0,0,0,1,1/)
! Specify which modes are soluble
      MODESOL=(/1,1,1,1,0,0,0/)
! Component names
      COMPONENT_NAMES=                                                  &
     & (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      COMPONENT_CHOICE=(/0,0,0,0,1,0/) ! ***all cpts on except dust***
! Components that are soluble
      SOLUBLE_CHOICE=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      COMPONENT_MODE(1,1:NCP)=(/1,0,0,0,0,1/)       !allowed in nuc_sol
      COMPONENT_MODE(2,1:NCP)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      COMPONENT_MODE(3,1:NCP)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      COMPONENT_MODE(4,1:NCP)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      COMPONENT_MODE(5,1:NCP)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      COMPONENT_MODE(6,1:NCP)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      COMPONENT_MODE(7,1:NCP)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!! set dplim34 here to be 500nm to match value found from bin-mode comparison
!!      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,1.0e-6,1.0e-8,1.0e-7,1.0e-6/)
!!      ddplim1=(/1.0e-8,1.0e-7,1.0e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)
      DDPLIM0=(/1.0E-9_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-8_JPRB,1.0E-7_JPRB,1.0E-6_JPRB/)
      DDPLIM1=(/1.0E-8_JPRB,1.0E-7_JPRB,0.5E-6_JPRB,1.0E-5_JPRB,1.0E-7_JPRB,1.0E-6_JPRB,1.0E-5_JPRB/)

! Specify fixed geometric standard deviation for each mode
      SIGMAG=(/1.59_JPRB,1.59_JPRB,1.40_JPRB,2.0_JPRB,1.59_JPRB,1.59_JPRB,2.0_JPRB/) ! to match M7 but sigacc=1.4
!!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO JMODE=1,NMODES
        X(JMODE)=EXP(4.5_JPRB*LOG(SIGMAG(JMODE))*LOG(SIGMAG(JMODE)))
      ENDDO
!
! Specify threshold for ND (per cc) below which don't do calculations
!!      NUM_EPS=(/1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB,1.0e-5_JPRB/)
      !NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-4,1.0e-3,1.0e-3,1.0e-4/)
      NUM_EPS=(/1.0E-1,1.0E-2,1.0E-3,1.0E-4,1.0E-2,1.0E-3,1.0E-4/)
!!      NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3/)
!
      DO JMODE=1,NMODES
       DDPMID(JMODE)=                                                   &
        & EXP(0.5_JPRB*(LOG(DDPLIM0(JMODE))+LOG(DDPLIM1(JMODE))))
       MMID(JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPMID(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MLO (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM0(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
       MHI (JMODE)=                                                     &
        & (PPI/6.0_JPRB)*(DDPLIM1(JMODE)**3.0_JPRB)*(RHOSUL*AVC/MMSUL)*X(JMODE)
      ENDDO

! Initial fractions of mass in each mode among components
      MFRAC_0(1,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !nucln. soluble
      MFRAC_0(2,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken soluble
      MFRAC_0(3,1:NCP)=(/1.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !accum. soluble
      MFRAC_0(4,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/) !coarse soluble
      MFRAC_0(5,1:NCP)=(/0.0_JPRB,0.5_JPRB,0.5_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB/) !Aitken insoluble
      MFRAC_0(6,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !accum. insoluble
      MFRAC_0(7,1:NCP)=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB/) !coarse insoluble

! Modes resulting when two modes coagulate
      COAG_MODE(1,1:NMODES)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      COAG_MODE(2,1:NMODES)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(3,1:NMODES)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(4,1:NMODES)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(5,1:NMODES)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(6,1:NMODES)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      COAG_MODE(7,1:NMODES)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      MM=(/0.098_JPRB,0.012_JPRB,0.0168_JPRB,0.05844_JPRB,0.100_JPRB,0.0168_JPRB/)
!          h2so4  bc     oc    nacl   dust  sec_org
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 C-H ratio)
! assume secondary organic species has mm=mm_oc
! Mass density of components (kg m^-3)
      RHOCOMP=(/1769.0_JPRB,1500.0_JPRB,1500.0_JPRB,1600.0_JPRB,2650.0_JPRB,1500.0_JPRB/)
! number of dissociating ions in soluble components
      NO_IONS=(/3.0_JPRB,0.0_JPRB,0.0_JPRB,2.0_JPRB,0.0_JPRB,0.0_JPRB/)
!
      FRACBCEM=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/)
      FRACOCEM=(/0.0_JPRB,0.0_JPRB,0.0_JPRB,0.0_JPRB,1.0_JPRB,0.0_JPRB,0.0_JPRB/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
!
! Set logical variables
      MODE=(MODE_CHOICE > 0)
      COMPONENT=.FALSE.
      SOLUBLE=.FALSE.
      DO JMODE=1,NMODES
        DO JCP=1,NCP
          IF(((COMPONENT_MODE(JMODE,JCP) == 1).AND.                     &
           &    (COMPONENT_CHOICE(JCP) == 1)).AND.                      &
           &    (MODE_CHOICE(JMODE) == 1)) THEN
             COMPONENT(JMODE,JCP)=.TRUE.
          ENDIF
          IF(SOLUBLE_CHOICE(JCP) == 1) SOLUBLE(JCP)=.TRUE.
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('UKCA_MODE_SETUP:UKCA_MODE_DUONLY_2MODE',1,ZHOOK_HANDLE)

      END SUBROUTINE UKCA_MODE_DUONLY_2MODE

      END MODULE UKCA_MODE_SETUP
