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

MODULE YOMTRAJ

!     Purpose.
!     --------
!       Variables to control trajectory

!     Author.
!     -------
!        Y. Tremolet *ECMWF*

!     Modifications.
!     --------------
!        Original : 2002-12-03
!        Y.Tremolet: 14-01-04 Rename a few things and add comments
!        Y.Tremolet    18-Mar-2004 Add lprttraj option for diagnostics
!        K.Yessad  : 27-Jun-2007 Add missing comments
!        F. Vana  28-Nov-2013 : Redesigned trajectory handling
!        M. Janiskova  03-May-2012: fields for surface scheme
!        M/ Fisher 9 July 2015: Moved BCKG into YOMSPJB and made TRAJEC and BCKG the same type.
!        E. Holm  22 July 2016: Hacked BACKGR04-10, need to do neater!
!        F. Vana  21 Nov 2017 : Option LSLDP_CURV and fixed type of IGFLL95
!        F. Vana   July 2018: RK4 scheme for trajectory research.
!        F. Vana  11-Jul-2019: Option LRHS_CURV
!        F. Vana  11-Sep-2020: SLAVEPP and prognostic cloud
!        F. Vana  22-Sep-2021: Specific CAMS trajectory for LETRAJPT
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRM
USE SPECTRAL_FIELDS_MOD, ONLY : ASSIGNMENT(=), SPECTRAL_FIELD

IMPLICIT NONE
SAVE
PUBLIC

LOGICAL :: LTRAJHR              ! Use interpolated trajectory
LOGICAL :: LTRAJHR_ALTI         ! Use interpolated trajectory (upper air flds)
LOGICAL :: LTRAJHR_SURF         ! Use interpolated trajectory (surface flds)
LOGICAL :: LTRAJGP              ! Store trajectory in grid-point space
LOGICAL :: LTRAJSAVE = .FALSE.  ! Saving trajectory
LOGICAL :: LTRAJSLAG = .FALSE.  ! Semi-Lagrangian trajectory is used
LOGICAL :: LTRAJPHYS = .FALSE.  ! Physics trajectory is used
LOGICAL :: LTRAJCST  = .FALSE.  ! Constants in trajectory (needs tidying-up)
LOGICAL :: LTRAJALLOC= .FALSE.  ! Trajectory storage is allocated
LOGICAL :: LPRTTRAJ  = .FALSE.  ! Prints more info from trajectory routines
LOGICAL :: LTRAJRESET= .FALSE.  ! Recomputing and changing trajectory
LOGICAL :: LREADGPTRAJ          ! Write/read main trajectory in grid-point space
LOGICAL :: LSNOWTRAJCONS       ! Use conservative interpolation for snow multilayer

REAL(KIND=JPRB) :: TSTEP_TRAJ       ! Time interval between saved trajectory steps
INTEGER(KIND=JPIM) :: NSMAX_TRAJ    ! Spectral truncation for trajectory


! * NSMAX_BACKGR00 to NSMAX_BACKGR10: Spectral truncation for background state
!   for the various inner loops.
INTEGER(KIND=JPIM) :: NSMAX_BACKGR00, NSMAX_BACKGR01, NSMAX_BACKGR02, &
 & NSMAX_BACKGR03, NSMAX_BACKGR04, NSMAX_BACKGR05, NSMAX_BACKGR06, &
 & NSMAX_BACKGR07, NSMAX_BACKGR08, NSMAX_BACKGR09, NSMAX_BACKGR10

! * NSPTRAJ          : number of spectral fields stored in the trajectory
! * NGPTRAJ          : number of grid-point fields stored in the trajectory
! * NTRAJP           : number of fields in buffer TRAJ_PHYS
! * NGP5             : number of fields in TRAJEC%SRFC
! * NTRAJ_CST        : number of fields in TRAJEC%CST
! * NSTEPTRAJ        : number of timesteps where the trajectory must be stored
! * MSTART           : MSTEPTRAJW(NSTART) == 1 (by definition)
! * MSTEPTRAJW       : numbering for timesteps where trajectory is written
! * MSTEPTRAJR       : numbering for timesteps where trajectory is read
! * MIOTRAJMAIN      : ???
INTEGER(KIND=JPIM) :: NSPTRAJ
INTEGER(KIND=JPIM) :: NGPTRAJ
INTEGER(KIND=JPIM) :: NTRAJP    !  should be removed
INTEGER(KIND=JPIM) :: NGP5
INTEGER(KIND=JPIM) :: NTRAJ_CST
INTEGER(KIND=JPIM) :: NSTEPTRAJ
INTEGER(KIND=JPIM) :: MSTART
INTEGER(KIND=JPIM), ALLOCATABLE :: MSTEPTRAJW(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: MSTEPTRAJR(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: MIOTRAJMAIN(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: MIOTRAJSURF(:)
INTEGER(KIND=JPIM), PARAMETER :: MKINDTRAJ = JPRM

! For no packing, replace previous line by
! INTEGER, PARAMETER :: MKINDTRAJ = JPRB
! For packing, replace previous line by
! INTEGER, PARAMETER :: MKINDTRAJ = JPRM

! Trajectory types

! Trajectory used in physics from NL to TL/AD models
TYPE TRAJ_PHYS_TYPE
  LOGICAL :: LASTCHUNK
  ! Types used by ECMWF (defined as single precission)
  REAL(KIND=MKINDTRAJ), DIMENSION(:,:), POINTER :: PTT95=>NULL(), PQT95=>NULL(), PUT95=>NULL(), PVT95=>NULL(), &
  &    PTET5=>NULL(), PTEQ5=>NULL(), PTEU5=>NULL(), PTEV5=>NULL(), PVVEL5=>NULL(), PGFLT9_LRCH45=>NULL()   ! NFLEV
  REAL(KIND=JPRB), DIMENSION(:,:), POINTER ::  PLT95=>NULL(), PIT95=>NULL(), PTEL5=>NULL(), PTEI5=>NULL() ! NFLEV cloud variables
  REAL(KIND=MKINDTRAJ), DIMENSION(:,:), POINTER :: PEMTED5=>NULL(), PTRSOL5=>NULL()  ! NFLEV+1
  REAL(KIND=MKINDTRAJ), DIMENSION(:,:,:), POINTER :: PGFLT9_GHG5=>NULL(), PTEGFL_GHG5=>NULL(), &
                      & PGFLT9_AERO5=>NULL(), PTEGFL_AERO5=>NULL(), PGFLT9_CHEM5=>NULL(), PTEGFL_CHEM5=>NULL()
  REAL(KIND=MKINDTRAJ), DIMENSION(:,:), POINTER :: PTSA95=>NULL(), PWSA95=>NULL(), PSNS95=>NULL(),PRSN95=>NULL(), PTIA95=>NULL()
  REAL(KIND=MKINDTRAJ), DIMENSION(:), POINTER   :: PSPT95=>NULL(), PSPT15=>NULL(), PWL95=>NULL(), &
  &    PTL95=>NULL(), PNTOP5=>NULL(), PNBAS5=>NULL(), PACPR5=>NULL(), PACCPR5=>NULL(), PGZ0M5=>NULL(), PLZ0H5=>NULL()
  REAL(KIND=MKINDTRAJ), DIMENSION(:,:), POINTER :: PUSTRTI5=>NULL(), PVSTRTI5=>NULL(), PAHFSTI5=>NULL(), &
  &    PEVAPTI5=>NULL(), PTSKTI5=>NULL()  ! NTILES
  REAL(KIND=JPRB)     , DIMENSION(:,:), POINTER      :: RR2=>NULL()   ! rainggauges

  ! Types used by Meteo-France physics (double precission at the moment)
  !  In order to be distinguished from ECMWF arrays in single precission, they are appendixed by 'MF')
  REAL(KIND=JPRB), DIMENSION(:,:), POINTER  :: PUT9MF5=>NULL(),PVT9MF5=>NULL(),PTT9MF5=>NULL(),PQT9MF5=>NULL(), &
  &    PQLT9MF5=>NULL(),PQIT9MF5=>NULL()
  REAL(KIND=JPRB), DIMENSION(:)  , POINTER  :: PSP9MF5=>NULL()
  REAL(KIND=JPRB), DIMENSION(:,:), POINTER  :: PKTROVMF=>NULL(),PKUROVMF=>NULL(),PRAPGWDMF=>NULL() ! NFLEV+1
  REAL(KIND=JPRB), DIMENSION(:,:), POINTER  :: PQLMF=>NULL(),PQIMF=>NULL(),PQRMF=>NULL(),PQSMF=>NULL(),PLIQCVPPKFMF=>NULL(), &
  &    PNEBCVPPKFMF=>NULL()
  REAL(KIND=JPRB), DIMENSION(:)  , POINTER  :: PCDROVMF=>NULL(),PCHROVMF=>NULL()
  REAL(KIND=JPRB), DIMENSION(:)  , POINTER  :: PQSSMF=>NULL(),PTSMF=>NULL(),PSNSMF=>NULL()  ! surface quantities trajs
  REAL(KIND=JPRB), DIMENSION(:)  , POINTER  :: PQSS1MF=>NULL(),PTS1MF=>NULL()  ! surface quantities used in next time-step
END TYPE TRAJ_PHYS_TYPE

! Trajectory to handle semi-prognostic cloud (between SL and physics)
TYPE CLD_TRAJ_TYPE
  REAL(KIND=JPRB), DIMENSION(:,:), POINTER  :: YA=>NULL(),YL=>NULL(),YI=>NULL()
END TYPE CLD_TRAJ_TYPE

! Trajectory used in physics to communicate between TL and AD models without time dimension
TYPE TRAJ_RAINGG_TLAD_TYPE
  LOGICAL :: LASTCHUNK
  REAL(KIND=JPRB), DIMENSION(:,:) , POINTER :: RR=>NULL()
END TYPE TRAJ_RAINGG_TLAD_TYPE

! Trajectory used in TL and AD physics to save computation in the AD model. 
!   (It must remain the same precision for both models, i.e. KIND parameter
!    has to be the same as in TL model (being currently JPRB)!)
TYPE TRAJ_PHYS_TLAD_TYPE
  LOGICAL :: LASTCHUNK
  ! ECMWF fields
  REAL (KIND=JPRB), DIMENSION (:,:), POINTER :: PEMTED5=>NULL(), PTRSOL5=>NULL()   ! NFLEVG+1
  REAL (KIND=JPRB), DIMENSION (:,:), POINTER :: PA5=>NULL(), PL5=>NULL(), PI5=>NULL(), &
   & PLUDE5=>NULL(), PLU5=>NULL(), PMFU5=>NULL(), PMFD5=>NULL()
  REAL (KIND=JPRB), DIMENSION (:,:), POINTER :: PTENUWM5=>NULL(), PTENVWM5=>NULL()
  REAL (KIND=JPRB), DIMENSION (:,:), POINTER :: PTENTVD5=>NULL(), PTENQVD5=>NULL(), PTENUVD5=>NULL(), &
   & PTENVVD5=>NULL(), PTENTGW5=>NULL(), PTENUGW5=>NULL(), PTENVGW5=>NULL(), PTENTCO5=>NULL(), PTENQCO5=>NULL(), &
   & PTENUCO5=>NULL(), PTENVCO5=>NULL(), PTENTSL5=>NULL(), PTENQSL5=>NULL()
  REAL (KIND=JPRB), DIMENSION (:,:), POINTER :: PDIFTQ5=>NULL(), PDIFTS5=>NULL()   ! NFLEVG+1
  REAL (KIND=JPRB), DIMENSION (:,:), POINTER :: PTENTCL5=>NULL()
  REAL (KIND=JPRB), DIMENSION (:,:), POINTER :: PFPLCL5=>NULL(), PFPLCN5=>NULL()   ! NFLEVG+1
  REAL (KIND=JPRB), DIMENSION (:), POINTER :: PUCFL5=>NULL(),PVCFL5=>NULL(),PTCFL5=>NULL(),PQCFL5=>NULL()
  !  surface scheme
  REAL (KIND=JPRB), DIMENSION (:,:), POINTER :: PFRSOTI5=>NULL(),PAHFSTI5=>NULL(),PEVAPTI5=>NULL()
  REAL (KIND=JPRB), DIMENSION (:), POINTER :: PEVAPSNW5=>NULL(),PFPLSL5=>NULL(),PFPLSN5=>NULL(),PFRTH5=>NULL()
  ! CAMS tendencies for LETRAJP
  REAL (KIND=JPRB), DIMENSION (:,:,:), POINTER :: TENDC55D=>NULL()
  ! LSLPHY (including SL trajectory)
  REAL (KIND=JPRB),DIMENSION(:,:), POINTER :: PHTENU5=>NULL(),PHTENV5=>NULL(),PHTENT5=>NULL(),PHTENQ5=>NULL()
  REAL (KIND=JPRB),DIMENSION(:,:), POINTER :: PUP95=>NULL(),PVP95=>NULL()
  REAL (KIND=JPRB),DIMENSION(:,:), POINTER :: PUSLP95=>NULL(),PVSLP95=>NULL(),PZSLP95=>NULL(), &
   & PTSLP95=>NULL(),PQSLP95=>NULL()
  ! MF fields
  REAL(KIND=JPRB), DIMENSION(:,:), POINTER  :: PCVGQMF=>NULL(),PTENDUMF=>NULL(),PTENDVMF=>NULL(),PTENDHMF=>NULL(),&
   & PTENDQMF=>NULL(), PTENDQLMF=>NULL(), PTENDQIMF=>NULL()
  REAL(KIND=JPRB), DIMENSION(:,:), POINTER  :: PQLMF=>NULL(),PQIMF=>NULL(),PQRMF=>NULL(),PQSMF=>NULL(),PLIQCVPPKFMF=>NULL(),&
   & PNEBCVPPKFMF=>NULL()
  REAL(KIND=JPRB), DIMENSION(:,:), POINTER  :: PDIFCQMF=>NULL(),PDIFCSMF=>NULL(),PDIFTQMF=>NULL(),PDIFTSMF=>NULL(), &
   & PFCCQLMF=>NULL(),PFCCQNMF=>NULL(),PFCSQLMF=>NULL(),PFCSQNMF=>NULL(),PFPLCLMF=>NULL(),PFPLCNMF=>NULL(),PFPLSLMF=>NULL(),&
   & PFPLSNMF=>NULL(),PFRSOMF=>NULL(),PFRTHMF=>NULL(), PSTRCUMF=>NULL(),PSTRCVMF=>NULL(),PSTRDUMF=>NULL(),PSTRDVMF=>NULL(),&
   & PSTRTUMF=>NULL(),PSTRTVMF=>NULL(),PSTRMUMF=>NULL(),PSTRMVMF=>NULL(),PFRMHMF=>NULL() ! NFLEV+1
  REAL(KIND=JPRB), DIMENSION(:,:), POINTER  :: PKTROVMF=>NULL(),PKUROVMF=>NULL(),PXTROVMF=>NULL(),&
   & PXUROVMF=>NULL(),PRAPGWDMF=>NULL() ! NFLEV+1
  REAL(KIND=JPRB), DIMENSION(:)  , POINTER  :: PCDROVMF=>NULL(),PCHROVMF=>NULL(),PXDROVMF=>NULL(),PXHROVMF=>NULL(),&
   & PTRIGCONVMF=>NULL()
END TYPE TRAJ_PHYS_TLAD_TYPE

! Just preserving the current way using the surface_fields_mix tools.
TYPE TRAJ_CST_TYPE
  LOGICAL :: LASTCHUNK
  REAL(KIND=MKINDTRAJ), DIMENSION(:,:), POINTER :: MIX_CST=>NULL()
END TYPE TRAJ_CST_TYPE

TYPE TRAJ_SRFC_TYPE
  LOGICAL :: LASTCHUNK
  REAL (KIND=MKINDTRAJ), DIMENSION(:,:), POINTER :: MIX_SRFC=>NULL()
END TYPE TRAJ_SRFC_TYPE

TYPE TRAJ_SLAG_TYPE
  LOGICAL :: LASTCHUNK
  REAL(KIND=MKINDTRAJ), DIMENSION(:,:), POINTER :: PUL95=>NULL(), PVL95=>NULL(), PTL95=>NULL(), &
   &  PCL95=>NULL(), PUL05=>NULL(), PVL05=>NULL(), PTL05=>NULL(), PUT15=>NULL(), PVT15=>NULL(), &
   &  PURL05=>NULL(), PVRL05=>NULL(), PZRL05=>NULL(), PWRL05=>NULL(), PZL95=>NULL(), PZL05=>NULL(), &
   &  PURL005=>NULL(), PVRL005=>NULL(), PZRL005=>NULL(), PWRL005=>NULL(), &
   &  PURL5=>NULL(), PVRL5=>NULL(), PWRL5=>NULL(), PZRL5=>NULL(), &
   &  PKAPPA5=>NULL(), PKAPPAT5=>NULL(), PTT15=>NULL(), PQT15=>NULL(), PWRL95=>NULL(), &
   &  PSTDDISU5=>NULL(), PSTDDISV5=>NULL(), PSTDDISW5=>NULL()
  REAL(KIND=MKINDTRAJ), DIMENSION(:,:,:), POINTER :: PSAVEDP5=>NULL()
  REAL(KIND=JPRB), DIMENSION(:,:,:), POINTER :: PGFLL95=>NULL() ! better to be JPRB
  INTEGER(KIND=JPIM), DIMENSION(:), POINTER :: IGFLL95=>NULL() ! flag identifying trajectory content
  REAL(KIND=MKINDTRAJ), DIMENSION(:), POINTER :: PX95=>NULL()
  ! quantities required by 2D model only
  REAL(KIND=MKINDTRAJ), DIMENSION(:), POINTER :: PUT95=>NULL(), PVT95=>NULL(), PSPNLT95=>NULL()       ! NTRSLTYPE=1
  REAL(KIND=MKINDTRAJ), DIMENSION(:), POINTER :: PUSI5=>NULL(), PVSI5=>NULL(), PSPL95=>NULL(), PSPL05=>NULL() ! NTRSLTYPE=2
END TYPE TRAJ_SLAG_TYPE

TYPE TRAJ_MAIN_TYPE
  REAL(KIND=MKINDTRAJ), DIMENSION (:,:,:), POINTER  :: GFL=>NULL(), GMV=>NULL()
  REAL(KIND=MKINDTRAJ), DIMENSION (:,:)  , POINTER  :: GMVS=>NULL()
END TYPE TRAJ_MAIN_TYPE

! Super-trajectory structure
TYPE TRAJ_TYPE
  TYPE (SPECTRAL_FIELD)       , POINTER :: SPEC(:)=>NULL()
  TYPE (TRAJ_SRFC_TYPE)       , POINTER :: SRFC(:)=>NULL()
  TYPE (TRAJ_CST_TYPE)        , POINTER :: CST(:)=>NULL()
  TYPE (TRAJ_SLAG_TYPE)       , POINTER :: SLAG(:)=>NULL()
  TYPE (TRAJ_PHYS_TYPE)       , POINTER :: PHYS(:)=>NULL()
  TYPE (TRAJ_PHYS_TLAD_TYPE)  , POINTER :: PHYS_TLAD(:)=>NULL()
  TYPE (TRAJ_RAINGG_TLAD_TYPE), POINTER :: RAINGG_TLAD(:)=>NULL()
  TYPE (TRAJ_MAIN_TYPE)       , POINTER :: MAIN(:)=>NULL()
END TYPE TRAJ_TYPE

!-----------------------------------------------------------------------
! * TRAJ%<x> contains the actual trajectory
!   where <x> = SPEC for spectral data, SRFC for surface data, CST for constants,
!    SLAG for semi-Lagrangian data, PHYS for physics data, ZERO (cf. PHYS
!    but for some particuliar instants), GMV for GMV variables, GMVS for GMVS
!    variables, GFL for GFL variables.
!   TRAJ_SPEC_TMP0 is a copy of TRAJ_SPEC at the initial instant.


TYPE (SPECTRAL_FIELD)                 :: TRAJ_SPEC_TMP0

! Main model trajectory (containing all the other structures)
TYPE (TRAJ_TYPE)            , POINTER :: TRAJEC(:) =>NULL()

!-----------------------------------------------------------------------

! * TRAJ_GRIB, MAIN_GRIB, BACKGR_GRIB: information for grib header
! * MTYPE_[X]_TRAJ: grib codes ([X] = SURF for surface, SLAG for semi-Lag data,
!   PHYS for physics data, MAIN3 for upper air 3D variables, MAIN2 for
!   2D variables).
INTEGER(KIND=JPIM) :: MTRAJ_GRIB
INTEGER(KIND=JPIM) :: MTRAJ_GRIB2
INTEGER(KIND=JPIM) :: MAIN_GRIB
INTEGER(KIND=JPIM), ALLOCATABLE :: MBACKGR_GRIB(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: MTYPE_SURF_TRAJ(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: MTYPE_MAIN3_TRAJ(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: MTYPE_MAIN2_TRAJ(:)

!-----------------------------------------------------------------------

END MODULE YOMTRAJ


