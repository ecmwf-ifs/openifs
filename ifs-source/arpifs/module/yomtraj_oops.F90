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

MODULE YOMTRAJ_OOPS

!     Purpose.
!     --------
!       Variables to control trajectory for OOPS

!     Author.
!     -------
!        O. Marsden *ECMWF*
!        Copied from original YOMTRAJ 
!
!     Modifications.
!     --------------
!        Original : September 2016
!     ------------------------------------------------------------------

USE PARKIND1           , ONLY : JPIM, JPRB, JPRM
USE SPECTRAL_FIELDS_MOD, ONLY : ASSIGNMENT(=), SPECTRAL_FIELD
USE YOMTRAJ
USE SURFACE_FIELDS_MIX , ONLY : TSURF

IMPLICIT NONE
SAVE
PUBLIC


LOGICAL :: LTRAJALLOC_oops= .FALSE.  ! Trajectory storage is allocated
LOGICAL :: LTRAJSAVE_oops = .FALSE.  ! Saving trajectory
LOGICAL :: LTRAJCST_oops  = .FALSE.  ! Constants in trajectory (needs tidying-up)
LOGICAL :: LTRAJSLAG_oops = .FALSE.  ! Semi-Lagrangian trajectory is used
LOGICAL :: LTRAJPHYS_oops = .FALSE.  ! Physics trajectory is used
LOGICAL :: LTRAJRESET_oops= .FALSE.  ! Recomputing and changing trajectory


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
INTEGER(KIND=JPIM) :: NSPTRAJ_oops
INTEGER(KIND=JPIM) :: NGPTRAJ_oops
INTEGER(KIND=JPIM) :: NTRAJP_oops    !  should be removed
! INTEGER(KIND=JPIM) :: NGP5_oops
! INTEGER(KIND=JPIM) :: NTRAJ_CST_oops
INTEGER(KIND=JPIM), ALLOCATABLE :: MIOTRAJMAIN_oops(:)  !! olivier addition
INTEGER(KIND=JPIM), ALLOCATABLE :: MSTEPTRAJ_oops(:)

! For no packing, replace previous line by
! INTEGER, PARAMETER :: MKINDTRAJ = JPRB
! For packing, replace previous line by
! INTEGER, PARAMETER :: MKINDTRAJ = JPRM

! Constant trajectory (should not exist...)
REAL(KIND=MKINDTRAJ), ALLOCATABLE :: CONSTANT_TRAJ_GLOBAL(:,:,:)


! * TRAJ_GRIB, MAIN_GRIB, BACKGR_GRIB: information for grib header
! * MTYPE_[X]_TRAJ: grib codes ([X] = SURF for surface, SLAG for semi-Lag data,
!   PHYS for physics data, MAIN3 for upper air 3D variables, MAIN2 for
!   2D variables).
INTEGER(KIND=JPIM) :: MTRAJ_GRIB_oops
INTEGER(KIND=JPIM) :: MAIN_GRIB_oops
INTEGER(KIND=JPIM), ALLOCATABLE :: MBACKGR_GRIB_oops(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: MTYPE_SURF_TRAJ_oops(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: MTYPE_MAIN3_TRAJ_oops(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: MTYPE_MAIN2_TRAJ_oops(:)



!-----------------------------------------------------------------------



! Trajectory types


TYPE TRAJ_MAIN_TYPE_oops
  REAL(KIND=MKINDTRAJ), POINTER :: GFL(:,:,:,:) => NULL()   !! (1:NPROMA, 1:NFLEVG, 1:YGFL%NDIM5, 1:TDIM%NGPBLKS)
  REAL(KIND=MKINDTRAJ), POINTER :: GMV(:,:,:,:) => NULL()   !! (1:NPROMA, 1:NFLEVG, 1:TGMV5%YT5%NDIM, 1:TDIM%NGPBLKS)
  REAL(KIND=MKINDTRAJ), POINTER :: GMVS(:,:,:)  => NULL()   !! (1:NPROMA, 1:TGMV5%YT5%NDIMS, 1:TDIM%NGPBLKS))
  TYPE(SPECTRAL_FIELD) :: SPEC
END TYPE TRAJ_MAIN_TYPE_oops

TYPE TRAJ_SRFC_TYPE_oops
  REAL(KIND=MKINDTRAJ), ALLOCATABLE :: TRAJSURF(:,:,:)
  TYPE(TSURF)                       :: Y_SURF ! Only metadata
END TYPE TRAJ_SRFC_TYPE_oops


TYPE TRAJ_TYPE_oops
  TYPE(TRAJ_MAIN_TYPE_oops)    :: MAIN
  TYPE(TRAJ_SRFC_TYPE_oops)    :: SURFACE
  TYPE(TRAJ_SLAG_TYPE),     ALLOCATABLE :: SLAG(:) 
  TYPE(TRAJ_PHYS_TYPE),     ALLOCATABLE :: PHYS(:) 
  TYPE(TRAJ_PHYS_TLAD_TYPE),ALLOCATABLE :: PHYS_TLAD(:)  
  TYPE(TRAJ_SRFC_TYPE)     ,ALLOCATABLE :: SRFC(:)
  TYPE(TRAJ_CST_TYPE)      ,ALLOCATABLE :: CST(:)
END TYPE TRAJ_TYPE_oops



!-----------------------------------------------------------------------

END MODULE YOMTRAJ_OOPS
