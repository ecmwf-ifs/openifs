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

MODULE YOMFP_SERV

!**** *YOMFP_SERV*  - FP server definition

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO-FRANCE*
!      Original : 15-04-2016


USE PARKIND1, ONLY : JPIM
USE MPL_MPIF, ONLY : MPI_COMM_NULL
USE YOMIO_SERV, ONLY : IO_SERV
USE YOMFP_SERV_DINF, ONLY : FP_SERV_DINF

IMPLICIT NONE

INTEGER (KIND=JPIM), PARAMETER :: NFP_SERV_SYN_TYPE_CLO = 1, &
                                & NFP_SERV_SYN_TYPE_STO = 2

TYPE FP_SERV_DISTR
! Grid-point distribution
! Sort points received by procs from remote model (NGPTOTMX_R, NPROC_RG)
  INTEGER (KIND=JPIM), POINTER :: ISORTR2L (:,:)   => NULL ()
! Remote processors this task communicates with
  INTEGER (KIND=JPIM), POINTER :: IPROCR2L (:)     => NULL ()
! Number of points to exchange with remote processors
  INTEGER (KIND=JPIM), POINTER :: ISORTCNT (:)     => NULL ()
  INTEGER (KIND=JPIM), POINTER :: ISORTOFF (:)     => NULL ()
  TYPE (FP_SERV_DINF), POINTER :: YFSDINF  (:)     => NULL ()
  INTEGER (KIND=JPIM), POINTER :: IRANKSET_R (:,:) => NULL () ! (IPROC, 1) -> IPROC, (ISETW, ISETV) -> IPROC
  INTEGER (KIND=JPIM), POINTER :: ISETLEV_R (:)    => NULL () ! (ILEV) -> 1        , (ILEV) -> ISETV        
  INTEGER (KIND=JPIM), POINTER :: ISETLEV_L (:)    => NULL () ! (ILEV) -> 1        , (ILEV) -> ISETV        
  INTEGER (KIND=JPIM), POINTER :: ISETV_R (:)      => NULL () ! List of remote V-sets we need to talk with  
! Number of remote processors we need to talk to
  INTEGER (KIND=JPIM) :: NPROC_R  = 0
END TYPE FP_SERV_DISTR


TYPE FP_SERV

! Fp server has sent some data
  LOGICAL :: LSENTFLD = .FALSE.

! Fp client writes historic files

  LOGICAL :: LFP_CLIENT_WRITE = .FALSE.

! Fp client computes filtering matrices and sends them to the Fp server

  LOGICAL :: LFP_SERVER_FPMTS = .TRUE.

! Global communicator
  INTEGER (KIND=JPIM) :: MYPROC   = 0
  INTEGER (KIND=JPIM) :: NPROC    = 0
  INTEGER (KIND=JPIM) :: NCOMM    = MPI_COMM_NULL

! Local communicator
  INTEGER (KIND=JPIM) :: MYPROC_L = 0
  INTEGER (KIND=JPIM) :: NPROC_L  = 0
  INTEGER (KIND=JPIM) :: NCOMM_L  = MPI_COMM_NULL

! Remote communicator
  INTEGER (KIND=JPIM) :: NPROC_R  = 0


  LOGICAL :: LFP_SERVER = .FALSE.
  LOGICAL :: LFP_CLIENT = .FALSE.
  INTEGER (KIND=JPIM) :: NPROC_WR = 0 
  INTEGER (KIND=JPIM) :: NPROC_FP = 0

  INTEGER (KIND=JPIM), POINTER :: MYPROCS_L (:)      => NULL ()
  INTEGER (KIND=JPIM), POINTER :: MYPROCS_R (:)      => NULL ()

  TYPE (FP_SERV_DISTR) :: YGP_DISTR
  TYPE (FP_SERV_DISTR) :: YSP_DISTR

! For Fullpos server
  TYPE (IO_SERV)               :: YIOSS
! For Fullpos client
  TYPE (IO_SERV),      POINTER :: YIOS_GP_CL  (:,:)  => NULL ()
  TYPE (IO_SERV),      POINTER :: YIOS_SP_CL  (:,:)  => NULL ()

END TYPE FP_SERV

TYPE (FP_SERV) :: FP_SERV_C001

SAVE

END MODULE YOMFP_SERV

