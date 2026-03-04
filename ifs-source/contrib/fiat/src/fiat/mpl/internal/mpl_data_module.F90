! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_DATA_MODULE

!     Author.
!     -------
!        D.Dent, M.Hamrud     ECMWF

!     Modifications.
!     --------------
!        Original: 2000-09-01

!     ------------------------------------------------------------------

!      variables controlling the execution of MPL

!  MPL_METHOD   : buffering type
!  MPL_MBX_SIZE : size of application mailbox, (bytes)
!                 used when MPL_METHOD=JP_BLOCKING_BUFFERED
!  MPL_COMM     : default communicator in use
!  MPL_COMM_OML : communicators for messages between corresponding OML-threads
!  MPL_UNIT     : Fortran I/O unit for messages (default=6)
!  MPL_ERRUNIT  : Fortran I/O unit for error messages (default=0)
!  MPL_OUTPUT   : controls contents of Output (see mpl_init_mod.F90 for values/default)
!  MPL_RANK     : rank of the process within MPL_COMM_OML(1)
!  MPL_NUMPROC  : number of processes in MPL_COMM_OML(1)
!  MPL_IDS      : array of processor numbers
!  LUSEHLMPI    : always use high level MPI calls (collective comm.)
!  LINITMPI_VIA_MPL : true if MPI has been initialized from within MPL_INIT()
!  LTHSAFEMPI   : Thread safe MPI, if .TRUE. (default)

USE MPL_MPIF , ONLY : MPI_COMM_WORLD
USE EC_PARKIND  ,ONLY : JPIM

IMPLICIT NONE

PRIVATE :: JPIM

SAVE

PUBLIC 

INTEGER(KIND=JPIM) :: MPL_METHOD, MPL_MBX_SIZE, MPL_UNIT=6, MPL_OUTPUT=1
INTEGER(KIND=JPIM) :: MPL_RANK=0,MPL_NUMPROC = -1,MPL_ERRUNIT=0
INTEGER(KIND=JPIM),ALLOCATABLE :: MPL_IDS(:)
INTEGER(KIND=JPIM) :: MPL_COMM
INTEGER(KIND=JPIM),ALLOCATABLE :: MPL_COMM_OML(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: MPL_OPPONENT(:)
INTEGER(KIND=JPIM) :: MPL_NCPU_PER_NODE=1
INTEGER(KIND=JPIM) :: MPL_MAX_TASK_PER_NODE
INTEGER(KIND=JPIM),ALLOCATABLE :: MPL_TASK_PER_NODE(:)
INTEGER(KIND=JPIM) :: MPL_NNODES
LOGICAL :: LFULLNODES
INTEGER(KIND=JPIM) :: MPL_MYNODE=0
INTEGER(KIND=JPIM),ALLOCATABLE :: MPL_NODE(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: MPL_NODE_TASKS(:,:)
!INTEGER_M,ALLOCATABLE :: MPL_ATTACHED_BUFFER(:)
!   needs to ge a TARGET for coexistence with MPE
INTEGER(KIND=JPIM),ALLOCATABLE,TARGET :: MPL_ATTACHED_BUFFER(:)
LOGICAL :: LUSEHLMPI
LOGICAL :: LINITMPI_VIA_MPL = .FALSE.
LOGICAL :: LTHSAFEMPI = .TRUE.
INTEGER(KIND=JPIM),PARAMETER :: JP_ATTACHED_BUFFER_BYTES    = JPIM
INTEGER(KIND=JPIM),PARAMETER :: JP_BLOCKING_STANDARD        = 1
INTEGER(KIND=JPIM),PARAMETER :: JP_BLOCKING_BUFFERED        = 2
INTEGER(KIND=JPIM),PARAMETER :: JP_BLOCKING_SYNCHRONOUS     = 3
INTEGER(KIND=JPIM),PARAMETER :: JP_BLOCKING_READY           = 4
INTEGER(KIND=JPIM),PARAMETER :: JP_NON_BLOCKING_STANDARD    = 5
INTEGER(KIND=JPIM),PARAMETER :: JP_NON_BLOCKING_BUFFERED    = 6
INTEGER(KIND=JPIM),PARAMETER :: JP_NON_BLOCKING_SYNCHRONOUS = 7
INTEGER(KIND=JPIM),PARAMETER :: JP_NON_BLOCKING_READY       = 8
LOGICAL :: LMPLUSERCOMM = .FALSE.
INTEGER(KIND=JPIM) :: MPLUSERCOMM = -1
INTEGER(KIND=JPIM) :: MPL_SEND_COUNT, MPL_SEND_BYTES
INTEGER(KIND=JPIM) :: MPL_RECV_COUNT, MPL_RECV_BYTES
INTEGER(KIND=JPIM) :: MPL_WORLD_RANK = -1
INTEGER(KIND=JPIM) :: MPL_WORLD_SIZE =  0

END MODULE MPL_DATA_MODULE
