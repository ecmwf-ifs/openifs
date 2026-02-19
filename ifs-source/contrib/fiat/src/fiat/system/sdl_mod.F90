! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SDL_MOD

!    Interface between user applications and system-dependent intrinsic
!    routines, provided by the computer vendors.

!    All routines which wish to call these routines must contain:
!    USE SDL_MOD

! Author :
! ------
!   11-Apr-2005 R. El Khatib  *METEO-FRANCE*
!   26-Apr-2006 S.T.Saarinen  Dr.Hook trace, calls to EC_RAISE, Intel/ifort traceback

USE EC_PARKIND  ,ONLY : JPIM
USE OML_MOD   ,ONLY : OML_MY_THREAD
USE MPL_MPIF
IMPLICIT NONE

SAVE

PRIVATE

INTEGER, PARAMETER :: SIGABRT = 6 ! Hardcoded

PUBLIC :: SDL_SRLABORT, SDL_DISABORT, SDL_TRACEBACK

CONTAINS

!-----------------------------------------------------------------------------
SUBROUTINE SDL_TRACEBACK(KTID)
USE YOMHOOK, ONLY:  DR_HOOK_CALLTREE

! Purpose :
! -------
!   Traceback

!   KTID : thread 

INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KTID

INTEGER(KIND=JPIM) ITID
CHARACTER(LEN=80) :: CLTRBK
INTEGER(KIND=JPIM) :: IERROR,IPROC
LOGICAL :: LMPI_INITIALIZED

IPROC=1
CALL MPI_INITIALIZED(LMPI_INITIALIZED,IERROR) ! always thread safe, see standard !
IF( LMPI_INITIALIZED ) THEN
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,IPROC,IERROR) ! always thread safe, see standard !
  IPROC = IPROC+1 ! 1-based in IFS context
ENDIF

IF (PRESENT(KTID)) THEN
  ITID = KTID
ELSE
  ITID = OML_MY_THREAD()
ENDIF

WRITE(0,'(A,I0,A,I0,A)') 'SDL_TRACEBACK [PROC=',IPROC,',THRD=',ITID,'] ...'
CALL DR_HOOK_CALLTREE(ITID)
#if defined(__INTEL_COMPILER)
CALL INTEL_TRBK()  ! runs LINUX_TRBK as well inside with environment EC_LINUX_TRBK=1 -- See gentrbk.F90
#else
CALL LINUX_TRBK()
CALL GDB_TRBK()    ! needs environment GNUDEBUGGER=1 -- See linuxtrbk.c
CALL DBX_TRBK()    ! needs environment DBXDEBUGGER=1 -- See linuxtrbk.c
#endif
WRITE(0,'(A,I0,A,I0,A)') 'SDL_TRACEBACK [PROC=',IPROC,',THRD=',ITID,'] ... DONE'

END SUBROUTINE SDL_TRACEBACK

!-----------------------------------------------------------------------------
SUBROUTINE SDL_SRLABORT

! Purpose :
! -------
!   To abort in serial environment

CALL EC_RAISE(SIGABRT)
STOP 'SDL_SRLABORT'

END SUBROUTINE SDL_SRLABORT
!-----------------------------------------------------------------------------
SUBROUTINE SDL_DISABORT()

! Purpose :
! -------
!   To abort in distributed environment

USE YOMHOOK, ONLY : LHOOK

INTEGER(KIND=JPIM) :: IRETURN_CODE,IERROR
CHARACTER(LEN=80) :: CLJOBID
CHARACTER(LEN=80) :: CLTRBK

#if defined(__INTEL_COMPILER)
! Intel compiler seems to hang in MPI_ABORT -- on all but the failing task(s)
! ... when linux trbk is used. REK
IF (LHOOK) THEN
  CALL GET_ENVIRONMENT_VARIABLE("EC_LINUX_TRBK",CLTRBK)
  IF (CLTRBK=='1') THEN
    CALL GET_ENVIRONMENT_VARIABLE("SLURM_JOBID",CLJOBID)
    IF (CLJOBID /= ' ') THEN
      CALL SYSTEM("set -x; sleep 10; scancel --signal=TERM "//trim(CLJOBID)//" &")
    ENDIF
  ENDIF
ENDIF
#endif

IRETURN_CODE=SIGABRT
CALL MPI_ABORT(MPI_COMM_WORLD,IRETURN_CODE,IERROR) ! Tracked by the supervisor/process-damager (manager) -- KCOMM /= MPI_COMM_WORLD may hang as sub-communicator

CALL EC_RAISE(SIGABRT) ! In case ever ends up here
STOP 'SDL_DISABORT'

END SUBROUTINE SDL_DISABORT
!-----------------------------------------------------------------------------

END MODULE SDL_MOD
