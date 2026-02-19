! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_ARG_MOD

!**** MPL_GETARG : A substitute for GET_COMMAND_ARGUMENT (formerly GETARG) for MPL applications
!     MPL_IARGC  : A substitute for function COMMAND_ARGUMENT_COUNT (formerly IARGC) for MPL applications

!     Purpose.
!     --------
!     MPL-task#1 calls GET_COMMAND_ARGUMENT until COMMAND_ARGUMENT_COUNT() arguments read
!     or until the argument is a terminating argument
!     Then arguments are passed on to other processors
!     If MPL has not been initialized, it will be done now.

!**   Interface.
!     ----------
!        CALL MPL_GETARG(KARG, CDARG)

!        Input required arguments :
!        -------------------------
!           KARG     -  The argument number requested (INTEGER(4))
!                       Range : [ 0 .. MPL_IARGC() ]

!        Output required arguments :
!        ---------------------------
!           CDARG    -  Return argument value (CHARACTER(LEN=*))
!
!**   Interface.
!     ----------
!        INUM_ARGS = MPL_IARGC()
!
!     where INUM_ARGS is INTEGER(4)

!     Author.
!     -------
!        S.Saarinen, G.Mozdzynski  ECMWF

!     Modifications.
!     --------------
!        Original: 2006-03-15

USE EC_PARKIND  ,ONLY : JPIM
USE MPL_MPIF
USE MPL_DATA_MODULE, ONLY : MPL_NUMPROC,LINITMPI_VIA_MPL,LMPLUSERCOMM,MPLUSERCOMM

IMPLICIT NONE

PRIVATE

CHARACTER(LEN=10), SAVE :: CL_TERMINATE = '-^' ! terminating argument

INTEGER(KIND=JPIM), PARAMETER :: JP_ARGLEN = 1024
CHARACTER(LEN=JP_ARGLEN), ALLOCATABLE, SAVE :: CL_ARGS(:)
INTEGER(KIND=JPIM), SAVE :: N_ARGS = -1

PUBLIC :: MPL_GETARG
PUBLIC :: MPL_IARGC
PUBLIC :: MPL_ARG_SET_CL_TERMINATE
PUBLIC :: MPL_ARG_GET_CL_TERMINATE

CONTAINS

SUBROUTINE MPL_ARG_SET_CL_TERMINATE(CDTERM)
CHARACTER(LEN=*), INTENT(IN) :: CDTERM
CL_TERMINATE = CDTERM
END SUBROUTINE MPL_ARG_SET_CL_TERMINATE

SUBROUTINE MPL_ARG_GET_CL_TERMINATE(CDTERM)
CHARACTER(LEN=*), INTENT(OUT) :: CDTERM
CDTERM = CL_TERMINATE
END SUBROUTINE MPL_ARG_GET_CL_TERMINATE

SUBROUTINE INIT_ARGS()

#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_INITIALIZED => MPI_INITIALIZED8, MPI_COMM_SIZE => MPI_COMM_SIZE8, &
    MPI_COMM_RANK => MPI_COMM_RANK8, MPI_BCAST => MPI_BCAST8, &
    MPI_INIT => MPI_INIT8
#endif

INTEGER(KIND=JPIM) :: IARGS
INTEGER(KIND=JPIM) :: IERROR, IROOT, ICOUNT
INTEGER(KIND=JPIM) :: IRANK, INUMPROC, IRET, J
INTEGER(KIND=JPIM) :: IARGC_C
CHARACTER(LEN=LEN(CL_TERMINATE)) :: ENV_CL_TERMINATE
CHARACTER(LEN=JP_ARGLEN) :: CLARG0
LOGICAL LLCARGS
INTEGER LLINIT
INTEGER(KIND=JPIM) :: ICOMM

IF (N_ARGS == -1) THEN
  IF (MPL_NUMPROC == -1) THEN
    ! This is complicated, but I hope it works:
    ! MPI has not yet been initialized, when this routines was called.
    ! Initialize MPI, but NOT via MPL_INIT to avoid recursion in MPL_IARGC()
    ! However, must pretend that MPL_INIT has actually initialized it, but 
    ! MPL_NUMPROC will not be set
    CALL MPI_INITIALIZED(LLINIT,IRET)
    IF (LLINIT == 0) THEN
      CALL MPI_INIT(IERROR)
      LINITMPI_VIA_MPL = .TRUE.
      CALL EC_MPI_ATEXIT() ! ifsaux/support/endian.c: to make sure MPI_FINALIZE gets called
    ENDIF
  ENDIF

  ! If LMPLUSERCOMM is not set use MPI_COMM_WORLD
  IF (LMPLUSERCOMM) THEN
    ICOMM = MPLUSERCOMM
  ELSE
    ICOMM = MPI_COMM_WORLD
  ENDIF

  CALL MPI_COMM_SIZE(ICOMM,INUMPROC,IERROR)
  CALL MPI_COMM_RANK(ICOMM,IRANK,IERROR)
  IRANK=IRANK+1

  IF (IRANK == 1 .OR. INUMPROC == 1) THEN
    CALL GET_ENVIRONMENT_VARIABLE('MPL_CL_TERMINATE',ENV_CL_TERMINATE)
    IF (ENV_CL_TERMINATE /= ' ') CL_TERMINATE = ENV_CL_TERMINATE
    IARGS = COMMAND_ARGUMENT_COUNT()
    LLCARGS = (IARGS < 0) ! Should be true for non-F90 main programs
    IF (LLCARGS) THEN
      IARGS = IARGC_C()
      LLCARGS = (IARGS >= 0)
      CALL GETARG_C(0,CLARG0) ! The executable name (see ifsaux/support/cargs.c)
    ELSE
      CALL PUTARG_INFO(IARGS, TRIM(CL_TERMINATE)) ! (see ifsaux/support/cargs.c)
      CALL GET_COMMAND_ARGUMENT(0,CLARG0)         ! The executable name (normal F2003 way)
      CALL PUTARG_C(0,TRIM(CLARG0))               ! (see ifsaux/support/cargs.c)
    ENDIF
    IF (IARGS < 0) IARGS = 0
    ALLOCATE(CL_ARGS(0:IARGS))
    N_ARGS = 0
    CL_ARGS(0) = CLARG0
    DO J=1,IARGS ! Other args (repeat until end of loop or terminating argument found)
      IF (LLCARGS) THEN
        CALL GETARG_C(J,CL_ARGS(J))
      ELSE
        CALL GET_COMMAND_ARGUMENT(J,CL_ARGS(J))
        CALL PUTARG_C(J,TRIM(CL_ARGS(J)))
      ENDIF
      IF (CL_ARGS(J) == CL_TERMINATE) EXIT
      N_ARGS = N_ARGS + 1
    ENDDO
  ENDIF

  IF (INUMPROC > 1) THEN
    IROOT = 0
    IARGS = 0
    IF (IRANK == 1) IARGS = N_ARGS
    ! The following broadcast does not use "mailbox" nor attached buffer, both potentially yet to be allocated
    CALL MPI_BCAST(IARGS,1,MPI_INTEGER,IROOT,ICOMM,IERROR)
    ICOUNT = JP_ARGLEN
    IF (IRANK > 1) ALLOCATE(CL_ARGS(0:IARGS))
    IF (IRANK > 1) CALL PUTARG_INFO(IARGS, TRIM(CL_TERMINATE))
    DO J=0,IARGS
     ! The following broadcast does not use "mailbox" nor attached buffer, both potentially yet to be allocated
      CALL MPI_BCAST(CL_ARGS(J),ICOUNT,MPI_BYTE,IROOT,ICOMM,IERROR)
      IF (IRANK > 1) CALL PUTARG_C(J,TRIM(CL_ARGS(J)))
    ENDDO
    IF (IRANK > 1) N_ARGS = IARGS
  ENDIF
ENDIF
END SUBROUTINE INIT_ARGS

SUBROUTINE MPL_GETARG(KARG, CDARG)
INTEGER(KIND=JPIM), INTENT(IN) :: KARG
CHARACTER(LEN=*), INTENT(OUT)  :: CDARG
IF (N_ARGS == -1) CALL INIT_ARGS()
IF (KARG >= 0 .AND. KARG <= N_ARGS) THEN
  CDARG = CL_ARGS(KARG)
ELSE
  CDARG = ' '
ENDIF
END SUBROUTINE MPL_GETARG

FUNCTION MPL_IARGC() RESULT(IRET)
INTEGER(KIND=JPIM) :: IRET
IF (N_ARGS == -1) CALL INIT_ARGS()
IRET = N_ARGS
END FUNCTION MPL_IARGC

END MODULE MPL_ARG_MOD
