! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE ABOR1FL(CDFILE, KLINENUM, CDTEXT)
  !! Abort that prints file, line, and message
  !! Tracebacks will be printed if possible
  !! All processes will be terminated in parallel MPI context

USE EC_PARKIND  ,ONLY : JPIM
USE EC_LUN   ,ONLY : NULOUT, NULERR
USE MPL_MODULE    ,ONLY : MPL_ABORT
USE OML_MOD       ,ONLY : OML_MY_THREAD
IMPLICIT NONE

CHARACTER(LEN=*),   INTENT(IN) :: CDFILE
INTEGER(KIND=JPIM), INTENT(IN) :: KLINENUM
CHARACTER(LEN=*),   INTENT(IN) :: CDTEXT

IF (LEN(CDFILE) > 0 .AND. KLINENUM > 0) THEN
  IF (NULOUT >= 0 .AND. NULOUT /= 6 .AND. NULOUT /= NULERR) THEN
    WRITE(NULOUT,'(A,I0,A,I0,A,A,A,I0,A,A)') 'ABOR1     [PROC=',MYPROC(),',THRD=',OML_MY_THREAD(),'] from [',CDFILE,' +',KLINENUM,'] : ', CDTEXT
  ENDIF
  IF (NULERR >= 0) THEN
    WRITE(NULERR,'(A,I0,A,I0,A,A,A,I0,A,A)') 'ABOR1     [PROC=',MYPROC(),',THRD=',OML_MY_THREAD(),'] from [',CDFILE,' +',KLINENUM,'] : ', CDTEXT
  ENDIF
ELSE
  IF (NULOUT >= 0 .AND. NULOUT /= 6 .AND. NULOUT /= NULERR) THEN
    WRITE(NULOUT,'(A,I0,A,I0,A,A)') 'ABOR1     [PROC=',MYPROC(),',THRD=',OML_MY_THREAD(),'] : ', CDTEXT
  ENDIF
  IF (NULERR >= 0 ) THEN
    WRITE(NULERR,'(A,I0,A,I0,A,A)') 'ABOR1     [PROC=',MYPROC(),',THRD=',OML_MY_THREAD(),'] : ', CDTEXT
  ENDIF
ENDIF

IF (NULOUT >= 0) THEN
  CALL EC_FLUSH(NULOUT)
  IF (NULOUT /= 0 .AND. NULOUT /= 6) CLOSE(NULOUT)
ENDIF
CALL EC_FLUSH(NULERR)

IF(LEN(CDTEXT) <= 512) THEN
  CALL MPL_ABORT(CDTEXT)
ELSE
  CALL MPL_ABORT
ENDIF

CONTAINS

FUNCTION MYPROC() RESULT(IPROC)
  USE MPL_MPIF
  IMPLICIT NONE
  INTEGER(KIND=JPIM) :: IERROR,IPROC
  LOGICAL :: LMPI_INITIALIZED
  IPROC = 1
  CALL MPI_INITIALIZED(LMPI_INITIALIZED,IERROR) ! always thread safe, see standard !
  IF( LMPI_INITIALIZED ) THEN
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,IPROC,IERROR) ! always thread safe, see standard !
    IPROC = IPROC+1 ! 1-based in IFS context
  ENDIF
END FUNCTION

END SUBROUTINE ABOR1FL


SUBROUTINE ABOR1(CDTEXT)
  !! Abort that prints message without file and line number
  !! Delegates to ABOR1FL
  !! Tracebacks will be printed if possible
  !! All processes will be terminated in parallel MPI context

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: CDTEXT
CALL ABOR1FL("",0,CDTEXT)
END SUBROUTINE ABOR1


SUBROUTINE ABOR1_EXCEPTION_HANDLER()
  !! This routine, when registered as the fckit exception handler, will be called
  !! whenever any C++ exception is thrown. The exception is intercepted and can
  !! be inquired through the variable FCKIT_EXCEPTION.
  !! An exception can also be thrown within Fortran:
  !!    CALL FCKIT_EXCEPTION%ABORT("I have my reasons")

#ifdef WITH_FCKIT
USE FCKIT_MODULE, ONLY : FCKIT_EXCEPTION
IF( FCKIT_EXCEPTION%LOCATION%IS_SET() ) then 
  CALL ABOR1FL( FCKIT_EXCEPTION%LOCATION%FILE(), FCKIT_EXCEPTION%LOCATION%LINE(), FCKIT_EXCEPTION%WHAT() )
ELSE
  CALL ABOR1( FCKIT_EXCEPTION%WHAT() )
ENDIF
#else
CALL ABOR1( "An unknown exception is handled via ABOR1_EXCEPTION_HANDLER. Compile fiat with fckit to get a better error message" )
#endif
END SUBROUTINE

SUBROUTINE SET_ABOR1_EXCEPTION_HANDLER() BIND(C,NAME="set_abor1_exception_handler")
  !! This routine registers ABOR1 as the C++ exception handler, will be called
  !! whenever any C++ exception is thrown and not caught.
  !! The exception is intercepted and can be inquired through
  !! the variable FCKIT_EXCEPTION.
  !! An exception can also be thrown within Fortran:
  !!    CALL FCKIT_EXCEPTION%ABORT("I have my reasons")

#ifdef WITH_FCKIT
USE FCKIT_MODULE, ONLY : FCKIT_EXCEPTION, FCKIT_EXCEPTION_HANDLER
EXTERNAL :: ABOR1_EXCEPTION_HANDLER
PROCEDURE(FCKIT_EXCEPTION_HANDLER), POINTER :: FUNPTR
FUNPTR => ABOR1_EXCEPTION_HANDLER
CALL FCKIT_EXCEPTION%SET_HANDLER( FUNPTR )
#endif
END SUBROUTINE
