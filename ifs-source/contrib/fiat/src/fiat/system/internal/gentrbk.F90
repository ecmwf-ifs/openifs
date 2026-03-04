! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!-- Generic traceback calls here

SUBROUTINE GENTRBK_DUMMY
END SUBROUTINE GENTRBK_DUMMY

SUBROUTINE INTEL_TRBK()
#ifdef __INTEL_COMPILER
USE IFCORE
#endif
USE MPL_DATA_MODULE, ONLY : MPL_RANK
LOGICAL :: DONE_TRACEBACK = .FALSE.
INTEGER :: MYPROC,MYTHREAD
CHARACTER(LEN=512) :: CLTRBK
CHARACTER(LEN=512) :: MESSAGE
#ifdef _OPENMP
INTEGER,EXTERNAL :: OMP_GET_THREAD_NUM
#endif

IF(DONE_TRACEBACK) THEN
  WRITE(0,*) "INTEL_TRBK already called"
  RETURN
ENDIF

MYPROC=MPL_RANK
#ifdef _OPENMP
MYTHREAD=OMP_GET_THREAD_NUM() + 1
#else
MYTHREAD=1
#endif

#ifdef __INTEL_COMPILER
WRITE(MESSAGE,'(A,I4,A,I2,A)') &
  &           "Process ",MYPROC," thread ",MYTHREAD, &
  &           " calling tracebackqq from intel_trbk()"
CALL TRACEBACKQQ(TRIM(MESSAGE), USER_EXIT_CODE=-1)
CALL GET_ENVIRONMENT_VARIABLE("EC_LINUX_TRBK",CLTRBK)
#else
  CLTRBK = '1'
#endif
IF( CLTRBK == '1' ) THEN
  WRITE(0,*) "Process ",MYPROC," thread ",MYTHREAD, &
    &        " calling linux_trbk from intel_trbk()"
  CALL LINUX_TRBK() ! See linuxtrbk.c
ENDIF
DONE_TRACEBACK=.TRUE.
END SUBROUTINE INTEL_TRBK
