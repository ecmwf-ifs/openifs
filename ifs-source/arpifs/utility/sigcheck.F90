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

SUBROUTINE SIGCHECK(LDSIGSTOP,LDSIGREST)

!**** *SIGCHECK* - Signal check

!     Purpose.
!     --------
!        Check for external signals requesting creation of restart files

!**   Interface.
!     ----------
!        *CALL* *SIGCHECK(...)*

!        Explicit arguments : (output)
!        --------------------  
!                LDSIGSTOP: set .TRUE. if stop signal received
!                LDSIGREST: set .TRUE. if restart file signal received

!        Implicit arguments :  MYPROC
!        --------------------

!     Method.
!     -------      Signals are unblocked.
!                  Master process checks flags.
!                  They may be set by signal processor (see IFSSIG)
!                  Values of flags are sent to all other processes

!     Externals.   IFSSIGB, MPL_SEND, MPL_RECV, MPL_BARRIER
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      David Dent  *ECMWF*
!      Original : 96-07-29

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NPROC, NPRCIDS, MYPROC
USE YOMTAG   , ONLY : MTAGSIG
USE YOMRES   , ONLY : NFLSTOP, NFLREST
USE YOMLUN   , ONLY : NULOUT
USE MPL_MODULE

IMPLICIT NONE

LOGICAL           ,INTENT(OUT)   :: LDSIGSTOP 
LOGICAL           ,INTENT(OUT)   :: LDSIGREST 
INTEGER(KIND=JPIM) :: IBUF(2)

INTEGER(KIND=JPIM) ::  IPROC, ITAG
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SIGCHECK',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!         unblock signals

CALL IFSSIGB

!         check for arrival of signals

LDSIGSTOP=.FALSE.

CALL GSTATS(658,0)
IF(MYPROC == 1) THEN
!        master process checks flags
  IF(NFLREST == 1) THEN
    WRITE(NULOUT,'(A)')' SIGNAL RECEIVED: WRITE RESTART FILES'
  ENDIF
  IF(NFLSTOP == 1) THEN
    WRITE(NULOUT,'(A)')' SIGNAL RECEIVED: STOP'
  ENDIF
  IBUF(1)=NFLSTOP
  IBUF(2)=NFLREST
ENDIF

ITAG=MTAGSIG
IPROC=1
CALL MPL_BROADCAST(IBUF,KTAG=ITAG,KROOT=NPRCIDS(IPROC), &
 & CDSTRING='SIGCHECK:')
NFLSTOP=IBUF(1)
NFLREST=IBUF(2)
IF(NFLREST == 1) THEN
  LDSIGREST=.TRUE.
ENDIF
IF(NFLSTOP == 1) THEN
  LDSIGSTOP=.TRUE.
ENDIF

CALL GSTATS(658,1)
!       reset flag in case execution continues
NFLREST=0
IF(NPROC > 1) THEN
  CALL GSTATS(718,0)
  CALL MPL_BARRIER(CDSTRING='SIGCHECK:')
  CALL GSTATS(718,1)
ENDIF

IF (LHOOK) CALL DR_HOOK('SIGCHECK',1,ZHOOK_HANDLE)
END SUBROUTINE SIGCHECK
