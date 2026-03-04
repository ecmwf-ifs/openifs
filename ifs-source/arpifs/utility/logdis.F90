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

SUBROUTINE LOGDIS(KSTEP,PMODELTIME,PTIME1,PTIME2,LDTIMER,LDDATER,CDSTEP)

!**** *LOGDIS*  - Display the current state of the model execution on stderr and the model listing

!     Purpose.
!     --------
!       LOGFILE DISPLAY :  To display the current state of the model execution

!**   Interface.
!     ----------
!        *CALL* *LOGDIS

!        Explicit arguments :
!        --------------------
!        KSTEP : Model step
!        PMODELTIME : time of the model
!        PTIME1 : cpu start time
!        PTIME2 : cpu end time
!        LDTIMER : .T. to display the time
!        LDDATER : .T. to display the date
!        CDSTEP : 'step' string

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      R. El Khatib *Meteo-France*
!      Original : 02-Apr-2015 from CNT4.

! Modifications
! -------------
!       R. El Khatib 28-Jun-2018 'step' string
! End Modifications
!      ----------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT   ,NULERR
USE YOMMP0   , ONLY : MYPROC

!      ----------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KSTEP 
REAL(KIND=JPRB), INTENT(IN) :: PMODELTIME
REAL(KIND=JPRD), INTENT(IN) :: PTIME1
REAL(KIND=JPRD), INTENT(IN) :: PTIME2
LOGICAL, INTENT(IN) :: LDTIMER
LOGICAL, INTENT(IN) :: LDDATER
CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: CDSTEP

CHARACTER (LEN = 46) ::  CLDAYF
CHARACTER (LEN = 10) ::  CLTIMEOD,CLDAT(3)
CHARACTER (LEN = 8)  ::  CLSTEP

INTEGER(KIND=JPIM) :: IVALUES(8)
INTEGER(KIND=JPIM) :: IHOUR, IMIN

REAL(KIND=JPRB) :: ZSEC, ZDTIME

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "meminfo.intfb.h"

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LOGDIS',0,ZHOOK_HANDLE)
!      ----------------------------------------------------------------

CALL DATE_AND_TIME(CLDAT(1),CLTIMEOD,CLDAT(3),IVALUES)
ZSEC   = PMODELTIME
IHOUR  = INT( ZSEC/3600._JPRB )
IMIN   = INT( (ZSEC-REAL(IHOUR,JPRB)*3600._JPRB)/60._JPRB )
ZDTIME=PTIME2-PTIME1
IF (PRESENT(CDSTEP)) THEN
  CLSTEP=CDSTEP
ELSE
  CLSTEP='STEP'
ENDIF
WRITE (UNIT=CLDAYF,&
 & FMT='(1X,A,'':'',A,'':'',A,'' '',A,I5,'' H='',I4,'':'',I2.2,'' +CPU='',F7.3)')&
 & CLTIMEOD(1:2),CLTIMEOD(3:4),CLTIMEOD(5:6),TRIM(CLSTEP),KSTEP,IHOUR,IMIN,ZDTIME
IF (MYPROC == 1) THEN
  WRITE(NULERR,*) CLDAYF
  CALL MEMINFO(NULERR,KSTEP)
ENDIF
IF (.NOT.LDTIMER) THEN
  ZDTIME=0._JPRB
ENDIF
IF (.NOT.LDDATER) THEN
  CLTIMEOD(1:10)='0000000000'
ENDIF
WRITE (UNIT=CLDAYF,&
 & FMT='(1X,A,'':'',A,'':'',A,'' '',A,I5,'' H='',I4,'':'',I2.2,'' +CPU='',F7.3)')&
 & CLTIMEOD(1:2),CLTIMEOD(3:4),CLTIMEOD(5:6),TRIM(CLSTEP),KSTEP,IHOUR,IMIN,ZDTIME
WRITE (UNIT=NULOUT, FMT='(A42)') CLDAYF

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LOGDIS',1,ZHOOK_HANDLE)
END SUBROUTINE LOGDIS
