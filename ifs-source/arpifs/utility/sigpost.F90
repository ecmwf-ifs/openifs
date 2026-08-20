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

SUBROUTINE SIGPOST(KTIME)

!**** *SIGPOST*  - Post events to signal completion of I/Os

!     Purpose.
!     --------
!       Post  events to signal completion of I/Os

!**   Interface.
!     ----------
!        *CALL* *SIGPOST

!        Explicit arguments :
!        --------------------
!        KTIME : time of the model (in s.)

!        Implicit arguments :
!        --------------------
!        None

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
!      R. El Khatib *Meteo-France
!      Original : 02-Apr-2015 from CNT4.

! Modifications
! -------------
!      R. El Khatib 10-Dec-2015 KSTEP in argument (OOPS)
!      R. El Khatib : 23-Aug-2016 fullpos-arpege stamp file moved away 
! End Modifications
!      ----------------------------------------------------------------

USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN       , ONLY : NULOUT, NULERR
USE YOMCT0       , ONLY : LSMSSIG  , CMETER
USE YOMMP0       , ONLY : MYPROC
USE ECFLOW_LIGHT , ONLY : ECFLOW_LIGHT_UPDATE_METER

!      ----------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KTIME

CHARACTER (LEN = 256) ::  CLMULTIO

CHARACTER (LEN = 19), PARAMETER :: CL_CPENV_MULTIO='MULTIO_NOTIFY_FLUSH'

INTEGER(KIND=JPIM) :: ICPLEN, IPPTR

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

CHARACTER(LEN = 64) :: METER_NAME
INTEGER(KIND=JPIM)  :: METER_VALUE
INTEGER(KIND=JPIM)  :: ERROR
CHARACTER(LEN = 16) :: ERROR_STR

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SIGPOST',0,ZHOOK_HANDLE)
!      ----------------------------------------------------------------

CALL GET_ENVIRONMENT_VARIABLE(CL_CPENV_MULTIO, CLMULTIO, LENGTH=ICPLEN)
IF(ICPLEN==0) CLMULTIO = 'NOMULTIO'
! CLMULTIO="                                             "
! CALL UTIL_CGETENV(CL_CPENV_MULTIO, 'NOMULTIO', CLMULTIO, ICPLEN)
IF (ICPLEN > 0.AND.CLMULTIO(1:8) /= 'NOMULTIO' ) THEN
  IF (CLMULTIO(1:8) /= 'JUSTBARR' ) THEN
     IF(MYPROC == 1) THEN
        IPPTR=INT(REAL(KTIME,JPRB)/3600._JPRB)
        CALL IMULTIO_NOTIFY_STEP(IPPTR)
     ENDIF
  ENDIF
ENDIF

!*       3.15   Signal SMS event for completion of post_processing

IF(LSMSSIG) THEN
  CALL GSTATS(1940,0)
  IF(MYPROC == 1) THEN
    IPPTR=INT(REAL(KTIME,JPRB)/3600._JPRB)

    METER_NAME = "step"
    METER_VALUE = IPPTR

    WRITE(NULERR,FMT='("Setting task meter """,A,""" to value """,I0,"""")') TRIM(METER_NAME),METER_VALUE
    ERROR = ECFLOW_LIGHT_UPDATE_METER(METER_NAME, METER_VALUE)

    WRITE(ERROR_STR, *) ERROR
    WRITE(NULERR,FMT='("Update task meter finished, with result """,A,"""")') ERROR_STR
  ENDIF
  CALL GSTATS(1940,1)
ENDIF

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SIGPOST',1,ZHOOK_HANDLE)
END SUBROUTINE SIGPOST
