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

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMCT0   , ONLY : LSMSSIG  , CMETER
USE YOMMP0   , ONLY : MYPROC

!      ----------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KTIME

CHARACTER (LEN = 40) ::  CLSETEV
CHARACTER (LEN = 256) ::  CLSMSNAME,CLECFNAME,CLMULTIO

CHARACTER (LEN = 7), PARAMETER :: CL_CPENV='SMSNAME'
CHARACTER (LEN = 8), PARAMETER :: CL_CPENV_ECF='ECF_NAME'
CHARACTER (LEN = 19), PARAMETER :: CL_CPENV_MULTIO='MULTIO_NOTIFY_FLUSH'

INTEGER(KIND=JPIM) :: ICPLEN,ICPLEN_ECF, IPPTR, ISTAT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SIGPOST',0,ZHOOK_HANDLE)
!      ----------------------------------------------------------------
  !OIFS/FC-only - Call to uti_cgetenv removed because routine is part of odb, may cause issues


!*       3.15   Signal SMS event for completion of post_processing

IF(LSMSSIG) THEN
  CALL GSTATS(1940,0)
  IF(MYPROC == 1) THEN
    IPPTR=INT(REAL(KTIME,JPRB)/3600._JPRB)
    WRITE(CLSETEV,' (A25,'' step '',I8,''&'') ') CMETER,IPPTR
    CLSMSNAME="                                             "
    CLECFNAME="                                             "
    CALL GET_ENVIRONMENT_VARIABLE(NAME=CL_CPENV, VALUE=CLSMSNAME, LENGTH=ICPLEN, STATUS=ISTAT)
    IF (ISTAT /= 0) THEN
      CLSMSNAME='NOSMS'
      ICPLEN=5
    ENDIF
    CALL GET_ENVIRONMENT_VARIABLE(NAME=CL_CPENV_ECF, VALUE=CLECFNAME, LENGTH=ICPLEN_ECF, STATUS=ISTAT)
    IF (ISTAT /= 0) THEN
      CLECFNAME='NOECF'
      ICPLEN=5
    ENDIF
    IF ((ICPLEN > 0.AND.CLSMSNAME(1:5) /= 'NOSMS') .OR.  &
       &(ICPLEN_ECF > 0.AND.CLECFNAME(1:5) /= 'NOECF')    ) THEN
      CALL SYSTEM(CLSETEV)
      WRITE(UNIT=NULOUT,FMT='(A25,I8,'' posted '')') CMETER,IPPTR
    ELSE
      WRITE(UNIT=NULOUT,&
       & FMT='(A25,I8,&
       & '' not posted  because neither SMSNAME nor ECF_NAME  is defined. LEN='',&
       & I3,A16,A16)') CMETER,IPPTR, ICPLEN, CLSMSNAME, CLECFNAME
    ENDIF
  ENDIF
  CALL GSTATS(1940,1)
ENDIF

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SIGPOST',1,ZHOOK_HANDLE)
END SUBROUTINE SIGPOST
