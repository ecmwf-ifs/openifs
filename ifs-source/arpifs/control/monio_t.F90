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

SUBROUTINE MONIO_T(KSTART,YDRIP,K___TS,KN1___,KNFR___,KN___TS,KN___TSMIN,LDCOND,LDACTIVE)

!   Purpose:
!   --------
!    ??????

!   Author:
!   -------
!    ??????

!   Modifications:
!   --------------
!      K. Yessad (July 2014): Move some variables.
!      R. El Khatib  04-Jun-2018 refactor suct1 against monio
!--------------------------------------------------------------------------------

USE PARKIND1, ONLY : JPRB, JPIM
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMRIP  , ONLY : TRIP

!--------------------------------------------------------------------------------

IMPLICIT NONE

INTEGER (KIND=JPIM), INTENT (IN)           :: KSTART           ! lower bound of K___TS
TYPE(TRIP)         , INTENT (IN)           :: YDRIP
INTEGER (KIND=JPIM), INTENT (INOUT)        :: K___TS (KSTART:) !  KPOSTS
INTEGER (KIND=JPIM), INTENT (IN)           :: KN1___           !  N1POS
INTEGER (KIND=JPIM), INTENT (IN)           :: KNFR___          !  NFRPOS
INTEGER (KIND=JPIM), INTENT (IN)           :: KN___TS (0:)     !  NPOSTS
INTEGER (KIND=JPIM), INTENT (IN), OPTIONAL :: KN___TSMIN (0:)  !  NPOSTSMIN
LOGICAL,             INTENT (IN), OPTIONAL :: LDCOND
LOGICAL,             INTENT (IN), OPTIONAL :: LDACTIVE

!--------------------------------------------------------------------------------

REAL (KIND=JPRB), PARAMETER :: ZUNIT = 3600._JPRB

INTEGER (KIND=JPIM) :: JS, ISS, I1, IEND
INTEGER (KIND=JPIM) :: ILS
REAL (KIND=JPRB) :: ZISS, ZISSX, ZMIN
REAL (KIND=JPRB) :: ZTSTEP, ZLS
LOGICAL          :: LLCOND, LLACTIVE

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MONIO_T',0,ZHOOK_HANDLE)
ASSOCIATE(NSTART=>YDRIP%NSTART,NSTOP=>YDRIP%NSTOP, TSTEP=>YDRIP%TSTEP)
!--------------------------------------------------------------------------------

IEND=UBOUND(K___TS,1)

K___TS(:)=0

LLCOND = .TRUE.
IF (PRESENT (LDCOND)) LLCOND = LDCOND

LLACTIVE = .TRUE.
IF (PRESENT (LDACTIVE)) LLACTIVE = LDACTIVE

IF (LLACTIVE) THEN
  I1=KN1___
ELSE
  I1=0
ENDIF

ZTSTEP = ABS (TSTEP)
IF(ZTSTEP < TINY (ZLS))THEN
  ZTSTEP = TINY (ZLS)
ENDIF
ZISSX = REAL (HUGE (ILS), JPRB)

IF(KN___TS(0) >= 1)THEN
  DO JS=1,KN___TS(0)
    ISS=KN___TS(JS)*KNFR___
    IF(ISS >= NSTART.AND.ISS <= NSTOP .AND. KN___TS(JS) >= KSTART .AND. KN___TS(JS) <= IEND )THEN
      K___TS(KN___TS(JS))=I1
    ENDIF
  ENDDO                  
ELSEIF(KN___TS(0) <= -1)THEN
  DO JS=1,-KN___TS(0)
    IF (PRESENT (KN___TSMIN)) THEN
      ZMIN=REAL (KN___TSMIN(JS), JPRB)
    ELSE
      ZMIN=0._JPRB
    ENDIF
    ZISS=ZUNIT*REAL((-KN___TS(JS)+ZMIN/60._JPRB),JPRB)*REAL(KNFR___,JPRB)/ZTSTEP
    ZISS=MIN(ZISSX,ZISS)
    ISS=NINT(ZISS)
    IF(ISS >= NSTART.AND.ISS <= NSTOP .AND. ISS/MAX(1,KNFR___) >= KSTART .AND. ISS/MAX(1,KNFR___) <= IEND )THEN
      K___TS(ISS/MAX(1,KNFR___))=I1
    ENDIF
  ENDDO
ELSEIF (LLCOND) THEN
  DO JS=NSTART,NSTOP
    IF(MOD(JS,MAX(KNFR___,1)) == 0 .AND. JS/MAX(1,KNFR___) >= KSTART .AND. JS/MAX(1,KNFR___) <= IEND )THEN
      K___TS(JS/MAX(1,KNFR___))=I1
    ENDIF
  ENDDO
ENDIF

!--------------------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('MONIO_T',1,ZHOOK_HANDLE)
END SUBROUTINE MONIO_T

