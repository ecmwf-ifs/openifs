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

MODULE DATETIME_TMP_MOD

! Temporary datetime derived type until the one from OOPS or fckit can be used
! in the IFS. We will not attempt to re-implement all functionality here.

USE PARKIND1, ONLY : JPIM, JPRB, JPRD
USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

PRIVATE
PUBLIC DATETIME_TMP, SET_DATETIME, GET_DATETIME, SETRIP0, GETRIP0, &
     & SETRIPSTEP, GETRIPSTEP

TYPE DATETIME_TMP
  INTEGER(KIND=JPIM) :: NDATE  ! date in the form YYYYMMDD
  INTEGER(KIND=JPIM) :: NSECS  ! time in seconds (e.g. for 12h, 43200)
END TYPE DATETIME_TMP

!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------

SUBROUTINE SET_DATETIME(dt, yy, mm, dd, hh, nn, ss)
IMPLICIT NONE
TYPE(DATETIME_TMP), intent(inout) :: dt
INTEGER(KIND=JPIM), intent(in) :: yy, mm, dd, hh, nn, ss
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DATETIME_TMP_MOD:SET_DATETIME',0,ZHOOK_HANDLE)

dt%ndate = yy * 10000 + mm * 100 + dd
dt%nsecs = hh * 3600 + nn * 60 + ss

IF (LHOOK) CALL DR_HOOK('DATETIME_TMP_MOD:SET_DATETIME',1,ZHOOK_HANDLE)
END SUBROUTINE SET_DATETIME

!-----------------------------------------------------------------------

SUBROUTINE GET_DATETIME(dt, yy, mm, dd, hh, nn, ss)
IMPLICIT NONE
TYPE(DATETIME_TMP), intent(in) :: dt
INTEGER(KIND=JPIM), intent(out) :: yy, mm, dd, hh, nn, ss
INTEGER(KIND=JPIM) :: ii
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#include "fcttim.func.h"

IF (LHOOK) CALL DR_HOOK('DATETIME_TMP_MOD:GET_DATETIME',0,ZHOOK_HANDLE)

dd=NDD(dt%ndate)
mm=NMM(dt%ndate)
yy=NCCAA(dt%ndate)

hh = dt%nsecs / 3600
ii = dt%nsecs - 3600 * hh
nn = ii / 60
ss = ii - 60 * nn

IF (LHOOK) CALL DR_HOOK('DATETIME_TMP_MOD:GET_DATETIME',1,ZHOOK_HANDLE)
END SUBROUTINE GET_DATETIME

!-----------------------------------------------------------------------

SUBROUTINE SETRIP0(ydt0)
USE YOMCST , ONLY : RDAY
USE YOMRIP0, ONLY : NINDAT, NSSSSS, RTIMST
USE YOMCT3,  ONLY : NSTEP

IMPLICIT NONE
TYPE(DATETIME_TMP), intent(in) :: ydt0

INTEGER(KIND=JPIM) :: IA, ID, IM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#include "fcttim.func.h"

IF (LHOOK) CALL DR_HOOK('DATETIME_TMP_MOD:SETRIP0',0,ZHOOK_HANDLE)

NINDAT=ydt0%NDATE
NSSSSS=ydt0%NSECS
ID=NDD(NINDAT)
IM=NMM(NINDAT)
IA=NCCAA(NINDAT)
RTIMST=RTIME(IA,IM,ID,NSSSSS,RDAY)
NSTEP=0

IF (LHOOK) CALL DR_HOOK('DATETIME_TMP_MOD:SETRIP0',1,ZHOOK_HANDLE)
END SUBROUTINE SETRIP0

!-----------------------------------------------------------------------

SUBROUTINE GETRIP0(ydt0)
USE YOMRIP0, ONLY : NINDAT, NSSSSS

IMPLICIT NONE
TYPE(DATETIME_TMP), intent(inout) :: ydt0
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#include "fcttim.func.h"

IF (LHOOK) CALL DR_HOOK('DATETIME_TMP_MOD:GETRIP0',0,ZHOOK_HANDLE)

ydt0%NDATE=NINDAT
ydt0%NSECS=NSSSSS

IF (LHOOK) CALL DR_HOOK('DATETIME_TMP_MOD:GETRIP0',1,ZHOOK_HANDLE)
END SUBROUTINE GETRIP0

!-----------------------------------------------------------------------

SUBROUTINE SETRIPSTEP(PTSTEP,YDT0,KSECS)
USE YOMCT3, ONLY : NSTEP

IMPLICIT NONE
REAL(KIND=JPRB)   , INTENT(IN) :: PTSTEP
TYPE(DATETIME_TMP), INTENT(IN) :: YDT0
INTEGER(KIND=JPIM), INTENT(IN) :: KSECS

REAL(KIND=JPRB) :: ZSECS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#include "fcttim.func.h"

IF (LHOOK) CALL DR_HOOK('DATETIME_TMP_MOD:SETRIPSTEP',0,ZHOOK_HANDLE)

CALL SETRIP0(YDT0)
ZSECS = REAL(KSECS,JPRB)
NSTEP = NINT(ZSECS/PTSTEP)

IF (LHOOK) CALL DR_HOOK('DATETIME_TMP_MOD:SETRIPSTEP',1,ZHOOK_HANDLE)
END SUBROUTINE SETRIPSTEP

!-----------------------------------------------------------------------

SUBROUTINE GETRIPSTEP(PTSTEP,YDT0,KSECS)
USE YOMCT3, ONLY : NSTEP

IMPLICIT NONE
REAL(KIND=JPRB)   , INTENT(IN)    :: PTSTEP
TYPE(DATETIME_TMP), INTENT(INOUT) :: YDT0
INTEGER(KIND=JPIM), INTENT(OUT)   :: KSECS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#include "fcttim.func.h"

IF (LHOOK) CALL DR_HOOK('DATETIME_TMP_MOD:GETRIPSTEP',0,ZHOOK_HANDLE)

CALL GETRIP0(YDT0)
KSECS=NINT(PTSTEP)*NSTEP

IF (LHOOK) CALL DR_HOOK('DATETIME_TMP_MOD:GETRIPSTEP',1,ZHOOK_HANDLE)
END SUBROUTINE GETRIPSTEP

!-----------------------------------------------------------------------

END MODULE DATETIME_TMP_MOD
