! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

PROGRAM MASTER1C

!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE OML_MOD   ,ONLY : OML_INIT
!USE YOMERRTRAP

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
#include "cnt1c.intfb.h"
!     ------------------------------------------------------------------

! DR_HOOK doesn't work with SCM (the relevant code is not complete there).
LHOOK=.false.

IF (LHOOK) CALL DR_HOOK('MASTER1C',0,ZHOOK_HANDLE)


!     ------------------------------------------------------------------
!     Error trapping from master.F90.
!     ------------------------------------------------------------------

!CALL SET_ERR_TRAP


!     ------------------------------------------------------------------
!     Call model.
!     ------------------------------------------------------------------

CALL OML_INIT()
CALL CNT1C


IF (LHOOK) CALL DR_HOOK('MASTER1C',1,ZHOOK_HANDLE)
STOP
END PROGRAM MASTER1C



