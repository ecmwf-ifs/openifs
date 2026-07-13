MODULE ABORT_SURF_MOD
CONTAINS
SUBROUTINE ABORT_SURF(CDTEXT)
 
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE MPL_MODULE

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

IMPLICIT NONE

CHARACTER(LEN=*) :: CDTEXT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ABORT_SURF_MOD:ABORT_SURF',0,ZHOOK_HANDLE)

!this is temporary /mats
CALL ABOR1(CDTEXT)
IF (LHOOK) CALL DR_HOOK('ABORT_SURF_MOD:ABORT_SURF',1,ZHOOK_HANDLE)

END SUBROUTINE ABORT_SURF
END MODULE ABORT_SURF_MOD
