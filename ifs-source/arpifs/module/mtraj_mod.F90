!
! Copyright 2011 ECMWF
! 
! This software was developed at ECMWF for evaluation
! and may be used for academic and research purposes only.
! The software is provided as is without any warranty.
! 
! This software can be used, copied and modified but not
! redistributed or sold. This notice must be reproduced
! on each copy made.
!

!> Handle model configuration for the IFS model

MODULE MTRAJ_MOD

USE YOMGFL, ONLY : TGFL
USE YOMGMV, ONLY : TGMV
IMPLICIT NONE
PRIVATE

PUBLIC :: MTRAJ, MTRAJ_CREATE, MTRAJ_DELETE

! ------------------------------------------------------------------------------

TYPE :: MTRAJ
  TYPE(TGMV) :: YRGMV5
  TYPE(TGFL) :: YRGFL5
END TYPE MTRAJ

! ------------------------------------------------------------------------------
CONTAINS
! ------------------------------------------------------------------------------

SUBROUTINE MTRAJ_CREATE(SELF)
USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

TYPE(MTRAJ), INTENT(INOUT) :: SELF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MTRAJ_MOD:MTRAJ_CREATE',0,ZHOOK_HANDLE)

!WRITE(0,*)'MTRAJ_MOD:MTRAJ_CREATE NOT IMPLEMENTED'

IF (LHOOK) CALL DR_HOOK('MTRAJ_MOD:MTRAJ_CREATE',1,ZHOOK_HANDLE)

END SUBROUTINE MTRAJ_CREATE

! ------------------------------------------------------------------------------

SUBROUTINE MTRAJ_DELETE(SELF)
USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE
TYPE(MTRAJ), INTENT(INOUT) :: SELF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MTRAJ_MOD:MTRAJ_DELETE',0,ZHOOK_HANDLE)

!WRITE(0,*)'MTRAJ_MOD:MTRAJ_DELETE NOT IMPLEMENTED'

IF (LHOOK) CALL DR_HOOK('MTRAJ_MOD:MTRAJ_DELETE',1,ZHOOK_HANDLE)

END SUBROUTINE MTRAJ_DELETE

! ------------------------------------------------------------------------------

END MODULE MTRAJ_MOD
