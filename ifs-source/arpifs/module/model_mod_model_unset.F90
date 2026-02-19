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

SUBMODULE (MODEL_MOD) MODEL_MOD_MODEL_UNSET
IMPLICIT NONE

CONTAINS

MODULE SUBROUTINE MODEL_UNSET(SELF)
USE PARKIND1,  ONLY : JPRB, JPIM
USE YOMHOOK,   ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE
TYPE(MODEL), INTENT(IN) :: SELF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODEL_MOD:MODEL_UNSET',0,ZHOOK_HANDLE)

IF (LHOOK) CALL DR_HOOK('MODEL_MOD:MODEL_UNSET',1,ZHOOK_HANDLE)
END SUBROUTINE MODEL_UNSET

END SUBMODULE
