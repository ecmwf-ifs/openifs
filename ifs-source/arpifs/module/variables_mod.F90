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

!> Handle variables for the IFS model

MODULE VARIABLES_MOD
use parkind1, only : jpim

IMPLICIT NONE
PRIVATE

PUBLIC :: VARIABLES, VARIABLES_CREATE, VARIABLES_DELETE, VARIABLES_CLONE, &
        & HAS_MODEL_FIELDS, HAS_OLD_FIELDS, IS_LINEAR

! ------------------------------------------------------------------------------

TYPE :: VARIABLES
  LOGICAL :: linit = .false.
  LOGICAL :: lctrl = .false.
  LOGICAL :: linear = .false.
  integer(kind=jpim), pointer :: fieldids(:)=>null()
END TYPE VARIABLES

! ------------------------------------------------------------------------------

interface variables_create
  module procedure variables_create_old, variables_create_new
end interface

! ------------------------------------------------------------------------------

#include "abor1.intfb.h"

! ------------------------------------------------------------------------------
CONTAINS
! ------------------------------------------------------------------------------

SUBROUTINE VARIABLES_CREATE_OLD(SELF, LDCTL, LDLIN)
IMPLICIT NONE
TYPE(VARIABLES), INTENT(INOUT) :: SELF
LOGICAL, INTENT(IN) :: LDCTL
LOGICAL, OPTIONAL, INTENT(IN) :: LDLIN
SELF%LCTRL = LDCTL
IF (PRESENT(LDLIN)) SELF%LINEAR=LDLIN
SELF%LINIT = .true.
END SUBROUTINE VARIABLES_CREATE_OLD

! ------------------------------------------------------------------------------

SUBROUTINE VARIABLES_CREATE_NEW(SELF, KIDS)
IMPLICIT NONE
TYPE(VARIABLES), INTENT(INOUT) :: SELF
integer(kind=jpim), intent(in) :: kids(:)

allocate(self%fieldids(size(kids)))
self%fieldids(:)=kids(:)

self%lctrl = .true.
self%linit = .true.
END SUBROUTINE VARIABLES_CREATE_NEW

! ------------------------------------------------------------------------------

SUBROUTINE VARIABLES_CLONE(SELF, OTHER)
IMPLICIT NONE
TYPE(VARIABLES), INTENT(INOUT) :: SELF
TYPE(VARIABLES), INTENT(IN) :: OTHER
IF (.not.other%linit) call abor1('VARIABLES_CLONE: input not initialized')
IF (ASSOCIATED(other%fieldids)) THEN
  ALLOCATE(self%fieldids(size(other%fieldids)))
  self%fieldids(:) = other%fieldids(:)
ENDIF
self%lctrl = other%lctrl
self%linear = other%linear
self%linit = .true.
END SUBROUTINE VARIABLES_CLONE

! ------------------------------------------------------------------------------

SUBROUTINE VARIABLES_DELETE(SELF)
IMPLICIT NONE
TYPE(VARIABLES), INTENT(INOUT) :: SELF
IF (ASSOCIATED(self%fieldids)) DEALLOCATE(self%fieldids)
self%linit = .false.
END SUBROUTINE VARIABLES_DELETE

! ------------------------------------------------------------------------------

LOGICAL FUNCTION HAS_MODEL_FIELDS(SELF)
IMPLICIT NONE
TYPE(VARIABLES), INTENT(IN) :: SELF
IF (.not.self%linit) call abor1('VARIABLES:HAS_MODEL_FIELDS: input not initialized')
HAS_MODEL_FIELDS = .NOT.self%lctrl
END FUNCTION HAS_MODEL_FIELDS

! ------------------------------------------------------------------------------

LOGICAL FUNCTION HAS_OLD_FIELDS(SELF)
IMPLICIT NONE
TYPE(VARIABLES), INTENT(IN) :: SELF
IF (.not.self%linit) call abor1('VARIABLES:HAS_OLD_FIELDS: input not initialized')
HAS_OLD_FIELDS = .NOT.ASSOCIATED(self%fieldids)
END FUNCTION HAS_OLD_FIELDS

! ------------------------------------------------------------------------------

LOGICAL FUNCTION IS_LINEAR(SELF)
IMPLICIT NONE
TYPE(VARIABLES), INTENT(IN) :: SELF
IF (.not.self%linit) call abor1('VARIABLES:IS_LINEAR: input not initialized')
IS_LINEAR = self%linear
END FUNCTION IS_LINEAR

! ------------------------------------------------------------------------------

END MODULE VARIABLES_MOD
