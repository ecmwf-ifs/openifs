! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE DR_HOOK_STACKCHECK_MOD

! Used by dr_hook_util to monitor thread stack usage 

USE EC_PARKIND  ,ONLY : JPIB

IMPLICIT NONE

SAVE

PRIVATE :: JPIB

PUBLIC

INTEGER(KIND=JPIB), ALLOCATABLE :: ISAVE(:) 
INTEGER(KIND=JPIB), ALLOCATABLE :: IMAXSTACK(:) 
LOGICAL,            ALLOCATABLE :: LL_THREAD_FIRST(:)

END MODULE DR_HOOK_STACKCHECK_MOD

