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

MODULE YOMVCGL

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE CONTROL_VECTORS_BASE_MIX , ONLY : CONTROL_VECTOR

IMPLICIT NONE
SAVE

!     ------------------------------------------------------------------

!*     Global arrays for intermediate results of the forecast error calculation
!      or used in the CONGRAD minimization (preconditioner).

! YVCGLPC: eigenvectors (from an earlier minimization)
!          that are used to construct the preconditioner.
! RCGLPC : eigenvalues (from an earlier minimization)
!          that are used to construct the preconditioner.
! NVCGLPC: the number of eigenpairs used to form the preconditioner.
! YVCGLEV: eigenvectors for the current minimization.
! RCGLEV : eigenvalues for the current minimization.
! NVCGLEV: the number of eigenpairs for the current minimization.
! NSAVEPC: Number of PC to be saved (for diagnostics).
! NSAVEEV: Number of EV to be saved (for diagnostics).

#ifdef RS6K
! xlf90 bug
TYPE(CONTROL_VECTOR), ALLOCATABLE, DIMENSION(:), SAVE :: YVCGLPC
TYPE(CONTROL_VECTOR), ALLOCATABLE, DIMENSION(:), SAVE :: YVCGLEV
#else
TYPE(CONTROL_VECTOR), ALLOCATABLE, DIMENSION(:) :: YVCGLPC
TYPE(CONTROL_VECTOR), ALLOCATABLE, DIMENSION(:) :: YVCGLEV
#endif
REAL(KIND=JPRB), ALLOCATABLE :: RCGLPC(:)
INTEGER(KIND=JPIM) :: NVCGLPC
INTEGER(KIND=JPIM) :: NSAVEPC
REAL(KIND=JPRB), ALLOCATABLE :: RCGLEV(:)
INTEGER(KIND=JPIM) :: NVCGLEV
INTEGER(KIND=JPIM) :: NSAVEEV

!     ------------------------------------------------------------------

END MODULE YOMVCGL
