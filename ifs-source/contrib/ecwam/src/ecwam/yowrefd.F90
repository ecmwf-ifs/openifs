! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWREFD

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *REFDOT* - DEPTH AND CURRENT PART OF THETA DOT.

      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:)   :: THDD
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:)   :: THDC
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:,:) :: SDOT

      LOGICAL :: LLUPDTTD

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *THDD*      REAL      DEPTH GRADIENT PART OF THETA DOT.
!      *THDC*      REAL      CURRENT GRADIENT PART OF THETA DOT.
!      *SDOT*      SIGMA DOT ARRAY

!      *LLUPDTTD*  LOGICAL  IF TRUE THETA DOT ARRAYS NEED UPDATING.

! ----------------------------------------------------------------------
      END MODULE YOWREFD
