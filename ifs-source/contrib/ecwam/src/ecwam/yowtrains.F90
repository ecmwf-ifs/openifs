! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWTRAINS

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    SWELL TRAINS PARAMETERS.

      REAL(KIND=JWRB), ALLOCATABLE :: EMTRAIN(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: THTRAIN(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: PMTRAIN(:,:)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *ETRAIN*   REAL      TOTAL ENERGY.
!      *THTRAIN*  REAL      MEAN DIRECTION.
!      *P1TRAIN*  REAL      MEAN PERIOD (-1).
! ----------------------------------------------------------------------
      END MODULE YOWTRAINS
