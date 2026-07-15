! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWCOER

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *COERS* - CONTROLS INFORMATION OF ERS OUTPUT.          

      INTEGER(KIND=JWIM)              :: NERS
      INTEGER(KIND=JWIM)              :: IDELERS 
      INTEGER(KIND=JWIM)              :: IERS 
      INTEGER(KIND=JWIM), ALLOCATABLE :: IJERS(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: IGERS(:) 

      CHARACTER(LEN=14)               :: CDTERS 

      END MODULE YOWCOER
