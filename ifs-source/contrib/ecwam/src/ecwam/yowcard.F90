! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWCARD

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *WAMCARD* - TEXT OF THE WAMINFO FILE.

      INTEGER(KIND=JWIM), PARAMETER :: JPCL=17
      CHARACTER(LEN=72) :: CARD(JPCL) 
      END MODULE YOWCARD
