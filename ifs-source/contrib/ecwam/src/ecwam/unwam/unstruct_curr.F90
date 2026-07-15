! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE UNSTRUCT_CURR
  !
  ! Rule of conversion of angle in WAM.
  ! theta_{WAM} = 90 - theta_{trigonometric}
  !
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      IMPLICIT NONE
      REAL(KIND=JWRU), ALLOCATABLE :: CURTXY (:,:)

END MODULE
