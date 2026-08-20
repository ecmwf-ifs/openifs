! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWMEAN

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWDRVTYPE, ONLY : INTGT_PARAM_FIELDS

      IMPLICIT NONE

!*    **  INTEGRATED PARAMETERS OF A BLOCK.

      TYPE(INTGT_PARAM_FIELDS) :: INTFLDS

      END MODULE YOWMEAN
