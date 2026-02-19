! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

MODULE VARBC_PRED
  ! Header and Type definition removed for OpenIFS and forecast only
  ! Retain declaration of JPREDNAME since required for CLASS_VARBC Type
  ! declaration
  USE PARKIND1   , ONLY : JPIM

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC JPREDNAME

  INTEGER(KIND=JPIM), PARAMETER :: JPREDNAME = 21
!-----------------------------------------------------------------------
END MODULE VARBC_PRED
