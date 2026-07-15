! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWCINP

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *CINP* USER INPUT: AREAS TO BE CHANGED

      INTEGER(KIND=JWIM)              :: NOUT 
      REAL(KIND=JWRB), ALLOCATABLE    :: XOUTW(:) 
      REAL(KIND=JWRB), ALLOCATABLE    :: XOUTS(:) 
      REAL(KIND=JWRB), ALLOCATABLE    :: XOUTE(:) 
      REAL(KIND=JWRB), ALLOCATABLE    :: XOUTN(:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: NOUTD(:) 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *NOUT*      INTEGER   NUMBER OF AREAS TO BE ADJUSTED.
!      *XOUTW*     REAL      WESTERN-MOST LONG OF AREA TO BE CHANGED.
!      *XOUTE*     REAL      EASTERN-MOST LONG OF AREA TO BE CHANGED.
!      *XOUTS*     REAL      SOUTHERN-MOST LAT OF AREA TO BE CHANGED.
!      *XOUTN*     REAL      NORTHERN-MOST LAT OF AREA TO BE CHANGED.
!      *NOUTD*     INTEGER   DEPTH IN AREA IN METRES -999 FOR LAND.

! ----------------------------------------------------------------------
      END MODULE YOWCINP
