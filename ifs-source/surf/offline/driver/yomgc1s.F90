MODULE YOMGC1S
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE
SAVE

!     ------------------------------------------------------------------
!*    Grid point boundaries

! -   VARIABLES

!       RCORI   -  Coriolis parameter 
!       GEMU    -  SIN of latitude on the real earth
!       GSQM2   -  COS of latitude on the real earth
!       GELAM   -  longitude on the real earth
!       GELAT   -  latitude on the real earth
!       GECLO   -  COS of longitude on the real earth
!       GESLO   -  SIN of longitude on the real earth
!       LMASK   -  Logical Mask Array of original input field

REAL(KIND=JPRB),ALLOCATABLE :: RCORI(:)
REAL(KIND=JPRB),ALLOCATABLE :: GEMU(:)
REAL(KIND=JPRB),ALLOCATABLE :: GSQM2(:)
REAL(KIND=JPRB),ALLOCATABLE :: GELAM(:)
REAL(KIND=JPRB),ALLOCATABLE :: GELAT(:)
REAL(KIND=JPRB),ALLOCATABLE :: GECLO(:)
REAL(KIND=JPRB),ALLOCATABLE :: GESLO(:)
LOGICAL,ALLOCATABLE :: LMASK(:)

!      ----------------------------------------------------------------
END MODULE YOMGC1S
