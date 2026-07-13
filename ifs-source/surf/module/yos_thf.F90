MODULE YOS_THF
 
USE PARKIND1  ,ONLY : JPIM     ,JPRB

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

IMPLICIT NONE

SAVE

REAL(KIND=JPRB) :: RHOH2O    ! DENSITY OF LIQUID WATER.   (RATM/100.)
REAL(KIND=JPRB) :: RVTMP2    ! RVTMP2=RCPV/RCPD-1.
REAL(KIND=JPRB) :: R2ES      ! Constants for computation of esat
REAL(KIND=JPRB) :: R3LES     ! Constants for computation of esat
REAL(KIND=JPRB) :: R3IES     ! Constants for computation of esat
REAL(KIND=JPRB) :: R4LES     ! Constants for computation of esat
REAL(KIND=JPRB) :: R4IES     ! Constants for computation of esat
REAL(KIND=JPRB) :: R5LES     ! Constants for computation of esat
REAL(KIND=JPRB) :: R5IES     ! Constants for computation of esat

END MODULE YOS_THF
