! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE YOM_GRID_BICONSERV

USE PARKIND1 , ONLY : JPRB

IMPLICIT NONE

SAVE
! ----------------------------------------------------------------------------

! Global surface ln pressure fields for conserving interpolation of trajectory
! and increments. 
! * RGPPRS_HR    : grid-point high resolution (outer loop) lnps
! * RGPPRS_LR    : grid-point low resolution (inner loop) lnps

REAL(KIND=JPRB),ALLOCATABLE :: RGPPRS_HR (:)
REAL(KIND=JPRB),ALLOCATABLE :: RGPPRS_LR (:)

! ----------------------------------------------------------------------------
END MODULE YOM_GRID_BICONSERV
