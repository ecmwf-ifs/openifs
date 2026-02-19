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

MODULE YOMOROG

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Orography: structure TOROG

! OROG       : grid-point surface orography "Phi_s".
! OROGL      : zonal component of "grad Phi_s".
! OROGM      : meridian component of "grad Phi_s".
! OROGLL,OROGMM,OROGLM: second-order reduced derivatives of surface orography.

TYPE TOROG
  REAL(KIND=JPRB), ALLOCATABLE :: OROG(:)
  REAL(KIND=JPRB), ALLOCATABLE :: OROGL(:)
  REAL(KIND=JPRB), ALLOCATABLE :: OROGM(:)
  REAL(KIND=JPRB), ALLOCATABLE :: OROGLL(:)
  REAL(KIND=JPRB), ALLOCATABLE :: OROGMM(:)
  REAL(KIND=JPRB), ALLOCATABLE :: OROGLM(:)
END TYPE TOROG

TYPE TOROG_BLOCKED
  REAL(KIND=JPRB), ALLOCATABLE :: OROG(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: OROGL(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: OROGM(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: OROGLL(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: OROGMM(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: OROGLM(:,:)
END TYPE TOROG_BLOCKED

! ------------------------------------------------------------------

END MODULE YOMOROG
