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

MODULE YOMSP

USE PARKIND1, ONLY : JPRB
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Spectral arrays
!     SP(X) for (X)

! -   "EXTRA SPECTRAL ARRAYS FOR SAVING PURPOSES"

TYPE(SPECTRAL_FIELD) :: SPAINI

!     ------------------------------------------------------------------

!*    Spectral arrays spatially filtered
!     SP(X) for (X)

! -   STATE VARIABLES OF THE MODEL

REAL(KIND=JPRB), POINTER :: SPVOR_FLT(:,:) => NULL()! Vorticity (rot)  related filtered.
REAL(KIND=JPRB), POINTER :: SPDIV_FLT(:,:) => NULL()! Divergence (div) related filtered.
REAL(KIND=JPRB), POINTER :: SPUB_FLT(:)    => NULL()! mean value of gradpsl
REAL(KIND=JPRB), POINTER :: SPVB_FLT(:)    => NULL()! mean value of gradpsm
END MODULE YOMSP
