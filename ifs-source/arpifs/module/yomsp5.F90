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

MODULE YOMSP5

USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD

IMPLICIT NONE
SAVE

!     ------------------------------------------------------------------

!*    Spectral arrays
!     SP5(X) for trajectory of (X),

! -   TRAJECTORY

TYPE(SPECTRAL_FIELD) :: SPA5

END MODULE YOMSP5
