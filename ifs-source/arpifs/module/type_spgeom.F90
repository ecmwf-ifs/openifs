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

MODULE TYPE_SPGEOM

! Module for spectral geometry structures.

USE PARKIND1 , ONLY : JPIM, JPRB

!     ------------------------------------------------------------------

IMPLICIT NONE
PRIVATE

SAVE

! GMR     : coefficients for spectral multiplication by GM.
! SCGMAP  : coefficients for multiplication by (GM**2) in spectral space (global model).
! ESCGMAP : coefficients for multiplication by (GM**2) in spectral space (LAM model).

TYPE, PUBLIC :: TSPGEOM
REAL(KIND=JPRB),ALLOCATABLE :: GMR(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SCGMAP(:,:)
REAL(KIND=JPRB) :: ESCGMAP(3)
CONTAINS
FINAL :: TSPGEOM_FINAL
END TYPE TSPGEOM

!!TYPE(TSPGEOM), POINTER,PUBLIC :: YSPGEOM => NULL()

CONTAINS

SUBROUTINE TSPGEOM_FINAL(THIS)
  TYPE(TSPGEOM) :: THIS
  ! If we don't add this, we may get a internal compiler error
END SUBROUTINE TSPGEOM_FINAL

END MODULE TYPE_SPGEOM
