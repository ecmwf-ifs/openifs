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

MODULE YOMNEMO

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Control variables for the NEMO model

! NEMOCSTEP : Current NEMO time step
! NEMONSTEP : Number of NEMO steps between each coupling

INTEGER(KIND=JPIM) :: NEMOCSTEP
INTEGER(KIND=JPIM) :: NEMONSTEP

!     ------------------------------------------------------------------
END MODULE YOMNEMO

