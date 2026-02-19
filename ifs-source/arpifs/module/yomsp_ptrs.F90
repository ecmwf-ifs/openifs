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

MODULE YOMSP_PTRS

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

! Pointers (offsets) for individual fields in spectral arrays

INTEGER(KIND=JPIM) :: MSPTR_T
INTEGER(KIND=JPIM) :: MSPTR_SP

END MODULE YOMSP_PTRS
