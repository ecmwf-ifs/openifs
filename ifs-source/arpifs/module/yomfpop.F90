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

MODULE YOMFPOP

USE YOMFPIOS  , ONLY : TNAMFPIOS
USE TYPE_FAOPH, ONLY : TFAOPH
USE TYPE_FPOFN, ONLY : TFPOFN

IMPLICIT NONE

SAVE

! I/O handling for Fullpos

TYPE TFPIOH

! General I/O control parameters
TYPE(TNAMFPIOS) :: YNAMFPIOS

! FA I/O parameters
TYPE (TFAOPH), ALLOCATABLE :: YFPOPH(:)

! Filenames
TYPE (TFPOFN), ALLOCATABLE :: YFPOFN(:)

END TYPE TFPIOH

END MODULE YOMFPOP
