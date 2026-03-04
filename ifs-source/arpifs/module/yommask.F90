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

MODULE YOMMASK

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

! ----------------------------------------------------------------------

! Module for On demand Semi-Lagrangian Masks

! * NFIXSFLD: description of fixed fields (fields not sent/received)
!   in on-demand communications.
INTEGER(KIND=JPIM) :: NFIXSFLD(2) 

! ----------------------------------------------------------------------

END MODULE YOMMASK

