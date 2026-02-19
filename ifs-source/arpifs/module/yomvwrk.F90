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

MODULE YOMVWRK

USE PARKIND1, ONLY : JPIM

IMPLICIT NONE

SAVE

!   -----------------------------------------------------------------

!    NTRSLTYPE       : trajectory storage method -
!                       =0 no storage (only works for LELTRA=T)
!                       =1 minimum storage
!                       =2 more storage, less computation
!       ** N.B. Only NTRSLTYPE=2 is supported in the 3-d model **

INTEGER(KIND=JPIM) :: NTRSLTYPE

!   ----------------------------------------------------------------
END MODULE YOMVWRK
