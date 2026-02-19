! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

MODULE PARKIND_WAVE
!
!     *** Define usual kinds ***
!
USE PARKIND1, ONLY : JPIM, JPRB, JPRD

IMPLICIT NONE
SAVE
!
!     Integer Kinds
!     -------------
!
INTEGER, PARAMETER :: JWIM = JPIM

!
!     Real Kinds
!     ----------
!
INTEGER, PARAMETER :: JWRB = JPRB
INTEGER, PARAMETER :: JWRU = JPRD 
!

END MODULE PARKIND_WAVE
