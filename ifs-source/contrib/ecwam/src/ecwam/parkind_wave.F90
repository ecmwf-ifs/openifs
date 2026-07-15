! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PARKIND_WAVE
!
!     *** Define usual kinds ***
!
USE PARKIND1, ONLY : JPIM, JPRB, JPRD

IMPLICIT NONE
SAVE
PRIVATE
!
!     Integer Kinds
!     -------------
!
INTEGER, PARAMETER, PUBLIC :: JWIM = JPIM

!
!     Real Kinds
!     ----------
!
INTEGER, PARAMETER, PUBLIC :: JWRB = JPRB
INTEGER, PARAMETER, PUBLIC :: JWRU = JPRD
#ifdef PARKIND1_SINGLE_NEMO
INTEGER, PARAMETER, PUBLIC :: JWRO = SELECTED_REAL_KIND(6,37)
#else
INTEGER, PARAMETER, PUBLIC :: JWRO = SELECTED_REAL_KIND(13,300)
#endif
!

END MODULE PARKIND_WAVE
