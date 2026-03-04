! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE EC_PARKIND
!
!     *** Define usual kinds for strong typing ***
!
IMPLICIT NONE
SAVE
!
!     Integer Kinds
!     -------------
!
INTEGER, PARAMETER :: JPIM = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: JPIB = SELECTED_INT_KIND(12)

!
!     Real Kinds
!     ----------
!
INTEGER, PARAMETER :: JPRM = SELECTED_REAL_KIND(6,37)
INTEGER, PARAMETER :: JPRD = SELECTED_REAL_KIND(13,300)


END MODULE EC_PARKIND
