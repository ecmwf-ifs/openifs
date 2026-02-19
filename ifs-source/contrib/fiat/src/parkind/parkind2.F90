! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PARKIND2
!
!     *** Define huge kinds for strong typing ***
!
IMPLICIT NONE
SAVE
!
!     Integer Kinds
!     -------------
!
INTEGER, PARAMETER :: JPIH = SELECTED_INT_KIND(18)
!
!     Real Kinds
!     ----------
!
#ifdef REALHUGE
INTEGER, PARAMETER :: JPRH = SELECTED_REAL_KIND(31,291)
#else
INTEGER, PARAMETER :: JPRH = SELECTED_REAL_KIND(13,300)
#endif
!
END MODULE PARKIND2
