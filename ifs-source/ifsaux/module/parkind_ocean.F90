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

MODULE PARKIND_OCEAN
!
!     *** Define usual kinds ***
!
  
IMPLICIT NONE
SAVE
!
!     Integer Kinds
!     -------------
!


!
!     Real Kinds
!     ----------
!
#ifdef PARKIND1_SINGLE_NEMO
INTEGER, PARAMETER :: JPRO = SELECTED_REAL_KIND(6,37)
#else
INTEGER, PARAMETER :: JPRO = SELECTED_REAL_KIND(13,300)
#endif
!

END MODULE PARKIND_OCEAN
