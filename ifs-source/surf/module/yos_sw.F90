MODULE YOS_SW

USE PARKIND1  ,ONLY : JPIM     ,JPRB

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

IMPLICIT NONE

SAVE

!       ----------------------------------------------------------------
!*    ** *YOS_SW* - COEFFICIENTS FOR SHORTWAVE RADIATION TRANSFER
!       ----------------------------------------------------------------

TYPE :: TSW
REAL(KIND=JPRB), ALLOCATABLE :: RSUN(:)      ! SOLAR FRACTION IN SPECTRAL INTERVALS
REAL(KIND=JPRB), ALLOCATABLE :: RALBICE_AR(:,:) ! monthly sea-ice albedo 
                                       ! in SW spectral intervals 
!                            for Antarctica
REAL(KIND=JPRB), ALLOCATABLE :: RALBICE_AN(:,:) ! for Arctic
!                     1 - ocean (flat response)
!                     2 - sea-ice (as snow or ice depending on month)
!                     3 - wet skin (flat)
!                     4 - vegetation low and snow-free (BR, 1986)
!                     5 - snow on low vegetation       (Warren, 1982)
!                     6 - vegetation high and snow-free (BR, 1986)
!                     7 - snow under high vegetation   (Warren, 1982)
!                     8 - bare soil (Briegleb, Ramanathan, 1986)
!     -----------------------------------------------------------------
END TYPE TSW

END MODULE YOS_SW
