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


MODULE PARERR

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!*    PARERR - OBS. ERRORS VARIOUS PARAMETERS

!        D. VASILJEVIC    ECMWF    19/09/94

!     NAME                    MEANING
!     ----                    -------
!     JPMXESLV    MAX. NO. OF ERROR STANDARD LEVELS
!     JPMXESLY    MAX. NO. OF ERROR STANDARD LAYERS
!     JPMXEARE    MAX. NO. OF ERROR AREAS / RADIOSONDE ERROR PROFILES
!     JPMXRSTP    MAX. NO. OF RADIOSONDE TYPES

INTEGER(KIND=JPIM), PARAMETER :: JPMXESLV=15
INTEGER(KIND=JPIM), PARAMETER :: JPMXESLY=14
INTEGER(KIND=JPIM), PARAMETER :: JPMXEARE=6
INTEGER(KIND=JPIM), PARAMETER :: JPMXRSTP=200

!-----------------------------------------------------------------------

END MODULE PARERR
