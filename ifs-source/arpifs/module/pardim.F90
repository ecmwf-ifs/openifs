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

MODULE PARDIM

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!     JPMXLE : MAXIMUM NUMBER OF LEVELS
!     JPMXGL : MAXIMUM NUMBER OF GAUSSIAN LATITUDES 
!     JPSLWIDE: maximum allowed number of rows for interpolation halos.
!     JPFPPYX : Maximum number of fields in catalogue cilipdy
!     JPNPPM  : Number of interpolation methods in post-processing

INTEGER(KIND=JPIM), PARAMETER :: JPMXLE=200
INTEGER(KIND=JPIM), PARAMETER :: JPMXGL=5120
INTEGER(KIND=JPIM), PARAMETER :: JPSLWIDE=48
INTEGER(KIND=JPIM), PARAMETER :: JPFPPYX=17
INTEGER(KIND=JPIM), PARAMETER :: JPNULNAM=4
INTEGER(KIND=JPIM), PARAMETER :: JPNPPM=4

!     ------------------------------------------------------------------
END MODULE PARDIM
