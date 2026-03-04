! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

MODULE YOMGP1C9

!    Grid point array at time t in the surface physics

! - 1D (PROGNOSTIC QUANTITIES) .

! -   Splitting of the array GP9

!    SNS*       : MASS OF SNOW PER UNIT SURFACE.
!    ASN*       : SNOW ALBEDO
!    RSN*       : SNOW DENSITY
!    TSN*       : SNOW TEMPERATURE
!    TSA*       : MULTI-LAYER SOIL TEMPERATURE
!    WSA*       : MULTI-LAYER SOIL WETNESS
!    TIA*       : MULTI-LAYER ICE TEMPERATURE
!    WL*        : SKIN RESERVOIR WATER CONTENT
!    TL*        : SKIN (BRIGHTNESS) TEMPERATURE
!    EXTRP*     : EXTRA MULTI-LAYER FIELDS
!    XTRP2*     : EXTRA FIELDS

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE PARDIM1C , ONLY : JPGPP
!USE YOMPHYDS , ONLY : YRPHYDS

IMPLICIT NONE

SAVE

INTEGER(KIND=JPIM), PARAMETER :: ICSS = 4
INTEGER(KIND=JPIM), PARAMETER :: KLON = 1
REAL(KIND=JPRB) :: GP9(JPGPP)
REAL(KIND=JPRB) :: SNS9(KLON)
REAL(KIND=JPRB) :: ASN9(KLON)
REAL(KIND=JPRB) :: RSN9(KLON)
REAL(KIND=JPRB) :: TSN9(KLON)
REAL(KIND=JPRB) :: TSA9(ICSS)
REAL(KIND=JPRB) :: WSA9(ICSS)
REAL(KIND=JPRB) :: TIA9(ICSS)
REAL(KIND=JPRB) :: WL9(KLON)
REAL(KIND=JPRB) :: TL9(KLON)
REAL(KIND=JPRB) :: EXTRP9(50,20)
REAL(KIND=JPRB) :: XTRP29(50)
!REAL(KIND=JPRB) :: EXTRP9(YRPHYDS%JPCXP,YRPHYDS%JPVXP)
!REAL(KIND=JPRB) :: XTRP29(YRPHYDS%JPVXP2)

!     ------------------------------------------------------------------

END MODULE YOMGP1C9
