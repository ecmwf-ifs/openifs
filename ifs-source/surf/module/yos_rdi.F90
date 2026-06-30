MODULE YOS_RDI
  
USE PARKIND1  ,ONLY : JPIM     ,JPRB

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

IMPLICIT NONE

! Maximum number of spectral intervals to store longwave emissivity
! for different surface types
INTEGER(KIND=JPIM), PARAMETER :: NMAXLWEMISS = 16

SAVE

!       ----------------------------------------------------------------
!*    ** *YOERDI* - COEFFICIENTS WITHIN RADIATION INTERFACE
!       ----------------------------------------------------------------

TYPE :: TRDI
REAL(KIND=JPRB) :: RALBSFO   ! DEFAULT SNOW ALBEDO IN THE PRESENCE OF FOREST
! Snow albedo for 0-20 high vegetation types and 2 spectral intervals (1=UV/Vis, 2=Near-IR)
REAL(KIND=JPRB) :: RALB_SNOW_FOREST(0:20,2)
REAL(KIND=JPRB) :: RALBSEAD  ! OPEN SEA ALBEDO FOR DIFFUSE RADIATION
REAL(KIND=JPRB) :: REPALB    ! SECURITY TO AVOID ZERO ALBEDOS.
REAL(KIND=JPRB) :: REMISS_DESERT(NMAXLWEMISS) ! LONGWAVE EMISSIVITY OF DESERT LAND SURFACE
REAL(KIND=JPRB) :: REMISS_LAND  (NMAXLWEMISS) ! LONGWAVE EMISSIVITY OF LAND
REAL(KIND=JPRB) :: REMISS_SNOW  (NMAXLWEMISS) ! LONGWAVE EMISSIVITY OF SNOW
REAL(KIND=JPRB) :: REMISS_SEA   (NMAXLWEMISS) ! LONGWAVE EMISSIVITY OF SEA
! APPROX WEIGHT OF EMISSIVITY SPECTRAL INTERVALS FOR BROADBAND EMISSIVITY DIAGNOSTIC
REAL(KIND=JPRB) :: REMISS_WEIGHT(NMAXLWEMISS)
! APPROX WEIGHT OF EMISSIVITY SPECTRAL INTERVALS (1) OUTSIDE INFRARED
! WINDOW AND (2) INSIDE INFRARED WINDOW (800-1250 CM-1)
REAL(KIND=JPRB) :: REMISS_OLD_WEIGHT(NMAXLWEMISS,2)

END TYPE TRDI

END MODULE YOS_RDI
