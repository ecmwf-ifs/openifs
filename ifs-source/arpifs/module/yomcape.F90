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

MODULE YOMCAPE

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!*
!     ------------------------------------------------------------------

!     VARIABLES TO CONTROL CAPE COMPUTATION IN FULLPOS:

!        NCAPEITER : NUMBER OF ITERATIONS IN THE NEWTON LOOPS. Default
!                    is the same as NBITER in YOMPHY.

!        NETAPES : NUMBER OF INTERMEDIATE-LAYERS USED FOR CALCULATION OF
!                  VERTICAL ASCENT BETWEEN TWO MODEL PRESSURE LEVELS.
!                  Default value is 2.

!        GCAPEPSD : RATIO OF p/ps OF LAYER ABOVE THE GROUND
!                   IN CASE OF TYPE 2: IN WHICH MOST UNSTABLE PARCEL IS 
!                                       SEARCHED FOR (CAPE Pressure Search Depth)
!                   IN CASE OF TYPE 6: DEFINES DEPTH OF LAYER FOR MIXING THE PARCEL TO BE LIFTED
!                                      (ML CAPE computation)
!
!        NCAPEPSD : DEPTH of GCAPEPSD recalculated homogeneously

!        GCAPERET: FRACTION OF THE CONDENSATES WHICH IS RETAINED, 
!                  I.E. WHICH DOES NOT PRECIPITATE.
!            IF GCAPERET=1. ==> REVERSIBLE MOIST ASCENT.
!                      IT IS ASSUMED THAT ALL THE PARCEL'S CONDENSED
!                      WATER IS RETAINED, THUS CLOUD CONDENSATES 
!                      REDUCE THE BUOYANCY.
!            IF GCAPERET=0. ==> "IRREVERSIBLE" (PSEUDO-ADIABATIC) MOIST ASCENT.
!                       CLOUD CONDENSATES PRECIPITATE  INSTANTANEOUSLY
!                       AND THUS DO NOT AFFECT THE BUOYANCY.
!            GCAPERET CAN BE USED WITH VALUES BETWEEN 0. AND 1..
!            Default value is 0.
!-------------------------------------------------

INTEGER(KIND=JPIM) :: NCAPEITER
INTEGER(KIND=JPIM) :: NETAPES
INTEGER(KIND=JPIM) :: NCAPEPSD

REAL(KIND=JPRB) :: GCAPERET
REAL(KIND=JPRB) :: GCAPEPSD

!     ------------------------------------------------------------------
END MODULE YOMCAPE
