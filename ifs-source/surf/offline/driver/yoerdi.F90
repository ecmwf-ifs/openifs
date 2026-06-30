MODULE YOERDI
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERDI* - COEFFICIENTS WITHIN RADIATION INTERFACE
!     -----------------------------------------------------------------

REAL(KIND=JPRB) :: RRAE
REAL(KIND=JPRB) :: RSUNDUR
REAL(KIND=JPRB) :: RCARDI
REAL(KIND=JPRB) :: RCH4
REAL(KIND=JPRB) :: RN2O
REAL(KIND=JPRB) :: RO3
REAL(KIND=JPRB) :: RCFC11
REAL(KIND=JPRB) :: RCFC12
REAL(KIND=JPRB) :: REPCLC
REAL(KIND=JPRB) :: REPH2O
REAL(KIND=JPRB) :: RCCO2, RCCH4, RCN2O, RCCFC11, RCCFC12
REAL(KIND=JPRB) :: RHVAR(151,5), RFVAR(101,3,5), RINCSOL(154), RSOLINC

!        * E.C.M.W.F. PHYSICS PACKAGE *

!     Original  J.-J. MORCRETTE       E.C.M.W.F.      89/07/14
!     Modified  P. Viterbo    99/03/26    Surface tiling
!     Modified  P. Viterbo    24/05/2004  surf library
!     Modified JJMorcrette    2005/01/19  GHG and Solar constant variability

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! RRAE   : EFFECT OF EARTH'S CURVATURE ON COSINE SOLAR ZENITH ANGLE
! RSUNDUR: MINIMUM DIRECT SOLAR FOR COMPUTING SOLAR DURATION
! RCARDI : SPECIFIC ATMOSPHERIC CONTENT IN CO2
! REPCLC : SECURITY TO AVOID ZERO OR ONE CLOUD COVERS
! REPH2O : SECURITY TO AVOID WATER VAPOUR CONTENT IN A LAYER
!          TO BE MORE THAN THE RESPECTIVE VALUE AT SATURATION.
!     -----------------------------------------------------------------
END MODULE YOERDI
