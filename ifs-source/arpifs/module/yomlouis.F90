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

MODULE YOMLOUIS

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

TYPE TLOUIS
!*
!  ---------------------------------------------------------------------
!  This comdeck contains the fit of R in QNSE scheme
!  ---------------------------------------------------------------------

!R at neutrality
REAL(KIND=JPRB) :: RLOUIS_S0=0.136_JPRB
REAL(KIND=JPRB) :: PLOUIS_S0=0.130_JPRB

!fit of R(Ri) for QNSE
!unstable stratification
REAL(KIND=JPRB) :: RLOUIS_GU1=-2.67_JPRB
REAL(KIND=JPRB) :: RLOUIS_GU2=-6.15_JPRB
REAL(KIND=JPRB) :: PLOUIS_GU1=-1.13_JPRB
REAL(KIND=JPRB) :: PLOUIS_GU2=-3.80_JPRB

!stable stratification
REAL(KIND=JPRB) :: RLOUIS_GS1=1.0_JPRB
REAL(KIND=JPRB) :: RLOUIS_GS2=2.08_JPRB
REAL(KIND=JPRB) :: RLOUIS_GS3=1.0_JPRB
REAL(KIND=JPRB) :: RLOUIS_GS4=2.7_JPRB
REAL(KIND=JPRB) :: PLOUIS_GS1=1.0_JPRB
REAL(KIND=JPRB) :: PLOUIS_GS2=1.93_JPRB
REAL(KIND=JPRB) :: PLOUIS_GS3=1.1_JPRB
REAL(KIND=JPRB) :: PLOUIS_GS4=3.55_JPRB

END TYPE TLOUIS

TYPE(TLOUIS), POINTER :: YRLOUIS => NULL()

END MODULE YOMLOUIS
