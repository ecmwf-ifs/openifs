! (C) Copyright 1989- ECMWF.
! (C) Copyright 1989- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOERDU

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOERDU* - CONTROL, PARAMETERS AND SECURITY IN RADIATION
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: NUAER
INTEGER(KIND=JPIM) :: NTRAER
INTEGER(KIND=JPIM) :: NIMP
INTEGER(KIND=JPIM) :: NOUT

REAL(KIND=JPRB) :: R10E
REAL(KIND=JPRB) :: REPLOG
REAL(KIND=JPRB) :: REPSC
REAL(KIND=JPRB) :: REPSCA
REAL(KIND=JPRB) :: REPSCO
REAL(KIND=JPRB) :: REPSCQ
REAL(KIND=JPRB) :: REPSCT
REAL(KIND=JPRB) :: REPSCW
REAL(KIND=JPRB) :: DIFF

!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! NUAER  : INTEGER   NUMBER OF ABSORBER AMOUNTS        W OR W/O AEROSOLS
! NTRAER : INTEGER   NUMBER OF TRANSMISSION FUNCTIONS  W OR W/O AEROSOLS
! NIMP   : INTEGER   INDEX FOR EXTRA PRINTS WITHIN RADIATION CODE
! NOUT   : INTEGER   UNIT NUMBER FOR THE EXTRA PRINTS
! RCDAY  : REAL
! CCO2   : REAL      CONVERSION COEFFICIENT FOR CO2 IN S.W. CODE
! CH2O   : REAL      CONVERSION COEFFICIENT FOR H2O IN S.W. CODE
! R10E   : REAL      DECIMAL/NATURAL LOG.FACTOR
! DIFF   : REAL      DIFFUSIVITY FACTOR
!-SECURITY THRESHOLDS
! REPLOG : REAL      SEC. EPSILON FOR ABS.AMOUNT IN LAPLACE TRANSFORM
! REPSC  : REAL      SEC. EPSILON FOR CLOUD COVER
! REPSCO : REAL      SEC. EPSILON FOR OZONE AMOUNT
! REPSCQ : REAL      SEC. EPSILON FOR WATER VAPOR
! REPSCT : REAL      SEC. EPSILON FOR SHORTWAVE OPTICAL THICKNESS
! REPSCW : REAL      SEC. EPSILON FOR CLOUD LIQUID WATER PATH

!     -----------------------------------------------------------------
END MODULE YOERDU
