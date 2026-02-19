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

MODULE YOERRTO13

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO13* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 13
!     BAND 13:  2080-2250 cm-1 (low - H2O,N2O; high - nothing)
!     ABozzo 201306 updated to rrtmg v4.85
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO13 = 16

REAL(KIND=JPRB) :: FRACREFAO(NO13,9)
REAL(KIND=JPRB) , DIMENSION(NO13) :: FRACREFBO

REAL(KIND=JPRB) :: KAO(9,5,13,NO13)
REAL(KIND=JPRD) :: KAO_D(9,5,13,NO13)
REAL(KIND=JPRB) :: SELFREFO(10,NO13)
REAL(KIND=JPRB) :: KAO_MCO2(9,19,NO13)
REAL(KIND=JPRB) :: KAO_MCO(9,19,NO13)
REAL(KIND=JPRB) :: KBO_MO3(19,NO13)
REAL(KIND=JPRB) :: FORREFO(4,NO13)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL
!     -----------------------------------------------------------------
END MODULE YOERRTO13
