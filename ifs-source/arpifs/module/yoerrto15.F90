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

MODULE YOERRTO15

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO15* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 15
!     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; high - nothing)
!     ABozzo 2001306 updated to rrtmg v4.85
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!                              (high - nothing)
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO15 = 16

REAL(KIND=JPRB) :: FRACREFAO(NO15,9)

REAL(KIND=JPRB) :: KAO(9,5,13,NO15)
REAL(KIND=JPRD) :: KAO_D(9,5,13,NO15)
REAL(KIND=JPRB) :: SELFREFO(10,NO15)
REAL(KIND=JPRB) :: FORREFO(4,NO15)
REAL(KIND=JPRB) :: KAO_MN2(9,19,NO15)
!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL 
!     -----------------------------------------------------------------
END MODULE YOERRTO15
