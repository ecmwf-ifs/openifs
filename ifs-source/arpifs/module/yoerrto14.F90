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

MODULE YOERRTO14

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO14* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 14
!     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)
!     ABozzo 201306 updated to rrtmg v4.85
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO14 = 16

REAL(KIND=JPRB) , DIMENSION(NO14) :: FRACREFAO
REAL(KIND=JPRB) , DIMENSION(NO14) :: FRACREFBO

REAL(KIND=JPRB) :: KAO(5,13,NO14)
REAL(KIND=JPRB) :: KBO(5,13:59,NO14)
REAL(KIND=JPRD) :: KAO_D(5,13,NO14)
REAL(KIND=JPRD) :: KBO_D(5,13:59,NO14)
REAL(KIND=JPRB) :: SELFREFO(10,NO14)
REAL(KIND=JPRB) :: FORREFO(4,NO14)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO14
