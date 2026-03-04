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

MODULE YOERRTO4

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO4* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 4
!     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)
!     201306 ABozzo updated to rrtmg v4.85
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO4  = 16

REAL(KIND=JPRB) :: FRACREFAO(NO4,9)  ,FRACREFBO(NO4,5) 
REAL(KIND=JPRB) :: KAO(9,5,13,NO4)
REAL(KIND=JPRB) :: KBO(5,5,13:59,NO4)
REAL(KIND=JPRD) :: KAO_D(9,5,13,NO4)
REAL(KIND=JPRD) :: KBO_D(5,5,13:59,NO4)
REAL(KIND=JPRB) :: SELFREFO(10,NO4),FORREFO(4,NO4)


!     ------------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO4
