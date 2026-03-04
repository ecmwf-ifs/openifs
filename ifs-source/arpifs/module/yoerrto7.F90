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

MODULE YOERRTO7

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO7* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 7
!     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)
!     ABozzo updated to rrtmg v4.85
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO7  = 16

REAL(KIND=JPRB) :: FRACREFAO(NO7,9)

REAL(KIND=JPRB) , DIMENSION(NO7) :: FRACREFBO
REAL(KIND=JPRB) :: KAO(9,5,13,NO7)
REAL(KIND=JPRB) :: KBO(5,13:59,NO7)
REAL(KIND=JPRD) :: KAO_D(9,5,13,NO7)
REAL(KIND=JPRD) :: KBO_D(5,13:59,NO7)
REAL(KIND=JPRB) :: SELFREFO(10,NO7)
REAL(KIND=JPRB) :: KAO_MCO2(9,19,NO7)
REAL(KIND=JPRB) :: KBO_MCO2(19,NO7)
REAL(KIND=JPRB) :: FORREFO(4,NO7)



!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL  
!     -----------------------------------------------------------------
END MODULE YOERRTO7
