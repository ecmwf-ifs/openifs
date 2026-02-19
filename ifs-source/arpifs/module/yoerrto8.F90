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

MODULE YOERRTO8

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO8* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 8
!     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)
!     ABozzo 201306 updated to rrtmg v4.85
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO8  = 16

REAL(KIND=JPRB) , DIMENSION(NO8) :: FRACREFAO
REAL(KIND=JPRB) , DIMENSION(NO8) :: FRACREFBO
REAL(KIND=JPRB) , DIMENSION(NO8) :: CFC12O
REAL(KIND=JPRB) , DIMENSION(NO8) :: CFC22ADJO

REAL(KIND=JPRB) :: KAO(5,13,NO8)
REAL(KIND=JPRB) :: KBO(5,13:59,NO8)
REAL(KIND=JPRD) :: KAO_D(5,13,NO8)
REAL(KIND=JPRD) :: KBO_D(5,13:59,NO8)
REAL(KIND=JPRB) :: SELFREFO(10,NO8)
REAL(KIND=JPRB) :: KAO_MCO2(19,NO8)
REAL(KIND=JPRB) :: KAO_MN2O(19,NO8)
REAL(KIND=JPRB) :: KAO_MO3(19,NO8)
REAL(KIND=JPRB) :: KBO_MCO2(19,NO8)
REAL(KIND=JPRB) :: KBO_MN2O(19,NO8)
REAL(KIND=JPRB) :: FORREFO(4,NO8)



!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSCO2A : REAL     
! ABSCO2B : REAL     
! ABSN2OA : REAL     
! ABSN2OB : REAL 
! CFC12   : REAL     
! CFC22ADJ: REAL     
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO8
