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

MODULE YOERRTO6

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO6* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 6
!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)
!     ABozzo 201306 update to rrtmg v4.85
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO6  = 16

REAL(KIND=JPRB) , DIMENSION(NO6) :: FRACREFAO

REAL(KIND=JPRB) , DIMENSION(NO6) :: CFC11ADJO
REAL(KIND=JPRB) , DIMENSION(NO6) :: CFC12O


REAL(KIND=JPRB) :: KAO(5,13,NO6)
REAL(KIND=JPRD) :: KAO_D(5,13,NO6)
REAL(KIND=JPRB) :: SELFREFO(10,NO6)
REAL(KIND=JPRB) :: KAO_MCO2(19,NO6)
REAL(KIND=JPRB) :: FORREFO(4,NO6)



!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO6
