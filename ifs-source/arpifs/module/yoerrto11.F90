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

MODULE YOERRTO11

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO11* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 11
!     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)
!     ABozzo 201306 updated to rrtmg v4.85
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO11 = 16

REAL(KIND=JPRB) , DIMENSION(NO11) :: FRACREFAO
REAL(KIND=JPRB) , DIMENSION(NO11) :: FRACREFBO

REAL(KIND=JPRB) :: KAO(5,13,NO11)
REAL(KIND=JPRB) :: KBO(5,13:59,NO11)
REAL(KIND=JPRD) :: KAO_D(5,13,NO11)
REAL(KIND=JPRD) :: KBO_D(5,13:59,NO11)
REAL(KIND=JPRB) :: KAO_MO2(19,NO11)
REAL(KIND=JPRB) :: KBO_MO2(19,NO11)
REAL(KIND=JPRB) :: SELFREFO(10,NO11)
REAL(KIND=JPRB) :: FORREFO(4,NO11)

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
END MODULE YOERRTO11
