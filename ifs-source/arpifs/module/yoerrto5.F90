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

MODULE YOERRTO5

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO5* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 5
!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)
!     ABozzo 201306 updated to rrtmg v4.85
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO5  = 16

REAL(KIND=JPRB) :: FRACREFAO(NO5,9) ,FRACREFBO(NO5,5)

REAL(KIND=JPRB) , DIMENSION(NO5) :: CCL4O

REAL(KIND=JPRB) :: KAO_MO3(9,19,NO5)
REAL(KIND=JPRB) :: KAO(9,5,13,NO5)
REAL(KIND=JPRB) :: KBO(5,5,13:59,NO5)
REAL(KIND=JPRD) :: KAO_D(9,5,13,NO5)
REAL(KIND=JPRD) :: KBO_D(5,5,13:59,NO5)
REAL(KIND=JPRB) :: SELFREFO(10,NO5)
REAL(KIND=JPRB) :: FORREFO(4,NO5)


!     -----------------------------------------------------------------
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
END MODULE YOERRTO5
