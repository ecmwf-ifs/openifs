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

MODULE YOERRTA8

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA8* - RRTM COEFFICIENTS FOR INTERVAL 8
!     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG8  = 8

REAL(KIND=JPRB) , DIMENSION(NG8) :: FRACREFA
REAL(KIND=JPRB) , DIMENSION(NG8) :: FRACREFB
REAL(KIND=JPRB) , DIMENSION(NG8) :: CFC12
REAL(KIND=JPRB) , DIMENSION(NG8) :: CFC22ADJ

REAL(KIND=JPRB) :: KA(5,13,NG8)    ,ABSA(65,NG8)
REAL(KIND=JPRB) :: KB(5,13:59,NG8) ,ABSB(235,NG8)
REAL(KIND=JPRB) :: KA_MCO2(19,NG8)
REAL(KIND=JPRB) :: KA_MN2O(19,NG8)
REAL(KIND=JPRB) :: KA_MO3(19,NG8)
REAL(KIND=JPRB) :: KB_MCO2(19,NG8)
REAL(KIND=JPRB) :: KB_MN2O(19,NG8)
REAL(KIND=JPRB) :: SELFREF(10,NG8)
REAL(KIND=JPRB) :: FORREF(4,NG8)


EQUIVALENCE (KA(1,1,1),ABSA(1,1)),(KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL     absorption coefficient of major absorber for M reference tropospheric 
!                    pressures and N reference tropospheric temperatures 
! ABSB    : REAL     absorption coefficient of secondary absorber for M reference stratospheric
!                    pressures and N reference stratospheric temperatures 
! ABSCO2A : REAL     absorption coefficient for CO2 for M reference tropospheric 
!                    pressures and N reference tropospheric temperatures 
! ABSCO2B : REAL     absorption coefficient for CO2 for M reference stratospheric 
!                    pressures and N reference stratospheric temperatures     
! ABSN2OA : REAL     absorption coefficient for N2O for M reference tropospheric 
!                    pressures and N reference tropospheric temperatures     
! ABSN2OB : REAL     absorption coefficient for N2O for M reference stratospheric 
!                    pressures and N reference stratospheric temperatures 
! CFC12   : REAL     absorption coefficient for CFC-12
! CFC22ADJ: REAL     absorption coefficient for CFC-22 (adjusted)
! FRACREFA: REAL     distance from r and T reference tabulated points (troposphere)
! FRACREFB: REAL     distance from r and T reference tabulated points (stratosphere)
! H2OREF  : REAL     reference profile for H2O
! KA      : REAL     absorption coefficient of major absorber (equiv. to ABSA)   
! KB      : REAL     absorption coefficient of secondary absorber (equiv. to ABSB)   
! N2OREF  : REAL     reference profile for N2O
! O3REF   : REAL     reference profile for O3
! SELFREF : REAL     self broadening coefficient for water vapour
!     -----------------------------------------------------------------
END MODULE YOERRTA8
