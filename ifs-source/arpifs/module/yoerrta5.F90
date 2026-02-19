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

MODULE YOERRTA5

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA5* - RRTM COEFFICIENTS FOR INTERVAL 5
!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG5  = 16

REAL(KIND=JPRB) :: FRACREFA(NG5,9) ,FRACREFB(NG5,5)

REAL(KIND=JPRB) , DIMENSION(NG5) :: CCL4

REAL(KIND=JPRB) :: KA(9,5,13,NG5)   ,ABSA(585,NG5)
REAL(KIND=JPRB) :: KB(5,5,13:59,NG5),ABSB(1175,NG5)
REAL(KIND=JPRB) :: KA_MO3(9,19,NG5)
REAL(KIND=JPRB) :: SELFREF(10,NG5)
REAL(KIND=JPRB) :: FORREF(4,NG5)

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL     absorption coefficient of major absorber for M reference tropospheric 
!                    pressures and N reference tropospheric temperatures 
! ABSB    : REAL     absorption coefficient of secondary absorber for M reference stratospheric
!                    pressures and N reference stratospheric temperatures 
! CCL4    : REAL     absorption coefficient for CCl4
! FRACREFA: REAL     distance from r and T reference tabulated points (troposphere)
! FRACREFB: REAL     distance from r and T reference tabulated points (stratosphere)
! KA      : REAL     absorption coefficient of major absorber (equiv. to ABSA)   
! KB      : REAL     absorption coefficient of secondary absorber (equiv. to ABSB)   
! SELFREF : REAL     self broadening coefficient for water vapour
! STRRATn : REAL     weighting factors for the transition between tropospheric 
!                    and stratospheric computations
!     -----------------------------------------------------------------
END MODULE YOERRTA5
