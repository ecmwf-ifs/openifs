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

MODULE YOERRTA12

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA12* - RRTM COEFFICIENTS FOR INTERVAL 12
!     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG12 = 8

REAL(KIND=JPRB) :: FRACREFA(NG12,9)
REAL(KIND=JPRB) :: KA(9,5,13,NG12) ,ABSA(585,NG12)
REAL(KIND=JPRB) :: SELFREF(10,NG12)
REAL(KIND=JPRB) :: FORREF(4,NG12)

REAL(KIND=JPRB) :: STRRAT

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL     absorption coefficient of major absorber for M reference tropospheric 
!                    pressures and N reference tropospheric temperatures 
! FRACREFA: REAL     distance from r and T reference tabulated points (troposphere)
! KA      : REAL     absorption coefficient of major absorber (equiv. to ABSA)   
! SELFREF : REAL     self broadening coefficient for water vapour
! STRRAT  : REAL     weighting factors for the transition between tropospheric 
!                    and stratospheric computations
!     -----------------------------------------------------------------
END MODULE YOERRTA12
