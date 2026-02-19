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

MODULE YOERRTA11

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA11* - RRTM COEFFICIENTS FOR INTERVAL 11
!     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG11 = 8

REAL(KIND=JPRB) , DIMENSION(NG11) :: FRACREFA
REAL(KIND=JPRB) , DIMENSION(NG11) :: FRACREFB

REAL(KIND=JPRB) :: KA(5,13,NG11)   , ABSA(65,NG11)
REAL(KIND=JPRB) :: KB(5,13:59,NG11), ABSB(235,NG11)
REAL(KIND=JPRB) :: KA_MO2(19,NG11)
REAL(KIND=JPRB) :: KB_MO2(19,NG11)
REAL(KIND=JPRB) :: SELFREF(10,NG11)
REAL(KIND=JPRB) :: FORREF(4,NG11)

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
! FRACREFA: REAL     distance from r and T reference tabulated points (troposphere)
! FRACREFB: REAL     distance from r and T reference tabulated points (stratosphere)
! KA      : REAL     absorption coefficient of major absorber (equiv. to ABSA)   
! KB      : REAL     absorption coefficient of secondary absorber (equiv. to ABSB)   
! SELFREF : REAL     self broadening coefficient for water vapour
!     -----------------------------------------------------------------
END MODULE YOERRTA11
