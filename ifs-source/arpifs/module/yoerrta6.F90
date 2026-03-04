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

MODULE YOERRTA6

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA6* - RRTM COEFFICIENTS FOR INTERVAL 6
!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)
!     ABozzo 201306 updaten to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG6  = 8

REAL(KIND=JPRB) , DIMENSION(NG6) :: FRACREFA

REAL(KIND=JPRB) , DIMENSION(NG6) :: CFC11ADJ
REAL(KIND=JPRB) , DIMENSION(NG6) :: CFC12


REAL(KIND=JPRB) :: KA(5,13,NG6),ABSA(65,NG6)
REAL(KIND=JPRB) :: SELFREF(10,NG6)
REAL(KIND=JPRB) :: KA_MCO2(19,NG6)
REAL(KIND=JPRB) :: FORREF(4,NG6)

EQUIVALENCE (KA(1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSCO2  : REAL     absorption coefficient for CO2
! ABSA    : REAL     absorption coefficient of major absorber for M reference tropospheric 
!                    pressures and N reference tropospheric temperatures 
! CFC11ADJ: REAL     absorption coefficient for CFC-11 (adjusted)
! CFC12   : REAL     absorption coefficient for CFC-12
! FRACREFA: REAL     distance from r and T reference tabulated points (troposphere)
! KA      : REAL     absorption coefficient of major absorber (equiv. to ABSA)   
! SELFREF : REAL     self broadening coefficient for water vapour
!     -----------------------------------------------------------------
END MODULE YOERRTA6
