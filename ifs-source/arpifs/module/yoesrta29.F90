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

MODULE YOESRTA29

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA29* - SRTM COEFFICIENTS FOR INTERVAL 29
!     BAND 29:  820-2600 cm-1 (low - H2O; high - CO2)
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPG = 16, NG29 = 16

REAL(KIND=JPRB) :: KA(5,13,JPG)   
REAL(KIND=JPRB) :: KB(5,13:59,JPG)
REAL(KIND=JPRD) :: KA_D(5,13,JPG)   
REAL(KIND=JPRD) :: KB_D(5,13:59,JPG)
REAL(KIND=JPRB) :: SELFREF(10,JPG),FORREF(4,JPG)
REAL(KIND=JPRB) :: SFLUXREF(JPG)  ,ABSH2O(JPG)  , ABSCO2(JPG)
REAL(KIND=JPRB) :: RAYL
INTEGER(KIND=JPIM) :: LAYREFFR

REAL(KIND=JPRB) :: KAC(5,13,NG29)   ,ABSA(65,NG29)
REAL(KIND=JPRB) :: KBC(5,13:59,NG29),ABSB(235,NG29)
REAL(KIND=JPRB) :: SELFREFC(10,NG29),FORREFC(4,NG29)
REAL(KIND=JPRB) :: SFLUXREFC(NG29)  ,ABSH2OC(NG29)  , ABSCO2C(NG29)

!EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
EQUIVALENCE (KAC(1,1,1),ABSA(1,1)), (KBC(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL     absorption coefficient of major absorber
! KB      : REAL     absorption coefficient of secondary absorber
! SELFREF : REAL     self brodening coefficient for water vapour
! FORREF  : REAL     foreign broadening coefficient for water vapour
! SFLUXREF: REAL     Incident solar radiation in the spectral interval
! ABSH2O  : REAL     line absorption coefficient for H2O
! ABSCO2  : REAL     line absorption coefficient for CO2
! RAYL    : REAL     Rayleigh scattering parameter
! LAYREFFR: INTEGER  reference level for the transition
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! ABSH2OC : REAL     Reduced g-point array for ABSH2O
! ABSCO2C : REAL     Reduced g-point array for ABSCO2
!     -----------------------------------------------------------------
END MODULE YOESRTA29

