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

MODULE YOESRTA24

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA24* - SRTM COEFFICIENTS FOR INTERVAL 24
!     BAND 24: 12850-16000 cm-1 (low - H2O,O2; high - O2)
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPG = 16, NG24 = 16

REAL(KIND=JPRB) :: KA(9,5,13,JPG) 
REAL(KIND=JPRB) :: KB(5,13:59,JPG)
REAL(KIND=JPRD) :: KA_D(9,5,13,JPG) 
REAL(KIND=JPRD) :: KB_D(5,13:59,JPG)
REAL(KIND=JPRB) :: SELFREF(10,JPG),FORREF(3,JPG)
REAL(KIND=JPRB) :: SFLUXREF(JPG,9)
REAL(KIND=JPRB) :: ABSO3A(JPG), ABSO3B(JPG), RAYLA(JPG,9), RAYLB(JPG)
REAL(KIND=JPRB) :: STRRAT
INTEGER(KIND=JPIM) :: LAYREFFR

REAL(KIND=JPRB) :: KAC(9,5,13,NG24) ,ABSA(585,NG24)
REAL(KIND=JPRB) :: KBC(5,13:59,NG24),ABSB(235,NG24)
REAL(KIND=JPRB) :: SELFREFC(10,NG24),FORREFC(3,NG24)
REAL(KIND=JPRB) :: SFLUXREFC(NG24,9)
REAL(KIND=JPRB) :: ABSO3AC(NG24), ABSO3BC(NG24), RAYLAC(NG24,9), RAYLBC(NG24)

!EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
EQUIVALENCE (KAC(1,1,1,1),ABSA(1,1)), (KBC(1,13,1),ABSB(1,1))

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
! ABSO3A  : REAL     O3 absorption coefficient in first part of band
! ABSO3B  : REAL     O3 absorption coefficient in second part of band
! RAYLA   : REAL     Rayleigh scattering parameter in first part of band
! RAYLB   : REAL     Rayleigh scattering parameter in second part of band   
! STRRAT  : REAL     weighting factor for the transition between tropospheric 
!                    and stratospheric computations
! LAYREFFR: INTEGER  reference level for the transition
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! ABSO3AC : REAL     Reduced g-point array for ABSO3A
! ABSO3BC : REAL     Reduced g-point array for ABSO3B
! RAYLAC  : REAL     Reduced g-point array for RAYLA
! RAYLBC  : REAL     Reduced g-point array for RAYLB
!     -----------------------------------------------------------------
END MODULE YOESRTA24

