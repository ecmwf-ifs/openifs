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

MODULE YOESRTA16

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA16* - SRTM COEFFICIENTS FOR INTERVAL 16
!     BAND 16:  2600-3250 cm-1 (low - H2O,CH4; high - CH4)
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPG=16, NG16 = 16, NGS15 = 0

REAL(KIND=JPRB) :: KA(9,5,13,JPG) 
REAL(KIND=JPRB) :: KB(5,13:59,JPG)
REAL(KIND=JPRD) :: KA_D(9,5,13,JPG) 
REAL(KIND=JPRD) :: KB_D(5,13:59,JPG)
REAL(KIND=JPRB) :: SELFREF(10,JPG),FORREF(3,JPG)
REAL(KIND=JPRB) :: SFLUXREF(JPG)
REAL(KIND=JPRB) :: RAYL            ,STRRAT1
INTEGER(KIND=JPIM) :: LAYREFFR

REAL(KIND=JPRB) :: KAC(9,5,13,NG16),ABSA(585,NG16)
REAL(KIND=JPRB) :: KBC(5,13:59,NG16),ABSB(235,NG16)
REAL(KIND=JPRB) :: SELFREFC(10,NG16),FORREFC(3,NG16)
REAL(KIND=JPRB) :: SFLUXREFC(NG16)

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
! RAYL    : REAL     Rayleigh scattering parameter
! STRRAT1 : REAL     weighting factor for the transition between tropospheric 
!                    and stratospheric computations
! LAYREFFR: INTEGER  reference level for the transition
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
!     -----------------------------------------------------------------
END MODULE YOESRTA16

