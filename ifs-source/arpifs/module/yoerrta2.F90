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

MODULE YOERRTA2

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA2* - RRTM COEFFICIENTS FOR INTERVAL 2
!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)
! ABozzo May 2013 updated to the last rrtmg
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG2  = 12

!     The ith set of reference fractions are from the ith reference
!     pressure level.
REAL(KIND=JPRB) :: FRACREFA(NG2), FRACREFB(NG2)
REAL(KIND=JPRB) :: KA(5,13,NG2)   , ABSA(65,NG2)
REAL(KIND=JPRB) :: KB(5,13:59,NG2), ABSB(235,NG2)
REAL(KIND=JPRB) :: SELFREF(10,NG2), FORREF(4,NG2)

EQUIVALENCE (KA(1,1,1),ABSA(1,1)),(KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL     absorption coefficient of major absorber for M reference tropospheric 
!                    pressures and N reference tropospheric temperatures 
! ABSB    : REAL     absorption coefficient of secondary absorber for M reference stratospheric
!                    pressures and N reference stratospheric temperatures 
! FRACREFA: REAL     distance from r and T reference tabulated points (troposphere)
! FRACREFB: REAL     distance from r and T reference tabulated points (stratosphere)
! FORREF  : REAL     foreign broadening coefficient for water vapour
! KA      : REAL     absorption coefficient of major absorber (equiv. to ABSA)   
! KB      : REAL     absorption coefficient of secondary absorber (equiv. to ABSB)   
! REFPARAM: REAL     reference water vapour mixing ratio for use in Planck function interpolation
! SELFREF : REAL     self broadening coefficient for water vapour
!     -----------------------------------------------------------------
END MODULE YOERRTA2

