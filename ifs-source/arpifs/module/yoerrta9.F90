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

MODULE YOERRTA9

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA9* - RRTM COEFFICIENTS FOR INTERVAL 9
!     BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG9  = 12

REAL(KIND=JPRB) :: FRACREFA(NG9,9)

REAL(KIND=JPRB) , DIMENSION(NG9) :: FRACREFB

REAL(KIND=JPRB) :: KA(9,5,13,NG9) ,ABSA(585,NG9)
REAL(KIND=JPRB) :: KB(5,13:59,NG9) ,ABSB(235,NG9)
REAL(KIND=JPRB) :: KA_MN2O(9,19,NG9)
REAL(KIND=JPRB) :: KB_MN2O(19,NG9)
REAL(KIND=JPRB) :: SELFREF(10,NG9)
REAL(KIND=JPRB) :: FORREF(4,NG9)

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL     absorption coefficient of major absorber for M reference tropospheric 
!                    pressures and N reference tropospheric temperatures 
! ABSB    : REAL     absorption coefficient of secondary absorber for M reference stratospheric
!                    pressures and N reference stratospheric temperatures 
! ABSN2O  : REAL     absorption coefficient for N2O
! CH4REF  : REAL     reference profile for CH4
! ETAREF  : REAL     reference eta profile
! FRACREFA: REAL     distance from r and T reference tabulated points (troposphere)
! FRACREFB: REAL     distance from r and T reference tabulated points (stratosphere)
! H2OREF  : REAL     reference profile for H2O
! KA      : REAL     absorption coefficient of major absorber (equiv. to ABSA)   
! KB      : REAL     absorption coefficient of secondary absorber (equiv. to ABSB)   
! N2OREF  : REAL     reference profile for N2O
! SELFREF : REAL     self broadening coefficient for water vapour
! STRRAT  : REAL     weighting factors for the transition between tropospheric 
!                    and stratospheric computations
!     -----------------------------------------------------------------
END MODULE YOERRTA9
