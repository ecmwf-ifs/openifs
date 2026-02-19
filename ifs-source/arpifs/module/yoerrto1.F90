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

MODULE YOERRTO1

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO1* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 1
!     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
!     ABozzo may 2013 update to rrtmg v4.85
!     band 1:  10-350 cm-1
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO1  = 16

REAL(KIND=JPRB) :: FRACREFAO(NO1)  , FRACREFBO(NO1)
REAL(KIND=JPRB) :: KAO(5,13,NO1)
REAL(KIND=JPRB) :: KBO(5,13:59,NO1)
REAL(KIND=JPRD) :: KAO_D(5,13,NO1)
REAL(KIND=JPRD) :: KBO_D(5,13:59,NO1)
REAL(KIND=JPRB) :: KAO_MN2(19,NO1) , KBO_MN2(19,NO1)
REAL(KIND=JPRB) :: SELFREFO(10,NO1), FORREFO(4,NO1)


!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
!FRACREFAO: REAL    
!FRACREFBO: REAL
! FORREFO : REAL
! KAO     : REAL     
! KBO     : REAL     
! SELFREFO: REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO1
