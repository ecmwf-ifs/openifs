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

MODULE YOERRTO3

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO3* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 3
!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)
!      ABozzo 200130517 updated to rrtmg_lw_v4.85:
!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
!                           (high key - h2o,co2; high minor - n2o)
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO3  = 16

REAL(KIND=JPRB) :: FRACREFAO(NO3,9) ,FRACREFBO(NO3,5)

REAL(KIND=JPRB) :: KAO_MN2O(9,19,NO3), KBO_MN2O(5,19,NO3)
REAL(KIND=JPRB) :: KAO(9,5,13,NO3)
REAL(KIND=JPRB) :: KBO(5,5,13:59,NO3)
REAL(KIND=JPRD) :: KAO_D(9,5,13,NO3)
REAL(KIND=JPRD) :: KBO_D(5,5,13:59,NO3)
REAL(KIND=JPRB) :: SELFREFO(10,NO3),FORREFO(4,NO3)


!     ------------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSN2OAO: REAL
! ABSN2OBO: REAL
!FRACREFAO: REAL    
!FRACREFBO: REAL
! KAO     : REAL     
! KBO     : REAL     
! SELFREFO: REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO3
