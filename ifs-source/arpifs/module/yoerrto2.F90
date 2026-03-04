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

MODULE YOERRTO2

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO2* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 2
!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)
! ABozzo May 2013 updated to rrtmg v4.85
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO2  = 16

!     The ith set of reference fractions are from the ith reference
!     pressure level.
REAL(KIND=JPRB) :: FRACREFAO(NO2), FRACREFBO(NO2)
REAL(KIND=JPRB) :: KAO(5,13,NO2)
REAL(KIND=JPRB) :: KBO(5,13:59,NO2)
REAL(KIND=JPRD) :: KAO_D(5,13,NO2)
REAL(KIND=JPRD) :: KBO_D(5,13:59,NO2)
REAL(KIND=JPRB) :: SELFREFO(10,NO2) , FORREFO(4,NO2)



!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
!FRACREFAO: REAL    
!FRACREFBO: REAL
! KAO     : REAL     
! KBO     : REAL     
! SELFREFO: REAL
! FORREFO : REAL  
!     -----------------------------------------------------------------
END MODULE YOERRTO2
