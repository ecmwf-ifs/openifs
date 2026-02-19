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

MODULE YOERRTRWT

USE PARKIND1, ONLY : JPRB
!USE PARRRTM,  ONLY : JPGPT, JPG, JPBAND
USE PARRRTM,  ONLY : JPG, JPBAND, JPGMAX

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

!    -------------------------------------------------------------------

REAL(KIND=JPRB) :: FREFA  (JPGMAX,13)
REAL(KIND=JPRB) :: FREFB  (JPGMAX,6)
REAL(KIND=JPRB) :: FREFADF(JPGMAX,13)
REAL(KIND=JPRB) :: FREFBDF(JPGMAX,6)
REAL(KIND=JPRB) :: RWGT   (JPG*JPBAND)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! FREFA  :  REAL   :
! FREFB  :  REAL   :
! FREFADF:  REAL   :
! FREFBDF:  REAL   :
! RWT    :  REAL   : 
!    -------------------------------------------------------------------
END MODULE YOERRTRWT
