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

MODULE YOERRTFTR

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!USE PARRRTM, ONLY : JPBAND, JPG, JPGMAX
!USE YOERRTM, ONLY : JPGPT

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

!    -------------------------------------------------------------------

INTEGER(KIND=JPIM) :: NGC(16)
INTEGER(KIND=JPIM) :: NGS(16)

INTEGER(KIND=JPIM) :: NGN(256)
INTEGER(KIND=JPIM) :: NGB(256)

INTEGER(KIND=JPIM) :: NGM(256)
REAL(KIND=JPRB) ::    WT(16)

!INTEGER(KIND=JPIM) :: NGC(JPBAND)
!INTEGER(KIND=JPIM) :: NGS(JPBAND)
!
!INTEGER(KIND=JPIM) :: NGN(JPGMAX)
!INTEGER(KIND=JPIM) :: NGB(JPGMAX)
!
!INTEGER(KIND=JPIM) :: NGM(JPG*JPBAND)
!REAL(KIND=JPRB) ::    WT(JPG)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
!  NGC   : INTEGER :
!  NGS   : INTEGER :
!  NGN   : INTEGER :
!  NGB   : INTEGER :
!  NGM   : INTEGER :
!  WT    : REAL    :
!    -------------------------------------------------------------------
END MODULE YOERRTFTR

