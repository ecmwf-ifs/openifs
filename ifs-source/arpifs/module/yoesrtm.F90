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

MODULE YOESRTM

USE PARKIND1  ,ONLY : JPIM
USE PARSRTM   ,ONLY : JPGMAX

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!     Parameters relevant to AER's RRTM-SW radiation scheme: Part 2

!     20110610 JJMorcrette

!     Modified to allow possibilities of different g-point numbers.  
!     ------------------------------------------------------------------

!INTEGER(KIND=JPIM) :: JPGPT
!INTEGER(KIND=JPIM) :: JPGSW

!INTEGER(KIND=JPIM) :: NG16
!INTEGER(KIND=JPIM) :: NG17
!INTEGER(KIND=JPIM) :: NG18
!INTEGER(KIND=JPIM) :: NG19
!INTEGER(KIND=JPIM) :: NG20
!INTEGER(KIND=JPIM) :: NG21
!INTEGER(KIND=JPIM) :: NG22
!INTEGER(KIND=JPIM) :: NG23
!INTEGER(KIND=JPIM) :: NG24
!INTEGER(KIND=JPIM) :: NG25
!INTEGER(KIND=JPIM) :: NG26
!INTEGER(KIND=JPIM) :: NG27
!INTEGER(KIND=JPIM) :: NG28
!INTEGER(KIND=JPIM) :: NG29

!-- NGnn : number of g-points in bands nn=16 to 29
!- as used for the (operational) 112 g-points version of RRTM_SW
INTEGER(KIND=JPIM), PARAMETER :: JPGPT  = 112
INTEGER(KIND=JPIM), PARAMETER :: JPGSW  = 112

INTEGER(KIND=JPIM), PARAMETER :: NG16 = 6
INTEGER(KIND=JPIM), PARAMETER :: NG17 = 12
INTEGER(KIND=JPIM), PARAMETER :: NG18 = 8
INTEGER(KIND=JPIM), PARAMETER :: NG19 = 8
INTEGER(KIND=JPIM), PARAMETER :: NG20 = 10
INTEGER(KIND=JPIM), PARAMETER :: NG21 = 10
INTEGER(KIND=JPIM), PARAMETER :: NG22 = 2
INTEGER(KIND=JPIM), PARAMETER :: NG23 = 10
INTEGER(KIND=JPIM), PARAMETER :: NG24 = 8
INTEGER(KIND=JPIM), PARAMETER :: NG25 = 6
INTEGER(KIND=JPIM), PARAMETER :: NG26 = 6
INTEGER(KIND=JPIM), PARAMETER :: NG27 = 8
INTEGER(KIND=JPIM), PARAMETER :: NG28 = 6
INTEGER(KIND=JPIM), PARAMETER :: NG29 = 12


INTEGER(KIND=JPIM) :: NGN(JPGMAX), NGBSW(JPGMAX)

!     ------------------------------------------------------------------
END MODULE YOESRTM

