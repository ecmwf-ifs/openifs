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

MODULE PARSRTM

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!     Parameters relevant to AER's RRTM-SW radiation scheme

!     030224  JJMorcrette

!     Modified for g-point reduction from 224 to 112.  
!     Swap code below to restore 224 g-point set. 
!     Mar2004 MJIacono, AER
!     20110322 JJMorcrette : additional comments
!     20110603 JJMorcrette reduced number of g-points
!     ------------------------------------------------------------------

!-- basic spectral information unrelated to number of g-points
! JPG     : INTEGER : maximum number of g-points in a given spectral band
! JPBAND  : INTEGER : total number of spectral bands 
! JPSW    : INTEGER : total number of shortwave spectral bands
! JPB1    : INTEGER : starting index of shortwave spectrum
! JPB2    : INTEGER : end index of shortwave spectrum

INTEGER(KIND=JPIM), PARAMETER :: JPG    = 16
INTEGER(KIND=JPIM), PARAMETER :: JPBAND = 29
INTEGER(KIND=JPIM), PARAMETER :: JPSW   = 14
INTEGER(KIND=JPIM), PARAMETER :: JPB1   = 16
INTEGER(KIND=JPIM), PARAMETER :: JPB2   = 29
INTEGER(KIND=JPIM), PARAMETER :: JPGMAX = 224

!-- other information that could be relevant for RRTM_SW
!-- NB: The following parameters are unused within the ECMWF IFS. 
!       They relate to the description of the optical properties 
!       in the original cloud model embedded in RRTM_SW
!INTEGER(KIND=JPIM), PARAMETER :: JMCMU  = 32
!INTEGER(KIND=JPIM), PARAMETER :: JMUMU  = 32
!INTEGER(KIND=JPIM), PARAMETER :: JMPHI  = 3
!INTEGER(KIND=JPIM), PARAMETER :: JMXANG = 4
!INTEGER(KIND=JPIM), PARAMETER :: JMXSTR = 16

!-- original spectral grid before spectral averaging
!-- original from AER, Inc with 224 g-points
INTEGER(KIND=JPIM), PARAMETER :: NGS16 = 0
INTEGER(KIND=JPIM), PARAMETER :: NGS17 = 16
INTEGER(KIND=JPIM), PARAMETER :: NGS18 = 32
INTEGER(KIND=JPIM), PARAMETER :: NGS19 = 48
INTEGER(KIND=JPIM), PARAMETER :: NGS20 = 64
INTEGER(KIND=JPIM), PARAMETER :: NGS21 = 80
INTEGER(KIND=JPIM), PARAMETER :: NGS22 = 96
INTEGER(KIND=JPIM), PARAMETER :: NGS23 = 112
INTEGER(KIND=JPIM), PARAMETER :: NGS24 = 128
INTEGER(KIND=JPIM), PARAMETER :: NGS25 = 144
INTEGER(KIND=JPIM), PARAMETER :: NGS26 = 160
INTEGER(KIND=JPIM), PARAMETER :: NGS27 = 176
INTEGER(KIND=JPIM), PARAMETER :: NGS28 = 192
INTEGER(KIND=JPIM), PARAMETER :: NGS29 = 208

!-------------------------------------------------------------------------------
!-- NGnn : number of g-points in bands nn=16 to 29
!- as used in the Ng g-points version of RRTM_SW
!-------------------------------------------------------------------------------
!-- configuration with 14 spectral intervals
!   and a total of 56 g-points (14xvariable number)

!INTEGER(KIND=JPIM), PARAMETER :: JPGPT  = 56
!
!INTEGER(KIND=JPIM), PARAMETER :: NG16 = 3
!INTEGER(KIND=JPIM), PARAMETER :: NG17 = 6
!INTEGER(KIND=JPIM), PARAMETER :: NG18 = 4
!INTEGER(KIND=JPIM), PARAMETER :: NG19 = 4
!INTEGER(KIND=JPIM), PARAMETER :: NG20 = 5
!INTEGER(KIND=JPIM), PARAMETER :: NG21 = 5
!INTEGER(KIND=JPIM), PARAMETER :: NG22 = 1
!INTEGER(KIND=JPIM), PARAMETER :: NG23 = 5
!INTEGER(KIND=JPIM), PARAMETER :: NG24 = 4
!INTEGER(KIND=JPIM), PARAMETER :: NG25 = 3
!INTEGER(KIND=JPIM), PARAMETER :: NG26 = 3
!INTEGER(KIND=JPIM), PARAMETER :: NG27 = 4
!INTEGER(KIND=JPIM), PARAMETER :: NG28 = 3
!INTEGER(KIND=JPIM), PARAMETER :: NG29 = 6
!-------------------------------------------------------------------------------
!-- configuration with 14 spectral intervals
!   and a total of 112 g-points (14xvariable number)
!
!INTEGER(KIND=JPIM), PARAMETER :: JPGPT  = 112
!
!INTEGER(KIND=JPIM), PARAMETER :: NG16 = 6
!INTEGER(KIND=JPIM), PARAMETER :: NG17 = 12
!INTEGER(KIND=JPIM), PARAMETER :: NG18 = 8
!INTEGER(KIND=JPIM), PARAMETER :: NG19 = 8
!INTEGER(KIND=JPIM), PARAMETER :: NG20 = 10
!INTEGER(KIND=JPIM), PARAMETER :: NG21 = 10
!INTEGER(KIND=JPIM), PARAMETER :: NG22 = 2
!INTEGER(KIND=JPIM), PARAMETER :: NG23 = 10
!INTEGER(KIND=JPIM), PARAMETER :: NG24 = 8
!INTEGER(KIND=JPIM), PARAMETER :: NG25 = 6
!INTEGER(KIND=JPIM), PARAMETER :: NG26 = 6
!INTEGER(KIND=JPIM), PARAMETER :: NG27 = 8
!INTEGER(KIND=JPIM), PARAMETER :: NG28 = 6
!INTEGER(KIND=JPIM), PARAMETER :: NG29 = 12

!-------------------------------------------------------------------------------
!-- configuration with 14 spectral intervals 
!   and a total of 224 g-points (14x16)
! 
!INTEGER(KIND=JPIM), PARAMETER :: JPGPT  = 224

!INTEGER(KIND=JPIM), PARAMETER :: NG16 = 16
!INTEGER(KIND=JPIM), PARAMETER :: NG17 = 16
!INTEGER(KIND=JPIM), PARAMETER :: NG18 = 16
!INTEGER(KIND=JPIM), PARAMETER :: NG19 = 16
!INTEGER(KIND=JPIM), PARAMETER :: NG20 = 16
!INTEGER(KIND=JPIM), PARAMETER :: NG21 = 16
!INTEGER(KIND=JPIM), PARAMETER :: NG22 = 16
!INTEGER(KIND=JPIM), PARAMETER :: NG23 = 16
!INTEGER(KIND=JPIM), PARAMETER :: NG24 = 16
!INTEGER(KIND=JPIM), PARAMETER :: NG25 = 16
!INTEGER(KIND=JPIM), PARAMETER :: NG26 = 16
!INTEGER(KIND=JPIM), PARAMETER :: NG27 = 16
!INTEGER(KIND=JPIM), PARAMETER :: NG28 = 16
!INTEGER(KIND=JPIM), PARAMETER :: NG29 = 16

!     ------------------------------------------------------------------
END MODULE PARSRTM

