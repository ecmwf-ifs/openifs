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

MODULE YOMVAREPS

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE
SAVE

!     -----------------------------------------------------------------
!*    ** *YOMVAREPS* - CONTROL PARAMETERS FOR VAREPS CONFIGURATION
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NAFVAREPSMAX=47

INTEGER(KIND=JPIM) :: NFCHO_TRUNC_INI
INTEGER(KIND=JPIM) :: NFCLENGTH_INI
INTEGER(KIND=JPIM) :: NST_FCLENGTH_INI
INTEGER(KIND=JPIM) :: NFCHO_TRUNC
INTEGER(KIND=JPIM) :: NST_TRUNC_INI
INTEGER(KIND=JPIM) :: NST_TRUNC
INTEGER(KIND=JPIM) :: NAFVAREPS
INTEGER(KIND=JPIM) :: NAFVAREPSGC(NAFVAREPSMAX)
LOGICAL :: LVAREPS

!     ------------------------------------------------------------------

!*     Variables for controlling VAREPS

!      Original

!      Author: Roberto Buizza - 10 March 2005

!      LVAREPS          : .T. when running with variable resolution
!      NFCHO_TRUNC_INI  : forecast step used to define the ICs (ie NFCHO_TRUNC of previous VAREPS LEG)
!      NFCLENGTH_INI    : length of the forecast of the previous VAREPS LEG
!      NST_FCLENGTH_INI : NFCLENGTH_INI in number of time steps of current LEG
!      NFCHO_TRUNC      : forecast hour when resolution is reduced
!      NST_TRUNC_INI    : NFCHO_TRUNC_INI in number of time steps
!      NST_TRUNC        : NFCHO_TRUNC in number of time steps of current LEG                     
!      NAFVAREPSMAX     : maximum number of accumulated fields to be re-set if VAREPS LEG gt 1
!      NAFVAREPS        : number of accumulated fields to be re-set if VAREPS LEG gt 1
!      NAFVAREPSGC      : grid codes of accumulated fields to be re-set

!      Schematic of a 3-LEG system with overlap periods:

!      LEG1:   !................!..............!
!              ^                ^              ^
!        NFCHO_TRUNC_INI   NFCHO_TRUNC       FCLENGTH
!      LEG2:                    !..............!..................!........!
!                               ^              ^                  ^        ^
!         in hours .. NFCHO_TRUNC_INI     NFCLENGTH_INI    NFCHO_TRUNC  FCLENGTH      
!         in num t-steps ..  (NST=0)   (NST_NFCLENGTH_INI)  (NST_TRUNC)
!      LEG3:                                                      !........!................!
!                                                                 ^        ^                ^
!         in hours ..                                  NFCHO_TRUNC_INI NFCLENGTH_INI      FCLENGTH   
!         in num t-steps ..                               (NST=0)     (NST_NFCLENGTH_INI)  

!      Schematic of a 3-LEG system T399(0-7d)+T255(6d-15d)+T159(14d-30d) with 1-day overlap periods
!      between each legs (ie the T255 starts at 6d and not 7d, and the T159 starts at 14d)

!      LEG1:   !0............ 144 ......... 168!
!      LEG2:              !144(144+0) .. 168(144+24) .... 336(144+192) .. 360(144+216)!
!      LEG3:                                             !336(336+0)   .. 360(336+24) .. 720(336+384)!
!              ^                ^              ^                  ^        ^                ^

!      Then the T-variables should be set:
!      - LEG1 (0-7d):
!        - NFCHO_TRUNC_INI=0
!        - NFCLENGTH_INI  =0  
!        - NFCHO_TRUNC    =144 (=6d)
!        - FCLENGTH       =168 (=7d)
!      - LEG2 (6d-15d):
!        - NFCHO_TRUNC_INI=144 (=0+144=NFCHO_TRUNC_INI(LEG1)+NFCHO_TRUNC(LEG1))
!        - NFCLENGTH_INI  =168 (=0+168=NFCHO_TRUNC_INI(LEG1)+FCLENGTH(LEG1))
!        - NFCHO_TRUNC    =192 (=8d=14d-6d)
!        - FCLENGTH       =216 (=9d=15d-6d)
!      - LEG3 (14d-30d):
!        - NFCHO_TRUNC_INI=336 (=144+192=NFCHO_TRUNC_INI(LEG2)+NFCHO_TRUNC(LEG2))
!        - NFCLENGTH_INI  =360 (=144+216=NFCHO_TRUNC_INI(LEG2)+FCLENGTH(LEG2))
!        - NFCHO_TRUNC    =0   (since there is no other leg)
!        - FCLENGTH       =384 (=16d)


END MODULE YOMVAREPS
