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

MODULE TOVSCV_BGC_MOD
USE PARKIND1   , ONLY : JPIM, JPRB, JPRD
USE YOMVAR     , ONLY : LTOVSCV, LCLDSINK
USE YOMLUN     , ONLY : NULOUT, NULERR
USE YOMANCS    , ONLY : RMDI
USE YOMHOOK    , ONLY : LHOOK, DR_HOOK, JPHOOK
USE TOVSCV_BASE_MOD
USE TOVSCV_MOD
IMPLICIT NONE
PRIVATE
TYPE,EXTENDS(TOVSCV_BASE),PUBLIC :: TOVSCV_BGC
!     TOVSCVBG: TOVS control variable, background, not used by OOPS
!     TOVSCVER: TOVS control variable, standard deviation of error in bg (sqrt of diagonal of B)
  REAL(KIND=JPRB),ALLOCATABLE :: TOVSCVER(:,:)
  REAL(KIND=JPRB),ALLOCATABLE :: TOVSCVBG(:,:)
END TYPE TOVSCV_BGC


END MODULE TOVSCV_BGC_MOD
