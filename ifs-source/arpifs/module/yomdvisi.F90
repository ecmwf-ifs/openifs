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

MODULE YOMDVISI

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE TDVISI
!*
!     ------------------------------------------------------------------
!     COEFFICIENTS EXTINCTION FOR VISIBILITY
!     
REAL(KIND=JPRB) :: HVISI
REAL(KIND=JPRB) :: COEF_CM1
REAL(KIND=JPRB) :: COEF_CM2
REAL(KIND=JPRB) :: COEF_CM3
REAL(KIND=JPRB) :: COEF_CM4
REAL(KIND=JPRB) :: COEF_RM1
REAL(KIND=JPRB) :: COEF_RM2
REAL(KIND=JPRB) :: COEF_IM1
REAL(KIND=JPRB) :: COEF_IM2
REAL(KIND=JPRB) :: COEF_SM1
REAL(KIND=JPRB) :: COEF_SM2
REAL(KIND=JPRB) :: COEF_GM1
REAL(KIND=JPRB) :: COEF_GM2

END TYPE TDVISI

!     ------------------------------------------------------------------
END MODULE YOMDVISI
