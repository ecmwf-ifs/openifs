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

MODULE YOMANA

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!*
!      YOMANA - ANALYSIS PARAMETERS

!        D. VASILJEVIC   ECMWF     30/5/90

!     NAME    TYPE     MEANING
!     ----    ----     -------

!     NANDAT   I   ANALYSIS DATE, IN FORM YYMMDD
!     NANTIM   I   ANALYSIS TIME, IN FORM HHMMSS
!     NANMIN   I   ANALYSIS TIME, IN MINUTES FROM 00Z
INTEGER(KIND=JPIM) :: NANDAT
INTEGER(KIND=JPIM) :: NANTIM
INTEGER(KIND=JPIM) :: NANMIN

!     ------------------------------------------------------------------

END MODULE YOMANA
