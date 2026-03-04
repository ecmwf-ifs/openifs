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

MODULE YOMDIAGVAR
! 
! Diagnostics of 4DVAR for print-out
!
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE



TYPE DIAGVAR
  INTEGER(KIND=JPIM)  :: NOBS,NCOST
  REAL(KIND=JPRB)     :: CONDNUM,JO,JB,JC,JQ,JP,JH,JCVARBC
END TYPE DIAGVAR 

TYPE(DIAGVAR) :: DIAG_4DVAR
END MODULE YOMDIAGVAR
