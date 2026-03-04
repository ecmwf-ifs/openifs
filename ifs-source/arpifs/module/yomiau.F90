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

MODULE YOMIAU

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

! Logical to activate IAU
LOGICAL :: LIAU

! Integer
INTEGER(KIND=JPIM) :: NSTARTIAU ! fist time step of the IAU
INTEGER(KIND=JPIM) :: NSTOPIAU  ! last time step of the IAU

! Real
REAL(KIND=JPRB) :: TSTARTIAU ! fist time (sec.) of the IAU
REAL(KIND=JPRB) :: TSTOPIAU  ! last time (sec.) of the IAU
REAL(KIND=JPRB) :: ALPHAIAU  ! part of the total increment added during the execution
                             ! then increment added at each time step is ALPHAIAU*(xa-xb)/(NSTARTIAU-NSTOPIAU) 

END MODULE YOMIAU
