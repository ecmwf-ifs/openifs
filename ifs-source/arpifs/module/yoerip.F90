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
MODULE YOERIP

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Real time related variables for ECMWF radiation: updated UPDTIER and UPDTIM

!     RCODECM: COSINE OF THE DECLINATION
!     RSIDECM:   SINE OF THE DECLINATION

!     RCOVSRM: COSINE OF TRUE SOLAR TIME
!     RSIVSRM:   SINE OF TRUE SOLAR TIME

!     RSOVRM:  Solar time in seconds at centre of radiation timestep
!              (GMT + equation of time)
!     RWSOVRM: As RSOVRM but in radians

TYPE :: TERIP
REAL(KIND=JPRB) :: RCODECM
REAL(KIND=JPRB) :: RSIDECM
REAL(KIND=JPRB) :: RCOVSRM
REAL(KIND=JPRB) :: RSIVSRM
REAL(KIND=JPRB) :: RSOVRM
REAL(KIND=JPRB) :: RWSOVRM


END TYPE TERIP

!!TYPE(TERIP), POINTER :: YRERIP => NULL()

!     ------------------------------------------------------------------
END MODULE YOERIP
