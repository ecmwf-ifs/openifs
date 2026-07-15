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

INTERFACE
SUBROUTINE SEND_OASIS3_SFX(YDGEOMETRY,PTIMEC,PTSTEP)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1, ONLY : JPRB, JPIM
IMPLICIT NONE
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
REAL(KIND=JPRB),INTENT(IN)    :: PTIMEC    ! time of atmospheric model
REAL(KIND=JPRB),INTENT(IN)    :: PTSTEP    ! time-step of atmospheric model
END SUBROUTINE SEND_OASIS3_SFX
END INTERFACE
