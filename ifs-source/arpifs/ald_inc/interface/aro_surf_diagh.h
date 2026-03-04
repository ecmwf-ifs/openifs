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
SUBROUTINE ARO_SURF_DIAGH(YDGEOMETRY,YDMODEL,YDGFL,YDSURF,YDSPEC,YDCFU,YDXFU,YDRIP)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL   , ONLY : MODEL
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMGFL , ONLY : TGFL
USE YOMRIP , ONLY : TRIP
USE YOMCFU , ONLY : TCFU
USE YOMXFU , ONLY : TXFU
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD
TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL)         ,INTENT(INOUT) :: YDMODEL
TYPE(TGFL), INTENT(INOUT) :: YDGFL
TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(SPECTRAL_FIELD), INTENT(INOUT) :: YDSPEC
TYPE(TCFU), INTENT(INOUT) :: YDCFU
TYPE(TXFU), INTENT(INOUT) :: YDXFU
TYPE(TRIP), INTENT(IN) :: YDRIP
END SUBROUTINE ARO_SURF_DIAGH
END INTERFACE

