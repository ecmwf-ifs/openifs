! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUMODERR(YDGEOMETRY,YDGMV,YDML_GCONF)
USE MODEL_GENERAL_CONF_MOD, ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMGMV, ONLY : TGMV
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TGMV), INTENT(INOUT) :: YDGMV
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT) :: YDML_GCONF
call abor1("SUMODERR should never be called with this build configuration - EXITING")
END SUBROUTINE SUMODERR
