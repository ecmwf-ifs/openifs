! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUMODERRCONFPTRS(YDMODERRCONF,YDSPMODERR,YDGPMODERR,YDGEOMETRY,YDML_GCONF,YDGMV)
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMGMV, ONLY : TGMV
USE YOMMODERRCONF, ONLY : TMODERR_CONF
USE GRIDPOINT_FIELDS_MIX, ONLY : GRIDPOINT_FIELD
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD
TYPE(TMODERR_CONF), INTENT(INOUT) :: YDMODERRCONF
TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSPMODERR(:)
TYPE(GRIDPOINT_FIELD), INTENT(IN) :: YDGPMODERR(:)
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN) :: YDML_GCONF
TYPE(TGMV), INTENT(IN) :: YDGMV
call abor1("SUMODERRCONFPTRS should never be called with this build configuration - EXITING")
END SUBROUTINE SUMODERRCONFPTRS
