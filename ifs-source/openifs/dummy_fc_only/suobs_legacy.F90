! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUOBS_LEGACY(YDGEOMETRY,YDEPHY,YDEAERLID,YDML_GCONF,YDPHY2,YDODB)
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOEAERLID , ONLY : TEAERLID
USE YOMPHY2 , ONLY : TPHY2
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOEPHY , ONLY : TEPHY
USE YOMDB
USE DBASE_MOD , ONLY : DBASE
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(TEAERLID) ,INTENT(INOUT) :: YDEAERLID
TYPE(TEPHY) ,INTENT(INOUT) :: YDEPHY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT) :: YDML_GCONF
TYPE(TPHY2) ,INTENT(INOUT) :: YDPHY2
CLASS(DBASE) ,INTENT(INOUT) :: YDODB
call abor1("SUOBS_LEGACY should never be called with this build configuration - EXITING")
END SUBROUTINE SUOBS_LEGACY
