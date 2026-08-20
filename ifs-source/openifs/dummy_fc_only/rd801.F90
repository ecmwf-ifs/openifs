! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE RD801(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYN,YDDYNA,YDML_LBC,YDSP)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMDYN , ONLY : TDYN
USE YOMDYNA , ONLY : TDYNA
USE YEMLBC_MODEL , ONLY : TELBC_MODEL
USE YOMGFL , ONLY : TGFL
USE SPECTRAL_FIELDS_MOD , ONLY : SPECTRAL_FIELD
TYPE(GEOMETRY) ,INTENT(INOUT) :: YDGEOMETRY
TYPE(TGFL) ,INTENT(INOUT) :: YDGFL
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TDYN) ,INTENT(INOUT) :: YDDYN
TYPE(TDYNA) ,INTENT(INOUT) :: YDDYNA
TYPE(TELBC_MODEL) ,INTENT(INOUT) :: YDML_LBC
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSP
call abor1("RD801 should never be called with this build configuration - EXITING")
END SUBROUTINE RD801
