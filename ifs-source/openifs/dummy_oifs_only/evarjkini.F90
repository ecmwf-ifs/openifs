! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE EVARJKINI(YDGEOMETRY,YDGFL,YDGFL5,YDDYN,YDDYNA,YDML_GCONF,YDML_LBC,YDSP)
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMDYN , ONLY : TDYN
USE YOMDYNA , ONLY : TDYNA
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGFL , ONLY : TGFL
USE YEMLBC_MODEL, ONLY : TELBC_MODEL
USE SPECTRAL_FIELDS_DATA, ONLY: SPECTRAL_FIELD
TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TGFL), INTENT(INOUT) :: YDGFL
TYPE(TGFL), INTENT(INOUT) :: YDGFL5
TYPE(TDYN) ,INTENT(INOUT) :: YDDYN
TYPE(TDYNA) ,INTENT(INOUT) :: YDDYNA
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TELBC_MODEL), INTENT(INOUT) :: YDML_LBC
TYPE(SPECTRAL_FIELD), INTENT(INOUT) :: YDSP
call abor1("EVARJKINI should never be called with this build configuration - EXITING")
END SUBROUTINE EVARJKINI
