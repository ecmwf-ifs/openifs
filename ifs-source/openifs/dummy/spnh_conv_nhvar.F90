! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SPNH_CONV_NHVAR(YDGEOMETRY,YDML_GCONF,YDDYNA,LDMODEL_TO_FILE,YDSP)
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SPECTRAL_FIELDS_MOD , ONLY : SPECTRAL_FIELD
USE YOMDYNA , ONLY : TDYNA
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN) :: YDML_GCONF
TYPE(TDYNA) ,INTENT(IN) :: YDDYNA
LOGICAL ,INTENT(IN) :: LDMODEL_TO_FILE
TYPE(SPECTRAL_FIELD) ,INTENT(INOUT) :: YDSP
call abor1("spnh_conv_nhvar.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SPNH_CONV_NHVAR
