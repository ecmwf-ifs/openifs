! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ADD_MODERR_TL(YDGEOMETRY,YDFIELDS,YDML_GCONF,YDDYNA,YDMODERRCONF,YDSPERR,YDGPERR,KSTEP)
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE FIELDS_MOD , ONLY : FIELDS
use parkind1 , only:&
 & jpim
USE YOMDYNA , ONLY : TDYNA
use gridpoint_fields_mix , only:&
 & gridpoint_field
use spectral_fields_mod , only:&
 & spectral_field
USE YOMMODERRCONF , ONLY : TMODERR_CONF
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(FIELDS) ,INTENT(INOUT) :: YDFIELDS
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TDYNA) ,INTENT(IN) :: YDDYNA
TYPE(TMODERR_CONF) ,INTENT(IN) :: YDMODERRCONF
TYPE(SPECTRAL_FIELD) ,INTENT(IN) :: YDSPERR(:)
TYPE(GRIDPOINT_FIELD),INTENT(IN) :: YDGPERR(:)
INTEGER(KIND=JPIM) ,INTENT(IN) :: KSTEP
call abor1("ADD_MODERR_TL should never be called with this build configuration - EXITING")
END SUBROUTINE ADD_MODERR_TL
