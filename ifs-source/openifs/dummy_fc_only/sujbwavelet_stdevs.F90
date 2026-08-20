! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUJBWAVELET_STDEVS(YDGEOMETRY,YDML_GCONF,YDCHEM,YD_JB_STRUCT)
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMJG, ONLY : TYPE_JB_STRUCT
USE YOMCHEM, ONLY : TCHEM
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT) :: YDML_GCONF
TYPE(TCHEM) ,INTENT(INOUT) :: YDCHEM
TYPE(TYPE_JB_STRUCT) ,INTENT(INOUT) :: YD_JB_STRUCT
call abor1("SUJBWAVELET_STDEVS should never be called with this build configuration - EXITING")
END SUBROUTINE SUJBWAVELET_STDEVS
