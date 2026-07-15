! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

 SUBROUTINE CHEM_INIT_TLAD(YDGEOMETRY,YDML_GCONF,YDML_CHEM)
USE MODEL_CHEM_MOD , ONLY : MODEL_CHEM_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(MODEL_CHEM_TYPE),INTENT(INOUT):: YDML_CHEM
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
call abor1("CHEM_INIT_TLAD should never be called with this build configuration - EXITING")
END SUBROUTINE CHEM_INIT_TLAD
