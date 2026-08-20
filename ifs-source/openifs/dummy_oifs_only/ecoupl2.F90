! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ECOUPL2(YDGEOMETRY,YDMODEL,YDFIELDS,LD_DFISTEP)
USE TYPE_MODEL, ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE FIELDS_MOD, ONLY: FIELDS
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(MODEL) ,INTENT(IN) :: YDMODEL
TYPE(FIELDS) ,INTENT(INOUT) :: YDFIELDS
LOGICAL ,INTENT(IN) :: LD_DFISTEP
call abor1("ECOUPL2 should never be called with this build configuration - EXITING")
END SUBROUTINE ECOUPL2
