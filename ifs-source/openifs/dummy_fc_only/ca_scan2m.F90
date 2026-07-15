! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CA_SCAN2M(YDGEOMETRY,YDFIELDS,YDMODEL,YDODB)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE FIELDS_MOD , ONLY : FIELDS
USE TYPE_MODEL , ONLY : MODEL
USE DBASE_MOD , ONLY : DBASE
CLASS(DBASE), INTENT(INOUT) :: YDODB
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(FIELDS) ,INTENT(INOUT) :: YDFIELDS
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
call abor1("CA_SCAN2M should never be called with this build configuration - EXITING")
END SUBROUTINE CA_SCAN2M
