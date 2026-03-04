! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ELSRW(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,LDIFI,LDREADCOU)
USE TYPE_MODEL , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE FIELDS_MOD , ONLY : FIELDS
USE MTRAJ_MOD , ONLY : MTRAJ
TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(FIELDS), INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ), INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
LOGICAL ,INTENT(IN) :: LDIFI
LOGICAL, OPTIONAL ,INTENT(OUT) :: LDREADCOU
call abor1("elsrw.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ELSRW
