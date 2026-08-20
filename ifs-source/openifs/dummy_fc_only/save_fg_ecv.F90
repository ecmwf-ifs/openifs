! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SAVE_FG_ECV(YDGEOMETRY,YDFIELDS)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE FIELD_CONTAINER_GP_MOD, ONLY : FIELD_CONTAINER_GP
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(FIELD_CONTAINER_GP) , INTENT(INOUT) :: YDFIELDS
call abor1("SAVE_FG_ECV should never be called with this build configuration - EXITING")
END SUBROUTINE SAVE_FG_ECV
