! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE READ_FG_ECV(YDGEOMETRY,YDFIELDS,LDBCK)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE FIELDS_MOD , ONLY : FIELDS
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(FIELDS) ,INTENT(INOUT) :: YDFIELDS
LOGICAL, OPTIONAL ,INTENT(IN) :: LDBCK
call abor1("read_fg_ecv.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE READ_FG_ECV
