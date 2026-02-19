! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUEINFCE(YDGEOMETRY,YD_JB_STRUCT)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMJG, ONLY : TYPE_JB_STRUCT
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TYPE_JB_STRUCT), INTENT(INOUT) :: YD_JB_STRUCT
call abor1("sueinfce.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUEINFCE
