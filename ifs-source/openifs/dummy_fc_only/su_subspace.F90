! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SU_SUBSPACE(YDGEOMETRY,YDDIMF)
USE YOMDIMF , ONLY : TDIMF
USE GEOMETRY_MOD, ONLY : GEOMETRY
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TDIMF) ,INTENT(INOUT):: YDDIMF
call abor1("SU_SUBSPACE should never be called with this build configuration - EXITING")
END SUBROUTINE SU_SUBSPACE
