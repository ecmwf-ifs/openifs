! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUEHEG(YDGEOMETRY,YDDYN,YDEDYN,YDRIP)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYN , ONLY : TDYN
USE YEMDYN , ONLY : TEDYN
USE YOMRIP , ONLY : TRIP
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TDYN) ,INTENT(INOUT):: YDDYN
TYPE(TEDYN) ,INTENT(INOUT):: YDEDYN
TYPE(TRIP) ,INTENT(INOUT):: YDRIP
call abor1("sueheg.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUEHEG
