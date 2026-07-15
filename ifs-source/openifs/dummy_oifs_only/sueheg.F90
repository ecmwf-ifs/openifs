! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUEHEG(YDGEOMETRY,YDDYN,YDDYNA,YDEDYN,YDRIP)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYN , ONLY : TDYN
USE YOMDYNA , ONLY : TDYNA
USE YEMDYN , ONLY : TEDYN
USE YOMRIP , ONLY : TRIP
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TDYN) ,INTENT(INOUT):: YDDYN
TYPE(TDYNA) ,INTENT(INOUT):: YDDYNA
TYPE(TEDYN) ,INTENT(INOUT):: YDEDYN
TYPE(TRIP) ,INTENT(INOUT):: YDRIP
call abor1("SUEHEG should never be called with this build configuration - EXITING")
END SUBROUTINE SUEHEG
