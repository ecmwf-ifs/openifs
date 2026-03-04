! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUNHSI_TESTCONV(YDGEOMETRY,YDRIP,YDDYN,YDDYNA)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYN , ONLY : TDYN
USE YOMDYNA , ONLY : TDYNA
USE YOMRIP , ONLY : TRIP
TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TDYN) ,INTENT(INOUT) :: YDDYN
TYPE(TDYNA) ,INTENT(INOUT) :: YDDYNA
TYPE(TRIP) ,INTENT(INOUT) :: YDRIP
call abor1("sunhsi_testconv.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUNHSI_TESTCONV
