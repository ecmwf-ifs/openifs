! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE PERTOBS(YDRIP,YDODB,YDGEOMETRY)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMRIP , ONLY : TRIP
USE YOMDB
USE RANDOM_NUMBERS_MIX
USE DBASE_MOD, ONLY : DBASE
TYPE(TRIP) ,INTENT(INOUT) :: YDRIP
CLASS(DBASE), INTENT(INOUT) :: YDODB
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
call abor1("PERTOBS should never be called with this build configuration - EXITING")
END SUBROUTINE PERTOBS
