! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUOBSB(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDRIP,YDVARBC,YDODB,YDTCV5,YDTCV_BGC)
USE YOMRIP , ONLY : TRIP
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE FIELDS_MOD , ONLY : FIELDS
USE MTRAJ_MOD , ONLY : MTRAJ
USE VARBC_CLASS,ONLY: CLASS_VARBC
USE TOVSCV_MOD , ONLY : TOVSCV
USE TOVSCV_BGC_MOD ,ONLY : TOVSCV_BGC
USE DBASE_MOD, ONLY : DBASE
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(FIELDS) ,INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ) ,INTENT(IN) :: YDMTRAJ
TYPE(TRIP) ,INTENT(INOUT) :: YDRIP
TYPE(CLASS_VARBC),INTENT(INOUT) :: YDVARBC
CLASS(DBASE) ,INTENT(INOUT) :: YDODB
TYPE(TOVSCV),OPTIONAL,INTENT(INOUT) :: YDTCV5
TYPE(TOVSCV_BGC),OPTIONAL,INTENT(INOUT) :: YDTCV_BGC
call abor1("SUOBSB should never be called with this build configuration - EXITING")
END SUBROUTINE SUOBSB
