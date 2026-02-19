! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUECGES(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YD_BG,YD_JB_STRUCT)
USE TYPE_MODEL , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
use fields_mod , only:&
 & fields
USE MTRAJ_MOD , ONLY : MTRAJ
USE YOMJG , ONLY : TYPE_JB_STRUCT
use yomtraj , only:&
 & traj_type
TYPE(GEOMETRY) ,INTENT(INOUT) :: YDGEOMETRY
TYPE(FIELDS) ,INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ) ,INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
TYPE(TRAJ_TYPE) ,INTENT(INOUT) :: YD_BG
TYPE(TYPE_JB_STRUCT),INTENT(INOUT) :: YD_JB_STRUCT
call abor1("suecges.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUECGES
