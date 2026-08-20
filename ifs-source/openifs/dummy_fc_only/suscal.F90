! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUSCAL(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,STRUCT,YDVARBC,YD_JB_STRUCT)
USE TYPE_MODEL , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
use fields_mod , only:&
 & fields
USE MTRAJ_MOD , ONLY : MTRAJ
use yomcva , only:&
 & scalp_struct_type
USE YOMJG , ONLY : TYPE_JB_STRUCT
USE VARBC_CLASS,ONLY: CLASS_VARBC
TYPE(GEOMETRY) ,INTENT(INOUT) :: YDGEOMETRY
TYPE(FIELDS) ,INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ) ,INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
TYPE(CLASS_VARBC) ,INTENT(INOUT) :: YDVARBC
TYPE(TYPE_JB_STRUCT),INTENT(INOUT) :: YD_JB_STRUCT
TYPE(SCALP_STRUCT_TYPE), POINTER, INTENT(INOUT) :: STRUCT
call abor1("SUSCAL should never be called with this build configuration - EXITING")
END SUBROUTINE SUSCAL
