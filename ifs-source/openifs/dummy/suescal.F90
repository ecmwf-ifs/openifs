! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUESCAL(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,STRUCT,YDVARBC,YD_JB_STRUCT,YDTCV)
USE TYPE_MODEL , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
use fields_mod , only:&
 & fields
USE MTRAJ_MOD , ONLY : MTRAJ
use yomcva , only:&
 & scalp_struct_type
USE VARBC_CLASS, ONLY : CLASS_VARBC
USE YOMJG , ONLY : TYPE_JB_STRUCT
USE TOVSCV_MOD, ONLY : TOVSCV
USE CONTROL_VECTORS_COMM_MOD
USE CONTROL_VECTORS_MOD
TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(FIELDS), INTENT(INOUT) :: YDFIELDS
TYPE(MODEL), INTENT(INOUT) :: YDMODEL
TYPE(MTRAJ), INTENT(INOUT) :: YDMTRAJ
TYPE(CLASS_VARBC), INTENT(INOUT) :: YDVARBC
TYPE(TYPE_JB_STRUCT), INTENT(INOUT) :: YD_JB_STRUCT
TYPE(TOVSCV),OPTIONAL, INTENT(INOUT) :: YDTCV
TYPE(SCALP_STRUCT_TYPE), POINTER :: STRUCT
call abor1("suescal.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUESCAL
