! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE INCLI0(YDGEOMETRY,YDSURF,YDGFL,YDEPHY,YDML_PHY_MF,YDMCC)
USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE YOEPHY , ONLY : TEPHY
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMMCC , ONLY : TMCC
USE YOM_YGFL , ONLY : TYPE_GFLD
TYPE(GEOMETRY),INTENT(INOUT) :: YDGEOMETRY
TYPE(TSURF) ,INTENT(INOUT) :: YDSURF
TYPE(TYPE_GFLD) ,INTENT(INOUT):: YDGFL
TYPE(TEPHY) ,INTENT(INOUT) :: YDEPHY
TYPE(TMCC) ,INTENT(INOUT) :: YDMCC
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
call abor1("incli0.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE INCLI0
