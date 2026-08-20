! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUALLT(YDGEOMETRY,YDMTRAJ,YDGMV,YDSURF,YDMODEL)
USE TYPE_MODEL , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE MTRAJ_MOD , ONLY : MTRAJ
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMGMV , ONLY : TGMV
TYPE(GEOMETRY),INTENT(IN) :: YDGEOMETRY
TYPE(TGMV) ,INTENT(INOUT) :: YDGMV
TYPE(MTRAJ) ,INTENT(INOUT) :: YDMTRAJ
TYPE(TSURF) ,INTENT(INOUT) :: YDSURF
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
call abor1("SUALLT should never be called with this build configuration - EXITING")
END SUBROUTINE SUALLT
