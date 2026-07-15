! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUJB(YDGEOMETRY,YGFL,YDJBS,YDCVA_DATA)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMJG , ONLY : TYPE_JB_STRUCT
USE YOMCVA , ONLY : CVA_DATA_TYPE
USE YOM_YGFL , ONLY : TYPE_GFLD
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(TYPE_GFLD) ,INTENT(INOUT) :: YGFL
TYPE(TYPE_JB_STRUCT),INTENT(INOUT) :: YDJBS
TYPE(CVA_DATA_TYPE) ,INTENT(INOUT) :: YDCVA_DATA
call abor1("SUJB should never be called with this build configuration - EXITING")
END SUBROUTINE SUJB
