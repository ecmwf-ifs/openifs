! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUALGES(YDGEOMETRY,YD_JB_STRUCT,YDCVA_DATA)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMJG , ONLY : TYPE_JB_STRUCT
USE YOMCVA , ONLY : CVA_DATA_TYPE
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(TYPE_JB_STRUCT),INTENT(INOUT) :: YD_JB_STRUCT
TYPE(CVA_DATA_TYPE) ,INTENT(IN) :: YDCVA_DATA
call abor1("SUALGES should never be called with this build configuration - EXITING")
END SUBROUTINE SUALGES
