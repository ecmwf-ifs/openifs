! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUALCTV(YDGEOMETRY,YDCVA_STRUCT, YDCVA_DATA, YD_JB_STRUCT)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMCVA , ONLY : CVA_DATA_TYPE
USE YOMJG , ONLY : TYPE_JB_STRUCT
use control_vectors_data_mix, only:&
 & control_vector_data_struct
TYPE(GEOMETRY),INTENT(IN) :: YDGEOMETRY
TYPE(CONTROL_VECTOR_DATA_STRUCT), POINTER, INTENT(IN) :: YDCVA_STRUCT
TYPE(CVA_DATA_TYPE), POINTER, INTENT(INOUT) :: YDCVA_DATA
TYPE(TYPE_JB_STRUCT), INTENT(INOUT) :: YD_JB_STRUCT
call abor1("SUALCTV should never be called with this build configuration - EXITING")
END SUBROUTINE SUALCTV
