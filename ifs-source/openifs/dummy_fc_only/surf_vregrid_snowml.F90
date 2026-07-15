! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SURF_VREGRID_SNOWML(YDGEOMTRY,YSURF,YDSURF,FIELD_INOUT,KFIDS)
use parkind1 , only:&
 & jpim
USE FIELD_CONTAINER_GP_MOD, ONLY : FIELD_CONTAINER_GP
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE ISO_C_BINDING
USE SRFSN_VGRID_MOD
USE SRFSN_REGRID_MOD
TYPE(GEOMETRY) , INTENT(IN) :: YDGEOMTRY
TYPE(C_PTR) , INTENT(IN) :: YSURF
TYPE(TSURF) , INTENT(IN) :: YDSURF
TYPE(FIELD_CONTAINER_GP) , INTENT(INOUT) :: FIELD_INOUT
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFIDS(:)
call abor1("SURF_VREGRID_SNOWML should never be called with this build configuration - EXITING")
END SUBROUTINE SURF_VREGRID_SNOWML
