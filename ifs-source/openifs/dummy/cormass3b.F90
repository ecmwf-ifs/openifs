! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CORMASS3B(YDGEOMETRY,YDEDYN,YDSPEC)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YEMDYN , ONLY : TEDYN
use spectral_fields_mod, only:&
 & spectral_field
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(TEDYN) ,INTENT(INOUT) :: YDEDYN
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSPEC
call abor1("cormass3b.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE CORMASS3B
