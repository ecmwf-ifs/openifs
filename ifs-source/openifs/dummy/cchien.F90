! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CCHIEN(YDGEOMETRY,CDNOMA,KFILE,KINF)
USE GEOMETRY_MOD , ONLY : GEOMETRY
use parkind1 , only:&
 & jpim
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
CHARACTER(LEN=*), INTENT(IN) :: CDNOMA
INTEGER(KIND=JPIM), INTENT(IN) :: KFILE
INTEGER(KIND=JPIM), INTENT(INOUT) :: KINF
call abor1("cchien.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE CCHIEN
