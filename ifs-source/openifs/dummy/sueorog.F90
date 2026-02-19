! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUEOROG(YDGEOMETRY,PSPOR,POROG,POROGL,POROGM,POROGLL,POROGLM,POROGMM)
USE GEOMETRY_MOD , ONLY : GEOMETRY
use parkind1 , only:&
 & jprb
TYPE(GEOMETRY) ,INTENT(INOUT) :: YDGEOMETRY
REAL(KIND=JPRB) ,INTENT(IN) :: PSPOR(YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) ,INTENT(OUT) :: POROG (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) ,INTENT(OUT) :: POROGL (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) ,INTENT(OUT) :: POROGM (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) ,INTENT(OUT) :: POROGLL(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) ,INTENT(OUT) :: POROGMM(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) ,INTENT(OUT) :: POROGLM(YDGEOMETRY%YRGEM%NGPTOT)
call abor1("sueorog.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUEOROG
