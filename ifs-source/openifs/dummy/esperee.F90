! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ESPEREE(YDGEOMETRY,KFLSUR,KFIELD,PSPEC,PREEL,PREELL,PREELM)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 ,ONLY : JPIM ,JPRB
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN) :: KFLSUR
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELD
REAL(KIND=JPRB),INTENT(IN) :: PSPEC(KFLSUR,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB),INTENT(OUT) :: PREEL(YDGEOMETRY%YRGEM%NGPTOT,KFIELD,1)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PREELL(YDGEOMETRY%YRGEM%NGPTOT,KFIELD,1)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PREELM(YDGEOMETRY%YRGEM%NGPTOT,KFIELD,1)
call abor1("esperee.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ESPEREE
