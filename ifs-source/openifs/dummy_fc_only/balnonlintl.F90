! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BALNONLINTL(YDGEOMETRY,PVOR,PZZP,PGP7, &
 & KOFF_VOR, KOFF_DIV, KOFF_U, KOFF_V, KOFF_EW_U, KOFF_EW_V)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
REAL(KIND=JPRB) ,INTENT(INOUT) :: PVOR(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) ,INTENT(OUT) :: PZZP(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) ,INTENT(IN) :: PGP7(YDGEOMETRY%YRDIM%NPROMA,9*YDGEOMETRY%YRDIMV%NFLEVG, &
 &(YDGEOMETRY%YRGEM%NGPTOT+YDGEOMETRY%YRDIM%NPROMA-1)/YDGEOMETRY%YRDIM%NPROMA)
INTEGER(KIND=JPIM), INTENT(IN) :: KOFF_VOR, KOFF_DIV, KOFF_U, KOFF_V, KOFF_EW_U, KOFF_EW_V
call abor1("BALNONLINTL should never be called with this build configuration - EXITING")
END SUBROUTINE BALNONLINTL
