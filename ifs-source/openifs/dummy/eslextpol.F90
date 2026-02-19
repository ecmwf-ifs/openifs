! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ESLEXTPOL(YDGEOMETRY,YDSL,KFLDSLB,KFIXFLD,KTYP,PB1A,KMASK2)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE EINT_MOD , ONLY : SL_STRUCT
TYPE(GEOMETRY) , INTENT(IN) :: YDGEOMETRY
TYPE(SL_STRUCT), INTENT(IN) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN) :: KFLDSLB
INTEGER(KIND=JPIM),INTENT(IN) :: KFIXFLD(2)
INTEGER(KIND=JPIM),INTENT(IN) :: KTYP
REAL(KIND=JPRB) ,INTENT(INOUT) :: PB1A(YDSL%NASLB1,KFLDSLB,1)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KMASK2(YDSL%NASLB1+YDGEOMETRY%YRDIM%NSTENCILWIDE*2)
call abor1("eslextpol.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ESLEXTPOL
