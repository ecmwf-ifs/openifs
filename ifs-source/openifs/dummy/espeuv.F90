! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ESPEUV(YDGEOMETRY,PETPSI,PDIKHI,PSPUB,PSPVB,PU,PV,KFLSUR,KFIELD,KINP)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 ,ONLY : JPIM ,JPRB
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN) :: KFLSUR
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELD
REAL(KIND=JPRB) ,INTENT(IN) :: PETPSI(KFLSUR,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) ,INTENT(IN) :: PDIKHI(KFLSUR,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) ,INTENT(IN) :: PSPUB(KFIELD)
REAL(KIND=JPRB) ,INTENT(IN) :: PSPVB(KFIELD)
REAL(KIND=JPRB) ,INTENT(OUT) :: PU(YDGEOMETRY%YRGEM%NGPTOT,KFIELD)
REAL(KIND=JPRB) ,INTENT(OUT) :: PV(YDGEOMETRY%YRGEM%NGPTOT,KFIELD)
INTEGER(KIND=JPIM),INTENT(IN) :: KINP
call abor1("espeuv.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ESPEUV
