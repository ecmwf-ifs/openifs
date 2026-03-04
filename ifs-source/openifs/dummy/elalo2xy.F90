! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ELALO2XY(KGP,PLON0,PLAT0,PLON1,PLAT1,PLONC,PLATC,LDMRT,&
 & PGELAM,PGELAT,PX,PY,PGNORX,PGNORY)
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KGP
REAL(KIND=JPRB) ,INTENT(IN) :: PLON0
REAL(KIND=JPRB) ,INTENT(IN) :: PLAT0
REAL(KIND=JPRB) ,INTENT(IN) :: PLON1
REAL(KIND=JPRB) ,INTENT(IN) :: PLAT1
REAL(KIND=JPRB) ,INTENT(IN) :: PLONC
REAL(KIND=JPRB) ,INTENT(IN) :: PLATC
LOGICAL ,INTENT(IN) :: LDMRT
REAL(KIND=JPRB) ,DIMENSION(KGP) ,INTENT(IN) :: PGELAM
REAL(KIND=JPRB) ,DIMENSION(KGP) ,INTENT(IN) :: PGELAT
REAL(KIND=JPRB) ,DIMENSION(KGP) ,INTENT(OUT) :: PX
REAL(KIND=JPRB) ,DIMENSION(KGP) ,INTENT(OUT) :: PY
REAL(KIND=JPRB) ,DIMENSION(KGP) ,INTENT(OUT) :: PGNORX
REAL(KIND=JPRB) ,DIMENSION(KGP) ,INTENT(OUT) :: PGNORY
call abor1("elalo2xy.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ELALO2XY
