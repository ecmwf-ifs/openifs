! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE PPLTEMP(YDSTA,KPROMA,KST,KND,KFLEV,PGEO,PT,PXTEMP,LDTOP,PXGEO)
USE PARKIND1 , ONLY : JPIM, JPRB
use yomsta , only:&
 & tsta
TYPE(TSTA), INTENT(IN) :: YDSTA
INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN) :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KST
INTEGER(KIND=JPIM),INTENT(IN) :: KND
REAL(KIND=JPRB) ,INTENT(IN) :: PGEO(KPROMA,KFLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PT(KPROMA,KFLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PXTEMP
LOGICAL ,INTENT(IN) :: LDTOP
REAL(KIND=JPRB) ,INTENT(OUT) :: PXGEO(KPROMA)
call abor1("PPLTEMP should never be called with this build configuration - EXITING")
END SUBROUTINE PPLTEMP
