! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE PPCVIRT(KPROMA,KSTART,KPROF,KFLEV,PT,PQ,PQL,PQR,PQS,&
 & PQG,PQI,PTV)
USE PARKIND1 ,ONLY : JPIM, JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN) :: KSTART
INTEGER(KIND=JPIM),INTENT(IN) :: KPROF
INTEGER(KIND=JPIM),INTENT(IN) :: KFLEV
REAL(KIND=JPRB),INTENT(IN) :: PQ(KPROMA,KFLEV)
REAL(KIND=JPRB),INTENT(IN) :: PT(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL ,INTENT(IN) :: PQL(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL ,INTENT(IN) :: PQR(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL ,INTENT(IN) :: PQS(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL ,INTENT(IN) :: PQG(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL ,INTENT(IN) :: PQI(KPROMA,KFLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PTV(KPROMA,KFLEV)
call abor1("PPCVIRT should never be called with this build configuration - EXITING")
END SUBROUTINE PPCVIRT
