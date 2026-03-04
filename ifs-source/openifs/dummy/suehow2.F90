! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUEHOW2(KPROMA,KEND,PCR,PRGMSD,PGMSF,KBINL,LDML,&
 & PDLAT,PDLO,PDELY,PWXX,PWXY)
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN) :: KEND
REAL(KIND=JPRB) ,INTENT(IN) :: PCR(2)
REAL(KIND=JPRB) ,INTENT(IN) :: PRGMSD(KEND)
REAL(KIND=JPRB) ,INTENT(IN) :: PGMSF(KEND)
INTEGER(KIND=JPIM),INTENT(IN) :: KBINL
LOGICAL ,INTENT(IN) :: LDML
REAL(KIND=JPRB) ,INTENT(IN) :: PDLAT(4,KPROMA)
REAL(KIND=JPRB) ,INTENT(IN) :: PDLO(0:3,4,KPROMA)
REAL(KIND=JPRB) ,INTENT(IN) :: PDELY
REAL(KIND=JPRB) ,INTENT(OUT) :: PWXX(KPROMA,KBINL)
REAL(KIND=JPRB) ,INTENT(OUT) :: PWXY(KPROMA,5:16)
call abor1("suehow2.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUEHOW2
