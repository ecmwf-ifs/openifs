! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE POAERO(PZ,PT,PUW,PVW,PTROP,PTROT,PUJET,PVJET,PJET,&
 & PNIVP,KDIM,KBEGIN,KEND,KMX,KPOL2)
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KDIM
INTEGER(KIND=JPIM),INTENT(IN) :: KMX
REAL(KIND=JPRB) ,INTENT(IN) :: PZ(KDIM,KMX)
REAL(KIND=JPRB) ,INTENT(IN) :: PT(KDIM,KMX)
REAL(KIND=JPRB) ,INTENT(IN) :: PUW(KDIM,KMX)
REAL(KIND=JPRB) ,INTENT(IN) :: PVW(KDIM,KMX)
REAL(KIND=JPRB) ,INTENT(OUT) :: PTROP(KDIM)
REAL(KIND=JPRB) ,INTENT(OUT) :: PTROT(KDIM)
REAL(KIND=JPRB) ,INTENT(OUT) :: PUJET(KDIM)
REAL(KIND=JPRB) ,INTENT(OUT) :: PVJET(KDIM)
REAL(KIND=JPRB) ,INTENT(OUT) :: PJET(KDIM)
REAL(KIND=JPRB) ,INTENT(IN) :: PNIVP(KDIM,KMX)
INTEGER(KIND=JPIM),INTENT(IN) :: KBEGIN
INTEGER(KIND=JPIM),INTENT(IN) :: KEND
INTEGER(KIND=JPIM),INTENT(OUT) :: KPOL2
call abor1("POAERO should never be called with this build configuration - EXITING")
END SUBROUTINE POAERO
