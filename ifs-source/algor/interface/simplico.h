! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

INTERFACE
SUBROUTINE SIMPLICO(KM,KSMAX,KFLEV,KFLSUR,PALPHA,PDENIM,&
 & PFPLUS,PFMINUS,PSIVP,PRLAPDI,PBDT2,PY,PX) 
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KM
INTEGER(KIND=JPIM),INTENT(IN) :: KSMAX
INTEGER(KIND=JPIM),INTENT(IN) :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KFLSUR
REAL(KIND=JPRB) ,INTENT(IN) :: PALPHA(KM:KSMAX+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PDENIM(KM:KSMAX+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PFPLUS(KM:KSMAX+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PFMINUS(KM:KSMAX+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PSIVP(KFLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PRLAPDI(0:KSMAX)
REAL(KIND=JPRB) ,INTENT(IN) :: PBDT2
REAL(KIND=JPRB) ,INTENT(INOUT) :: PY(KFLSUR,2,KM:KSMAX)
REAL(KIND=JPRB) ,INTENT(OUT) :: PX(KFLSUR,2,KM:KSMAX)
END SUBROUTINE SIMPLICO
END INTERFACE
