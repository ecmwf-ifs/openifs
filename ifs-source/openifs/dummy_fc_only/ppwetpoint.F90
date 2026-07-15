! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE PPWETPOINT(YDCST, YDPHY,KIDIA,KFDIA,KLON,PAPRS,PT,PQV,PQL,PQI,PWETPOINT)
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE YOMPHY , ONLY : TPHY
USE YOMCST , ONLY : TCST
TYPE (TCST), INTENT (IN) :: YDCST
TYPE(TPHY) ,INTENT(IN) :: YDPHY
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
REAL(KIND=JPRB) ,INTENT(IN) :: PAPRS(KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PT(KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PQV(KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PQL(KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PQI(KLON)
REAL(KIND=JPRB) ,INTENT(OUT) :: PWETPOINT(KLON)
call abor1("PPWETPOINT should never be called with this build configuration - EXITING")
END SUBROUTINE PPWETPOINT
