! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SIDD(YDGEOMETRY,YDDYN,KLEV,KLON,PDH,PDV,PRNH,PT,PSP,KNLON)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMDYN , ONLY : TDYN
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(TDYN) ,INTENT(IN) :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KNLON
REAL(KIND=JPRB) ,INTENT(INOUT) :: PDH(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) ,INTENT(OUT) :: PDV(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PRNH(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PT(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PSP(KNLON)
call abor1("sidd.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SIDD
