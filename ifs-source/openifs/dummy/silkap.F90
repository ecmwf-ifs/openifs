! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SILKAP(YDGEOMETRY,YDDYN,LDNHQE_C2,KLEV,KLON,KNLON,PIN,POUT,PMULFAC)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYN , ONLY : TDYN
USE PARKIND1 , ONLY : JPIM, JPRB
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(TDYN) ,INTENT(IN) :: YDDYN
LOGICAL ,INTENT(IN) :: LDNHQE_C2
INTEGER(KIND=JPIM) ,INTENT(IN) :: KLEV
INTEGER(KIND=JPIM) ,INTENT(IN) :: KLON
INTEGER(KIND=JPIM) ,INTENT(IN) :: KNLON
REAL(KIND=JPRB) ,INTENT(IN) :: PIN(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) ,INTENT(OUT) :: POUT(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PMULFAC
call abor1("silkap.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SILKAP
