! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SISEVE(YDGEOMETRY,YDDYN,YDDYNA,KLEV,KLON,PV1,PV2,KNLON,LD_LSTAR)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMDYN , ONLY : TDYN
USE YOMDYNA , ONLY : TDYNA
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(TDYN) ,INTENT(IN) :: YDDYN
TYPE(TDYNA) ,INTENT(IN) :: YDDYNA
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KNLON
REAL(KIND=JPRB) ,INTENT(IN) :: PV1(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) ,INTENT(OUT) :: PV2(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
LOGICAL,OPTIONAL ,INTENT(IN) :: LD_LSTAR
call abor1("SISEVE should never be called with this build configuration - EXITING")
END SUBROUTINE SISEVE
