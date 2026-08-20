! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GNHY(&
 & YDGEOMETRY,YDDYNA,KPROMA,KSTART,KEND,POROGL,POROGM,&
 & PSP,PUF,PVF,&
 & PNHY)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYNA , ONLY : TDYNA
USE PARKIND1 , ONLY : JPIM, JPRB
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TDYNA), INTENT(IN) :: YDDYNA
INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN) :: KSTART
INTEGER(KIND=JPIM),INTENT(IN) :: KEND
REAL(KIND=JPRB) ,INTENT(IN) :: POROGL (KPROMA)
REAL(KIND=JPRB) ,INTENT(IN) :: POROGM (KPROMA)
REAL(KIND=JPRB) ,INTENT(IN) :: PSP (KPROMA)
REAL(KIND=JPRB) ,INTENT(IN) :: PUF (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ,INTENT(IN) :: PVF (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ,INTENT(OUT) :: PNHY (KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
call abor1("GNHY should never be called with this build configuration - EXITING")
END SUBROUTINE GNHY
