! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUEFPG3(KFPBOYD,KFPDOM,YDFPUSERGEO,YDGEOMETRY,KGP,KGP_DEP,LDWIDER,PLA_DEP,PLO_DEP,LDMASK)
use parkind1 , only:&
 & jpim,&
 & jprb
USE TYPE_FPUSERGEO, ONLY : TFPUSERGEO
USE GEOMETRY_MOD , ONLY : GEOMETRY
INTEGER(KIND=JPIM),INTENT(IN) :: KFPBOYD
INTEGER(KIND=JPIM),INTENT(IN) :: KFPDOM
TYPE (TFPUSERGEO), INTENT(IN) :: YDFPUSERGEO(KFPDOM)
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN) :: KGP
INTEGER(KIND=JPIM),INTENT(IN) :: KGP_DEP
LOGICAL, INTENT(IN) :: LDWIDER
REAL(KIND=JPRB) ,INTENT(INOUT) :: PLA_DEP(KGP_DEP)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PLO_DEP(KGP_DEP)
LOGICAL ,INTENT(OUT) :: LDMASK(KGP)
call abor1("suefpg3.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUEFPG3
