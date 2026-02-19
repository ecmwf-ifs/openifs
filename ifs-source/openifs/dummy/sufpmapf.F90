! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUFPMAPF(KFPDOM,YDFPUSERGEO,KFPRGPG_DEP,PGM)
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE TYPE_FPUSERGEO, ONLY : TFPUSERGEO
INTEGER(KIND=JPIM), INTENT(IN) :: KFPDOM
TYPE (TFPUSERGEO), INTENT(IN) :: YDFPUSERGEO(KFPDOM)
INTEGER(KIND=JPIM), INTENT(IN) :: KFPRGPG_DEP
REAL(KIND=JPRB), INTENT(INOUT) :: PGM(KFPRGPG_DEP)
call abor1("sufpmapf.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUFPMAPF
