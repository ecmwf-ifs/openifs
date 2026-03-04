! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ESPNORMB(YDLAP,YDLEP,YDDIM,YDEDIM,PX,KLEV,PSN,PSM,PMET,KFLEV,KFLSUR,LDNWAVE)
USE YOMDIM , ONLY : TDIM
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE YOMLAP , ONLY : TLAP
USE YEMLAP , ONLY : TLEP
USE YEMDIM , ONLY : TEDIM
TYPE(TLAP) , INTENT(IN) :: YDLAP
TYPE(TLEP) , INTENT(IN) :: YDLEP
TYPE(TDIM) , INTENT(IN) :: YDDIM
TYPE(TEDIM), INTENT(IN) :: YDEDIM
INTEGER(KIND=JPIM),INTENT(IN) :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KFLSUR
REAL(KIND=JPRB) ,INTENT(IN) :: PX(KFLSUR,YDDIM%NSPEC2)
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
REAL(KIND=JPRB) ,INTENT(OUT) :: PSM(KFLEV,YDDIM%NUMP)
REAL(KIND=JPRB) ,INTENT(IN) :: PMET(YDDIM%NSPECG)
LOGICAL ,INTENT(IN) :: LDNWAVE
REAL(KIND=JPRB) :: PSN(KFLEV,YDDIM%NUMP)
call abor1("espnormb.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ESPNORMB
