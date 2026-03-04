! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ACSOLW ( YDPHY1,KIDIA,KFDIA,KLON,&
 & PARG,PD2,PLSM,PIVEG,PSAB,&
 & LDHMT,&
 & PWFC,PWPMX,PWSAT,PWSMX,PWWILT)
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE YOMPHY1 , ONLY : TPHY1
TYPE(TPHY1) ,INTENT(IN) :: YDPHY1
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
REAL(KIND=JPRB) ,INTENT(IN) :: PARG(KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PD2(KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PLSM(KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PIVEG(KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PSAB(KLON)
LOGICAL ,INTENT(IN) :: LDHMT
REAL(KIND=JPRB) ,INTENT(INOUT) :: PWFC(KLON)
REAL(KIND=JPRB) ,INTENT(OUT) :: PWPMX(KLON)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PWSAT(KLON)
REAL(KIND=JPRB) ,INTENT(OUT) :: PWSMX(KLON)
REAL(KIND=JPRB) ,INTENT(OUT) :: PWWILT(KLON)
call abor1("acsolw.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ACSOLW
