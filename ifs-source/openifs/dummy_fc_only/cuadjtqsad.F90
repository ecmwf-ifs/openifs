! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUADJTQSAD &
 &(KIDIA, KFDIA, KLON, KLEV,&
 & KK,&
 & PSP5, PT5, PQ5,&
 & PSP, PT, PQ, LDFLAG, KCALL)
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KK
REAL(KIND=JPRB) ,INTENT(IN) :: PSP5(KLON)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PT5(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PQ5(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PSP(KLON)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PT(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PQ(KLON,KLEV)
LOGICAL ,INTENT(IN) :: LDFLAG(KLON)
INTEGER(KIND=JPIM),INTENT(IN) :: KCALL
call abor1("CUADJTQSAD should never be called with this build configuration - EXITING")
END SUBROUTINE CUADJTQSAD
