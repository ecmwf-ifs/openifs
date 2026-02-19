! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

INTERFACE
SUBROUTINE EGGMLT (PRPI, KDLUX, KDLUN, KDGUX, KDGUN, KULOUT,&
 & KPRINT, PRPK, PLON0U, PLON1U, PLON2U, KSOTRP, PLAT1R, PLAT2R,&
 & PHSUD,PBETA)  

!--------------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
!--------------------------------------------------------------------------
IMPLICIT NONE
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPI 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPRINT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPK 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLON0U 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLON1U 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLON2U 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTRP 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAT1R 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAT2R 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHSUD 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBETA 
!--------------------------------------------------------------------------
END SUBROUTINE EGGMLT
END INTERFACE
