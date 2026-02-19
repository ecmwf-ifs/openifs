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
SUBROUTINE RIEN(CDNAMC,KTYPTR,PSLAPO,PLOCEN,&
 & PCODIL,KTRONC,KDGL,KNXLON,KNLOPA,KNOZPA,PSINLA,&
 & KHTYP,KFLEV,PREF,PVALH,PVBH,KQUAD,&
 & KDGSA,KDGEN,PEPS,LDFICP,KULOUT)  
!     ------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
!     ------------------------------------------------------------------
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(INOUT) :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGSA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGEN 
CHARACTER(LEN=16) ,INTENT(IN)    :: CDNAMC
INTEGER(KIND=JPIM),INTENT(INOUT) :: KTYPTR 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSLAPO 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLOCEN 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCODIL 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KTRONC 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KDGL 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KNXLON 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KNLOPA(KDGSA:KDGEN) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KNOZPA(KDGSA:KDGEN) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSINLA(KDGSA:KDGEN) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KHTYP 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PREF 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVALH(0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVBH(0:KFLEV) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KQUAD 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEPS 
LOGICAL           ,INTENT(INOUT) :: LDFICP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
!     ------------------------------------------------------------------
END SUBROUTINE RIEN

END INTERFACE
