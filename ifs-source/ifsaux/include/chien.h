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
SUBROUTINE CHIEN(CDNAMC,KTYPTR,PSLAPO,PLOCEN,&
 & PCODIL,KTRONC,KDGL,KNXLON,KNLOPA,KNOZPA,&
 & KHTYP,KFLEV,PREF,PVALH,PVBH,KQUAD,KINF,&
 & KDGSA,KDGEN,PEPS,LDFICP,KULOUT)  
!     ------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
!     ------------------------------------------------------------------
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)  :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGSA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGEN 
CHARACTER(LEN=16) ,INTENT(IN)  :: CDNAMC
INTEGER(KIND=JPIM),INTENT(IN)  :: KTYPTR 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PSLAPO 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PLOCEN 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PCODIL 
INTEGER(KIND=JPIM),INTENT(IN)  :: KTRONC 
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGL 
INTEGER(KIND=JPIM),INTENT(IN)  :: KNXLON 
INTEGER(KIND=JPIM),INTENT(IN)  :: KNLOPA(KDGSA:KDGEN) 
INTEGER(KIND=JPIM),INTENT(IN)  :: KNOZPA(KDGSA:KDGEN) 
INTEGER(KIND=JPIM),INTENT(IN)  :: KHTYP 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PREF 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PVALH(0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PVBH(0:KFLEV) 
INTEGER(KIND=JPIM),INTENT(IN)  :: KQUAD 
INTEGER(KIND=JPIM),INTENT(IN)  :: KINF 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PEPS 
LOGICAL           ,INTENT(OUT) :: LDFICP 
INTEGER(KIND=JPIM),INTENT(IN)  :: KULOUT 
!     ------------------------------------------------------------------
END SUBROUTINE CHIEN

END INTERFACE
