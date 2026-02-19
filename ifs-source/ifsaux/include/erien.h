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
SUBROUTINE ERIEN(CDNAMC,KTYPTR,LDMAP,&
 & KTRONC,KDGL,KNXLON,KNLOPA,PSINLA,&
 & KFLEV,PREF,PVALH,PVBH,&
 & PEPS,KULOUT)  
!     ------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
IMPLICIT NONE
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM), PARAMETER :: JPXGEO=18
INTEGER(KIND=JPIM), PARAMETER :: JPXPAH=8
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM),INTENT(INOUT) :: KFLEV 
CHARACTER(LEN=16),INTENT(IN)     :: CDNAMC
INTEGER(KIND=JPIM),INTENT(INOUT) :: KTYPTR 
LOGICAL           ,INTENT(INOUT) :: LDMAP 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KTRONC 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KDGL 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KNXLON 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KNLOPA(JPXPAH) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSINLA(JPXGEO) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PREF 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVALH(0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVBH(0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEPS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT
!     ------------------------------------------------------------------
END SUBROUTINE ERIEN
END INTERFACE
