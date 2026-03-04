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

SUBROUTINE SUEHOX1(KSLWIDE,KDGSA,KDGEN,KFRSTLOFF,KPROMA,KEND,PLAT,PLON,PDELX,&
                  &PDELY,KDLUN,KDLUX,KDGUN,KDGUX,KLA,KLO)

!**** *SUEHOX1*  - Compute the coordinates of points of a column
!                  in the halo near a target interpolation point
!                  of Fullpos.

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014

!     Description.
!     ------------
!     This routine is very similar to SUEHOW1; all arguments have 
!     the same meaning as those of SUEHOW1, except that KLO holds 
!     coordinates of a full column of length 2*KSLWIDE of the
!     halo.

!         KSLWIDE  : width of halo part of interpolation buffer


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KSLWIDE
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGSA
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGEN
INTEGER(KIND=JPIM),INTENT(IN)    :: KFRSTLOFF
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT(KEND)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLON(KEND)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELY
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUN
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUX
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUN
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUX
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLA(KPROMA)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLO(KPROMA,2*KSLWIDE)

REAL(KIND=JPRB) :: ZYIN (KDGSA:KDGEN)

INTEGER(KIND=JPIM) :: ILA, ILAG, JLEN, JY  
INTEGER(KIND=JPIM) :: ILO

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUEHOX1',0,ZHOOK_HANDLE)

DO JY = KDGSA,KDGEN
  ZYIN(JY) = PDELY*REAL(JY+KFRSTLOFF-KDGUN,JPRB)
ENDDO

DO JLEN = 1,KEND
  ILAG = INT(PLAT(JLEN)/PDELY) + 1
  ILAG = MIN(MAX(ILAG,KDGUN),KDGUX)
  ILA  = ILAG - KFRSTLOFF
  KLA(JLEN) = ILA
ENDDO

DO JLEN = 1,KEND
  ILO=INT(PLON(JLEN)/PDELX)+1
  ILO=MIN(MAX(ILO,KDLUN),KDLUX)
  KLO (JLEN,:) = ILO-1
ENDDO

IF (LHOOK) CALL DR_HOOK('SUEHOX1',1,ZHOOK_HANDLE)
END SUBROUTINE SUEHOX1

