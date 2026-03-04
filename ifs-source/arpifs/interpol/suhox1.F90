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

SUBROUTINE SUHOX1(KSLWIDE,KDGSA,KDGEN,KFRSTLOFF,KPROMA,KEND,P4JP,&
                 &PI,KLOEN,PLATIN,PLAT,PLON,KLA,KLO)

!**** *SEHOX1*  - Compute the coordinates of points of a column
!                 in the halo near a target interpolation point
!                 of Fullpos.

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014

!     Description.
!     ------------
!     This routine is very similar to SUHOW1; all arguments have 
!     the same meaning as those of SUHOW1, except that KLO holds 
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: P4JP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI
INTEGER(KIND=JPIM),INTENT(IN)    :: KLOEN(KDGSA+KFRSTLOFF:KDGEN+KFRSTLOFF)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLATIN(KDGSA+KFRSTLOFF:KDGEN+KFRSTLOFF)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT(KEND)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLON(KEND)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLA(KPROMA)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLO(KPROMA,2*KSLWIDE)

REAL(KIND=JPRB) :: ZLSDEPI(KDGSA:KDGEN)

INTEGER(KIND=JPIM) :: IGLGLO, ILA, IZLAT, IZLATG, JGL, JLEN, JLO

REAL(KIND=JPRB) :: ZDEPI, ZLOO
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUHOX1',0,ZHOOK_HANDLE)

ZDEPI=2.0_JPRB*PI

DO JGL=KDGSA,KDGEN
  IGLGLO=JGL+KFRSTLOFF
  ZLSDEPI(JGL)=REAL(KLOEN(IGLGLO),JPRB)/ZDEPI
ENDDO

DO JLEN = 1,KEND

  IZLAT =INT(P4JP*(0.5_JPRB*PI-PLAT(JLEN))+0.75_JPRB)-KFRSTLOFF
  IZLATG=IZLAT+KFRSTLOFF
  ILA   =IZLAT+NINT(SIGN(0.5_JPRB,PLATIN(IZLATG)-PLAT(JLEN))-1.5_JPRB)
  KLA(JLEN)=ILA+1

  ZLOO = MOD(PLON(JLEN),ZDEPI)

  DO JLO = 1, 2*KSLWIDE
    KLO (JLEN,JLO)=INT(ZLOO*ZLSDEPI(ILA+(JLO-KSLWIDE)))
  ENDDO

ENDDO

IF (LHOOK) CALL DR_HOOK('SUHOX1',1,ZHOOK_HANDLE)
END SUBROUTINE SUHOX1

