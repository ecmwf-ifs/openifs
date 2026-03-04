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

SUBROUTINE SUHOW1(KDGSA,KDGEN,KFRSTLOFF,KPROMA,KEND, &
 & P4JP,PI,KLOEN,PLATIN,KBINL,LDML,PLAT,PLON,&
 & KLA,KLON,KLONN,KLOS,KLOSS,PDLAT,PDLO)  

!**** *SUHOW1* - Computes the coordinates of nearest points of an
!****            interpolation point for FULL-POS. 

!     Purpose.
!     -------
!       Determines the interpolation grid for 16 points
!       interpolations for FULL-POS fields.
!       Numbering of the points (I is the interpolation point):
!                     13       5       6      14

!                      7       1       2       8
!                                 (I)
!                      9       3       4      10

!                     15      11      12      16
!       - computation of the latitude and the longitude of the
!         points numbered 13, 7, 9, 15.
!       - computation of linear weights.
!       Storage of coordinates and weights.

!       For distributed memory SUHOW1 is executed locally on each processor
!       (iproca,iprocb).
!       Abbreviation DM-global means 'variable global for all processors'.
!       Abbreviation DM-local means 'variable local for each processor
!        iproca and iprocb'.

!**   Interface.
!**   ---------
!**   CALL SUHOW1(...)
!**
!        Explicit arguments :
!        --------------------

!     INPUT:
!      KDGSA     : Lower DM-locally numbered latitude 'jgl' involved (DM-local).
!      KDGEN     : Upper DM-locally numbered latitude 'jgl' involved (DM-local).
!      KFRSTLOFF : Difference 'DM-global numbering of latitudes' -
!                  'DM-local numbering of latitudes'
!                  KFRSTLOFF=0 if not distributed memory (DM-local).
!      KPROMA    : Length of output vector, to which data are interpolated (DM-local).
!      KEND      : Actual length of output vector, to which data are interpolated (DM-local).
!      P4JP      : Approximative inverse of the differences of latitudes (DM-global).
!      PI        : Constant PI (DM-global).
!      KLOEN     : KLOEN(jglg): number of longitudes for the DM-globally 
!                  numbered latitude 'jglg'.
!                  (DM-global, but only elements KDGSA+KFRSTLOFF to 
!                  KDGEN+KFRSTLOFF are provided here).
!      PLATIN    : Gaussian/Lobatto latitudes.
!                  (DM-global, but only elements KDGSA+KFRSTLOFF to
!                  KDGEN+KFRSTLOFF are provided here).
!      KBINL     : Options for choice of interpolations (DM-global):
!                  = 0 ==> 12 points interpolation
!                  = 1 ==> bilinear interpolation
!      LDML      : Multilinear interpolations (DM-global).
!      PLAT      : Latitudes of output grid, in radians (DM-local).
!      PLON      : Longitudes of output grid, in radians (DM-local).

!     OUTPUT (all of them are DM-local):
!      KLA       : Latitude index of the point numbered 1.
!      KLON      : Longitude index of the point numbered 7.
!      KLONN     : Longitude index of the point numbered 13.
!      KLOS      : Longitude index of the point numbered 9.
!      KLOSS     : Longitude index of the point numbered 15.
!      PDLAT     : Meridian linear weights.
!      PDLO      : Zonal linear weights.

!        Implicit arguments :
!        --------------------
!        none

!     Method.
!     -------
!        See documentation about FULL-POS.

!     Externals.
!     ----------
!        No external.

!     Reference.
!     ----------

!     Author.
!     -------
!     Author : K. YESSAD after routines SUHOW and LASCAW.
!     Date   : MARCH 1997.

!     Modifications
!     -------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad (oct 2010): multi-linear interpolations.
!     ---------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ---------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KDGSA
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGEN
INTEGER(KIND=JPIM),INTENT(IN)    :: KFRSTLOFF
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
REAL(KIND=JPRB)   ,INTENT(IN)    :: P4JP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI
INTEGER(KIND=JPIM),INTENT(IN)    :: KLOEN(KDGSA+KFRSTLOFF:KDGEN+KFRSTLOFF)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLATIN(KDGSA+KFRSTLOFF:KDGEN+KFRSTLOFF)
INTEGER(KIND=JPIM),INTENT(IN)    :: KBINL
LOGICAL           ,INTENT(IN)    :: LDML
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT(KEND)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLON(KEND)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLA(KPROMA)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLON(KPROMA)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLONN(KPROMA)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLOS(KPROMA)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLOSS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLAT(4,KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLO(0:3,4,KPROMA)

!     ---------------------------------------------------------------------

REAL(KIND=JPRB) :: ZLSDEPI(KDGSA:KDGEN)

INTEGER(KIND=JPIM) :: IGLGLO, ILA, ILAG, ILO, ILO1, ILO2, ILO3,&
 & IZLAT, IZLATG, JGL, JLEN  

REAL(KIND=JPRB) :: ZDEPI, ZLO, ZLO1, ZLO2, ZLO3, ZLOO, ZTHIRD
REAL(KIND=JPRB) :: ZDLAT1, ZDLAT4
REAL(KIND=JPRB) :: ZDLO1_1, ZDLO2_1, ZDLO1_4, ZDLO2_4
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ---------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUHOW1',0,ZHOOK_HANDLE)
!     ---------------------------------------------------------------------

!*       0. INITIAL SETTING.
!           ----------------

ZDEPI=2.0_JPRB*PI
IF (LDML) ZTHIRD=1.0_JPRB/3.0_JPRB

!     Numbering of latitudes is DM-local in ZLSDEPI.
DO JGL=KDGSA,KDGEN
  IGLGLO=JGL+KFRSTLOFF
  ZLSDEPI(JGL)=REAL(KLOEN(IGLGLO),JPRB)/ZDEPI
ENDDO

!     ---------------------------------------------------------------------

!*       1. COORDINATES.
!           ------------

!        Latitudes of points 5,1,3,11 are always contiguous.

DO JLEN = 1,KEND

  IZLAT =INT(P4JP*(0.5_JPRB*PI-PLAT(JLEN))+0.75_JPRB)-KFRSTLOFF
  IZLATG=IZLAT+KFRSTLOFF
  ILA   =IZLAT+NINT(SIGN(0.5_JPRB,PLATIN(IZLATG)-PLAT(JLEN))-1.5_JPRB)
  ILAG  =ILA+KFRSTLOFF
  KLA(JLEN)=ILA+1
  ZDLAT1=(PLAT(JLEN)-PLATIN(ILAG+1))/(PLATIN(ILAG+2)-PLATIN(ILAG+1))
  PDLAT(1,JLEN)=ZDLAT1
  IF (LDML) THEN
    ZDLAT4=(PLAT(JLEN)-PLATIN(ILAG))/(PLATIN(ILAG+3)-PLATIN(ILAG))
    PDLAT(2,JLEN)=ZDLAT4
    PDLAT(3,JLEN)=ZDLAT1
    PDLAT(4,JLEN)=ZDLAT4
  ENDIF

  ZLOO = MOD(PLON(JLEN),ZDEPI)

  ZLO   =ZLOO*ZLSDEPI(ILA  )
  ILO   =INT(ZLO )
  KLONN(JLEN)=ILO

  ZLO1  =ZLOO*ZLSDEPI(ILA+1)
  ILO1  =INT(ZLO1)
  KLON(JLEN)=ILO1

  ZLO2  =ZLOO*ZLSDEPI(ILA+2)
  ILO2  =INT(ZLO2)
  KLOS(JLEN)=ILO2

  ZLO3  =ZLOO*ZLSDEPI(ILA+3)
  ILO3  =INT(ZLO3)
  KLOSS(JLEN)=ILO3

  IF(KBINL == 0) PDLO(0,1,JLEN)=ZLO -REAL(ILO ,JPRB)
  ZDLO1_1=ZLO1-REAL(ILO1,JPRB)
  PDLO(1,1,JLEN)=ZDLO1_1
  ZDLO2_1=ZLO2-REAL(ILO2,JPRB)
  PDLO(2,1,JLEN)=ZDLO2_1
  IF(KBINL == 0) PDLO(3,1,JLEN)=ZLO3-REAL(ILO3,JPRB)

  IF (LDML) THEN
    ZDLO1_4=ZTHIRD*(ZLO -REAL(ILO ,JPRB)+1.0_JPRB)
    ZDLO2_4=ZTHIRD*(ZLO3-REAL(ILO3,JPRB)+1.0_JPRB)
    PDLO(1,2,JLEN)=ZDLO1_1
    PDLO(2,2,JLEN)=ZDLO2_1
    PDLO(1,3,JLEN)=ZDLO1_4
    PDLO(2,3,JLEN)=ZDLO2_4
    PDLO(1,4,JLEN)=ZDLO1_4
    PDLO(2,4,JLEN)=ZDLO2_4
  ENDIF

ENDDO

!     ---------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUHOW1',1,ZHOOK_HANDLE)
END SUBROUTINE SUHOW1
