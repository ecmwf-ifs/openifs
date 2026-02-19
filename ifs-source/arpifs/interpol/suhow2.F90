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

SUBROUTINE SUHOW2(KDGSA,KDGEN,KFRSTLOFF,KPROMA,KEND,PCR,PLATIN,PRGMSD,PGMSF,KBINL,LDML,KLA,&
 & PDLAT,PDLO,PWXX,PWXY)  

!**** *SUHOW2* - Compute the weights of nearest points.
!****
!****            No land/sea/ice mask is involved here.

!     Purpose.
!     -------
!       Computes the weights for 12 points interpolations or
!       bilinear interpolations for FULL-POS fields.
!       Numbering of the points (I is the interpolation point):
!                     13       5       6      14

!                      7       1       2       8
!                                 (I)
!                      9       3       4      10

!                     15      11      12      16
!       Points number 1 to 12 are used if 12 points interpolations.
!       Points number 1 to 4 are used if bilinear interpolations.
!       All points are used if multi-linear interpolations.
!       Storage of weights.

!       For distributed memory SUHOW2 is executed locally on each processor
!       (iproca,iprocb).
!       Abbreviation DM-global means 'variable global for all processors'.
!       Abbreviation DM-local means 'variable local for each processor
!        iproca and iprocb'.

!**   Interface.
!**   ---------
!**   CALL SUHOW2(...)
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
!      PCR       : Coefficients to compute ZR.
!      PLATIN    : Gaussian/Lobatto latitudes.
!                  (DM-global, but only elements KDGSA+KFRSTLOFF to
!                  KDGEN+KFRSTLOFF are provided here).
!      PRGMSD    : 1/(departure local mesh-size) (DM-local).
!      PGMSF     : Target local mesh-size (DM-local).
!      KBINL     : Options for choice of interpolations (DM-global).
!                  = 12 ==> 12 points interpolation
!                  =  4 ==> bilinear interpolation
!      LDML      : Multilinear interpolations (DM-global).
!      KLA       : Latitude index of the point numbered 1 (DM-local).
!      PDLAT     : Meridian linear weights (DM-local).
!      PDLO      : Zonal linear weights.

!     OUTPUT:
!      PWXX      : Weights array (DM-local).
!      PWXY      : Additional weights array for multi-linear interpolations (DM-local).

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCR(2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLATIN(KDGSA+KFRSTLOFF:KDGEN+KFRSTLOFF)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRGMSD(KEND)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMSF(KEND)
INTEGER(KIND=JPIM),INTENT(IN)    :: KBINL
LOGICAL           ,INTENT(IN)    :: LDML
INTEGER(KIND=JPIM),INTENT(IN)    :: KLA(KEND)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLAT(4,KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLO(0:3,4,KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWXX(KPROMA,KBINL)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWXY(KPROMA,5:16)

!     ---------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ILA1G, ILA2G, ILA3G, ILAG, JLEN, JBINL

REAL(KIND=JPRB) :: ZCXN, ZCXNN, ZCXS, ZCXSS, ZCY, ZDLAT, ZDLO1,&
 & ZDLO2, ZDY, ZDY10, ZDY21, ZDY32, ZWXN, ZWXN0, &
 & ZWXN1, ZWXN2, ZWXN3, ZWXNN, ZWXS, ZWXS0, &
 & ZWXS1, ZWXS2, ZWXS3, ZWXSS, ZWY, ZWY0, ZWY1, ZWY2, ZWY3
REAL(KIND=JPRB) :: ZR(KPROMA), ZRW(KPROMA,4)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ---------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUHOW2',0,ZHOOK_HANDLE)
!     ---------------------------------------------------------------------

!*       1. RESOLUTIONS RATIO.
!           ------------------

IF (LDML) THEN
  DO JLEN = 1,KEND
    ZR(JLEN)=PCR(1)*PRGMSD(JLEN)*PGMSF(JLEN)+PCR(2)
    ZRW(JLEN,4)=0.25_JPRB*MIN(1.0_JPRB,MAX(0.0_JPRB,ZR(JLEN)-2.0_JPRB))
    ZRW(JLEN,2)=(1.0_JPRB/3.0_JPRB)*(1.0_JPRB-ZRW(JLEN,4))*MIN(1.0_JPRB,MAX(0.0_JPRB,ZR(JLEN)-1.0_JPRB))
    ZRW(JLEN,3)=ZRW(JLEN,2)
    ZRW(JLEN,1)=1.0_JPRB-ZRW(JLEN,2)-ZRW(JLEN,3)-ZRW(JLEN,4)
  ENDDO
ENDIF

!     ---------------------------------------------------------------------

!*       2. WEIGHTS FOR BILINEAR INTERPOLATION.
!           -----------------------------------

IF (KBINL == 4) THEN

  DO JLEN = 1,KEND

    ZWY=PDLAT(1,JLEN)
    ZCY=1.0_JPRB-PDLAT(1,JLEN)

    ZWXN=PDLO(1,1,JLEN)
    ZCXN=1.0_JPRB-PDLO(1,1,JLEN)

    ZWXS=PDLO(2,1,JLEN)
    ZCXS=1.0_JPRB-PDLO(2,1,JLEN)

    PWXX(JLEN,1) = ZCXN*ZCY
    PWXX(JLEN,2) = ZWXN*ZCY
    PWXX(JLEN,3) = ZCXS*ZWY
    PWXX(JLEN,4) = ZWXS*ZWY

  ENDDO

ENDIF

!     ---------------------------------------------------------------------

!*       3. WEIGHTS FOR 12 POINTS INTERPOLATION.
!           ------------------------------------

IF (KBINL == 12) THEN

  DO JLEN = 1,KEND

    ILA1G=KLA(JLEN)+KFRSTLOFF
    ILAG =ILA1G-1
    ILA2G=ILA1G+1
    ILA3G=ILA1G+2

    ZDLAT = PDLAT(1,JLEN)
    ZDY   = ZDLAT*(PLATIN(ILA2G)-PLATIN(ILA1G))
    ZDY10 = PLATIN(ILA1G)-PLATIN(ILAG)
    ZDY21 = PLATIN(ILA2G)-PLATIN(ILA1G)
    ZDY32 = PLATIN(ILA3G)-PLATIN(ILA2G)

    ! * L - polynom in x- direction; parallel of point 1

    ZDLO1 = PDLO(1,1,JLEN)
    ZWXN3 = ((1.0_JPRB+ZDLO1)*(ZDLO1)*(ZDLO1-1.0_JPRB))/6._JPRB
    ZWXN2 = - ((1.0_JPRB+ZDLO1)*(ZDLO1)*(ZDLO1-2.0_JPRB))/2.0_JPRB
    ZWXN1 = ((1.0_JPRB+ZDLO1)*(ZDLO1-1.0_JPRB)*(ZDLO1-2.0_JPRB))/2.0_JPRB
    ZWXN0 = 1.0_JPRB-ZWXN1-ZWXN2-ZWXN3

    ! * L - polynom in x- direction; parallel of point 3

    ZDLO2 = PDLO(2,1,JLEN)
    ZWXS3 = ((1.0_JPRB+ZDLO2)*(ZDLO2)*(ZDLO2-1.0_JPRB))/6._JPRB
    ZWXS2 = - ((1.0_JPRB+ZDLO2)*(ZDLO2)*(ZDLO2-2.0_JPRB))/2.0_JPRB
    ZWXS1 = ((1.0_JPRB+ZDLO2)*(ZDLO2-1.0_JPRB)*(ZDLO2-2.0_JPRB))/2.0_JPRB
    ZWXS0 = 1.0_JPRB-ZWXS1-ZWXS2-ZWXS3

    ! * L - polynom in y- direction

    ZWY3 = ((ZDY+ZDY10)*(ZDY)*(ZDY-ZDY21))&
     & /((ZDY10+ZDY21+ZDY32)*(ZDY21+ZDY32)*(ZDY32))  
    ZWY2 = ((ZDY+ZDY10)*(ZDY)*(ZDY-ZDY21-ZDY32))&
     & /((ZDY10+ZDY21)*(ZDY21)*(-ZDY32))  
    ZWY1 = ((ZDY+ZDY10)*(ZDY-ZDY21)*(ZDY-ZDY21-ZDY32))&
     & /((ZDY10)*(-ZDY21)*(-ZDY21-ZDY32))  
    ZWY0 = 1.0_JPRB-ZWY1-ZWY2-ZWY3

    ! * Linear parts for extreme rows.

    ZWXNN=PDLO(0,1,JLEN)
    ZCXNN=1.0_JPRB-PDLO(0,1,JLEN)

    ZWXSS=PDLO(3,1,JLEN)
    ZCXSS=1.0_JPRB-PDLO(3,1,JLEN)

    ! * Weights for 12 points interpolation.

    PWXX(JLEN,1)  = ZWXN1*ZWY1
    PWXX(JLEN,2)  = ZWXN2*ZWY1
    PWXX(JLEN,3)  = ZWXS1*ZWY2
    PWXX(JLEN,4)  = ZWXS2*ZWY2
    PWXX(JLEN,5)  = ZCXNN*ZWY0
    PWXX(JLEN,6)  = ZWXNN*ZWY0
    PWXX(JLEN,7)  = ZWXN0*ZWY1
    PWXX(JLEN,8)  = ZWXN3*ZWY1
    PWXX(JLEN,9)  = ZWXS0*ZWY2
    PWXX(JLEN,10) = ZWXS3*ZWY2
    PWXX(JLEN,11) = ZCXSS*ZWY3
    PWXX(JLEN,12) = ZWXSS*ZWY3

  ENDDO

ENDIF

!     ---------------------------------------------------------------------

!*       4. ADDITIONAL LINEAR WEIGHTS FOR MULTI-LINEAR INTERPOLATION.
!           ---------------------------------------------------------

IF (LDML) THEN
  ! * points number 1, 2, 3, 4 (and 5 to 12 if cubic interpolations):
  !   multiply by ZRW(.,1).
  DO JBINL = 1,KBINL
    DO JLEN = 1,KEND
      PWXX(JLEN,JBINL)=ZRW(JLEN,1)*PWXX(JLEN,JBINL)
    ENDDO
  ENDDO

  ! * points number 5, 6, 11, 12:
  DO JLEN = 1,KEND
    ZWY=PDLAT(2,JLEN)
    ZCY=1.0_JPRB-PDLAT(2,JLEN)

    ZWXN=PDLO(1,2,JLEN)
    ZCXN=1.0_JPRB-PDLO(1,2,JLEN)

    ZWXS=PDLO(2,2,JLEN)
    ZCXS=1.0_JPRB-PDLO(2,2,JLEN)

    PWXY(JLEN, 5)=ZRW(JLEN,2)*ZCXN*ZCY
    PWXY(JLEN, 6)=ZRW(JLEN,2)*ZWXN*ZCY
    PWXY(JLEN,11)=ZRW(JLEN,2)*ZCXS*ZWY
    PWXY(JLEN,12)=ZRW(JLEN,2)*ZWXS*ZWY
  ENDDO

  ! * points number 7, 8, 9, 10:
  DO JLEN = 1,KEND
    ZWY=PDLAT(3,JLEN)
    ZCY=1.0_JPRB-PDLAT(3,JLEN)

    ZWXN=PDLO(1,3,JLEN)
    ZCXN=1.0_JPRB-PDLO(1,3,JLEN)

    ZWXS=PDLO(2,3,JLEN)
    ZCXS=1.0_JPRB-PDLO(2,3,JLEN)

    PWXY(JLEN, 7)=ZRW(JLEN,3)*ZCXN*ZCY
    PWXY(JLEN, 8)=ZRW(JLEN,3)*ZWXN*ZCY
    PWXY(JLEN, 9)=ZRW(JLEN,3)*ZCXS*ZWY
    PWXY(JLEN,10)=ZRW(JLEN,3)*ZWXS*ZWY
  ENDDO

  ! * points number 13, 14, 15, 16:
  DO JLEN = 1,KEND
    ZWY=PDLAT(4,JLEN)
    ZCY=1.0_JPRB-PDLAT(4,JLEN)

    ZWXN=PDLO(1,4,JLEN)
    ZCXN=1.0_JPRB-PDLO(1,4,JLEN)

    ZWXS=PDLO(2,4,JLEN)
    ZCXS=1.0_JPRB-PDLO(2,4,JLEN)

    PWXY(JLEN,13)=ZRW(JLEN,4)*ZCXN*ZCY
    PWXY(JLEN,14)=ZRW(JLEN,4)*ZWXN*ZCY
    PWXY(JLEN,15)=ZRW(JLEN,4)*ZCXS*ZWY
    PWXY(JLEN,16)=ZRW(JLEN,4)*ZWXS*ZWY
  ENDDO
ELSE
  ! PWXY is not used later and filled with zeros.
  PWXY(:,:)=0._JPRB
ENDIF

!     ---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUHOW2',1,ZHOOK_HANDLE)
END SUBROUTINE SUHOW2
