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

SUBROUTINE SUHOWLSM(KDGSA,KDGEN,KFRSTLOFF,KPROMA,KEND,KAFPB,PCR,PLATIN,PRGMSD,PGMSF,&
 & KLSIM,KBINL,LDML,KLA,KL0,PDLAT,PDLO,PWXX,PWXY,&
 & PTMERGL,PTSIN,PLSMIN,PTSOUT,PLSMOUT,LDPSL)

!**** *SUHOWLSM* - Compute the weights of nearest points.
!****
!****            Possibility to use the land/sea/ice mask

!     Purpose.
!     -------
!       Computes the weights for 12 points interpolations or
!       bilinear interpolations for FULL-POS fields.
!       Same as SUHOW2 but possibility of using the land/sea/ice mask.
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

!       For distributed memory SUHOWLSM is executed locally on each processor
!       (iproca,iprocb).
!       Abbreviation DM-global means 'variable global for all processors'.
!       Abbreviation DM-local means 'variable local for each processor
!        iproca and iprocb'.

!**   Interface.
!**   ---------
!**   CALL SUHOWLSM(...)
!**
!        Explicit arguments :
!        --------------------

!     INPUT:
!      KDGSA     : Lower DM-locally numbered latitude 'jgl' involved (DM-local).
!      KDGEN     : Upper DM-locally numbered latitude 'jgl' involved (DM-local).
!      KFRSTLOFF : Difference 'DM-global numbering of latitudes' -
!                  'DM-local numbering of latitudes'
!                  KFRSTLOFF=0 if not distributed memory (DM-local).
!      KPROMA    : Length of output vector, to which data are interpolated (DM-local)
!      KEND      : actual length of output vector
!      KAFPB     : horizontal dimension of a field in the buffer (DM-local).
!      PCR       : Coefficients to compute ZR.
!      PLATIN    : Gaussian/Lobatto latitudes.
!                  (DM-global, but only elements KDGSA+KFRSTLOFF to
!                  KDGEN+KFRSTLOFF are provided here).
!      PRGMSD    : 1/(departure local mesh-size) (DM-local).
!      PGMSF     : Target local mesh-size (DM-local).
!      PTMERGL   : Melting temperature of floating ice (DM-global).
!      PTSIN     : Interpolation buffer containing surface temperature (DM-local).
!      PLSMIN    : Interpolation buffer containing land/sea mask (DM-local).
!      PTSOUT    : Interpolated surface temperature (DM-local).
!      PLSMOUT   : Interpolated surface land/sea mask (DM-local).
!      KLSIM     : Options of land/sea/ice mask (DM-global):
!                  = 0 ==> no use of l/s/i mask
!                  = 1 ==> use of l/s mask
!                  = 2 ==> use of l/s/i mask
!                  = 3 ==> use of general masks
!      KBINL     : Options for choice of interpolations (DM-global):
!                  = 12 ==> 12 points interpolation
!                  =  4 ==> bilinear interpolation
!      LDML      : Multilinear interpolations (DM-global).
!      KLA       : Latitude index of the point numbered 1 (DM-local).
!      KL0       : index of the four western point (numbered 13, 7, 9, 15)
!                  of the 16 points interpolation grid (DM-local).
!      PDLAT     : Meridian linear weights (DM-local).
!      PDLO      : Zonal linear weights.
!      LDPSL     : For the 4-point interpolation weights calculation, which
!      points should be taken into account : LDPSL (i) == T implies that the
!      i-th point will be used to compute the weights

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
!     Author : K. YESSAD after routine SUHOW.
!     Date   : MARCH 1997

!     Modifications
!     -------------
!      R. El Khatib 01-06-14 LFPLAKE => NFPLAKE
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad (oct 2010): multi-linear interpolations.
!      P. Marguinaud (oct 2013): Add LDPSL, make some arguments optional,
!      cleaning using inlined functions, KLSIM==3 for general masks
!     ---------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE PARFPOS   , ONLY : JPROW

!     ---------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KDGSA
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGEN
INTEGER(KIND=JPIM),INTENT(IN)    :: KFRSTLOFF
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
INTEGER(KIND=JPIM),INTENT(IN)    :: KAFPB
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCR(2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLATIN(KDGSA+KFRSTLOFF:KDGEN+KFRSTLOFF)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRGMSD(KEND)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMSF(KEND)
INTEGER(KIND=JPIM),INTENT(IN)    :: KLSIM
INTEGER(KIND=JPIM),INTENT(IN)    :: KBINL
LOGICAL           ,INTENT(IN)    :: LDML
INTEGER(KIND=JPIM),INTENT(IN)    :: KLA(KEND)
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(KPROMA,2*JPROW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLAT(4,KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLO(0:3,4,KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWXX(KPROMA,KBINL)
REAL(KIND=JPRB), OPTIONAL  ,INTENT(OUT) :: PWXY(KPROMA,5:16)
REAL(KIND=JPRB), OPTIONAL  ,INTENT(IN)  :: PTMERGL
REAL(KIND=JPRB), OPTIONAL  ,INTENT(IN)  :: PTSIN(KAFPB)
REAL(KIND=JPRB), OPTIONAL  ,INTENT(IN)  :: PLSMIN(KAFPB)
REAL(KIND=JPRB), OPTIONAL  ,INTENT(IN)  :: PTSOUT(KPROMA)
REAL(KIND=JPRB), OPTIONAL  ,INTENT(IN)  :: PLSMOUT(KPROMA)
LOGICAL,         OPTIONAL  ,INTENT(IN)  :: LDPSL (4)

!     ---------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ILA1G, ILA2G, ILA3G, ILAG, JLEN, JBINL
INTEGER(KIND=JPIM) :: IMKOUT
INTEGER(KIND=JPIM) :: IMKIN1, IMKIN2, IMKIN3, IMKIN4, IMKIN5, IMKIN6, &
                    & IMKIN7, IMKIN8, IMKIN9, IMKIN10, IMKIN11, IMKIN12
INTEGER(KIND=JPIM) :: JR
INTEGER(KIND=JPIM) :: ILOL, ILOR, ILAN, ILAS, INW, INE, ISW, ISE
REAL(KIND=JPRB) :: ZCXN, ZCXNN, ZCXS, ZCXSS, ZCY, ZDLAT, &
 & ZDLO1, ZDLO2, ZDY, ZDY10, ZDY21, ZDY32, &
 & ZEPSIL, ZWLSI1, ZWLSI10, ZWLSI11, ZWLSI12, &
 & ZWLSI13, ZWLSI14, ZWLSI15, ZWLSI16, ZWLSI17, &
 & ZWLSI2, ZWLSI3, ZWLSI4, ZWLSI5, ZWLSI6, ZWLSI7, &
 & ZWLSI8, ZWLSI9, ZWXN, ZWXN0, ZWXN1, ZWXN2, &
 & ZWXN3, ZWXNN, ZWXS, ZWXS0, ZWXS1, ZWXS2, &
 & ZWXS3, ZWXSS, ZWY, ZWY0, ZWY1, ZWY2, ZWY3  
REAL(KIND=JPRB) :: ZR(KPROMA), ZRW(KPROMA,4)
REAL(KIND=JPRB) :: ZSUM
LOGICAL :: LLDI

#include "suhowlsm.decl.h"

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ---------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUHOWLSM',0,ZHOOK_HANDLE)
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

ZEPSIL=1.E-5_JPRB

!     ---------------------------------------------------------------------

!*       2. WEIGHTS FOR BILINEAR INTERPOLATION.
!           -----------------------------------

IF (KBINL == 4) THEN

  DO JLEN = 1,KEND

    ! * Preliminary initialisation of ZWLSI1 to ZWLSI4, ZWLSI14, ZWLSI15, ZWLSI17.

    ZWLSI1=1._JPRB
    ZWLSI2=1._JPRB
    ZWLSI3=1._JPRB
    ZWLSI4=1._JPRB

    ZWLSI14=0._JPRB
    ZWLSI15=0._JPRB
    ZWLSI17=0._JPRB

    ! * Modification of ZWLSI1 to ZWLSI4, ZWLSI14, ZWLSI15,
    !   ZWLSI17 according to value of the land/sea/ice mask.

    IF ( KLSIM  /=  0 ) THEN

      IMKOUT = JMKOUT ()

      IF(JMKIN (1) /= IMKOUT) ZWLSI1=0.0_JPRB
      IF(JMKIN (2) /= IMKOUT) ZWLSI2=0.0_JPRB
      IF(JMKIN (3) /= IMKOUT) ZWLSI3=0.0_JPRB
      IF(JMKIN (4) /= IMKOUT) ZWLSI4=0.0_JPRB

    ENDIF

    LLDI = (KLSIM /= 3) .OR. (ZWLSI1 + ZWLSI2 + ZWLSI3 + ZWLSI4) > 0._JPRB
  
    IF (.NOT. LLDI) THEN
      PWXX (JLEN,:) = 0._JPRB
      CYCLE
    ENDIF

    ZWLSI14=MIN(ZWLSI1+ZWLSI2,1.0_JPRB)
    ZWLSI15=MIN(ZWLSI3+ZWLSI4,1.0_JPRB)
    ZWLSI17=MIN(ZWLSI14+ZWLSI15,1.0_JPRB)

    ! * Definite values for ZWLSI1 to ZWLSI4, ZWLSI14 and ZWLSI15.

    ZWLSI1=1.0_JPRB-INT(ZWLSI14+ZEPSIL)*(1.0_JPRB-ZWLSI1)
    ZWLSI2=1.0_JPRB-INT(ZWLSI14+ZEPSIL)*(1.0_JPRB-ZWLSI2)
    ZWLSI3=1.0_JPRB-INT(ZWLSI15+ZEPSIL)*(1.0_JPRB-ZWLSI3)
    ZWLSI4=1.0_JPRB-INT(ZWLSI15+ZEPSIL)*(1.0_JPRB-ZWLSI4)

    ZWLSI14=1.0_JPRB-INT(ZWLSI17+ZEPSIL)*(1.0_JPRB-ZWLSI14)
    ZWLSI15=1.0_JPRB-INT(ZWLSI17+ZEPSIL)*(1.0_JPRB-ZWLSI15)

    ! * Weights for bilinear interpolation.

    ZWY=ZWLSI15*(1.0_JPRB+ZWLSI14*(PDLAT(1,JLEN)-1.0_JPRB))
    ZCY=1.0_JPRB-ZWY

    ZWXN=ZWLSI2*(1.0_JPRB+ZWLSI1*(PDLO(1,1,JLEN)-1.0_JPRB))
    ZCXN=1.0_JPRB-ZWXN

    ZWXS=ZWLSI4*(1.0_JPRB+ZWLSI3*(PDLO(2,1,JLEN)-1.0_JPRB))
    ZCXS=1.0_JPRB-ZWXS

    PWXX(JLEN,1)  = ZCXN*ZCY
    PWXX(JLEN,2)  = ZWXN*ZCY
    PWXX(JLEN,3)  = ZCXS*ZWY
    PWXX(JLEN,4)  = ZWXS*ZWY

    ! Force some weights to zero; they should already be zero, but sometimes
    ! they are 1E-16

    IF (KLSIM == 3) THEN
      IF(IMKIN1 /= IMKOUT) PWXX(JLEN,1)=0.0_JPRB
      IF(IMKIN2 /= IMKOUT) PWXX(JLEN,2)=0.0_JPRB
      IF(IMKIN3 /= IMKOUT) PWXX(JLEN,3)=0.0_JPRB
      IF(IMKIN4 /= IMKOUT) PWXX(JLEN,4)=0.0_JPRB
      ZSUM = SUM (PWXX(JLEN,1:4))
      IF (ZSUM > 0._JPRB) THEN
        PWXX(JLEN,1:4) = PWXX(JLEN,1:4) / ZSUM
      ENDIF
    ENDIF

  ENDDO

ENDIF

!     ---------------------------------------------------------------------

!*       3. WEIGHTS FOR 12 POINTS INTERPOLATION.
!           ------------------------------------

IF (KBINL == 12) THEN

  DO JLEN = 1,KEND

    ! * Preliminary initialisation of ZWLSI1 to ZWLSI17.

    ZWLSI1=1._JPRB
    ZWLSI2=1._JPRB
    ZWLSI3=1._JPRB
    ZWLSI4=1._JPRB
    ZWLSI5=1._JPRB
    ZWLSI6=1._JPRB
    ZWLSI7=1._JPRB
    ZWLSI8=1._JPRB
    ZWLSI9=1._JPRB
    ZWLSI10=1._JPRB
    ZWLSI11=1._JPRB
    ZWLSI12=1._JPRB

    ZWLSI13=0._JPRB
    ZWLSI14=0._JPRB
    ZWLSI15=0._JPRB
    ZWLSI16=0._JPRB
    ZWLSI17=0._JPRB

    ! * Modification of ZWLSI1 to ZWLSI17
    !   according to value of the land/sea/ice mask.

    IF ( KLSIM  /=  0 ) THEN

      IMKOUT = JMKOUT ()

      IF(JMKIN ( 1) /= IMKOUT) ZWLSI1=0.0_JPRB
      IF(JMKIN ( 2) /= IMKOUT) ZWLSI2=0.0_JPRB
      IF(JMKIN ( 3) /= IMKOUT) ZWLSI3=0.0_JPRB
      IF(JMKIN ( 4) /= IMKOUT) ZWLSI4=0.0_JPRB
      IF(JMKIN ( 5) /= IMKOUT) ZWLSI5=0.0_JPRB
      IF(JMKIN ( 6) /= IMKOUT) ZWLSI6=0.0_JPRB
      IF(JMKIN ( 7) /= IMKOUT) ZWLSI7=0.0_JPRB
      IF(JMKIN ( 8) /= IMKOUT) ZWLSI8=0.0_JPRB
      IF(JMKIN ( 9) /= IMKOUT) ZWLSI9=0.0_JPRB
      IF(JMKIN (10) /= IMKOUT) ZWLSI10=0.0_JPRB
      IF(JMKIN (11) /= IMKOUT) ZWLSI11=0.0_JPRB
      IF(JMKIN (12) /= IMKOUT) ZWLSI12=0.0_JPRB

    ENDIF

    LLDI = (KLSIM /= 3) .OR.  (ZWLSI1 + ZWLSI2 + ZWLSI3 + ZWLSI4 + ZWLSI5 + ZWLSI6 &
                           & + ZWLSI7 + ZWLSI8 + ZWLSI9 + ZWLSI10 + ZWLSI11 + ZWLSI12) > 0._JPRB

    IF (.NOT. LLDI) THEN
      PWXX(JLEN,:)  = 0._JPRB
      CYCLE
    ENDIF

    ZWLSI13=MIN(ZWLSI5+ZWLSI6,1.0_JPRB)
    ZWLSI14=MIN(ZWLSI7+ZWLSI1+ZWLSI2+ZWLSI8,1.0_JPRB)
    ZWLSI15=MIN(ZWLSI9+ZWLSI3+ZWLSI4+ZWLSI10,1.0_JPRB)
    ZWLSI16=MIN(ZWLSI11+ZWLSI12,1.0_JPRB)
    ZWLSI17=MIN(ZWLSI13+ZWLSI14+ZWLSI15+ZWLSI16,1.0_JPRB)

    ! * Definite values for ZWLSI1 to ZWLSI16.

    ZWLSI1=1.0_JPRB-INT(ZWLSI14+ZEPSIL)*(1.0_JPRB-ZWLSI1)
    ZWLSI2=1.0_JPRB-INT(ZWLSI14+ZEPSIL)*(1.0_JPRB-ZWLSI2)
    ZWLSI3=1.0_JPRB-INT(ZWLSI15+ZEPSIL)*(1.0_JPRB-ZWLSI3)
    ZWLSI4=1.0_JPRB-INT(ZWLSI15+ZEPSIL)*(1.0_JPRB-ZWLSI4)
    ZWLSI5=1.0_JPRB-INT(ZWLSI13+ZEPSIL)*(1.0_JPRB-ZWLSI5)
    ZWLSI6=1.0_JPRB-INT(ZWLSI13+ZEPSIL)*(1.0_JPRB-ZWLSI6)
    ZWLSI7=1.0_JPRB-INT(ZWLSI14+ZEPSIL)*(1.0_JPRB-ZWLSI7)
    ZWLSI8=1.0_JPRB-INT(ZWLSI14+ZEPSIL)*(1.0_JPRB-ZWLSI8)
    ZWLSI9=1.0_JPRB-INT(ZWLSI15+ZEPSIL)*(1.0_JPRB-ZWLSI9)
    ZWLSI10=1.0_JPRB-INT(ZWLSI15+ZEPSIL)*(1.0_JPRB-ZWLSI10)
    ZWLSI11=1.0_JPRB-INT(ZWLSI16+ZEPSIL)*(1.0_JPRB-ZWLSI11)
    ZWLSI12=1.0_JPRB-INT(ZWLSI16+ZEPSIL)*(1.0_JPRB-ZWLSI12)

    ZWLSI13=1.0_JPRB-INT(ZWLSI17+ZEPSIL)*(1.0_JPRB-ZWLSI13)
    ZWLSI14=1.0_JPRB-INT(ZWLSI17+ZEPSIL)*(1.0_JPRB-ZWLSI14)
    ZWLSI15=1.0_JPRB-INT(ZWLSI17+ZEPSIL)*(1.0_JPRB-ZWLSI15)
    ZWLSI16=1.0_JPRB-INT(ZWLSI17+ZEPSIL)*(1.0_JPRB-ZWLSI16)

    ! * Weights for 12 points interpolation.

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
    ZWXN3 = ZWLSI8*((1.0_JPRB+ZWLSI7*ZDLO1)*(1.0_JPRB+ZWLSI1*(ZDLO1-1.0_JPRB))&
     & *(1.0_JPRB+ZWLSI2*(ZDLO1-2.0_JPRB))) &
     & /((1.0_JPRB+2.0_JPRB*ZWLSI7)*(1.0_JPRB+ZWLSI1))  
    ZWXN2 = ZWLSI2*((1.0_JPRB+ZWLSI7*ZDLO1)*(1.0_JPRB+ZWLSI1*(ZDLO1-1.0_JPRB))&
     & *(1.0_JPRB+ZWLSI8*(ZDLO1-3._JPRB))) &
     & /((1.0_JPRB+ZWLSI7)*(1.0_JPRB-2.0_JPRB*ZWLSI8))  
    ZWXN1 = ZWLSI1*((1.0_JPRB+ZWLSI7*ZDLO1)*(1.0_JPRB+ZWLSI2*(ZDLO1-2.0_JPRB))&
     & *(1.0_JPRB+ZWLSI8*(ZDLO1-3._JPRB))) &
     & /((1.0_JPRB-2.0_JPRB*ZWLSI2)*(1.0_JPRB-3._JPRB*ZWLSI8))  
    ZWXN0 = 1.0_JPRB-ZWXN1-ZWXN2-ZWXN3

    ! * L - polynom in x- direction; parallel of point 3

    ZDLO2 = PDLO(2,1,JLEN)
    ZWXS3 = ZWLSI10*((1.0_JPRB+ZWLSI9*ZDLO2) &
     & *(1.0_JPRB+ZWLSI3*(ZDLO2-1.0_JPRB))&
     & *(1.0_JPRB+ZWLSI4*(ZDLO2-2.0_JPRB))) &
     & /((1.0_JPRB+2.0_JPRB*ZWLSI9)*(1.0_JPRB+ZWLSI3))  
    ZWXS2 = ZWLSI4*((1.0_JPRB+ZWLSI9*ZDLO2)*(1.0_JPRB+ZWLSI3*(ZDLO2-1.0_JPRB))&
     & *(1.0_JPRB+ZWLSI10*(ZDLO2-3._JPRB))) &
     & /((1.0_JPRB+ZWLSI9)*(1.0_JPRB-2.0_JPRB*ZWLSI10))  
    ZWXS1 = ZWLSI3*((1.0_JPRB+ZWLSI9*ZDLO2)*(1.0_JPRB+ZWLSI4*(ZDLO2-2.0_JPRB))&
     & *(1.0_JPRB+ZWLSI10*(ZDLO2-3._JPRB))) &
     & /((1.0_JPRB-2.0_JPRB*ZWLSI4)*(1.0_JPRB-3._JPRB*ZWLSI10))  
    ZWXS0 = 1.0_JPRB-ZWXS1-ZWXS2-ZWXS3

    ! * L - polynom in y- direction

    ZWY3 = ZWLSI16*((1.0_JPRB+ZWLSI13*(ZDY+ZDY10-1.0_JPRB)) &
     & *(1.0_JPRB+ZWLSI14*(ZDY-1.0_JPRB)) &
     & *(1.0_JPRB+ZWLSI15*(ZDY-ZDY21-1.0_JPRB)))/((1.0_JPRB+ZWLSI13 &
     & *(ZDY10+ZDY21+ZDY32-1.0_JPRB))*(1.0_JPRB+ZWLSI14 &
     & *(ZDY21+ZDY32-1.0_JPRB))*(1.0_JPRB+ZWLSI15*(ZDY32-1.0_JPRB)))  

    ZWY2 = ZWLSI15*((1.0_JPRB+ZWLSI13*(ZDY+ZDY10-1.0_JPRB)) &
     & *(1.0_JPRB+ZWLSI14*(ZDY-1.0_JPRB)) &
     & *(1.0_JPRB+ZWLSI16*(ZDY-ZDY21-ZDY32-1.0_JPRB)))/((1.0_JPRB+ZWLSI13 &
     & *(ZDY10+ZDY21-1.0_JPRB))*(1.0_JPRB+ZWLSI14*(ZDY21-1.0_JPRB)) &
     & *(1.0_JPRB+ZWLSI16*(-ZDY32-1.0_JPRB)))  

    ZWY1 = ZWLSI14*((1.0_JPRB+ZWLSI13*(ZDY+ZDY10-1.0_JPRB)) &
     & *(1.0_JPRB+ZWLSI15*(ZDY-ZDY21-1.0_JPRB))&
     & *(1.0_JPRB+ZWLSI16*(ZDY-ZDY21-ZDY32-1.0_JPRB)))/((1.0_JPRB+ZWLSI13 &
     & *(ZDY10-1.0_JPRB))*(1.0_JPRB+ZWLSI15*(-ZDY21-1.0_JPRB)) &
     & *(1.0_JPRB+ZWLSI16*(-ZDY21-ZDY32-1.0_JPRB)))  

    ZWY0 = 1.0_JPRB-ZWY1-ZWY2-ZWY3

    ! * Linear parts for extreme rows.

    ZWXNN=ZWLSI6*(1.0_JPRB+ZWLSI5*(PDLO(0,1,JLEN)-1.0_JPRB))
    ZCXNN=1.0_JPRB-ZWXNN

    ZWXSS=ZWLSI12*(1.0_JPRB+ZWLSI11*(PDLO(3,1,JLEN)-1.0_JPRB))
    ZCXSS=1.0_JPRB-ZWXSS

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

    ! Force some weights to zero; they should already be zero, but sometimes
    ! they are 1E-16

    IF (KLSIM == 3) THEN
      IF(IMKIN1  /= IMKOUT) PWXX(JLEN,1) =0.0_JPRB
      IF(IMKIN2  /= IMKOUT) PWXX(JLEN,2) =0.0_JPRB
      IF(IMKIN3  /= IMKOUT) PWXX(JLEN,3) =0.0_JPRB
      IF(IMKIN4  /= IMKOUT) PWXX(JLEN,4) =0.0_JPRB
      IF(IMKIN5  /= IMKOUT) PWXX(JLEN,5) =0.0_JPRB
      IF(IMKIN6  /= IMKOUT) PWXX(JLEN,6) =0.0_JPRB
      IF(IMKIN7  /= IMKOUT) PWXX(JLEN,7) =0.0_JPRB
      IF(IMKIN8  /= IMKOUT) PWXX(JLEN,8) =0.0_JPRB
      IF(IMKIN9  /= IMKOUT) PWXX(JLEN,9) =0.0_JPRB
      IF(IMKIN10 /= IMKOUT) PWXX(JLEN,10)=0.0_JPRB
      IF(IMKIN11 /= IMKOUT) PWXX(JLEN,11)=0.0_JPRB
      IF(IMKIN12 /= IMKOUT) PWXX(JLEN,12)=0.0_JPRB
      ZSUM = SUM (PWXX(JLEN,1:12))
      IF (ZSUM > 0._JPRB) THEN
        PWXX(JLEN,1:12) = PWXX(JLEN,1:12) / ZSUM
      ENDIF

    ENDIF

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

  ! * additional linear weights for points 5 to 16 (fill PWXY).
  DO JR=2,4

    SELECT CASE (JR)
    CASE(2)
      ! * points number 5, 6, 11, 12:
      ILOL=1
      ILOR=2
      ILAN=0
      ILAS=3
      INW=5
      INE=6
      ISW=11
      ISE=12
    CASE(3)
      ! * points number 7, 8, 9, 10:
      ILOL=0
      ILOR=3
      ILAN=1
      ILAS=2
      INW=7
      INE=8
      ISW=9
      ISE=10
    CASE(4)
      ! * points number 13, 14, 15, 16:
      ILOL=0
      ILOR=3
      ILAN=0
      ILAS=3
      INW=13
      INE=14
      ISW=15
      ISE=16
    END SELECT

    DO JLEN = 1,KEND

      ! * Preliminary initialisation of ZWLSI1 to ZWLSI4, ZWLSI14, ZWLSI15, ZWLSI17.

      ZWLSI1=1._JPRB
      ZWLSI2=1._JPRB
      ZWLSI3=1._JPRB
      ZWLSI4=1._JPRB

      ZWLSI14=0._JPRB
      ZWLSI15=0._JPRB
      ZWLSI17=0._JPRB

      ! * Modification of ZWLSI1 to ZWLSI4, ZWLSI14, ZWLSI15,
      !   ZWLSI17 according to value of the land/sea/ice mask.

      IF ( KLSIM  /=  0 ) THEN

        IMKOUT = JMKOUT ()

        IF(JMKIN (KIN1=ILAN, KIN2=ILOL) /= IMKOUT) ZWLSI1=0.0_JPRB
        IF(JMKIN (KIN1=ILAN, KIN2=ILOR) /= IMKOUT) ZWLSI2=0.0_JPRB
        IF(JMKIN (KIN1=ILAS, KIN2=ILOL) /= IMKOUT) ZWLSI3=0.0_JPRB
        IF(JMKIN (KIN1=ILAS, KIN2=ILOR) /= IMKOUT) ZWLSI4=0.0_JPRB

      ENDIF

      ZWLSI14=MIN(ZWLSI1+ZWLSI2,1.0_JPRB)
      ZWLSI15=MIN(ZWLSI3+ZWLSI4,1.0_JPRB)
      ZWLSI17=MIN(ZWLSI14+ZWLSI15,1.0_JPRB)

      ! * Definite values for ZWLSI1 to ZWLSI4, ZWLSI14 and ZWLSI15.

      ZWLSI1=1.0_JPRB-INT(ZWLSI14+ZEPSIL)*(1.0_JPRB-ZWLSI1)
      ZWLSI2=1.0_JPRB-INT(ZWLSI14+ZEPSIL)*(1.0_JPRB-ZWLSI2)
      ZWLSI3=1.0_JPRB-INT(ZWLSI15+ZEPSIL)*(1.0_JPRB-ZWLSI3)
      ZWLSI4=1.0_JPRB-INT(ZWLSI15+ZEPSIL)*(1.0_JPRB-ZWLSI4)

      ZWLSI14=1.0_JPRB-INT(ZWLSI17+ZEPSIL)*(1.0_JPRB-ZWLSI14)
      ZWLSI15=1.0_JPRB-INT(ZWLSI17+ZEPSIL)*(1.0_JPRB-ZWLSI15)

      ! * Weights for bilinear interpolation.

      ZWY=ZWLSI15*(1.0_JPRB+ZWLSI14*(PDLAT(JR,JLEN)-1.0_JPRB))
      ZCY=1.0_JPRB-ZWY

      ZWXN=ZWLSI2*(1.0_JPRB+ZWLSI1*(PDLO(1,JR,JLEN)-1.0_JPRB))
      ZCXN=1.0_JPRB-ZWXN

      ZWXS=ZWLSI4*(1.0_JPRB+ZWLSI3*(PDLO(2,JR,JLEN)-1.0_JPRB))
      ZCXS=1.0_JPRB-ZWXS

      PWXY(JLEN,INW)=ZRW(JLEN,JR)*ZCXN*ZCY
      PWXY(JLEN,INE)=ZRW(JLEN,JR)*ZWXN*ZCY
      PWXY(JLEN,ISW)=ZRW(JLEN,JR)*ZCXS*ZWY
      PWXY(JLEN,ISE)=ZRW(JLEN,JR)*ZWXS*ZWY
    ENDDO
  ENDDO

ELSEIF (PRESENT (PWXY)) THEN
  ! PWXY is not used later and filled with zeros.
  PWXY(:,:)=0._JPRB
ENDIF

!     ---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUHOWLSM',1,ZHOOK_HANDLE)

CONTAINS

#include "suhowlsm.func.h"

END SUBROUTINE SUHOWLSM

