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

!ocl list_copy(32,IADDR)
!ocl list_copy(32,PLSDEPI)
!ocl list_copy(32,PLATI)
!ocl list_copy(32,PIPI0)
!ocl list_copy(32,PIPI1)
!ocl list_copy(32,PIPI2)
SUBROUTINE RDSCAW(YDSL,KPROMB,KSTART,KPROF,&
 & KSTABUF,KWIS,KFRSTLOFF,&
 & P4JP,PIS2,PLSDEPI,PLATI,&
 & PIPI,PLON,PLAT,&
 & PDLAT,PCLA,PDLO,PCLO,KL0)

!**** *RDSCAW  -  Interpolator for radiation scheme:
!                 Storage of Coordinates And Weights.

!     Purpose.
!     --------
!       Same kind of action as LASCAW but for the interpolator used
!        in the radiation scheme.
!       Determines the interpolation grid:
!       - computation of the latitude and the longitude of the
!         point situated at the upper left corner of the 16 points
!         square, and of the interpolation point.
!       - computation of weights.
!       Storage of coordinates and weights.

!**   Interface.
!     ----------
!        *CALL* *RDSCAW( ... )

!        Explicit arguments :
!        --------------------

!        INPUT:
!          YDSL    - SL_STRUCT definition
!          KPROMB  - horizontal dimension for interpolation point quantities.
!          KSTART  - first element of arrays where computations are performed.
!          KPROF   - depth of work.
!          KSTABUF - for a latitude IGL, KSTABUF(IGL) is the
!                    address of the element corresponding to
!                    (ILON=0,IGL) in the KPROMB arrays.
!          KWIS    - kind of interpolation.
!          KFRSTLOFF   - offset for first lat of own a-set in grid-point space,
!                        see NFRSTLOFF
!          P4JP    - Approximative inverse of the differences of latitudes.
!          PIS2    - PI / 2
!          PLSDEPI - (Number of points by latitude) / (2 * PI) .
!          PLATI   - latitude on the computational sphere.
!          PIPI   - coefficients for the bicubic interpolations.
!          PLON    - Interpolation point longitude on the computational sphere.
!          PLAT    - Interpolation point latitude on the computational sphere.

!        OUTPUT:
!          PDLAT   - distance for horizontal linear interpolations
!                    in latitude
!          PCLA    - weights for horizontal cubic interpolations
!                    in latitude
!          PDLO    - distances for horizontal linear interpolations
!                    in longitude (latitude rows 0, 1, 2, 3)
!          PCLO    - weights for horizontal cubic interpolations in
!                    longitude (latitude rows 1, 2)
!          KL0     - index of the four western points
!                    of the 16 points interpolation grid.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!      See documentation about ECMWF physics, and semi-Lagrangian interpolator.

!     Externals.
!     ----------

!        No external.
!        Called by RADINTG

!     Reference.
!     ----------

!     Author.
!     -------
!        G.Mozdzynski (based on lascaw)

!     Modifications.
!     --------------
!        Original : November 2001
!        Modified 02-10-01 GMozdzynski support on-demand radiation comms
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        R. El Khatib 28-Jul-2006 Porting to NEC
!        F. Vana  21-Aug-2008 - precomputation of longitudal interpolation 
!                             - generalized support for memory conflicts
!        G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!        G.Mozdzynski (Aug 2011): support higher order interpolation
!        G. Mozdzynski (May 2012): further cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE EINT_MOD , ONLY : SL_STRUCT,JPDUP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(SL_STRUCT),   INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMB 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KWIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFRSTLOFF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P4JP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLATI(YDSL%NDGSAH:YDSL%NDGENH) 
REAL(KIND=JPRD)   ,INTENT(IN)    :: PIPI(YDSL%NDGSAH:YDSL%NDGENH,1:3) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLON(KPROMB) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT(KPROMB) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLAT(KPROMB) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLA(KPROMB,3) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLO(KPROMB,0:3) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLO(KPROMB,3,2)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KL0(KPROMB,0:3) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IADDR(YDSL%NDGSAH:YDSL%NDGENH)

INTEGER(KIND=JPIM) :: ILA, ILO, ILO1, ILO2, ILO3, IZLAT, JLAT, JROF

REAL(KIND=JPRB) :: PD, ZDA, ZDB, ZDC, ZDD, ZLO, ZLO1, ZLO2, ZLO3

REAL(KIND=JPRB),PARAMETER :: Z6_R=1.0_JPRB/6.0_JPRB

REAL(KIND=JPRB) :: FLAG1, FLAG2, FLAG3, ZEPS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! Duplicata of some dummy or local arrays for optimisation on NEC platform.
INTEGER(KIND=JPIM) :: J_, JK_
INTEGER(KIND=JPIM) :: IADDR_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB) :: ZLSDEPI_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB) :: ZRLATI_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRD) :: ZRIPI0_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRD) :: ZRIPI1_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRD) :: ZRIPI2_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)

!     ------------------------------------------------------------------

!     FUNCTIONS:
! weights for cubic Lagrange interpolation (regular nodes)
FLAG1(PD)= 0.5_JPRB*(PD+1.0_JPRB)   *(PD-1.0_JPRB)*(PD-2.0_JPRB)
FLAG2(PD)=-0.5_JPRB*(PD+1.0_JPRB)*PD              *(PD-2.0_JPRB)
FLAG3(PD)= Z6_R    *(PD+1.0_JPRB)*PD*(PD-1.0_JPRB)

!     ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RDSCAW',0,ZHOOK_HANDLE)

!     ----------------------------------------------------------------

!*       0.    PRELIMINARY INITIALISATIONS.
!              ----------------------------
ZEPS=100.0_JPRB*TINY(1.0_JPRB)

DO JLAT=YDSL%NDGSAH,YDSL%NDGENH
  IADDR(JLAT)=KSTABUF(JLAT)
ENDDO
DO J_ = 1, JPDUP
  DO JK_ = YDSL%NDGSAH, YDSL%NDGENH
    IADDR_(J_,JK_)=IADDR(JK_)
  ENDDO
ENDDO
!CDIR SHORTLOOP
DO J_ = 1, JPDUP
  ZLSDEPI_(J_,YDSL%NDGSAH:YDSL%NDGENH)=PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)
  ZRLATI_(J_,YDSL%NDGSAH:YDSL%NDGENH)=PLATI(YDSL%NDGSAH:YDSL%NDGENH)
  ZRIPI0_(J_,YDSL%NDGSAH:YDSL%NDGENH)=PIPI(YDSL%NDGSAH:YDSL%NDGENH,1)
  ZRIPI1_(J_,YDSL%NDGSAH:YDSL%NDGENH)=PIPI(YDSL%NDGSAH:YDSL%NDGENH,2)
  ZRIPI2_(J_,YDSL%NDGSAH:YDSL%NDGENH)=PIPI(YDSL%NDGSAH:YDSL%NDGENH,3)
ENDDO


!        Coordinates and weights for bilinear interpolations.

IF (KWIS == 201) THEN

!CDIR NODEP
  DO JROF=KSTART,KPROF

    J_ = MOD(JROF+1-KSTART,JPDUP)+1

    IZLAT =INT(P4JP*(PIS2-PLAT(JROF))+0.75_JPRB+ZEPS)-KFRSTLOFF
    ILA   =IZLAT+NINT(SIGN(0.5_JPRB,ZRLATI_(J_,IZLAT)-PLAT(JROF)+ZEPS)-1.5_JPRB)

    PDLAT(JROF)=(PLAT(JROF)-ZRLATI_(J_,ILA+1))&
     & /(ZRLATI_(J_,ILA+2)-ZRLATI_(J_,ILA+1)+ZEPS)  

    ZLO1  =PLON(JROF)*ZLSDEPI_(J_,ILA+1)
    ILO1  =INT(ZLO1)
    PDLO(JROF,1)=ZLO1-REAL(ILO1,JPRB)
    ZLO2  =PLON(JROF)*ZLSDEPI_(J_,ILA+2)
    ILO2  =INT(ZLO2)
    PDLO(JROF,2)=ZLO2-REAL(ILO2,JPRB)

    KL0(JROF,1)=IADDR_(J_,ILA+1)+YDSL%NSLEXT(ILO1,ILA+1)
    KL0(JROF,2)=IADDR_(J_,ILA+2)+YDSL%NSLEXT(ILO2,ILA+2)

    IF( YDSL%LSLONDEM.AND..NOT.YDSL%LSLONDEM_ACTIVE )THEN
      YDSL%MASK_SL2(KL0(JROF,1)  )=1
      YDSL%MASK_SL2(KL0(JROF,1)+1)=1
      YDSL%MASK_SL2(KL0(JROF,1)+2)=1
      YDSL%MASK_SL2(KL0(JROF,1)+3)=1
      YDSL%MASK_SL2(KL0(JROF,2)  )=1
      YDSL%MASK_SL2(KL0(JROF,2)+1)=1
      YDSL%MASK_SL2(KL0(JROF,2)+2)=1
      YDSL%MASK_SL2(KL0(JROF,2)+3)=1
    ENDIF

  ENDDO

ENDIF

!        Coordinates and weights for 12 points interpolations.

IF (KWIS == 203) THEN

!CDIR NODEP
  DO JROF=KSTART,KPROF

    J_ = MOD(JROF+1-KSTART,JPDUP)+1

    IZLAT =INT(P4JP*(PIS2-PLAT(JROF))+0.75_JPRB+ZEPS)-KFRSTLOFF
    ILA   =IZLAT+NINT(SIGN(0.5_JPRB,ZRLATI_(J_,IZLAT)-PLAT(JROF)+ZEPS)-1.5_JPRB)

    PDLAT(JROF)=(PLAT(JROF)-ZRLATI_(J_,ILA+1))&
     & /(ZRLATI_(J_,ILA+2)-ZRLATI_(J_,ILA+1)+ZEPS)  

    ZDA   =PLAT(JROF)-ZRLATI_(J_,ILA)
    ZDB   =PLAT(JROF)-ZRLATI_(J_,ILA+1)
    ZDC   =PLAT(JROF)-ZRLATI_(J_,ILA+2)
    ZDD   =PLAT(JROF)-ZRLATI_(J_,ILA+3)
    
    PCLA(JROF,1)=REAL((ZDA*ZDC)*ZDD,JPRD)*ZRIPI0_(J_,ILA+1)
    PCLA(JROF,2)=REAL((ZDA*ZDB)*ZDD,JPRD)*ZRIPI1_(J_,ILA+1)
    PCLA(JROF,3)=REAL((ZDA*ZDB)*ZDC,JPRD)*ZRIPI2_(J_,ILA+1)

    ZLO   =PLON(JROF)*ZLSDEPI_(J_,ILA  )
    ILO   =INT(ZLO )
    PDLO(JROF,0)=ZLO -REAL(ILO ,JPRB)
    ZLO1  =PLON(JROF)*ZLSDEPI_(J_,ILA+1)
    ILO1  =INT(ZLO1)
    PDLO(JROF,1)=ZLO1-REAL(ILO1,JPRB)
    ZLO2  =PLON(JROF)*ZLSDEPI_(J_,ILA+2)
    ILO2  =INT(ZLO2)
    PDLO(JROF,2)=ZLO2-REAL(ILO2,JPRB)
    ZLO3  =PLON(JROF)*ZLSDEPI_(J_,ILA+3)
    ILO3  =INT(ZLO3)
    PDLO(JROF,3)=ZLO3-REAL(ILO3,JPRB)

    PCLO(JROF,1,1)= FLAG1(PDLO(JROF,1))
    PCLO(JROF,2,1)= FLAG2(PDLO(JROF,1))
    PCLO(JROF,3,1)= FLAG3(PDLO(JROF,1))
    PCLO(JROF,1,2)= FLAG1(PDLO(JROF,2))
    PCLO(JROF,2,2)= FLAG2(PDLO(JROF,2))
    PCLO(JROF,3,2)= FLAG3(PDLO(JROF,2))

    KL0(JROF,0)=IADDR_(J_,ILA  )+YDSL%NSLEXT(ILO ,ILA  )
    KL0(JROF,1)=IADDR_(J_,ILA+1)+YDSL%NSLEXT(ILO1,ILA+1)
    KL0(JROF,2)=IADDR_(J_,ILA+2)+YDSL%NSLEXT(ILO2,ILA+2)
    KL0(JROF,3)=IADDR_(J_,ILA+3)+YDSL%NSLEXT(ILO3,ILA+3)

    IF( YDSL%LSLONDEM.AND..NOT.YDSL%LSLONDEM_ACTIVE )THEN
      YDSL%MASK_SL2(KL0(JROF,0)  )=1
      YDSL%MASK_SL2(KL0(JROF,0)+1)=1
      YDSL%MASK_SL2(KL0(JROF,0)+2)=1
      YDSL%MASK_SL2(KL0(JROF,0)+3)=1
      YDSL%MASK_SL2(KL0(JROF,1)  )=1
      YDSL%MASK_SL2(KL0(JROF,1)+1)=1
      YDSL%MASK_SL2(KL0(JROF,1)+2)=1
      YDSL%MASK_SL2(KL0(JROF,1)+3)=1
      YDSL%MASK_SL2(KL0(JROF,2)  )=1
      YDSL%MASK_SL2(KL0(JROF,2)+1)=1
      YDSL%MASK_SL2(KL0(JROF,2)+2)=1
      YDSL%MASK_SL2(KL0(JROF,2)+3)=1
      YDSL%MASK_SL2(KL0(JROF,3)  )=1
      YDSL%MASK_SL2(KL0(JROF,3)+1)=1
      YDSL%MASK_SL2(KL0(JROF,3)+2)=1
      YDSL%MASK_SL2(KL0(JROF,3)+3)=1
    ENDIF

  ENDDO

ENDIF

!     ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RDSCAW',1,ZHOOK_HANDLE)
END SUBROUTINE RDSCAW
