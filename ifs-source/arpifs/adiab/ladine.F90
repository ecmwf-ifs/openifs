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

SUBROUTINE LADINE(YDGEOMETRY,YDRIP,YDDYN,YDDYNA,KSTART,KPROF,YDSL,YDSLINT,KPROMA,KFLDN,KFLDX,&
 & KIBL,&
 & PURL0,PVRL0,PURL,PVRL,&
 & PUL0,PVL0,PCL0,PUL9,PVL9,PCL9,&
 & PUSI,PVSI,&
 & PUT1,PVT1,PSPT1 )  

!**** *LADINE* - semi-LAgrangian scheme:
!                Interface subroutine for interpolations and lagrangian
!                trends. (Programme INterface d'Ensemble).
!                2D Shallow water.

!     Purpose.
!     --------
!           Grid point calculations in dynamics.

!           Computes the medium point M and origin point O.
!           Does interpolations at O.
!           Computes lagrangian trends.

!**   Interface.
!     ----------
!        *CALL* *LADINE(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KSTART  - first element of arrays where
!                    computations are performed.
!          KPROF   - depth of work.
!          YDSL    - SL_STRUCT definition.
!          KPROMA  - horizontal dimension.
!          KFLDN   - number of the first field for semi-lag variables.
!          KFLDX   - number of the last field for semi-lag variables.
!          KIBL    - index into YRGSGEOM/YRCSGEOM instances in YDGEOMETRY
!          PURL0   - U-wind, to be interpolated.
!          PVRL0   - V-wind, to be interpolated.
!          PURL    - U-wind for research of medium point.
!          PVRL    - V-wind for research of medium point.
!          PUL0    - Quantity to interpolate linearly at O in the
!                    U-wind equation if NWLAG=3.
!          PVL0    - Quantity to interpolate linearly at O in the
!                    V-wind equation if NWLAG=3.
!          PCL0    - Quantity to interpolate linearly at O in the
!                    continuity equation if NVLAG=3.
!          PUL9    - Quantity to interpolate at O in the U-wind equation.
!          PVL9    - Quantity to interpolate at O in the V-wind equation.
!          PCL9    - Quantity to interpolate at O in the continuity equation.
!          PUSI    - Part of the t semi-implicit term evaluated at the
!                    final point in the U-wind equation.
!          PVSI    - Part of the t semi-implicit term evaluated at the
!                    final point in the V-wind equation.

!        INPUT AND OUTPUT:
!          PUT1    - t+dt U-wind.
!          PVT1    - t+dt V-wind.
!          PSPT1   - t+dt equivalent height.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!        K. YESSAD after the subroutine LAG2 written by
!        Maurice IMBARD   METEO FRANCE/EERM/CRMD
!        Original : FEBRUARY 1992

!     Modifications.
!     --------------
!   Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!   R. El Khatib : 01-08-07 Pruning options
!   M.Hamrud      01-Oct-2003 CY28 Cleaning
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   J. Masek  30-Jun-2008  Dataflow for new SLHD interpolators.
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   K. Yessad (Jan 2011): LARCHE2 instead of LARCHE; remove LRPLANE features.
!   G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!   G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM, TCSGEOM and TCSGLEG
!   G. Mozdzynski (May 2012): further cleaning
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN       , ONLY : NULOUT
USE YOMCT0       , ONLY : LRPLANE
USE YOMDYNA      , ONLY : TDYNA
USE YOMDYN       , ONLY : TDYN
USE YOMSLINT     , ONLY : TSLINT
USE YOMCST       , ONLY : RPI, ROMEGA, RA
USE YOMRIP       , ONLY : TRIP
USE EINT_MOD     , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
TYPE(TSLINT)      ,INTENT(IN)    :: YDSLINT
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(IN)    :: PURL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PURL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUL9(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVL9(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCL9(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUT1(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVT1(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPT1(KPROMA) 

!     ------------------------------------------------------------------

! - trigonometric arrays.
REAL(KIND=JPRB) :: ZSINCO(KPROMA)
REAL(KIND=JPRB) :: ZCOSCO(KPROMA)
REAL(KIND=JPRB) :: ZSINLA(KPROMA)
REAL(KIND=JPRB) :: ZCOPHI(KPROMA)
REAL(KIND=JPRB) :: ZP    (KPROMA)
REAL(KIND=JPRB) :: ZQ    (KPROMA)
REAL(KIND=JPRB) :: ZR    (KPROMA)
REAL(KIND=JPRB) :: ZS    (KPROMA)
REAL(KIND=JPRB) :: ZLON(KPROMA)
REAL(KIND=JPRB) :: ZLAT(KPROMA)
REAL(KIND=JPRB) :: ZLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)
! - interpolated quantities arrays.
REAL(KIND=JPRB) :: ZU9(KPROMA)
REAL(KIND=JPRB) :: ZV9(KPROMA)
REAL(KIND=JPRB) :: ZC9(KPROMA)
REAL(KIND=JPRB) :: ZZUZ9(KPROMA)         ! To allow YRDYN%L2TLFF option
REAL(KIND=JPRB) :: ZZVZ9(KPROMA)         ! To allow YRDYN%L2TLFF option

INTEGER(KIND=JPIM) :: ISTALAT(YDSL%NDGSAH:YDSL%NDGENH)

INTEGER(KIND=JPIM) :: IGLGLO, JGL, JROF

LOGICAL :: LLSETTLST

REAL(KIND=JPRB) :: ZCOPH, ZDEPI,&
 & ZDSTRET, ZPHI, ZPIS2, ZPP, ZPP1, ZPP9, ZQQ,&
 & ZQQ1, ZQQ9, ZSINX, ZUT1, ZVT1  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "lacone.intfb.h"
#include "lainor2.intfb.h"
#include "larche2.intfb.h"
#include "larmes2.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LADINE',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL),YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(L2TLFF=>YDDYN%L2TLFF, &
 & LADVF=>YDDYN%LADVF, LQMHP=>YDDYN%LQMHP, LQMHW=>YDDYN%LQMHW, &
 & NITMP=>YDDYN%NITMP, NVLAG=>YDDYN%NVLAG, NWLAG=>YDDYN%NWLAG, &
 & VMAX1=>YDDYN%VMAX1, VMAX2=>YDDYN%VMAX2, NSTTYP=>YDGEOMETRY%YRGEM%NSTTYP, &
 & R4JP=>YDGEOMETRY%YRGEM%R4JP, RC2M1=>YDGEOMETRY%YRGEM%RC2M1, RC2P1=>YDGEOMETRY%YRGEM%RC2P1, &
 & RLOCEN=>YDGEOMETRY%YRGEM%RLOCEN, RMUCEN=>YDGEOMETRY%YRGEM%RMUCEN, RSTRET=>YDGEOMETRY%YRGEM%RSTRET, &
 & RDTS22=>YDRIP%RDTS22, &
 & RDTS62=>YDRIP%RDTS62, RDTSA=>YDRIP%RDTSA, RTDT=>YDRIP%RTDT, &
 & YDGEM=>YDGEOMETRY%YRGEM)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS:
!              ----------------------------

ZDEPI  =2.0_JPRB*RPI
ZPIS2  =0.5_JPRB*RPI

ZDSTRET=2.0_JPRB*RSTRET
DO JGL=MAX(YDSL%NDGSAG,YDSL%NDGSAL+YDSL%NFRSTLOFF-YDSL%NSLWIDE)-YDSL%NFRSTLOFF,&
   & MIN(YDSL%NDGENG,YDSL%NDGENL+YDSL%NFRSTLOFF+YDSL%NSLWIDE)-YDSL%NFRSTLOFF  
  IGLGLO=JGL+YDSL%NFRSTLOFF
  ZLSDEPI(JGL)=REAL(YDSL%NLOENG(IGLGLO),JPRB)/ZDEPI
ENDDO

DO JGL=MAX(YDSL%NDGSAH,LBOUND(YDSL%NSLOFF,1)),MIN(YDSL%NDGENH,UBOUND(YDSL%NSLOFF,1))
  ISTALAT(JGL)=YDSL%NSLOFF(JGL)
ENDDO

!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF THE MEDIUM POINT AND
!              INTERPOLATIONS AT THIS POINT.
!              ------------------------------------

LLSETTLST=YDDYNA%LSETTLST.AND.(.NOT.YDDYNA%LELTRA)
IF (LRPLANE) THEN
  CALL ABOR1('LADINE: ELARMES2 NOT YET CODED ')
ELSE
  CALL LARMES2(YDGEOMETRY,YDDYN,YDDYNA,YDSL,YDSLINT,KPROMA,KSTART,KPROF,KFLDN,KFLDX,&
   & NULOUT,ISTALAT,&
   & NVLAG,NWLAG,NITMP,VMAX1,VMAX2,&
   & LRPLANE,LLSETTLST,YDDYNA%LELTRA,&
   & PURL0,PVRL0,PURL,PVRL,&
   & NSTTYP,ZDSTRET,RC2M1,RC2P1,R4JP,RPI,ZDEPI,ZPIS2,&
   & RDTSA,RDTS62,RDTS22,&
   & RLOCEN,RMUCEN,ZLSDEPI,YDCSGLEG%RLATI(YDSL%NDGSAH:),&
   & KIBL,&
   & ZLON,ZLAT,&
   & ZCOSCO,ZSINCO,ZSINLA,ZCOPHI,&
   & ZR,ZS)
ENDIF

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF THE ORIGIN POINT AND
!              INTERPOLATIONS AT THIS POINT.
!              ------------------------------------

CALL LAINOR2(YDGEOMETRY,YDDYN,YDDYNA,YDSL,YDSLINT,KPROMA,KSTART,KPROF,KFLDN,KFLDX,&
 & ISTALAT,&
 & NVLAG,NWLAG,&
 & LRPLANE,LQMHW,LQMHP,&
 & PUL9,PVL9,PCL9,PUL0,PVL0,PCL0,&
 & KIBL,&
 & NSTTYP,ZDSTRET,RC2M1,RC2P1,R4JP,RPI,ZDEPI,ZPIS2,&
 & RLOCEN,RMUCEN,ZLSDEPI,YDCSGLEG%RLATI(YDSL%NDGSAH:),&
 & ZCOSCO,ZSINCO,ZSINLA,ZCOPHI,&
 & ZLON,ZLAT,&
 & ZP,ZQ,&
 & ZU9,ZV9,ZC9,&
 & ZZUZ9,ZZVZ9)

!     ------------------------------------------------------------------

!*       4.    COMPUTATION OF LAGRANGIAN TRENDS.
!              ---------------------------------

!*       4.1   MOMENTUM EQUATION.

DO JROF=KSTART,KPROF
  ZPP=ZP(JROF)*YDGSGEOM%GNORDM(JROF)-ZQ(JROF)*YDGSGEOM%GNORDL(JROF)
  ZQQ=ZQ(JROF)*YDGSGEOM%GNORDM(JROF)+ZP(JROF)*YDGSGEOM%GNORDL(JROF)
  PUT1(JROF) = PUT1(JROF) + ZPP*ZU9(JROF) + ZQQ*ZV9(JROF)
  PVT1(JROF) = PVT1(JROF) - ZQQ*ZU9(JROF) + ZPP*ZV9(JROF)
ENDDO

IF (LADVF) THEN

  IF (LRPLANE) THEN
    ! not implemented in the 2D model.
  ELSE
    CALL LACONE(KPROMA,KSTART,KPROF,1,NSTTYP,&
     & ROMEGA,ZDSTRET,RC2M1,RC2P1,RMUCEN,&
     & RA,RPI,ZLON,ZLAT,ZU9,ZV9)  
  ENDIF

  DO JROF=KSTART,KPROF
    ZPP=ZP(JROF)*YDGSGEOM%GNORDM(JROF)-ZQ(JROF)*YDGSGEOM%GNORDL(JROF)
    ZQQ=ZQ(JROF)*YDGSGEOM%GNORDM(JROF)+ZP(JROF)*YDGSGEOM%GNORDL(JROF)
    PUT1(JROF) = PUT1(JROF) + ZPP*ZU9(JROF) + ZQQ*ZV9(JROF)
    PVT1(JROF) = PVT1(JROF) - ZQQ*ZU9(JROF) + ZPP*ZV9(JROF)
  ENDDO

ENDIF

!     Refined treatment of Coriolis term in 2TL scheme (or in 3TL scheme).
!     Works for the following configurations:
!      - 3TL SL scheme, LADVF=.T. or 2TL SL scheme, LADVF=.T.
!      - RW2TLFF=1.
!     Contrary to the 3D model where (PUSI, PVSI) are already substracted
!      in the not lagged calculations (because of t+dt physics), (PUSI,PVSI)
!      has to be substracted here.

IF (L2TLFF) THEN

  IF (LRPLANE) THEN
    ! not implemented in the 2D model.
  ELSE

    ! * compute new origin point:
    DO JROF=KSTART,KPROF
      ZPP=ZP(JROF)*YDGSGEOM%GNORDM(JROF)-ZQ(JROF)*YDGSGEOM%GNORDL(JROF)
      ZQQ=ZQ(JROF)*YDGSGEOM%GNORDM(JROF)+ZP(JROF)*YDGSGEOM%GNORDL(JROF)
      ZPP9=ZP(JROF)
      ZQQ9=ZQ(JROF)
      ZPP1=YDGSGEOM%GNORDM(JROF)
      ZQQ1=-YDGSGEOM%GNORDL(JROF)
      ! * new estimate of trajectory wind:
      ZUT1=0.5_JPRB*(&
       & ZPP1*(PUT1(JROF)-PUSI(JROF))+ZQQ1*(PVT1(JROF)-PVSI(JROF))&
       & +ZPP9*ZZUZ9(JROF)+ZQQ9*ZZVZ9(JROF))  
      ZVT1=0.5_JPRB*(&
       & -ZQQ1*(PUT1(JROF)-PUSI(JROF))+ZPP1*(PVT1(JROF)-PVSI(JROF))&
       & -ZQQ9*ZZUZ9(JROF)+ZPP9*ZZVZ9(JROF))  
      ! * subtract old estimate of advected Coriolis term:
      PUT1(JROF)=PUT1(JROF)-ZPP*ZU9(JROF)-ZQQ*ZV9(JROF)
      PVT1(JROF)=PVT1(JROF)+ZQQ*ZU9(JROF)-ZPP*ZV9(JROF)
      ! * compute new origin point:
      ZPHI=0.5_JPRB*RTDT*SQRT(ZUT1*ZUT1+ZVT1*ZVT1)/RA
      ZCOPH=1.0_JPRB-0.5_JPRB*ZPHI*ZPHI
      ZSINX=ZCOPH*RTDT*(1.0_JPRB-ZPHI*ZPHI/6.0_JPRB)/RA
      ZSINLA(JROF)=YDGSGEOM%GEMU(JROF)*(2.0_JPRB*ZCOPH*ZCOPH-1.0_JPRB)&
       & -ZVT1*YDGSGEOM%GSQM2(JROF)*ZSINX  
      ZCOSCO(JROF)=YDGSGEOM%GSQM2(JROF)*(2.0_JPRB*ZCOPH*ZCOPH-1.0_JPRB)&
       & +ZVT1*YDGSGEOM%GEMU(JROF)*ZSINX  
      ZSINCO(JROF)=-ZUT1*ZSINX
      ZCOPHI(JROF)=2.0_JPRB*ZCOPH*ZCOPH-1.0_JPRB
    ENDDO
    ! * compute rotation matrix and Coriolis term at the new origin point.
    CALL LARCHE2(KPROMA,KSTART,KPROF,NSTTYP,ZDSTRET,RC2M1,RC2P1,&
     & RPI,ZDEPI,RLOCEN,RMUCEN,YDGSGEOM,YDCSGEOM,&
     & ZCOSCO,ZSINCO,ZSINLA,ZCOPHI,1,ZLON,ZLAT,ZP,ZQ)  
    CALL LACONE(KPROMA,KSTART,KPROF,1,NSTTYP,&
     & ROMEGA,ZDSTRET,RC2M1,RC2P1,RMUCEN,&
     & RA,RPI,ZLON,ZLAT,ZU9,ZV9)  
    ! * add new estimate of advected Coriolis term:
    DO JROF=KSTART,KPROF
      ZPP=ZP(JROF)*YDGSGEOM%GNORDM(JROF)-ZQ(JROF)*YDGSGEOM%GNORDL(JROF)
      ZQQ=ZQ(JROF)*YDGSGEOM%GNORDM(JROF)+ZP(JROF)*YDGSGEOM%GNORDL(JROF)
      PUT1(JROF) = PUT1(JROF) + ZPP*ZU9(JROF) + ZQQ*ZV9(JROF)
      PVT1(JROF) = PVT1(JROF) - ZQQ*ZU9(JROF) + ZPP*ZV9(JROF)
    ENDDO

  ENDIF

ENDIF

!*       4.2   CONTINUITY EQUATION.

DO JROF=KSTART,KPROF
  PSPT1(JROF)=PSPT1(JROF)+ZC9(JROF)
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LADINE',1,ZHOOK_HANDLE)
END SUBROUTINE LADINE
