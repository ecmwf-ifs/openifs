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

SUBROUTINE LARCIN2(YDGEOMETRY,YDDYN,YDDYNA,YDSL,YDSLINT,KPROMA,KSTART,KPROF,KFLDN,KFLDX,&
 & KSTABUF,&
 & KWISA,KXLAG,KROT,LDPLANE,&
 & LDQMHW,LDQMHP,&
 & KSTTYP,PDSTRET,PC2M1,PC2P1,P4JP,PI,PDEPI,PIS2,&
 & PLOCEN,PMUCEN,PLSDEPI,PLATI,&
 & KIBL,&
 & PCOSCO,PSINCO,PSINLA,PCOPHI,&
 & PURL0,PVRL0,&
 & PUSL,PVSL,PCSL,PUSL2,PVSL2,PCSL2,&
 & PLON,PLAT,PO,PQ,&
 & PUF0,PVF0,&
 & PUF,PVF,PCF,&
 & PUFZ,PVFZ)

!**** *LARCIN2 -  semi-LAgrangian scheme:
!                 Research of the Coordinates (of the medium or origin
!                 point) and INterpolations (2-d model).

!     Purpose.
!     --------
!       Computes the longitude and latitude of the interpolation
!       point from its cartesian coordinates.
!       Then computes the vector displacement matrix
!                        I po pq I
!                        I       I
!                        I-pq po I
!       from the interpolation point to the grid point.
!       At last determines the interpolation grid:
!       - computation of the latitude and the longitude of the
!         point situated at the upper left corner of the 16 points square,
!         and of the interpolation point.
!       - interpolations.

!**   Interface.
!     ----------
!        *CALL* *LARCIN2(......)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          YDSL    - SL_STRUCT definition
!          KPROMA  - horizontal dimension.
!          KSTART  - first element of arrays where
!                    computations are performed.
!          KPROF   - depth of work.
!          KFLDN   - number of the first field.
!          KFLDX   - number of the last field.
!          KSTABUF - for a latitude IGL, KSTABUF(IGL) is the
!                    address of the element corresponding to
!                    (ILON=1,IGL) in the KPROMA arrays.
!          KWISA   - switch :
!                    KWISA=1: interpolations for wind to find trajectory.
!                    KWISA=3: origin point interp.
!          KXLAG   - array containing switches concerning discretisations
!                    of equations (NVLAG,NTLAG,NWLAG,NSVDLAG,NSPDLAG).
!                    Only KXLAG(1) and KXLAG(4) are used in the 2D model.
!          KROT    - KROT=1: computation of the elements po and pq
!                    of the wind displacement matrix.
!                    KROT=0: no computation.
!          LDPLANE - switch: .T. = plane geometry; .F. = spherical geometry.
!          LDQMH[X]- Use quasi-monotone interpolation in the horizontal
!                    for field X (Wind, eq. height).
!          KSTTYP  - 1: Not tilted pole;  2: Tilted pole.
!          PDSTRET - 2*c (where c is the stretching factor).
!          PC2M1   - c*c-1.
!          PC2P1   - c*c+1.
!          P4JP    - Approximative inverse of the differences of latitudes.
!          PI      - number PI
!          PDEPI   - 2 * PI
!          PIS2    - PI / 2
!          PLOCEN  - geographical longitude of the stretching pole.
!          PMUCEN  - sine of the geographical latitude of the stretching pole.
!          PLSDEPI - (Number of points by latitude) / (2 * PI) .
!          PLATI   - latitude on the computational sphere.
!          KIBL    - index into YRGSGEOM/YRCSGEOM instances in YDGEOMETRY
!          PCOSCO  - cos(Longitude-Longitude(grid-point))*cos(Latitude)
!                    of the interpolation point (geographical longitude
!                    and latitude).
!                    ! if LDPLANE: x - coordinate (fractional system).
!          PSINCO  - sin(Longitude-Longitude(grid-point))*cos(Latitude)
!                    of the interpolation point (geographical longitude
!                    and latitude).
!                    ! if LDPLANE: y - coordinate (fractional system).
!          PSINLA  - sine of the interpolation point geographical latitude.
!          PCOPHI  - cosine of the geographical angle between the
!                    interpolation point and the grid-point.
!          PURL0   - U-component of the wind.
!          PVRL0   - V-component of the wind.
!          PUSL    - Quantity to interpolate at O in the U-wind equation.
!          PVSL    - Quantity to interpolate at O in the V-wind equation.
!          PCSL    - Quantity to interpolate at O in the continuity equation.
!          PUSL2   - Second quantity to interpolate at O
!                    (if NWLAG=3) in the U-component of the wind equation.
!          PVSL2   - Second quantity to interpolate at O
!                    (if NWLAG=3) in the V-component of the wind equation.
!          PCSL2   - Second quantity to interpolate at O
!                    (if NVLAG=3) in the continuity equation.

!        OUTPUT:
!          PLON    - computational sphere longitude of interpolation point.
!          PLAT    - computational sphere latitude of interpolation point.
!          PO      - first element of the wind displacement matrix.
!          PQ      - second element of the wind displacement matrix.
!          PUF0    - Interpolated U-wind at the medium point.
!          PVF0    - Interpolated V-wind at the medium point.
!          PUF     - Interpolated quantity at O in the U-wind equation.
!          PVF     - Interpolated quantity at O in the V-wind equation.
!          PCF     - Interpolated quantity at O in the continuity equation.
!          PUFZ    - Interpolated quantity at O in the U-wind equation
!                    in the SL2TL scheme for case L2TLFF=.T.
!          PVFZ    - Interpolated quantity at O in the V-wind equation
!                    in the SL2TL scheme for case L2TLFF=.T.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!      C. Temperton (ECMWF), by simplification of LARCIN
!      Original : 98/01/19

!     Modifications.
!     --------------
!   Modified 09-Jan-2007 by F. Vana: update of arguments for LASCAW
!   Modified 28-Aug-2007 by F. Vana: update of arguments for LASCAW
!   30-Jun-2008 J. Masek    Dataflow for new SLHD interpolators.
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   K. Yessad Nov 2008: interpolation routines: merge QM with not-QM version.
!   K. Yessad (Aug 2009): use RIPI, RSLD
!   K. Yessad (Jan 2011): LARCHE2 instead of LARCHE; remove LRPLANE features.
!   F. Vana  21-Feb-2011:  update of arguments for LASCAW
!   G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!   G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!   G. Mozdzynski (May 2012): further cleaning
!   F. Vana  13-Feb-2014  update of arguments for LASCAW
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   S. Malardel (Nov 2013): COMAD weights for SL interpolations
!   K. Yessad (March 2017): simplify level numbering in interpolator.
!   F. Vana    21-Nov-2017: Option LHOISLT
!   F. Vana  18-Jul-2019:  SLVF
! End Modifications
!-------------------------------------------------------------------------------

USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : TDYNA
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMSLINT     , ONLY : TSLINT
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE EINT_MOD     , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
TYPE(TSLINT)      ,INTENT(IN)    :: YDSLINT
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KWISA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KXLAG(6) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KROT 
LOGICAL           ,INTENT(IN)    :: LDPLANE 
LOGICAL           ,INTENT(IN)    :: LDQMHW 
LOGICAL           ,INTENT(IN)    :: LDQMHP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTTYP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDSTRET
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC2M1 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC2P1
REAL(KIND=JPRB)   ,INTENT(IN)    :: P4JP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDEPI 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLOCEN
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMUCEN
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLATI(YDSL%NDGSAH:YDSL%NDGENH) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOSCO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSINCO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSINLA(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOPHI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PURL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSL(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSL(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSL(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSL2(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSL2(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSL2(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLON(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLAT(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQ(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUF0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVF0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUF(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVF(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCF(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUFZ(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVFZ(KPROMA) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZDLAT(KPROMA)
REAL(KIND=JPRB) :: ZCLA (KPROMA,1,3)
REAL(KIND=JPRB) :: ZDLO (KPROMA,1,0:3)
REAL(KIND=JPRB) :: ZCLO (KPROMA,1,3,2)
REAL(KIND=JPRB) :: ZUF  (KPROMA)
REAL(KIND=JPRB) :: ZVF  (KPROMA)
REAL(KIND=JPRB) :: ZTF  (KPROMA)

INTEGER(KIND=JPIM) :: IL0 (KPROMA,1,0:3)
INTEGER(KIND=JPIM) :: ILH0(KPROMA,1,0:3)

INTEGER(KIND=JPIM) :: ISPLTHOI, IDIMK
INTEGER(KIND=JPIM) :: IHOR, IVLAG, IWIS, IWLAG, JROF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! unused arguments in call to LASCAW
! a) input - array dimensions omitted
INTEGER(KIND=JPIM) :: IUN_VAUTF(1)
REAL(KIND=JPRB) :: ZUN_SLD(1,3)
REAL(KIND=JPRB) :: ZUN_SLDW(1)
REAL(KIND=JPRB) :: ZUN_3DTW(1)
REAL(KIND=JPRB) :: ZUN_LEV(1)
REAL(KIND=JPRB) :: ZUN_VETAF(1)
REAL(KIND=JPRB) :: ZUN_VCUICO(1)
REAL(KIND=JPRB) :: ZUN_VSLD(1)
REAL(KIND=JPRB) :: ZUN_VSLDW(1)
REAL(KIND=JPRB) :: ZUN_VSLVF(1)
REAL(KIND=JPRB) :: ZUN_KAPPA(1),ZUN_KAPPAT(1),ZUN_KAPPAM(1),ZUN_KAPPAH(1)
REAL(KIND=JPRB) :: ZUN_MADU(1),ZUN_MADV(1),ZUN_MADW(1)
REAL(KIND=JPRB) :: ZUN_GAMMA_WENO(1)
! b) output - array dimensions kept for safety
INTEGER(KIND=JPIM) :: IUN_LEV(KPROMA,1)
INTEGER(KIND=JPIM) :: IUN_NOWENO(KPROMA,1)
REAL(KIND=JPRB) :: ZUN_DLAMAD(KPROMA)
REAL(KIND=JPRB) :: ZUN_CLASLD(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_CLAMAD(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_DLOMAD (KPROMA,1,0:3)
REAL(KIND=JPRB) :: ZUN_CLOSLD(KPROMA,1,3,2)
REAL(KIND=JPRB) :: ZUN_CLOMAD(KPROMA,1,3,2)
REAL(KIND=JPRB) :: ZUN_CLASLT(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_CLOSLT(KPROMA,1,3,2)
REAL(KIND=JPRB) :: ZUN_DVER(KPROMA,1)
REAL(KIND=JPRB) :: ZUN_DVERMAD(KPROMA,1)
REAL(KIND=JPRB) :: ZUN_VINTW(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_VINTWSLD(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_VINTWMAD(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_VINTWSLT(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_VINTWSLVF(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_VINTWS(KPROMA,1,1:4)
REAL(KIND=JPRB) :: ZUN_VDERW(KPROMA,1,0,0)  ! case KHVI=0
REAL(KIND=JPRB) :: ZUN_HVW(KPROMA,1,0)      ! case KHVI=0
REAL(KIND=JPRB) :: ZUN_CW(KPROMA,1,3)

LOGICAL         :: LLSLHD,LLSLHDQUAD,LLSLHD_OLD,LL3DTURB
LOGICAL         :: LLCOMAD, LLCOMADH, LLCOMADV
!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "laiddi.intfb.h"
#include "laidli.intfb.h"
#include "larche2.intfb.h"
#include "lascaw.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LARCIN2',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

! Switching off the SLHD/COMAD related weights computation (not coded for 2D model)
LLSLHD=.FALSE.
LLSLHDQUAD=.FALSE.
LLSLHD_OLD=.FALSE.
LLCOMAD=.FALSE.
LLCOMADH=.FALSE.
LLCOMADV=.FALSE.

ISPLTHOI=0
LL3DTURB=.FALSE.
IDIMK=1

IHOR=0

! * I[X]LAG:
IVLAG=KXLAG(1)
IWLAG=KXLAG(4)

! * Input variable IWIS for LASCAW.
IF (KWISA == 1) THEN
  ! * trajectory research.
  IF (YDDYNA%LHOISLT) THEN
    IWIS=202
  ELSE
    IWIS=201
  ENDIF
ELSEIF (KWISA == 3) THEN
  ! * origin point interpolations.
  IWIS=203
ELSE
  IWIS=-999
  CALL ABOR1('LARCIN2: WRONG VALUE FOR IWIS')
ENDIF

!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF LAT LON OF THE INTERPOLATION POINT.
!              IF KROT=1 COMPUTATION OF THE WIND DISPLACEMENT MATRIX
!              FROM THE INTERPOLATION POINT TO THE FINAL POINT.
!              ( T FOR LATITUDE THETA, L FOR LONGITUDE LAMBDA).
!              PO = ( 1 / (1+cos(PHI)) )
!                  *( cos(TG)*cos(T) + (1+sin(TG)*sin(T))*cos(L-LG) )
!              PQ = (-1 / (1+cos(PHI)) )
!                  *( sin(TG)+sin(T) ) * sin(L-LG)

!     ------------------------------------------------------------------

IF (LDPLANE) THEN
  ! not implemented in the 2D model
ELSE

  CALL LARCHE2(KPROMA,KSTART,KPROF,&
   & KSTTYP,PDSTRET,PC2M1,PC2P1,PI,PDEPI,&
   & PLOCEN,PMUCEN,YDGEOMETRY%YRGSGEOM(KIBL),YDGEOMETRY%YRCSGEOM(KIBL),&
   & PCOSCO,PSINCO,PSINLA,PCOPHI,&
   & KROT,PLON,PLAT,PO,PQ)  

  CALL LASCAW(&
   ! --- INPUT ----------------------------------------------------------------
   & YDSLINT%YRVSPLIP,YDSL,YDGEOMETRY%YRCVER%LREGETA,YDDYNA, &
   & KPROMA,IDIMK,KSTART,KPROF,1,&
   & KFLDN,KSTABUF,IWIS,IHOR,1,0,&
   & LLSLHD,LLSLHDQUAD,LLSLHD_OLD,YDDYN%LSLHDHEAT,LL3DTURB,&
   & LLCOMAD,LLCOMADH,LLCOMADV,ISPLTHOI,&
   & P4JP,PIS2,PLSDEPI,PLATI,&
   & YDSLINT%YRHSLMER%RIPI(YDSL%NDGSAH:YDSL%NDGENH,1:3),&
   & ZUN_SLD,ZUN_SLDW,ZUN_3DTW,&
   & PLON,PLAT,ZUN_LEV,&
   & ZUN_VETAF,IUN_VAUTF,&
   & ZUN_VCUICO,ZUN_VSLD,ZUN_VSLDW,ZUN_VSLVF,ZUN_GAMMA_WENO,1,1.0_JPRB,&
   & ZUN_KAPPA,ZUN_KAPPAT,ZUN_KAPPAM,ZUN_KAPPAH,&
   & ZUN_MADU,ZUN_MADV,ZUN_MADW,&
   ! --- OUTPUT ---------------------------------------------------------------
   & ZDLAT,ZUN_DLAMAD,ZCLA,ZUN_CLASLD,ZUN_CLASLT,ZUN_CLAMAD,&
   & ZDLO,ZUN_DLOMAD,ZCLO,ZUN_CLOSLD,ZUN_CLOSLT,ZUN_CLOMAD,&
   & IL0,ILH0,IUN_LEV,IUN_NOWENO,ZUN_CW,&
   & ZUN_DVER,ZUN_DVERMAD,ZUN_VINTW,ZUN_VINTWSLD,ZUN_VINTWSLT,&
   & ZUN_VINTWSLVF,ZUN_VINTWMAD,ZUN_VINTWS,&
   & ZUN_VDERW,ZUN_HVW)

ENDIF

!     ------------------------------------------------------------------

!*       3.    INTERPOLATIONS OF WIND FOR TRAJECTORY RESEARCH.
!              -----------------------------------------------

IF (KWISA == 1) THEN

! * Bilinear interpolations.
  CALL LAIDLI(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
   & KFLDX,ZDLAT,ZDLO(1,1,1),IL0(1,1,1),&
   & PURL0,PUF0)  
  CALL LAIDLI(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
   & KFLDX,ZDLAT,ZDLO(1,1,1),IL0(1,1,1),&
   & PVRL0,PVF0)  

ENDIF

!     ------------------------------------------------------------------

!*       5.    ORIGIN POINT INTERPOLATIONS.
!              ----------------------------

IF (KWISA == 3) THEN

  IF(.NOT.LDQMHW) THEN
    CALL LAIDDI(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,0,ZCLA,ZDLO,ZCLO,IL0(1,1,0),&
     & PUSL,PUF)  
    CALL LAIDDI(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,0,ZCLA,ZDLO,ZCLO,IL0(1,1,0),&
     & PVSL,PVF)  
  ELSE
    CALL LAIDDI(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,1,ZCLA,ZDLO,ZCLO,IL0(1,1,0),&
     & PUSL,PUF)  
    CALL LAIDDI(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,1,ZCLA,ZDLO,ZCLO,IL0(1,1,0),&
     & PVSL,PVF)  
  ENDIF
  IF(.NOT.LDQMHP) THEN
    CALL LAIDDI(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,0,ZCLA,ZDLO,ZCLO,IL0(1,1,0),&
     & PCSL,PCF)  
  ELSE
    CALL LAIDDI(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,1,ZCLA,ZDLO,ZCLO,IL0(1,1,0),&
     & PCSL,PCF)  
  ENDIF

  IF (IWLAG == 3) THEN
!   Save PUF,PVF in case of refined treatment of
!   Coriolis term in 2TL scheme:
    PUFZ(KSTART:KPROF)=PUF(KSTART:KPROF)
    PVFZ(KSTART:KPROF)=PVF(KSTART:KPROF)
    CALL LAIDLI(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,ZDLAT,ZDLO(1,1,1),IL0(1,1,1),&
     & PUSL2,ZUF)  
    CALL LAIDLI(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,ZDLAT,ZDLO(1,1,1),IL0(1,1,1),&
     & PVSL2,ZVF)  
    DO JROF=KSTART,KPROF
      PUF(JROF)=PUF(JROF)+ZUF(JROF)
      PVF(JROF)=PVF(JROF)+ZVF(JROF)
    ENDDO
  ENDIF

  IF (IVLAG == 3) THEN
    CALL LAIDLI(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,ZDLAT,ZDLO(1,1,1),IL0(1,1,1),&
     & PCSL2,ZTF)  
    DO JROF=KSTART,KPROF
      PCF(JROF)=PCF(JROF)+ZTF(JROF)
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LARCIN2',1,ZHOOK_HANDLE)
END SUBROUTINE LARCIN2

