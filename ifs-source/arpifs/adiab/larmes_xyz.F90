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

#ifdef RS6K
@PROCESS NOCHECK
#endif
!option! -O nomove
SUBROUTINE LARMES_XYZ(YDGEOMETRY,YDRIP,YDML_DYN,KST,KPROF,&
 & YDSL,KSTABUF,PB1,PB2,PWRL9,&
 & PLSDEPI,KIBL,KVSEPC,KVSEPL,&
 & PSAVEDP,PSCO,PLEV,PCCO,PUF,PVF,PWF,&
 & KL0,KLH0,KLEV)  

!----compiled for Cray with -hcontiguous----

!**** *LARMES_XYZ - semi-LAgrangian scheme:
!       calculates DP in Cartesian co-ordinate system XYZ

!     Purpose.
!     --------

!     The computation of the location of the interpolation point of
!     the lagrangian trajectory is performed by SETTLS iterative
!     scheme expressed in a Cartesian/eta co-ordinate system (X,Y,Z,eta)
!     where X,Y,Z are coordinates of points on the surface of the sphere while
!     eta is the vertical coordinate. It computes lat/lon without need to call
!     LARCHE
!     

!**   Interface.
!     ----------
!        *CALL* *LARMES_XYZ(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST      - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          YDSL     - SL_STRUCT definition
!          KSTABUF  - for a latitude IGL, KSTABUF(IGL) is the
!                     address of the element corresponding to
!                     (ILON=1,IGL) in the NPROMA arrays.
!          PB1      - SLBUF1-buffer for interpolations.
!          PB2      - SLBUF2-buf to communicate info from non lag. to lag. dyn.
!          PWRL9    - Etadot vertical velocity at t-dt
!          PLSDEPI  - (Number of points by latitude) / (2 * PI) .
!          KIBL     - index into YRCSGEOM/YRGSGEOM instances in YDGEOMETRY

!        INPUT/OUTPUT:
!          KVSEPC   - vertical separation (used in S/L adjoint, cubic interp.)
!          KVSEPL   - vertical separation (used in S/L adjoint, linear interp.)
!          PSAVEDP  - departure point coordinates (lon,lat,eta) from previous timestep
!
!          PSCO     - information about geographic position of interpol. point. (not used
!                     here but provided for compatibility with larcina)
!          PLEV     - vertical coordinate of the interpolation point.
!          PCCO     - information about comput. space position of interpol. point.
!          PUF      - U-comp of wind necessary to
!                     find the position of the origin point,
!                     in a local repere linked to computational sphere.
!          PVF      - V-comp of wind necessary to
!                     find the position of the origin point,
!                     in a local repere linked to computational sphere.
!          PWF      - interpolated etadot used to find the vertical displacement.
!          KL0      - index of the four western points
!                     of the 16 points interpolation grid.
!          KLH0     - second value of index of the four western points
!                     of the 16 points interpolation grid if needed.
!          KLEV     - lower level of the vertical interpolation
!                     grid needed for vertical interpolations.


!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See tech memo XXX.

!     Externals.
!     ----------
!        Calls  LARCINA.
!        Is called by LAPINEA (3D model)

!     Reference.
!     ----------

!     Author.
!     -------
!     M. Diamantakis

!     Modifications.
!     --------------
! End Modifications
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST             , ONLY : RPI
USE YOMLUN             , ONLY : NULERR
USE YOMRIP             , ONLY : TRIP
USE YOMCT3             , ONLY : NSTEP
USE EINT_MOD           , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB1(YDSL%NASLB1,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWRL9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
INTEGER(KIND=JPIM),INTENT(INOUT) :: KVSEPC
INTEGER(KIND=JPIM),INTENT(INOUT) :: KVSEPL
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSAVEDP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,3)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KLH0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ILEVEXP, IROFEXP
INTEGER(KIND=JPIM) :: IHVI, IROT, ITIP
INTEGER(KIND=JPIM) :: JLEV, JROF, JITER
INTEGER(KIND=JPIM) :: ISTEPSAVE
INTEGER(KIND=JPIM) :: ISTESB(YDGEOMETRY%YRDIM%NPROMA)
INTEGER(KIND=JPIM) :: ISTEST(YDGEOMETRY%YRDIM%NPROMA)
INTEGER(KIND=JPIM) :: ISTEST0(YDGEOMETRY%YRDIM%NPROMA)
INTEGER(KIND=JPIM) :: IDEP (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: INOWENO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

LOGICAL :: LLINTV, LLSLHD, LLSLHDQUAD
LOGICAL :: LL_SLDPSAVE_ACTIVE, LL_PCDP3

REAL(KIND=JPRB) :: ZDT2, ZEW, ZEWX, ZNOR, ZNORX, ZVMAX1, ZVMAX2
REAL(KIND=JPRB) :: ZVETAON, ZVETAOX, ZLEVB, ZLEVO, ZLEVT
REAL(KIND=JPRD) :: ZRDP, ZXDP, ZYDP, ZZDP
REAL(KIND=JPRB) :: ZLON2, ZPI2
REAL(KIND=JPRB) :: ZSINLAT0, ZCOSLAT0, ZSINDLON, ZCOSDLON
REAL(KIND=JPRB) :: ZQ11, ZQ12, ZQ21, ZQ22, ZQ33
REAL(KIND=JPRB) :: ZUN_KAPPA(1),ZUN_KAPPAT(1),ZUN_KAPPAM(1),ZUN_KAPPAH(1),ZWFSM(1)

REAL(KIND=JPRB) :: ZWF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZLSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM)
REAL(KIND=JPRB) :: ZRSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "larcina.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LARMES_XYZ',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM,  &
 & YDMP=>YDGEOMETRY%YRMP,  YDSTA=>YDGEOMETRY%YRSTA, YDVETA=>YDGEOMETRY%YRVETA,        &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL), YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL),          &
 & YDDYN=>YDML_DYN%YRDYN,YDDYNA=>YDML_DYN%YRDYNA, &
 & YDPTRSLB1=>YDML_DYN%YRPTRSLB1,YDPTRSLB2=>YDML_DYN%YRPTRSLB2, &
 & YDTCCO=>YDML_DYN%YYTCCO)

ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & LFINDVSEP=>YDDYN%LFINDVSEP, NITMP=>YDDYN%NITMP, &
 & NPC_ITER=>YDDYN%NCURRENT_ITER, &
 & VETAON=>YDDYN%VETAON, VETAOX=>YDDYN%VETAOX, LRHS_CURV=>YDDYN%LRHS_CURV, &
 & VMAX1=>YDDYN%VMAX1, VMAX2=>YDDYN%VMAX2,LSLDP_SAVE=>YDDYN%LSLDP_SAVE, &
 & MSLB1UR0=>YDPTRSLB1%MSLB1UR0, MSLB1VR0=>YDPTRSLB1%MSLB1VR0, &
 & MSLB1ZR0=>YDPTRSLB1%MSLB1ZR0, MSLB1WR0=>YDPTRSLB1%MSLB1WR0, &
 & NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & MSLB2STDDISU=>YDPTRSLB2%MSLB2STDDISU, MSLB2STDDISV=>YDPTRSLB2%MSLB2STDDISV, &
 & MSLB2STDDISW=>YDPTRSLB2%MSLB2STDDISW, MSLB2URL=>YDPTRSLB2%MSLB2URL, &
 & MSLB2VRL=>YDPTRSLB2%MSLB2VRL, MSLB2WRL=>YDPTRSLB2%MSLB2WRL,         & 
 & MSLB2ZRL=>YDPTRSLB2%MSLB2ZRL, NFLDSLB2=>YDPTRSLB2%NFLDSLB2,         &
 & RTDT=>YDRIP%RTDT)
!------------------------------------------------------------------
!             PRELIMINARY INITIALISATIONS AND TESTS.
!------------------------------------------------------------------
ZPI2=2.0_JPRB*RPI
ZDT2=0.5_JPRB*RTDT
ZLEVT=YDVETA%VETAF(0)
ZLEVB=YDVETA%VETAF(NFLEVG+1)
ZVETAON=(1.0_JPRB-VETAON)*YDVETA%VETAH(0)+VETAON*YDVETA%VETAF(1)
ZVETAOX=(1.0_JPRB-VETAOX)*YDVETA%VETAH(NFLEVG)+VETAOX*YDVETA%VETAF(NFLEVG)
ISTEST0(:)=0
ISTEST(:)=0
ISTESB(:)=0
ISTEPSAVE=1
IF (LSLDP_SAVE .AND. NSTEP<ISTEPSAVE) THEN
  LL_SLDPSAVE_ACTIVE=.FALSE.
ELSE
  LL_SLDPSAVE_ACTIVE=LSLDP_SAVE
ENDIF
! If full PC scheme is used always use saved DP for eta coordinate
LL_PCDP3=LL_SLDPSAVE_ACTIVE.AND.YDDYNA%LPC_FULL.AND.(.NOT.YDDYNA%LPC_CHEAP).AND.NPC_ITER>0
!*       Flags for LARCINA call
!
IHVI=0
IROT=0
ITIP=1
LLINTV=.TRUE.
LLSLHD=.FALSE.
LLSLHDQUAD=.FALSE.

!*       1.1   Test that wind is not too strong.
ZNOR=0.0_JPRB
ZEW=0.0_JPRB
ZVMAX1=VMAX1*VMAX1
ZVMAX2=VMAX2*VMAX2

DO JLEV=1,NFLEVG
  DO JROF=KST,KPROF
    ZNOR=MAX(PB2(JROF,MSLB2VRL+JLEV-1)*PB2(JROF,MSLB2VRL+JLEV-1),ZNOR)
    ZEW =MAX(PB2(JROF,MSLB2URL+JLEV-1)*PB2(JROF,MSLB2URL+JLEV-1),ZEW )
  ENDDO
ENDDO

IF (ZNOR > ZVMAX1) WRITE(NULERR,*) ' MAX V WIND=',SQRT(ZNOR)
IF (ZEW  > ZVMAX1) WRITE(NULERR,*) ' MAX U WIND=',SQRT(ZEW)

IF (ZNOR > ZVMAX2) THEN
  ZNOR=0.0_JPRB
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZNORX=PB2(JROF,MSLB2VRL+JLEV-1)*PB2(JROF,MSLB2VRL+JLEV-1)
      ZNOR=MAX(ZNORX,ZNOR)
      IF(ZNOR == ZNORX) THEN
        ILEVEXP=JLEV
        IROFEXP=JROF
      ENDIF
    ENDDO
  ENDDO
  WRITE(NULERR,*) ' V WIND =',SQRT(ZNOR),' IS TOO STRONG, EXPLOSION.'
  WRITE(NULERR,*) ' LEVEL= ',ILEVEXP,' POINT= ',IROFEXP
  WRITE(NULERR,*) ' LON  = ',ACOS(YDCSGEOM%RCOLON(IROFEXP))*180._JPRB/RPI,' degrees'
  WRITE(NULERR,*) ' LAT  = ',ASIN(YDGSGEOM%GEMU(IROFEXP))*180._JPRB/RPI,' degrees'
  CALL ABOR1(' !V WIND TOO STRONG, EXPLOSION!!!')
ENDIF

IF (ZEW > ZVMAX2) THEN
  ZEW=0.0_JPRB
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZEWX=PB2(JROF,MSLB2URL+JLEV-1)*PB2(JROF,MSLB2URL+JLEV-1)
      ZEW=MAX(ZEWX,ZEW)
      IF(ZEW == ZEWX) THEN
        ILEVEXP=JLEV
        IROFEXP=JROF
      ENDIF
    ENDDO
  ENDDO
  WRITE(NULERR,*) ' U WIND =',SQRT(ZEW),' IS TOO STRONG, EXPLOSION.'
  WRITE(NULERR,*) ' LEVEL= ',ILEVEXP,' POINT= ',IROFEXP
  WRITE(NULERR,*) ' LON  = ',ACOS(YDCSGEOM%RCOLON(IROFEXP))*180._JPRB/RPI,' degrees'
  WRITE(NULERR,*) ' LAT  = ',ASIN(YDGSGEOM%GEMU(IROFEXP))*180._JPRB/RPI,' degrees'
  CALL ABOR1(' !U WIND TOO STRONG, EXPLOSION!!!')
ENDIF

!------------------------------------------------------------------------------------
! Starting values for SETTLS iter in cartesian coordinates X,Y,Z
!------------------------------------------------------------------------------------
IF (LL_SLDPSAVE_ACTIVE) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      PCCO(JROF,JLEV,YDTCCO%M_RLON)=PSAVEDP(JROF,JLEV,1)
      PCCO(JROF,JLEV,YDTCCO%M_RLAT)=PSAVEDP(JROF,JLEV,2)
      IF ( LL_PCDP3.OR.PWRL9(JROF,JLEV)*PB2(JROF,MSLB2WRL+JLEV-1) > 0.0_JPRB ) THEN
        PLEV(JROF,JLEV)=PSAVEDP(JROF,JLEV,3)
      ELSE ! use 1st order Euler predictor if time oscilation detected
        ZLEVO=YDVETA%VETAF(JLEV)-RTDT*PB2(JROF,MSLB2WRL+JLEV-1)
        PLEV(JROF,JLEV)=MIN(ZVETAOX,MAX(ZVETAON,ZLEVO))
      ENDIF
    ENDDO
  ENDDO
ELSE ! 
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZXDP  = YDGSGEOM%GEXCO(JROF) - RTDT*PB2(JROF,MSLB2URL+JLEV-1)
      ZYDP  = YDGSGEOM%GEYCO(JROF) - RTDT*PB2(JROF,MSLB2VRL+JLEV-1)
      ZZDP  = YDGSGEOM%GEZCO(JROF) - RTDT*PB2(JROF,MSLB2ZRL+JLEV-1)
      ZRDP  = ZZDP/SQRT(ZXDP*ZXDP+ZYDP*ZYDP+ZZDP*ZZDP)
! Compute corresponding spherical coordinates of DP: lambda_d, theta_d
      ZLON2 = ATAN2(ZYDP,ZXDP)
      PCCO(JROF,JLEV,YDTCCO%M_RLON)=(0.5_JPRB-SIGN(0.5_JPRB,ZLON2))*ZPI2+ZLON2
      PCCO(JROF,JLEV,YDTCCO%M_RLAT)=ASIN(ZRDP)
      ZLEVO = YDVETA%VETAF(JLEV)   - RTDT*PB2(JROF,MSLB2WRL+JLEV-1)
      PLEV(JROF,JLEV)=MIN(ZVETAOX,MAX(ZVETAON,ZLEVO))
    ENDDO
  ENDDO
ENDIF

!------------------------------------------------------------------------------------
! Iterate to find the DP
!------------------------------------------------------------------------------------
DO JITER=2, NITMP
  CALL LARCINA(YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,IHVI,KSTABUF,LFINDVSEP,LLSLHD,LLSLHDQUAD,LLINTV,&
    & ITIP,IROT,.FALSE.,PLSDEPI,&
    & KIBL,PSCO,PLEV,&
    & ZUN_KAPPA,ZUN_KAPPAT,ZUN_KAPPAM,ZUN_KAPPAH,&
    & PB2(1,MSLB2STDDISU),PB2(1,MSLB2STDDISV),PB2(1,MSLB2STDDISW),&
    & PB1(1,MSLB1UR0),PB1(1,MSLB1VR0),PB1(1,MSLB1ZR0),PB1(1,MSLB1WR0),&
    & KVSEPC,KVSEPL,PCCO,&
    & PUF,PVF,ZZF,ZWF,ZWFSM,&
    & KL0,KLH0,KLEV,ZLSCAW,ZRSCAW,IDEP,INOWENO)
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZXDP  = YDGSGEOM%GEXCO(JROF) - ZDT2*(PUF(JROF,JLEV)+PB2(JROF,MSLB2URL+JLEV-1))
      ZYDP  = YDGSGEOM%GEYCO(JROF) - ZDT2*(PVF(JROF,JLEV)+PB2(JROF,MSLB2VRL+JLEV-1))
      ZZDP  = YDGSGEOM%GEZCO(JROF) - ZDT2*(ZZF(JROF,JLEV)+PB2(JROF,MSLB2ZRL+JLEV-1))
      ZRDP  = ZZDP/SQRT(ZXDP*ZXDP+ZYDP*ZYDP+ZZDP*ZZDP)
! Compute corresponding latitude,longitude of the DP
      ZLON2 = ATAN2(ZYDP,ZXDP)
      PCCO(JROF,JLEV,YDTCCO%M_RLON)=(0.5_JPRB-SIGN(0.5_JPRB,ZLON2))*ZPI2+ZLON2
      PCCO(JROF,JLEV,YDTCCO%M_RLAT)=ASIN(ZRDP)
! Compute vertical coordinate of DP eta_d
      ZLEVO = YDVETA%VETAF(JLEV) - ZDT2*(ZWF(JROF,JLEV)+PB2(JROF,MSLB2WRL+JLEV-1))
      PLEV(JROF,JLEV)=MIN(ZVETAOX,MAX(ZVETAON,ZLEVO))
      IF (JITER==NITMP) THEN
        ISTEST(JROF)=ISTEST(JROF)-MIN(0,MAX(-1,NINT(ZLEVO-ZLEVT-0.5_JPRB)))
        ISTESB(JROF)=ISTESB(JROF)-MIN(0,MAX(-1,NINT(ZLEVB-ZLEVO-0.5_JPRB)))
      ENDIF
    ENDDO
  ENDDO
ENDDO

IF (YDDYNA%LSLDIA) PWF(KST:KPROF,1:NFLEVG)=ZWF(KST:KPROF,1:NFLEVG)

!----------------------------------------------------------------------------------
! Compute rotation matrix elements according to Staniforth et al, QJRMS 2010 
! equivalent to Temperton et al, QJRMS 2001
!----------------------------------------------------------------------------------
IF (.NOT.LRHS_CURV) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZSINLAT0=SIN(PCCO(JROF,JLEV,YDTCCO%M_RLAT))
      ZCOSLAT0=COS(PCCO(JROF,JLEV,YDTCCO%M_RLAT))
      ZSINDLON=SIN(YDGSGEOM%GELAM(JROF)-PCCO(JROF,JLEV,YDTCCO%M_RLON))
      ZCOSDLON=COS(YDGSGEOM%GELAM(JROF)-PCCO(JROF,JLEV,YDTCCO%M_RLON))
      ZQ11=ZCOSDLON
      ZQ12=ZSINLAT0*ZSINDLON
      ZQ21=-YDGSGEOM%GEMU(JROF)*ZSINDLON
      ZQ22=YDGSGEOM%GSQM2(JROF)*ZCOSLAT0+YDGSGEOM%GEMU(JROF)*ZSINLAT0*ZCOSDLON
      ZQ33=YDGSGEOM%GEMU(JROF)*ZSINLAT0+YDGSGEOM%GSQM2(JROF)*ZCOSLAT0*ZCOSDLON
      PCCO(JROF,JLEV,YDTCCO%M_RQX)=(ZQ11+ZQ22)/(1.0_JPRB+ZQ33)
      PCCO(JROF,JLEV,YDTCCO%M_RQY)=(ZQ12-ZQ21)/(1.0_JPRB+ZQ33)
    ENDDO
  ENDDO
ENDIF

!----------------------------------------------------------------------------------
! Save DP at current timestep
!----------------------------------------------------------------------------------
IF (LSLDP_SAVE) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF  
      PSAVEDP(JROF,JLEV,1)=PCCO(JROF,JLEV,YDTCCO%M_RLON)
      PSAVEDP(JROF,JLEV,2)=PCCO(JROF,JLEV,YDTCCO%M_RLAT)
      PSAVEDP(JROF,JLEV,3)=PLEV(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!----------------------------------------------------------------------------------
! Report trajectories outside bounds
!----------------------------------------------------------------------------------
IF( ANY( ISTEST /= ISTEST0 )) THEN
  WRITE(NULERR,'(A,I6,A)') ' SMILAG TRAJECTORY OUT OF ATM ',SUM(ISTEST(KST:KPROF)),' TIMES.'
  ! print statistics
  IF (YDDYNA%LRPRSLTRJ) THEN
    DO JROF=KST,KPROF
      IF( ISTEST(JROF) /= 0 ) THEN
        WRITE(NULERR,*) ' POINT= ',JROF,' MAX ETADOT VERTICAL VEL.= ',MAXVAL(ZWF(JROF,:))
        WRITE(NULERR,*) ' LON  = ',ACOS(YDCSGEOM%RCOLON(JROF))*180._JPRB/RPI,' degrees'
        WRITE(NULERR,*) ' LAT  = ',ASIN(YDGSGEOM%GEMU(JROF))*180._JPRB/RPI,' degrees'
      ENDIF
    ENDDO
  ENDIF
ENDIF

IF( ANY( ISTESB /= ISTEST0 )) THEN
  WRITE(NULERR,'(A,I6,A)') ' SMILAG TRAJECTORY UNDERGROUND ',SUM(ISTESB(KST:KPROF)),' TIMES.'
  ! print statistics
  IF (YDDYNA%LRPRSLTRJ) THEN
    DO JROF=KST,KPROF
      IF(ISTESB(JROF) /= 0 ) THEN
        WRITE(NULERR,*) ' POINT= ',JROF,' MAX ETADOT VERTICAL VEL.= ',MAXVAL(ZWF(JROF,:))
        WRITE(NULERR,*) ' LON  = ',ACOS(YDCSGEOM%RCOLON(JROF))*180._JPRB/RPI,' degrees'
        WRITE(NULERR,*) ' LAT  = ',ASIN(YDGSGEOM%GEMU(JROF))*180._JPRB/RPI,' degrees'
      ENDIF
    ENDDO
  ENDIF
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LARMES_XYZ',1,ZHOOK_HANDLE)
END SUBROUTINE LARMES_XYZ
