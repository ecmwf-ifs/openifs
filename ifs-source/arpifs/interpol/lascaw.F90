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

!option! -O extendreorder
SUBROUTINE LASCAW(&
 ! --- INPUT -------------------------------------------------
 & YDVSPLIP,YDSL,LDREGETA,YDDYNA,KPROMB,KDIMK,KST,KPROF,KFLEV,&
 & KFLDN,KSTABUF,KWIS,KHOR,KWENO,KHVI,&
 & LDSLHD,LDSLHDQUAD,LDSLHD_OLD,LDSLHDHEAT,LD3DTURB,&
 & LDCOMAD,LDCOMADH,LDCOMADV,KSPLTHOI,&
 & P4JP,PIS2,PLSDEPI,PLATI,&
 & PIPI,PSLD,PSLDW,P3DTW,&
 & PLON,PLAT,PLEV,&
 & PVETA,KVAUT,&
 & PVCUICO,PVSLD,PVSLDW,PVSLVF,PGAMMA_WENO,KRLEVX,PVRLEVX,&
 & PKAPPA,PKAPPAT,PKAPPAM,PKAPPAH,&
 & PSTDDISU,PSTDDISV,PSTDDISW,&
 ! --- OUTPUT ------------------------------------------------
 & PDLAT,PDLAMAD,PCLA,PCLASLD,PCLASLT,PCLAMAD,&
 & PDLO ,PDLOMAD,PCLO,PCLOSLD,PCLOSLT,PCLOMAD,&
 & KL0,KLH0,KLEV,KNOWENO,PCW,&
 & PDVER,PDVERMAD,PVINTW,PVINTWSLD,PVINTWSLT,PVINTWSLVF,PVINTWMAD,PVINTWS,&
 & PVDERW,PHVW&
 & )  

!     ------------------------------------------------------------------

!**** *LASCAW  -  Externalisable interpolator:
!                 Storage of Coordinates And Weights.
!                 Spherical geometry version (Gaussian grids)

!     Purpose.
!     --------
!       Determines the interpolation grid:
!       - computation of the coordinates of the
!         point situated at the upper left corner of the 16 points
!         square, and of the interpolation point.
!       - computation of weights.
!       Storage of coordinates and weights.

!       Note that this routine should not know if levels are half levels or full levels;
!       this information must remain in the caller.

!**   Interface.
!     ----------
!        *CALL* *LASCAW( ... )

!        Explicit arguments :
!        --------------------

!        INPUT:
!          YDSL    - SL_STRUCT definition
!          KPROMB  - horizontal dimension for interpolation point quantities.
!          KDIMK   - last dimension for some non-linear weights.
!          KST     - first element of arrays where computations are performed.
!          KPROF   - depth of work.
!          KFLEV   - vertical dimension.
!          KFLDN   - number of the first field.
!          KSTABUF - for a latitude IGL, KSTABUF(IGL) is the
!                    address of the element corresponding to
!                    (ILON=0,IGL) in the KPROMB arrays.
!          KWIS    - kind of interpolation.
!          KHOR    - 0: Horizontal interpolation for each level
!                       of a 3D variable.
!                    1: Interpolation for all origin points corresponding
!                       to a final point of a 2D variable.
!          KWENO   - 1/3 = off/on WENO: extra dimension used for WENO
!          KHVI    - 1/0: filling weights arrays PVDERW and PHVW is necessary/not necessary.
!          LDSLHD  - key activating SLHD weights precomputation
!          LDSLHDQUAD - key activating quadratic weights precomputation
!          LDSLHD_OLD - use old SLHD interpolator
!          LDSLHDHEAT   - If true, the triggering function for heat variables differs from the one for momentum variables
!          LD3DTURB- key activating 3D turbulence weights precomputation
!          LDCOMAD -  key activating COMAD weight computation
!          LDCOMADH-  key activating hor. COMAD
!          LDCOMADV-  key activating ver. COMAD
!          KSPLTHOI- controls additional weights precomputation
!          P4JP    - Approximative inverse of the differences of latitudes.
!          PIS2    - PI / 2
!          PLSDEPI - (Number of points by latitude) / (2 * PI) .
!          PLATI   - latitude on the computational sphere.
!          PIPI    - coefficients for the bicubic interpolations.
!          PSLD    - auxiliary quantity for SLHD interpolation in latitude
!          PSLDW   - weights for SLHD Laplacian smoother in latitude
!          P3DTW   - weights for 3D turb. Laplacian smoother in latitude
!          PLON    - Interpolation point longitude on the computational sphere.
!          PLAT    - Interpolation point latitude on the computational sphere.
!          PLEV    - vertical coordinate of the interpolation point.
!          PVETA   - Values of ETA.
!          KVAUT   - Help table for vertical box search: gives the number
!                    of the layer immediately above eta.
!          PVCUICO - Denominators of the vertical cubic interpolation coefficients
!          PVSLD   - auxiliary quantities for vertical SLHD interpolation
!          PVSLDW  - weights for SLHD vertical Laplacian smoother
!          PVSLVF  - weights for SLVF vertical Laplacian smoother
!          PGAMMA_WENO - weights for vertical WENO interpolation
!          KRLEVX  - Dimension of KVAUT
!          PVRLEVX - REAL(KRLEVX).
!          PKAPPA  - kappa function ("coefficient of SLHD") based on the
!                    rescaled horizontal deformation of the flow evaluated
!                    at instant "t" for the final point F
!          PKAPPAT - PKAPPA for heat variable
!          PKAPPAM - horizontal exchange coefficient for momentum in 3D turb. 
!          PKAPPAH - horizontal exchange coefficient for heat in 3D turb.
!          PSTDDISU- COMAD correction coefficient based on estimated flow deformation
!                    along zonal direction of the trajectory but computed
!                    with wind derivatives at instant "t" at the final point F
!          PSTDDISV- COMAD correction coefficient based on estimated flow deformation
!                    along meridional direction of the trajectory but computed
!                    with wind derivatives at instant "t" at the final point F
!          PSTDDISW- COMAD correction coefficient based on estimated flow deformation
!                    along vertical direction of the trajectory but computed
!                    with wind derivatives at instant "t" at the final point F
!        OUTPUT:
!          PDLAT     - distance for horizontal linear interpolations in latitude
!          PDLAMAD   - PDLAT, COMAD case
!          PCLA      - weights for horizontal cubic interpolations in latitude
!          PCLASLD   - weights for horizontal cubic interpolations in latitude, SLHD case
!          PCLAMAD   - PCLA, COMAD case
!          PCLASLT   - weights for horizontal cubic interpolations in latitude, SLHD case on T
!          PDLO      - distances for horizontal linear interpolations
!                      in longitude (latitude rows 0, 1, 2, 3)
!          PDLOMAD   - PDLO, COMAD case
!          PCLO      - weights for horizontal cubic interpolations in
!                      longitude (latitude rows 1, 2)
!          PCLOSLD   - weights for horizontal cubic interpolations in
!                      longitude, SLHD case (latitude rows 1, 2)
!          PCLOMAD   - PCLO, COMAD case
!          PCLOSLT   - weights for horizontal cubic interpolations in
!                      longitude, SLHD case (latitude rows 1, 2) on T
!          KL0       - index of the four western points
!                      of the 16 points interpolation grid.
!          KLH0      - second value of index of the four western points
!                      of the 16 points interpolation grid if needed.
!          KLEV      - lower level of the vertical interpolation
!                      grid needed for vertical interpolations.
!          KNOWENO   - specific boundary treatment for WENO
!          PCW       - C_k weights for the vertical WENO interpolation.
!          PDVER     - distance for vertical linear interpolation
!          PDVERMAD  - PDVER, COMAD case
!          PVINTW    - vertical cubic interpolation weights
!          PVINTWSLD - vertical cubic interpolation weights, SLHD case
!          PVINTWMAD - PVINTW, COMAD case
!          PVINTWSLT - vertical cubic interpolation weights, SLHD case on T
!          PVINTWSLVF -  vertical cubic interpolation weights, SLVF case on T
!          PVINTWS   - Vertical spline interpolation weights.
!          PVDERW    - weights to compute vertical derivatives necessary for
!                      Hermite cubic vertical interpolation.
!          PHVW      - Hermite vertical cubic interpolation weights.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!      K. YESSAD, after the subroutine LAGINT1
!      written by Maurice IMBARD, Alain CRAPLET and Michel ROCHAS
!      METEO-FRANCE, CNRM/GMAP.
!      Original : MARCH 1992.

!     Modifications.
!     --------------
!      F. Vana 28-Aug-2007  removed 4-points cubic spline interpolation
!      07-Nov-2007 J. Masek   New weights for SLHD interpolators.
!      F. Vana 26-Aug-2008  vectorization support
!      K. Yessad (Dec 2008): remove useless dummy arguments
!      K. Yessad (Feb 2009): split loops, rewrite in a shorter way.
!      R. El Khatib 07-08-2009 Optimisation directive for NEC
!      K. Yessad (Aug 2009): use RIPI, RSLD
!      F. Vana 22-Feb-2011: Extra weights for horiz. turbulence and phy. diff
!      G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      G.Mozdzynski (Aug 2011): support higher order interpolation
!      G. Mozdzynski (May 2012): further cleaning
!      F. Vana 13-feb-2014 SLHD weights for heat variables
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      S. Malardel (Nov 2013): COMAD weights for SL interpolations
!      K. Yessad (March 2017): simplify level numbering.
!      F. Vana, P. Smolikova & A. Craciun (Aug-2017): high order traj research & WENO
!      F. Vana  20-Feb-2019: quintic vertical interpolation 
!      F. Vana  18-Jul-2019: SLVF & cleaning
! End Modifications
!     ------------------------------------------------------------------

USE YOMVSPLIP , ONLY : TVSPLIP
USE PARKIND1  , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK

! arp/ifs dependencies to be solved later.
USE YOMDYNA   , ONLY : TDYNA
USE YOMMP0    , ONLY : NPROC
USE EINT_MOD  , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVSPLIP)     ,INTENT(IN)    :: YDVSPLIP
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
LOGICAL           ,INTENT(IN)    :: LDREGETA
LOGICAL           ,INTENT(IN)    :: LDSLHDHEAT
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMB 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIMK
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KWIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KHOR 
INTEGER(KIND=JPIM),INTENT(IN)    :: KWENO
INTEGER(KIND=JPIM),INTENT(IN)    :: KHVI 
LOGICAL           ,INTENT(IN)    :: LDSLHD
LOGICAL           ,INTENT(IN)    :: LDSLHDQUAD
LOGICAL           ,INTENT(IN)    :: LDSLHD_OLD
LOGICAL           ,INTENT(IN)    :: LD3DTURB
LOGICAL           ,INTENT(IN)    :: LDCOMAD
LOGICAL           ,INTENT(IN)    :: LDCOMADH
LOGICAL           ,INTENT(IN)    :: LDCOMADV
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPLTHOI
INTEGER(KIND=JPIM),INTENT(IN)    :: KRLEVX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P4JP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLATI(YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRD)   ,INTENT(IN)    :: PIPI(YDSL%NDGSAH:YDSL%NDGENH,3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLD(YDSL%NDGSAH:YDSL%NDGENH,3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLDW(3,3,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P3DTW(3,3,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLON(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLEV(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVETA(0:KFLEV+1) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVAUT(0:KRLEVX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCUICO(4,0:KFLEV-1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSLD(3,0:KFLEV-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSLDW(3,3,0:KFLEV-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSLVF(3,3,0:KFLEV-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAMMA_WENO(KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRLEVX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPPA(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPPAT(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPPAM(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPPAH(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTDDISU(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTDDISV(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTDDISW(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLAT(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLAMAD(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLA(KPROMB,KFLEV,3,KDIMK) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLASLD(KPROMB,KFLEV,3,KDIMK)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLASLT(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLAMAD(KPROMB,KFLEV,3,KDIMK)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLO(KPROMB,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLOMAD(KPROMB,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLO(KPROMB,KFLEV,3,2,KDIMK)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLOSLD(KPROMB,KFLEV,3,2,KDIMK)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLOSLT(KPROMB,KFLEV,3,2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLOMAD(KPROMB,KFLEV,3,2,KDIMK)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KL0(KPROMB,KFLEV,0:3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLH0(KPROMB,KFLEV,0:3) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLEV(KPROMB,KFLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KNOWENO(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCW(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDVER(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDVERMAD(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTW(KPROMB,KFLEV,3*KWENO) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWSLD(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWMAD(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWSLT(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWSLVF(KPROMB,KFLEV,3*KWENO)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWS(KPROMB,KFLEV,1:4) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDERW(KPROMB,KFLEV,2*KHVI,2*KHVI) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHVW(KPROMB,KFLEV,4*KHVI)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IADDR(YDSL%NDGSAH:YDSL%NDGENH)

INTEGER(KIND=JPIM) :: IFLVM2,ILEV,ILEVV,JLAT,JLEV,JROF,IBCLIM

INTEGER(KIND=JPIM) :: ILAV(KST:KPROF,KFLEV)
INTEGER(KIND=JPIM) :: IZLATV(KST:KPROF,KFLEV)

INTEGER(KIND=JPIM) :: JJ,JK
INTEGER(KIND=JPIM) :: ILO(KPROMB), ILO1(KPROMB), ILO2(KPROMB), ILO3(KPROMB)
INTEGER(KIND=JPIM) :: ILOIK(KST:KPROF,KFLEV), ILO1IK(KST:KPROF,KFLEV)
INTEGER(KIND=JPIM) :: ILO2IK(KST:KPROF,KFLEV), ILO3IK(KST:KPROF,KFLEV)
INTEGER(KIND=JPIM) :: IILA(KPROMB,KFLEV)
INTEGER(KIND=JPIM) :: IILEV(KPROMB,KFLEV)

REAL(KIND=JPRB) :: PD, ZD2, ZDA, ZDB, ZDC, ZDD,&
 & ZDAMAD, ZDBMAD, ZDCMAD, ZDDMAD,&
 & ZDEN1, ZDEN2, ZDVER, ZFAC, ZLO, ZLO1, ZLO2, ZLO3, ZNUM, ZEPS
REAL(KIND=JPRB) :: ZSLHDKMINH,ZSLHDKMINV
REAL(KIND=JPRB), PARAMETER :: ZSLHDKMINV_WENO=0._JPRB   ! WENO only runs with Lagrangian cubic!!!


REAL(KIND=JPRB) :: ZKHTURB(KPROMB,KFLEV,KDIMK)
LOGICAL         :: LLT_SLHD(4),LLT_PHYS(4),LLSLHD,LLSLHDQUAD,LLSLHD_OLD,LL3DTURB,LLSLVF
LOGICAL         :: LLCOMAD,LLCOMADH,LLCOMADV
REAL(KIND=JPRB) :: ZZWH(KPROMB,3,KFLEV)
REAL(KIND=JPRB) :: ZZWHMAD(KPROMB,3,KFLEV)
REAL(KIND=JPRB) :: ZCLA(KPROMB,KFLEV,3) 
REAL(KIND=JPRB) :: ZCLO(KPROMB,KFLEV,3,2)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------
! functions

REAL(KIND=JPRB) :: FHLO1, FHLO2, FHLO3, FHLO4

! auxiliary functions for Hermite cubic interpolation
FHLO1(PD)= (1.0_JPRB-PD)*(1.0_JPRB-PD)*(1.0_JPRB+2.0_JPRB*PD)
FHLO2(PD)= PD*PD*(3._JPRB-2.0_JPRB*PD)
FHLO3(PD)= PD*(1.0_JPRB-PD)*(1.0_JPRB-PD)
FHLO4(PD)=-PD*PD*(1.0_JPRB-PD)

!     ------------------------------------------------------------------

#include "lascaw_cla.intfb.h"
#include "lascaw_clo.intfb.h"
#include "lascaw_vintw.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LASCAW',0,ZHOOK_HANDLE)
ASSOCIATE(RFVV=>YDVSPLIP%RFVV,&
 & LSLHD=>YDDYNA%LSLHD,LSLHDQUAD=>YDDYNA%LSLHDQUAD,SLHDKMIN=>YDDYNA%SLHDKMIN, &
 & HOISLTH=>YDDYNA%HOISLTH,HOISLTV=>YDDYNA%HOISLTV,LHOISLT=>YDDYNA%LHOISLT, &
 & LSLTVWENO=>YDDYNA%LSLTVWENO,LSLVF=>YDDYNA%LSLVF)
!     ------------------------------------------------------------------

!*       0.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

ZEPS=100.0_JPRB*TINY(1.0_JPRB)

! cases relevant for SLHD scheme (switches LDSLHD, LDSLHDQUAD are
! deactivated during computation of medium points in LAPINEA)
LLSLHD=LDSLHD.AND.(KWIS==103.OR.KWIS==104.OR.KWIS==105.OR.KWIS==106.OR.KWIS==203)
LLSLHDQUAD=(LDSLHDQUAD.AND.(KWIS==103.OR.KWIS==104.OR.KWIS==105.OR.KWIS==106.OR.KWIS==203)).OR. &
   & ((KWIS==102.OR.KWIS==202).AND.LHOISLT.AND.((HOISLTV/=0_JPRB).OR.(HOISLTH/=0_JPRB)))
LL3DTURB=LD3DTURB.AND.(KWIS==103.OR.KWIS==104.OR.KWIS==105.OR.KWIS==106)

! switch for old SLHD scheme
LLSLHD_OLD=LLSLHD.AND.LDSLHD_OLD

LLT_SLHD(1)=LLSLHD
LLT_SLHD(2)=LLSLHDQUAD
LLT_SLHD(3)=LLSLHD_OLD
LLT_SLHD(4)=.FALSE.

! switch for SLVF
LLSLVF=(KWIS /= 202).AND.LSLVF

! cases relevant for COMAD scheme (switches LDCOMADH and LDCOMADV  are
! deactivated during computation of interpolation points in LAPINEA)
LLCOMAD =LDCOMAD.AND.(KWIS==103.OR.KWIS==104.OR.KWIS==105.OR.KWIS==203.OR.KWIS==106)
LLCOMADH=LLCOMAD.AND.LDCOMADH
LLCOMADV=LLCOMAD.AND.LDCOMADV

! switches for interpolation of physics 
! It holds the same value for every iteration step (in ICI scheme). 
LLT_PHYS(1)=LSLHD
LLT_PHYS(2)=LSLHDQUAD
LLT_PHYS(3)=LDSLHD_OLD
LLT_PHYS(4)=.FALSE.

ZSLHDKMINH=SLHDKMIN
ZSLHDKMINV=SLHDKMIN

! Modify previous defaults for high order trajectory research
IF (KWIS==102.OR.KWIS==202) THEN
  ZSLHDKMINH=HOISLTH
  ZSLHDKMINV=HOISLTV
  ! Set LSLHDQUAD for the horizontal interpolation now
  LLT_SLHD(2)= (HOISLTH /= 0._JPRB)
ELSEIF (KWIS==106) THEN
  !Just make sure there is no change from default 3rd order Lagrangian cubic
  LLT_PHYS(2)=.FALSE.
ENDIF

DO JLAT=YDSL%NDGSAH,YDSL%NDGENH
  IADDR(JLAT)=KSTABUF(JLAT)+YDSL%NASLB1*(0-KFLDN)
ENDDO

IF ( LL3DTURB ) THEN
  ! copy from kappa
  ZKHTURB(KST:KPROF,1:KFLEV,2)=PKAPPAM(KST:KPROF,1:KFLEV)
  ZKHTURB(KST:KPROF,1:KFLEV,3)=PKAPPAH(KST:KPROF,1:KFLEV)
ENDIF

IF (KSPLTHOI == 1) THEN
  ! In this case the ZKHTURB is set to static mode with maximum diffusion.
  ZKHTURB(KST:KPROF,1:KFLEV,KDIMK)=1._JPRB
ENDIF

IFLVM2=KFLEV-2

!     ------------------------------------------------------------------

!*       1.    3D MODEL.
!              ---------

!*   distance for horizontal linear interpolations in latitude (common to all options)
#ifdef __INTEL_COMPILER
!$OMP SIMD PRIVATE(JLEV)
DO JROF=KST,KPROF
  DO JLEV=1,KFLEV
#else
DO JLEV=1,KFLEV
  DO JROF=KST,KPROF
#endif
    ! * Calculation of linear weights, KL0.
    IZLATV(JROF,JLEV)=INT(P4JP*(PIS2-PLAT(JROF,JLEV))+0.75_JPRB+ZEPS)-YDSL%NFRSTLOFF
    ILAV(JROF,JLEV)=IZLATV(JROF,JLEV)+NINT(SIGN(0.5_JPRB,PLATI(IZLATV(JROF,JLEV))-PLAT(JROF,JLEV)+ZEPS)-1.5_JPRB)
    PDLAT(JROF,JLEV)=(PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+1))&
         & /(PLATI(ILAV(JROF,JLEV)+2)-PLATI(ILAV(JROF,JLEV)+1)+ZEPS)  
  ENDDO
ENDDO

!        1.01  Coordinates and weights for trilinear interpolations.
IF (KWIS == 101) THEN
  ZFAC=PVRLEVX/(PVETA(KFLEV+1)-PVETA(0))
   !$OMP SIMD PRIVATE(ZLO,ZLO1,ZLO2,ZLO3,JLEV,ILEV)
   DO JROF=KST,KPROF
      DO JLEV=1,KFLEV
         ZLO   =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)  )
         ILOIK(JROF,JLEV)   =INT(ZLO)
         ILOIK(JROF,JLEV)=IADDR(ILAV(JROF,JLEV)  )+YDSL%NSLEXT(ILOIK(JROF,JLEV) ,ILAV(JROF,JLEV) )
         PDLO(JROF,JLEV,0)=ZLO -REAL(INT(ZLO) ,JPRB)
         ZLO1  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+1)
         ILO1IK(JROF,JLEV)  =INT(ZLO1)
         ILO1IK(JROF,JLEV)=IADDR(ILAV(JROF,JLEV)+1)+YDSL%NSLEXT(ILO1IK(JROF,JLEV),ILAV(JROF,JLEV)+1)
         PDLO(JROF,JLEV,1)=ZLO1-REAL(INT(ZLO1),JPRB)
         ZLO2  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+2)
         ILO2IK(JROF,JLEV)  =INT(ZLO2)
         ILO2IK(JROF,JLEV)=IADDR(ILAV(JROF,JLEV)+2)+YDSL%NSLEXT(ILO2IK(JROF,JLEV),ILAV(JROF,JLEV)+2)
         PDLO(JROF,JLEV,2)=ZLO2-REAL(INT(ZLO2),JPRB)
         ZLO3  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+3)
         ILO3IK(JROF,JLEV)  =INT(ZLO3)
         ILO3IK(JROF,JLEV)=IADDR(ILAV(JROF,JLEV)+3)+YDSL%NSLEXT(ILO3IK(JROF,JLEV),ILAV(JROF,JLEV)+3)
         PDLO(JROF,JLEV,3)=ZLO3-REAL(INT(ZLO3),JPRB)
      ENDDO
   ENDDO
  ! NOTE: Loop merged with previous loop
  ! DO JLEV=1,KFLEV
!DIR$ PREFERVECTOR
  !   DO JROF=KST,KPROF
  !     ILOIK(JROF,JLEV)=IADDR(ILAV(JROF,JLEV)  )+YDSL%NSLEXT(ILOIK(JROF,JLEV) ,ILAV(JROF,JLEV) )
  !     ILO1IK(JROF,JLEV)=IADDR(ILAV(JROF,JLEV)+1)+YDSL%NSLEXT(ILO1IK(JROF,JLEV),ILAV(JROF,JLEV)+1)
  !     ILO2IK(JROF,JLEV)=IADDR(ILAV(JROF,JLEV)+2)+YDSL%NSLEXT(ILO2IK(JROF,JLEV),ILAV(JROF,JLEV)+2)
  !     ILO3IK(JROF,JLEV)=IADDR(ILAV(JROF,JLEV)+3)+YDSL%NSLEXT(ILO3IK(JROF,JLEV),ILAV(JROF,JLEV)+3)
  !   ENDDO
  ! ENDDO

  !$OMP SIMD PRIVATE(JLEV,ILEV)
  DO JROF=KST,KPROF
     DO JLEV=1,KFLEV
        ILEV  =KVAUT(INT(PLEV(JROF,JLEV)*ZFAC))-1
        IF(ILEV < IFLVM2.AND.&
             & (PLEV(JROF,JLEV)-PVETA(ILEV+2)) > 0.0_JPRB) ILEV=ILEV+1  
        KLEV(JROF,JLEV)=ILEV
        KL0(JROF,JLEV,0)=ILOIK(JROF,JLEV)+YDSL%NASLB1*ILEV
        KL0(JROF,JLEV,1)=ILO1IK(JROF,JLEV)+YDSL%NASLB1*ILEV
        KL0(JROF,JLEV,2)=ILO2IK(JROF,JLEV)+YDSL%NASLB1*ILEV
        KL0(JROF,JLEV,3)=ILO3IK(JROF,JLEV)+YDSL%NASLB1*ILEV
     ENDDO
     ! NOTE: Loop partially combined with the previous loop
     DO JLEV=1,KFLEV
        ILEV=KLEV(JROF,JLEV)
        PDVER(JROF,JLEV)=(PLEV(JROF,JLEV)-PVETA(ILEV+1))/&
             & (PVETA(ILEV+2)-PVETA(ILEV+1))
     ENDDO
  ENDDO
#ifdef __INTEL_COMPILER
    ! * Mask calculation for on-demand communications:
  IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
     !$OMP SIMD PRIVATE(JROF,JJ,JK)
     DO JLEV=1,KFLEV
        DO JROF=KST,KPROF
           JJ=ILOIK (JROF,JLEV)
           DO JK=JJ,JJ+3
              YDSL%MASK_SL2(JK)=1
           ENDDO
        ENDDO
     ENDDO
     !$OMP SIMD PRIVATE(JROF,JJ,JK)
     DO JLEV=1,KFLEV
        DO JROF=KST,KPROF
           JJ=ILO1IK(JROF,JLEV)
           DO JK=JJ,JJ+3
              YDSL%MASK_SL2(JK)=1
           ENDDO
        ENDDO
     ENDDO
     !$OMP SIMD PRIVATE(JROF,JJ,JK)
     DO JLEV=1,KFLEV
        DO JROF=KST,KPROF
           JJ=ILO2IK(JROF,JLEV)
           DO JK=JJ,JJ+3
              YDSL%MASK_SL2(JK)=1
           ENDDO
        ENDDO
     ENDDO
     !$OMP SIMD PRIVATE(JROF,JJ,JK)
     DO JLEV=1,KFLEV
        DO JROF=KST,KPROF
           JJ=ILO3IK(JROF,JLEV)
           DO JK=JJ,JJ+3
              YDSL%MASK_SL2(JK)=1
           ENDDO
        ENDDO
     ENDDO
  ENDIF
#else
    ! * Mask calculation for on-demand communications:
  IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        JJ=ILOIK (JROF,JLEV)
        YDSL%MASK_SL2(JJ:JJ+3)=1
        JJ=ILO1IK(JROF,JLEV)
        YDSL%MASK_SL2(JJ:JJ+3)=1
        JJ=ILO2IK(JROF,JLEV)
        YDSL%MASK_SL2(JJ:JJ+3)=1
        JJ=ILO3IK(JROF,JLEV)
        YDSL%MASK_SL2(JJ:JJ+3)=1
      ENDDO
    ENDDO
  ENDIF
#endif

ENDIF

!        1.03  Coordinates and weights for ( horizontal 12 points
!              + vertical cubic + 32 points interpolations ) or
!              ( horizontal 12 points + 32 points interpolations ).
!              Optionally, Hermite cubic, cubic B-spline or WENO vertical interpolations weights are computed.

IF (KWIS == 102 .OR. KWIS == 103 .OR. KWIS == 104 .OR. KWIS == 105 .OR. KWIS == 106) THEN

ZFAC=PVRLEVX/(PVETA(KFLEV+1)-PVETA(0))

IF(.NOT.LLCOMADH.AND..NOT.LLCOMADV)THEN


  DO JLEV=1,KFLEV
    ! * Calculation of linear weights, KL0, KLH0.
!CDIR NODEP
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF

      ! meridional interpolation: linear weights and input for cubic weights (LASCAW_CLA)
      ! general case 
      ZDA   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV))
      ZDB   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+1)
      ZDC   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+2)
      ZDD   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+3)
      ZZWH(JROF,1,JLEV)=(ZDA*ZDC)*ZDD*PIPI(ILAV(JROF,JLEV)+1,1)
      ZZWH(JROF,2,JLEV)=(ZDA*ZDB)*ZDD*PIPI(ILAV(JROF,JLEV)+1,2)
      ZZWH(JROF,3,JLEV)=(ZDA*ZDB)*ZDC*PIPI(ILAV(JROF,JLEV)+1,3)
    ENDDO
    IILA(KST:KPROF,JLEV)=ILAV(KST:KPROF,JLEV)
    PDLAMAD(KST:KPROF,JLEV)=PDLAT(KST:KPROF,JLEV)
    ZZWHMAD(KST:KPROF,1:3,JLEV)=ZZWH(KST:KPROF,1:3,JLEV)
  ENDDO

  DO JLEV=1,KFLEV
!CDIR NODEP
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      ZLO   =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)  )
      ILOIK(JROF,JLEV)   =INT(ZLO)
      PDLO(JROF,JLEV,0)=ZLO -REAL(INT(ZLO) ,JPRB)
      ZLO1  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+1)
      ILO1IK(JROF,JLEV)  =INT(ZLO1)
      PDLO(JROF,JLEV,1)=ZLO1-REAL(INT(ZLO1),JPRB)
      ZLO2  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+2)
      ILO2IK(JROF,JLEV)  =INT(ZLO2)
      PDLO(JROF,JLEV,2)=ZLO2-REAL(INT(ZLO2),JPRB)
      ZLO3  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+3)
      ILO3IK(JROF,JLEV)  =INT(ZLO3)
      PDLO(JROF,JLEV,3)=ZLO3-REAL(INT(ZLO3),JPRB)
    ENDDO
  ENDDO
  DO JLEV=1,KFLEV
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      ILOIK(JROF,JLEV)=IADDR(ILAV(JROF,JLEV)  )+YDSL%NSLEXT(ILOIK(JROF,JLEV) ,ILAV(JROF,JLEV) )
      ILO1IK(JROF,JLEV)=IADDR(ILAV(JROF,JLEV)+1)+YDSL%NSLEXT(ILO1IK(JROF,JLEV),ILAV(JROF,JLEV)+1)
      ILO2IK(JROF,JLEV)=IADDR(ILAV(JROF,JLEV)+2)+YDSL%NSLEXT(ILO2IK(JROF,JLEV),ILAV(JROF,JLEV)+2)
      ILO3IK(JROF,JLEV)=IADDR(ILAV(JROF,JLEV)+3)+YDSL%NSLEXT(ILO3IK(JROF,JLEV),ILAV(JROF,JLEV)+3)
    ENDDO
  ENDDO

  DO JLEV=1,KFLEV
    PDLOMAD(KST:KPROF,JLEV,0)=PDLO(KST:KPROF,JLEV,0)
    PDLOMAD(KST:KPROF,JLEV,1)=PDLO(KST:KPROF,JLEV,1)
    PDLOMAD(KST:KPROF,JLEV,2)=PDLO(KST:KPROF,JLEV,2)
    PDLOMAD(KST:KPROF,JLEV,3)=PDLO(KST:KPROF,JLEV,3)
  ENDDO

  DO JLEV=1,KFLEV
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      ! vertical interpolation: linear weights
      ! the cubic weight computation are done in 
      ! LASCAW_VINTW (including terms for grid irregularity)
      ILEVV=KVAUT(INT(PLEV(JROF,JLEV)*ZFAC))-1
      IF(ILEVV < IFLVM2.AND.&
       & (PLEV(JROF,JLEV)-PVETA(ILEVV+2)) > 0.0_JPRB) ILEVV=ILEVV+1  
      KLEV(JROF,JLEV)=ILEVV
      ! general case
      PDVER(JROF,JLEV)=(PLEV(JROF,JLEV)-PVETA(ILEVV+1))/&
       & (PVETA(ILEVV+2)-PVETA(ILEVV+1))  
      PDVERMAD(JROF,JLEV)=PDVER(JROF,JLEV)
    ENDDO
  ENDDO

!CDIR NODEP
  DO JLEV=1,KFLEV
    IF (KHOR == 0) ILEV=JLEV
    IF (KHOR == 1) ILEV=KFLDN
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      KL0(JROF,JLEV,0)=ILOIK(JROF,JLEV)+YDSL%NASLB1*KLEV(JROF,JLEV)
      KL0(JROF,JLEV,1)=ILO1IK(JROF,JLEV)+YDSL%NASLB1*KLEV(JROF,JLEV)
      KL0(JROF,JLEV,2)=ILO2IK(JROF,JLEV)+YDSL%NASLB1*KLEV(JROF,JLEV)
      KL0(JROF,JLEV,3)=ILO3IK(JROF,JLEV)+YDSL%NASLB1*KLEV(JROF,JLEV)

      KLH0(JROF,JLEV,0)=ILOIK(JROF,JLEV)+YDSL%NASLB1*ILEV
      KLH0(JROF,JLEV,1)=ILO1IK(JROF,JLEV)+YDSL%NASLB1*ILEV
      KLH0(JROF,JLEV,2)=ILO2IK(JROF,JLEV)+YDSL%NASLB1*ILEV
      KLH0(JROF,JLEV,3)=ILO3IK(JROF,JLEV)+YDSL%NASLB1*ILEV

    ENDDO
  ENDDO

    ! * Mask calculation for on-demand communications:
  IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        JJ=ILOIK (JROF,JLEV)
        YDSL%MASK_SL2(JJ:JJ+3)=1
        JJ=ILO1IK(JROF,JLEV)
        YDSL%MASK_SL2(JJ:JJ+3)=1
        JJ=ILO2IK(JROF,JLEV)
        YDSL%MASK_SL2(JJ:JJ+3)=1
        JJ=ILO3IK(JROF,JLEV)
        YDSL%MASK_SL2(JJ:JJ+3)=1
      ENDDO
    ENDDO
  ENDIF

ELSE

  DO JLEV=1,KFLEV
    IF (KHOR == 0) ILEV=JLEV
    IF (KHOR == 1) ILEV=KFLDN

    ! * Calculation of linear weights, KL0, KLH0.
!CDIR NODEP
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF

      ! meridional interpolation: linear weights and input for cubic weights (LASCAW_CLA)
      ! general case 
      ZDA   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV))
      ZDB   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+1)
      ZDC   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+2)
      ZDD   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+3)
      IILA(JROF,JLEV)=ILAV(JROF,JLEV)
      ZZWH(JROF,1,JLEV)=(ZDA*ZDC)*ZDD*PIPI(ILAV(JROF,JLEV)+1,1)
      ZZWH(JROF,2,JLEV)=(ZDA*ZDB)*ZDD*PIPI(ILAV(JROF,JLEV)+1,2)
      ZZWH(JROF,3,JLEV)=(ZDA*ZDB)*ZDC*PIPI(ILAV(JROF,JLEV)+1,3)

      ! COMAD meridional interpolation 
      IF (LLCOMADH) THEN
        PDLAMAD(JROF,JLEV)=PDLAT(JROF,JLEV)*PSTDDISV(JROF,JLEV) + 0.5_JPRB * (1._JPRB-PSTDDISV(JROF,JLEV))
        ZDAMAD = ZDA-ZDB
        ZDDMAD = ZDD-ZDC
        ZDBMAD = ZDB*PSTDDISV(JROF,JLEV) +0.5_JPRB*(PLATI(ILAV(JROF,JLEV)+2)-PLATI(ILAV(JROF,JLEV)+1)) *&
         & (1._JPRB-PSTDDISV(JROF,JLEV))
        ZDCMAD = ZDC*PSTDDISV(JROF,JLEV) -0.5_JPRB*(PLATI(ILAV(JROF,JLEV)+2)-PLATI(ILAV(JROF,JLEV)+1)) *&
         & (1._JPRB-PSTDDISV(JROF,JLEV))
        ZDAMAD = ZDAMAD + ZDBMAD
        ZDDMAD = ZDDMAD + ZDCMAD
        IILA(JROF,JLEV)=ILAV(JROF,JLEV)
        ZZWHMAD(JROF,1,JLEV)=(ZDAMAD*ZDCMAD)*ZDDMAD*PIPI(ILAV(JROF,JLEV)+1,1)
        ZZWHMAD(JROF,2,JLEV)=(ZDAMAD*ZDBMAD)*ZDDMAD*PIPI(ILAV(JROF,JLEV)+1,2)
        ZZWHMAD(JROF,3,JLEV)=(ZDAMAD*ZDBMAD)*ZDCMAD*PIPI(ILAV(JROF,JLEV)+1,3)
      ELSE
        PDLAMAD(JROF,JLEV)=PDLAT(JROF,JLEV)
        ZZWHMAD(JROF,1,JLEV)=ZZWH(JROF,1,JLEV)
        ZZWHMAD(JROF,2,JLEV)=ZZWH(JROF,2,JLEV)
        ZZWHMAD(JROF,3,JLEV)=ZZWH(JROF,3,JLEV)
      ENDIF

      ! zonal interpolation: linear weights for 4 lat. cicles
      ! as the grid is regular in the zonal direction,
      ! the cubic weight computation does not need 
      ! other input than linear weights (LASCAW_CLO)
      ! general case
      ZLO   =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)  )
      ILO(JROF)   =INT(ZLO )
      PDLO(JROF,JLEV,0)=ZLO -REAL(ILO(JROF) ,JPRB)
      ZLO1  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+1)
      ILO1(JROF)  =INT(ZLO1)
      PDLO(JROF,JLEV,1)=ZLO1-REAL(ILO1(JROF),JPRB)
      ZLO2  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+2)
      ILO2(JROF)  =INT(ZLO2)
      PDLO(JROF,JLEV,2)=ZLO2-REAL(ILO2(JROF),JPRB)
      ZLO3  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+3)
      ILO3(JROF)  =INT(ZLO3)
      PDLO(JROF,JLEV,3)=ZLO3-REAL(ILO3(JROF),JPRB)

      ! COMAD zonal interpolation 
      IF (LLCOMADH) THEN
        PDLOMAD(JROF,JLEV,0)=PDLO(JROF,JLEV,0)*PSTDDISU(JROF,JLEV) + 0.5_JPRB * (1._JPRB-PSTDDISU(JROF,JLEV))
        PDLOMAD(JROF,JLEV,1)=PDLO(JROF,JLEV,1)*PSTDDISU(JROF,JLEV) + 0.5_JPRB * (1._JPRB-PSTDDISU(JROF,JLEV))
        PDLOMAD(JROF,JLEV,2)=PDLO(JROF,JLEV,2)*PSTDDISU(JROF,JLEV) + 0.5_JPRB * (1._JPRB-PSTDDISU(JROF,JLEV))
        PDLOMAD(JROF,JLEV,3)=PDLO(JROF,JLEV,3)*PSTDDISU(JROF,JLEV) + 0.5_JPRB * (1._JPRB-PSTDDISU(JROF,JLEV))
      ELSE
        PDLOMAD(JROF,JLEV,0)=PDLO(JROF,JLEV,0)
        PDLOMAD(JROF,JLEV,1)=PDLO(JROF,JLEV,1)
        PDLOMAD(JROF,JLEV,2)=PDLO(JROF,JLEV,2)
        PDLOMAD(JROF,JLEV,3)=PDLO(JROF,JLEV,3)
      ENDIF

      ! vertical interpolation: linear weights
      ! the cubic weight computation are done in 
      ! LASCAW_VINTW (including terms for grid irregularity)
      ILEVV=KVAUT(INT(PLEV(JROF,JLEV)*ZFAC))-1
      IF(ILEVV < IFLVM2.AND.&
       & (PLEV(JROF,JLEV)-PVETA(ILEVV+2)) > 0.0_JPRB) ILEVV=ILEVV+1  
      KLEV(JROF,JLEV)=ILEVV
      ! general case
      PDVER(JROF,JLEV)=(PLEV(JROF,JLEV)-PVETA(ILEVV+1))/&
       & (PVETA(ILEVV+2)-PVETA(ILEVV+1))  
      ! COMAD vertical interpolation 
      IF (LLCOMADV) THEN
        PDVERMAD(JROF,JLEV)=PDVER(JROF,JLEV)*PSTDDISW(JROF,JLEV) + 0.5_JPRB * (1._JPRB-PSTDDISW(JROF,JLEV))
      ELSE
        PDVERMAD(JROF,JLEV)=PDVER(JROF,JLEV)
      ENDIF

      ILO(JROF)=IADDR(ILAV(JROF,JLEV)  )+YDSL%NSLEXT(ILO(JROF) ,ILAV(JROF,JLEV)  )
      ILO1(JROF)=IADDR(ILAV(JROF,JLEV)+1)+YDSL%NSLEXT(ILO1(JROF),ILAV(JROF,JLEV)+1)
      ILO2(JROF)=IADDR(ILAV(JROF,JLEV)+2)+YDSL%NSLEXT(ILO2(JROF),ILAV(JROF,JLEV)+2)
      ILO3(JROF)=IADDR(ILAV(JROF,JLEV)+3)+YDSL%NSLEXT(ILO3(JROF),ILAV(JROF,JLEV)+3)

      KL0(JROF,JLEV,0)=ILO(JROF)+YDSL%NASLB1*ILEVV
      KL0(JROF,JLEV,1)=ILO1(JROF)+YDSL%NASLB1*ILEVV
      KL0(JROF,JLEV,2)=ILO2(JROF)+YDSL%NASLB1*ILEVV
      KL0(JROF,JLEV,3)=ILO3(JROF)+YDSL%NASLB1*ILEVV

      KLH0(JROF,JLEV,0)=ILO(JROF)+YDSL%NASLB1*ILEV
      KLH0(JROF,JLEV,1)=ILO1(JROF)+YDSL%NASLB1*ILEV
      KLH0(JROF,JLEV,2)=ILO2(JROF)+YDSL%NASLB1*ILEV
      KLH0(JROF,JLEV,3)=ILO3(JROF)+YDSL%NASLB1*ILEV

    ENDDO
  ENDDO

    ! * Mask calculation for on-demand communications:
  IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        JJ=ILO (JROF)
        YDSL%MASK_SL2(JJ:JJ+3)=1
        JJ=ILO1(JROF)
        YDSL%MASK_SL2(JJ:JJ+3)=1
        JJ=ILO2(JROF)
        YDSL%MASK_SL2(JJ:JJ+3)=1
        JJ=ILO3(JROF)
        YDSL%MASK_SL2(JJ:JJ+3)=1
      ENDDO
    ENDDO
  ENDIF

ENDIF



  IF (LLSLHD.AND.LDSLHDHEAT) THEN
    ! Computes the weights for heat fields affected by SLHD
    !  all the rest is recomputed once again bellow.
    LLT_SLHD(4)=.FALSE.

    ! * Calculation of PCLA and PCLASLD:
    CALL LASCAW_CLA(YDSL,YDDYNA,KFLEV,&
     & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,&
     & IILA,ZZWH,ZZWHMAD,PDLAT,PDLAMAD,PKAPPAT,ZKHTURB(:,:,1),&
     & PSLD,PSLDW,P3DTW,&
     & PCLA(:,:,:,1),PCLAMAD(:,:,:,1),PCLASLT(:,:,:))

    ! * Calculation of PCLO and PCLOSLD:
    CALL LASCAW_CLO(YDDYNA,KFLEV,&
     & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLO(:,:,1),PDLOMAD(:,:,1),&
     & PKAPPAT,ZKHTURB(:,:,1),&
     & PCLO(:,:,:,1,1),PCLOMAD(:,:,:,1,1),PCLOSLT(:,:,:,1))
    CALL LASCAW_CLO(YDDYNA,KFLEV,&
     & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLO(:,:,2),PDLOMAD(:,:,2),&
     & PKAPPAT,ZKHTURB(:,:,1),&
     & PCLO(:,:,:,2,1),PCLOMAD(:,:,:,2,1),PCLOSLT(:,:,:,2))

  ENDIF

  ! Loop over all horiz. weights for 3Dturb (computing in addition
  !                         two sets with and without SLHD)
  DO JJ=1,KDIMK
    IF ((KSPLTHOI == 1).AND.(JJ == KDIMK)) THEN
      ! Bit specific case computing diffusive weights for physical tendencies.
      ! In this case SLHD weights are of no use.

      ! * Calculation of PCLA and PCLASLD:
      CALL LASCAW_CLA(YDSL,YDDYNA,KFLEV,&
       & KPROMB,KST,KPROF,LLT_PHYS,ZSLHDKMINH,&
       & IILA,ZZWH,ZZWHMAD,PDLAT,PDLAMAD,ZKHTURB(:,:,JJ),ZKHTURB(:,:,JJ),&
       & PSLD,PSLDW,P3DTW,&
       & ZCLA(:,:,:),PCLAMAD(:,:,:,JJ),PCLA(:,:,:,JJ))

      ! * Calculation of PCLO and PCLOSLD:
      CALL LASCAW_CLO(YDDYNA,KFLEV,&
       & KPROMB,KST,KPROF,LLT_PHYS,ZSLHDKMINH,PDLO(:,:,1),PDLOMAD(:,:,1),&
       & ZKHTURB(:,:,JJ),ZKHTURB(:,:,JJ),&
       & ZCLO(:,:,:,1),PCLOMAD(:,:,:,1,JJ),PCLO(:,:,:,1,JJ))
      CALL LASCAW_CLO(YDDYNA,KFLEV,&
       & KPROMB,KST,KPROF,LLT_PHYS,ZSLHDKMINH,PDLO(:,:,2),PDLOMAD(:,:,2),&
       & ZKHTURB(:,:,JJ),ZKHTURB(:,:,JJ),&
       & ZCLO(:,:,:,2),PCLOMAD(:,:,:,2,JJ),PCLO(:,:,:,2,JJ))

    ELSE

      IF (JJ == 1) THEN
        LLT_SLHD(4)=.FALSE.
      ELSE
        LLT_SLHD(4)=LL3DTURB
      ENDIF

      ! * Calculation of PCLA, PCLAMAD and PCLASLD:
      CALL LASCAW_CLA(YDSL,YDDYNA,KFLEV,&
       & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,&
       & IILA,ZZWH,ZZWHMAD,PDLAT,PDLAMAD,PKAPPA,ZKHTURB(:,:,JJ),&
       & PSLD,PSLDW,P3DTW,&
       & PCLA(:,:,:,JJ),PCLAMAD(:,:,:,JJ),PCLASLD(:,:,:,JJ))

      ! * Calculation of PCLO and PCLOSLD for central lat 1 and 2
      ! (linear int. only for lat 0 and 3)
      CALL LASCAW_CLO(YDDYNA,KFLEV,&
       & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLO(:,:,1),PDLOMAD(:,:,1),&
       & PKAPPA,ZKHTURB(:,:,JJ),&
       & PCLO(:,:,:,1,JJ),PCLOMAD(:,:,:,1,JJ),PCLOSLD(:,:,:,1,JJ))
      CALL LASCAW_CLO(YDDYNA,KFLEV,&
       & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLO(:,:,2),PDLOMAD(:,:,2),&
       & PKAPPA,ZKHTURB(:,:,JJ),&
       & PCLO(:,:,:,2,JJ),PCLOMAD(:,:,:,2,JJ),PCLOSLD(:,:,:,2,JJ))

    ENDIF

  ENDDO

  ! * Calculation of PVINTW and PVINTWSLD:

  ! Do update of quadratic weight for vertical interpolation
  IF (KWIS == 102) LLT_SLHD(2) = (HOISLTV /= 0._JPRB)

  IF (((KWIS == 102).AND.LSLTVWENO) .OR. (KWIS == 106)) THEN

    ! Set value for boundary offset
    IF (LDREGETA) THEN
      IBCLIM=0
    ELSE
      IBCLIM=1
    ENDIF

    ! WENO computation
    KNOWENO(KST:KPROF,1:KFLEV) = 0
    DO JJ=3,1,-1

      SELECT CASE (JJ)
        CASE (1)
          IILEV(KST:KPROF,1:KFLEV)=KLEV(KST:KPROF,1:KFLEV)
        CASE (2)
          IILEV(KST:KPROF,1:KFLEV)=MIN(IFLVM2-IBCLIM,KLEV(KST:KPROF,1:KFLEV)+1) 
          KNOWENO(KST:KPROF,1:KFLEV)=KNOWENO(KST:KPROF,1:KFLEV) &
           & + IILEV(KST:KPROF,1:KFLEV)-KLEV(KST:KPROF,1:KFLEV) - 1
        CASE (3)
          ! can't be  KSLEV-1 as there is no half level on -1
          IILEV(KST:KPROF,1:KFLEV)=MAX(IBCLIM,KLEV(KST:KPROF,1:KFLEV)-1)
          KNOWENO(KST:KPROF,1:KFLEV)=KNOWENO(KST:KPROF,1:KFLEV) &
           & + IILEV(KST:KPROF,1:KFLEV)-KLEV(KST:KPROF,1:KFLEV) + 1
        CASE DEFAULT
          CALL ABOR1(' LASCAW: WENO PROBLEM')
      END SELECT    
        
      CALL LASCAW_VINTW(YDDYNA,&
       & KPROMB,KFLEV,KST,KPROF,LLCOMADV,LLT_SLHD,LLSLVF,LDSLHDHEAT,ZSLHDKMINV_WENO,IILEV,&
       & PLEV,PDVER,PDVERMAD,PSTDDISW,PKAPPA,PKAPPAT,PVETA,PVCUICO,PVSLD,PVSLDW,PVSLVF,&
       & PVINTW(:,:,3*(JJ-1)+1),PVINTWMAD,PVINTWSLD,PVINTWSLT,&
       & PVINTWSLVF(:,:,3*(JJ-1)+1))

    ENDDO

    ! make sure it only keeps -1,0,+1 values
    IF ((MAXVAL(KNOWENO(KST:KPROF,1:KFLEV)) > 1+IBCLIM) .OR. &
     &  (MINVAL(KNOWENO(KST:KPROF,1:KFLEV)) <-1-IBCLIM))     &
     &  CALL ABOR1(' LASCAW: Something strange is happenig about level shifts.')

    ! C_k functions 
    IF (LDREGETA) THEN
      DO JLEV=1,KFLEV
        DO JROF=KST,KPROF
          ! regular mesh (LREGETA=.t. case)
          ! Note: This code doesn't seem to work for irregular vertical spacing.
          !       Hence the smart LWENOBC code is only allowed with LREGETA=t.
          PCW(JROF,JLEV,1)=0.10_JPRB*(2.0_JPRB+PDVER(JROF,JLEV))*(3.0_JPRB-PDVER(JROF,JLEV))   ! central
          PCW(JROF,JLEV,2)=0.05_JPRB*(2.0_JPRB+PDVER(JROF,JLEV))*(1.0_JPRB+PDVER(JROF,JLEV))   ! lower
          PCW(JROF,JLEV,3)=0.05_JPRB*(2.0_JPRB-PDVER(JROF,JLEV))*(3.0_JPRB-PDVER(JROF,JLEV))   ! upper
        ENDDO
      ENDDO
    ELSE
      DO JLEV=1,KFLEV
        DO JROF=KST,KPROF
          ! general form
          ILEVV=KLEV(JROF,JLEV)
          IF ((ILEVV > 1) .AND. (ILEVV < KFLEV-3)) THEN
            PCW(JROF,JLEV,1)=PGAMMA_WENO(ILEVV,1) &
             & *(PLEV(JROF,JLEV)-PVETA(ILEVV-1))*(PLEV(JROF,JLEV)-PVETA(ILEVV+4))  ! central
            PCW(JROF,JLEV,2)=PGAMMA_WENO(ILEVV,2) &
             & *(PLEV(JROF,JLEV)-PVETA(ILEVV-1))*(PLEV(JROF,JLEV)-PVETA(ILEVV  ))  ! lower
            PCW(JROF,JLEV,3)=PGAMMA_WENO(ILEVV,3) &
             & *(PLEV(JROF,JLEV)-PVETA(ILEVV+3))*(PLEV(JROF,JLEV)-PVETA(ILEVV+4))  ! upper
          ENDIF
        ENDDO
      ENDDO
    ENDIF

  ELSE

    ! All the other cases but WENO
    CALL LASCAW_VINTW(YDDYNA,&
     & KPROMB,KFLEV,KST,KPROF,LLCOMADV,LLT_SLHD,LLSLVF,LDSLHDHEAT,ZSLHDKMINV,KLEV,&
     & PLEV,PDVER,PDVERMAD,PSTDDISW,PKAPPA,PKAPPAT,PVETA,PVCUICO,PVSLD,PVSLDW,PVSLVF,&
     & PVINTW,PVINTWMAD,PVINTWSLD,PVINTWSLT,PVINTWSLVF)

  ENDIF

  IF (KWIS == 104) THEN
    DO JLEV=1,KFLEV
      ! * Calculation of PHVW:
      DO JROF=KST,KPROF
        ZDVER=PDVER(JROF,JLEV)
        PHVW(JROF,JLEV,1)=FHLO1(ZDVER)
        PHVW(JROF,JLEV,2)=FHLO2(ZDVER)
        PHVW(JROF,JLEV,3)=FHLO3(ZDVER)
        PHVW(JROF,JLEV,4)=FHLO4(ZDVER)
      ENDDO
      ! * Calculation of PVDERW:
      DO JROF=KST,KPROF
        ILEVV=KLEV(JROF,JLEV)
        ZNUM=PVETA(ILEVV+2)-PVETA(ILEVV+1)
        ZDEN1=0.5_JPRB*(PVETA(ILEVV+2)-PVETA(ILEVV))
        ZDEN2=0.5_JPRB*(PVETA(ILEVV+3)-PVETA(ILEVV+1))
        IF(ILEVV >= 1.AND.ILEVV <= KFLEV-3) THEN
          PVDERW(JROF,JLEV,1,1)=0.5_JPRB*ZNUM/ZDEN1
          PVDERW(JROF,JLEV,2,1)=0.5_JPRB*ZNUM/ZDEN1
          PVDERW(JROF,JLEV,1,2)=0.5_JPRB*ZNUM/ZDEN2
          PVDERW(JROF,JLEV,2,2)=0.5_JPRB*ZNUM/ZDEN2
        ELSEIF (ILEVV == 0) THEN
          PVDERW(JROF,JLEV,1,1)=0.0_JPRB
          PVDERW(JROF,JLEV,2,1)=ZNUM/ZDEN1
          PVDERW(JROF,JLEV,1,2)=0.5_JPRB*ZNUM/ZDEN2
          PVDERW(JROF,JLEV,2,2)=0.5_JPRB*ZNUM/ZDEN2
        ELSEIF (ILEVV == KFLEV-2) THEN
          PVDERW(JROF,JLEV,1,1)=0.5_JPRB*ZNUM/ZDEN1
          PVDERW(JROF,JLEV,2,1)=0.5_JPRB*ZNUM/ZDEN1
          PVDERW(JROF,JLEV,1,2)=ZNUM/ZDEN2
          PVDERW(JROF,JLEV,2,2)=0.0_JPRB
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  IF (KWIS == 105) THEN
    ! * Calculation of PVINTWS (weights for cubic spline interpolation).
    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        ILEVV=KLEV(JROF,JLEV)
        ZD2=PLEV(JROF,JLEV)-PVETA(ILEVV+1)
        PVINTWS(JROF,JLEV,1)=RFVV(4,ILEVV  ,1)+ZD2*( RFVV(4,ILEVV  ,2) +&
         & ZD2*(RFVV(4,ILEVV   ,3) + ZD2*RFVV(4,ILEVV  ,4) ) )  
        PVINTWS(JROF,JLEV,2)=RFVV(3,ILEVV+1,1)+ZD2*( RFVV(3,ILEVV+1,2) +&
         & ZD2*( RFVV(3,ILEVV+1,3) + ZD2*RFVV(3,ILEVV+1,4) ) )  
        PVINTWS(JROF,JLEV,3)=RFVV(2,ILEVV+2,1)+ZD2*( RFVV(2,ILEVV+2,2) +&
         & ZD2*( RFVV(2,ILEVV+2,3) + ZD2*RFVV(2,ILEVV+2,4) ) )  
        PVINTWS(JROF,JLEV,4)=RFVV(1,ILEVV+3,1)+ZD2*( RFVV(1,ILEVV+3,2) +&
         & ZD2*( RFVV(1,ILEVV+3,3) + ZD2*RFVV(1,ILEVV+3,4) ) )  
      ENDDO
    ENDDO
  ENDIF

ENDIF

!     ----------------------------------------------------------------

!*       2.    2D MODEL AND CASES IN THE 3D MODEL WHERE ONLY
!              2D INTERPOLATIONS ARE NEEDED.
!              ---------------------------------------------

!        2.01  Coordinates and weights for bilinear interpolations.

IF (KWIS == 201) THEN

  DO JLEV=1,KFLEV
    IF (KHOR == 0) ILEV=JLEV
    IF (KHOR == 1) ILEV=KFLDN

    ! * Calculation of linear weights, KL0.
!CDIR NODEP
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF

      ZLO1  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+1)
      ILO1(JROF)  =INT(ZLO1)
      PDLO(JROF,JLEV,1)=ZLO1-REAL(ILO1(JROF),JPRB)
      ZLO2  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+2)
      ILO2(JROF)  =INT(ZLO2)
      PDLO(JROF,JLEV,2)=ZLO2-REAL(ILO2(JROF),JPRB)

      KL0(JROF,JLEV,1)=IADDR(ILAV(JROF,JLEV)+1)+YDSL%NSLEXT(ILO1(JROF),ILAV(JROF,JLEV)+1)+YDSL%NASLB1*ILEV
      KL0(JROF,JLEV,2)=IADDR(ILAV(JROF,JLEV)+2)+YDSL%NSLEXT(ILO2(JROF),ILAV(JROF,JLEV)+2)+YDSL%NASLB1*ILEV

    ENDDO

    ! * Mask calculation for on-demand communications:
    IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
      DO JROF=KST,KPROF
        YDSL%MASK_SL2(KL0(JROF,JLEV,1):KL0(JROF,JLEV,1)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2):KL0(JROF,JLEV,2)+3)=1
      ENDDO
    ENDIF

  ENDDO

ENDIF

!        2.03  Coordinates and weights for 12 points interpolations.

IF (KWIS == 202 .OR. KWIS == 203) THEN

  DO JLEV=1,KFLEV
    IF (KHOR == 0) ILEV=JLEV
    IF (KHOR == 1) ILEV=KFLDN

    ! * Calculation of linear weights, KL0.
!CDIR NODEP
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF

      ! meridional interpolation: linear weights and input for cubic weights (LASCAW_CLA)
      ! general case 
      PDLAT(JROF,JLEV)=(PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+1))&
       & /(PLATI(ILAV(JROF,JLEV)+2)-PLATI(ILAV(JROF,JLEV)+1))  
      ZDA   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV))
      ZDB   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+1)
      ZDC   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+2)
      ZDD   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+3)

      IILA(JROF,JLEV)=ILAV(JROF,JLEV)
      ZZWH(JROF,1,JLEV)=(ZDA*ZDC)*ZDD*PIPI(ILAV(JROF,JLEV)+1,1)
      ZZWH(JROF,2,JLEV)=(ZDA*ZDB)*ZDD*PIPI(ILAV(JROF,JLEV)+1,2)
      ZZWH(JROF,3,JLEV)=(ZDA*ZDB)*ZDC*PIPI(ILAV(JROF,JLEV)+1,3)

      ! COMAD meridional interpolation 
      IF (LLCOMADH) THEN
        PDLAMAD(JROF,JLEV)=PDLAT(JROF,JLEV)*PSTDDISV(JROF,JLEV) + 0.5_JPRB * (1._JPRB-PSTDDISV(JROF,JLEV))
        ZDAMAD = ZDA-ZDB
        ZDDMAD = ZDD-ZDC
        ZDBMAD = ZDB*PSTDDISV(JROF,JLEV) +0.5_JPRB*(PLATI(ILAV(JROF,JLEV)+2)-PLATI(ILAV(JROF,JLEV)+1)) *&
         & (1._JPRB-PSTDDISV(JROF,JLEV))
        ZDCMAD = ZDC*PSTDDISV(JROF,JLEV) -0.5_JPRB*(PLATI(ILAV(JROF,JLEV)+2)-PLATI(ILAV(JROF,JLEV)+1)) *&
         & (1._JPRB-PSTDDISV(JROF,JLEV))
        ZDAMAD = ZDAMAD + ZDBMAD
        ZDDMAD = ZDDMAD + ZDCMAD
        IILA(JROF,JLEV)=ILAV(JROF,JLEV)
        ZZWHMAD(JROF,1,JLEV)=(ZDAMAD*ZDCMAD)*ZDDMAD*PIPI(ILAV(JROF,JLEV)+1,1)
        ZZWHMAD(JROF,2,JLEV)=(ZDAMAD*ZDBMAD)*ZDDMAD*PIPI(ILAV(JROF,JLEV)+1,2)
        ZZWHMAD(JROF,3,JLEV)=(ZDAMAD*ZDBMAD)*ZDCMAD*PIPI(ILAV(JROF,JLEV)+1,3)
      ELSE
        PDLAMAD(JROF,JLEV)=PDLAT(JROF,JLEV)
        ZZWHMAD(JROF,1,JLEV)=ZZWH(JROF,1,JLEV)
        ZZWHMAD(JROF,2,JLEV)=ZZWH(JROF,2,JLEV)
        ZZWHMAD(JROF,3,JLEV)=ZZWH(JROF,3,JLEV)
      ENDIF

      ! zonal interpolation: linear weights for 4 lat. cicles
      ! as the grid is regular in the zonal direction,
      ! the cubic weight computation does not need 
      ! other input than linear weights (LASCAW_CLO)
      ! general case
      ZLO   =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)  )
      ILO(JROF)   =INT(ZLO )
      PDLO(JROF,JLEV,0)=ZLO -REAL(ILO(JROF) ,JPRB)
      ZLO1  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+1)
      ILO1(JROF)  =INT(ZLO1)
      PDLO(JROF,JLEV,1)=ZLO1-REAL(ILO1(JROF),JPRB)
      ZLO2  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+2)
      ILO2(JROF)  =INT(ZLO2)
      PDLO(JROF,JLEV,2)=ZLO2-REAL(ILO2(JROF),JPRB)
      ZLO3  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+3)
      ILO3(JROF)  =INT(ZLO3)
      PDLO(JROF,JLEV,3)=ZLO3-REAL(ILO3(JROF),JPRB)

      ! COMAD zonal interpolation 
      IF (LLCOMADH) THEN
        PDLOMAD(JROF,JLEV,0)=PDLO(JROF,JLEV,0)*PSTDDISU(JROF,JLEV) + 0.5_JPRB * (1._JPRB-PSTDDISU(JROF,JLEV))
        PDLOMAD(JROF,JLEV,1)=PDLO(JROF,JLEV,1)*PSTDDISU(JROF,JLEV) + 0.5_JPRB * (1._JPRB-PSTDDISU(JROF,JLEV))
        PDLOMAD(JROF,JLEV,2)=PDLO(JROF,JLEV,2)*PSTDDISU(JROF,JLEV) + 0.5_JPRB * (1._JPRB-PSTDDISU(JROF,JLEV))
        PDLOMAD(JROF,JLEV,3)=PDLO(JROF,JLEV,3)*PSTDDISU(JROF,JLEV) + 0.5_JPRB * (1._JPRB-PSTDDISU(JROF,JLEV))
      ELSE
        PDLOMAD(JROF,JLEV,0)= PDLO(JROF,JLEV,0)
        PDLOMAD(JROF,JLEV,1)= PDLO(JROF,JLEV,1)
        PDLOMAD(JROF,JLEV,2)= PDLO(JROF,JLEV,2)
        PDLOMAD(JROF,JLEV,3)= PDLO(JROF,JLEV,3)
      ENDIF

      ILO(JROF)=IADDR(ILAV(JROF,JLEV)  )+YDSL%NSLEXT(ILO(JROF) ,ILAV(JROF,JLEV)  )
      ILO1(JROF)=IADDR(ILAV(JROF,JLEV)+1)+YDSL%NSLEXT(ILO1(JROF),ILAV(JROF,JLEV)+1)
      ILO2(JROF)=IADDR(ILAV(JROF,JLEV)+2)+YDSL%NSLEXT(ILO2(JROF),ILAV(JROF,JLEV)+2)
      ILO3(JROF)=IADDR(ILAV(JROF,JLEV)+3)+YDSL%NSLEXT(ILO3(JROF),ILAV(JROF,JLEV)+3)

      KL0(JROF,JLEV,0)=ILO(JROF)+YDSL%NASLB1*ILEV
      KL0(JROF,JLEV,1)=ILO1(JROF)+YDSL%NASLB1*ILEV
      KL0(JROF,JLEV,2)=ILO2(JROF)+YDSL%NASLB1*ILEV
      KL0(JROF,JLEV,3)=ILO3(JROF)+YDSL%NASLB1*ILEV

    ENDDO

    ! * Mask calculation for on-demand communications:
    IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
      DO JROF=KST,KPROF
        YDSL%MASK_SL2(KL0(JROF,JLEV,0):KL0(JROF,JLEV,0)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1):KL0(JROF,JLEV,1)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2):KL0(JROF,JLEV,2)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,3):KL0(JROF,JLEV,3)+3)=1
      ENDDO
    ENDIF

  ENDDO

  ! * Calculation of PCLA and PCLASLD:
  CALL LASCAW_CLA(YDSL,YDDYNA,KFLEV,&
   & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,&
   & IILA,ZZWH,ZZWHMAD,PDLAT,PDLAMAD,PKAPPA,ZKHTURB(:,:,1),&
   & PSLD,PSLDW,P3DTW,&
   & PCLA(:,:,:,1),PCLAMAD(:,:,:,1),PCLASLD(:,:,:,1))

  ! * Calculation of PCLO and PCLOSLD:
  CALL LASCAW_CLO(YDDYNA,KFLEV,&
   & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLO(:,:,1),PDLOMAD(:,:,1),&
   & PKAPPA,ZKHTURB(:,:,1),&
   & PCLO(:,:,:,1,1),PCLOMAD(:,:,:,1,1),PCLOSLD(:,:,:,1,1))
  CALL LASCAW_CLO(YDDYNA,KFLEV,&
   & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLO(:,:,2),PDLOMAD(:,:,2),&
   & PKAPPA,ZKHTURB(:,:,1),&
   & PCLO(:,:,:,2,1),PCLOMAD(:,:,:,2,1),PCLOSLD(:,:,:,2,1))

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('LASCAW',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAW
