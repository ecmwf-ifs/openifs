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
!      H Petithomme (Dec 2020): optimisation and simplification, COMAD bugfix
!      R. El Khatib 03-Nov-2021 Reduce cache misses on MASK_SL2 + vectorization
!      F. Vana  14-Mar-2024: fixed CY49R1
! End Modifications
!     ------------------------------------------------------------------

USE YOMVSPLIP , ONLY : TVSPLIP
USE PARKIND1  , ONLY : JPIA,JPIM, JPRB, JPRD
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
INTEGER(KIND=JPIM),INTENT(OUT) :: KL0(KPROMB,KFLEV,0:3)
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

INTEGER(KIND=JPIA) :: ILEV64,IL
INTEGER(KIND=JPIM) :: ILEV,IOFF,JOFF,INP
INTEGER(KIND=JPIM) :: IADDR(YDSL%NDGSAH:YDSL%NDGENH)
INTEGER(KIND=JPIM) :: IFLVM2,JLAT,JLEV,JROF,IBCLIM,JJ,IZLATV
INTEGER(KIND=JPIM) :: ILAV(KPROMB,KFLEV)
INTEGER(KIND=JPIM),ALLOCATABLE :: IILEV(:,:)
INTEGER(KIND=JPIM),TARGET :: ILO(4*KPROMB)
INTEGER(KIND=JPIM),CONTIGUOUS,POINTER :: ILO1(:),ILO2(:),ILO3(:)
INTEGER(KIND=JPIM) :: ILO0IK(KST:KPROF,KFLEV), ILO1IK(KST:KPROF,KFLEV)
INTEGER(KIND=JPIM) :: ILO2IK(KST:KPROF,KFLEV), ILO3IK(KST:KPROF,KFLEV)

REAL(KIND=JPRB) :: ZDAMAD, ZDBMAD, ZDCMAD, ZDDMAD

REAL(KIND=JPRB) :: ZKHTURB(KPROMB,KFLEV,KDIMK)
LOGICAL         :: LLT_SLHD(4),LLT_PHYS(4),LLSLHD,LLSLHDQUAD,LLSLHD_OLD,LL3DTURB,LLSLVF
LOGICAL         :: LLCOMAD,LLCOMADH,LLCOMADV
REAL(KIND=JPRB) :: PD,ZDA,ZDB,ZDC,ZDD,ZDVER,ZFAC,ZLO,ZLO1,ZLO2,ZLO3,ZEPS
REAL(KIND=JPRB) :: ZSLHDKMINH,ZSLHDKMINV,ZW1,ZW2,Z1,Z2,Z3,Z4
REAL(KIND=JPRB), PARAMETER :: ZSLHDKMINV_WENO=0._JPRB   ! WENO only runs with Lagrangian cubic!!!

REAL(KIND=JPRB) :: ZCLA(KPROMB,KFLEV,3),ZCLO(KPROMB,KFLEV,3)
REAL(KIND=JPRB) :: ZDUM(KPROMB,KFLEV,3) ! output field to be ditched

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
#include "lascaw_claturb.intfb.h"
#include "lascaw_cloturb.intfb.h"
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
LLCOMAD =LDCOMAD.AND.(KWIS==103.OR.KWIS==104.OR.KWIS==105.OR.KWIS==106.OR.KWIS==203)
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
  IADDR(JLAT)=KSTABUF(JLAT)-YDSL%NASLB1*KFLDN
ENDDO

IF (LL3DTURB) THEN
  ! copy from kappa
  ZKHTURB(KST:KPROF,1:KFLEV,2)=PKAPPAM(KST:KPROF,1:KFLEV)
  ZKHTURB(KST:KPROF,1:KFLEV,3)=PKAPPAH(KST:KPROF,1:KFLEV)
ENDIF

! In this case ZKHTURB is used as KAPPA and is set to static mode with maximum diffusion.
IF (KSPLTHOI == 1) ZKHTURB(KST:KPROF,1:KFLEV,KDIMK)=1._JPRB

IFLVM2=KFLEV-2

! optim: general indexing of longitudinal positions on latitudes using pointers on
! adjacent positions so that ilo is contiguously filled over kst:krpof+3*np
! please bear in mind that pointers do overlap (by kst-1), but not final positions in ilo
INP = KPROF-KST+1
ILO1 => ILO(INP+1:)
ILO2 => ILO(2*INP+1:)
ILO3 => ILO(3*INP+1:)

!     ------------------------------------------------------------------

!*       1.    3D MODEL.
!              ---------

!*   distance for horizontal linear interpolations in latitude (common to all options)

! optim: enhanced vectorization by way of intermediate variable (lowers gathers)
DO JLEV=1,KFLEV
  DO JROF=KST,KPROF
    ! * Calculation of linear weights
    IZLATV=INT(P4JP*(PIS2-PLAT(JROF,JLEV))+0.75_JPRB+ZEPS)-YDSL%NFRSTLOFF
    ILAV(JROF,JLEV)=IZLATV+NINT(SIGN(0.5_JPRB,PLATI(IZLATV)-PLAT(JROF,JLEV)+ZEPS)-1.5_JPRB)

    Z1 = PLATI(ILAV(JROF,JLEV)+1)
    PDLAT(JROF,JLEV)=(PLAT(JROF,JLEV)-Z1)/(PLATI(ILAV(JROF,JLEV)+2)-Z1+ZEPS)
  ENDDO
ENDDO

!        1.01  Coordinates and weights for trilinear interpolations.
!              horizontal 12 points + vertical cubic + 32 points interpolations
!              Optionally, Hermite cubic, cubic B-spline or WENO vertical interpolations weights are computed.

! note: IOFF varies or not, depending on KHOR
IOFF=YDSL%NASLB1*KFLDN

IF (101 <= KWIS.AND.KWIS <= 106) THEN
  ZFAC=PVRLEVX/(PVETA(KFLEV+1)-PVETA(0))

  DO JLEV=1,KFLEV
    ! vertical interpolation: linear weights
    ! the cubic weight computation are done in LASCAW_VINTW
    !CDIR NODEP
    !DIR$ PREFERVECTOR
    !DIR$ IVDEP
    DO JROF=KST,KPROF
      KLEV(JROF,JLEV)=KVAUT(INT(PLEV(JROF,JLEV)*ZFAC))-1
      IF (KLEV(JROF,JLEV) < IFLVM2 .AND. PLEV(JROF,JLEV) > PVETA(KLEV(JROF,JLEV)+2)) THEN
        KLEV(JROF,JLEV)=KLEV(JROF,JLEV)+1
      ENDIF
      Z1 = PVETA(KLEV(JROF,JLEV)+1)
      PDVER(JROF,JLEV)=(PLEV(JROF,JLEV)-Z1)/(PVETA(KLEV(JROF,JLEV)+2)-Z1)
    ! * Calculation of linear weights, KL0, KLH0.
    ! zonal interpolation: linear weights for 4 lat. circles
    ! as the grid is regular in the zonal direction, the cubic weight computation does not need
    ! other input than linear weights (LASCAW_CLO)
      ZLO   =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV))
      ILO(JROF)   =INT(ZLO )
      PDLO(JROF,JLEV,0)=ZLO -ILO(JROF)
      ZLO1  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+1)
      ILO1(JROF)  =INT(ZLO1)
      PDLO(JROF,JLEV,1)=ZLO1-ILO1(JROF)
      ZLO2  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+2)
      ILO2(JROF)  =INT(ZLO2)
      PDLO(JROF,JLEV,2)=ZLO2-ILO2(JROF)
      ZLO3  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+3)
      ILO3(JROF)  =INT(ZLO3)
      PDLO(JROF,JLEV,3)=ZLO3-ILO3(JROF)
    ENDDO

    ! optim: does not vectorize, use 64-bit indexing (IL)
    DO JROF=KST,KPROF
      IL=ILAV(JROF,JLEV)
      ILO(JROF)=IADDR(ILAV(JROF,JLEV))+YDSL%NSLEXT(ILO(JROF),IL)
      ILO1(JROF)=IADDR(ILAV(JROF,JLEV)+1)+YDSL%NSLEXT(ILO1(JROF),IL+1)
      ILO2(JROF)=IADDR(ILAV(JROF,JLEV)+2)+YDSL%NSLEXT(ILO2(JROF),IL+2)
      ILO3(JROF)=IADDR(ILAV(JROF,JLEV)+3)+YDSL%NSLEXT(ILO3(JROF),IL+3)
    ENDDO

    DO JROF=KST,KPROF
      ILO0IK(JROF,JLEV)=ILO (JROF)
      ILO1IK(JROF,JLEV)=ILO1(JROF)
      ILO2IK(JROF,JLEV)=ILO2(JROF)
      ILO3IK(JROF,JLEV)=ILO3(JROF)
    ENDDO

    IF(KWIS == 101) THEN
      ! note: case 101 only set linear weights, so only requires central lats
      DO JROF=KST,KPROF
        JOFF=YDSL%NASLB1*KLEV(JROF,JLEV)
        KL0(JROF,JLEV,1)=ILO1(JROF)+JOFF
        KL0(JROF,JLEV,2)=ILO2(JROF)+JOFF
      ENDDO
    ELSE
      DO JROF=KST,KPROF
        JOFF=YDSL%NASLB1*KLEV(JROF,JLEV)
        KL0(JROF,JLEV,0)=ILO(JROF)+JOFF
        KL0(JROF,JLEV,1)=ILO1(JROF)+JOFF
        KL0(JROF,JLEV,2)=ILO2(JROF)+JOFF
        KL0(JROF,JLEV,3)=ILO3(JROF)+JOFF
      ENDDO

      ! note: mind the offsets IOFF/JOFF
      IF(KHOR == 0) IOFF=YDSL%NASLB1*JLEV
      KLH0(KST:KPROF,JLEV,0)=ILO(KST:KPROF)+IOFF
      KLH0(KST:KPROF,JLEV,1)=ILO1(KST:KPROF)+IOFF
      KLH0(KST:KPROF,JLEV,2)=ILO2(KST:KPROF)+IOFF
      KLH0(KST:KPROF,JLEV,3)=ILO3(KST:KPROF)+IOFF
    ENDIF
  ENDDO

  ! * Mask calculation for on-demand communications:
  IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE) THEN
    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        YDSL%MASK_SL2(ILO1IK(JROF,JLEV):ILO1IK(JROF,JLEV)+3)=1
      ENDDO
    ENDDO
    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        YDSL%MASK_SL2(ILO2IK(JROF,JLEV):ILO2IK(JROF,JLEV)+3)=1
      ENDDO
    ENDDO
    IF (KWIS /= 101) THEN
      DO JLEV=1,KFLEV
        DO JROF=KST,KPROF
          YDSL%MASK_SL2(ILO0IK(JROF,JLEV):ILO0IK(JROF,JLEV)+3)=1
        ENDDO
      ENDDO
      DO JLEV=1,KFLEV
        DO JROF=KST,KPROF
          YDSL%MASK_SL2(ILO3IK(JROF,JLEV):ILO3IK(JROF,JLEV)+3)=1
        ENDDO
      ENDDO
    ENDIF
  ENDIF
ENDIF

IF (102 <= KWIS.AND.KWIS <= 106) THEN
  DO JLEV=1,KFLEV
    DO JROF=KST,KPROF
      ! meridional interpolation: linear weights and input for cubic weights (LASCAW_CLA)
      ! general case 
      ZDA   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV))
      ZDB   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+1)
      ZDC   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+2)
      ZDD   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+3)
      Z2 = ZDA*ZDB
      ! here we need to keep this original array separated from the weight 
      ! PCLA as that might be modified by every call of LASCAW_CLA
      ZCLA(JROF,JLEV,1)=(ZDA*ZDC)*ZDD*PIPI(ILAV(JROF,JLEV)+1,1)
      ZCLA(JROF,JLEV,2)=Z2*ZDD*PIPI(ILAV(JROF,JLEV)+1,2)
      ZCLA(JROF,JLEV,3)=Z2*ZDC*PIPI(ILAV(JROF,JLEV)+1,3)

      IF (LLCOMADH) THEN
        Z1 = 1._JPRB/PSTDDISV(JROF,JLEV)
        Z4 = 0.5_JPRB*(PLATI(ILAV(JROF,JLEV)+2)-PLATI(ILAV(JROF,JLEV)+1))*(Z1-1._JPRB)
        ZDAMAD = (ZDA-ZDB)*Z1
        ZDDMAD = (ZDD-ZDC)*Z1
        ZDBMAD = ZDB+Z4
        ZDCMAD = ZDC-Z4
        ZDAMAD = ZDAMAD + ZDBMAD
        ZDDMAD = ZDDMAD + ZDCMAD
        Z2 = ZDAMAD*ZDBMAD
        Z3 = PSTDDISV(JROF,JLEV)**3
        PCLAMAD(JROF,JLEV,1,1)=ZDAMAD*ZDCMAD*ZDDMAD*PIPI(ILAV(JROF,JLEV)+1,1)*Z3
        PCLAMAD(JROF,JLEV,2,1)=Z2*ZDDMAD*PIPI(ILAV(JROF,JLEV)+1,2)*Z3
        PCLAMAD(JROF,JLEV,3,1)=Z2*ZDCMAD*PIPI(ILAV(JROF,JLEV)+1,3)*Z3
      ELSE
        PCLAMAD(JROF,JLEV,1,1)=ZCLA(JROF,JLEV,1)
        PCLAMAD(JROF,JLEV,2,1)=ZCLA(JROF,JLEV,2)
        PCLAMAD(JROF,JLEV,3,1)=ZCLA(JROF,JLEV,3)
      ENDIF
    ENDDO

    ! COMAD zonal linear weights
    IF (LLCOMADH) THEN
      DO JROF=KST,KPROF
          PDLAMAD(JROF,JLEV)=0.5_JPRB+(PDLAT(JROF,JLEV)-0.5_JPRB)*PSTDDISV(JROF,JLEV)
          PDLOMAD(JROF,JLEV,0)=0.5_JPRB+(PDLO(JROF,JLEV,0)-0.5_JPRB)*PSTDDISU(JROF,JLEV)
          PDLOMAD(JROF,JLEV,1)=0.5_JPRB+(PDLO(JROF,JLEV,1)-0.5_JPRB)*PSTDDISU(JROF,JLEV)
          PDLOMAD(JROF,JLEV,2)=0.5_JPRB+(PDLO(JROF,JLEV,2)-0.5_JPRB)*PSTDDISU(JROF,JLEV)
          PDLOMAD(JROF,JLEV,3)=0.5_JPRB+(PDLO(JROF,JLEV,3)-0.5_JPRB)*PSTDDISU(JROF,JLEV)
      ENDDO
    ELSE
      DO JROF=KST,KPROF
        PDLAMAD(JROF,JLEV)=PDLAT(JROF,JLEV)
        PDLOMAD(JROF,JLEV,0)=PDLO(JROF,JLEV,0)
        PDLOMAD(JROF,JLEV,1)=PDLO(JROF,JLEV,1)
        PDLOMAD(JROF,JLEV,2)=PDLO(JROF,JLEV,2)
        PDLOMAD(JROF,JLEV,3)=PDLO(JROF,JLEV,3)
      ENDDO
    ENDIF

    ! COMAD vertical interpolation 
    IF (LLCOMADV) THEN
      DO JROF=KST,KPROF
        PDVERMAD(JROF,JLEV)=0.5_JPRB+(PDVER(JROF,JLEV)-0.5_JPRB)*PSTDDISW(JROF,JLEV)
      ENDDO
    ELSE
      DO JROF=KST,KPROF
        PDVERMAD(JROF,JLEV)=PDVER(JROF,JLEV)
      ENDDO
    ENDIF
  ENDDO

  IF (LLSLHD.AND.LDSLHDHEAT) THEN
    ! Computes the weights for heat fields affected by SLHD
    !  all the rest is recomputed once again below.

    ! * Calculation of PCLA and PCLASLD:
    CALL LASCAW_CLA(YDSL,YDDYNA,KFLEV,KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,&
     & ILAV,ZCLA,PDLAT,PKAPPAT,&
     & PSLD,PSLDW,PCLA(:,:,:,1),PCLASLT)

    ! * Calculation of PCLO and PCLOSLD:
    CALL LASCAW_CLO(YDDYNA,KFLEV,&
     & KPROMB,KST,KPROF,LLCOMADH,LLT_SLHD,ZSLHDKMINH,&
     & PDLO(:,:,1),PDLOMAD(:,:,1),PKAPPAT,&
     & PCLO(:,:,:,1,1),PCLOMAD(:,:,:,1,1),PCLOSLT(:,:,:,1))
    CALL LASCAW_CLO(YDDYNA,KFLEV,&
     & KPROMB,KST,KPROF,LLCOMADH,LLT_SLHD,ZSLHDKMINH,&
     & PDLO(:,:,2),PDLOMAD(:,:,2),PKAPPAT,&
     & PCLO(:,:,:,2,1),PCLOMAD(:,:,:,2,1),PCLOSLT(:,:,:,2))
  ENDIF

  ! Loop over all horiz. weights for 3Dturb (computing in addition
  !                         two sets with and without SLHD)
  DO JJ=1,KDIMK
    ! * Calculation of PCLA and PCLASLD:
    IF(JJ > 1) THEN
      PCLAMAD(KST:KPROF,1:KFLEV,1:3,JJ) = PCLAMAD(KST:KPROF,1:KFLEV,1:3,1)
    ENDIF

    IF ((KSPLTHOI == 1).AND.(JJ == KDIMK)) THEN
      ! Bit specific case computing diffusive weights for physical tendencies.
      ! In this case SLHD weights are of no use.

      ! warning: PCLA must not be modified in LASCAW_CLA as PCLA, but as PCLASLD
      CALL LASCAW_CLA(YDSL,YDDYNA,KFLEV,KPROMB,KST,KPROF,LLT_PHYS,ZSLHDKMINH, &
       & ILAV,ZCLA,PDLAT,&
       & ZKHTURB(:,:,JJ),PSLD,PSLDW,ZDUM,PCLA(:,:,:,JJ))

      ! * Calculation of PCLO (as PCLOSLD) and PCLOMAD:
      CALL LASCAW_CLO(YDDYNA,KFLEV,&
       & KPROMB,KST,KPROF,LLCOMADH,LLT_PHYS,ZSLHDKMINH,&
       & PDLO(:,:,1),PDLOMAD(:,:,1),ZKHTURB(:,:,JJ),&
       & ZCLO,PCLOMAD(:,:,:,1,JJ),PCLO(:,:,:,1,JJ))
      CALL LASCAW_CLO(YDDYNA,KFLEV,&
       & KPROMB,KST,KPROF,LLCOMADH,LLT_PHYS,ZSLHDKMINH,&
       & PDLO(:,:,2),PDLOMAD(:,:,2),ZKHTURB(:,:,JJ),&
       & ZCLO,PCLOMAD(:,:,:,2,JJ),PCLO(:,:,:,2,JJ))
    ELSE
      ! * Calculation of PCLA, PCLAMAD and PCLASLD:
      CALL LASCAW_CLA(YDSL,YDDYNA,KFLEV,KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,&
        & ILAV,ZCLA,PDLAT,PKAPPA,&
       & PSLD,PSLDW,PCLA(:,:,:,JJ),PCLASLD(:,:,:,JJ))

      ! * Calculation of PCLO and PCLOSLD for central lat 1 and 2
      ! (linear int. only for lat 0 and 3)
      CALL LASCAW_CLO(YDDYNA,KFLEV,&
       & KPROMB,KST,KPROF,LLCOMADH,LLT_SLHD,ZSLHDKMINH,&
       & PDLO(:,:,1),PDLOMAD(:,:,1),PKAPPA,&
       & PCLO(:,:,:,1,JJ),PCLOMAD(:,:,:,1,JJ),PCLOSLD(:,:,:,1,JJ))
      CALL LASCAW_CLO(YDDYNA,KFLEV,&
       & KPROMB,KST,KPROF,LLCOMADH,LLT_SLHD,ZSLHDKMINH,&
       & PDLO(:,:,2),PDLOMAD(:,:,2),PKAPPA,&
       & PCLO(:,:,:,2,JJ),PCLOMAD(:,:,:,2,JJ),PCLOSLD(:,:,:,2,JJ))

      IF (JJ > 1.AND.LL3DTURB) THEN
        CALL LASCAW_CLATURB(YDSL,KFLEV,KPROMB,KST,KPROF,ILAV,ZKHTURB(:,:,JJ),P3DTW,&
         & PCLA(:,:,:,JJ),PCLASLD(:,:,:,JJ))
        CALL LASCAW_CLOTURB(KFLEV,KPROMB,KST,KPROF,ZKHTURB(:,:,JJ),PCLO(:,:,:,1,JJ),&
         & PCLOSLD(:,:,:,1,JJ))
        CALL LASCAW_CLOTURB(KFLEV,KPROMB,KST,KPROF,ZKHTURB(:,:,JJ),PCLO(:,:,:,2,JJ),&
         & PCLOSLD(:,:,:,2,JJ))
      ENDIF
    ENDIF
  ENDDO

  ! * Calculation of PVINTW and PVINTWSLD:

  ! Do update of quadratic weight for vertical interpolation
  IF (KWIS == 102) LLT_SLHD(2) = HOISLTV /= 0._JPRB

  IF (KWIS == 102.AND.LSLTVWENO.OR.KWIS == 106) THEN
    ALLOCATE(IILEV(KPROMB,KFLEV))

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
          DO JLEV=1,KFLEV
            !dir$ ivdep
            DO JROF=KST,KPROF
              KNOWENO(JROF,JLEV)=KNOWENO(JROF,JLEV)+IILEV(JROF,JLEV)-KLEV(JROF,JLEV)-1
            ENDDO
          ENDDO
        CASE (3)
          ! can't be  KSLEV-1 as there is no half level on -1
          IILEV(KST:KPROF,1:KFLEV)=MAX(IBCLIM,KLEV(KST:KPROF,1:KFLEV)-1)
          DO JLEV=1,KFLEV
            !dir$ ivdep
            DO JROF=KST,KPROF
              KNOWENO(JROF,JLEV)=KNOWENO(JROF,JLEV)+IILEV(JROF,JLEV)-KLEV(JROF,JLEV)+1
            ENDDO
          ENDDO
        CASE DEFAULT
          CALL ABOR1(' LASCAW: WENO PROBLEM')
      END SELECT    
        
      CALL LASCAW_VINTW(YDDYNA,&
       & KPROMB,KFLEV,KST,KPROF,LLCOMADV,LLT_SLHD,LLSLVF,LDSLHDHEAT,ZSLHDKMINV_WENO,IILEV,&
       & PLEV,PDVER,PDVERMAD,PSTDDISW,PKAPPA,PKAPPAT,PVETA,PVCUICO,PVSLD,PVSLDW,PVSLVF,&
       & PVINTW(:,:,3*(JJ-1)+1),PVINTWMAD,PVINTWSLD,PVINTWSLT,&
       & PVINTWSLVF(:,:,3*(JJ-1)+1))

    ENDDO

    DEALLOCATE(IILEV)

    ! make sure it only keeps -1,0,+1 values
    IF (ANY(ABS(KNOWENO(KST:KPROF,1:KFLEV)) > 1+IBCLIM))&
     &  CALL ABOR1(' LASCAW: Something strange is happening about level shifts.')

    ! C_k functions 
    IF (LDREGETA) THEN
      DO JLEV=1,KFLEV
        !dir$ ivdep
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
        ! optim: does not vectorize (yet ?), but use of interm. variables lowers gathers
        DO JROF=KST,KPROF
          ! general form
          ILEV=KLEV(JROF,JLEV)
          IF (ILEV > 1.AND.ILEV < KFLEV-3) THEN
            Z1 = PLEV(JROF,JLEV)-PVETA(ILEV-1)
            Z4 = PLEV(JROF,JLEV)-PVETA(ILEV+4)
            PCW(JROF,JLEV,1)=PGAMMA_WENO(ILEV,1)*Z1*Z4 ! central
            PCW(JROF,JLEV,2)=PGAMMA_WENO(ILEV,2)*Z1*(PLEV(JROF,JLEV)-PVETA(ILEV)) ! lower
            PCW(JROF,JLEV,3)=PGAMMA_WENO(ILEV,3)*(PLEV(JROF,JLEV)-PVETA(ILEV+3))*Z4 ! upper
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
      ! optim: now vectorizes with help of intermediate variables (lowers gathers)
      DO JROF=KST,KPROF
        ILEV=KLEV(JROF,JLEV)
        Z1 = PVETA(ILEV+1)
        Z2 = PVETA(ILEV+2)
        ZW1=(Z2-Z1)/(Z2-PVETA(ILEV))
        ZW2=(Z2-Z1)/(PVETA(ILEV+3)-Z1)
        IF(ILEV >= 1.AND.ILEV <= KFLEV-3) THEN
          PVDERW(JROF,JLEV,1,1)=ZW1
          PVDERW(JROF,JLEV,2,1)=ZW1
          PVDERW(JROF,JLEV,1,2)=ZW2
          PVDERW(JROF,JLEV,2,2)=ZW2
        ELSEIF (ILEV == 0) THEN
          PVDERW(JROF,JLEV,1,1)=0.0_JPRB
          PVDERW(JROF,JLEV,2,1)=2._jprb*ZW1
          PVDERW(JROF,JLEV,1,2)=ZW2
          PVDERW(JROF,JLEV,2,2)=ZW2
        ELSEIF (ILEV == KFLEV-2) THEN
          PVDERW(JROF,JLEV,1,1)=ZW1
          PVDERW(JROF,JLEV,2,1)=ZW1
          PVDERW(JROF,JLEV,1,2)=2._jprb*ZW2
          PVDERW(JROF,JLEV,2,2)=0.0_JPRB
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  IF (KWIS == 105) THEN
    ! * Calculation of PVINTWS (weights for cubic spline interpolation).
    DO JLEV=1,KFLEV
      ! optim: does not vectorize
      DO JROF=KST,KPROF
        ILEV64=KLEV(JROF,JLEV)
        Z1=PLEV(JROF,JLEV)-PVETA(ILEV64+1)
        PVINTWS(JROF,JLEV,1)=RFVV(4,ILEV64,1)+Z1*(RFVV(4,ILEV64,2) +&
         & Z1*(RFVV(4,ILEV64,3) + Z1*RFVV(4,ILEV64,4)))
        PVINTWS(JROF,JLEV,2)=RFVV(3,ILEV64+1,1)+Z1*(RFVV(3,ILEV64+1,2) +&
         & Z1*(RFVV(3,ILEV64+1,3) + Z1*RFVV(3,ILEV64+1,4)))
        PVINTWS(JROF,JLEV,3)=RFVV(2,ILEV64+2,1)+Z1*(RFVV(2,ILEV64+2,2) +&
         & Z1*(RFVV(2,ILEV64+2,3) + Z1*RFVV(2,ILEV64+2,4)))
        PVINTWS(JROF,JLEV,4)=RFVV(1,ILEV64+3,1)+Z1*(RFVV(1,ILEV64+3,2) +&
         & Z1*(RFVV(1,ILEV64+3,3) + Z1*RFVV(1,ILEV64+3,4)))
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
  IF (KHOR == 1) IOFF=0

  DO JLEV=1,KFLEV
    IF (KHOR == 0) IOFF=YDSL%NASLB1*JLEV

    ! * Calculation of linear weights, KL0.
!CDIR NODEP
!DIR$ PREFERVECTOR
    !dir$ ivdep
    DO JROF=KST,KPROF
      ZLO1=PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+1)
      ILO1(JROF)=INT(ZLO1)
      PDLO(JROF,JLEV,1)=ZLO1-ILO1(JROF)
      ZLO2=PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+2)
      ILO2(JROF)=INT(ZLO2)
      PDLO(JROF,JLEV,2)=ZLO2-ILO2(JROF)
    ENDDO

    ! optim: does not vectorize, use of IL (kind JPIA) since nslext is 64-bit indexed
!CDIR NODEP
!DIR$ PREFERVECTOR
    !dir$ ivdep
    DO JROF=KST,KPROF
      IL=ILAV(JROF,JLEV)
      KL0(JROF,JLEV,1)=IADDR(ILAV(JROF,JLEV)+1)+YDSL%NSLEXT(ILO1(JROF),IL+1)+IOFF
      KL0(JROF,JLEV,2)=IADDR(ILAV(JROF,JLEV)+2)+YDSL%NSLEXT(ILO2(JROF),IL+2)+IOFF
    ENDDO

  ENDDO

  ! * Mask calculation for on-demand communications (using KL0)
  IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        YDSL%MASK_SL2(KL0(JROF,JLEV,1):KL0(JROF,JLEV,1)+3)=1
      ENDDO
    ENDDO
    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        YDSL%MASK_SL2(KL0(JROF,JLEV,2):KL0(JROF,JLEV,2)+3)=1
      ENDDO
    ENDDO
  ENDIF
ENDIF

!        2.03  Coordinates and weights for 12 points interpolations.

IF (KWIS == 202 .OR. KWIS == 203) THEN
  IF (KHOR == 1) IOFF=0

  DO JLEV=1,KFLEV
    IF (KHOR == 0) IOFF=YDSL%NASLB1*JLEV

    ! * Calculation of linear weights, KL0.
!CDIR NODEP
!DIR$ PREFERVECTOR
    !dir$ ivdep
    DO JROF=KST,KPROF
      ! zonal interpolation: linear weights for 4 lat. cicles
      ! as the grid is regular in the zonal direction,
      ! the cubic weight computation does not need 
      ! other input than linear weights (LASCAW_CLO)
      ! general case
      ZLO   =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV))
      ILO(JROF)=INT(ZLO)
      PDLO(JROF,JLEV,0)=ZLO-ILO(JROF)
      ZLO1  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+1)
      ILO1(JROF) =INT(ZLO1)
      PDLO(JROF,JLEV,1)=ZLO1-ILO1(JROF)
      ZLO2  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+2)
      ILO2(JROF) =INT(ZLO2)
      PDLO(JROF,JLEV,2)=ZLO2-ILO2(JROF)
      ZLO3  =PLON(JROF,JLEV)*PLSDEPI(ILAV(JROF,JLEV)+3)
      ILO3(JROF)=INT(ZLO3)
      PDLO(JROF,JLEV,3)=ZLO3-ILO3(JROF)
    ENDDO

    ! optim: does not vectorize, use of IL (kind JPIA) since nslext is 64-bit indexed
!CDIR NODEP
!DIR$ PREFERVECTOR
    !dir$ ivdep
    DO JROF=KST,KPROF
      IL=ILAV(JROF,JLEV)
      KL0(JROF,JLEV,0)=IADDR(ILAV(JROF,JLEV))+YDSL%NSLEXT(ILO(JROF),IL)+IOFF
      KL0(JROF,JLEV,1)=IADDR(ILAV(JROF,JLEV)+1)+YDSL%NSLEXT(ILO1(JROF),IL+1)+IOFF
      KL0(JROF,JLEV,2)=IADDR(ILAV(JROF,JLEV)+2)+YDSL%NSLEXT(ILO2(JROF),IL+2)+IOFF
      KL0(JROF,JLEV,3)=IADDR(ILAV(JROF,JLEV)+3)+YDSL%NSLEXT(ILO3(JROF),IL+3)+IOFF
    ENDDO

  ENDDO

  ! * Mask calculation for on-demand communications:
  IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        YDSL%MASK_SL2(KL0(JROF,JLEV,0):KL0(JROF,JLEV,0)+3)=1
      ENDDO
    ENDDO
    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        YDSL%MASK_SL2(KL0(JROF,JLEV,1):KL0(JROF,JLEV,1)+3)=1
      ENDDO
    ENDDO
    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        YDSL%MASK_SL2(KL0(JROF,JLEV,2):KL0(JROF,JLEV,2)+3)=1
      ENDDO
    ENDDO
    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        YDSL%MASK_SL2(KL0(JROF,JLEV,3):KL0(JROF,JLEV,3)+3)=1
      ENDDO
    ENDDO
  ENDIF

  DO JLEV=1,KFLEV
    !dir$ ivdep
    DO JROF=KST,KPROF
      ! meridional interpolation: linear weights and input for cubic weights (LASCAW_CLA)
      ! general case 
      Z1 = PLATI(ILAV(JROF,JLEV)+1)
      Z2 = PLATI(ILAV(JROF,JLEV)+2)
      PDLAT(JROF,JLEV)=(PLAT(JROF,JLEV)-Z1)/(Z2-Z1)  
      ZDA   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV))
      ZDB   =PLAT(JROF,JLEV)-Z1
      ZDC   =PLAT(JROF,JLEV)-Z2
      ZDD   =PLAT(JROF,JLEV)-PLATI(ILAV(JROF,JLEV)+3)
      Z2 = ZDA*ZDB
      ZCLA(JROF,JLEV,1)=(ZDA*ZDC)*ZDD*PIPI(ILAV(JROF,JLEV)+1,1)
      ZCLA(JROF,JLEV,2)=Z2*ZDD*PIPI(ILAV(JROF,JLEV)+1,2)
      ZCLA(JROF,JLEV,3)=Z2*ZDC*PIPI(ILAV(JROF,JLEV)+1,3)

      ! COMAD meridional interpolation 
      IF (LLCOMADH) THEN
        PDLAMAD(JROF,JLEV)=PDLAT(JROF,JLEV)*PSTDDISV(JROF,JLEV)+0.5_JPRB*(1._JPRB-PSTDDISV(JROF,JLEV))
        Z1 = 1._JPRB/PSTDDISV(JROF,JLEV)
        Z4 = 0.5_JPRB*(Z2-Z1)*(Z1-1._JPRB)
        ZDA = (ZDA-ZDB)*Z1
        ZDD = (ZDD-ZDC)*Z1
        ZDB = ZDB +Z4
        ZDC = ZDC -Z4
        ZDA = ZDA + ZDB
        ZDD = ZDD + ZDC
        Z2 = ZDA*ZDB
        Z3 = PSTDDISV(JROF,JLEV)**3
        PCLAMAD(JROF,JLEV,1,1)=ZDA*ZDC*ZDD*PIPI(ILAV(JROF,JLEV)+1,1)*Z3
        PCLAMAD(JROF,JLEV,2,1)=Z2*ZDD*PIPI(ILAV(JROF,JLEV)+1,2)*Z3
        PCLAMAD(JROF,JLEV,3,1)=Z2*ZDC*PIPI(ILAV(JROF,JLEV)+1,3)*Z3

        PDLOMAD(JROF,JLEV,0)=PDLO(JROF,JLEV,0)*PSTDDISU(JROF,JLEV)+0.5_JPRB*(1._JPRB-PSTDDISU(JROF,JLEV))
        PDLOMAD(JROF,JLEV,1)=PDLO(JROF,JLEV,1)*PSTDDISU(JROF,JLEV)+0.5_JPRB*(1._JPRB-PSTDDISU(JROF,JLEV))
        PDLOMAD(JROF,JLEV,2)=PDLO(JROF,JLEV,2)*PSTDDISU(JROF,JLEV)+0.5_JPRB*(1._JPRB-PSTDDISU(JROF,JLEV))
        PDLOMAD(JROF,JLEV,3)=PDLO(JROF,JLEV,3)*PSTDDISU(JROF,JLEV)+0.5_JPRB*(1._JPRB-PSTDDISU(JROF,JLEV))
      ELSE
        PDLAMAD(JROF,JLEV)=PDLAT(JROF,JLEV)
        PCLAMAD(JROF,JLEV,1,1)=ZCLA(JROF,JLEV,1)
        PCLAMAD(JROF,JLEV,2,1)=ZCLA(JROF,JLEV,2)
        PCLAMAD(JROF,JLEV,3,1)=ZCLA(JROF,JLEV,3)
        PDLOMAD(JROF,JLEV,0)= PDLO(JROF,JLEV,0)
        PDLOMAD(JROF,JLEV,1)= PDLO(JROF,JLEV,1)
        PDLOMAD(JROF,JLEV,2)= PDLO(JROF,JLEV,2)
        PDLOMAD(JROF,JLEV,3)= PDLO(JROF,JLEV,3)
      ENDIF
    ENDDO
  ENDDO

  ! * Calculation of PCLA and PCLASLD:
  CALL LASCAW_CLA(YDSL,YDDYNA,KFLEV,KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,&
   & ILAV,ZCLA,PDLAT,PKAPPA,&
   & PSLD,PSLDW,PCLA(:,:,:,1),PCLASLD(:,:,:,1))

  ! * Calculation of PCLO and PCLOSLD:
  CALL LASCAW_CLO(YDDYNA,KFLEV,&
   & KPROMB,KST,KPROF,LLCOMADH,LLT_SLHD,ZSLHDKMINH,&
   & PDLO(:,:,1),PDLOMAD(:,:,1),PKAPPA,&
   & PCLO(:,:,:,1,1),PCLOMAD(:,:,:,1,1),PCLOSLD(:,:,:,1,1))
  CALL LASCAW_CLO(YDDYNA,KFLEV,&
   & KPROMB,KST,KPROF,LLCOMADH,LLT_SLHD,ZSLHDKMINH,&
   & PDLO(:,:,2),PDLOMAD(:,:,2),PKAPPA,&
   & PCLO(:,:,:,2,1),PCLOMAD(:,:,:,2,1),PCLOSLD(:,:,:,2,1))
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('LASCAW',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAW
