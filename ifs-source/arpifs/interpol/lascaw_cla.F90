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

SUBROUTINE LASCAW_CLA(YDSL,YDDYNA,KFLEV,&
 & KPROM,KST,KPROF,LDT,PSLHDKMIN,&
 & KILA,PWH,PWHMAD,PDLAT,PDLAMAD,PKAPPA,PKHTURB,&
 & PSLD,PSLDW,P3DTW,&
 & PCLA,PCLAMAD,PCLASLD)

!     ------------------------------------------------------------------

!**** *LASCAW_CLA  -  Weights for semi-LAgrangian interpolator:
!                     Computes PCLA, PCLAMAD and PCLASLD for one layer
!                     (high-order zonal weights)
!      Spherical geometry only.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LASCAW_CLA( ... )

!        Explicit arguments :
!        --------------------

!        INPUT:
!          YDSL     - SL_STRUCT definition
!          KFLEV    - Vertical dimension
!          KPROM    - horizontal dimension.
!          KST      - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          LDT      - key for SLHD and horizontal turbulence.
!          KILA     - cf. ILA in LASCAW.
!          PSLHDKMIN - either HOISLH or SLHDKMIN
!          PWH      - cf. ZZWH in LASCAW.
!          PWHMAD   - cf. ZZWHMAD in LASCAW .
!          PDLAT    - distance for horizontal linear interpolations in latitude
!          PDLAMAD  - PDLAT for COMAD
!          PKAPPA   - kappa function ("coefficient of SLHD").
!          PKHTURB  - horizontal exchange coefficients for 3D turbulence
!          PSLD     - auxiliary quantity for SLHD interpolation in latitude
!          PSLDW    - weights for SLHD Laplacian smoother in latitude
!          P3DTW    - weights for 3D turbulence Laplacian smoother in latitude

!        OUTPUT:
!          PCLA    - weights for horizontal cubic interpolations in latitude.
!          PCLAMAD - cf. PCLA, COMAD case.
!          PCLASLD - cf. PCLA, SLHD case.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------

!        No external.
!        Called by LASCAW.

!     Reference.
!     ----------

!     Author.
!     -------
!        K. YESSAD, after former LASCAW code (JAN 2009).
!        METEO-FRANCE, CNRM/GMAP.

!     Modifications.
!     --------------
!      F. Vana  22-Feb-2011: horizontal turbulence and phys tendencies diff
!      G. Mozdzynski (May 2012): further cleaning
!      S. Malardel (Nov 2013): COMAD weights for SL interpolations
!      F. Vana 13-feb-2014 SLHD weights for heat variables
!      K. Yessad (March 2017): simplify level numbering.
!      F. Vana    21-Nov-2017: Option LHOISLT
!      F. Vana  18-Jul-2019: SLVF + cleaning
!      F. Vana  30-Oct-2019: More precision to LSLHDQUAD
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK 

! arp/ifs dependencies to be solved later.
USE YOMDYNA  , ONLY : TDYNA
USE EINT_MOD , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(SL_STRUCT),    INTENT(IN)  :: YDSL
TYPE(TDYNA),        INTENT(IN)  :: YDDYNA
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLEV
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROM
INTEGER(KIND=JPIM), INTENT(IN)  :: KST
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROF
LOGICAL           , INTENT(IN)  :: LDT(4)
INTEGER(KIND=JPIM), INTENT(IN)  :: KILA(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLHDKMIN
REAL(KIND=JPRB)   , INTENT(IN)  :: PWH(KPROM,3,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PWHMAD(KPROM,3,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PDLAT(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PDLAMAD(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKAPPA(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKHTURB(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLD(YDSL%NDGSAH:YDSL%NDGENH,3)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLDW(3,3,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   , INTENT(IN)  :: P3DTW(3,3,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLA(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLAMAD(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLASLD(KPROM,KFLEV,3)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,ILA,JLEV
REAL(KIND=JPRB) :: ZWA1,ZWA2,ZWA3,ZWD1,ZWD2,ZWD3,ZWH1,ZWH2,ZWH3
REAL(KIND=JPRB) :: ZWDS1,ZWDS2,ZWDS3,ZWL1,ZWL2,ZWL3,ZSQR
REAL(KIND=JPRB) :: ZREPSH,ZSIGN,ZSLHDKMIN
LOGICAL :: LLSLHD,LLSLHDQUAD,LLSLHD_OLD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_CLA',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

LLSLHD=LDT(1)
LLSLHDQUAD=LDT(2)
LLSLHD_OLD=LDT(3)

! * Calculation of PCLA and PCLASLD:
DO JLEV=1,KFLEV
!CDIR NODEP
DO JROF=KST,KPROF

  ILA=KILA(JROF,JLEV)

  ZWH1=PWH(JROF,1,JLEV)
  ZWH2=PWH(JROF,2,JLEV)
  ZWH3=PWH(JROF,3,JLEV)
  IF (LLSLHDQUAD) THEN
    ZSQR=PDLAT(JROF,JLEV)*(1._JPRB-PDLAT(JROF,JLEV))
    ZWL1=ZSQR*PSLD(ILA+1,1) + (1._JPRB-PDLAT(JROF,JLEV))
    ZWL2=ZSQR*PSLD(ILA+1,2) + PDLAT(JROF,JLEV)
    ZWL3=ZSQR*PSLD(ILA+1,3)
  ELSEIF (LLSLHD_OLD) THEN
    ZWL2=PDLAT(JROF,JLEV)
    ZWL1=1.0_JPRB-ZWL2
    ZWL3=0.0_JPRB
  ENDIF

  IF (LLSLHD) THEN
    ZSIGN=SIGN(0.5_JPRB,PKAPPA(JROF,JLEV))
    ZSLHDKMIN=(0.5_JPRB+ZSIGN)*PSLHDKMIN - (ZSIGN-0.5_JPRB)*YDDYNA%SLHDKREF
    ZWA1=ZWH1+ZSLHDKMIN*(ZWL1-ZWH1)
    ZWA2=ZWH2+ZSLHDKMIN*(ZWL2-ZWH2)
    ZWA3=ZWH3+ZSLHDKMIN*(ZWL3-ZWH3)
    ZWD1=ZWH1+YDDYNA%SLHDKMAX*(ZWL1-ZWH1)
    ZWD2=ZWH2+YDDYNA%SLHDKMAX*(ZWL2-ZWH2)
    ZWD3=ZWH3+YDDYNA%SLHDKMAX*(ZWL3-ZWH3)
    ZWDS1=PSLDW(1,1,ILA+1)*ZWD1+PSLDW(1,2,ILA+1)*ZWD2+&
     & PSLDW(1,3,ILA+1)*ZWD3
    ZWDS2=PSLDW(2,1,ILA+1)*ZWD1+PSLDW(2,2,ILA+1)*ZWD2+&
     & PSLDW(2,3,ILA+1)*ZWD3
    ZWDS3=PSLDW(3,1,ILA+1)*ZWD1+PSLDW(3,2,ILA+1)*ZWD2+&
     & PSLDW(3,3,ILA+1)*ZWD3
    PCLA(JROF,JLEV,1)=ZWA1
    PCLA(JROF,JLEV,2)=ZWA2
    PCLA(JROF,JLEV,3)=ZWA3
    PCLASLD(JROF,JLEV,1)=ZWA1+ABS(PKAPPA(JROF,JLEV))*(ZWDS1-ZWA1)
    PCLASLD(JROF,JLEV,2)=ZWA2+ABS(PKAPPA(JROF,JLEV))*(ZWDS2-ZWA2)
    PCLASLD(JROF,JLEV,3)=ZWA3+ABS(PKAPPA(JROF,JLEV))*(ZWDS3-ZWA3)
  ELSEIF (LLSLHDQUAD) THEN
    ZWA1=ZWH1+PSLHDKMIN*(ZWL1-ZWH1)
    ZWA2=ZWH2+PSLHDKMIN*(ZWL2-ZWH2)
    ZWA3=ZWH3+PSLHDKMIN*(ZWL3-ZWH3)
    PCLA(JROF,JLEV,1)=ZWA1
    PCLA(JROF,JLEV,2)=ZWA2
    PCLA(JROF,JLEV,3)=ZWA3
    PCLASLD(JROF,JLEV,1)=ZWA1
    PCLASLD(JROF,JLEV,2)=ZWA2
    PCLASLD(JROF,JLEV,3)=ZWA3
  ELSE
    PCLA(JROF,JLEV,1)=ZWH1
    PCLA(JROF,JLEV,2)=ZWH2
    PCLA(JROF,JLEV,3)=ZWH3
    PCLASLD(JROF,JLEV,1)=ZWH1
    PCLASLD(JROF,JLEV,2)=ZWH2
    PCLASLD(JROF,JLEV,3)=ZWH3
  ENDIF
  PCLAMAD(JROF,JLEV,1)=PWHMAD(JROF,1,JLEV)
  PCLAMAD(JROF,JLEV,2)=PWHMAD(JROF,2,JLEV)
  PCLAMAD(JROF,JLEV,3)=PWHMAD(JROF,3,JLEV)
ENDDO
ENDDO

! In case of 3D turbulence apply in addition the horizontal Laplacian to
!  both PCLO and PCLOSLD:
IF (LDT(4)) THEN
  ZREPSH=2.0_JPRB

  DO JLEV=1,KFLEV
  !CDIR NODEP
  DO JROF=KST,KPROF
    ILA=KILA(JROF,JLEV)
    ZWDS1=P3DTW(1,1,ILA+1)*PCLA(JROF,JLEV,1)&
     &   +P3DTW(1,2,ILA+1)*PCLA(JROF,JLEV,2)&
     &   +P3DTW(1,3,ILA+1)*PCLA(JROF,JLEV,3)
    ZWDS2=P3DTW(2,1,ILA+1)*PCLA(JROF,JLEV,1)&
     &   +P3DTW(2,2,ILA+1)*PCLA(JROF,JLEV,2)&
     &   +P3DTW(2,3,ILA+1)*PCLA(JROF,JLEV,3)
    ZWDS3=P3DTW(3,1,ILA+1)*PCLA(JROF,JLEV,1)&
     &   +P3DTW(3,2,ILA+1)*PCLA(JROF,JLEV,2)&
     &   +P3DTW(3,3,ILA+1)*PCLA(JROF,JLEV,3)
    PCLA(JROF,JLEV,1)=PCLA(JROF,JLEV,1)&
     &  +PKHTURB(JROF,JLEV)*ZREPSH*(ZWDS1-PCLA(JROF,JLEV,1))
    PCLA(JROF,JLEV,2)=PCLA(JROF,JLEV,2)&
     &  +PKHTURB(JROF,JLEV)*ZREPSH*(ZWDS2-PCLA(JROF,JLEV,2))
    PCLA(JROF,JLEV,3)=PCLA(JROF,JLEV,3)&
     &  +PKHTURB(JROF,JLEV)*ZREPSH*(ZWDS3-PCLA(JROF,JLEV,3))
    ZWDS1=P3DTW(1,1,ILA+1)*PCLASLD(JROF,JLEV,1)&
     &   +P3DTW(1,2,ILA+1)*PCLASLD(JROF,JLEV,2)&
     &   +P3DTW(1,3,ILA+1)*PCLASLD(JROF,JLEV,3)
    ZWDS2=P3DTW(2,1,ILA+1)*PCLASLD(JROF,JLEV,1)&
     &   +P3DTW(2,2,ILA+1)*PCLASLD(JROF,JLEV,2)&
     &   +P3DTW(2,3,ILA+1)*PCLASLD(JROF,JLEV,3)
    ZWDS3=P3DTW(3,1,ILA+1)*PCLASLD(JROF,JLEV,1)&
     &   +P3DTW(3,2,ILA+1)*PCLASLD(JROF,JLEV,2)&
     &   +P3DTW(3,3,ILA+1)*PCLASLD(JROF,JLEV,3)
    PCLASLD(JROF,JLEV,1)=PCLASLD(JROF,JLEV,1)+PKHTURB(JROF,JLEV)*ZREPSH&
     &   *(ZWDS1-PCLASLD(JROF,JLEV,1))
    PCLASLD(JROF,JLEV,2)=PCLASLD(JROF,JLEV,2)+PKHTURB(JROF,JLEV)*ZREPSH&
     &   *(ZWDS2-PCLASLD(JROF,JLEV,2))
    PCLASLD(JROF,JLEV,3)=PCLASLD(JROF,JLEV,3)+PKHTURB(JROF,JLEV)*ZREPSH&
     &   *(ZWDS3-PCLASLD(JROF,JLEV,3))
  ENDDO
  ENDDO
ENDIF
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_CLA',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAW_CLA
