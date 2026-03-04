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

SUBROUTINE LASCAW_CLO(YDDYNA,KFLEV,&
 & KPROM,KST,KPROF,LDT,PSLHDKMIN,PDLO,PDLOMAD,PKAPPA,PKHTURB,PCLO,PCLOMAD,PCLOSLD)

!     ------------------------------------------------------------------

!**** *LASCAW_CLO  -  Weights for semi-LAgrangian interpolator:
!                     Computes PCLO, PCLOMAD and PCLOSLD for one layer
!                     (high-order meridian weights)

!      Can be also called for zonal high-order weights if plane geometry
!       (no need to code a specific ELASCAW_CLA).

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LASCAW_CLO( ... )

!        Explicit arguments :
!        --------------------

!        * INPUT:
!        KFLEV    - Vertical dimension
!        KPROM    - horizontal dimension.
!        KST      - first element of arrays where computations are performed.
!        KPROF    - depth of work.
!        LDT      - keys for SLHD and horizontal turbulence.
!        PSLHDKMIN - either HOISLH or SLHDKMIN
!        PDLO     - distances for horizontal linear interpolations in longitude.
!        PDLOMAD  -  PDLO for COMAD
!        PKAPPA   - kappa function ("coefficient of SLHD").
!        PKHTURB  - horizontal exchange coefficients for 3D turbulence.

!        * OUTPUT:
!        PCLO     - weights for horizontal cubic interpolations in longitude.
!        PCLOMAD  - cf. PCLO, COMAD case.
!        PCLOSLD  - cf. PCLO, SLHD case.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------

!        No external.
!        Called by some (E)LASCAW.. routines.

!     Reference.
!     ----------

!     Author.
!     -------
!        K. YESSAD, after former LASCAW code (JAN 2009).
!        METEO-FRANCE, CNRM/GMAP.

!     Modifications.
!     --------------
!     F. Vana 22-Feb-2011: Horiz. turbulence and diff on phys. tendencies
!     S. Malardel (Nov 2013): COMAD weights for SL interpolations
!     F. Vana 13-feb-2014 SLHD weights for heat variables
!     K. Yessad (March 2017): simplify level numbering.
!     F. Vana    21-Nov-2017: Option LHOISLT
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK 

! arp/ifs dependencies to be solved later.
USE YOMDYNA  , ONLY : TDYNA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDYNA)       , INTENT(IN)  :: YDDYNA
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLEV
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROM
INTEGER(KIND=JPIM), INTENT(IN)  :: KST
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROF
LOGICAL           , INTENT(IN)  :: LDT(4)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLHDKMIN
REAL(KIND=JPRB)   , INTENT(IN)  :: PDLO(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PDLOMAD(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKAPPA(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKHTURB(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLO(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLOMAD(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLOSLD(KPROM,KFLEV,3)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,JLEV
REAL(KIND=JPRB) :: ZWA1,ZWA2,ZWA3,ZWD1,ZWD2,ZWD3,ZWH1,ZWH2,ZWH3
REAL(KIND=JPRB) :: ZWHMAD1,ZWHMAD2,ZWHMAD3
REAL(KIND=JPRB) :: ZWDS1,ZWDS2,ZWDS3,ZWL1,ZWL2,ZWL3
REAL(KIND=JPRB) :: Z1M2EPSH, ZSIGN, ZSLHDKMIN
LOGICAL :: LLSLHD,LLSLHDQUAD,LLSLHD_OLD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

REAL(KIND=JPRB),PARAMETER :: PP6_R=1.0_JPRB/6.0_JPRB
REAL(KIND=JPRB) :: PD
REAL(KIND=JPRB) :: FLAG1, FLAG2, FLAG3
REAL(KIND=JPRB) :: FQUAD1, FQUAD2, FQUAD3

! weights for cubic Lagrange interpolation (regular nodes)
FLAG1(PD)= 0.5_JPRB*(PD+1.0_JPRB)   *(PD-1.0_JPRB)*(PD-2.0_JPRB)
FLAG2(PD)=-0.5_JPRB*(PD+1.0_JPRB)*PD              *(PD-2.0_JPRB)
FLAG3(PD)= PP6_R   *(PD+1.0_JPRB)*PD*(PD-1.0_JPRB)

! weights for quadratic SLHD interpolation (regular nodes)
FQUAD1(PD)=(1.0_JPRB-PD)*(1.0_JPRB+0.25_JPRB*PD)
FQUAD2(PD)=PD*(1.25_JPRB-0.25_JPRB*PD)
FQUAD3(PD)=0.25_JPRB*PD*(PD-1.0_JPRB)

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_CLO',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

LLSLHD=LDT(1)
LLSLHDQUAD=LDT(2)
LLSLHD_OLD=LDT(3)

! * Auxiliary quantity for Laplacian smoother:
Z1M2EPSH=1.0_JPRB-2.0_JPRB*YDDYNA%SLHDEPSH

! * Calculation of PCLO and PCLOSLD:
DO JLEV=1,KFLEV
!CDIR NODEP
DO JROF=KST,KPROF

  ZWH1=FLAG1(PDLO(JROF,JLEV))
  ZWH2=FLAG2(PDLO(JROF,JLEV))
  ZWH3=FLAG3(PDLO(JROF,JLEV))
  ZWHMAD1=FLAG1(PDLOMAD(JROF,JLEV))
  ZWHMAD2=FLAG2(PDLOMAD(JROF,JLEV))
  ZWHMAD3=FLAG3(PDLOMAD(JROF,JLEV))
  IF (LLSLHDQUAD) THEN
    ZWL1=FQUAD1(PDLO(JROF,JLEV))
    ZWL2=FQUAD2(PDLO(JROF,JLEV))
    ZWL3=FQUAD3(PDLO(JROF,JLEV))
  ELSEIF (LLSLHD_OLD) THEN
    ZWL2=PDLO(JROF,JLEV)
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
    ZWDS1=Z1M2EPSH*ZWD1+YDDYNA%SLHDEPSH*ZWD2
    ZWDS2=YDDYNA%SLHDEPSH*ZWD1+Z1M2EPSH*ZWD2
    ZWDS3=YDDYNA%SLHDEPSH*ZWD2+ZWD3
    PCLO(JROF,JLEV,1)=ZWA1
    PCLO(JROF,JLEV,2)=ZWA2
    PCLO(JROF,JLEV,3)=ZWA3
    PCLOSLD(JROF,JLEV,1)=ZWA1+ABS(PKAPPA(JROF,JLEV))*(ZWDS1-ZWA1)
    PCLOSLD(JROF,JLEV,2)=ZWA2+ABS(PKAPPA(JROF,JLEV))*(ZWDS2-ZWA2)
    PCLOSLD(JROF,JLEV,3)=ZWA3+ABS(PKAPPA(JROF,JLEV))*(ZWDS3-ZWA3)
  ELSEIF (LLSLHDQUAD) THEN
    ZWA1=ZWH1+PSLHDKMIN*(ZWL1-ZWH1)
    ZWA2=ZWH2+PSLHDKMIN*(ZWL2-ZWH2)
    ZWA3=ZWH3+PSLHDKMIN*(ZWL3-ZWH3)
    PCLO(JROF,JLEV,1)=ZWA1
    PCLO(JROF,JLEV,2)=ZWA2
    PCLO(JROF,JLEV,3)=ZWA3
    PCLOSLD(JROF,JLEV,1)=ZWA1
    PCLOSLD(JROF,JLEV,2)=ZWA2
    PCLOSLD(JROF,JLEV,3)=ZWA3
  ELSE
    PCLO(JROF,JLEV,1)=ZWH1
    PCLO(JROF,JLEV,2)=ZWH2
    PCLO(JROF,JLEV,3)=ZWH3
    PCLOSLD(JROF,JLEV,1)=ZWH1
    PCLOSLD(JROF,JLEV,2)=ZWH2
    PCLOSLD(JROF,JLEV,3)=ZWH3
  ENDIF
  PCLOMAD(JROF,JLEV,1)=ZWHMAD1
  PCLOMAD(JROF,JLEV,2)=ZWHMAD2
  PCLOMAD(JROF,JLEV,3)=ZWHMAD3
ENDDO
ENDDO

! In case of 3D turbulence apply also the horizontal Laplacian to
!  both PCLO and PCLOSLD:
IF (LDT(4)) THEN
  DO JLEV=1,KFLEV
  !CDIR NODEP
  DO JROF=KST,KPROF
    ZWDS1=(1._JPRB-2._JPRB*PKHTURB(JROF,JLEV))*PCLO(JROF,JLEV,1)&
      &  +PKHTURB(JROF,JLEV)*PCLO(JROF,JLEV,2)
    ZWDS2=(1._JPRB-2._JPRB*PKHTURB(JROF,JLEV))*PCLO(JROF,JLEV,2)&
      &  +PKHTURB(JROF,JLEV)*PCLO(JROF,JLEV,1)
    ZWDS3=PKHTURB(JROF,JLEV)*PCLO(JROF,JLEV,2)+PCLO(JROF,JLEV,3)
    PCLO(JROF,JLEV,1)=ZWDS1
    PCLO(JROF,JLEV,2)=ZWDS2
    PCLO(JROF,JLEV,3)=ZWDS3
    ZWDS1=(1._JPRB-2._JPRB*PKHTURB(JROF,JLEV))*PCLOSLD(JROF,JLEV,1)&
      &  +PKHTURB(JROF,JLEV)*PCLOSLD(JROF,JLEV,2)
    ZWDS2=(1._JPRB-2._JPRB*PKHTURB(JROF,JLEV))*PCLOSLD(JROF,JLEV,2)&
      &  +PKHTURB(JROF,JLEV)*PCLOSLD(JROF,JLEV,1)
    ZWDS3=PKHTURB(JROF,JLEV)*PCLOSLD(JROF,JLEV,2)+PCLOSLD(JROF,JLEV,3)
    PCLOSLD(JROF,JLEV,1)=ZWDS1
    PCLOSLD(JROF,JLEV,2)=ZWDS2
    PCLOSLD(JROF,JLEV,3)=ZWDS3
  ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_CLO',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAW_CLO
