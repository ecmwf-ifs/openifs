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

SUBROUTINE LASCAW_VINTW(YDDYNA,KPROM,KFLEV,KST,KPROF,&
 & LDCOMADV,LDT_SLHD,LDSLVF,LDSLHDHEAT,PSLHDKMIN,&
 & KLEV,PLEV,PDVER,PDVERMAD,PSTDDISW,PKAPPA,PKAPPAT,PVETA,PVCUICO,PVSLD,PVSLDW,PVSLVF,&
 & PVINTW,PVINTWMAD,PVINTWSLD,PVINTWSLT,PVINTWSLVF)

!     ------------------------------------------------------------------

!**** *LASCAW_VINTW  -  Weights for semi-LAgrangian interpolator:
!                       Computes PVINTW, PVINTWMAD and PVINTWSLD/T for one layer
!                       (high-order vertical weights)

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LASCAW_VINTW( ... )

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROM    - horizontal dimension.
!          KFLEV    - vertical dimension.
!          KST      - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          LDCOMADV - key for COMAD.
!          LDT_SLHD - key for SLHD.
!          LDSLVF   - key for SLVF.
!          PSLHDKMIN - either HOISLV or SLHDKMIN
!          KLEV     - lower level of the vertical interpolation
!                     grid needed for vertical interpolations.
!          PLEV     - vertical coordinate of the interpolation point.
!          PDVER    - distance for vertical linear interpolation.
!          PDVERMAD - PDVER for COMAD
!          PSTDDISW - STDDISW correction coef. for vertical COMAD
!          PKAPPA   - kappa function ("coefficient of SLHD").
!          PKAPPAT  - kappa function ("coefficient of SLHD") on T.
!          PVETA    - Values of ETA.
!          PVCUICO  - Denominators of the vertical cubic interpolation coef.
!          PVSLD    - auxiliary quantities for vertical SLHD interpolation.
!          PVSLDW   - weights for SLHD vertical Laplacian smoother.
!          PVSLVF   - weights for SLVF vertical Laplacian smoother.

!        OUTPUT:
!          PVINTW    - vertical cubic interpolation weights.
!          PVINTWMAD - vertical cubic interpolation weights, COMAD case.
!          PVINTWSLD - vertical cubic interpolation weights, SLHD case.
!          PVINTWSLT - vertical cubic interpolation weights, SLHD case on T.
!          PVINTWSLVF - vertical cubic interpolation weights, SLVF case on T.

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
!        S. Malardel (Nov 2013): COMAD weights for SL interpolations
!        F. Vana 13-feb-2014 SLHD weights for heat variables
!        K. Yessad (March 2017): simplify level numbering.
!        F. Vana    21-Nov-2017: Option LHOISLT
!        F. Vana  18-Jul-2019: SLVF + cleaning
!        F. Vana  30-Oct-2019: More precision to LSLHDQUAD
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK 

! arp/ifs dependencies to be solved later.
USE YOMDYNA  , ONLY : TDYNA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDYNA),        INTENT(IN)  :: YDDYNA
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROM
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLEV
INTEGER(KIND=JPIM), INTENT(IN)  :: KST
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROF
LOGICAL           , INTENT(IN)  :: LDCOMADV
LOGICAL           , INTENT(IN)  :: LDT_SLHD(3)
LOGICAL           , INTENT(IN)  :: LDSLVF
LOGICAL           , INTENT(IN)  :: LDSLHDHEAT
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLHDKMIN
INTEGER(KIND=JPIM), INTENT(IN)  :: KLEV(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PLEV(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PDVER(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PDVERMAD(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKAPPA(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSTDDISW(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKAPPAT(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PVETA(0:KFLEV+1)
REAL(KIND=JPRB)   , INTENT(IN)  :: PVCUICO(4,0:KFLEV-1)
REAL(KIND=JPRB)   , INTENT(IN)  :: PVSLD(3,0:KFLEV-1)
REAL(KIND=JPRB)   , INTENT(IN)  :: PVSLDW(3,3,0:KFLEV-1)
REAL(KIND=JPRB)   , INTENT(IN)  :: PVSLVF(3,3,0:KFLEV-1)
REAL(KIND=JPRB)   , INTENT(OUT) :: PVINTW(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PVINTWMAD(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PVINTWSLD(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PVINTWSLT(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PVINTWSLVF(KPROM,KFLEV,3)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,ILEVV,JLEV
REAL(KIND=JPRB) :: ZWA1,ZWA2,ZWA3,ZWD1,ZWD2,ZWD3,ZWH1,ZWH2,ZWH3
REAL(KIND=JPRB) :: ZWHMAD1,ZWHMAD2,ZWHMAD3
REAL(KIND=JPRB) :: ZWDS1,ZWDS2,ZWDS3,ZWL1,ZWL2,ZWL3,ZSQR
REAL(KIND=JPRB) :: ZD0,ZD1,ZD2,ZD3,ZSIGN,ZSLHDKMIN
LOGICAL :: LLSLHD,LLSLHDQUAD,LLSLHD_OLD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_VINTW',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

LLSLHD=LDT_SLHD(1)
LLSLHDQUAD=LDT_SLHD(2)
LLSLHD_OLD=LDT_SLHD(3)

! * Calculation of PVINTW and PVINTWSLD/T:
DO JLEV=1,KFLEV
!CDIR NODEP
DO JROF=KST,KPROF

  ILEVV=KLEV(JROF,JLEV)
  IF (ILEVV >= 1.AND.ILEVV <= KFLEV-3) THEN
    ZD0=PLEV(JROF,JLEV)-PVETA(ILEVV  )
    ZD1=PLEV(JROF,JLEV)-PVETA(ILEVV+1)
    ZD2=PLEV(JROF,JLEV)-PVETA(ILEVV+2)
    ZD3=PLEV(JROF,JLEV)-PVETA(ILEVV+3)
    ZWH1=ZD0    *ZD2*ZD3*PVCUICO(2,ILEVV)
    ZWH2=ZD0*ZD1    *ZD3*PVCUICO(3,ILEVV)
    ZWH3=ZD0*ZD1*ZD2    *PVCUICO(4,ILEVV)
    IF (LDCOMADV) THEN
      ZD0 = (ZD0-ZD1)
      ZD3 = (ZD3-ZD2)
      ZD1 = ZD1*PSTDDISW(JROF,JLEV)&
       & +0.5_JPRB*(PVETA(ILEVV+2)-PVETA(ILEVV+1))*(1._JPRB-PSTDDISW(JROF,JLEV))
      ZD2 = ZD2*PSTDDISW(JROF,JLEV)&
       & -0.5_JPRB*(PVETA(ILEVV+2)-PVETA(ILEVV+1))*(1._JPRB-PSTDDISW(JROF,JLEV))
      ZD0 = ZD0 + ZD1
      ZD3 = ZD3 + ZD2
      ZWHMAD1=ZD0     *ZD2*ZD3*PVCUICO(2,ILEVV)
      ZWHMAD2=ZD0*ZD1     *ZD3*PVCUICO(3,ILEVV)
      ZWHMAD3=ZD0*ZD1*ZD2     *PVCUICO(4,ILEVV)
    ELSE
      ZWHMAD1=ZWH1
      ZWHMAD2=ZWH2
      ZWHMAD3=ZWH3
    ENDIF

    PVINTWMAD(JROF,JLEV,1)=ZWHMAD1
    PVINTWMAD(JROF,JLEV,2)=ZWHMAD2
    PVINTWMAD(JROF,JLEV,3)=ZWHMAD3

    IF (LLSLHDQUAD) THEN
      ZSQR=PDVER(JROF,JLEV)*(1._JPRB-PDVER(JROF,JLEV))
      ZWL1=ZSQR*PVSLD(1,ILEVV) + (1._JPRB-PDVER(JROF,JLEV))
      ZWL2=ZSQR*PVSLD(2,ILEVV) + PDVER(JROF,JLEV)
      ZWL3=ZSQR*PVSLD(3,ILEVV)
    ELSEIF (LLSLHD_OLD) THEN
      ZWL2=PDVER(JROF,JLEV)
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
      ZWDS1=PVSLDW(1,1,ILEVV)*ZWD1+PVSLDW(1,2,ILEVV)*ZWD2+&
       &    PVSLDW(1,3,ILEVV)*ZWD3
      ZWDS2=PVSLDW(2,1,ILEVV)*ZWD1+PVSLDW(2,2,ILEVV)*ZWD2+&
       &    PVSLDW(2,3,ILEVV)*ZWD3
      ZWDS3=PVSLDW(3,1,ILEVV)*ZWD1+PVSLDW(3,2,ILEVV)*ZWD2+&
       &    PVSLDW(3,3,ILEVV)*ZWD3
      PVINTW(JROF,JLEV,1)=ZWA1
      PVINTW(JROF,JLEV,2)=ZWA2
      PVINTW(JROF,JLEV,3)=ZWA3
      PVINTWSLD(JROF,JLEV,1)=ZWA1+ABS(PKAPPA(JROF,JLEV))*(ZWDS1-ZWA1)
      PVINTWSLD(JROF,JLEV,2)=ZWA2+ABS(PKAPPA(JROF,JLEV))*(ZWDS2-ZWA2)
      PVINTWSLD(JROF,JLEV,3)=ZWA3+ABS(PKAPPA(JROF,JLEV))*(ZWDS3-ZWA3)
      IF (LDSLHDHEAT) THEN
        PVINTWSLT(JROF,JLEV,1)=ZWA1+ABS(PKAPPAT(JROF,JLEV))*(ZWDS1-ZWA1)
        PVINTWSLT(JROF,JLEV,2)=ZWA2+ABS(PKAPPAT(JROF,JLEV))*(ZWDS2-ZWA2)
        PVINTWSLT(JROF,JLEV,3)=ZWA3+ABS(PKAPPAT(JROF,JLEV))*(ZWDS3-ZWA3)
      ENDIF
      IF (LDSLVF) THEN
        PVINTWSLVF(JROF,JLEV,1)=PVSLVF(1,1,ILEVV)*ZWA1+PVSLVF(1,2,ILEVV)*ZWA2+&
         &                      PVSLVF(1,3,ILEVV)*ZWA3
        PVINTWSLVF(JROF,JLEV,2)=PVSLVF(2,1,ILEVV)*ZWA1+PVSLVF(2,2,ILEVV)*ZWA2+&
         &                      PVSLVF(2,3,ILEVV)*ZWA3
        PVINTWSLVF(JROF,JLEV,3)=PVSLVF(3,1,ILEVV)*ZWA1+PVSLVF(3,2,ILEVV)*ZWA2+&
         &                      PVSLVF(3,3,ILEVV)*ZWA3
      ENDIF
    ELSEIF (LLSLHDQUAD) THEN
      ZWA1=ZWH1+PSLHDKMIN*(ZWL1-ZWH1)
      ZWA2=ZWH2+PSLHDKMIN*(ZWL2-ZWH2)
      ZWA3=ZWH3+PSLHDKMIN*(ZWL3-ZWH3)
      PVINTW(JROF,JLEV,1)=ZWA1
      PVINTW(JROF,JLEV,2)=ZWA2
      PVINTW(JROF,JLEV,3)=ZWA3
      PVINTWSLD(JROF,JLEV,1)=ZWA1
      PVINTWSLD(JROF,JLEV,2)=ZWA2
      PVINTWSLD(JROF,JLEV,3)=ZWA3
     !IF (LDSLHDHEAT) THEN
     !  PVINTWSLT(JROF,JLEV,1)=ZWA1
     !  PVINTWSLT(JROF,JLEV,2)=ZWA2
     !  PVINTWSLT(JROF,JLEV,3)=ZWA3
     !ENDIF
      IF (LDSLVF) THEN
        PVINTWSLVF(JROF,JLEV,1)=PVSLVF(1,1,ILEVV)*ZWA1+PVSLVF(1,2,ILEVV)*ZWA2+&
         &                      PVSLVF(1,3,ILEVV)*ZWA3
        PVINTWSLVF(JROF,JLEV,2)=PVSLVF(2,1,ILEVV)*ZWA1+PVSLVF(2,2,ILEVV)*ZWA2+&
         &                      PVSLVF(2,3,ILEVV)*ZWA3
        PVINTWSLVF(JROF,JLEV,3)=PVSLVF(3,1,ILEVV)*ZWA1+PVSLVF(3,2,ILEVV)*ZWA2+&
         &                      PVSLVF(3,3,ILEVV)*ZWA3
      ENDIF
    ELSE
      PVINTW(JROF,JLEV,1)=ZWH1
      PVINTW(JROF,JLEV,2)=ZWH2
      PVINTW(JROF,JLEV,3)=ZWH3
      PVINTWSLD(JROF,JLEV,1)=ZWH1
      PVINTWSLD(JROF,JLEV,2)=ZWH2
      PVINTWSLD(JROF,JLEV,3)=ZWH3
     !IF (LDSLHDHEAT) THEN
     !  PVINTWSLT(JROF,JLEV,1)=ZWH1
     !  PVINTWSLT(JROF,JLEV,2)=ZWH2
     !  PVINTWSLT(JROF,JLEV,3)=ZWH3
     !ENDIF
      IF (LDSLVF) THEN
        PVINTWSLVF(JROF,JLEV,1)=PVSLVF(1,1,ILEVV)*ZWH1+PVSLVF(1,2,ILEVV)*ZWH2+&
         &                      PVSLVF(1,3,ILEVV)*ZWH3
        PVINTWSLVF(JROF,JLEV,2)=PVSLVF(2,1,ILEVV)*ZWH1+PVSLVF(2,2,ILEVV)*ZWH2+&
         &                      PVSLVF(2,3,ILEVV)*ZWH3
        PVINTWSLVF(JROF,JLEV,3)=PVSLVF(3,1,ILEVV)*ZWH1+PVSLVF(3,2,ILEVV)*ZWH2+&
         &                      PVSLVF(3,3,ILEVV)*ZWH3
      ENDIF
    ENDIF
  ELSE
    PVINTW(JROF,JLEV,1)=1.0_JPRB-PDVER(JROF,JLEV)
    PVINTW(JROF,JLEV,2)=PDVER(JROF,JLEV)
    PVINTW(JROF,JLEV,3)=0.0_JPRB
    IF (LLSLHD) THEN
      PVINTWSLD(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
      PVINTWSLD(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
      PVINTWSLD(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)
      IF (LDSLHDHEAT) THEN
        PVINTWSLT(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
        PVINTWSLT(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
        PVINTWSLT(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)
      ENDIF
    ENDIF
    IF (LDSLVF) THEN
      PVINTWSLVF(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
      PVINTWSLVF(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
      PVINTWSLVF(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)
    ENDIF
    PVINTWMAD(JROF,JLEV,1)=1.0_JPRB-PDVERMAD(JROF,JLEV)
    PVINTWMAD(JROF,JLEV,2)=PDVERMAD(JROF,JLEV)
    PVINTWMAD(JROF,JLEV,3)=0.0_JPRB
  ENDIF
ENDDO
ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_VINTW',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAW_VINTW
