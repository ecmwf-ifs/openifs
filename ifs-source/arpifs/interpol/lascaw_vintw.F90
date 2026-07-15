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

SUBROUTINE LASCAW_VINTW(YDDYNA,KPROMA,KFLEV,KST,KEND,&
 & LDCOMAD,LDT_SLHD,LDSLVF,LDSLHDHEAT,PSLHDKMIN,&
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
!          KPROMA    - horizontal dimension.
!          KFLEV    - vertical dimension.
!          KST      - first element of arrays where computations are performed.
!          KEND    - depth of work.
!          LDCOMAD  - key for COMAD.
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
!        H Petithomme (March 2021): optimization (test hoisting, access lowering)
!        R. El Khatib 01-Jun-2022 Remove JPDUP
!        F. Vana  14-Mar-2024: Fixed CY49R1
!        R. El Khatib 27-Mar-2024 Vectorization safety
!        R. El Khatib 27-Mar-2024 Additional cleanings to simplify the phasing to cycle 50
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK 

! arp/ifs dependencies to be solved later.
USE YOMDYNA  , ONLY : TDYNA
IMPLICIT NONE

TYPE(TDYNA),        INTENT(IN)  :: YDDYNA
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLEV
INTEGER(KIND=JPIM), INTENT(IN)  :: KST
INTEGER(KIND=JPIM), INTENT(IN)  :: KEND
LOGICAL           , INTENT(IN)  :: LDCOMAD
LOGICAL           , INTENT(IN)  :: LDT_SLHD(3)
LOGICAL           , INTENT(IN)  :: LDSLVF
LOGICAL           , INTENT(IN)  :: LDSLHDHEAT
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLHDKMIN
INTEGER(KIND=JPIM), INTENT(IN)  :: KLEV(KPROMA,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PLEV(KPROMA,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PDVER(KPROMA,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PDVERMAD(KPROMA,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKAPPA(KPROMA,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSTDDISW(KPROMA,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKAPPAT(KPROMA,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PVETA(0:KFLEV+1)
REAL(KIND=JPRB)   , INTENT(IN)  :: PVCUICO(4,0:KFLEV-1)
REAL(KIND=JPRB)   , INTENT(IN)  :: PVSLD(3,0:KFLEV-1)
REAL(KIND=JPRB)   , INTENT(IN)  :: PVSLDW(3,3,0:KFLEV-1)
REAL(KIND=JPRB)   , INTENT(IN)  :: PVSLVF(3,3,0:KFLEV-1)
REAL(KIND=JPRB)   , INTENT(OUT) :: PVINTW(KPROMA,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PVINTWMAD(KPROMA,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PVINTWSLD(KPROMA,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PVINTWSLT(KPROMA,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PVINTWSLVF(KPROMA,KFLEV,3)

INTEGER(KIND=JPIM) :: JROF,ILEV,JLEV
LOGICAL :: LLSLHD,LLSLHDQUAD
REAL(KIND=JPRB) :: ZWH1(KPROMA),ZWH2(KPROMA),ZWH3(KPROMA),ZWDS1(KPROMA),ZWDS2(KPROMA),ZWDS3(KPROMA)
REAL(KIND=JPRB) :: ZWD1,ZWD2,ZWD3,ZWL1,ZWL2,ZWL3,ZSQR
REAL(KIND=JPRB) :: ZD0,ZD1,ZD2,ZD3,ZDV,Z1,Z2,ZC2,ZC3,ZC4,ZSLHDKMIN,ZKAP,ZKAPT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('LASCAW_VINTW',0,ZHOOK_HANDLE)

LLSLHD=LDT_SLHD(1)
LLSLHDQUAD=LDT_SLHD(2)

! * Calculation of PVINTW and PVINTWSLD/T:
  ! optim: manual hoisting of SLHD and COMAD options helps compiler optimizer
IF (LLSLHD) THEN
    ! branch SLHD, old or quad scheme + varying SLD/SLT coefs
    ! optim: does not vectorize
  DO JLEV=1,KFLEV
    !CDIR NODEP
    DO JROF=KST,KEND
      ILEV=KLEV(JROF,JLEV)
      IF (ILEV >= 1.AND.ILEV <= KFLEV-3) THEN
        Z1=PVETA(ILEV+1)
        Z2=PVETA(ILEV+2)
        ZD0=PLEV(JROF,JLEV)-PVETA(ILEV)
        ZD1=PLEV(JROF,JLEV)-Z1
        ZD2=PLEV(JROF,JLEV)-Z2
        ZD3=PLEV(JROF,JLEV)-PVETA(ILEV+3)

        ZC2=PVCUICO(2,ILEV)
        ZC3=PVCUICO(3,ILEV)
        ZC4=PVCUICO(4,ILEV)

        ZDV=ZD0*ZD1
        ZWH1(JROF)=ZD0*ZD2*ZD3*ZC2
        ZWH2(JROF)=ZDV*ZD3*ZC3
        ZWH3(JROF)=ZDV*ZD2*ZC4

        ! warning: LSLHD and LSLHDQUAD can both be true, so no test simplification
        ! in this branch LLSLHD, alternative means old (linear) SLHD scheme
        ZWL2=PDVER(JROF,JLEV)
        ZWL1=1.0_JPRB-ZWL2

        IF (LLSLHDQUAD) THEN
          ZSQR=ZWL1*ZWL2
          ZWL1=ZSQR*PVSLD(1,ILEV) + ZWL1
          ZWL2=ZSQR*PVSLD(2,ILEV) + ZWL2
          ZWL3=ZSQR*PVSLD(3,ILEV)
        ELSE
          ZWL3=0.0_JPRB
        ENDIF

        ZWD1=ZWH1(JROF)+YDDYNA%SLHDKMAX*(ZWL1-ZWH1(JROF))
        ZWD2=ZWH2(JROF)+YDDYNA%SLHDKMAX*(ZWL2-ZWH2(JROF))
        ZWD3=ZWH3(JROF)+YDDYNA%SLHDKMAX*(ZWL3-ZWH3(JROF))
        ZWDS1(JROF)=PVSLDW(1,1,ILEV)*ZWD1+PVSLDW(1,2,ILEV)*ZWD2+PVSLDW(1,3,ILEV)*ZWD3
        ZWDS2(JROF)=PVSLDW(2,1,ILEV)*ZWD1+PVSLDW(2,2,ILEV)*ZWD2+PVSLDW(2,3,ILEV)*ZWD3
        ZWDS3(JROF)=PVSLDW(3,1,ILEV)*ZWD1+PVSLDW(3,2,ILEV)*ZWD2+PVSLDW(3,3,ILEV)*ZWD3

        ! optim: no use of sign/abs functions anymore
        IF (PKAPPA(JROF,JLEV) < 0._JPRB) THEN
          ZKAP=-PKAPPA(JROF,JLEV)
          ZSLHDKMIN=YDDYNA%SLHDKREF
        ELSE
          ZKAP=PKAPPA(JROF,JLEV)
          ZSLHDKMIN=PSLHDKMIN
        ENDIF

        ZWH1(JROF)=ZWH1(JROF)+ZSLHDKMIN*(ZWL1-ZWH1(JROF))
        ZWH2(JROF)=ZWH2(JROF)+ZSLHDKMIN*(ZWL2-ZWH2(JROF))
        ZWH3(JROF)=ZWH3(JROF)+ZSLHDKMIN*(ZWL3-ZWH3(JROF))
        PVINTW(JROF,JLEV,1)=ZWH1(JROF)
        PVINTW(JROF,JLEV,2)=ZWH2(JROF)
        PVINTW(JROF,JLEV,3)=ZWH3(JROF)
        PVINTWSLD(JROF,JLEV,1)=ZWH1(JROF)+ZKAP*(ZWDS1(JROF)-ZWH1(JROF))
        PVINTWSLD(JROF,JLEV,2)=ZWH2(JROF)+ZKAP*(ZWDS2(JROF)-ZWH2(JROF))
        PVINTWSLD(JROF,JLEV,3)=ZWH3(JROF)+ZKAP*(ZWDS3(JROF)-ZWH3(JROF))

       ! When LSLHD=true it must be always LDCOMAD=false.
       !IF (LDCOMAD) THEN
       !  ZDV = 0.5_JPRB*(Z2-Z1)*(1._JPRB-PSTDDISW(JROF,JLEV))
       !  ZD1 = ZD1*PSTDDISW(JROF,JLEV)+ZDV
       !  ZD2 = ZD2*PSTDDISW(JROF,JLEV)-ZDV
       !  ZD0 = Z1-PVETA(ILEV)+ZD1
       !  ZD3 = Z2-PVETA(ILEV+3)+ZD2
       !ENDIF

        PVINTWMAD(JROF,JLEV,1) = ZD0*ZD2*ZD3*ZC2
        PVINTWMAD(JROF,JLEV,2) = ZD0*ZD1*ZD3*ZC3
        PVINTWMAD(JROF,JLEV,3) = ZD0*ZD1*ZD2*ZC4

      ELSE

        PVINTW(JROF,JLEV,1)=1.0_JPRB-PDVER(JROF,JLEV)
        PVINTW(JROF,JLEV,2)=PDVER(JROF,JLEV)
        PVINTW(JROF,JLEV,3)=0.0_JPRB
        PVINTWSLD(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
        PVINTWSLD(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
        PVINTWSLD(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)

        PVINTWMAD(JROF,JLEV,1) = 1.0_JPRB-PDVERMAD(JROF,JLEV)
        PVINTWMAD(JROF,JLEV,2) = PDVERMAD(JROF,JLEV)
        PVINTWMAD(JROF,JLEV,3) = 0.0_JPRB

      ENDIF

    ENDDO

    IF (LDSLHDHEAT) THEN
      DO JROF=KST,KEND
        ILEV=KLEV(JROF,JLEV)
        IF (ILEV >= 1.AND.ILEV <= KFLEV-3) THEN
          ZKAPT = ABS(PKAPPAT(JROF,JLEV))
          PVINTWSLT(JROF,JLEV,1)=ZWH1(JROF)+ZKAPT*(ZWDS1(JROF)-ZWH1(JROF))
          PVINTWSLT(JROF,JLEV,2)=ZWH2(JROF)+ZKAPT*(ZWDS2(JROF)-ZWH2(JROF))
          PVINTWSLT(JROF,JLEV,3)=ZWH3(JROF)+ZKAPT*(ZWDS3(JROF)-ZWH3(JROF))
        ELSE
          PVINTWSLT(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
          PVINTWSLT(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
          PVINTWSLT(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)
        ENDIF
      ENDDO
    ENDIF
    IF (LDSLVF) THEN
      DO JROF=KST,KEND
        ILEV=KLEV(JROF,JLEV)
        IF (ILEV >= 1.AND.ILEV <= KFLEV-3) THEN
          PVINTWSLVF(JROF,JLEV,1)=PVSLVF(1,1,ILEV)*PVINTW(JROF,JLEV,1)+PVSLVF(1,2,ILEV)*PVINTW(JROF,JLEV,2)+&
           &                      PVSLVF(1,3,ILEV)*PVINTW(JROF,JLEV,3)
          PVINTWSLVF(JROF,JLEV,2)=PVSLVF(2,1,ILEV)*PVINTW(JROF,JLEV,1)+PVSLVF(2,2,ILEV)*PVINTW(JROF,JLEV,2)+&
           &                      PVSLVF(2,3,ILEV)*PVINTW(JROF,JLEV,3)
          PVINTWSLVF(JROF,JLEV,3)=PVSLVF(3,1,ILEV)*PVINTW(JROF,JLEV,1)+PVSLVF(3,2,ILEV)*PVINTW(JROF,JLEV,2)+&
           &                      PVSLVF(3,3,ILEV)*PVINTW(JROF,JLEV,3)
        ELSE
          PVINTWSLVF(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
          PVINTWSLVF(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
          PVINTWSLVF(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)
        ENDIF
      ENDDO
    ENDIF

  ENDDO

ELSE  ! SLHD=.false.

  ! branch no SLHD or LLSLHDQUAD + constant SLD/SLT coefs
  ! optim: does not vectorize
  DO JLEV=1,KFLEV
    !CDIR NODEP
    DO JROF=KST,KEND
      ILEV=KLEV(JROF,JLEV)
      IF (ILEV >= 1.AND.ILEV <= KFLEV-3) THEN
        Z1=PVETA(ILEV+1)
        Z2=PVETA(ILEV+2)
        ZD0=PLEV(JROF,JLEV)-PVETA(ILEV)
        ZD1=PLEV(JROF,JLEV)-Z1
        ZD2=PLEV(JROF,JLEV)-Z2
        ZD3=PLEV(JROF,JLEV)-PVETA(ILEV+3)

        ZC2=PVCUICO(2,ILEV)
        ZC3=PVCUICO(3,ILEV)
        ZC4=PVCUICO(4,ILEV)

        ZDV=ZD0*ZD1
        ZWH1(JROF)=ZD0*ZD2*ZD3*ZC2
        ZWH2(JROF)=ZDV*ZD3*ZC3
        ZWH3(JROF)=ZDV*ZD2*ZC4

        IF (LLSLHDQUAD) THEN
          ZWL2=PDVER(JROF,JLEV)
          ZWL1=1.0_JPRB-ZWL2
          ZSQR=ZWL1*ZWL2
          ZWL1=ZSQR*PVSLD(1,ILEV) + ZWL1
          ZWL2=ZSQR*PVSLD(2,ILEV) + ZWL2
          ZWL3=ZSQR*PVSLD(3,ILEV)

          ZWH1(JROF)=ZWH1(JROF)+PSLHDKMIN*(ZWL1-ZWH1(JROF))
          ZWH2(JROF)=ZWH2(JROF)+PSLHDKMIN*(ZWL2-ZWH2(JROF))
          ZWH3(JROF)=ZWH3(JROF)+PSLHDKMIN*(ZWL3-ZWH3(JROF))
        ENDIF

        PVINTW(JROF,JLEV,1)=ZWH1(JROF)
        PVINTW(JROF,JLEV,2)=ZWH2(JROF)
        PVINTW(JROF,JLEV,3)=ZWH3(JROF)

        IF (LDCOMAD) THEN
          ZDV = 0.5_JPRB*(Z2-Z1)*(1._JPRB-PSTDDISW(JROF,JLEV))
          ZD1 = ZD1*PSTDDISW(JROF,JLEV)+ZDV
          ZD2 = ZD2*PSTDDISW(JROF,JLEV)-ZDV
          ZD0 = PVETA(ILEV+1)-PVETA(ILEV)  + ZD1
          ZD3 = PVETA(ILEV+2)-PVETA(ILEV+3)+ ZD2
          PVINTWMAD(JROF,JLEV,1) = ZD0*ZD2*ZD3*ZC2
          PVINTWMAD(JROF,JLEV,2) = ZD0*ZD1*ZD3*ZC3
          PVINTWMAD(JROF,JLEV,3) = ZD0*ZD1*ZD2*ZC4
        ELSE
          PVINTWMAD (JROF,JLEV,1)=PVINTW (JROF,JLEV,1)
          PVINTWMAD (JROF,JLEV,2)=PVINTW (JROF,JLEV,2)
          PVINTWMAD (JROF,JLEV,3)=PVINTW (JROF,JLEV,3)
        ENDIF

      ELSE

        PVINTW(JROF,JLEV,1)=1.0_JPRB-PDVER(JROF,JLEV)
        PVINTW(JROF,JLEV,2)=PDVER(JROF,JLEV)
        PVINTW(JROF,JLEV,3)=0.0_JPRB
        PVINTWMAD(JROF,JLEV,1) = 1.0_JPRB-PDVERMAD(JROF,JLEV)
        PVINTWMAD(JROF,JLEV,2) = PDVERMAD(JROF,JLEV)
        PVINTWMAD(JROF,JLEV,3) = 0.0_JPRB

      ENDIF
    ENDDO

    IF (LDSLVF) THEN
      DO JROF=KST,KEND
        ILEV=KLEV(JROF,JLEV)
        IF (ILEV >= 1.AND.ILEV <= KFLEV-3) THEN
          PVINTWSLVF(JROF,JLEV,1)=PVSLVF(1,1,ILEV)*ZWH1(JROF)+PVSLVF(1,2,ILEV)*ZWH2(JROF)+&
           &                      PVSLVF(1,3,ILEV)*ZWH3(JROF)
          PVINTWSLVF(JROF,JLEV,2)=PVSLVF(2,1,ILEV)*ZWH1(JROF)+PVSLVF(2,2,ILEV)*ZWH2(JROF)+&
           &                      PVSLVF(2,3,ILEV)*ZWH3(JROF)
          PVINTWSLVF(JROF,JLEV,3)=PVSLVF(3,1,ILEV)*ZWH1(JROF)+PVSLVF(3,2,ILEV)*ZWH2(JROF)+&
           &                      PVSLVF(3,3,ILEV)*ZWH3(JROF)
        ELSE
          PVINTWSLVF(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
          PVINTWSLVF(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
          PVINTWSLVF(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)
        ENDIF
      ENDDO
    ENDIF

    !CDIR NODEP
    DO JROF=KST,KEND
      PVINTWSLD(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
      PVINTWSLD(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
      PVINTWSLD(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)
    ENDDO

    IF (LDSLHDHEAT) THEN
      !CDIR NODEP
      DO JROF=KST,KEND
        PVINTWSLT(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
        PVINTWSLT(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
        PVINTWSLT(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)
      ENDDO
    ENDIF
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('LASCAW_VINTW',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAW_VINTW
