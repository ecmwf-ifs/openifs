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

SUBROUTINE LATTES(YDGEOMETRY,YDGMV,YDRIP,YDML_DYN,KST,KPROF,PDTS2,PBDT,PESGP,PESGM,&
 & KIBL,POROGL,POROGM,&
 & PSDIV0,PSDIV9,PSDVBC,PRES0,PGMVS,PMIXNL,&
 & PGMV,PGMVT1S,PB1,PB2)

!**** *LATTES*   Semi-Lagrangian scheme.
!                Computation of the t and t-dt useful quantities
!                 at grid-points. Equations for bi-dimensional
!                 variables (continuity equation).

!     Purpose.
!     --------
!        * This subroutine computes the equation quantities to be
!          interpolated at each grid-point of the colocation grid.
!          The terms considered here are the explicit terms and
!          the explicit part of the semi-implicit scheme (terms
!          previously computed in LASSIE or LANHSI).
!          Equations considered here are equations for bi-dimensional
!          variables: continuity equation.
!        * Remark 1: when an alternate averaging is used for linear terms
!          in the 2TL SL scheme, the first timestep is treated differently
!          (first order uncentering), no uncentering is applied to the
!          total term ("cursive A") and the term saved in P[X]NLT9
!          is [ (Delta t/2) ("cursive A" - (1 + xidt) beta "cursive B") ]
!          instead of [ (Delta t/2) ("cursive A" - beta "cursive B") ].
!        * Remark 2: for lsettls=true, uncentering is applied to
!          the 'stable' extrapolation if vesl > 0 to avoid instability
!          in the momentum equation.
!        * Remark 3: for PC schemes:
!          - this routine is called for nsiter=0.
!          - this routine is called for the predictor of lpc_full.
!          - this routine is called for the corrector of lpc_full.

!**   Interface.
!     ----------
!        *CALL* *LATTES(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST         - first element of work.
!          KPROF       - depth of work.
!          PDTS2       - 0.5*"pdt", where "pdt" =
!                        time step for the first time-integration step of
!                        a leap-frog scheme or all time-integration steps of
!                        a two-time level scheme; 2*time step for the following
!                        time-integration steps of a leap-frog scheme.
!          PBDT        - zbt if semi-implicit scheme with unreduced
!                        divergence, zbt*(c**2/GM**2) if semi-implicit
!                        scheme with reduced divergence, where zbt
!                        is equal to PDTS2*BETADT (BETADT is in YOMDYN) or zero
!                        according to the value of configuration.
!          PESGP       - (1 + uncentering factor).
!          PESGM       - (1 - uncentering factor).
!          KIBL        - index into YROROG instance in YDGEOMETRY
!          POROGL      - zonal component of "grad(surf orography)"
!          POROGM      - meridian component of "grad(surf orography)"
!          PSDIV0      - SI term at time t for continuity equation (Nu*D).
!          PSDIV9      - SI term at time t-dt for continuity equation (Nu*D).
!          PSDVBC      - vertical integral of divergence computed in "gpcty",
!                        including the "lrubc" and "delta m=1" contributions
!                        of (etadot d prehyd/d eta).
!          PRES0       - hydrostatic pressure at half levels at t.
!          PGMVS       - GMVS variables at t-dt and t.
!          PMIXNL      - extrapolation control variable for mixed NESC/SETTLS scheme.

!        INPUT/OUTPUT:
!          PGMV        - GMV variables at t-dt and t.
!          PGMVT1S     - GMVS variables at t+dt.
!          PB1         - "SLBUF1" buffer for interpolations.
!          PB2         - "SLBUF2" buffer.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none
!           Called by LACDYN.

!     Reference.
!     ----------
!             Arpege documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!      K. YESSAD (METEO FRANCE/CNRM/GMAP) after old parts 3.3 et 3.4 of LACDYN. 
!      Original : AUGUST 1995.

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Aug 2008): simplify XIDT treatment with PC + cleanings
!   K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (Oct 2013): allow NESC without ICI-PC scheme.
!   K. Yessad (July 2014): Move some variables, rename some variables.
!   K. Yessad (Dec 2016): Prune obsolete options.
!   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   J. Vivoda (July 2018): mixed NESC/SETTLS scheme.
! End Modifications
!------------------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE YOMGMV             , ONLY : TGMV
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST             , ONLY : RD
USE YOMCT3             , ONLY : NSTEP
USE YOMSTA             , ONLY : RTSUR
USE YOMRIP             , ONLY : TRIP

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)           ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)               ,INTENT(INOUT) :: YDGMV
TYPE(TRIP)               ,INTENT(IN)    :: YDRIP
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN)    :: YDML_DYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBDT(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSDIV0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSDIV9(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSDVBC(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRES0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMIXNL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1S(YDGEOMETRY%YRDIM%NPROMA,YDGMV%YT1%NDIMS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
!------------------------------------------------------------------------------
REAL(KIND=JPRB) :: ZMOY1SP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZMOYSPNL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZSPNLT0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZSPNLT1(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZSPNLT_FE(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZSPTB(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZSPT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZVWEI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZCMSLP, ZXIGP
REAL(KIND=JPRB) :: ZXIDT0, ZXIDT9
REAL(KIND=JPRB) :: ZNESC, ZSETTLS, ZMIXNL

LOGICAL :: LLCT, LLCTC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------

#include "verdisint.intfb.h"

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LATTES',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, &
 & YDOROG=>YDGEOMETRY%YROROG(KIBL), YDCVER=>YDGEOMETRY%YRCVER,YDDYN=>YDML_DYN%YRDYN,YDDYNA=>YDML_DYN%YRDYNA, & 
 & YDPTRSLB1=>YDML_DYN%YRPTRSLB1, &
 & YDPTRSLB2=>YDML_DYN%YRPTRSLB2)

ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, NVLAG=>YDDYN%NVLAG, &
 & RCMSLP0=>YDDYN%RCMSLP0, XIDT=>YDDYN%XIDT, &
 & NFOST=>YDRIP%NFOST, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT1=>YDGMV%YT1, YT9=>YDGMV%YT9, &
 & MSLB1C9=>YDPTRSLB1%MSLB1C9, &
 & MSLB1C9_NL=>YDPTRSLB1%MSLB1C9_NL, &
 & MSLB1SP9=>YDPTRSLB1%MSLB1SP9, NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & MSLB2SPSI=>YDPTRSLB2%MSLB2SPSI, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2)
!------------------------------------------------------------------------------

!############################################
! 1. AUXILIARITIES
!############################################

LLCT = YDDYNA%LPC_FULL .AND. NCURRENT_ITER > 0 ! corrector for LPC_FULL
LLCTC = YDDYNA%LPC_CHEAP .AND. NCURRENT_ITER > 0

ZXIDT0=1.0_JPRB+XIDT
ZXIDT9=1.0_JPRB+XIDT
ZXIGP=1.0_JPRB+XIDT

!------------------------------------------------------------------------------

!############################################
! 2. NONLINEAR MODEL
!############################################

ZCMSLP=RCMSLP0/(RD*RTSUR)
DO JROF=KST,KPROF
  ZSPTB(JROF)=PSDVBC(JROF,NFLEVG)/PRES0(JROF,NFLEVG)
ENDDO

! * [vertical weight]_l = [delta B]_l
DO JLEV=1,NFLEVG
  ZVWEI(KST:KPROF,JLEV)=YDVAB%VDELB(JLEV)
ENDDO

! linear model at time t

DO JROF=KST,KPROF
  PB2(JROF,MSLB2SPSI)=PBDT(JROF)*PSDIV0(JROF)
  ZSPTB(JROF)=-PDTS2*ZSPTB(JROF)
ENDDO

! nonlinear model at time t

DO JLEV=1,NFLEVG
  DO JROF=KST,KPROF
    ZMOY1SP(JROF,JLEV)=PDTS2*&
     & (PGMV(JROF,JLEV,YT0%MU)*(PGMVS(JROF,YT0%MSPL)+ZCMSLP*POROGL(JROF))&
     & +PGMV(JROF,JLEV,YT0%MV)*(PGMVS(JROF,YT0%MSPM)+ZCMSLP*POROGM(JROF)))
  ENDDO
  DO JROF=KST,KPROF
    ZMOY1SP(JROF,JLEV)=ZMOY1SP(JROF,JLEV)+ZSPTB(JROF)
  ENDDO
ENDDO

!*       *   Continuity equation.

IF (YDDYNA%LTWOTL) THEN

  IF(.NOT.LLCT)THEN

    !############################################
    ! Predictor for LPC_FULL
    ! or case nsiter=0.
    !############################################

    !---------------------------------------------------------------
    ! 3D buffers 
    !---------------------------------------------------------------
    DO JLEV=1,NFLEVG

      ! nonlinear residual time t
      DO JROF=KST,KPROF
        ZSPNLT0(JROF)=ZMOY1SP(JROF,JLEV)+ZXIDT0*PB2(JROF,MSLB2SPSI)
      ENDDO

      ! Fill PB1(.,MSLB1C9),PGMVT1S(.,YT1%MSP),ZSPNLT_FE:
      IF (NVLAG == 2 .OR. NVLAG == 3) THEN

        IF ( NSTEP <= NFOST .OR. YDDYNA%LNESC ) THEN

          DO JROF=KST,KPROF
            PB1(JROF,MSLB1C9+JLEV-NFLSA)=PB1(JROF,MSLB1C9+JLEV-NFLSA)&
             & +PESGM*ZMOY1SP(JROF,JLEV)
          ENDDO
          IF (YDCVER%LVERTFE) THEN
            DO JROF=KST,KPROF
              ZSPNLT_FE(JROF,JLEV)=ZSPNLT0(JROF)&
               & *ZVWEI(JROF,JLEV)*YDVETA%VFE_RDETAH(JLEV)*PESGP
            ENDDO
          ELSE
            DO JROF=KST,KPROF
              PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)&
               & +PESGP*ZVWEI(JROF,JLEV)*ZSPNLT0(JROF)
            ENDDO
          ENDIF
 
        ELSEIF (YDDYNA%LSETTLS) THEN

          IF (YDDYNA%LPC_CHEAP) THEN
            DO JROF=KST,KPROF
              ZNESC  =PESGM*ZMOY1SP(JROF,JLEV)
              ZSETTLS=ZSPNLT0(JROF)-PGMV(JROF,JLEV,YT9%MSPNL)
              ZMIXNL =PMIXNL(JROF,JLEV)
              PB1(JROF,MSLB1C9+JLEV-NFLSA)=PB1(JROF,MSLB1C9+JLEV-NFLSA)+ZNESC
              PB1(JROF,MSLB1C9_NL+JLEV-NFLSA)=ZMIXNL*ZSETTLS
            ENDDO
          ELSE
            DO JROF=KST,KPROF
              ZNESC  =PESGM*ZMOY1SP(JROF,JLEV)
              ZSETTLS=ZSPNLT0(JROF)-PGMV(JROF,JLEV,YT9%MSPNL)
              ZMIXNL =PMIXNL(JROF,JLEV)
              PB1(JROF,MSLB1C9+JLEV-NFLSA)=PB1(JROF,MSLB1C9+JLEV-NFLSA)&
               & +ZNESC+ZMIXNL*ZSETTLS
            ENDDO
          ENDIF
          IF( YDCVER%LVERTFE )THEN
            DO JROF=KST,KPROF
              ZSPNLT_FE(JROF,JLEV)=&
               & ZSPNLT0(JROF)*PESGP*ZVWEI(JROF,JLEV)*YDVETA%VFE_RDETAH(JLEV)
            ENDDO
          ELSE
            DO JROF=KST,KPROF
              PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)&
               & +PESGP*ZVWEI(JROF,JLEV)*ZSPNLT0(JROF)
            ENDDO
          ENDIF

        ELSE

          ! remaining case: ldsettls=false, lnesc=false.
          DO JROF=KST,KPROF
            PB1(JROF,MSLB1C9+JLEV-NFLSA)=PB1(JROF,MSLB1C9+JLEV-NFLSA)&
             & +PESGM*ZMOY1SP(JROF,JLEV)&
             & +0.5_JPRB*PESGM*(ZSPNLT0(JROF)-PGMV(JROF,JLEV,YT9%MSPNL))
          ENDDO 
          IF (YDCVER%LVERTFE) THEN
            DO JROF=KST,KPROF
              ZSPNLT_FE(JROF,JLEV)=&
               & (1.5_JPRB*ZSPNLT0(JROF)-0.5_JPRB*PGMV(JROF,JLEV,YT9%MSPNL))&
               & *ZVWEI(JROF,JLEV)*YDVETA%VFE_RDETAH(JLEV)*PESGP
            ENDDO
          ELSE
            DO JROF=KST,KPROF
              PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)&
               & +PESGP*ZVWEI(JROF,JLEV)&
               & *(1.5_JPRB*ZSPNLT0(JROF)-0.5_JPRB*PGMV(JROF,JLEV,YT9%MSPNL))
            ENDDO
          ENDIF

        ENDIF

      ENDIF

      ! Save quantities for corrector step
      IF( YDDYNA%LPC_FULL )THEN
        ! save nonlinear model at time t
        PGMV(KST:KPROF,JLEV,YT9%MCSPNL)=ZMOY1SP(KST:KPROF,JLEV)
      ENDIF

      IF( .NOT.YDDYNA%LNESC )THEN
        ! save of nonlinear residual at time t
        ! to be used as nonlinear residual at time t-dt next time step
        DO JROF=KST,KPROF
          PGMV(JROF,JLEV,YT9%MSPNL)=ZMOY1SP(JROF,JLEV)+ZXIDT9*PB2(JROF,MSLB2SPSI)
        ENDDO
      ENDIF

    ENDDO

    !------------------------------------------------------------------------
    ! 3D additional actions
    !------------------------------------------------------------------------
    IF (YDCVER%LVERTFE) THEN
      IF(NVLAG == 2 .OR. NVLAG == 3)THEN
        ZSPNLT_FE(KST:KPROF,0)=0.0_JPRB
        ZSPNLT_FE(KST:KPROF,NFLEVG+1)=0.0_JPRB
        CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',NPROMA,KST,KPROF,NFLEVG,ZSPNLT_FE,ZSPT1)

        DO JROF=KST,KPROF
          PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)+ZSPT1(JROF,NFLEVG+1)
        ENDDO
      ENDIF
    ENDIF

    !------------------------------------------------------------------------
    ! 2D buffers
    !------------------------------------------------------------------------
    IF(NVLAG == 2 .OR. NVLAG == 3)THEN
      DO JROF=KST,KPROF
        PB1(JROF,MSLB1SP9)=PB1(JROF,MSLB1SP9)+PGMVS(JROF,YT0%MSP)&
         & +ZCMSLP*YDOROG%OROG(JROF)
        PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)-ZCMSLP*YDOROG%OROG(JROF)
      ENDDO
    ENDIF

  ELSE

    !############################################
    ! Corrector for LPC_FULL
    !############################################
    !---------------------------------------------------------------
    ! 3D buffers 
    !---------------------------------------------------------------
    DO JLEV=1,NFLEVG

      ! nonlinear residual time t

      DO JROF=KST,KPROF
        ZSPNLT1(JROF)=ZMOY1SP(JROF,JLEV)+ZXIDT9*PB2(JROF,MSLB2SPSI)
      ENDDO

      IF (.NOT.LLCTC) THEN
        DO JROF=KST,KPROF
          PB1(JROF,MSLB1C9+JLEV-NFLSA)=PB1(JROF,MSLB1C9+JLEV-NFLSA)&
           & +PESGM*PGMV(JROF,JLEV,YT9%MCSPNL)
        ENDDO
      ENDIF
      IF (YDCVER%LVERTFE) THEN
        DO JROF=KST,KPROF
          ZSPNLT_FE(JROF,JLEV)=ZSPNLT1(JROF)&
           & *ZVWEI(JROF,JLEV)*YDVETA%VFE_RDETAH(JLEV)*PESGP
        ENDDO
      ELSE
        DO JROF=KST,KPROF
          PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)&
           & +ZVWEI(JROF,JLEV)*ZSPNLT1(JROF)*PESGP
        ENDDO
      ENDIF

    ENDDO

    !------------------------------------------------------------------------
    ! 3D additional actions
    !------------------------------------------------------------------------
    IF (YDCVER%LVERTFE) THEN
      ZSPNLT_FE(KST:KPROF,0)=0.0_JPRB
      ZSPNLT_FE(KST:KPROF,NFLEVG+1)=0.0_JPRB
      CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',NPROMA,KST,KPROF,NFLEVG,ZSPNLT_FE,ZSPT1)

      DO JROF=KST,KPROF
        PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)+ZSPT1(JROF,NFLEVG+1)
      ENDDO
    ENDIF

    !------------------------------------------------------------------------
    ! 2D buffers
    !------------------------------------------------------------------------
    IF (.NOT.LLCTC) THEN
      DO JROF=KST,KPROF
        PB1(JROF,MSLB1SP9)=PB1(JROF,MSLB1SP9)+PGMVS(JROF,YT9%MSP)&
         & +ZCMSLP*YDOROG%OROG(JROF)
      ENDDO
    ENDIF
    DO JROF=KST,KPROF
      PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)-ZCMSLP*YDOROG%OROG(JROF)
    ENDDO

  ENDIF

  ! add preliminary quantity for LAPGHY physics

  IF( XIDT > 0.0_JPRB )THEN
    DO JROF=KST,KPROF
      PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)-ZXIGP*PB2(JROF,MSLB2SPSI)
      PB2(JROF,MSLB2SPSI)=ZXIGP*PB2(JROF,MSLB2SPSI)
    ENDDO
  ELSE
    DO JROF=KST,KPROF
      PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)-PESGP*PB2(JROF,MSLB2SPSI)
      PB2(JROF,MSLB2SPSI)=PESGP*PB2(JROF,MSLB2SPSI)
    ENDDO
  ENDIF

ELSE

  IF(NVLAG == 2 .OR. NVLAG == 3)THEN

    IF(YDCVER%LVERTFE) THEN
      DO JLEV=1,NFLEVG
        DO JROF=KST,KPROF
          ZMOYSPNL(JROF)=ZMOY1SP(JROF,JLEV)+PB2(JROF,MSLB2SPSI)
          PB1(JROF,MSLB1C9+JLEV-NFLSA)=PB1(JROF,MSLB1C9+JLEV-NFLSA)&
           & +PESGM*ZMOYSPNL(JROF)-PESGM*PBDT(JROF)*PSDIV9(JROF)
          ZSPNLT_FE(JROF,JLEV)=PESGP*ZVWEI(JROF,JLEV)*ZMOYSPNL(JROF)
        ENDDO
      ENDDO
      ZSPNLT_FE(KST:KPROF,0)=0.0_JPRB
      ZSPNLT_FE(KST:KPROF,NFLEVG+1)=0.0_JPRB
      CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',NPROMA,KST,KPROF,NFLEVG,ZSPNLT_FE,ZSPT1)

      DO JROF=KST,KPROF
        PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)+ZSPT1(JROF,NFLEVG+1)
      ENDDO
    ELSE
      DO JLEV=1,NFLEVG
        DO JROF=KST,KPROF
          ZMOYSPNL(JROF)=ZMOY1SP(JROF,JLEV)+PB2(JROF,MSLB2SPSI)
          PB1(JROF,MSLB1C9+JLEV-NFLSA)=PB1(JROF,MSLB1C9+JLEV-NFLSA)&
           & +PESGM*ZMOYSPNL(JROF)-PESGM*PBDT(JROF)*PSDIV9(JROF)
          PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)&
           & +PESGP*ZVWEI(JROF,JLEV)*ZMOYSPNL(JROF)
        ENDDO
      ENDDO
    ENDIF

    DO JROF=KST,KPROF
      PB1(JROF,MSLB1SP9)=PB1(JROF,MSLB1SP9)+PGMVS(JROF,YT9%MSP)&
       & +ZCMSLP*YDOROG%OROG(JROF)
      PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)-ZCMSLP*YDOROG%OROG(JROF)&
       & -PESGP*PB2(JROF,MSLB2SPSI)
      PB2(JROF,MSLB2SPSI)=PESGP*PB2(JROF,MSLB2SPSI)
    ENDDO

  ENDIF

ENDIF

!------------------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTES',1,ZHOOK_HANDLE)
END SUBROUTINE LATTES

