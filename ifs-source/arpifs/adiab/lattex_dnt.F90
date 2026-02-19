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

SUBROUTINE LATTEX_DNT(YDGEOMETRY,YDLDDH,YDRIP,YDDYN,YDDYNA,KST,KPROF,&
 & LDSETTLS,KXLAG,PESGP,PESGM,PXT0,PXT9,PMOY1X,PMIXNL,&
 & PXSI,PXNLT9,PXT1,PXL0,PXL9,PXLF9,PCXNLT9,&
 & PSIDDHXT1,PSIDDHXT9,PSIDDHXL0,PXLF0,LDNESC)

!------------------------------------------------------------------------------
! LATTEX_DNT - Semi-Lagrangian scheme.
!              Computation of the t and t-dt useful quantities
!              at grid-points. Equations for tri-dimensional
!              variables for 2TL scheme.

! Purpose
! -------

! Interface
! ---------
!   CALL LATTEX_DNT(..)

! Explicit arguments :
! --------------------

! * INPUT:
!        KST       - first element of work.
!        KPROF     - depth of work.
!        LDSETTLS  - .T./.F.: Stable/Conventional extrapolations for SL2TL.
!        KXLAG     - type of SL discretisation 
!        PESGP     - (1 + uncentering factor).
!        PESGM     - (1 - uncentering factor).
!        PXT0      - prognostic variable time t (predictor or SI),
!                    preliminary t+dt (corrector)
!        PXT9      - prognostic variable time t (corrector),
!                    not used for predictor or SI
!        PMOY1X    - full nonlinear model at time t [(Delta t/2) "cursive" A]
!        PMIXNL    - extrapolation control variable for mixed NESC/SETTLS scheme

! * INPUT/OUTPUT:
!        PXSI      - semi-implicit linear model at time t
!                    [- (Delta t/2) "cursive" B]
!        PXNLT9    - buffer used during predictor resp. SI time step
!        PXT1      - t+dt term and other final point terms 
!        PXL0      - second SL quantity to be interpolated (linearly) at O if KXLAG>3
!        PXLF0     - third SL quantity to be interpolated (linearly) at O if KXLAG==3
!        PXL9      - SL quantity to be interpolated at O 
!                    (if NSPLTHOI=1 with diffusive interpolation)
!        PXLF9     - SL quantity to be interpolated at O with high order
!                    interpolation (only meaningfull with NSPLTHOI /= 0)
!        PCXNLT9   - buffer used during corrector
!        PSIDDHX.. - buffers for DDH (linear terms).

! * OPTIONAL INPUT:
!        LDNESC    - cf. LNESC in YOMDYNA.

! Externals
! ---------
!   none
!   Called by LATTEX.

! Method
! ------

! Reference
! ---------
!   Arpege documentation about semi-Lagrangian scheme.

! Author
! ------
!      Mar-2002 J. VIVODA (SHMI/CHMI/LACE) - rationalization of LATTEX

! Modifications
! -------------
!   12-Oct-2002 J. Masek   PC bugfix
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   09-Jun-2004 J. Masek   NH cleaning (LPC_NOTR, LDFULLIMP)
!   01-Jul-2004 K. Yessad  Make clearer the tests for PC scheme.
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Aug 2008): simplify XIDT treatment with PC + cleanings
!   K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!   F. Vana  15-Oct-2009: option NSPLTHOI
!   F. Vana  22-Feb-2011: option K[X]LAG=4
!   K. Yessad (Oct 2013): allow NESC without ICI-PC scheme.
!   K. Yessad (July 2014): Move some variables, rename some variables.
!   J. Vivoda (July 2018): mixed NESC/SETTLS scheme.
! End Modifications
!------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYNA      , ONLY : TDYNA
USE YOMCT3       , ONLY : NSTEP
USE YOMDYN       , ONLY : TDYN
USE YOMRIP       , ONLY : TRIP
USE YOMLDDH      , ONLY : TLDDH

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
TYPE(TLDDH)       ,INTENT(IN)    :: YDLDDH
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
LOGICAL           ,INTENT(IN)    :: LDSETTLS
INTEGER(KIND=JPIM),INTENT(IN)    :: KXLAG
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGM
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMOY1X(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMIXNL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXSI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXNLT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXL9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXLF9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCXNLT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSIDDHXT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSIDDHXT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSIDDHXL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXLF0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
LOGICAL,OPTIONAL  ,INTENT(IN)    :: LDNESC

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV,JROF

!     * ZXNLT0 (resp. ZXNLT1) the non linear term at t (resp. t+dt).
REAL(KIND=JPRB) :: ZXNLT1(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZXNLT0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZXIDT0,ZXIDT9,ZXIGP
REAL(KIND=JPRB) :: ZNESC, ZSETTLS, ZMIXNL

LOGICAL :: LLCT,LLNESC

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LATTEX_DNT',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, NSPLTHOI=>YDDYN%NSPLTHOI, &
 & NFOST=>YDRIP%NFOST, &
 & XIDT=>YDDYN%XIDT, &
 & LRSIDDH=>YDLDDH%LRSIDDH)
!     ------------------------------------------------------------------

!*      1. AUXILIARY VARIABLES.
!       -----------------------

IF (PRESENT(LDNESC)) THEN
  LLNESC=LDNESC
ELSE
  LLNESC=YDDYNA%LNESC
ENDIF

LLCT = (YDDYNA%LPC_FULL .AND. NCURRENT_ITER > 0)

ZXIDT0=1.0_JPRB+XIDT
ZXIDT9=1.0_JPRB+XIDT
ZXIGP=1.0_JPRB+XIDT

!     ------------------------------------------------------------------

!*      2. MAIN CALCULATIONS.
!       ---------------------

IF( .NOT.LLCT )THEN

  !############################################
  ! 2.1 Predictor for LPC_FULL,
  !     or case NSITER=0.
  !############################################

  DO JLEV=1,NFLEVG

    ! nonlinear residual time t
    DO JROF=KST,KPROF
      ZXNLT0(JROF)=PMOY1X(JROF,JLEV)+ZXIDT0*PXSI(JROF,JLEV)
    ENDDO

    ! Fill PXL9,PXLF9,PXL0,PXT1.
    IF ((KXLAG == 2).AND.(NSPLTHOI /= 0)) THEN
      IF (NSTEP <= NFOST .OR. LLNESC) THEN
        DO JROF=KST,KPROF
          PXL9(JROF,JLEV)=PXL9(JROF,JLEV)+PXT0(JROF,JLEV)
          PXLF9(JROF,JLEV)=PXLF9(JROF,JLEV)+PESGM*PMOY1X(JROF,JLEV)
          PXT1(JROF,JLEV)=PXT1(JROF,JLEV)+PESGP*ZXNLT0(JROF)
        ENDDO
      ELSEIF (LDSETTLS) THEN
        DO JROF=KST,KPROF
          PXL9(JROF,JLEV)=PXL9(JROF,JLEV)+PXT0(JROF,JLEV)
          ZNESC   = PESGM*PMOY1X(JROF,JLEV)
          ZSETTLS = PESGM*PMOY1X(JROF,JLEV)+(ZXNLT0(JROF)-PXNLT9(JROF,JLEV)) 
          ZMIXNL  = PMIXNL(JROF,JLEV)
          PXLF9(JROF,JLEV)=PXLF9(JROF,JLEV)+ZMIXNL*ZSETTLS+(1.0_JPRB-ZMIXNL)*ZNESC
          PXT1(JROF,JLEV)=PXT1(JROF,JLEV)+PESGP*ZXNLT0(JROF)
        ENDDO
      ELSE
        ! * remaining case: ldsettls=false, lnesc=false.
        DO JROF=KST,KPROF
          PXL9(JROF,JLEV)=PXL9(JROF,JLEV)+PXT0(JROF,JLEV)
          PXLF9(JROF,JLEV)=PXLF9(JROF,JLEV)+PESGM*PMOY1X(JROF,JLEV)&
           & +0.5_JPRB*PESGM*(ZXNLT0(JROF)-PXNLT9(JROF,JLEV))
          PXT1(JROF,JLEV)=PXT1(JROF,JLEV)&
           & +PESGP*(1.5_JPRB*ZXNLT0(JROF)-0.5_JPRB*PXNLT9(JROF,JLEV))
        ENDDO
      ENDIF
    ELSEIF ((KXLAG == 2).AND.(NSPLTHOI == 0)) THEN
      IF (NSTEP <= NFOST .OR. LLNESC) THEN
        DO JROF=KST,KPROF
          PXL9(JROF,JLEV)=PXL9(JROF,JLEV)+PXT0(JROF,JLEV)&
           & +PESGM*PMOY1X(JROF,JLEV)
          PXT1(JROF,JLEV)=PXT1(JROF,JLEV)+PESGP*ZXNLT0(JROF)
        ENDDO
      ELSEIF (LDSETTLS) THEN
        DO JROF=KST,KPROF
          ZNESC   = PXT0(JROF,JLEV)+PESGM*PMOY1X(JROF,JLEV)
          ZSETTLS = PXT0(JROF,JLEV)+PESGM*PMOY1X(JROF,JLEV)+(ZXNLT0(JROF)-PXNLT9(JROF,JLEV))
          ZMIXNL  = PMIXNL(JROF,JLEV)
          PXL9(JROF,JLEV)=PXL9(JROF,JLEV)+ZMIXNL*ZSETTLS+(1.0_JPRB-ZMIXNL)*ZNESC
          PXT1(JROF,JLEV)=PXT1(JROF,JLEV)+PESGP*ZXNLT0(JROF)
        ENDDO
      ELSE
        ! * remaining case: ldsettls=false, lnesc=false.
        DO JROF=KST,KPROF
          PXL9(JROF,JLEV)=PXL9(JROF,JLEV)+PXT0(JROF,JLEV)&
           & +PESGM*PMOY1X(JROF,JLEV)&
           & +0.5_JPRB*PESGM*(ZXNLT0(JROF)-PXNLT9(JROF,JLEV))
          PXT1(JROF,JLEV)=PXT1(JROF,JLEV)&
           & +PESGP*(1.5_JPRB*ZXNLT0(JROF)-0.5_JPRB*PXNLT9(JROF,JLEV))
        ENDDO
      ENDIF
    ELSEIF (KXLAG >= 3) THEN
      IF (NSTEP <= NFOST .OR. LLNESC) THEN
        DO JROF=KST,KPROF
          PXL0(JROF,JLEV)=PXL0(JROF,JLEV)+PESGM*PMOY1X(JROF,JLEV)
          PXT1(JROF,JLEV)=PXT1(JROF,JLEV)+PESGP*ZXNLT0(JROF)
        ENDDO
      ELSEIF (LDSETTLS) THEN
        IF (YDDYNA%LPC_CHEAP) THEN
          DO JROF=KST,KPROF
            ZNESC    = PESGM*PMOY1X(JROF,JLEV)
            ZSETTLS  = ZXNLT0(JROF)-PXNLT9(JROF,JLEV)
            ZMIXNL   = PMIXNL(JROF,JLEV)
            PXL0 (JROF,JLEV)=PXL0(JROF,JLEV)+ZNESC
            PXLF0(JROF,JLEV)=ZMIXNL*ZSETTLS
            PXT1 (JROF,JLEV)=PXT1(JROF,JLEV)+PESGP*ZXNLT0(JROF)
          ENDDO
        ELSE
          DO JROF=KST,KPROF
            ZNESC   = PESGM*PMOY1X(JROF,JLEV)
            ZSETTLS = ZXNLT0(JROF)-PXNLT9(JROF,JLEV)
            ZMIXNL  = PMIXNL(JROF,JLEV)
            PXL0(JROF,JLEV)=PXL0(JROF,JLEV)+ZNESC+ZMIXNL*ZSETTLS
            PXT1(JROF,JLEV)=PXT1(JROF,JLEV)+PESGP*ZXNLT0(JROF)
          ENDDO
        ENDIF
      ELSE
        ! * remaining case: ldsettls=false, lnesc=false.
        DO JROF=KST,KPROF
          PXL0(JROF,JLEV)=PXL0(JROF,JLEV)+PESGM*PMOY1X(JROF,JLEV)&
           & +0.5_JPRB*PESGM*(ZXNLT0(JROF)-PXNLT9(JROF,JLEV))
          PXT1(JROF,JLEV)=PXT1(JROF,JLEV)&
           & +PESGP*(1.5_JPRB*ZXNLT0(JROF)-0.5_JPRB*PXNLT9(JROF,JLEV))
        ENDDO
      ENDIF
      DO JROF=KST,KPROF
        PXL9(JROF,JLEV)=PXL9(JROF,JLEV)+PXT0(JROF,JLEV)
      ENDDO
    ENDIF

    ! save quantities for corrector step
    IF( YDDYNA%LPC_FULL )THEN
      ! save nonlinear model at time t
      PCXNLT9(KST:KPROF,JLEV)=PMOY1X(KST:KPROF,JLEV)
    ENDIF

    IF( .NOT.LLNESC )THEN
      ! save of nonlinear residual at time t
      ! to be used as nonlinear residual at time t-dt next time step
      DO JROF=KST,KPROF
        PXNLT9(JROF,JLEV)=PMOY1X(JROF,JLEV)+ZXIDT9*PXSI(JROF,JLEV)
      ENDDO
    ENDIF

  ENDDO

ELSE

  !############################################
  ! 2.2 Corrector for LPC_FULL
  !############################################

  DO JLEV=1,NFLEVG

    ! nonlinear residual time t+dt (preliminary state from predictor)

    DO JROF=KST,KPROF
      ZXNLT1(JROF)=PMOY1X(JROF,JLEV)+ZXIDT9*PXSI(JROF,JLEV)
    ENDDO

    IF (KXLAG == 2) THEN
      IF (.NOT. YDDYNA%LPC_CHEAP) THEN
        IF (NSPLTHOI /= 0) THEN
          DO JROF=KST,KPROF
            PXL9(JROF,JLEV)=PXL9(JROF,JLEV)+PXT9(JROF,JLEV)
            PXLF9(JROF,JLEV)=PXLF9(JROF,JLEV)&
             & +PESGM*PCXNLT9(JROF,JLEV)
          ENDDO
        ELSE
          DO JROF=KST,KPROF
            PXL9(JROF,JLEV)=PXL9(JROF,JLEV)+PXT9(JROF,JLEV)&
             & +PESGM*PCXNLT9(JROF,JLEV)
          ENDDO
        ENDIF
      ENDIF
      DO JROF=KST,KPROF
        PXT1(JROF,JLEV)=PXT1(JROF,JLEV)+PESGP*ZXNLT1(JROF)
      ENDDO
    ELSEIF (KXLAG >= 3) THEN
      IF (.NOT. YDDYNA%LPC_CHEAP) THEN
        DO JROF=KST,KPROF
          PXL0(JROF,JLEV)=PXL0(JROF,JLEV)+PESGM*PCXNLT9(JROF,JLEV)
          PXL9(JROF,JLEV)=PXL9(JROF,JLEV)+PXT9(JROF,JLEV)
        ENDDO
      ENDIF
      DO JROF=KST,KPROF
        PXT1(JROF,JLEV)=PXT1(JROF,JLEV)+PESGP*ZXNLT1(JROF)
      ENDDO
    ENDIF

  ENDDO

ENDIF

!########################################################
! 2.3 Addition of preliminary quantity for LAGPHY physics
!########################################################

DO JLEV=1,NFLEVG
  IF(XIDT > 0.0_JPRB)THEN
    DO JROF=KST,KPROF
      PXT1(JROF,JLEV)=PXT1(JROF,JLEV)-ZXIGP*PXSI(JROF,JLEV)
      PXSI(JROF,JLEV)=ZXIGP*PXSI(JROF,JLEV)
    ENDDO
  ELSE
    DO JROF=KST,KPROF
      PXT1(JROF,JLEV)=PXT1(JROF,JLEV)-PESGP*PXSI(JROF,JLEV)
      PXSI(JROF,JLEV)=PESGP*PXSI(JROF,JLEV)
    ENDDO
  ENDIF
ENDDO

!########################################################
! 2.4  DDH computations for SI correction 
!########################################################

IF (LRSIDDH) THEN
  IF ( .NOT.LLCT ) THEN
    IF (KXLAG >= 3) THEN
      IF (NSTEP <= NFOST .OR. LLNESC) THEN
        PSIDDHXT1(KST:KPROF,1:NFLEVG)=PSIDDHXT1(KST:KPROF,1:NFLEVG)&
         & +PESGP*ZXIDT0*PXSI(KST:KPROF,1:NFLEVG)
      ELSEIF (LDSETTLS) THEN
        PSIDDHXL0(KST:KPROF,1:NFLEVG)=PSIDDHXL0(KST:KPROF,1:NFLEVG)&
         & +(ZXIDT0*PXSI(KST:KPROF,1:NFLEVG)-PSIDDHXT9(KST:KPROF,1:NFLEVG))
        PSIDDHXT1(KST:KPROF,1:NFLEVG)=PSIDDHXT1(KST:KPROF,1:NFLEVG)&
         &+PESGP*ZXIDT0*PXSI(KST:KPROF,1:NFLEVG)
      ELSE
        ! * remaining case: ldsettls=false, lnesc=false.
        PSIDDHXL0(KST:KPROF,1:NFLEVG)=PSIDDHXL0(KST:KPROF,1:NFLEVG)&
         & +0.5_JPRB*PESGM*(ZXIDT0*PXSI(KST:KPROF,1:NFLEVG)&
         & -PSIDDHXT9(KST:KPROF,1:NFLEVG))
        PSIDDHXT1(KST:KPROF,1:NFLEVG)=PSIDDHXT1(KST:KPROF,1:NFLEVG)&
         & +PESGP*(1.5_JPRB*ZXIDT0*PXSI(KST:KPROF,1:NFLEVG)&
         & -0.5_JPRB*PSIDDHXT9(KST:KPROF,1:NFLEVG))
      ENDIF
    ENDIF
    IF ( .NOT.LLNESC ) THEN
      ! save of semi-implicit linear term  at time t
      ! to be used at time t-dt next time step
      PSIDDHXT9(KST:KPROF,1:NFLEVG)=ZXIDT9*PXSI(KST:KPROF,1:NFLEVG)
    ENDIF
  ELSE
    IF (KXLAG >= 3) THEN
      PSIDDHXT1(KST:KPROF,1:NFLEVG)=PSIDDHXT1(KST:KPROF,1:NFLEVG)&
       & +PESGP*ZXIDT9*PXSI(KST:KPROF,1:NFLEVG)
    ENDIF
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTEX_DNT',1,ZHOOK_HANDLE)
END SUBROUTINE LATTEX_DNT

