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
SUBROUTINE CALL_SL(YDGEO,YDGMV,&
 ! --- INPUT ------------------------------------------------------------------
 & YDMODEL,PGFL,&
 ! --- INPUT/OUTPUT -----------------------------------------------------------
 & YDSL,KVSEPC,KVSEPL,PB1,PB2,PWRL9,PGMVT1,PGMVT1S,PGFLT1,PGFLPC,&
 ! --- OUTPUT -----------------------------------------------------------------
 & KL0,PLSCAW,PUP9,PVP9,PTP9,PGFLP9,PSPPTCLEARP9,KLH0,PRSCAW,PSCO,PCCO,PLSCAWH,PRSCAWH)

!**** *CALL_SL* - Semi-Lagrangian Interpolation

!     Purpose.
!     --------
!           Semi-Lagrangian Interpolation

!**   Interface.
!     ----------
!        *CALL* *CALL_SL*

!        Explicit arguments :
!        --------------------
!        INPUT:
!          PGFL        - grid-point fields

!        INPUT/OUTPUT:
!          YDSL        - SL_STRUCT definition
!          KVSEPC      - vertical separation (used in SL adjoint, cubic interp)
!          KVSEPL      - vertical separation (used in SL adjoint, linear intrp)
!          PB1         - SLBUF1-buffer for interpolations.
!          PB2         - SLBUF2-buf to communicate info from non lag to lag dyn.
!          PGMVT1      - t+dt GMV (upper air) variables.
!          PGMVT1S     - t+dt GMVS (surface) variables.
!          PGFLT1      - SL-buf to communicate info on unified_treatment
!                        grid-point fields
!          PGFLPC      - to store PC scheme quantities for GFL

!        OUTPUT:
!          KL0         - indices of the four western points of the 16 point
!                        interpolation grid
!          PLSCAW     -  linear weights (distances) for interpolations.
!          PUP9        - Interpolated quantity at O in U-wind eqn
!                        for split physics (if LSLPHY=.T. only).
!          PVP9        - Interpolated quantity at O in V-wind eqn
!                        for split physics (if LSLPHY=.T. only).
!          PTP9        - Interpolated quantity at O in T eqn
!                        for split physics (if LSLPHY=.T. only).
!          PGFLP9      - Interpolated quantity at O in unified_treatment
!                        grid-point fields
!          PSPPTCLEARP9- Interpolated clear sky radiation tendecy of T at O

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Called by GP_MODEL.

!     Author.
!     -------
!        Deborah Salmond 

! Modifications
! -------------
!   K. Yessad 07-02-2007: Splitting alti/surf for (gw) in NH+LGWADV.
!   K. Yessad 07-03-2007: Remove useless (gw)_surf interpolations in NH+LGWADV.
!   Modified 28-Aug-2007 F. Vana: removing 4 points splines arrays
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   26-Aug-2008 J. Masek and F. Vana   New dataflow for SLHD interpolators.
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   K. Yessad Dec 2008: merge SLCOMM+SLCOMM1 -> SLCOMM.
!   K. Yessad Dec 2008: merge the different (E)SLEXTPOL.. -> (E)SLEXTPOL.
!   K. Yessad (Aug 2009): always use root (QX,QY) for (p,q) variables names
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!   G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!   G. Mozdzynski (Aug 2011): support higher order interpolation
!   M. Diamantakis (June 2012): add code for quasi-monotone mass fixer 
!   G. Mozdzynski (May 2012): further cleaning
!   M. Diamantakis (Feb 2013): simplify mass fixers to avoid extra arg PCUB
!   T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): Move some variables.
!   R. El Khatib 27-Jul-2016 interface to eslextpol
!   K. Yessad (March 2017): simplify level numbering in interpolator.
!   F. Vana July 2018: RK4 scheme for trajectory research.
!   F. Vana  24-10-2018: Better control of SLPHYS interpolation
!   F. Vana 20-Feb-2019: quintic ivertical interpolation for RHS
!   R. El Khatib 22-May-2019 NGPTOT_CAP
!   M. Diamantakis Feb 2021: code for cartesian geometry option LSLDP_XYZ for departure point iterations
!   R. El Khatib 27-Oct-2021 optimize memory allocations for NH.
!   F. Vana and S. Lang January 2023: SPPT updated to new seq. physics order
! End Modifications
!     ------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0       , ONLY : LRPLANE
USE YOMMP0       , ONLY : NPROC, LSLDEBUG, LSYNC_CALSL
USE YOMMASK      , ONLY : NFIXSFLD
USE EINT_MOD     , ONLY : SL_STRUCT
USE MPL_MODULE   , ONLY : MPL_BARRIER

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEO
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM,YDGEO%YRDIM%NGPBLKS) 
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(INOUT) :: KVSEPC(YDGEO%YRDIM%NGPBLKS) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KVSEPL(YDGEO%YRDIM%NGPBLKS) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KL0(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,0:3,YDGEO%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(YDSL%NASLB1,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDGEO%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2,YDGEO%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWRL9(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDGEO%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDGMV%YT1%NDIM,YDGEO%YRDIM%NGPBLKS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1S(YDGEO%YRDIM%NPROMA,YDGMV%YT1%NDIMS,YDGEO%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLT1(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1,YDGEO%YRDIM%NGPBLKS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLPC(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIMPC,YDGEO%YRDIM%NGPBLKS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUP9(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDGEO%YRDIM%NGPBLKS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVP9(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDGEO%YRDIM%NGPBLKS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTP9(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDGEO%YRDIM%NGPBLKS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGFLP9(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIMSLP,YDGEO%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPPTCLEARP9(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDGEO%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLSCAW(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDMODEL%YRML_DYN%YYTLSCAW%NDIM,YDGEO%YRDIM%NGPBLKS)
! * Indices for interpolation of full-level data.
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLH0 (YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,0:3,YDGEO%YRDIM%NGPBLKS)
! * Interpolation weights and indices for interpolation of layer quantities
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRSCAW(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDMODEL%YRML_DYN%YYTRSCAW%NDIM,YDGEO%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSCO(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDMODEL%YRML_DYN%YYTSCO%NDIM,YDGEO%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCCO(YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDMODEL%YRML_DYN%YYTCCO%NDIM,YDGEO%YRDIM%NGPBLKS)
! * Interpolation weights and indices for interp. of interlayer quantities
! (neither in nor out : just allocated/deallocated before for a better computational performance.
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLSCAWH(YDGEO%YRDIM%NPROMNH,YDGEO%YRDIMV%NFLEVG,YDMODEL%YRML_DYN%YYTLSCAWH%NDIM,YDGEO%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRSCAWH(YDGEO%YRDIM%NPROMNH,YDGEO%YRDIMV%NFLEVG,YDMODEL%YRML_DYN%YYTRSCAWH%NDIM,YDGEO%YRDIM%NGPBLKS)

!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JSTGLO, IST, IEND, IBL

LOGICAL :: LLINC, LLSLPHY, LLDIAB, LLEXCLUDE
INTEGER(KIND=JPIM) :: IDUMARR(2)

!                                                       for interp. grid.
! * Indices for interpolation of half-level (interface) data.
INTEGER(KIND=JPIM) :: IL0H(YDGEO%YRDIM%NPROMNH,YDGEO%YRDIMV%NFLEVG,0:3,YDGEO%YRDIM%NGPBLKS)

! * interpolated quantities arrays.
REAL(KIND=JPRB) :: ZUF (YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDGEO%YRDIM%NGPBLKS)       ! cf. PUF in LAPINEA.
REAL(KIND=JPRB) :: ZVF (YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDGEO%YRDIM%NGPBLKS)       ! cf. PVF in LAPINEA.

! WENO indicator for the boundary
INTEGER(KIND=JPIM) :: INOWENO (YDGEO%YRDIM%NPROMA,YDGEO%YRDIMV%NFLEVG,YDGEO%YRDIM%NGPBLKS)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "eslextpol.intfb.h"
#include "lapinea.intfb.h"
#include "lapineb.intfb.h"
#include "slcomm.intfb.h"
#include "slcomm2a.intfb.h"
#include "slextpol.intfb.h"
#include "check_sl_struct.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CALL_SL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEO%YRDIM,YDDIMV=>YDGEO%YRDIMV,YDGEM=>YDGEO%YRGEM,  YDCSGEOM=>YDGEO%YRCSGEOM, &
 & YDGSGEOM=>YDGEO%YRGSGEOM, YDVAB=>YDGEO%YRVAB, YDVETA=>YDGEO%YRVETA, YDVFE=>YDGEO%YRVFE, &
 & YDSTA=>YDGEO%YRSTA, YDLAP=>YDGEO%YRLAP, &
 & YDCSGLEG=>YDGEO%YRCSGLEG, &
 & YDSPGEOM=>YDGEO%YSPGEOM,  &
 & YDDYN=>YDMODEL%YRML_DYN%YRDYN,YDDYNA=>YDMODEL%YRML_DYN%YRDYNA, &
 & YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,YDPTRSLB2=>YDMODEL%YRML_DYN%YRPTRSLB2, &
 & YDPTRSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1, YYTCCO=>YDMODEL%YRML_DYN%YYTCCO, &
 & YYTLSCAWH=>YDMODEL%YRML_DYN%YYTLSCAWH, YYTRSCAWH=>YDMODEL%YRML_DYN%YYTRSCAWH, &
 & YDECUCONVCA=>YDMODEL%YRML_PHY_EC%YRECUCONVCA,YGFL=>YDMODEL%YRML_GCONF%YGFL,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,  &
 & YDGPDDH=>YDMODEL%YRML_DIAG%YRGPDDH)

ASSOCIATE(NDIM=>YGFL%NDIM, NDIM1=>YGFL%NDIM1, NDIMPC=>YGFL%NDIMPC, &
 & NUMFLDS=>YGFL%NUMFLDS, YQ=>YGFL%YQ, &
 & NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, NPROMNH=>YDDIM%NPROMNH, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLSUL=>YDDIMV%NFLSUL, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, LSLDP_RK=>YDDYN%LSLDP_RK, &
 & LCA_ADVECT=>YDECUCONVCA%LCA_ADVECT, LCUCONV_CA=>YDECUCONVCA%LCUCONV_CA, &
 & RLATDEP=>YDECUCONVCA%RLATDEP, RLONDEP=>YDECUCONVCA%RLONDEP, RSAVEDP=>YDDYN%RSAVEDP, &
 & LAGPHY=>YDEPHY%LAGPHY, LEPHYS=>YDEPHY%LEPHYS, LSLPHY=>YDEPHY%LSLPHY,&
 & LMPHYS=>YDPHY%LMPHYS, &
 & NGPTOT=>YDGEM%NGPTOT, NGPTOT_CAP=>YDGEM%NGPTOT_CAP, &
 & YT1=>YDGMV%YT1, &
 & GFLTNDSL_DDH=>YDGPDDH%GFLTNDSL_DDH, GMVTNDSI_DDH=>YDGPDDH%GMVTNDSI_DDH, &
 & GMVTNDSL_DDH=>YDGPDDH%GMVTNDSL_DDH, &
 & MSLB1UR0=>YDPTRSLB1%MSLB1UR0, MSLB1WR0=>YDPTRSLB1%MSLB1WR0, MSLB1WR00=>YDPTRSLB1%MSLB1WR00, &
 & NFLDSLB1=>YDPTRSLB1%NFLDSLB1, RPARSL1=>YDPTRSLB1%RPARSL1, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS.
!              ---------------------------

! * Lagged physics is called.
LLDIAB = (LMPHYS.OR.LEPHYS) .AND. LAGPHY
LLSLPHY= LSLPHY.AND.LLDIAB
LLINC=.FALSE.
LLEXCLUDE=.TRUE.

!     ------------------------------------------------------------------

!*       2.    FIRST PART OF PROCESSOR COMMUNICATIONS WHEN NEEDED.
!              ---------------------------------------------------

IF( .NOT.(YDDYNA%LPC_CHEAP.AND.(NCURRENT_ITER > 0)) ) THEN

  CALL GSTATS(8,2)
  CALL GSTATS(51,0)

  IF (NPROC > 1) THEN

    IF(YDSL%LSLONDEM)THEN
      YDSL%LSLONDEM_ACTIVE=.TRUE.
      NFIXSFLD(1)=MSLB1UR0
      IF (LSLDP_RK) THEN
        NFIXSFLD(2)=MSLB1WR00+NFLEVG+NFLSUL
      ELSE
        NFIXSFLD(2)=MSLB1WR0+NFLEVG+NFLSUL
      ENDIF
    ENDIF

    IF (YDSL%LSLONDEM_ACTIVE) THEN
      CALL SLCOMM(YDSL,NFIXSFLD,NFLDSLB1,LLINC,1,PB1)
    ELSE
      CALL SLCOMM(YDSL,IDUMARR,NFLDSLB1,LLINC,0,PB1)
    ENDIF

  ENDIF

  IF (YDSL%LSLONDEM_ACTIVE) THEN
    YDSL%MASK_SL2=0
    IF (.NOT.LRPLANE) THEN
      CALL SLEXTPOL(YDGEO%YRDIM,YDSL,NFLDSLB1,NFIXSFLD,3,RPARSL1,PB1,KMASK2=YDSL%MASK_SL2)  
      IF( LSLDEBUG) CALL CHECK_SL_STRUCT(YDSL,'SLEXTPOL')
    ELSE
      CALL ESLEXTPOL(YDGEO,YDSL,NFLDSLB1,NFIXSFLD,3,PB1,KMASK2=YDSL%MASK_SL2)  
      IF (LSLDEBUG) THEN
        CALL ABOR1('CALL_SL: call to CHECK_SL_STRUCT not ready for plane geometry')
      ENDIF
    ENDIF
  ELSE
    IF (.NOT.LRPLANE) THEN
      CALL SLEXTPOL(YDGEO%YRDIM,YDSL,NFLDSLB1,IDUMARR,1,RPARSL1,PB1)  
      IF( LSLDEBUG) CALL CHECK_SL_STRUCT(YDSL,'SLEXTPOL')
    ELSE
      CALL ESLEXTPOL(YDGEO,YDSL,NFLDSLB1,IDUMARR,1,PB1)  
      IF (LSLDEBUG) THEN
        CALL ABOR1('CALL_SL: call to CHECK_SL_STRUCT not ready for plane geometry')
      ENDIF
    ENDIF
  ENDIF

  CALL GSTATS(51,1)
  CALL GSTATS(37,0)

ENDIF

!     ------------------------------------------------------------------

!*       3.    SL-TRAJECTORY RESEARCH, WEIGHTS AND INTERP. GRID CALCULATION.
!              -------------------------------------------------------------

IF( .NOT.(YDDYNA%LPC_CHEAP.AND.(NCURRENT_ITER > 0)) ) THEN

  CALL GSTATS(1004,0)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)&
!$OMP&PRIVATE(JSTGLO,IST,IEND,IBL)
  DO JSTGLO=1,NGPTOT,NPROMA
    IST=1
    IEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
    IBL=(JSTGLO-1)/NPROMA+1

    CALL LAPINEA(&
      ! --- INPUT --------------------------------------------------------------
      & YDGEO,YDMODEL%YRML_GCONF,YDMODEL%YRML_DYN,IST,IEND,YDSL,&
      & IBL,PB1(1,1),PB2(1,1,IBL),PWRL9(1,1,IBL),&
      ! --- INPUT/OUTPUT -------------------------------------------------------
      & KVSEPC(IBL),KVSEPL(IBL),&
      ! --- OUTPUT -------------------------------------------------------------
      & RSAVEDP(1,1,1,IBL),PCCO(1,1,1,IBL),ZUF(1,1,IBL),ZVF(1,1,IBL),&
      & KL0(1,1,0,IBL),KLH0(1,1,0,IBL),PLSCAW(1,1,1,IBL),PRSCAW(1,1,1,IBL),&
      & IL0H(1,1,0,IBL),PLSCAWH(1,1,1,IBL),PRSCAWH(1,1,1,IBL),&
      & PSCO(1,1,1,IBL),PGFLT1(1,1,1,IBL),INOWENO(1,1,IBL))

  ENDDO

!$OMP END PARALLEL DO
  CALL GSTATS(1004,1)

  IF (LCUCONV_CA .AND. LCA_ADVECT) THEN
    RLONDEP(1:NPROMA,1:NFLEVG,1:NGPBLKS)=PCCO(1:NPROMA,1:NFLEVG,YYTCCO%M_RLON,1:NGPBLKS)
    RLATDEP(1:NPROMA,1:NFLEVG,1:NGPBLKS)=PCCO(1:NPROMA,1:NFLEVG,YYTCCO%M_RLAT,1:NGPBLKS)
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       4.    SECOND PART OF PROCESSOR COMMUNICATIONS WHEN NEEDED.
!              ----------------------------------------------------

IF( .NOT.(YDDYNA%LPC_CHEAP.AND.(NCURRENT_ITER > 0)) .AND. (&
 & NPROC > 1 .AND. YDSL%LSLONDEM_ACTIVE)) THEN

  LLINC=.FALSE.
  CALL GSTATS(37,2)
  CALL GSTATS(89,0)

  CALL SLCOMM2A(YDGEO%YRDIM,YDSL,NFIXSFLD,LLEXCLUDE,&
   & YDSL%MASK_SL2,&
   & NFLDSLB1,LLINC,PB1)

  CALL GSTATS(89,1)
  CALL GSTATS(37,3)

  IF (.NOT.LRPLANE) THEN
    CALL SLEXTPOL(YDGEO%YRDIM,YDSL,NFLDSLB1,NFIXSFLD,2,RPARSL1,PB1)  
    IF( LSLDEBUG) CALL CHECK_SL_STRUCT(YDSL,'SLEXTPOL')
  ELSE
    CALL ESLEXTPOL(YDGEO,YDSL,NFLDSLB1,NFIXSFLD,2,PB1)  
    IF (LSLDEBUG) THEN
      CALL ABOR1('CALL_SL: call to CHECK_SL_STRUCT not ready for plane geometry')
    ENDIF
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       5.    INTERPOLATIONS.
!              ---------------

CALL GSTATS(1005,0)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JSTGLO,IST,IEND,IBL)
DO JSTGLO=1,NGPTOT_CAP,NPROMA
  IST=1
  IEND=MIN(NPROMA,NGPTOT_CAP-JSTGLO+1)
  IBL=(JSTGLO-1)/NPROMA+1

  CALL LAPINEB(YDGEO,YDGMV,&
   ! --- INPUT -------------------------------------------------------------
   & YDMODEL,NFLDSLB2,IST,IEND,YDSL,JSTGLO,LLSLPHY,&
   & IBL,INOWENO(1,1,IBL),&
   & KL0(1,1,0,IBL),KLH0(1,1,0,IBL),PLSCAW(1,1,1,IBL),PRSCAW(1,1,1,IBL),&
   & IL0H(1,1,0,IBL),PLSCAWH(1,1,1,IBL),PRSCAWH(1,1,1,IBL),&
   & ZUF(1,1,IBL),ZVF(1,1,IBL),PB1(1,1),PGFL(1,1,YQ%MP,IBL),&
   ! --- INPUT/OUTPUT ------------------------------------------------------
   & PSCO(1,1,1,IBL),PCCO(1,1,1,IBL),&
   & PB2(1,1,IBL),PGMVT1(1,1,1,IBL),PGMVT1S(1,1,IBL),PGFLT1(1,1,1,IBL),&
   & PGFLPC(:,:,:,IBL),&
   & GMVTNDSL_DDH(:,:,:,IBL),GFLTNDSL_DDH(:,:,:,IBL),GMVTNDSI_DDH(:,:,:,IBL),&
   ! --- OUTPUT ------------------------------------------------------------
   & PUP9(1,1,IBL),PVP9(1,1,IBL),PTP9(1,1,IBL),PGFLP9(:,:,:,IBL),&
   & PSPPTCLEARP9(1,1,IBL))

ENDDO

!$OMP END PARALLEL DO
CALL GSTATS(1005,1)

!     ------------------------------------------------------------------

!*       6.    FINAL ACTIONS AND DEALLOCATIONS.
!              ---------------------------------

! Synchronize - removes load imbalancing
IF (LSYNC_CALSL) THEN
   CALL GSTATS(37,2)
   CALL MPL_BARRIER(CDSTRING='CALL_SL')
   CALL GSTATS(37,3)
ENDIF

IF( .NOT.(YDDYNA%LPC_CHEAP.AND.(NCURRENT_ITER > 0)) ) THEN
  CALL GSTATS(37,1)
  CALL GSTATS(8,3)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('CALL_SL',1,ZHOOK_HANDLE)
END SUBROUTINE CALL_SL
