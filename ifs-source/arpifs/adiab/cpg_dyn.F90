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
SUBROUTINE CPG_DYN(YDGEOMETRY,YDGMV,&
 ! --- INPUT ----------------------------------------------------------
 & YDMODEL,KST,KEND,PBETADT,PDT,&
 & PSLHDA,PSLHDD0,&
 & KIBL,POROGL,POROGM,&
 & PRE0,PRDELP0,PUVH0,&
 & PHI0,PHIF0,PRCP0,&
 & PCTY0,PKENE0,&
 & PLSM,PTSOL,PQSOL,&
 & PATND,PDBBC,PRDPHI,PGWT0,PGWFT0,PGWT9,&
 & PNHXT0,PNHXT9,&
 & PGFL,&
 ! --- INPUT/OUTPUT ---------------------------------------------------
 & KSETTLOFF,PGFLPC,PGWS,PGMV,PGMVS,&
 & PB1,PB2,PGMVT1,PGMVT1S,PGFLT1,PGMVTNDSI,&
 & PWRL95,PSAVEDP5,PTRAJ_SLAG)

! -------------------------------------------------------------
!**** *CPG_DYN* - Grid point calculations: non lagged dynamics.

!     Purpose.
!     --------
!           Grid point calculations non lagged dynamics.
!           + simple physics.

!           Abbreviation "vwv" stands for "vertical wind variable".

!**   Interface.
!     ----------
!        *CALL* *CPG_DYN(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        KST       : first element of work.
!        KEND      : last element of work.
!        PBETADT   : BETADT or 0 according to configuration.
!        PDT       : For a leap-frog scheme (three time level scheme):
!                     'dt' at the first time-step, '2 dt' otherwise.
!                    For a 2TL SL scheme: timestep 'dt'.
!        PSLHDA    : Scaling factor of the deformation in f(d) function
!                    (including the model resolution correction)
!        PSLHDD0   : Threshold for deformation tensor enhancement
!        KIBL      : index into YRCSGEOM/YRCSGEOM instances in YDGEOMETRY
!        POROGL    : zonal component of "grad(surf orography)"
!        POROGM    : meridian component of "grad(surf orography)"
!        PRE0      : hydrostatic pressure "prehyd" at half levels at time t.
!        PRDELP0   : 1/(hydrostatic pressure depth of layers) at t.
!        PUVH0     : horizontal wind at time t at half levels.
!        PHI0      : geopotential height "gz" at t at half levels.
!        PHIF0     : geopotential height "gz" at t at full levels.
!        PRCP0     : contains "cp", "R" and "Kap=R/Cp" at t.
!        PCTY0     : contains vertical velocities, vertical integral of divergence at t.
!        PKENE0    : kinetic energy at t.
!        PLSM      : land-sea mask.
!        PTSOL     : surface temperature.
!        PQSOL     : surface specific humidity.
!        PATND     : adiabatic Lagrangian tendencies.
!        PDBBC     : [D (Gw)_surf / Dt]_adiab .
!        PRDPHI    : HYD: not used.
!                    NHEE: contains pre/(R T prehyd [Delta log(prehyd)]) at t.
!                    NHQE: contains 1/(R Tt [Delta log(prehyd)]) at t.
!                    "R" is the version of R (may be Rdry or Rmoist) used in the definition of vertical divergence "dver".
!        PGWT0     : [Gw] at t (LGWADV=T only; levels: see in CPG_GP).
!        PGWFT0    : [Gw] at full levels at t (LGWADV=T only; levels: see in CPG_GP)
!        PGWT9     : [Gw] at t-dt (LGWADV=T only; levels: see in CPG_GP).
!        PNHXT0    : "X" at full levels at t, diagnosed in CPG_GP.
!        PNHXT9    : "X" at full levels at t-dt, diagnosed in CPG_GP.
!        PGFL      : unified_treatment grid-point fields at t.

!     INPUT/OUTPUT:
!     -------------
!        PGFLPC    : unified_treatment grid-point fields at t (3TL PC only)
!        PGWS      : [Gw]_surf at t (LRDBBC only).
!        PGMV      : GMV variables at t-dt and t.
!        PGMVS     : GMVS variables at t-dt and t.
!        PB1       : "SLBUF1" buffer.
!        PB2       : "SLBUF2" buffer.
!        PGMVT1    : GMV variables at t+dt.
!        PGMVT1S   : GMVS variables at t+dt.
!        PGFLT1    : t+dt term for the unified_treatment of grid-point fields
!                    (as input contains the non lagged physics).
!        PGMVTNDSI : GMV: tendency due to linear terms (for DDH).
!        PTRAJ_SLAG: Stored SL quantities for TL/AD models

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE documentation about Eulerian and semi-Lagrangian schemes.

! Author
! ------
!   09-Aug-2001 K. YESSAD after part 7 of CPG.

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   09-Sep-2008 J. Masek  Dataflow for flow deformation along pressure levels
!   K. Yessad (Dec 2008): remove dummy CDLOCK + cleanings
!   F. Vana 13-Jan-2009: Storing KAPPA for TL/AD of SLHD
!   K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP.
!   K. Yessad (Nov 2009): prune lpc_old.
!   K. Yessad (Dec 2009): deep-layer NH dynamics for LGWADV=F.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   E. Bazile  19-01-2011 : Modified CP_FORCING input.
!   K. Yessad (Nov 2011): various contributions, call CP_FORCING moved in CPG_GP.
!   F. Vana    19-03-2012 : supression of storing not available values
!   K. Yessad (Nov 2012): simplify testings.
!   F. Vana  28-Nov-2013 : Redesigned trajectory handling
!   F. Vana    13-Feb-2014 Specific buffer for KAPPAT.
!   M. Diamantakis (Dec 2013): code for LSETTLSVF=true
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): Move some variables.
!   F. Vana    21-Nov-2017: Option LSLDP_CURV
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   F. Vana July 2018: RK4 scheme for trajectory research.
!   F. Vana 11-Jul-2019: Option LRHS_CURV
! End Modifications
!------------------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV   , ONLY : TGMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE INTDYN_MOD,ONLY : YYTHW0, YYTCTY0, YYTRCP0
USE YOMTRAJ  , ONLY : TRAJ_SLAG_TYPE,  LPRTTRAJ, LTRAJSAVE, LTRAJSLAG
USE YOMCT3   , ONLY : NSTEP
USE YOMLUN   , ONLY : NULOUT

!    ------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)       ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)           ,INTENT(INOUT) :: YDGMV
TYPE(MODEL)          ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KST
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KEND
INTEGER(KIND=JPIM)   ,INTENT(OUT)   :: KSETTLOFF(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PBETADT
REAL(KIND=JPRB)      ,INTENT(IN)    :: PDT
REAL(KIND=JPRB)      ,INTENT(IN)    :: PSLHDA(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PSLHDD0(YDGEOMETRY%YRDIM%NPROMA)
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KIBL
REAL(KIND=JPRB)      ,INTENT(IN)    :: POROGL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)      ,INTENT(IN)    :: POROGM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PRE0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PRDELP0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PUVH0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTHW0%NDIM)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PHI0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PHIF0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(OUT)   :: PRCP0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP0%NDIM)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PCTY0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTCTY0%NDIM)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PKENE0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PLSM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PTSOL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PQSOL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PATND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
 & YDMODEL%YRML_DYN%YRDYNA%YYTTND%NDIM)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PDBBC(YDGEOMETRY%YRDIM%NPROMNH)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PRDPHI(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PGWT0(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PGWFT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PGWT9(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PNHXT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PNHXT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGFLPC(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIMPC)
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGWS(YDGEOMETRY%YRDIM%NPROMNH)
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PB1(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGMVT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%YT1%NDIM)
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGMVT1S(YDGEOMETRY%YRDIM%NPROMA,YDGMV%YT1%NDIMS)
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1)
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGMVTNDSI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
 & YDMODEL%YRML_DIAG%YRMDDH%NDIMSIGMV)
TYPE (TRAJ_SLAG_TYPE),INTENT(INOUT) :: PTRAJ_SLAG
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PWRL95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PSAVEDP5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,3)
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JGFL, JJOFF
REAL(KIND=JPRB) :: ZTSTEP


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "cpeuldyn.intfb.h"
#include "lacdyn.intfb.h"




!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG_DYN',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDVAB=>YDGEOMETRY%YRVAB, &
 &  YDDYN=>YDMODEL%YRML_DYN%YRDYN,YDDYNA=>YDMODEL%YRML_DYN%YRDYNA, &
 & YDPTRSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1,YDPTRSLB2=>YDMODEL%YRML_DYN%YRPTRSLB2, &
 & YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH,YDRIP=>YDMODEL%YRML_GCONF%YRRIP, &
  & YGFL=>YDMODEL%YRML_GCONF%YGFL,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, &
  & YDPHLC=>YDMODEL%YRML_PHY_SLIN%YRPHLC,YDVDF=>YDMODEL%YRML_PHY_G%YRVDF)

ASSOCIATE(NDIM=>YGFL%NDIM, NDIM1=>YGFL%NDIM1, NDIMPC=>YGFL%NDIMPC, &
 & NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, YQ=>YGFL%YQ, &
 & NPROMA=>YDDIM%NPROMA, &
 & NPROMNH=>YDDIM%NPROMNH, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & L2TLFF=>YDDYN%L2TLFF, LADVF=>YDDYN%LADVF, LSETTLSVF=>YDDYN%LSETTLSVF, &
 & LSLHDHEAT=>YDDYN%LSLHDHEAT, LSLDP_CURV=>YDDYN%LSLDP_CURV, LRHS_CURV=>YDDYN%LRHS_CURV, &
 & LSLDP_XYZ=>YDDYN%LSLDP_XYZ,LSLDP_RK=>YDDYN%LSLDP_RK, LSLDP_SAVE=>YDDYN%LSLDP_SAVE, &
 & LEPHYS=>YDEPHY%LEPHYS, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT1=>YDGMV%YT1, YT9=>YDGMV%YT9, &
 & NDIMSIGMV=>YDMDDH%NDIMSIGMV, &
 & LSPHLC=>YDPHLC%LSPHLC, &
 & MSLB1C9=>YDPTRSLB1%MSLB1C9, MSLB1GFL9=>YDPTRSLB1%MSLB1GFL9, &
 & MSLB1SP9=>YDPTRSLB1%MSLB1SP9, MSLB1T0=>YDPTRSLB1%MSLB1T0, &
 & MSLB1T9=>YDPTRSLB1%MSLB1T9, MSLB1U0=>YDPTRSLB1%MSLB1U0, &
 & MSLB1U9=>YDPTRSLB1%MSLB1U9, MSLB1UR0=>YDPTRSLB1%MSLB1UR0, &
 & MSLB1V0=>YDPTRSLB1%MSLB1V0, MSLB1V9=>YDPTRSLB1%MSLB1V9, &
 & MSLB1Z0=>YDPTRSLB1%MSLB1Z0, MSLB1Z9=>YDPTRSLB1%MSLB1Z9, &
 & MSLB1VR0=>YDPTRSLB1%MSLB1VR0, MSLB1WR0=>YDPTRSLB1%MSLB1WR0, &
 & MSLB1ZR0=>YDPTRSLB1%MSLB1ZR0, NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & MSLB1UR00=>YDPTRSLB1%MSLB1UR00, MSLB1VR00=>YDPTRSLB1%MSLB1VR00, &
 & MSLB1ZR00=>YDPTRSLB1%MSLB1ZR00, MSLB1WR00=>YDPTRSLB1%MSLB1WR00, &
 & MSLB2KAPPA=>YDPTRSLB2%MSLB2KAPPA, MSLB2KAPPAT=>YDPTRSLB2%MSLB2KAPPAT, &
 & MSLB2STDDISU=>YDPTRSLB2%MSLB2STDDISU,MSLB2STDDISV=>YDPTRSLB2%MSLB2STDDISV, &
 & MSLB2STDDISW=>YDPTRSLB2%MSLB2STDDISW, &
 & MSLB2URL=>YDPTRSLB2%MSLB2URL, MSLB2VRL=>YDPTRSLB2%MSLB2VRL, &
 & MSLB2WRL=>YDPTRSLB2%MSLB2WRL, MSLB2ZRL=>YDPTRSLB2%MSLB2ZRL, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2, TSTEP=>YDRIP%TSTEP, &
 & NSTART=>YDRIP%NSTART, LMPHYS=>YDPHY%LMPHYS)
!     ------------------------------------------------------------------

!*       2.    Semi-lagrangian scheme
!              ----------------------

!* Cases nsiter=0, predictor and corrector for lpc_full.

IF (YDDYNA%LSLAG) THEN

!*       2.1   Semi-lagrangian dynamics (not lagged part).

  CALL LACDYN(YDGEOMETRY,YDGMV,&
   ! --- INPUT ----------------------------------------------------------------
   & YDMODEL,KST,KEND,PBETADT,PDT,PSLHDA,PSLHDD0,&
   & KIBL,POROGL,POROGM,PHIF0,&
   & PNHXT9,PNHXT0,PRE0,PRDELP0,PCTY0,PUVH0,&
   & PATND,PDBBC,PRDPHI,PGWT0,PGWFT0,PGWT9,PGFL,&
   ! --- INPUT/OUTPUT ---------------------------------------------------------
   & KSETTLOFF,PGMV(1,1,1),PGMVS(1,1),PB1(1,1),PB2(1,1),PGMVT1(1,1,1),&
   & PGMVT1S(1,1),PGWS,PGMVTNDSI(1,1,1),PWRL95)

!*       2.2   Simple physics.

  IF ((.NOT.(LMPHYS.OR.LEPHYS)).AND.LSPHLC) THEN
    IF (YDDYNA%LTWOTL) THEN
      CALL VDIFLCZ(YDVDF,YDVAB,YDPHLC,KST,KEND,NPROMA,NFLEVG,PDT,PLSM,&
       & PHI0,PHIF0,PRCP0(1,1,YYTRCP0%M_CP),PRCP0(1,1,YYTRCP0%M_R),PRDELP0,&
       & PGMV(1,1,YT0%MU),PGMV(1,1,YT0%MV),PGMV(1,1,YT0%MT),PRE0,&
       & PGMV(1,1,YT0%MU),PGMV(1,1,YT0%MV),PGMV(1,1,YT0%MT),&
       & PGMVT1(1,1,YT1%MU),PGMVT1(1,1,YT1%MV),PGMVT1(1,1,YT1%MT))
    ELSE
      CALL VDIFLCZ(YDVDF,YDVAB,YDPHLC,KST,KEND,NPROMA,NFLEVG,PDT,PLSM,&
       & PHI0,PHIF0,PRCP0(1,1,YYTRCP0%M_CP),PRCP0(1,1,YYTRCP0%M_R),PRDELP0,&
       & PGMV(1,1,YT0%MU),PGMV(1,1,YT0%MV),PGMV(1,1,YT0%MT),PRE0,&
       & PGMV(1,1,YT9%MU),PGMV(1,1,YT9%MV),PGMV(1,1,YT9%MT),&
       & PGMVT1(1,1,YT1%MU),PGMVT1(1,1,YT1%MV),PGMVT1(1,1,YT1%MT))
    ENDIF
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       3.    Eulerian scheme.
!              ----------------

!* Cases nsiter=0, predictor and corrector for lpc_full.

IF (.NOT.YDDYNA%LSLAG) THEN

!*       3.1   Dynamics with Eulerian advection treatment

  CALL CPEULDYN(YDGEOMETRY,YDGMV,&
   ! --- INPUT ----------------------------------------------------
   & YDEPHY,YGFL,YDMODEL%YRML_DYN,YDPHY,KST,KEND,PBETADT,PDT,&
   & KIBL,&
   & PGFL,PTSOL,PQSOL,PNHXT0,PKENE0,PCTY0,&
   & PRE0,PRDELP0,PATND,PNHXT9,&
   ! --- INPUT/OUTPUT ---------------------------------------------
   & PGFLT1,PGMVT1(1,1,1),PGMVT1S(1,1),PGFLPC,PGMV(1,1,1),PGMVS(1,1),&
   ! --- OUTPUT ---------------------------------------------------
   & PB2(1,1))

!*       3.2   Simple physics.

  IF ((.NOT.(LMPHYS.OR.LEPHYS)).AND.LSPHLC) THEN
    CALL VDIFLCZ(YDVDF,YDVAB,YDPHLC,KST,KEND,NPROMA,NFLEVG,PDT,PLSM,&
     & PHI0,PHIF0,PRCP0(1,1,YYTRCP0%M_CP),PRCP0(1,1,YYTRCP0%M_R),PRDELP0,&
     & PGMV(1,1,YT0%MU),PGMV(1,1,YT0%MV),PGMV(1,1,YT0%MT),PRE0,&
     & PGMV(1,1,YT9%MU),PGMV(1,1,YT9%MV),PGMV(1,1,YT9%MT),&
     & PGMVT1(1,1,YT1%MU),PGMVT1(1,1,YT1%MV),PGMVT1(1,1,YT1%MT))
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       4.   Stores some variables for AD/TL.
!             --------------------------------

IF (LTRAJSAVE .AND. LTRAJSLAG) THEN
  IF (TSTEP == 0.0_JPRB) THEN
    ZTSTEP=1.0_JPRB
  ELSE
    ZTSTEP=TSTEP
  ENDIF
  PTRAJ_SLAG%PUL95(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1U9+1-NFLSA:MSLB1U9+NFLEVG-NFLSA)
  PTRAJ_SLAG%PVL95(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1V9+1-NFLSA:MSLB1V9+NFLEVG-NFLSA)
  IF (LRHS_CURV) THEN
    PTRAJ_SLAG%PZL95(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1Z9+1-NFLSA:MSLB1Z9+NFLEVG-NFLSA)
  ENDIF
  PTRAJ_SLAG%PTL95(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1T9+1-NFLSA:MSLB1T9+NFLEVG-NFLSA)
  DO JGFL=1,NUMFLDS
    IF(YCOMP(JGFL)%LADV) THEN
      JJOFF=(NFLEN-NFLSA+1)*(YCOMP(JGFL)%MP_SL1-1)-NFLSA
      PTRAJ_SLAG%PGFLL95(KST:KEND,1:NFLEVG,YCOMP(JGFL)%MP_SL1) =&
         & PB1(KST:KEND,MSLB1GFL9+1+JJOFF:MSLB1GFL9+NFLEVG+JJOFF)
      ! store absolute position as a label for recognition in the TL/AD code
      PTRAJ_SLAG%IGFLL95(YCOMP(JGFL)%MP_SL1) = YCOMP(JGFL)%MP
    ENDIF
  ENDDO
  ! Some fields are tendencies which include a TSTEP factor which may
  ! be different when fields will be read in.
  PTRAJ_SLAG%PCL95(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1C9+1-NFLSA:MSLB1C9+NFLEVG-NFLSA)/ZTSTEP
  PTRAJ_SLAG%PUL05(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1U0+1-NFLSA:MSLB1U0+NFLEVG-NFLSA)/ZTSTEP
  PTRAJ_SLAG%PVL05(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1V0+1-NFLSA:MSLB1V0+NFLEVG-NFLSA)/ZTSTEP
  IF (LRHS_CURV) THEN
    PTRAJ_SLAG%PZL05(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1Z0+1-NFLSA:MSLB1Z0+NFLEVG-NFLSA)/ZTSTEP
  ENDIF
  PTRAJ_SLAG%PTL05(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1T0+1-NFLSA:MSLB1T0+NFLEVG-NFLSA)/ZTSTEP
  IF (LADVF.AND.L2TLFF) THEN
    PTRAJ_SLAG%PUT15(KST:KEND,1:NFLEVG)=PGMVT1(KST:KEND,1:NFLEVG,YT1%MU)
    PTRAJ_SLAG%PVT15(KST:KEND,1:NFLEVG)=PGMVT1(KST:KEND,1:NFLEVG,YT1%MV)
  ENDIF
  PTRAJ_SLAG%PURL05(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1UR0+1-NFLSA:MSLB1UR0+NFLEVG-NFLSA)
  PTRAJ_SLAG%PVRL05(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1VR0+1-NFLSA:MSLB1VR0+NFLEVG-NFLSA)
  IF (LSLDP_CURV) THEN
    PTRAJ_SLAG%PZRL05(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1ZR0+1-NFLSA:MSLB1ZR0+NFLEVG-NFLSA)
  ENDIF
  PTRAJ_SLAG%PWRL05(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1WR0+1-NFLSA:MSLB1WR0+NFLEVG-NFLSA)
  IF (LSLDP_RK) THEN
    PTRAJ_SLAG%PURL005(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1UR00+1-NFLSA:MSLB1UR00+NFLEVG-NFLSA)
    PTRAJ_SLAG%PVRL005(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1VR00+1-NFLSA:MSLB1VR00+NFLEVG-NFLSA)
    IF (LSLDP_CURV) THEN
      PTRAJ_SLAG%PZRL005(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1ZR00+1-NFLSA:MSLB1ZR00+NFLEVG-NFLSA)
    ENDIF
    PTRAJ_SLAG%PWRL005(KST:KEND,1:NFLEVG)=PB1(KST:KEND,MSLB1WR00+1-NFLSA:MSLB1WR00+NFLEVG-NFLSA)
  ENDIF
  PTRAJ_SLAG%PURL5(KST:KEND,1:NFLEVG)=PB2(KST:KEND,MSLB2URL:MSLB2URL+NFLEVG-1)
  PTRAJ_SLAG%PVRL5(KST:KEND,1:NFLEVG)=PB2(KST:KEND,MSLB2VRL:MSLB2VRL+NFLEVG-1)
  PTRAJ_SLAG%PWRL5(KST:KEND,1:NFLEVG)=PB2(KST:KEND,MSLB2WRL:MSLB2WRL+NFLEVG-1)
  IF (LSLDP_XYZ) PTRAJ_SLAG%PZRL5(KST:KEND,1:NFLEVG)=PB2(KST:KEND,MSLB2ZRL:MSLB2ZRL+NFLEVG-1)
  IF (YDDYNA%LSLHD.AND.(.NOT.YDDYNA%LSLHD_STATIC)) THEN
    PTRAJ_SLAG%PKAPPA5(KST:KEND,1:NFLEVG)=PB2(KST:KEND,MSLB2KAPPA:MSLB2KAPPA+NFLEVG-1)
    IF (LSLHDHEAT) THEN
      PTRAJ_SLAG%PKAPPAT5(KST:KEND,1:NFLEVG)=PB2(KST:KEND,MSLB2KAPPAT:MSLB2KAPPAT+NFLEVG-1)
    ENDIF
  ENDIF

! COMAD alpha parameter --> NEEDED for AD
  IF (YDDYNA%LCOMAD) THEN
    PTRAJ_SLAG%PSTDDISU5(KST:KEND,1:NFLEVG)=PB2(KST:KEND,MSLB2STDDISU:MSLB2STDDISU+NFLEVG-1)
    PTRAJ_SLAG%PSTDDISV5(KST:KEND,1:NFLEVG)=PB2(KST:KEND,MSLB2STDDISV:MSLB2STDDISV+NFLEVG-1)
    PTRAJ_SLAG%PSTDDISW5(KST:KEND,1:NFLEVG)=PB2(KST:KEND,MSLB2STDDISW:MSLB2STDDISW+NFLEVG-1)
  ENDIF

  IF (LSETTLSVF.OR.(LSLDP_XYZ.AND.LSLDP_SAVE)) THEN
    PTRAJ_SLAG%PWRL95(KST:KEND,1:NFLEVG)=PWRL95(KST:KEND,1:NFLEVG)
  ENDIF

  IF (NSTEP>NSTART.AND.LSLDP_XYZ.AND.LSLDP_SAVE) THEN
    PTRAJ_SLAG%PSAVEDP5(KST:KEND,1:NFLEVG,1:3)=PSAVEDP5(KST:KEND,1:NFLEVG,1:3)
  ENDIF

  PTRAJ_SLAG%PX95(KST:KEND)=PB1(KST:KEND,MSLB1SP9)

  IF (LPRTTRAJ.AND.PTRAJ_SLAG%LASTCHUNK) WRITE(NULOUT,*)'GREPTRAJ STORE TRAJ_SLAG in CPG_DYN'
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_DYN',1,ZHOOK_HANDLE)
END SUBROUTINE CPG_DYN
