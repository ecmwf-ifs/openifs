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

SUBROUTINE LACDYN(YDGEOMETRY,YDGMV,&
 ! --- INPUT ------------------------------------------------------------------
 & YDMODEL,KST,KPROF,PBETADT,PDT,PSLHDA,PSLHDD0,&
 & KIBL,POROGL,POROGM,PHIF0,&
 & PNHXT9,PNHXT0,&
 & PRES0,PRDELP0,PCTY0,PUVH0,&
 & PATND,PDBBC,PRDPHI,PGWT0,PGWFT0,PGWT9,PGFL,&
 ! --- INPUT/OUTPUT -----------------------------------------------------------
 & KSETTLOFF,PGMV,PGMVS,PB1,PB2,PGMVT1,PGMVT1S,PGWS,PGMVTNDSI,PWRL95)

!**** *LACDYN*   Semi-Lagrangian scheme.
!                Computation of the t and t-dt useful quantities at grid-points.

!     Purpose.
!     --------

!          Dynamic non-linear computations in grid-point space
!          for hydrostatic and NH primitive equations and SL scheme.

!          Additional remarks:
!          - notation "prehyds" is for hydrostatic surface pressure.
!          - for input and output upper air variables, values are at full levels
!            if no other information is provided.
!          - NH variables: "P variable" is the pressure departure variable
!            (pressure departure equation), i.e. ln(pre/prehyd) if npdvar=2;
!            abbreviation "vwv" stands for "vertical wind variable".
!          - for PC schemes:
!            * this routine is called for nsiter=0.
!            * this routine is called for the predictor of lpc_full.
!            * this routine is called for the corrector of lpc_full.

!          This subroutine fills the semi-Lagrangian buffers to be interpolated.

!**   Interface.
!     ----------
!        *CALL* *LACDYN(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST         - first element of work.
!          KPROF       - depth of work.
!          PBETADT     - BETADT or 0 according to configuration.
!          PDT         - For a leap-frog scheme (three time level scheme):
!                         'dt' at the first time-step, '2 dt' otherwise.
!                        For a 2TL SL scheme: timestep 'dt'.
!          PSLHDA      - Scaling factor of the deformation in f(d) function
!                        (including the model resolution correction)
!          PSLHDD0     - Threshold for deformation tensor enhancement
!          KIBL        - index into YDGSGEOM instance in YDGEOMETRY
!          POROGL      - zonal component of the orography gradient
!          POROGM      - merid component of the orography gradient
!          PHIF0       - geop. height at full levels
!          PNHXT9      - "X" at full levels at t-dt, diagnosed in CPG_GP.
!          PNHXT0      - "X" at full levels at t, diagnosed in CPG_GP.
!          PRES0       - hydrostatic pressure "prehyd" at half levels at t.
!          PRDELP0     - 1/(pressure depth of layers) at t.
!          PCTY0       - contains vertical velocities, vertical integral of divergence at t.
!          PUVH0       - horizontal wind at time t at half levels.
!          PATND       - adiabatic Lagrangian tendencies.
!          PDBBC       - [D (Gw)_surf / Dt]_adiab .
!          PRDPHI      - HYD: not used.
!                        NHEE: contains pre/(R T prehyd [Delta log(prehyd)]) at t.
!                        NHQE: contains 1/(R Tt [Delta log(prehyd)]) at t.
!                        "R" is the version of R (may be Rdry or Rmoist) used in the definition of vertical divergence "dver".
!          PGWT0       - [Gw] at t (LGWADV=T only; levels: see CPG_GP).
!          PGWFT0      - [Gw] at full-levels at t (LGWADV=T only; levels: see CPG_GP).
!          PGWT9       - [Gw] at t-dt (LGWADV=T only; levels: see CPG_GP).
!          PGFL        - unified_treatment grid-point fields

!        INPUT/OUTPUT:
!          KSETTLOFF   - counter for SETTLSTF=T (# of points new scheme activated at each lev);
!                        counter for LMIXETTLS (only for LSETTLSTF=F)
!          PGMV        - GMV variables at t-dt and t.
!          PGMVS       - GMVS variables at t-dt and t.
!          PB1         - "SLBUF1" buffer for interpolations.
!          PB2         - "SLBUF2" buffer.
!          PGMVT1      - GMV variables at t+dt.
!          PGMVT1S     - GMVS variables at t+dt.
!          PGWS        - [Gw]_surf at t (LRDBBC only).
!          PGMVTNDSI   - GMV: tendency due to linear terms (for DDH).
!          PWRL95      - store current timestep PGMV(,,YMEDOT%YT9) values before re-setting

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Called by CPG_DYN.

!     Reference.
!     ----------
!             Arpege documentation about semi-lagrangian scheme.

!     Author.
!     -------
!      K. YESSAD (METEO FRANCE/CNRM/GMAP) after routines
!        CPLGDY1 and LAGSIMP coded by Maurice IMBARD.
!      Original : FEBRUARY 1992.

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!   09-Sep-2008 J. Masek  Dataflow for flow deformation along pressure levels
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP.
!   K. Yessad (Nov 2009): prune lpc_old.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   K. Yessad (Dec 2011): various contributions.
!   M. Diamantakis (Feb 2014): code for LSETTLSVF=T
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   S. Malardel (Nov 2013): LATTE_STDDIS for COMAD corrections
!   K. Yessad (July 2014): Move some variables.
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
!   F. Vana    21-Nov-2017: Option LSLDP_CURV
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   J. Vivoda (July 2018): mixed NESC/SETTLS scheme.
!   F. Vana   14-Sep-2020: SLHD interacts with LSETTLSVF
! End Modifications
!------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE INTDYN_MOD   , ONLY : YYTHW0, YYTCTY0

!   ---------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(MODEL)       ,INTENT(IN)    :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KSETTLOFF(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBETADT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDA(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDD0(YDGEOMETRY%YRDIM%NPROMA)
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHIF0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHXT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHXT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRES0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCTY0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTCTY0%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUVH0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTHW0%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PATND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
 & YDMODEL%YRML_DYN%YRDYNA%YYTTND%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDBBC(YDGEOMETRY%YRDIM%NPROMNH)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDPHI(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWT0(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWFT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWT9(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%YT1%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1S(YDGEOMETRY%YRDIM%NPROMA,YDGMV%YT1%NDIMS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGWS(YDGEOMETRY%YRDIM%NPROMNH)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVTNDSI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_DIAG%YRMDDH%NDIMSIGMV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWRL95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
!     ------------------------------------------------------------------
! * computed in LASURE:
REAL(KIND=JPRB) :: ZBT, ZDTS2, ZESGM, ZESGP
REAL(KIND=JPRB) :: ZREDIV(YDGEOMETRY%YRDIM%NPROMA)
! * computed in LASSIE or LANHSI:
REAL(KIND=JPRB) :: ZBDT(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZGAGT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZGAGT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZGAGT9L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZGAGT9M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZLPD0(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZLPD9(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSDIV0(YDGEOMETRY%YRDIM%NPROMA),ZSDIV9(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZSPDS0 (YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZSPDS9(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTOD0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZTOD9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMIXNL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZALPHA(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

LOGICAL :: LL2TLFF1, LLGWADV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "lanhsi.intfb.h"
#include "lanhqesi.intfb.h"
#include "lanhsib.intfb.h"
#include "lassie.intfb.h"
#include "lasure.intfb.h"
#include "latte_bbc.intfb.h"
#include "latte_kappa.intfb.h"
#include "latte_nl.intfb.h"
#include "latte_stddis.intfb.h"
#include "lattes.intfb.h"
#include "lattex.intfb.h"
#include "lavabo.intfb.h"
#include "lavent.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LACDYN',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), &
 & YDDYN=>YDMODEL%YRML_DYN%YRDYN,YDDYNA=>YDMODEL%YRML_DYN%YRDYNA, &
 & YDPTRSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1, YDPTRSLB2=>YDMODEL%YRML_DYN%YRPTRSLB2, &
 & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
 & YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH, &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL)

ASSOCIATE(NDIM=>YGFL%NDIM, &
 & NPROMA=>YDDIM%NPROMA, &
 & NPROMNH=>YDDIM%NPROMNH, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT1=>YDGMV%YT1, &
 & NDIMSIGMV=>YDMDDH%NDIMSIGMV)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS:
!              ----------------------------

CALL LASURE(YDGEOMETRY,YDMODEL%YRML_PHY_EC%YREPHY,YDDYN,YDMODEL%YRML_DYN%YREDYN,YDMODEL%YRML_PHY_MF%YRPHY,KST,KPROF,PBETADT,PDT, &
 & KIBL,ZDTS2,ZBT,LL2TLFF1,ZBDT,ZREDIV,ZESGP,ZESGM)

!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF THE LINEAR TERMS FOR SEMI-IMPLICIT SCHEME.
!              ---------------------------------------------------------

IF (YDDYNA%LNHEE) THEN
  CALL LANHSI(YDGEOMETRY,YDGMV,YGFL,YDDYN,YDDYNA,KST,KPROF,YDDYNA%LGWADV,YDGSGEOM%RCORI,ZREDIV,PGMV,PGMVS,PGFL,&
   & ZSDIV0,ZSDIV9,ZTOD0,ZTOD9,ZGAGT0L,ZGAGT0M,ZGAGT9L,ZGAGT9M,&
   & ZSPDS0,ZSPDS9,ZLPD0,ZLPD9)
ELSEIF (YDDYNA%LNHQE) THEN
  CALL LANHQESI(YDGEOMETRY,YDGMV,YGFL,YDDYN,YDDYNA,KST,KPROF,YDDYNA%LGWADV,YDGSGEOM%RCORI,PGMV,PGMVS,PGFL,&
   & ZSDIV0,ZSDIV9,ZTOD0,ZTOD9,ZGAGT0L,ZGAGT0M,ZGAGT9L,ZGAGT9M,&
   & ZLPD0,ZLPD9)
ELSE
  CALL LASSIE(YDGEOMETRY,YDGMV,YGFL,YDDYN,YDDYNA,KST,KPROF,YDGSGEOM%RCORI,PGMV,PGMVS,PGFL,&
   & ZSDIV0,ZSDIV9,ZTOD0,ZTOD9,ZGAGT0L,ZGAGT0M,ZGAGT9L,ZGAGT9M)
ENDIF

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF THE WIND COMPONENTS NECESSARY FOR SL TRAJECTORY.
!              ---------------------------------------------------------------

IF ((NCURRENT_ITER == 0).OR. (&
 & NCURRENT_ITER > 0 .AND. YDDYNA%LPC_FULL .AND.(.NOT.YDDYNA%LPC_CHEAP) )) THEN

  CALL LAVENT(YDGEOMETRY,YDGMV,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_DYN,YDMODEL%YRML_PHY_MF%YRPARAR,&
   & KST,KPROF,KIBL,KSETTLOFF,LL2TLFF1, &
   & ZDTS2,PRDELP0,&
   & PCTY0(1,0,YYTCTY0%M_EVEL),PATND,PGMV,PB1,PB2,PWRL95,ZALPHA)

ENDIF

!     ------------------------------------------------------------------

!*       4.    COMPUTATION OF THE 3D-EQUATIONS RIGHT-HAND SIDE TERMS.
!              ------------------------------------------------------

!        4.1:  compute what kind of extrapolation will be used in each grid
!              point and level for LMIXETTLS=T:

IF (YDDYNA%LMIXETTLS) THEN
  CALL LATTE_NL(YDGEOMETRY,YDGMV,YDMODEL%YRML_DYN,YDMODEL%YRML_GCONF%YRRIP,&
   & KST,KPROF,ZDTS2,ZBDT,PGMV(1,1,YT0%MU),PGMV(1,1,YT0%MV),PCTY0(1,0,YYTCTY0%M_EVEL),&
   & PRDELP0,ZMIXNL,KSETTLOFF)
ELSE
  ZMIXNL(:,:) = 1.0_JPRB
ENDIF

!        4.2:  general case:

CALL LATTEX(YDGEOMETRY,YDGMV,YDMODEL%YRML_DIAG%YRLDDH,YDMDDH,YDMODEL%YRML_GCONF,YDMODEL%YRML_DYN,KST,KPROF,ZDTS2,ZBT,ZBDT, &
 & ZESGP,ZESGM,&
 & KIBL,POROGL,POROGM,PHIF0,&
 & ZGAGT0L,ZGAGT0M,ZTOD0,ZSPDS0,ZLPD0,&
 & ZGAGT9L,ZGAGT9M,ZTOD9,ZSPDS9,ZLPD9,&
 & PRDELP0,PCTY0(1,0,YYTCTY0%M_EVEL),PATND,&
 & PNHXT0,PGWT0,PGWFT0,PNHXT9,PGWT9,PGFL,ZMIXNL,&
 & PGMV,PGMVT1,PGMVTNDSI,PB1,PB2)

!        4.3:  Additional quantities required in the NH model
!              for option LRDBBC=T:

IF (YDDYNA%LNHDYN.AND.YDDYNA%LRDBBC.AND.(.NOT.YDDYNA%LGWADV)) THEN
  ! * Remarks:
  !   - for LGWADV=T this calculation is not required because
  !     in this case [Gw]_surf(t+dt) is computed by a diagnostic relationship
  !     in routine LAPINEB.
  CALL LATTE_BBC(YDGEOMETRY,YDGMV,YDMODEL%YRML_DYN,KST,KPROF,ZDTS2,ZESGP,ZESGM,PDBBC,PRDPHI,&
   & ZMIXNL,PGWS,PGMVS,PGMV,PB1,PB2)
ENDIF

!     ------------------------------------------------------------------

!*       5.    COMPUTATION OF THE 2D-EQUATIONS RIGHT-HAND SIDE TERMS.
!              ------------------------------------------------------

CALL LATTES(YDGEOMETRY,YDGMV,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_DYN,KST,KPROF,ZDTS2,ZBDT,ZESGP,ZESGM,&
 & KIBL,POROGL,POROGM,&
 & ZSDIV0,ZSDIV9,PCTY0(1,0,YYTCTY0%M_PSDVBC),PRES0,PGMVS,ZMIXNL,&
 & PGMV,PGMVT1S,PB1,PB2)

!     ------------------------------------------------------------------

!*       6.    COMPUTATION OF "KAPPA" AND STORE IT IN  "SLBUF2".
!              -----------------------------------------------------------

IF (YDDYNA%LSLHD) THEN
  CALL LATTE_KAPPA(YDGEOMETRY,YDGMV,YDDYN,YDDYNA,YDPTRSLB2,KST,KPROF,ZDTS2,PSLHDA,PSLHDD0,PGMV,PGMVS,&
   & PUVH0,PRES0(1,NFLEVG),PRDELP0,ZALPHA,PB2)
ENDIF

!     ------------------------------------------------------------------

!*       7.    COMPUTATION OF "STDDIS" for each direction
!              AND STORE THEM IN  "SLBUF2".
!              -----------------------------------------------------------

IF (YDDYNA%LCOMAD) THEN
  CALL LATTE_STDDIS(YDGEOMETRY,YDDYNA,YDGMV,YDPTRSLB2,KST,KPROF,ZDTS2,PGMV,PB2)
ENDIF

!     ------------------------------------------------------------------

!*       8.    UPPER AND LOWER LATERAL BOUNDARIES CONDITIONS.
!              ----------------------------------------------

IF ((NCURRENT_ITER == 0).OR. (&
 & NCURRENT_ITER > 0 .AND. YDDYNA%LPC_FULL .AND.(.NOT.YDDYNA%LPC_CHEAP) )) THEN
  CALL LAVABO(YDGEOMETRY,YDMODEL,KST,KPROF,LL2TLFF1,PB1)
ENDIF

!     ------------------------------------------------------------------

!*       8.    SI TERM COMPUTED WITH d3 OR d4 VARIABLE FOR GWADV 
!              -------------------------------------------------
!*    compute dt L X(t) (+ Xterm)
!*    this quantity is substracted later in CPGLAG

IF (YDDYNA%LGWADV) THEN

  LLGWADV=.FALSE.
  IF (YDDYNA%LNHEE) THEN
    CALL LANHSI(YDGEOMETRY,YDGMV,YGFL,YDDYN,YDDYNA,KST,KPROF,LLGWADV,YDGSGEOM%RCORI,ZREDIV,PGMV,PGMVS,PGFL,&
     & ZSDIV0,ZSDIV9,ZTOD0,ZTOD9,ZGAGT0L,ZGAGT0M,ZGAGT9L,ZGAGT9M,&
     & ZSPDS0,ZSPDS9,ZLPD0,ZLPD9)
  ELSEIF (YDDYNA%LNHQE) THEN
    CALL LANHQESI(YDGEOMETRY,YDGMV,YGFL,YDDYN,YDDYNA,KST,KPROF,LLGWADV,YDGSGEOM%RCORI,PGMV,PGMVS,PGFL,&
     & ZSDIV0,ZSDIV9,ZTOD0,ZTOD9,ZGAGT0L,ZGAGT0M,ZGAGT9L,ZGAGT9M,&
     & ZLPD0,ZLPD9)
  ENDIF

  CALL LANHSIB(YDGEOMETRY,YDGMV,YDMODEL%YRML_DYN,KST,KPROF,ZBT,ZBDT,ZESGP,&
   & ZGAGT0L,ZGAGT0M,ZTOD0,ZSPDS0,ZLPD0,ZSDIV0,PNHXT0,ZMIXNL,&
   & PGMV,PB1,PB2)

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LACDYN',1,ZHOOK_HANDLE)
END SUBROUTINE LACDYN
