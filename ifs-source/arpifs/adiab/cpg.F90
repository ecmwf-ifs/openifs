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
SUBROUTINE CPG(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDTMP, YDCPG_TND, YDCPG_SL1, YDCPG_SL2, &
& YDCPG_MISC, YDCPG_GPAR, YDCPG_PHY0, YDCPG_PHY9, YDMF_PHYS, YDCPG_DDH, YDCPG_DYN0, YDCPG_DYN9, &
& YDMF_PHYS_SURF, YDVARS, YDMODEL, YDFIELDS, YDSL, PGFLSLP, PSAVTEND, &
& PGMVTNDHD_DDH, PGFLTNDHD_DDH, PGPSDT2D, PSD_PF, PPGTERM, PGMV, PGMVS, PGFL, PGFLPC, PGFLPT, &
& PSD_XA, PGMVTNDSI_DDH, PHYS_MWAVE, KSETTLOFF, PB1, PGMVT1, PGMVT1S, PGFLT1, &
& PEXTRA, PGMVTNDSL_DDH, PGFLTNDSL_DDH, PTRAJ_PHYS, PTRAJ_SLAG, YDDDH, CDPART, &
& PSAVEDP)

!**** *CPG* - Grid point calculations.

!     Purpose.
!     --------
!           Grid point calculations.

!**   Interface.
!     ----------
!        *CALL* *CPG(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        LDCONFX      : T if this call to CPG is done only to do some diagnostics like CFU,XFU
!                       (LDCONFX=T <==> former CDCONF='X')
!        LD_DFISTEP      : 'D' -> DFI computations
!        LDFSTEP      : .T. if first step.
!        LDDIAB       : .T. if complete physics is activated.
!        LDSLPHY      : .T. if ECMWF split physics.
!        LDUSEPB1     : .T. if updating PB1
!        PDT          : For a leap-frog scheme (three time level scheme):
!                       'dt' at the first time-step, '2 dt' otherwise.
!                       For a 2TL SL scheme: timestep 'dt'.
!        PDTPHY       : timestep used in the physics.
!        PTE          : 1. or 0. according to different configurations.
!        PBETADT      : BETADT or 0. according to different configurations.
!        YDSL         : SL_STRUCT definition
!        PGFLSLP      : GFL array for use in semi-lagrangian physics
!        PSAVTEND     : tendencies buffer for split ECMWF physics.
!        PGMVTNDHD_DDH: tendencies of horizontal diffusion scheme for GMV.
!        PGFLTNDHD_DDH: tendencies of horizontal diffusion for spectrally treated GFL.
!        PGPSDT2D     : cf. YGPSDT%GP2D in YOMSPSDT (buffer for stochastic physics).

!     INPUT/OUTPUT:
!     -------------
!        PGMV         : upper air GMV variables at t and t-dt.
!        PGMVS        : surface GMV variables at t and t-dt.
!        PGFL         : unified_treatment grid-point fields at t
!        PGFLPC       : unified_treatment grid-point fields at t (3TL PC only)
!        PGFLPT       : tendency of X variable from phy.
!        PSD_PF       : precip fraction
!        PSD_XP       : precipitation type diagnostic
!        PPGTERM      : pressure gradient quantities, only used if LPGFSF=T
!        PGMVTNDSI_DDH: tendencies of semi-implicit scheme.
!        PGMU0        : COSINE OF SOLAR ZENITH ANGLE, APPROXIMATE ACTUAL VALUE
!                       linear T_e correction
!                       linear T_e correction
!        YDDDH        : diagnostic superstructure

!     OUTPUT:
!     -------
!        PB1          : "SLBUF1" buffer for interpolations in SL scheme.
!        PGMVT1       : upper air GMV variables buffer.
!        PGMVT1S      : surface GMV variables buffer.
!        PGFLT1       : GFL variables buffer.
!        PEXTRA       : additional quantity for diagnostics.
!        PGMVTNDSL_DDH: GMV(t+dt,F)-GMV(t or t-dt,O) for DDH

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE documentation vol 2 ch 1 and vol 3 ch 6

!     Author.
!     -------
!      Philippe Courtier  *DMN*
!      Original : 90-03-19

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!    4-Mar-2008 Y. Seity : Cleaning IR and WV similated sat pictures
!                            (replaced by Fullpos Calculations)
!   K. Yessad (Sep 2008): prune enhanced diffusion (lfrein).
!   18-Nov-2008 R.Brozkova: OpenMP bugfix
!   K. Yessad (Dec 2008): remove dummy CDLOCK + cleanings
!   18-May-2009 S. Riette : mean cls wind added
!   F. Vana  15-Oct-2009: NSPLTHOI option
!   15-October-2009 Y. Bouteloup : Store radiative cloud water and ice in GFL (YIRAD and YLRAD)
!   K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP.
!   K. Yessad (Nov 2009): prune lpc_old.
!   24-Feb-2010 S. Riette : max gust avec NXGSTPERIOD seconds
!   11-Mai-2010 F. Bouyssel : RINDX, RINDY in argument
!   L. Bengtsson-Sedlar & F. Vana  18-Feb-2011 : CA arrays for MF physics
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures + cleanings.
!   F. Vana  21-Feb-2011: horiz. turbulence, N[x]LAG=4, diffus. on phys. tend.
!   G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!   G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!   2011-09-07, J.M. Piriou and E. Bazile: open/close output LFA files for 1D model MUSC (LMUSCLFA).
!   K. Yessad (Fev 2012): various modifications.
!   R. El Khatib 16-Mar-2012 Automatic allocations + cleanings
!   F. Bouttier  Jul 2012: pass stochastic physics pattern ZGPSDT2D
!   M. Ahlgrimm  31-Oct-2011 add rain, snow and PEXTRA to DDH output
!   T. Wilhelmsson 29-Mar-2012 Remove LCPG_SPLIT option and move OpenMP loop to CPG_DRV
!   F. Vana  28-Nov-2013 : Redesigned trajectory handling.
!   M. Diamantakis Dec 2013: add KSETTLOFF - may be used when LSETTLSVF=T
!   M. Ahlgrimm Apr 2014: Add lake variables and precip fraction to DDH output
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   2013-11, J. Masek: Passing intermittency arrays for ACRANEB2.
!   K. Yessad (July 2014): Move some variables.
!   O. Marsden (May 2016): Replace NPROMA9 by NPROMA for dimensions, as NPROMA9 not compatible with some routines called
!   F.Taillefer (June 2016): add MF assim CLS arrays in CPG_DIA
!   K. Yessad (Dec 2016): Prune obsolete options.
!   2017-09, J. Masek: Shifted dimensioning of PGMU0.
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   R. El Khatib 05-Jun-2018 computation of periods moved from cnt4 (OOPS refactoring)
!   M. Diamantakis (June 2018): Extra arguments for multiple (Atlas) grid advection scheme
!   R. Brozkova (Sept 2018) : Dataflow for global normal irradiance and mean
!                             radiant temperature.
!   F. Vana  Oct-2018: More optimal trigger of SLPHY computation.
!   I. Etchevers (Jan 2019) : PSD_XP for precipitation type diagnostic
!   F. Vana 11-Jul-2019: Option LRHS_CURV
!   2019-09, M. Hrastinski: Dataflow for TKE and TTE terms in ALARO DDH (ZFTCNS)
!   I. Polichtchouk (June 2020): code for LPGFSF=T
!   F. Vana 21-Sep-2020: LPGREUSE
!   D. Nemec (July 2021): Add ZQGRAUPEL.
!   F. Vana and S. Lang January 2023: new seq. physics order (ICI and SPPT)
! End Modifications
!-------------------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
!
USE FIELDS_MOD         , ONLY : FIELDS
!
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_DYN_TYPE, CPG_PHY_TYPE, CPG_MISC_TYPE, &
                              & CPG_TND_TYPE, CPG_TMP_TYPE, CPG_DDH_TYPE, &
                              & CPG_SL1_TYPE, CPG_SL2_TYPE, CPG_GPAR_TYPE
USE CPG_OPTS_TYPE_MOD  , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,ONLY : MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE INTDYN_MOD         , ONLY : TPG_TYPE ,YYTHW0   ,YYTCTY0  ,&
 &                              YYTRCP0  ,YYTRCP9  ,YYTXYB0_PHY,YYTXYB9_PHY,YYTXYB0  ,YYTXYB9
USE EINT_MOD           , ONLY : SL_STRUCT
USE YOMTRAJ            , ONLY : TRAJ_PHYS_TYPE, TRAJ_SLAG_TYPE
USE DDH_MIX            , ONLY : TYP_DDH, SETDDH, CLEANDDH
USE YOE_PHYS_MWAVE     , ONLY : N_PHYS_MWAVE
USE SPP_MOD            , ONLY : TSPP_CONFIG, TSPP_DATA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_BNDS_TYPE),INTENT(IN)   :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE),INTENT(IN)   :: YDCPG_OPTS
TYPE(CPG_TMP_TYPE),INTENT(INOUT) :: YDTMP
TYPE(CPG_TND_TYPE),INTENT(INOUT) :: YDCPG_TND
TYPE(CPG_SL1_TYPE),INTENT(INOUT) :: YDCPG_SL1
TYPE(CPG_SL2_TYPE),INTENT(INOUT) :: YDCPG_SL2
TYPE(CPG_GPAR_TYPE),INTENT(INOUT):: YDCPG_GPAR
TYPE(CPG_MISC_TYPE),INTENT(INOUT):: YDCPG_MISC
TYPE(CPG_DDH_TYPE),INTENT(INOUT) :: YDCPG_DDH
TYPE(CPG_PHY_TYPE),INTENT(INOUT) :: YDCPG_PHY0
TYPE(CPG_PHY_TYPE),INTENT(INOUT) :: YDCPG_PHY9
TYPE(MF_PHYS_TYPE),INTENT(INOUT) :: YDMF_PHYS
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN0
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN9
TYPE(MF_PHYS_SURF_TYPE),INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
!
TYPE(FIELDS)      ,INTENT(INOUT) :: YDFIELDS
!
INTEGER(KIND=JPIM),INTENT(OUT)   :: KSETTLOFF(YDCPG_OPTS%KFLEVG)
TYPE(SL_STRUCT),   INTENT(IN)    :: YDSL
REAL(KIND=JPRB),   INTENT(IN)    :: PGFLSLP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIMSLP)
REAL(KIND=JPRB),   INTENT(IN)    :: PSAVTEND(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_PHY_G%YRSLPHY%NVTEND)
REAL(KIND=JPRB),   INTENT(IN)    :: PGMVTNDHD_DDH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,2+YDMODEL%YRML_GCONF%YRDIMF%NFTHER)
REAL(KIND=JPRB),   INTENT(IN)    :: PGFLTNDHD_DDH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,1)
REAL(KIND=JPRB),   INTENT(IN)    :: PGPSDT2D(YDCPG_OPTS%KLON,YDMODEL%YRML_SPPT%YGPSDT(1)%NG2D)
TYPE(TPG_TYPE),    INTENT(INOUT) :: PPGTERM
REAL(KIND=JPRB),   INTENT(INOUT) :: PSD_PF(YDCPG_OPTS%KLON,YDFIELDS%YRSURF%YSD_PFD%NLEVS,YDFIELDS%YRSURF%YSD_PFD%NDIM)
REAL(KIND=JPRB),   INTENT(INOUT) :: PGMV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDFIELDS%YRGMV%NDIMGMV)
REAL(KIND=JPRB),   INTENT(INOUT) :: PGMVS(YDCPG_OPTS%KLON,YDFIELDS%YRGMV%NDIMGMVS)
REAL(KIND=JPRB),   INTENT(INOUT) :: PGFL(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB),   INTENT(INOUT) :: PGFLPC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIMPC)
REAL(KIND=JPRB),   INTENT(INOUT) :: PGFLPT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIMPT)
REAL(KIND=JPRB),   INTENT(INOUT) :: PSD_XA(YDCPG_OPTS%KLON,YDFIELDS%YRSURF%YSD_XAD%NLEVS,YDFIELDS%YRSURF%YSD_XAD%NDIM)
REAL(KIND=JPRB),   INTENT(INOUT) :: PGMVTNDSI_DDH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_DIAG%YRMDDH%NDIMSIGMV)
REAL(KIND=JPRB),   INTENT(INOUT) :: PHYS_MWAVE(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,N_PHYS_MWAVE)
REAL(KIND=JPRB),   INTENT(OUT)   :: PB1(YDSL%NASLB1,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB),   INTENT(OUT)   :: PGMVT1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDFIELDS%YRGMV%YT1%NDIM)
REAL(KIND=JPRB),   INTENT(OUT)   :: PGMVT1S(YDCPG_OPTS%KLON,YDFIELDS%YRGMV%YT1%NDIMS)
REAL(KIND=JPRB),   INTENT(OUT)   :: PGFLT1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1)
REAL(KIND=JPRB),   INTENT(OUT)   :: PEXTRA(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_PHY_G%YRDPHY%NVEXTRDYN)
REAL(KIND=JPRB),   INTENT(OUT)   :: PGMVTNDSL_DDH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,2+YDMODEL%YRML_GCONF%YRDIMF%NFTHER)
REAL(KIND=JPRB),   INTENT(OUT)   :: PGFLTNDSL_DDH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
TYPE (TRAJ_PHYS_TYPE),INTENT(INOUT) :: PTRAJ_PHYS
TYPE (TRAJ_SLAG_TYPE),INTENT(INOUT) :: PTRAJ_SLAG
TYPE(TYP_DDH),        INTENT(INOUT) :: YDDDH
CHARACTER(LEN=*),  INTENT(IN) :: CDPART
REAL(KIND=JPRB),   INTENT(INOUT) :: PSAVEDP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,3)

INTEGER(KIND=JPIM) :: JFLD
INTEGER(KIND=JPIM) :: JLEV

REAL(KIND=JPRB) :: ZUCOMP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)  ! work arrays for wind transformation

INTEGER(KIND=JPIM) :: ISLB1GFL9,ISLB1T9,ISLB1V9,ISLB1U9,ISLB1VD9

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "cpdysldia.intfb.h"
#include "cpg_dia.intfb.h"
#include "cpg_dyn.intfb.h"
#include "cpg_end.intfb.h"
#include "cpg_gp.intfb.h"
#include "ec_phys_lslphy.intfb.h"
#include "gpiniddh.intfb.h"
!     ------------------------------------------------------------------

IF (CDPART == 'XXX') THEN
  IF (LHOOK) CALL DR_HOOK('CPG', 0, ZHOOK_HANDLE)
ELSE
  IF (LHOOK) CALL DR_HOOK('CPG.'//CDPART, 0, ZHOOK_HANDLE)
ENDIF

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, YDGSGEOM=>YDGEOMETRY%YRGSGEOM(YDCPG_BNDS%KBL), &
 & YDDYNA=>YDMODEL%YRML_DYN%YRDYNA, YDDYN=>YDMODEL%YRML_DYN%YRDYN, &
 & YDPTRSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1, YDPTRSLB2=>YDMODEL%YRML_DYN%YRPTRSLB2, &
 & YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL, YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH, YDECUCONVCA=>YDMODEL%YRML_PHY_EC%YRECUCONVCA, &
 & YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY, YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, YGFL=>YDMODEL%YRML_GCONF%YGFL, &
 & YDSLPHY=>YDMODEL%YRML_PHY_G%YRSLPHY, YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY)

ASSOCIATE(NPROMA=>YDDIM%NPROMA, NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, NSITER=>YDDYN%NSITER, SLHDA=>YDDYN%SLHDA, SLHDD0=>YDDYN%SLHDD0, &
 & RCUCONVCA=>YDECUCONVCA%RCUCONVCA, RNLCONVCA=>YDECUCONVCA%RNLCONVCA, LCUCONV_CA=>YDECUCONVCA%LCUCONV_CA, &
 & GM=>YDGSGEOM%GM, LSDDH=>YDLDDH%LSDDH, LFLEXDIA=>YDLDDH%LFLEXDIA, &
 & LDDH_OMP=>YDLDDH%LDDH_OMP, MSLB1GFLP9=>YDPTRSLB1%MSLB1GFLP9, MSLB1SP9=>YDPTRSLB1%MSLB1SP9, MSLB1TP9=>YDPTRSLB1%MSLB1TP9, &
 & MSLB1UP9=>YDPTRSLB1%MSLB1UP9, MSLB1VP9=>YDPTRSLB1%MSLB1VP9, MSLB1ZP9=>YDPTRSLB1%MSLB1ZP9, NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & MSLB1SPPTCLEARP9=>YDPTRSLB1%MSLB1SPPTCLEARP9, MSLB2VVEL=>YDPTRSLB2%MSLB2VVEL, &
 & YSD_VF=>YDFIELDS%YRSURF%YSD_VF, YSP_RR=>YDFIELDS%YRSURF%YSP_RR, LSIMPH=>YDSIMPHL%LSIMPH, LMPHYS=>YDPHY%LMPHYS, &
 & LRHS_CURV=>YDDYN%LRHS_CURV, &
 & KIDIA=>YDCPG_BNDS%KIDIA, KFDIA=>YDCPG_BNDS%KFDIA, KSTGLO=>YDCPG_BNDS%KSTGLO, KBL=>YDCPG_BNDS%KBL)
!     ------------------------------------------------------------------

IF (CDPART (1:1) == 'X') THEN

  !*       3.    READ BUFFERS, COMPUTE AUXILIARY QUANTITIES.
  !              -------------------------------------------

  YDCPG_GPAR%ZVIEW(:,:)=0.0_JPRB

  CALL CPG_GP(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDTMP%CPG_GP, YDCPG_TND, YDCPG_MISC, YDCPG_GPAR,  &
  & YDCPG_DYN0, YDCPG_DYN9, YDMF_PHYS_SURF, YDVARS, YDMODEL, YDFIELDS,  PPGTERM,                    &
  & PGFL, PGMVTNDSL_DDH, PGFLTNDSL_DDH, YDCPG_SL2%ZVIEW, PGFLT1, YDCPG_MISC%KOZO, YDDDH)

  !     ------------------------------------------------------------------

  !*       4.    PHYSICS + DIAGNOSTICS
  !              ---------------------

  !*       4.2    Sets-up some quantities for DDH.
  !               A subset of these quantities may be modified by the MF-physics.

  IF (LSDDH) THEN
    IF (NCURRENT_ITER == 0) &
      CALL GPINIDDH(YDGEOMETRY, YDMDDH, KIDIA, KFDIA, YDCPG_DDH%DDHI, &
       & YDCPG_MISC%DHSF, YDCPG_DDH%DHCV, YDCPG_MISC%NEB, YDCPG_MISC%CLCT, KSTGLO)
    IF (LFLEXDIA.AND.LDDH_OMP) THEN
      CALL SETDDH(YDMDDH, YDDDH, KIDIA, KFDIA, NPROMA, NFLEVG, NCURRENT_ITER, YDCPG_MISC%DHSF, &
       & YDCPG_DDH%DDHI, YDCPG_DDH%AUX3D, YDCPG_DDH%AUX2D, YDCPG_DDH%AUXSM)
    ENDIF
  ENDIF

  !*       4.1    Fill ZSLBUF1AU.

  IF (YDCPG_OPTS%LUSEPB1) THEN
    DO JFLD=1,NFLDSLB1
      YDCPG_SL1%ZVIEW(1:KFDIA,JFLD)=0.0_JPRB
    ENDDO
  ENDIF

  !     ------------------------------------------------------------------
  ! -- need probably a test on LEPHYS here.

  !*       4.3    ECMWF-PHYSICS.

  !*       4.3.1   Call ECMWF unlagged physics.

  ! CALL EC_PHYS(...) to plug here.

  !*       4.3.2   Call ECMWF split physics.

  !     ------------------------------------------------------------------

  IF (YDCPG_OPTS%LSLPHY) THEN
    CALL EC_PHYS_LSLPHY(&
      & YDGEOMETRY, YDMODEL, KIDIA, KFDIA, YDCPG_OPTS%LFSTEP, YDCPG_OPTS%ZDTPHY, &
      & YDCPG_SL1%ZVIEW(:,MSLB1UP9), YDCPG_SL1%ZVIEW(:,MSLB1VP9), &
      & YDCPG_SL1%ZVIEW(:,MSLB1TP9), YDCPG_SL1%ZVIEW(:,MSLB1GFLP9), &
      & YDCPG_SL1%ZVIEW(:,MSLB1SPPTCLEARP9), PSAVTEND, PGFLSLP)

  IF (LRHS_CURV) THEN
    DO JLEV=1,NFLEVG
      ZUCOMP(KIDIA:KFDIA,JLEV) = YDCPG_SL1%ZVIEW(KIDIA:KFDIA,MSLB1UP9+JLEV-NFLSA)
      YDCPG_SL1%ZVIEW(KIDIA:KFDIA,MSLB1ZP9+JLEV-NFLSA) = &
        & YDCPG_SL1%ZVIEW(KIDIA:KFDIA,MSLB1VP9+JLEV-NFLSA)*YDGSGEOM%GSQM2(KIDIA:KFDIA)
      YDCPG_SL1%ZVIEW(KIDIA:KFDIA,MSLB1UP9+JLEV-NFLSA) = &
        & -YDCPG_SL1%ZVIEW(KIDIA:KFDIA,MSLB1UP9+JLEV-NFLSA)*YDGSGEOM%GESLO(KIDIA:KFDIA) &
        & -YDCPG_SL1%ZVIEW(KIDIA:KFDIA,MSLB1VP9+JLEV-NFLSA)*YDGSGEOM%GECLO(KIDIA:KFDIA)*YDGSGEOM%GEMU(KIDIA:KFDIA)
      YDCPG_SL1%ZVIEW(KIDIA:KFDIA,MSLB1VP9+JLEV-NFLSA) = &
        & ZUCOMP(KIDIA:KFDIA,JLEV)*YDGSGEOM%GECLO(KIDIA:KFDIA) &
        & -YDCPG_SL1%ZVIEW(KIDIA:KFDIA,MSLB1VP9+JLEV-NFLSA)*YDGSGEOM%GESLO(KIDIA:KFDIA)*YDGSGEOM%GEMU(KIDIA:KFDIA)
    ENDDO
  ENDIF
  ENDIF

  IF (LCUCONV_CA) THEN
  ! Quick fix. The arrays should rather be dimensionned with (nproma,ngpblks). REK
    YDMF_PHYS%OUT%CUCONVCA(1:KFDIA)=RCUCONVCA(KSTGLO:KSTGLO+KFDIA-1)
    YDMF_PHYS%OUT%NLCONVCA(1:KFDIA)=RNLCONVCA(KSTGLO:KSTGLO+KFDIA-1)
  ENDIF
  !     ------------------------------------------------------------------

ENDIF


IF (CDPART (2:2) == 'X') THEN
  !     ------------------------------------------------------------------

  !*       4.5    DIAGNOSTICS.
  !*              - for dynamics, uses outputs of CPG_GP (MF and ECMWF).
  !*              - for physics, uses outputs of MF_PHYS (MF only).

  IF (NCURRENT_ITER == 0) THEN
    CALL CPG_DIA(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDCPG_TND, YDCPG_MISC, YDCPG_GPAR, YDMF_PHYS, &
      & YDCPG_DYN0, YDMF_PHYS_SURF, YDVARS, YDFIELDS%YRGMV, YDFIELDS%YRSURF, YDFIELDS%YRCFU, YDFIELDS%YRXFU, &
      & YDMODEL, YDCPG_DDH%DDHI, PGMV, PGFL, YDCPG_SL2%ZVIEW(:, MSLB2VVEL), &
      & PSD_XA, PSD_PF, PGMVTNDSI_DDH, PGMVTNDHD_DDH, PGFLTNDHD_DDH, YDCPG_DDH%DHCV, YDDDH, YDCPG_MISC%FTCNS)
  ENDIF

  IF (YDDYNA%LSLDIA .AND. (NCURRENT_ITER == NSITER)) THEN
    ! * Calculate and print SL dynamics diagnostics.
    CALL CPDYSLDIA(YDGEOMETRY, YDDPHY, YDMODEL%YRML_GCONF, NPROMA, YDCPG_OPTS%LFSTEP, KIDIA, KFDIA, &
      & NFLEVG, KBL, YDVARS%SP%T0, YDCPG_DYN0%XYB%DELP, YDCPG_DYN0%PREF, YDVARS%T%T0, PGFL, PEXTRA)
  ENDIF

  !     ------------------------------------------------------------------

  !*       5.    DYNAMICS.
  !              ---------

  !*       5.1   Call dynamics.

  IF (YSD_VF%YLSM%LSET) THEN
    YDCPG_MISC%LSM(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)
  ELSE
    YDCPG_MISC%LSM(KIDIA:KFDIA)=0.0_JPRB
  ENDIF
  IF (YSP_RR%YT%LSET) THEN
    YDCPG_MISC%TSOL(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSP_RR%PT_T0(KIDIA:KFDIA)
  ELSE
    YDCPG_MISC%TSOL(KIDIA:KFDIA)=YDVARS%T%T0(KIDIA:KFDIA,NFLEVG)
  ENDIF

  CALL CPG_DYN(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDCPG_TND, YDCPG_MISC, YDCPG_DYN0, YDCPG_DYN9, &
    & YDVARS, YDFIELDS%YRGMV, YDMODEL, SLHDA(KSTGLO), SLHDD0(KSTGLO), &
    & PGFL, KSETTLOFF, PGFLPC, PGMV, PGMVS, YDCPG_SL1%ZVIEW, YDCPG_SL2%ZVIEW, PGMVT1, PGMVT1S, PGFLT1, &
    & PGMVTNDSI_DDH, PSAVEDP, PTRAJ_SLAG)

  IF (YDCPG_OPTS%LUSEPB1) THEN
    CALL CPG_PB1 (YDCPG_BNDS, YDCPG_OPTS, YDMODEL, YDSL, PB1, YDCPG_SL1%ZVIEW)
  ENDIF

  !     ------------------------------------------------------------------

  !*       6.    FINAL PART OF NON-LAGGED GRID-POINT COMPUTATIONS.
  !              -------------------------------------------------

  CALL CPG_END(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDFIELDS%YRGMV, YDFIELDS%YRSURF, YDMODEL%YRML_PHY_EC%YREPHY, &
    & YDMODEL%YRML_GCONF, YDDYN, YDDYNA, YDMODEL%YRML_PHY_MF, YDMODEL%YRML_LBC%LTENC, YDDPHY, GM(1), &
    & YDCPG_MISC%QS, YDMF_PHYS_SURF%GSP_RR%PT_T9, PGFL, YDCPG_GPAR%ZVIEW, PGMV, PGMVS, PTRAJ_PHYS)

  !     ------------------------------------------------------------------

  !     ------------------------------------------------------------------

  !*       7.    DDH CLEANING.
  !              -------------------------------------------------

  IF (LSDDH.AND.LFLEXDIA.AND.LDDH_OMP) THEN
    CALL CLEANDDH(YDMODEL%YRML_DIAG%YRTDDH, YDDDH, KSTGLO, NCURRENT_ITER)
  ENDIF

ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE

IF (CDPART == 'XXX') THEN
  IF (LHOOK) CALL DR_HOOK('CPG', 1, ZHOOK_HANDLE)
ELSE
  IF (LHOOK) CALL DR_HOOK('CPG.'//CDPART, 1, ZHOOK_HANDLE)
ENDIF

END SUBROUTINE CPG
