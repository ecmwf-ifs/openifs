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
SUBROUTINE CPG_DIA(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDCPG_TND, YDCPG_MISC, YDCPG_GPAR, YDMF_PHYS, &
& YDCPG_DYN0, YDMF_PHYS_SURF, YDVARS, YDGMV, YDSURF, YDCFU, YDXFU, YDMODEL, KDDHI, &
& PGMV, PGFL, POMEGA, PEXT, PCOVPTOT, PGMVTNDSI, PGMVTNDHD, PGFLTNDHD, PDHCV,  &
& YDDDH, PFTCNS)

!**** *CPG_DIA* - Grid point calculations: diagnostics.

!     Purpose.
!     --------
!           Grid point calculations: diagnostics using non lagged dynamics
!           and non lagged physics (DDH, CFU, XFU).

!**   Interface.
!     ----------
!        *CALL* *CPG_DIA(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        KST       : first element of work.
!        KEND      : last element of work.
!        KSTC,KENDC: the same as KST,KEND but including zone C for ALADIN.
!        KSTGLO    : global offset.
!        KDDHI     : internal domain to which each point belongs (DDH).
!        LDCONFX   : (see in CPG)
!        LDLDIAB   : .T. if complete physics is activated.
!        LDLFSTEP  : .T. if first step.
!        LD_DFISTEP: 'D' -> DFI computations
!        KIBL      : index into YRGSGEOM/YRCSGEOM types in YDGEOMETRY
!        POROGL    : zonal component of "grad(surf orography)"
!        POROGM    : meridian component of "grad(surf orography)"
!        PGMV      : upper air GMV variables at time t and t-dt.
!        PGFL      : unified_treatment grid-point fields at t.
!        PRE0      : hydrostatic pressure "prehyd" at half levels at time t.
!        PRE0F     : hydrostatic pressure "prehyd" at full levels at time t.
!        PXYB0     : contains pressure depth, "delta", "alpha" at t.
!        PHI0      : geopotential height "gz" at t at half levels.
!        PHIF0     : geopotential height "gz" at t at full levels.
!        PHI0FL    : zonal component of "grad(gz)" at full levels.
!        PHI0FM    : meridian component of "grad(gz)" at full levels.
!        PRCP0     : contains "cp", "R" and "Kap=R/Cp" at t.
!        PCTY0     : contains vertical velocities, vertical integral of divergence at t.
!        POMEGA    : "omega" at full levels at t, including the "lrubc" and "delta m=1" effects.
!        PKENE0    : kinetic energy at t.
!        PRT0L     : zonal component of "grad(RT)" at full levels.
!        PRT0M     : meridian component of "grad(RT)" at full levels.
!        PRE0L     : zonal component of "grad(prehyds)" at full levels.
!        PRE0M     : meridian component of "grad(prehyds)" at full levels.
!        PNHPRE0F  : "pre" at full levels (time t).
!        PNHPRE0H  : "pre" at half levels (time t).
!        PQCHA0L   : zonal comp grad(log(pre/prehyd)).
!        PQCHA0M   : merid comp grad(log(pre/prehyd)).
!        PDHSF     : distribution of horizontal mean weights used for simplified radiation scheme.
!        ---------------------- output of aplpar ------------------------------
!        PALB      : model surface shortwave albedo.
!        PCAPE     : CAPE.
!        PCTOP     : pressure of top of convective cloud (diagnostic).
!        PCLCC     : convective cloud cover (diagnostic).
!        PCLCH     : high cloud cover (diagnostic).
!        PCLCL     : low cloud cover (diagnostic).
!        PCLCM     : medium cloud cover (diagnostic).
!        PCLCT     : total cloud cover (diagnostic).
!        PCLPH     : height (in meters) of the PBL.
!        PVEIN     : ventilation index in the PBL.
!        PCT       : thermical coefficient of soil-vegetation middle.
!        PDIFCQ    : convective flux of specific humidity (not rain/snow).
!        PDIFCQN   : convective flux of solid water (not rain/snow).
!        PDIFCQL   : convective flux of liquid water (not rain/snow).
!        PDIFCS    : convective flux of enthalpy (not rain/snow).
!        PDIFTQ    : turbulent flux (inc. "q" negative) of specific humidity.
!        PDIFTQN   : turbulent flux (inc. "q" negative) of solid water.
!        PDIFTQL   : turbulent flux (inc. "q" negative) of liquid water.
!        PDIFTS    : turbulent flux of enthalpy (or dry static energy).
!        PFCCQL    : convective condensation flux for liquid water.
!        PFCCQN    : convective condensation flux for ice.
!        PFCHSP    : heat flux from surface to deep soil.
!        PFCLL     : latent heat flux over liquid water (or wet soil).
!        PFCLN     : latent heat flux over snow (or ice).
!        PFCQNNG   : pseudo-flux of ice to correct for "qi"<0.
!        PFCQLNG   : pseudo-flux of liquid water to correct for "ql"<0.
!        PFCQNG    : pseudo-flux of water to correct for Q<0.
!        PFCS      : sensible heat flux at surface level.
!        PFCSQL    : stratiform condensation flux for liquid water.
!        PFCSQN    : stratiform condensation flux for ice.
!        PFEVL     : water vapour flux over liquid water (or wet soil).
!        PFEVN     : water vapour flux over snow (or ice) and frozen soil.
!        PFEVV     : evapotranspiration flux.
!        PFGEL     : freezing flux of soil water.
!        PFGELS    : freezing flux of soil water at surface level.
!        PFLWSP    : water flux from surface to deep soil.
!        PFONTE    : water flux corresponding to surface snow melt.
!        PFPLCL    : convective precipitation as rain.
!        PFPLCN    : convective precipitation as snow.
!        PFPLCG    : convective precipitation as graupel.
!        PFPLCH    : convective precipitation as hail.
!        PFPLSL    : stratiform precipitation as rain.
!        PFPLSN    : stratiform precipitation as snow.
!        PFPLSG    : stratiform precipitation as graupel.
!        PFPLSH    : stratiform precipitation as hail.
!        PMRT      : mean radiant temperature.
!        PFRMH     : mesospheric enthalpy flux.
!        PFRSO     : shortwave radiative flux.
!        PFRSOC    : shortwave clear sky radiative flux.
!        PFRSODS   : surface downwards solar flux.
!        PFRSOLU   : downward lunar flux at surface.
!        PFRSGNI   : surface global normal irradiance.
!        PFRSDNI   : surface direct normal irradiance.
!        PFRSOPS   : surface parallel solar flux.
!        PFRSOPT   : top parallel solar flux.
!        PFRTH     : longwave radiative flux.
!        PFRTHC    : longwave clear sky radiative flux.
!        PFRTHDS   : surface downwards IR flux.
!        PFTR      : transpiration flux.
!        PGZ0      : g*roughness lenght (current).
!        PGZ0H     : current g*thermal roughness lenght (if KVCLIV >= 8).
!        PNEB      : fractional cloudiness for radiation.
!        PQICE     : specific humidity of solid water for radiation.
!        PQLI      : specific humidity of liquid water for radiation.
!        PQRAIN    : specific humidity of rain for radiation.
!        PQSNOW    : specific humidity of snow for radiation.
!        PQGRAUPEL : specific humidity of graupel for radiation.
!        PQS       : specific humidity at surface level.
!        PEXT      : extra diagnostics from physics.
!        PRH       : relative humidity.
!        PRHCLS    : relative humidity at 2 meters (diagnostic).
!        PRUISL    : run-off flux out the interception water-tank.
!        PRUISP    : run-off flux in soil.
!        PRUISS    : run-off flux at surface level.
!        PSTRCU    : convective flux of momentum "U".
!        PSTRCV    : convective flux of momentum "V".
!        PSTRDU    : gravity wave drag flux "U".
!        PSTRDV    : gravity wave drag flux "V".
!        PSTRTU    : turbulent flux of momentum "U".
!        PSTRTV    : turbulent flux of momentum "V".
!        PUCLS     : U-component of wind at 10 meters (diagnostic).
!        PVCLS     : V-component of wind at 10 meters (diagnostic).
!        PNUCLS    : U-component of neutral wind at 10 meters (diagnostic).
!        PNVCLS    : V-component of neutral wind at 10 meters (diagnostic).
!        PTCLS     : temperature at 2 meters (diagnostic).
!        PQCLS     : specific humidity at 2 meters (diagnostic).
!        PTPWCLS   : wet bulb temperature at 2 meters (diagnostic).
!        PUGST     : U-component of gusts (diagnostic).
!        PVGST     : V-component of gusts (diagnostic).
!        PVISICLD  : visibility due to cloud (ice and liquid) water
!        PVISIHYRO : visibility due to precipitation
!        ---------------------- end of output of aplpar -----------------------
!        PMOCON    : moisture convergence.
!        PFDIS     : enthalpy flux due to dissipation of kinetic energy. 
!        PFHPCL    : liquid water convective condensation enthalpy flux.
!        PFHPCN    : snow convective condensation enthalpy flux.
!        PFHPSL    : liquid water stratiform condensation enthalpy flux.
!        PFHPSN    : snow stratiform condensation enthalpy flux.
!        PFHSCL    : sensible heat flux due to liquid convective precipitations
!        PFHSCN    : sensible heat flux due to snow convective precipitations.
!        PFHSSL    : sensible heat flux due to liquid stratiform precipitations
!        PFHSSN    : sensible heat flux due to snow stratiform precipitations.
!        PFPFPSL   : flux of liquid resolved precipitation: the generation term.
!        PFPFPSN   : flux of snow resolved precipitation: the generation term.
!        PFPFPCL   : flux of liquid conv. precipitation: the generation term.
!        PFPFPCN   : flux of solid conv. precipitation: the generation term.
!        PFPEVPSL  : resolved precipitation flux due to evaporation.
!        PFPEVPSN  : resolved precipitation flux due to sublimation of snow.
!        PFPEVPCL  : convective precipitation flux due to evaporation.
!        PFPEVPCN  : convective precipitation flux due to sublimation of snow.
!        PFCQRNG   : pseudo-flux of rain to correct for "qr"<0
!        PFCQNNG   : pseudo-flux of snow to correct for "qs"<0
!        PGMVTNDSI : tendencies of semi-implicit scheme.
!        PGMVTNDHD : tendencies of horizontal diffusion scheme for GMV.
!        PGFLTNDHD : tendencies of horizontal diffusion for spectrally treated GFL.
!        PATND     : adiabatic Lagrangian tendencies.

!        PCOVPTOT  : precip fraction 

!     INPUT/OUTPUT:
!     -------------
!        PDHCV     : result array containing flux divergences and fluxes.
!                    Entering this routine this array is assumed to be
!                    pre-initialised with zeros.
!        PDRNSHF   : Derivative of non solar surface fluxes
!        YDDDH     : diagnostic superstructure

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
!   13-Aug-2001 K. YESSAD after part 8 of CPG. 

! Modifications
! -------------
!   04-Mar-2008 Y. Seity  Cleaning MTS pictures (IR and WV) replaced by Fullpos
!   01-Oct-2008 O.riviere Introduction of new data flow for DDH
!   K. Yessad (Dec 2008): remove dummy CDLOCK + cleanings
!   18-May-2009 S. Riette : mean cls wind added
!   14-Oct-2009 A. Alias   ZUCLS, ZVCLS added in CUMCPLDM (E.Maisonnave)
!   04-Jan-2011 Y. Seity Add PDIAGH (AROME hail daignostic)
!   K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   2011-10-04 J.M. Piriou: do not call CPPHDDH if flexible DDH.
!   K. Yessad (Nov 2011): various modifications.
!   M. Ahlgrimm  31-Oct-2011 add rain, snow and PEXTRA to DDH output
!   M. Ahlgrimm Apr 2014: Add lake variables and precip fraction to DDH output
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): Move some variables, rename some variables.
!   R. El Khatib 07-Mar-2016 Pruning of ISP
!   Y. Seity (4-Sept-2017) add LDLRESET_GST* to replace KXXGSTTS
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   R. El Khatib 05-Jun-2018 computation of periods moved from cnt4 (OOPS refactoring)
!   R. Brozkova (Jul 2018): Dataflow for global normal irradiance and mean
!                           radiant temperature.
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!   I. Etchevers (Sept 2018) : Add visisbilites
!   I. Etchevers (Janv 2019) : Add precipitation type
!   M. Hrastinski (Sep-2019): Dataflow for TKE and TTE terms in ALARO DDH (PFTCNS)
!   D. Nemec (July 2021) : Add prognostic variables for graupel (PFHPSG, PFHPCG,
!                          PFHSSG, PFPFPSG, PFPEVPSG, PFPEVPCG, PFCQGNG, PQGRAUPEL)
! End Modifications
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_DYN_TYPE, CPG_MISC_TYPE, CPG_TND_TYPE, CPG_GPAR_TYPE
USE CPG_OPTS_TYPE_MOD   , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,ONLY : MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMCFU             , ONLY : TCFU
USE YOMXFU             , ONLY : TXFU
USE YOMGMV             , ONLY : TGMV
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0             , ONLY : LAROME
USE YOMCT3             , ONLY : NSTEP
USE YOMINI             , ONLY : LINITER
USE YOMLUN             , ONLY : NULOUT
USE INTDYN_MOD         , ONLY : YYTCTY0, YYTRCP0, YYTXYB0
USE DDH_MIX            , ONLY : NTOTFIELD, RDDH_FIELD, RESET_DDHFLEX, RDDHSURF_FIELD, TYP_DDH

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_BNDS_TYPE),INTENT(IN)   :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE),INTENT(IN)   :: YDCPG_OPTS
TYPE(CPG_TND_TYPE),INTENT(INOUT) :: YDCPG_TND
TYPE(CPG_MISC_TYPE),INTENT(INOUT):: YDCPG_MISC
TYPE(CPG_GPAR_TYPE),INTENT(INOUT):: YDCPG_GPAR
TYPE(MF_PHYS_TYPE),INTENT(INOUT) :: YDMF_PHYS
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN0
TYPE(MF_PHYS_SURF_TYPE),INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
TYPE(TCFU)        ,INTENT(INOUT) :: YDCFU
TYPE(TXFU)        ,INTENT(INOUT) :: YDXFU
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KDDHI(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POMEGA(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_PHY_G%YRDPHY%NVEXTR)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMVTNDSI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,2+YDMODEL%YRML_GCONF%YRDIMF%NFTHER)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMVTNDHD(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,2+YDMODEL%YRML_GCONF%YRDIMF%NFTHER)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFLTNDHD(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDHCV(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_DIAG%YRMDDH%NDHCVSUN) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCOVPTOT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,1)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFTCNS(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,6)
TYPE(TYP_DDH)     ,INTENT(INOUT) :: YDDDH 
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: ILEVCO, JLEV,  JROF, JCV
INTEGER (KIND=JPIM) :: IPOS,JEXT

LOGICAL :: LLDYN, LLPHY
LOGICAL :: LLCALL_CUMCPLDM,LLCALL_CPCFU,LLCALL_CPXFU

REAL(KIND=JPRB) :: ZDHCS(YDGEOMETRY%YRDIM%NPROMM,YDMODEL%YRML_DIAG%YRMDDH%NDHCSSU)
REAL(KIND=JPRB) :: ZDUM1(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZENTRA(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZENTRV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZNEB(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZQSAT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZRH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSTRCU(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSTRCV(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSTRDU(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSTRDV(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSTRTU(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSTRTV(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZUCLS(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZNUCLS(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZUGST(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZUZGEO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZVCLS(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZNVCLS(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZVGST(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZDPRECIPS(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS%NDTPREC)
REAL(KIND=JPRB) :: ZDPRECIPS2(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS%NDTPREC2)
REAL(KIND=JPRB) :: ZVMGEO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZFPLCL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZFPLCN(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZFPLCG(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZFPLSL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZFPLSN(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZFPLSG(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZQICE(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZQLI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZQRAIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZQSNOW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZQGRAUPEL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZEXT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_PHY_G%YRDPHY%NVEXTR)
REAL(KIND=JPRB) :: ZQS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZTSOL(YDGEOMETRY%YRDIM%NPROMA)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!#include "posinsub.intfb.h"

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "cpcfu.intfb.h"
#include "cpcuddh.intfb.h"
#include "cpcuddh_omp.intfb.h"
#include "cpdyddh.intfb.h"
#include "cpphddh.intfb.h"
#include "cpxfu.intfb.h"
#include "cumcpl.intfb.h"
#include "gprh.intfb.h"
#include "cpcls_assim.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG_DIA', 0, ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDVAB=>YDGEOMETRY%YRVAB, YDGSGEOM=>YDGEOMETRY%YRGSGEOM(YDCPG_BNDS%KBL), &
  & YDCVER=>YDGEOMETRY%YRCVER, YDCSGEOM=>YDGEOMETRY%YRCSGEOM(YDCPG_BNDS%KBL), &
  & YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY, YDSDDH=>YDMODEL%YRML_DIAG%YRSDDH, &
  & YDMCC=>YDMODEL%YRML_AOC%YRMCC, YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH, &
  & YDRIP=>YDMODEL%YRML_GCONF%YRRIP, YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY, YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL, &
  & YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, YGFL=>YDMODEL%YRML_GCONF%YGFL, YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
  & YDDIMF=>YDMODEL%YRML_GCONF%YRDIMF, YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR)

ASSOCIATE(NDIM=>YGFL%NDIM, YA=>YGFL%YA, YI=>YGFL%YI, YL=>YGFL%YL, YO3=>YGFL%YO3, YQ=>YGFL%YQ, YR=>YGFL%YR, &
  & YS=>YGFL%YS, YG=>YGFL%YG, LCUMFU=>YDCFU%LCUMFU, NPROMA=>YDDIM%NPROMA, NPROMM=>YDDIM%NPROMM, NPROMNH=>YDDIM%NPROMNH, &
  & NFTHER=>YDDIMF%NFTHER, NFLEVG=>YDDIMV%NFLEVG, NTSSG=>YDDPHY%NTSSG, NVEXTR=>YDDPHY%NVEXTR, NCEXTR=>YDDPHY%NCEXTR, &
  & LAGPHY=>YDEPHY%LAGPHY, LEPHYS=>YDEPHY%LEPHYS, NDIMGMV=>YDGMV%NDIMGMV, YT0=>YDGMV%YT0, LFLEXDIA=>YDLDDH%LFLEXDIA, &
  & LSDDH=>YDLDDH%LSDDH, LDDH_OMP=>YDLDDH%LDDH_OMP, NFRCPL=>YDMCC%NFRCPL, NDHCSSU=>YDMDDH%NDHCSSU, NDHCVSUN=>YDMDDH%NDHCVSUN, &
  & NSTOP=>YDRIP%NSTOP, LHDPAS=>YDSDDH%LHDPAS, NDTPREC=>YDPHY%YRDPRECIPS%NDTPREC, NDTPREC2=>YDPHY%YRDPRECIPS%NDTPREC2, &
  & YSD_VF=>YDSURF%YSD_VF, YSD_VFD=>YDSURF%YSD_VFD, YSP_RR=>YDSURF%YSP_RR, YSP_RRD=>YDSURF%YSP_RRD, YSP_SB=>YDSURF%YSP_SB, &
  & YSP_SBD=>YDSURF%YSP_SBD, YSP_SG=>YDSURF%YSP_SG, YSP_SGD=>YDSURF%YSP_SGD, YSD_XP=>YDSURF%YSD_XP, YSD_XP2=>YDSURF%YSD_XP2, &
  & LDPRECIPS=>YDPHY%LDPRECIPS, LDPRECIPS2=>YDPHY%LDPRECIPS2, LXFU=>YDXFU%LXFU, LDIRCLSMOD=>YDDPHY%LDIRCLSMOD, &
  & LTRAJPS=>YDSIMPHL%LTRAJPS, LSIMPH=>YDSIMPHL%LSIMPH, LMPHYS=>YDPHY%LMPHYS, LMPA=>YDARPHY%LMPA, MVQS=>YDPARAR%MVQS, &
  & MVTS=>YDPARAR%MVTS, LDIRSICMOD=>YDDPHY%LDIRSICMOD, &
  & KIDIA=>YDCPG_BNDS%KIDIA, KFDIA=>YDCPG_BNDS%KFDIA, KSTGLO=>YDCPG_BNDS%KSTGLO)
      
!     ------------------------------------------------------------------
IPOS=1
!CALL POSINSUB('CPG_DIA',IPOS)
!*       1.    DDH:
!              ----


IF ((NSTEP >= 0).AND.LSDDH.AND.(.NOT.LINITER).AND.(.NOT.YDCPG_OPTS%LCONFX)) THEN 

  IF (LMPA) THEN
    ZQICE(1:NPROMA,1:NFLEVG) = YDVARS%I%T0(1:NPROMA,1:NFLEVG)
    ZQLI(1:NPROMA,1:NFLEVG) = YDVARS%L%T0(1:NPROMA,1:NFLEVG)
    ZQRAIN(1:NPROMA,1:NFLEVG) = YDVARS%R%T0(1:NPROMA,1:NFLEVG)
    ZQSNOW(1:NPROMA,1:NFLEVG) = YDVARS%S%T0(1:NPROMA,1:NFLEVG)
    ZQGRAUPEL(1:NPROMA,1:NFLEVG) = YDVARS%G%T0(1:NPROMA,1:NFLEVG)
  ELSE
    IF(LMPHYS.OR.(LEPHYS.AND..NOT.LAGPHY).OR.LSIMPH) THEN
      ZFPLCL(1:NPROMA,0:NFLEVG) = YDMF_PHYS%OUT%FPLCL(1:NPROMA,0:NFLEVG)
      ZFPLCN(1:NPROMA,0:NFLEVG) = YDMF_PHYS%OUT%FPLCN(1:NPROMA,0:NFLEVG)
      ZFPLSL(1:NPROMA,0:NFLEVG) = YDMF_PHYS%OUT%FPLSL(1:NPROMA,0:NFLEVG)
      ZFPLSN(1:NPROMA,0:NFLEVG) = YDMF_PHYS%OUT%FPLSN(1:NPROMA,0:NFLEVG)
      ZFPLSG(1:NPROMA,0:NFLEVG) = YDMF_PHYS%OUT%FPLSG(1:NPROMA,0:NFLEVG)
      ZFPLCG(1:NPROMA,0:NFLEVG) = YDMF_PHYS%OUT%FPLCG(1:NPROMA,0:NFLEVG)
    ENDIF
    ZQICE(1:NPROMA,1:NFLEVG) = YDCPG_MISC%QICE(1:NPROMA,1:NFLEVG)
    ZQLI(1:NPROMA,1:NFLEVG) = YDCPG_MISC%QLI(1:NPROMA,1:NFLEVG)
    ZQRAIN(1:NPROMA,1:NFLEVG) = YDCPG_MISC%QRAIN(1:NPROMA,1:NFLEVG)
    ZQSNOW(1:NPROMA,1:NFLEVG) = YDCPG_MISC%QSNOW(1:NPROMA,1:NFLEVG)
    ZQGRAUPEL(1:NPROMA,1:NFLEVG) = YDCPG_MISC%QGRAUPEL(1:NPROMA,1:NFLEVG)
  ENDIF

  IF (LHDPAS) THEN
     IF (NVEXTR > 0) THEN
        IF( NCEXTR /= NFLEVG ) THEN
           CALL ABOR1(' CPG_DIA: NCEXTR /= NFLEVG')
        ENDIF
        DO JEXT=1,NVEXTR
           ZEXT(1:NPROMA,1:NFLEVG,JEXT)=PEXT(1:NPROMA,1:NFLEVG,JEXT)
        ENDDO
     ELSE
        ZEXT(1:NPROMA,1:NFLEVG,:)=0.0_JPRB
     ENDIF
  ENDIF 


  IF(LAGPHY.OR.(.NOT.YDCPG_OPTS%LDIAB)) THEN

    LLPHY=.FALSE.
    CALL GPRH(.FALSE., NPROMA, KIDIA, KFDIA, NFLEVG, 1.0_JPRB, 0.0_JPRB, YDVARS%Q%T0, &
      & YDVARS%T%T0, YDCPG_DYN0%PREF, ZQSAT, ZRH)
  
    ZDUM1(KIDIA:KFDIA)=0.0_JPRB
  
    IF(YA%LACTIVE) THEN
      ZNEB(KIDIA:KFDIA,1:NFLEVG)=YDVARS%A%T0(KIDIA:KFDIA,1:NFLEVG)
    ELSE
      IF (LMPA) THEN
        ZNEB(KIDIA:KFDIA,1:NFLEVG)=0.0_JPRB    ! ???
      ELSE
        ZNEB(KIDIA:KFDIA,1:NFLEVG)=YDCPG_MISC%NEB(KIDIA:KFDIA,1:NFLEVG)
      ENDIF
    ENDIF

  
    ! * COMP. OF DYNAMICAL VARIABLES AND FLUX/TENDENCIES
  
    CALL CPDYDDH(YDVAB, YDCVER, YDGMV, YDMODEL%YRML_GCONF, YDDPHY, YDMODEL%YRML_DIAG, YDPHY, YDMODEL%YRML_DYN%YRDYNA, &
      & NPROMA, NPROMNH, KIDIA, KFDIA, NFLEVG, YDGSGEOM, YDCSGEOM, YDCPG_DYN0%XYB%ZVIEW, YDCPG_DYN0%RCP%ZVIEW, &
      & YDCPG_DYN0%PHI, YDCPG_DYN0%PHIF, YDCPG_DYN0%PHIFL, YDCPG_DYN0%PHIFM, YDCPG_DYN0%KENE, YDCPG_DYN0%RTL, &
      & YDCPG_DYN0%RTM, YDCPG_DYN0%OROGL, YDCPG_DYN0%OROGM, YDCPG_DYN0%PREL, YDCPG_DYN0%PREM, YDCPG_DYN0%PREF, &
      & YDCPG_DYN0%PRE, YDCPG_DYN0%NHPREF, YDCPG_DYN0%NHPREH, YDCPG_DYN0%QCHAL, YDCPG_DYN0%QCHAM, YDCPG_DYN0%CTY%ZVIEW, &
      & POMEGA, ZFPLCL, ZFPLCN, ZFPLCG, ZFPLSL, ZFPLSN, ZFPLSG, ZDUM1, ZDUM1, ZRH, ZNEB, ZQLI, ZQICE, ZQRAIN, ZQSNOW, ZQGRAUPEL, &
      & ZEXT, PCOVPTOT, PGFL, PGMV, PGMVTNDSI, PGMVTNDHD, PGFLTNDHD, YDCPG_TND%ZVIEW, ZUZGEO, ZVMGEO, &
      & PDHCV, ZENTRA, ZENTRV, YDDDH)

  ELSE

    IF(.NOT.(LMPHYS.OR.LEPHYS)) THEN
      LLPHY=.FALSE.
      CALL GPRH(.FALSE., NPROMA, KIDIA, KFDIA, NFLEVG, 1.0_JPRB, 0.0_JPRB, YDVARS%Q%T0, &
        & YDVARS%T%T0, YDCPG_DYN0%PREF, ZQSAT, ZRH)
    ELSE
      IF (LMPA) THEN
        CALL GPRH(.FALSE., NPROMA, KIDIA, KFDIA, NFLEVG, 1.0_JPRB, 0.0_JPRB, YDVARS%Q%T0, &
          & YDVARS%T%T0, YDCPG_DYN0%PREF, ZQSAT, ZRH)
      ELSE
        ZRH(KIDIA:KFDIA,1:NFLEVG)=YDCPG_MISC%RH(KIDIA:KFDIA,1:NFLEVG)
      ENDIF
      LLPHY=.TRUE.
    ENDIF

    IF (LMPA) THEN
      ZQS(1:NPROMA) =  YDCPG_GPAR%ZVIEW(1:NPROMA,MVQS)
      ZTSOL(1:NPROMA) = YDCPG_GPAR%ZVIEW(1:NPROMA,MVTS)
    ELSE 
      ZQS(1:NPROMA) =  YDCPG_MISC%QS(1:NPROMA)
      ZTSOL(1:NPROMA) = YDMF_PHYS_SURF%GSP_RR%PT_T0(1:NPROMA)
    ENDIF

    ! * COMP. OF DYNAMICAL VARIABLES AND FLUX/TENDENCIES

!CALL POSINSUB('CPG_DIA',IPOS)
    CALL CPDYDDH(YDVAB, YDCVER, YDGMV, YDMODEL%YRML_GCONF, YDDPHY, YDMODEL%YRML_DIAG, YDPHY, YDMODEL%YRML_DYN%YRDYNA, &
      & NPROMA, NPROMNH, KIDIA, KFDIA, NFLEVG, YDGSGEOM, YDCSGEOM, YDCPG_DYN0%XYB%ZVIEW, YDCPG_DYN0%RCP%ZVIEW, &
      & YDCPG_DYN0%PHI, YDCPG_DYN0%PHIF, YDCPG_DYN0%PHIFL, YDCPG_DYN0%PHIFM, YDCPG_DYN0%KENE, YDCPG_DYN0%RTL, &
      & YDCPG_DYN0%RTM, YDCPG_DYN0%OROGL, YDCPG_DYN0%OROGM, YDCPG_DYN0%PREL, YDCPG_DYN0%PREM, YDCPG_DYN0%PREF, &
      & YDCPG_DYN0%PRE, YDCPG_DYN0%NHPREF, YDCPG_DYN0%NHPREH, YDCPG_DYN0%QCHAL, YDCPG_DYN0%QCHAM, YDCPG_DYN0%CTY%ZVIEW, &
      & POMEGA, ZFPLCL, ZFPLCN, ZFPLCG, ZFPLSL, ZFPLSN, ZFPLSG, ZQS, ZTSOL, ZRH, YDCPG_MISC%NEB, ZQLI, ZQICE, ZQRAIN, ZQGRAUPEL, &
      & ZQSNOW, ZEXT, PCOVPTOT, PGFL, PGMV, PGMVTNDSI, PGMVTNDHD, PGFLTNDHD, YDCPG_TND%ZVIEW, ZUZGEO, &
      & ZVMGEO, PDHCV, ZENTRA, ZENTRV, YDDDH)
  
!CALL POSINSUB('CPG_DIA',IPOS)
    ! * COMP. OF PHYS. FIELDS AND SURFACE FLUXES
  
    IF((LMPHYS.OR.LEPHYS).AND..NOT.LMPA.AND..NOT.LFLEXDIA) THEN
      CALL CPPHDDH (YDMODEL%YRML_DIAG, YDEPHY, YDRIP, YDMODEL%YRML_PHY_MF, NPROMA, KIDIA, KFDIA, NFLEVG, .TRUE., &
      & YDGSGEOM, ZUZGEO, ZVMGEO, YDMF_PHYS_SURF%GSP_SG%PF_T0, YDMF_PHYS_SURF%GSP_SB%PT_T0, &
      & YDMF_PHYS_SURF%GSP_RR%PT_T0, YDMF_PHYS_SURF%GSP_SB%PQ_T0, YDMF_PHYS_SURF%GSP_RR%PW_T0, YDVARS%T%T0, &
      & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS, &
      & YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLCG, &
      & YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLSG, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH, &
      & YDMF_PHYS%OUT%FRTHDS, YDMF_PHYS%OUT%FHPCL, YDMF_PHYS%OUT%FHPCN, YDMF_PHYS%OUT%FHPCG, YDMF_PHYS%OUT%FHPSL, &
      & YDMF_PHYS%OUT%FHPSN, YDMF_PHYS%OUT%FHPSG, YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, YDMF_PHYS%OUT%STRDU, &
      & YDMF_PHYS%OUT%STRDV, YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV, YDMF_PHYS%OUT%CT, YDMF_PHYS%OUT%FCHSP, &
      & YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FCS, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, &
      & YDMF_PHYS%OUT%FLWSP, YDMF_PHYS%OUT%FONTE, YDMF_PHYS%OUT%RUISP, YDMF_PHYS%OUT%RUISS, ZENTRA, &
      & ZENTRV, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCQN, &
      & YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FHSCL, YDMF_PHYS%OUT%FHSCN, YDMF_PHYS%OUT%FHSSL, &
      & YDMF_PHYS%OUT%FHSSN, YDMF_PHYS%OUT%FHSSG, YDMF_PHYS%OUT%FEPFP, YDMF_PHYS%OUT%FCMPCQ, YDMF_PHYS%OUT%FCMPSN, &
      & YDMF_PHYS%OUT%FCMPSL, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, &
      & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, &
      & YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, &
      & YDMF_PHYS%OUT%FPEVPCG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQGNG, YDCPG_DYN0%PHI, &
      & YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%QCLS, YDMF_PHYS%OUT%CLPH, &
      & YDMF_PHYS%OUT%ALB, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, PDHCV, ZDHCS, PFTCNS)
    ENDIF

    IF(LMPHYS.AND.LMPA) THEN

      ! Surface variables and fluxes set to zero.
      ! There are no surface variables for AROME in DDH, yet.
      IF(LLPHY) THEN
        DO JROF = 1, NPROMM
          DO JCV = 1, NDHCSSU
            ZDHCS(JROF,JCV)=0_JPRB
          ENDDO
        ENDDO
      ENDIF

    ENDIF

!CALL POSINSUB('CPG_DIA',IPOS)
  ENDIF
  
  ! * STORE/ACCUM. IN HORIZONTAL MEAN ARRAYS
  
  LLDYN=.TRUE.

  IF (LFLEXDIA) THEN
    IF (LDDH_OMP) THEN
      CALL CPCUDDH_OMP(YDMDDH, YDMODEL%YRML_DIAG%YRTDDH, NPROMA, KIDIA, KFDIA, NFLEVG, KDDHI, YDDDH, LLDYN)
    ELSE
      CALL CPCUDDH(YDMODEL%YRML_DIAG, YDRIP, NPROMA, KIDIA, KFDIA, NFLEVG, &
        & KDDHI, YDCPG_MISC%DHSF, RDDH_FIELD, RDDHSURF_FIELD, LLDYN, LLPHY, NTOTFIELD)
    ENDIF 
  ELSE 
    CALL CPCUDDH(YDMODEL%YRML_DIAG, YDRIP, NPROMA, KIDIA, KFDIA, NFLEVG, KDDHI, &
      & YDCPG_MISC%DHSF, PDHCV, ZDHCS, LLDYN, LLPHY, NDHCVSUN)
  ENDIF 
  
ENDIF

!     ------------------------------------------------------------------

!*       2.    CFU AND XFU:
!              ------------

!*    2.1   DEFINE THE CONDITIONS OF CALLING SOME ROUTINES.

LLCALL_CPCFU=&
 & LCUMFU.AND.(.NOT.LINITER).AND.(.NOT.YDCPG_OPTS%LCONFX).AND.(.NOT.LAGPHY)
LLCALL_CPXFU=LXFU.AND.YDCPG_OPTS%LDIAB.AND.(.NOT.LINITER).AND..NOT.YDCPG_OPTS%L_DFISTEP
LLCALL_CUMCPLDM= (&
 & NFRCPL /= 0).AND.(NFRCPL <= NSTOP).AND.(.NOT.YDCPG_OPTS%LCONFX).AND.(.NOT.LAGPHY).AND.(.NOT.LAROME)
! Check that if NPROMM /= NPROMA, all these LLCALL... keys are .F.:
IF ((LLCALL_CPCFU.OR.LLCALL_CPXFU.OR.LLCALL_CUMCPLDM)&
 & .AND.(NPROMM /= NPROMA)) THEN
  WRITE(NULOUT,*) ' CPG_DIA part 2: if NPROMM /= NPROMA,'
  WRITE(NULOUT,*) '  all the LLCALL_... keys must be equal to .F.'
  WRITE(NULOUT,*) '  Check the calculation of these keys.'
  WRITE(NULOUT,*) '  In particular, in adiabatic runs, CFU, XFU'
  WRITE(NULOUT,*) '  and cumul of coupled fields MUST NOT be called.'
  CALL ABOR1(' CPG_DIA: ABOR1 CALLED')
ENDIF

!*    2.2   DIVIDE SOME QUANTITIES LINKED TO (U,V) BY THE MAPPING FACTOR.

! * ZSTR.. arrays are used for CFU, XFU and in CUMCPL.
IF (LLCALL_CPCFU.OR.LLCALL_CPXFU.OR.LLCALL_CUMCPLDM) THEN
  DO JLEV=0,NFLEVG
    ZSTRCU(KIDIA:KFDIA,JLEV)=YDMF_PHYS%OUT%STRCU(KIDIA:KFDIA,JLEV)
    ZSTRCV(KIDIA:KFDIA,JLEV)=YDMF_PHYS%OUT%STRCV(KIDIA:KFDIA,JLEV)
    ZSTRDU(KIDIA:KFDIA,JLEV)=YDMF_PHYS%OUT%STRDU(KIDIA:KFDIA,JLEV)
    ZSTRDV(KIDIA:KFDIA,JLEV)=YDMF_PHYS%OUT%STRDV(KIDIA:KFDIA,JLEV)
    ZSTRTU(KIDIA:KFDIA,JLEV)=YDMF_PHYS%OUT%STRTU(KIDIA:KFDIA,JLEV)
    ZSTRTV(KIDIA:KFDIA,JLEV)=YDMF_PHYS%OUT%STRTV(KIDIA:KFDIA,JLEV)
  ENDDO
ENDIF

! * Z(U,V)CLS and Z(U,V)GST are used for XFU.
IF (LLCALL_CPXFU .OR. LDIRCLSMOD) THEN
  ZUCLS(KIDIA:KFDIA)=YDMF_PHYS%OUT%UCLS(KIDIA:KFDIA)
  ZVCLS(KIDIA:KFDIA)=YDMF_PHYS%OUT%VCLS(KIDIA:KFDIA)
  ZNUCLS(KIDIA:KFDIA)=YDMF_PHYS%OUT%NUCLS(KIDIA:KFDIA)
  ZNVCLS(KIDIA:KFDIA)=YDMF_PHYS%OUT%NVCLS(KIDIA:KFDIA)
  ZUGST(KIDIA:KFDIA)=YDMF_PHYS%OUT%UGST(KIDIA:KFDIA)
  ZVGST(KIDIA:KFDIA)=YDMF_PHYS%OUT%VGST(KIDIA:KFDIA)
ENDIF

! * Divide the arrays by the mapping factor for the cases
!   where the quantities to diagnose are the reduced ones.
!!! -- Test to be checked, seems suspicious -- !!!
IF (YDCPG_OPTS%LDIAB.OR.(LSIMPH.AND.(.NOT.LTRAJPS))) THEN  
  IF (LLCALL_CPCFU.OR.LLCALL_CPXFU.OR.LLCALL_CUMCPLDM) THEN
    DO JLEV=0,NFLEVG
      ZSTRCU(KIDIA:KFDIA,JLEV)=ZSTRCU(KIDIA:KFDIA,JLEV)/YDGSGEOM%GM(KIDIA:KFDIA)
      ZSTRCV(KIDIA:KFDIA,JLEV)=ZSTRCV(KIDIA:KFDIA,JLEV)/YDGSGEOM%GM(KIDIA:KFDIA)
      ZSTRDU(KIDIA:KFDIA,JLEV)=ZSTRDU(KIDIA:KFDIA,JLEV)/YDGSGEOM%GM(KIDIA:KFDIA)
      ZSTRDV(KIDIA:KFDIA,JLEV)=ZSTRDV(KIDIA:KFDIA,JLEV)/YDGSGEOM%GM(KIDIA:KFDIA)
      ZSTRTU(KIDIA:KFDIA,JLEV)=ZSTRTU(KIDIA:KFDIA,JLEV)/YDGSGEOM%GM(KIDIA:KFDIA)
      ZSTRTV(KIDIA:KFDIA,JLEV)=ZSTRTV(KIDIA:KFDIA,JLEV)/YDGSGEOM%GM(KIDIA:KFDIA)
    ENDDO
  ENDIF
  IF (LLCALL_CPXFU .OR. LDIRCLSMOD) THEN
    ZUCLS(KIDIA:KFDIA)=ZUCLS(KIDIA:KFDIA)/YDGSGEOM%GM(KIDIA:KFDIA)
    ZVCLS(KIDIA:KFDIA)=ZVCLS(KIDIA:KFDIA)/YDGSGEOM%GM(KIDIA:KFDIA)
    ZNUCLS(KIDIA:KFDIA)=ZNUCLS(KIDIA:KFDIA)/YDGSGEOM%GM(KIDIA:KFDIA)
    ZNVCLS(KIDIA:KFDIA)=ZNVCLS(KIDIA:KFDIA)/YDGSGEOM%GM(KIDIA:KFDIA)
    ZUGST(KIDIA:KFDIA)=ZUGST(KIDIA:KFDIA)/YDGSGEOM%GM(KIDIA:KFDIA)
    ZVGST(KIDIA:KFDIA)=ZVGST(KIDIA:KFDIA)/YDGSGEOM%GM(KIDIA:KFDIA)
  ENDIF
ENDIF

!*    2.3   CUMUL OF FLUX DIAGNOSTICS

IF(LLCALL_CPCFU)THEN  
  CALL CPCFU(YDGEOMETRY, YDCFU, YDRIP, KIDIA, KFDIA, YDCPG_BNDS%KIDIE, YDCPG_BNDS%KFDIE, KSTGLO, YDCPG_OPTS%LFSTEP, &
    & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS, &
    & YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLCL, &
    & YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLCG, YDMF_PHYS%OUT%FPLCH, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, &
    & YDMF_PHYS%OUT%FPLSG, YDMF_PHYS%OUT%FPLSH, YDMF_PHYS%OUT%FRMH, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRSOC, &
    & YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%FRTHC, YDMF_PHYS%OUT%FDIS(:, NFLEVG), ZSTRCU, ZSTRCV, ZSTRDU,&
    & ZSTRDV, ZSTRTU, ZSTRTV, YDCPG_DYN0%XYB%DELP, YDCPG_MISC%NEB, YDCPG_MISC%QICE, YDCPG_MISC%QLI, YDCPG_DYN0%PRE(1:,NFLEVG), &
    & YDMF_PHYS_SURF%GSP_SG%PF_T0, YDVARS%O3%T0, YDVARS%Q%T0, YDMF_PHYS_SURF%GSP_RR%PW_T0, YDMF_PHYS%OUT%FEVL, &
    & YDMF_PHYS%OUT%FEVN, YDCPG_MISC%CLCT, YDMF_PHYS%OUT%CLCH, YDMF_PHYS%OUT%CLCM, YDMF_PHYS%OUT%CLCL, &
    & YDMF_PHYS%OUT%CLCC, YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FCS, YDMF_PHYS%OUT%FONTE, &
    & YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FLWSP, YDMF_PHYS%OUT%RUISS, YDMF_PHYS%OUT%RUISP, YDMF_PHYS%OUT%FEVV, &
    & YDMF_PHYS%OUT%FTR, YDMF_PHYS%OUT%RUISL, YDMF_PHYS%OUT%FGEL, YDMF_PHYS%OUT%FGELS, YDMF_PHYS%OUT%FRSODS, &
    & YDMF_PHYS%OUT%FRSOPS, YDMF_PHYS%OUT%FRSOPT, YDMF_PHYS%OUT%FRSOLU, YDMF_PHYS%OUT%FRTHDS, YDMF_PHYS%OUT%FRSDNI, &
    & YDMF_PHYS%OUT%FRSGNI, YDMF_PHYS%OUT%FLASH)
ENDIF

IF(LLCALL_CUMCPLDM) THEN  
  ILEVCO=NFLEVG
  CALL CUMCPL(YDGEOMETRY%YRDIM, YDRIP, YDPHY, ILEVCO, KIDIA, KFDIA, KSTGLO, NTSSG, ZSTRTU(1, NFLEVG), ZSTRTV(1, NFLEVG), &
    & YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FCS, &
    & YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, YDMF_PHYS%OUT%FPLCL(:,NFLEVG), YDMF_PHYS%OUT%FPLCN(:,NFLEVG), &
    & YDMF_PHYS%OUT%FPLSL(:,NFLEVG), YDMF_PHYS%OUT%FPLSN(:, NFLEVG), YDMF_PHYS%OUT%RUISS, YDMF_PHYS%OUT%RUISP, &
    & YDMF_PHYS_SURF%GSP_RR%PT_T0, YDMF_PHYS%OUT%DRNSHF, YDMF_PHYS%OUT%ALB, ZUCLS, ZVCLS)
ENDIF

!*    2.5   INTERFACE FOR INSTANTANEOUS FLUXES

IF (LLCALL_CPXFU) THEN
  IF (LDPRECIPS) THEN
    ZDPRECIPS(KIDIA:KFDIA,1:NDTPREC)=YDMF_PHYS_SURF%GSD_XP%PPRECIP(KIDIA:KFDIA,1:NDTPREC)
  ENDIF
  IF (LDPRECIPS2) THEN
    ZDPRECIPS2(KIDIA:KFDIA,1:NDTPREC2)=YDMF_PHYS_SURF%GSD_XP2%PPRECIP2(KIDIA:KFDIA,1:NDTPREC2)
  ENDIF
  CALL CPXFU(YDGEOMETRY, YDXFU, YDPHY%YRDPRECIPS, KIDIA, KFDIA, YDCPG_BNDS%KIDIE, YDCPG_BNDS%KFDIE, KSTGLO, &
    & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%FPLCL, &
    & YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLCG, YDMF_PHYS%OUT%FPLCH, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, &
    & YDMF_PHYS%OUT%FPLSG, YDMF_PHYS%OUT%FPLSH, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDCPG_MISC%CLCT, &
    & ZSTRCU, ZSTRCV, ZSTRDU, ZSTRDV, ZSTRTU, ZSTRTV, YDCPG_MISC%NEB, YDCPG_MISC%QICE, YDCPG_MISC%QLI, &
    & ZUCLS, ZVCLS, ZNUCLS, ZNVCLS, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%QCLS, YDMF_PHYS%OUT%RHCLS, YDMF_PHYS%OUT%SIC, &
    & YDMF_PHYS%OUT%CLCH, YDMF_PHYS%OUT%CLCM, YDMF_PHYS%OUT%CLCL, YDMF_PHYS%OUT%CLCC, YDMF_PHYS%OUT%FCS, &
    & YDMF_PHYS%OUT%FONTE, YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FLWSP, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%RUISP, &
    & YDMF_PHYS%OUT%RUISS, YDMF_PHYS%OUT%CAPE, YDMF_PHYS%OUT%CTOP, YDMF_PHYS%OUT%MOCON, YDMF_PHYS%OUT%CLPH, &
    & YDMF_PHYS%OUT%VEIN, ZUGST, ZVGST, YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FTR, YDMF_PHYS%OUT%RUISL, &
    & YDMF_PHYS%OUT%DIAGH, YDMF_PHYS%OUT%MRT, YDMF_PHYS%OUT%VISICLD, YDMF_PHYS%OUT%VISIHYD, YDMF_PHYS%OUT%MXCLWC, &
    & YDMF_PHYS%OUT%TPWCLS, ZDPRECIPS, ZDPRECIPS2)
ENDIF

!*    2.6   INTERFACE FOR CLS FIELDS FOR ASSIMILATION

IF (LDIRCLSMOD) THEN
  IF (LDIRSICMOD) THEN
    CALL CPCLS_ASSIM(YDGEOMETRY, YDSURF, KIDIA, KFDIA, YDMF_PHYS_SURF%GSP_CL%PUCLS_T1, &
    & YDMF_PHYS_SURF%GSP_CL%PVCLS_T1, YDMF_PHYS_SURF%GSP_CL%PNUCLS_T1, YDMF_PHYS_SURF%GSP_CL%PNVCLS_T1, &
    & YDMF_PHYS_SURF%GSP_CL%PTCLS_T1, YDMF_PHYS_SURF%GSP_CL%PHUCLS_T1, ZUCLS, ZVCLS, ZNUCLS, ZNVCLS, &
    & YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%RHCLS, PXSIC=YDMF_PHYS%OUT%SIC)
  ELSE
    CALL CPCLS_ASSIM(YDGEOMETRY, YDSURF, KIDIA, KFDIA, YDMF_PHYS_SURF%GSP_CL%PUCLS_T1, &
    & YDMF_PHYS_SURF%GSP_CL%PVCLS_T1, YDMF_PHYS_SURF%GSP_CL%PNUCLS_T1, YDMF_PHYS_SURF%GSP_CL%PNVCLS_T1, &
    & YDMF_PHYS_SURF%GSP_CL%PTCLS_T1, YDMF_PHYS_SURF%GSP_CL%PHUCLS_T1, ZUCLS, ZVCLS, ZNUCLS, ZNVCLS, &
    & YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%RHCLS)
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       3.    RESETS DDH FLEXIBLE STRUCTURES:
!              -------------------------------

! resets ddh flexible structures if call only to phys

IF (YDCPG_OPTS%LCONFX.AND.LFLEXDIA.AND.(.NOT.LDDH_OMP)) CALL RESET_DDHFLEX

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_DIA', 1, ZHOOK_HANDLE)
END SUBROUTINE CPG_DIA
