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
SUBROUTINE CPG_GP_HYD(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDTMP, YDVARS, YDCPG_DYN0, YDCPG_DYN9, &
& YDMODEL, YDGMV, LDUVH, LDGW, LDLDIAB, LDMPA, POROG, PPGTERM, PGFL, PATND)

!**** *CPG_GP_HYD* - Grid point calculations at instants t and t-dt for hydrostatic model.

!     Purpose.
!     --------
!      Grid-point calculations at instants t and t-dt for hydrostatic model.
!      Computes some intermediate quantities.
!      Computes some Lagrangian equation RHS.

!**   Interface.
!     ----------
!        *CALL* *CPG_GP_HYD(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        KST       : first element of work.
!        KEND      : last element of work.
!        KSTGLO    : global offset.
!        LDUVH     : compute half level winds.
!        LDGW      : compute [gw].
!        LDLDIAB   : .T. if complete physics is activated and predictor step.
!        LDMPA     : AROME upper-air physics.
!        KIBL      : index into YRGSGEOM/YRCSGEOM types in YDGEOMETRY
!        POROG     : surface orography.
!        POROGL    : zonal component of "grad(surf orography)"
!        POROGM    : meridian component of "grad(surf orography)"
!        PRE0L     : zonal component of "grad prehyds" at t.
!        PRE0M     : meridian component of "grad prehyds" at t.

!     INPUT/OUTPUT:
!     -------------
!        PGFL      : unified_treatment grid-point fields at t
!        PGMV      : upper air GMV variables at time t and t-dt.
!        PPGTERM   : pressure gradient quantities, only used if LPGFSF=T
!        PRE0      : hydrostatic pressure "prehyd" at half levels at time t.
!        PRE9      : hydrostatic pressure "prehyd" at half levels at t-dt.

!     OUTPUT:
!     -------
!        PRE0F     : hydrostatic pressure "prehyd" at full levels at time t.
!        PXYB0     : contains pressure depth, "delta", "alpha" at t.
!        PUVH0     : horizontal wind at time t at half levels.
!        PHIF0     : geopotential height "gz" at t at full levels.
!        PRCP0     : contains "cp", "R" and "Kap=R/Cp" at t.
!        PCTY0     : contains vertical velocities, vertical integral of divergence at t.
!        PGWFT0    : [gw] at full levels at t.
!        PKENE0    : kinetic energy at t.
!        PRE9F     : hydrostatic pressure "prehyd" at full levels at time t-dt.
!        PXYB9     : contains pressure depth, "delta", "alpha" at t-dt.
!        PHI9      : geopotential height "gz" at t-dt at half levels.
!        PHIF9     : geopotential height "gz" at t-dt at full levels.
!        PRCP9     : contains "cp", "R" and "Kap=R/Cp" at t-dt.
!        PGWFT9    : [gw] at full levels at t-dt.
!        PATND     : adiabatic Lagrangian tendencies.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE documentation

!     Author.
!     -------
!        K. YESSAD, after part 3 of CPG_GP.
!        Original : March 2017

! Modifications
! -------------
!      K. Yessad (Feb 2018): remove deep-layer formulations.
!      K. Yessad (Apr 2018): introduce key L_RDRY_VD (ensure consistent definition of "dver" everywhere).
!      R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!      I. Polichtchouk (June 2020): add LPGFSF option
!      F. Vana 21-Sep-2020: LPGREUSE specific dataflow
!   P. Smolikova (Sep 2020): Remove obsolete calculations.
!   R. El Khatib 08-Jul-2022 Cleaning
!      F. Vana 20-Mar-2024: Fixed LPGFSF and LPGREUSE cleaned out
! End Modifications
!     ------------------------------------------------------------------

USE GEOMETRY_MOD, ONLY : GEOMETRY
USE FIELD_VARIABLES_MOD,ONLY : FIELD_VARIABLES
USE CPG_TYPE_MOD,ONLY : CPG_DYN_TYPE, CPG_GP_TMP_TYPE
USE CPG_OPTS_TYPE_MOD, ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE YOMGMV,       ONLY : TGMV
USE TYPE_MODEL,   ONLY : MODEL

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE SC2PRG_MOD, ONLY : SC2PRG
USE YOMCT0    , ONLY : LSFORC
USE YOMCST    , ONLY : RD, RV
USE YOMCT3    , ONLY : NSTEP
USE YOMLSFORC , ONLY : LSW_FRC, LSOMEGA_FRC
USE YOMGWDIAG , ONLY : UPDATE_GWDIAG,LGWDIAGS_ON
USE INTDYN_MOD, ONLY : TPG_TYPE, YYTHW0, YYTHW9,&
 & YYTCTY0, YYTXYBDER0, YYTRCP0, YYTRCP9, YYTXYB0, YYTXYB9

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_BNDS_TYPE),INTENT(IN)   :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE),INTENT(IN)   :: YDCPG_OPTS
TYPE(CPG_GP_TMP_TYPE),INTENT(INOUT) :: YDTMP
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN0
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN9
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
LOGICAL           ,INTENT(IN)    :: LDUVH
LOGICAL           ,INTENT(IN)    :: LDGW
LOGICAL           ,INTENT(IN)    :: LDLDIAB
LOGICAL           ,INTENT(IN)    :: LDMPA
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(YDGEOMETRY%YRDIM%NPROMA)
TYPE(TPG_TYPE)    ,INTENT(IN)    :: PPGTERM
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PATND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_DYN%YRDYNA%YYTTND%NDIM)
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV,JROF
LOGICAL :: LLDER,LLGWF0,LLGDWI

!    ==== MISCELLANEOUS ======================

REAL(KIND=JPRB), POINTER :: ZP1FORC(:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gphpre.intfb.h"
#include "gpgrxyb.intfb.h"
#include "gprcp.intfb.h"
#include "gprt.intfb.h"
#include "gpcty.intfb.h"
#include "gpcty_forc.intfb.h"
#include "gpgeo.intfb.h"
#include "gpgrgeo.intfb.h"
#include "gphlwi.intfb.h"
#include "gphluv.intfb.h"
#include "gpuvs.intfb.h"
#include "gpgrp.intfb.h"
#include "gpxx.intfb.h"
#include "gpgw.intfb.h"
#include "gp_tndlagadiab_uv.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG_GP_HYD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDGSGEOM=>YDGEOMETRY%YRGSGEOM(YDCPG_BNDS%KBL), &
 & YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, YGFL=>YDMODEL%YRML_GCONF%YGFL, &
 & YDCVER=>YDGEOMETRY%YRCVER, YDDYNA=>YDMODEL%YRML_DYN%YRDYNA, &
 & KIDIA=>YDCPG_BNDS%KIDIA, KFDIA=>YDCPG_BNDS%KFDIA, KSTGLO=>YDCPG_BNDS%KSTGLO)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, NFLEVG=>YDDIMV%NFLEVG, YRSPEC=>YGFL%YRSPEC, YFORC=>YGFL%YFORC,&
 & NCURRENT_ITER=>YDMODEL%YRML_DYN%YRDYN%NCURRENT_ITER, &
 & LSIMPH=>YDMODEL%YRML_PHY_MF%YRSIMPHL%LSIMPH)

CALL SC2PRG(1, YFORC(:)%MP, PGFL, ZP1FORC)

!     ------------------------------------------------------------------

!*       3.    INITIALISE AUXILIARY VARIABLES AND TENDENCIES.
!              ----------------------------------------------

!---------------------------------------------------
!      3.1         TIME t0 CALCULATIONS
!---------------------------------------------------

!*     3.1.1 COMPUTE PRE0, PXYB0, PRE0F.

CALL GPHPRE(NPROMA, NFLEVG, KIDIA, KFDIA, YDVAB, YDCVER, YDCPG_DYN0%PRE, PRESF=YDCPG_DYN0%PREF, &
& PDELP=YDCPG_DYN0%XYB%DELP, PLNPR=YDCPG_DYN0%XYB%LNPR, PRDELP=YDCPG_DYN0%XYB%RDELP, PALPH=YDCPG_DYN0%XYB%ALPH, &
& PRTGR=YDCPG_DYN0%XYB%RTGR, PRPRE=YDCPG_DYN0%XYB%RPRE, PRPP=YDCPG_DYN0%XYB%RPP)

!*     3.1.2 COMPUTE ZRPRE0F AND ZXYBDER0.

IF(YDCVER%LVERTFE) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KIDIA,KFDIA
      YDTMP%T0%RPREF(JROF,JLEV)=1.0_JPRB/YDCPG_DYN0%PREF(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF
CALL GPGRXYB(NPROMA, KIDIA, KFDIA, NFLEVG, .FALSE., YDVAB, YDCVER, YDCPG_DYN0%PREL, YDCPG_DYN0%PREM,  &
  & PDELP=YDCPG_DYN0%XYB%DELP, PLNPR=YDCPG_DYN0%XYB%LNPR, PRDELP=YDCPG_DYN0%XYB%RDELP, PALPH=YDCPG_DYN0%XYB%ALPH, &
  & PRTGR=YDCPG_DYN0%XYB%RTGR, PRPRE=YDCPG_DYN0%XYB%RPRE, PRPP=YDCPG_DYN0%XYB%RPP, PCOEFD_DER=YDTMP%T0%XYBDER%COEFD,  &
  & PLNPRL_DER=YDTMP%T0%XYBDER%LNPRL, PLNPRM_DER=YDTMP%T0%XYBDER%LNPRM, PCOEFA_DER=YDTMP%T0%XYBDER%COEFA, &
  & PCOEFAPL_DER=YDTMP%T0%XYBDER%COEFAPL, PALPHPLL_DER=YDTMP%T0%XYBDER%ALPHPLL, PALPHPLM_DER=YDTMP%T0%XYBDER%ALPHPLM, &
  & PALPHL_DER=YDTMP%T0%XYBDER%ALPHL, PALPHM_DER=YDTMP%T0%XYBDER%ALPHM)

!*     3.1.3 COMPUTE "R", "Cp" AND "kap=R/Cp".

CALL GPRCP(NPROMA, KIDIA, KFDIA, NFLEVG, PGFL=PGFL, PCP=YDCPG_DYN0%RCP%CP, PR=YDCPG_DYN0%RCP%R, &
& PKAP=YDCPG_DYN0%RCP%KAP)

!*     3.1.4 COMPUTES "RT" AND ITS GRADIENT.

IF(YRSPEC%LGP) THEN
  CALL GPRT(YDDYNA%LSPRT, NPROMA, KIDIA, KFDIA, NFLEVG, RD, RV, YDCPG_DYN0%RCP%R, YDVARS%T%T0, &
  & YDVARS%T%DL, YDVARS%T%DM, YDVARS%Q%DL, YDVARS%Q%DM, YDTMP%T0%RT, YDCPG_DYN0%RTL, YDCPG_DYN0%RTM, &
  & PRL=YDVARS%RSPEC%DL, PRM=YDVARS%RSPEC%DM)
ELSE
  CALL GPRT(YDDYNA%LSPRT, NPROMA, KIDIA, KFDIA, NFLEVG, RD, RV, YDCPG_DYN0%RCP%R, YDVARS%T%T0, &
  & YDVARS%T%DL, YDVARS%T%DM, YDVARS%Q%DL, YDVARS%Q%DM, YDTMP%T0%RT, YDCPG_DYN0%RTL, YDCPG_DYN0%RTM)
ENDIF

! ZRDT0 is the (RT) at "t" used in definition of "dver".
IF (YDDYNA%L_RDRY_VD) THEN
  YDTMP%T0%RDT(KIDIA:KFDIA,1:NFLEVG)=RD*YDVARS%T%T0(KIDIA:KFDIA,1:NFLEVG)
ELSE
  YDTMP%T0%RDT(KIDIA:KFDIA,1:NFLEVG)=YDTMP%T0%RT(KIDIA:KFDIA,1:NFLEVG)
ENDIF

!*     3.1.7c COMPUTATION OF SOME VERTICAL VELOCITIES.
!*           Solve continuity equation
!            Computation of the vertical velocities "etapt (d prehyd / d eta)"
!            "omega/prehyd" and "W".

!* Mass tendencies from physics as RHS of the continuity equation.
! Specific mass tendencies (unit: s-1) are computed
! in the last physics call (previous time step in IFS).
! They are stored in GFL YPHYCTY.
IF (YGFL%YPHYCTY%LACTIVE) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KIDIA,KFDIA
      YDTMP%T0%DPHYCTY(JROF,JLEV)=YDVARS%PHYCTY%T0(JROF,JLEV)
    ENDDO
  ENDDO
ELSE
  YDTMP%T0%DPHYCTY(:,:)=0.0_JPRB
ENDIF

CALL GPCTY(YDVFE, YDCVER, NPROMA, KIDIA, KFDIA, NFLEVG, YDDYNA%LRUBC, YDVAB, YDVETA, YDVARS%U%T0, &
& YDVARS%V%T0, YDVARS%DIV%T0, YDVARS%EDOT%T0, PSPL=YDCPG_DYN0%PREL, PSPM=YDCPG_DYN0%PREM, PRPREF=YDTMP%T0%RPREF, &
& PDPHYCTY=YDTMP%T0%DPHYCTY, PDELP=YDCPG_DYN0%XYB%DELP, PLNPR=YDCPG_DYN0%XYB%LNPR, PRDELP=YDCPG_DYN0%XYB%RDELP, &
& PALPH=YDCPG_DYN0%XYB%ALPH, PRTGR=YDCPG_DYN0%XYB%RTGR, PRPRE=YDCPG_DYN0%XYB%RPRE, PRPP=YDCPG_DYN0%XYB%RPP, &
& PEVEL=YDCPG_DYN0%CTY%EVEL, PVVEL=YDCPG_DYN0%CTY%VVEL, PPSDIV=YDCPG_DYN0%CTY%PSDIV, PPSDVBC=YDCPG_DYN0%CTY%PSDVBC, &
& PDIVDP=YDCPG_DYN0%CTY%DIVDP)
 

! Diagnostics of gravity wave noise reflected in surface pressure tendency
IF(LGWDIAGS_ON) THEN
  CALL UPDATE_GWDIAG(YDGEOMETRY, YDMODEL%YRML_GCONF%YRRIP, NPROMA, NFLEVG, KFDIA, NSTEP, KSTGLO, &
  & YDCPG_DYN0%PRE, YDCPG_DYN0%CTY%DIVDP(:, 1:))
ENDIF

! Overwrite PCTY0 for "EVEL" and "VVEL" with large scale forced values in case
! of 1D model (SCUM)
! (in SCUM, these are 0 in GPCTY)
IF ( LSFORC .AND. (LSW_FRC.OR.LSOMEGA_FRC) ) THEN
  CALL GPCTY_FORC(YDGEOMETRY, YDMODEL%YRML_GCONF, YDMODEL%YRML_PHY_MF%YRPHY2, KIDIA, KFDIA, &
  & ZP1FORC, YDCPG_DYN0%PREF, YDVARS%T%T0, YDCPG_DYN0%RCP%R, YDCPG_DYN0%CTY%EVEL, YDCPG_DYN0%CTY%VVEL(:,1:))
ENDIF

!*     3.1.8 COMPUTES THE GEOPOTENTIAL HEIGHT "gz" AND ITS GRADIENT.

!          if LPGFSF=T, these already computed in gptosptogp.F90 
!           make sure the two computations result in identity

  ! * "gz" at full levels and half levels.

YDCPG_DYN0%PHI(KIDIA:KFDIA,NFLEVG)=POROG(KIDIA:KFDIA)
CALL GPGEO(NPROMA, KIDIA, KFDIA, NFLEVG, YDCPG_DYN0%PHI, YDCPG_DYN0%PHIF, YDVARS%T%T0, &
& YDCPG_DYN0%RCP%R, YDCPG_DYN0%XYB%LNPR, YDCPG_DYN0%XYB%ALPH, YDGEOMETRY%YRVERT_GEOM)

  ! * "grad gz" at full levels and half levels.

CALL GPGRGEO(YDGEOMETRY, NPROMA, KIDIA, KFDIA, NFLEVG, YDTMP%T0%RT, YDCPG_DYN0%RTL, &
& YDCPG_DYN0%RTM, YDCPG_DYN0%XYB%LNPR, YDCPG_DYN0%XYB%ALPH, POROGL=YDCPG_DYN0%OROGL, POROGM=YDCPG_DYN0%OROGM, &
& PHIFL=YDCPG_DYN0%PHIFL, PHIFM=YDCPG_DYN0%PHIFM, PHIHL=YDTMP%T0%GPHL, PHIHM=YDTMP%T0%GPHM, PLNPRL_DER=YDTMP%T0%XYBDER%LNPRL, &
& PLNPRM_DER=YDTMP%T0%XYBDER%LNPRM, PALPHL_DER=YDTMP%T0%XYBDER%ALPHL, PALPHM_DER=YDTMP%T0%XYBDER%ALPHM)

!*     3.1.9 COMPUTES KINETIC ENERGY.

DO JLEV=1,NFLEVG
  DO JROF=KIDIA,KFDIA
    YDCPG_DYN0%KENE(JROF,JLEV)=0.5_JPRB*(YDVARS%U%T0(JROF,JLEV)**2+YDVARS%V%T0(JROF,JLEV)**2)
  ENDDO
ENDDO

!*     3.1.10 COMPUTES HALF-LEVEL WINDS.

IF (LDUVH) THEN
  CALL GPHLWI(YDGEOMETRY%YRDIMV, NPROMA, KIDIA, KFDIA, YDCPG_DYN0%XYB%LNPR, YDCPG_DYN0%XYB%ALPH, &
  & YDCPG_DYN0%UVH%WWI (1:, 1:))
  CALL GPHLUV(YDGEOMETRY%YRDIMV, NPROMA, KIDIA, KFDIA, YDVARS%U%T0, YDVARS%V%T0, &
  & PWWI=YDCPG_DYN0%UVH%WWI, PUH=YDCPG_DYN0%UVH%UH, PVH=YDCPG_DYN0%UVH%VH)
ENDIF
IF (LDGW.OR.YDDYNA%LSLHD) THEN
  LLDER=.FALSE.
  CALL GPUVS(NFLEVG, NPROMA, KIDIA, KFDIA, LLDER, YDVARS%U%T0, YDVARS%V%T0, YDTMP%T0%US, YDTMP%T0%VS)
ENDIF

!*     3.1.12 COMPUTES THE PRESSURE FORCE FOR THE RHS OF MOMENTUM EQN.

IF (YDDYNA%LPGFSF) THEN
  !  use the filtered pressure gradient force
  YDTMP%T0%PSGRTL(KIDIA:KFDIA,1:NFLEVG)=PPGTERM%PSGRTL(KIDIA:KFDIA,1:NFLEVG)
  YDTMP%T0%PSGRTM(KIDIA:KFDIA,1:NFLEVG)=PPGTERM%PSGRTM(KIDIA:KFDIA,1:NFLEVG)
ELSE
  CALL GPGRP(YDGEOMETRY, KIDIA, KFDIA, YDTMP%T0%RT, YDCPG_DYN0%RTL, YDCPG_DYN0%RTM, &
  & YDCPG_DYN0%PREL, YDCPG_DYN0%PREM, YDCPG_DYN0%XYB%RTGR, YDCPG_DYN0%XYB%ALPH, PHIHL=YDTMP%T0%GPHL, PHIHM=YDTMP%T0%GPHM, &
  & PHIFL=YDCPG_DYN0%PHIFL, PHIFM=YDCPG_DYN0%PHIFM, PSGRTL=YDTMP%T0%PSGRTL, PSGRTM=YDTMP%T0%PSGRTM, PALPHPLL_DER=YDTMP%T0%XYBDER%ALPHPLL, &
  & PALPHPLM_DER=YDTMP%T0%XYBDER%ALPHPLM)
ENDIF

! store for next time-step
IF( YDDYNA%LGRADSP ) THEN

  ! adjust field with filter result
  YDTMP%T0%SGRTL(KIDIA:KFDIA,1:NFLEVG)=YDTMP%T0%PSGRTL(KIDIA:KFDIA,1:NFLEVG)&
   & - (YDVARS%SGRTL%T9(KIDIA:KFDIA,1:NFLEVG) - YDVARS%SGRTL%T0(KIDIA:KFDIA,1:NFLEVG))
  YDTMP%T0%SGRTM(KIDIA:KFDIA,1:NFLEVG)=YDTMP%T0%PSGRTM(KIDIA:KFDIA,1:NFLEVG)&
   & - (YDVARS%SGRTM%T9(KIDIA:KFDIA,1:NFLEVG) - YDVARS%SGRTM%T0(KIDIA:KFDIA,1:NFLEVG))

  ! store for next time-step unfiltered fields
  YDVARS%SGRTL%T0(KIDIA:KFDIA,1:NFLEVG)=YDTMP%T0%PSGRTL(KIDIA:KFDIA,1:NFLEVG)
  YDVARS%SGRTM%T0(KIDIA:KFDIA,1:NFLEVG)=YDTMP%T0%PSGRTM(KIDIA:KFDIA,1:NFLEVG)

  ! set new values for current timestep
  YDTMP%T0%PSGRTL(KIDIA:KFDIA,1:NFLEVG)=YDTMP%T0%SGRTL(KIDIA:KFDIA,1:NFLEVG)
  YDTMP%T0%PSGRTM(KIDIA:KFDIA,1:NFLEVG)=YDTMP%T0%SGRTM(KIDIA:KFDIA,1:NFLEVG)

ENDIF

!*     3.1.13 COMPUTES TERM "NHX".

IF (LDGW) THEN
  CALL GPXX(YDGEOMETRY, NFLEVG, NPROMA, KIDIA, KFDIA, YDTMP%T0%GPHL, YDTMP%T0%GPHM,&
  & YDCPG_DYN0%PHIFL, YDCPG_DYN0%PHIFM, YDCPG_DYN0%XYB%LNPR, YDTMP%T0%RT, YDVARS%U%T0, YDVARS%V%T0, YDCPG_DYN0%UVH%UH, &
  & YDCPG_DYN0%UVH%VH, YDTMP%T0%NHXT)

!*     3.1.14 COMPUTES dver.

  IF (YDDYNA%L_RDRY_VD) THEN
    ! dver_hyd=-(R/Rd)*(D+X+(cv/cp)(omega/prehyd))
    DO JLEV=1,NFLEVG
      DO JROF=KIDIA,KFDIA
        YDTMP%T0%DVER(JROF,JLEV)=-(YDCPG_DYN0%RCP%R(JROF,JLEV)/RD)*(YDVARS%DIV%T0(JROF,JLEV)&
         & +YDTMP%T0%NHXT(JROF,JLEV)+(1.0_JPRB-YDCPG_DYN0%RCP%KAP(JROF,JLEV))*YDCPG_DYN0%CTY%VVEL(JROF,JLEV))
      ENDDO
    ENDDO
  ELSE
    ! dver_hyd=-(D+X+(cv/cp)(omega/prehyd))
    DO JLEV=1,NFLEVG
      DO JROF=KIDIA,KFDIA
        YDTMP%T0%DVER(JROF,JLEV)=-(YDVARS%DIV%T0(JROF,JLEV)&
         & +YDTMP%T0%NHXT(JROF,JLEV)+(1.0_JPRB-YDCPG_DYN0%RCP%KAP(JROF,JLEV))*YDCPG_DYN0%CTY%VVEL(JROF,JLEV))
      ENDDO
    ENDDO
  ENDIF

!*     3.1.16 COMPUTES term [gw].

! Caution: must use consistent definitions of RT in ZRDT0 and ZDVER0.
  LLGWF0=.TRUE.
  LLGDWI=.FALSE.
  CALL GPGW(YDGEOMETRY, YDDYNA%LNHDYN, NFLEVG, NPROMA, KIDIA, KFDIA, LLGWF0, LLGDWI, YDCPG_DYN0%OROGL, &
  & YDCPG_DYN0%OROGM, YDCPG_DYN0%XYB%LNPR, YDCPG_DYN0%XYB%ALPH, YDTMP%T0%US, YDTMP%T0%VS, YDTMP%T0%RDT, &
  & YDTMP%T0%DVER, YDTMP%T0%GWHT, YDCPG_DYN0%GWFT)
ENDIF

!*     3.1.20a COMPUTES [D V/Dt].

! Pressure gradient term + explicit Coriolis + Rayleigh friction
CALL GP_TNDLAGADIAB_UV(YDGEOMETRY, YDGMV, YDMODEL%YRML_PHY_EC%YREPHY, YDMODEL%YRML_DYN%YRDYN, KIDIA, KFDIA, &
 & YDGSGEOM%RCORI, YDGSGEOM%GNORDL, YDGSGEOM%GNORDM, YDTMP%T0%PSGRTL, YDTMP%T0%PSGRTM, YDVARS%U%T0, YDVARS%V%T0, &
 & PATND(1,1,YDDYNA%YYTTND%M_TNDU), PATND(1,1,YDDYNA%YYTTND%M_TNDV), &
 & PATND(1,1,YDDYNA%YYTTND%M_TNDU_NOC), PATND(1,1,YDDYNA%YYTTND%M_TNDV_NOC))

!*     3.1.24 COMPUTES [DT/Dt].

DO JLEV=1,NFLEVG
  DO JROF=KIDIA,KFDIA
    PATND(JROF,JLEV,YDDYNA%YYTTND%M_TNDT)=&
     & YDCPG_DYN0%RCP%KAP(JROF,JLEV)*YDVARS%T%T0(JROF,JLEV)*YDCPG_DYN0%CTY%VVEL(JROF,JLEV)
  ENDDO
ENDDO

!---------------------------------------------------
!      3.2         TIME t9 CALCULATIONS
!---------------------------------------------------

IF (.NOT.YDDYNA%LTWOTL .AND. NCURRENT_ITER == 0) THEN

!*     3.2.1 COMPUTE PRE9, PXYB9, PRE9F.

  CALL GPHPRE(NPROMA, NFLEVG, KIDIA, KFDIA, YDVAB, YDCVER, YDCPG_DYN9%PRE, PRESF=YDCPG_DYN9%PREF, &
  & PDELP=YDCPG_DYN9%XYB%DELP, PLNPR=YDCPG_DYN9%XYB%LNPR, PRDELP=YDCPG_DYN9%XYB%RDELP, PALPH=YDCPG_DYN9%XYB%ALPH, &
  & PRTGR=YDCPG_DYN9%XYB%RTGR, PRPRE=YDCPG_DYN9%XYB%RPRE, PRPP=YDCPG_DYN9%XYB%RPP)

!*     3.2.3 COMPUTE "R", "Cp" AND "Kap=R/Cp".

  CALL GPRCP(NPROMA, KIDIA, KFDIA, NFLEVG, PGFL=PGFL, KGFLTYP=9, PCP=YDCPG_DYN9%RCP%CP, &
  & PR=YDCPG_DYN9%RCP%R)

!*     3.2.4 COMPUTE "RT".

  ! ZR0T9 is used in calculation of NHX(t-dt).
  YDTMP%T9%R0T(KIDIA:KFDIA,1:NFLEVG)=YDVARS%T%T9(KIDIA:KFDIA,1:NFLEVG)*YDCPG_DYN0%RCP%R(KIDIA:KFDIA,1:NFLEVG)

  ! ZRDT9 is the (RT) at "t-dt" used in the definition of "dver":
  IF (YDDYNA%L_RDRY_VD) THEN
    YDTMP%T9%RDT(KIDIA:KFDIA,1:NFLEVG)=RD*YDVARS%T%T9(KIDIA:KFDIA,1:NFLEVG)
  ELSE
    ! use of R_moist at time t to be consistent with ZR0T9 and ZDVER9.
    DO JLEV=1,NFLEVG
      DO JROF=KIDIA,KFDIA
        YDTMP%T9%RDT(JROF,JLEV)=YDCPG_DYN0%RCP%R(JROF,JLEV)*YDVARS%T%T9(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF

!*     3.2.8 COMPUTES THE GEOPOTENTIAL HEIGHT "gz" AND ITS GRADIENT.

  ! * "gz" at full levels and half levels.
  IF(LDLDIAB.OR.LSIMPH) THEN
    YDCPG_DYN9%PHI(KIDIA:KFDIA,NFLEVG)=POROG(KIDIA:KFDIA)
    CALL GPGEO(NPROMA, KIDIA, KFDIA, NFLEVG, YDCPG_DYN9%PHI, YDCPG_DYN9%PHIF, YDVARS%T%T9, &
    & YDCPG_DYN9%RCP%R, YDCPG_DYN9%XYB%LNPR, YDCPG_DYN9%XYB%ALPH, YDGEOMETRY%YRVERT_GEOM)
  ENDIF

!*     3.2.10 COMPUTES HALF-LEVEL WINDS.

  IF (LDGW) THEN
    ! * same remark as for time t half-level winds.
    YDTMP%T9%UVH%WWI(KIDIA:KFDIA,1:NFLEVG)=YDCPG_DYN0%UVH%WWI(KIDIA:KFDIA,1:NFLEVG)
    CALL GPHLUV(YDGEOMETRY%YRDIMV, NPROMA, KIDIA, KFDIA, YDVARS%U%T9, YDVARS%V%T9, &
    & PWWI=YDTMP%T9%UVH%WWI, PUH=YDTMP%T9%UVH%UH, PVH=YDTMP%T9%UVH%VH)

    LLDER=.FALSE.
    CALL GPUVS(NFLEVG, NPROMA, KIDIA, KFDIA, LLDER, YDVARS%U%T9, YDVARS%V%T9, YDTMP%T9%US, &
    & YDTMP%T9%VS)

!*     3.2.13 COMPUTES TERM "NHX".

  ! caution: must use "R" at time t.
    CALL GPXX(YDGEOMETRY, NFLEVG, NPROMA, KIDIA, KFDIA, YDTMP%T0%GPHL, YDTMP%T0%GPHM, &
    & YDCPG_DYN0%PHIFL, YDCPG_DYN0%PHIFM, YDCPG_DYN9%XYB%LNPR, YDTMP%T9%R0T, YDVARS%U%T9, YDVARS%V%T9, &
    & YDTMP%T9%UVH%UH, YDTMP%T9%UVH%VH, YDTMP%T9%NHXT)

!*     3.2.14 COMPUTES dver.

    IF (YDDYNA%L_RDRY_VD) THEN
      ! dver_hyd=-(R/Rd)*(D+X+(cv/cp)(omega/prehyd))
      !  Caution: must use R, Kap and (omega/prehyd) at time t.
      DO JLEV=1,NFLEVG
        DO JROF=KIDIA,KFDIA
          YDTMP%T9%DVER(JROF,JLEV)=-(YDCPG_DYN0%RCP%R(JROF,JLEV)/RD)*(YDVARS%DIV%T9(JROF,JLEV)&
           & +YDTMP%T9%NHXT(JROF,JLEV)+(1.0_JPRB-YDCPG_DYN0%RCP%KAP(JROF,JLEV))*YDCPG_DYN0%CTY%VVEL(JROF,JLEV))
        ENDDO
      ENDDO
    ELSE
      ! dver_hyd=-(D+X+(cv/cp)(omega/prehyd))
      !  Caution: must use R, Kap and (omega/prehyd) at time t.
      DO JLEV=1,NFLEVG
        DO JROF=KIDIA,KFDIA
          YDTMP%T9%DVER(JROF,JLEV)=-(YDVARS%DIV%T9(JROF,JLEV)&
           & +YDTMP%T9%NHXT(JROF,JLEV)+(1.0_JPRB-YDCPG_DYN0%RCP%KAP(JROF,JLEV))*YDCPG_DYN0%CTY%VVEL(JROF,JLEV))
        ENDDO
      ENDDO
    ENDIF
  ENDIF

!*     3.2.16 COMPUTES term [gw].

  ! Caution: must use consistent definitions of RT in ZRDT9 and ZDVER9.
  LLGDWI=.FALSE.
  IF (LDGW) THEN
    CALL GPGW(YDGEOMETRY, YDDYNA%LNHDYN, NFLEVG, NPROMA, KIDIA, KFDIA, LDMPA, LLGDWI, YDCPG_DYN0%OROGL, &
    & YDCPG_DYN0%OROGM, YDCPG_DYN9%XYB%LNPR, YDCPG_DYN9%XYB%ALPH, YDTMP%T9%US, YDTMP%T9%VS, YDTMP%T9%RDT, &
    & YDTMP%T9%DVER, YDTMP%T9%GWHT, YDCPG_DYN9%GWFT)
  ENDIF

ENDIF

!     ------------------------------------------------------------------
! Cleaning
YDTMP%T0%GWHT(:,:)=0.0_JPRB
YDTMP%T9%GWHT(:,:)=0.0_JPRB

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_GP_HYD',1,ZHOOK_HANDLE)
END SUBROUTINE CPG_GP_HYD
