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
SUBROUTINE CPG_GP(YDGEOMETRY,YDMODEL,YDFIELDS,&
 !---------------------------------------------------------------------
 ! - INPUT .
 & LD_DFISTEP,KST,KEND,KSTC,KENDC,KBL,KSTGLO,LDLFSTEP,LDLDIAB,PDT,PTE,&
 & KIBL,&
 !---------------------------------------------------------------------
 ! - INPUT/OUTPUT .
 & PPGTERM,PGFL,PGMV,PGMVS,&
 & PSP_RR,PSP_SB,PSP_SG,PSD_VF,PSD_VV,&
 & PGMVTNDSL,PGFLTNDSL,&
 !---------------------------------------------------------------------
 ! - OUTPUT .
 & POROGL,POROGM,&
 & PRE0,PRE0L,PRE0M,PRE0F,PNHPRE0F,PNHPRE0H,PXYB0,PUVH0,&
 & PHIF0,PRCP0,PCTY0,PGWFT0,PKENE0,&
 & PRE9,PRE9F,PNHPRE9F,PNHPRE9H,PXYB9,&
 & PHI9,PHIF9,PRCP9,PGWFT9,&
 & PGPAR,PB2,&
 & PGMVT1,PGMVT1S,PGFLT1,&
 & PKOZO,PQS,PQICE,PQLI,PQRAIN,PQSNOW,&
 & PATND,PDBBC,PRDPHI,PGWT0,PGWT9,PGWS,PGWFL,PGWFM,&
 & PNHXT0,PNHXT9,PQCHA0L,PQCHA0M,PEXTRA,YDDDH)

!**** *CPG_GP* - Grid point calculations:
!                initial part of not lagged grid-point calculations.

!     Purpose.
!     --------
!           Grid point calculations:
!           initial part of not lagged grid-point calculations.
!           - get data in buffers.
!           - multiply p-order horizontal derivatives by M**p.
!           - second part of the temporal filter.
!           - grid-point calculations for nudging.
!           - calls some GP... routines to initialise some auxiliary variables.
!           - sets-up and PB2.

!           Abbreviation "vwv" stands for "vertical wind variable".

!**   Interface.
!     ----------
!        *CALL* *CPG_GP(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        LD_DFISTEP   : 'D' -> DFI computations
!        KST       : first element of work.
!        KEND      : last element of work.
!        KSTC,KENDC: the same as KST,KEND but including zone C for ALADIN.
!        KBL       : block number.
!        KSTGLO    : global offset.
!        LDLFSTEP  : .T.: first time-step?
!        LDLDIAB   : .T. if complete physics is activated and predictor step.
!        PDT       : For a leap-frog scheme (three time level scheme):
!                     'dt' at the first time-step, '2 dt' otherwise.
!                    For a 2TL SL scheme: timestep 'dt'.
!        PTE       : 1. or 0. according to different configurations.
!        KIBL      : index into YRGSGEOM/YRCSGEOM types in YDGEOMETRY 

!     INPUT/OUTPUT:
!     -------------
!        PGFL      : unified_treatment grid-point fields at t
!        PGMV      : upper air GMV variables at time t and t-dt.
!        PGMVS     : surface GMV variables at time t and t-dt.
!        PPGTERM   : pressure gradient quantities, only used if LPGFSF=T
!        PSP_RR    : surface prognostic quantities
!        PSP_SB    : soil prognostic quantities for the different reservoirs
!        PSP_SG    : surface snow prognostic quantities
!        PSD_VF    : climatological/geographical diagnostic fields
!        PSD_VV    : vegetation diagnostic fields
!        PGMVTNDSL : GMV(t+dt,F)-GMV(t or t-dt,O) for DDH
!        PGFLTNDSL : GFL(t+dt,F)-GFL(t or t-dt,O) for DDH

!     OUTPUT:
!     -------
!        POROGL    : zonal component of "grad(surf orography)"
!        POROGM    : meridian component of "grad(surf orography)"
!        PRE0      : hydrostatic pressure "prehyd" at half levels at time t.
!        PRE0L     : zonal component of "grad prehyds" at t.
!        PRE0M     : meridian component of "grad prehyds" at t.
!        PRE0F     : hydrostatic pressure "prehyd" at full levels at time t.
!        PNHPRE0F  : "pre" at full levels (time t).
!        PNHPRE0H  : "pre" at half levels (time t).
!        PXYB0     : contains pressure depth, "delta", "alpha" at t.
!        PUVH0     : horizontal wind at time t at half levels.
!        PHIF0     : geopotential height "gz" at t at full levels.
!        PRCP0     : contains "cp", "R" and "Kap=R/Cp" at t.
!        PCTY0     : contains vertical velocities, vertical integral of divergence at t.
!        PGWFT0    : [Gw] at full layers at t.
!        PKENE0    : kinetic energy at t.
!        PRE9      : hydrostatic pressure "prehyd" at half levels at t-dt.
!        PRE9F     : hydrostatic pressure "prehyd" at full levels at time t-dt.
!        PNHPRE9F  : "pre" at full levels (time t-dt).
!        PNHPRE9H  : "pre" at half levels (time t-dt).
!        PXYB9     : contains pressure depth, "delta", "alpha" at t-dt.
!        PHI9      : geopotential height "gz" at t-dt at half levels.
!        PHIF9     : geopotential height "gz" at t-dt at full levels.
!        PRCP9     : contains "cp", "R" and "Kap=R/Cp" at t-dt.
!        PGWFT9    : [Gw] at full levels at t-dt.
!        PGPAR     : surface fields for AROME.
!        PB2       : "SLB2" buffer.
!        PGMVT1    : upper air GMV variables at t+dt.
!        PGMVT1S   : surface GMV variables at t+dt.
!        PGFLT1    : GFL variables at t+dt.
!        PKOZO     : fields for photochemistery of ozon.
!        PQS       : specific humidity at surface level.
!        PQICE     : specific humidity of solid water for radiation.
!        PQLI      : specific humidity of liquid water for radiation.
!        PQRAIN    : specific humidity of rain for radiation.
!        PQSNOW    : specific humidity of snow for radiation.
!        PATND     : adiabatic Lagrangian tendencies.
!        PDBBC     : [D (Gw)_surf / Dt]_adiab.
!        PRDPHI    : HYD: not used.
!                    NHEE: contains pre/(R T prehyd [Delta log(prehyd)]) at t.
!                    NHQE: contains 1/(R Tt [Delta log(prehyd)]) at t.
!                    "R" is the version of R (may be Rdry or Rmoist) used in the definition of vertical divergence "dver".
!        PGWT0     : [Gw] at t (LGWADV=T only).
!        PGWT9     : [Gw] at t-dt (LGWADV=T only).
!                    for PGWT0 and PGWT9:
!                    * half level 0 to nflevg-1 values if LVFE_GW=F.
!                    * full level 1 to nflevg values if LVFE_GW=T.
!        PGWS      : [Gw]_surf at t (LRDBBC only).
!        PGWFL     : zonal comp grad(Gw) at full level at t.
!        PGWFM     : meridian comp grad(Gw) at full level at t.
!        PNHXT0    : term 'NHX' at t.
!        PNHXT9    : term 'NHX' at t-dt.
!        PQCHA0L   : zonal comp grad(log(pre/prehyd)).
!        PQCHA0M   : merid comp grad(log(pre/prehyd)).
!        PEXTRA    : additional quantity for diagnostics.
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
!        ARPEGE documentation vol 2 ch 1 and vol 3 ch 6

!     Author.
!     -------
!        K. YESSAD, after parts 1, 2 and 3 of old CPG. 
!        Original : 16-08-2001

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Sep 2008): prune enhanced diffusion (lfrein).
!   F. Bouyssel (Nov 2008): Removal of LPROCLD protection
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (March 2009): correct false comments for LRWSDLG=T
!   K. Yessad (Aug 2009): remove LPC_OLD in TL and AD codes.
!   K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP.
!   K. Yessad (Nov 2009): prune lpc_old.
!   K. Yessad (Dec 2009): LRWSDLW,LRWSDLR,LRWSDLG=T,T,T in NH model for LGWADV=F.
!   A. Alias  (Mar 2011) No nudging of SST when LSME=.T.
!   K. Yessad (Nov 2011): various contributions.
!   N. Wedi   (Nov 2011): add LGRADSP
!   M. Ahlgrimm  31-Oct-2011 add rain, snow and PEXTRA to DDH output
!   K. Yessad (Nov 2012): simplify testings.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!   K. Yessad (July 2014): Move some variables.
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   R. Roehrig (Sept 2018): add argument to CP_FORCING + add possibility to
!                           impose (time-evolving) surface pressure in MUSC
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!   I. Polichtchouk & F. Vana (Oct-2020): LPGFSF with LPGREUSE
!   I. Polichtchouk (Jul 2021): Add LSACC option
! End Modifications
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : GPPOPER
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE SC2PRG_MOD, ONLY : SC2PRG
USE YOMCT0    , ONLY : LSFORC
USE YOMCT3    , ONLY : NSTEP
USE YOMGPPB   , ONLY : GPARBUF
USE YOMNUD    , ONLY : NFNUDG, LNUDG
USE YOMLSFORC , ONLY : LSPS_FRC
USE YOMSNU    , ONLY : XPNUDG
USE INTDYN_MOD, ONLY : TPG_TYPE, YYTHW0, YYTHW9,&
 &                     YYTCTY0, YYTXYBDER0, YYTRCP0, YYTRCP9, YYTXYB0, YYTXYB9

USE DDH_MIX   , ONLY : TYP_DDH

! Added ---- 
USE FIELDS_MOD   ,ONLY : FIELDS
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
!---
TYPE(FIELDS)      ,INTENT(INOUT) :: YDFIELDS 
!---
LOGICAL           ,INTENT(IN)    :: LD_DFISTEP
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KENDC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTGLO 
LOGICAL           ,INTENT(IN)    :: LDLFSTEP 
LOGICAL           ,INTENT(IN)    :: LDLDIAB 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTE 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
TYPE(TPG_TYPE)    ,INTENT(INOUT) :: PPGTERM
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDFIELDS%YRGMV%NDIMGMV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDFIELDS%YRGMV%NDIMGMVS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSP_RR(YDGEOMETRY%YRDIM%NPROMA,YDFIELDS%YRSURF%YSP_RRD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSP_SB(YDGEOMETRY%YRDIM%NPROMA,YDFIELDS%YRSURF%YSP_SBD%NLEVS,YDFIELDS%YRSURF%YSP_SBD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSP_SG(YDGEOMETRY%YRDIM%NPROMA,YDFIELDS%YRSURF%YSP_SGD%NLEVS,YDFIELDS%YRSURF%YSP_SGD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_VF(YDGEOMETRY%YRDIM%NPROMA,YDFIELDS%YRSURF%YSD_VFD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_VV(YDGEOMETRY%YRDIM%NPROMA,YDFIELDS%YRSURF%YSD_VVD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVTNDSL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,2+YDMODEL%YRML_GCONF%YRDIMF%NFTHER)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLTNDSL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NUMFLDS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POROGL(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POROGM(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0L(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0M(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNHPRE0F(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNHPRE0H(YDGEOMETRY%YRDIM%NPROMNH,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXYB0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB0%NDIM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUVH0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTHW0%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PHIF0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRCP0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP0%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCTY0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTCTY0%NDIM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGWFT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKENE0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE9F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNHPRE9F(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNHPRE9H(YDGEOMETRY%YRDIM%NPROMNH,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXYB9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB9%NDIM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHI9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHIF9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRCP9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP9%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGWFT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGPAR(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_PHY_MF%YRPARAR%NGPAR+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGMVT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDFIELDS%YRGMV%YT1%NDIM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGMVT1S(YDGEOMETRY%YRDIM%NPROMA,YDFIELDS%YRGMV%YT1%NDIMS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKOZO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG*YDMODEL%YRML_PHY_G%YRDPHY%NVCLIS+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQS(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQICE(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQLI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQRAIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQSNOW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PATND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_DYN%YRDYNA%YYTTND%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDBBC(YDGEOMETRY%YRDIM%NPROMNH)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRDPHI(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGWT0(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGWT9(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGWS(YDGEOMETRY%YRDIM%NPROMNH)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGWFL(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGWFM(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNHXT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNHXT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQCHA0L(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQCHA0M(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEXTRA(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_PHY_G%YRDPHY%NVEXTRDYN)
TYPE(TYP_DDH)     ,INTENT(INOUT) :: YDDDH

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IFLAG, JGFL
!!! INTEGER(KIND=JPIM) :: INHFIELDS
LOGICAL :: LLSTR, LLGPXX
REAL(KIND=JPRB) :: ZEPS

REAL(KIND=JPRB), POINTER :: ZT0SP(:), ZT0SPD(:,:), ZT0SPL(:), ZT0SPM(:)
REAL(KIND=JPRB), POINTER :: ZT0SVD(:,:), ZT0T(:,:), ZT0U(:,:), ZT0V(:,:)
REAL(KIND=JPRB), POINTER :: ZT1NHX(:,:), ZT1SP(:), ZT1SPD(:,:), ZT1SVD(:,:)
REAL(KIND=JPRB), POINTER :: ZT1T(:,:), ZT1U(:,:), ZT1V(:,:), ZT9DIV(:,:)
REAL(KIND=JPRB), POINTER :: ZT9SP(:), ZT9SPD(:,:), ZT9SPL(:), ZT9SPM(:)
REAL(KIND=JPRB), POINTER :: ZT9SVD(:,:), ZT9T(:,:), ZT9U(:,:), ZT9V(:,:)
REAL(KIND=JPRB), POINTER :: ZVFALBF(:), ZVFEMISF(:), ZVFLSM(:), ZVFZ0F(:)
REAL(KIND=JPRB), POINTER :: ZVVD2(:), ZVVIVEG(:), ZVVZ0H(:), ZRRIC0(:)
REAL(KIND=JPRB), POINTER :: ZRRT0(:), ZRRW0(:), ZSBQ0(:,:), ZSBTL0(:,:)
REAL(KIND=JPRB), POINTER :: ZSBT0(:,:), ZSGA0(:,:), ZSGF0(:,:), ZSGR0(:,:)
REAL(KIND=JPRB), POINTER :: ZP1FORC(:,:), ZPQ(:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZPRE9L(YDGEOMETRY%YRDIM%NPROMA)            ! zonal comp of "grad prehyds" at t-dt
REAL(KIND=JPRB) :: ZPRE9M(YDGEOMETRY%YRDIM%NPROMA)            ! merid comp of "grad prehyds" at t-dt
REAL(KIND=JPRB) :: ZGDW0(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)     ! 'G dw' at full levels at t
REAL(KIND=JPRB) :: ZGDW9(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)     ! 'G dw' at full levels at t-dt
REAL(KIND=JPRB) :: ZGWHT0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)   ! 'Gw' at half levels at t
REAL(KIND=JPRB) :: ZGWHT9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)   ! 'Gw' at half levels at t-dt
REAL(KIND=JPRB) :: ZOROGLL(YDGEOMETRY%YRDIM%NPROMNH)          ! M**2 * YDOROG%OROGLL
REAL(KIND=JPRB) :: ZOROGMM(YDGEOMETRY%YRDIM%NPROMNH)          ! M**2 * YDOROG%OROGMM
REAL(KIND=JPRB) :: ZOROGLM(YDGEOMETRY%YRDIM%NPROMNH)          ! M**2 * YDOROG%OROGLM
REAL(KIND=JPRB) :: ZATND_Q(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! 1D model adiab Lagr tendency for "q"

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "cp_forcing.intfb.h"
#include "cp_forcing_ps.intfb.h"
#include "etenc.intfb.h"
#include "gpinislb.intfb.h"
#include "gpinozst.intfb.h"
#include "gpmpfc.intfb.h"
#include "gpnspng.intfb.h"
#include "gp_spv.intfb.h"
#include "gptf2.intfb.h"
#include "updsst.intfb.h"
#include "cpg_gp_hyd.intfb.h"
#include "cpg_gp_nhee.intfb.h"
#include "cpg_gp_nhqe.intfb.h"
#include "cpg_gp_sacc.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG_GP',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDOROG=>YDGEOMETRY%YROROG(KIBL), &
 & YDCVER=>YDGEOMETRY%YRCVER, &
 & YDPTRSLB2=>YDMODEL%YRML_DYN%YRPTRSLB2,YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR,YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH, &
 & YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY,YDDIMF=>YDMODEL%YRML_GCONF%YRDIMF,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
 & YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH,YGFL=>YDMODEL%YRML_GCONF%YGFL, &
 & YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY,YDDYN=>YDMODEL%YRML_DYN%YRDYN, &
 & YDDYNA=>YDMODEL%YRML_DYN%YRDYNA,YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL, &
 & YDPHLC=>YDMODEL%YRML_PHY_SLIN%YRPHLC)
ASSOCIATE(NDIM=>YGFL%NDIM, NDIM1=>YGFL%NDIM1, NUMFLDS=>YGFL%NUMFLDS, &
 & YCOMP=>YGFL%YCOMP, YFORC=>YGFL%YFORC, YQ=>YGFL%YQ, &
 & NPROMA=>YDDIM%NPROMA, &
 & NPROMNH=>YDDIM%NPROMNH, &
 & NFTHER=>YDDIMF%NFTHER, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NTOZ1D=>YDDPHY%NTOZ1D, NTOZ2D=>YDDPHY%NTOZ2D, NTOZ3D=>YDDPHY%NTOZ3D, &
 & NVCLIS=>YDDPHY%NVCLIS, NVEXTRDYN=>YDDPHY%NVEXTRDYN, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
 & LEO3CH=>YDEPHY%LEO3CH, &
 & RSTRET=>YDGEM%RSTRET, &
 & NDIMGMV=>YDFIELDS%YRGMV%NDIMGMV, NDIMGMVS=>YDFIELDS%YRGMV%NDIMGMVS, YT0=>YDFIELDS%YRGMV%YT0, &
 & YT1=>YDFIELDS%YRGMV%YT1, YT9=>YDFIELDS%YRGMV%YT9, &
 & LRSLDDH=>YDLDDH%LRSLDDH, &
 & MSLDDH_NHX=>YDMDDH%MSLDDH_NHX, MSLDDH_PD=>YDMDDH%MSLDDH_PD, &
 & MSLDDH_T=>YDMDDH%MSLDDH_T, MSLDDH_U=>YDMDDH%MSLDDH_U, &
 & MSLDDH_V=>YDMDDH%MSLDDH_V, MSLDDH_VD=>YDMDDH%MSLDDH_VD, &
 & LSPHLC=>YDPHLC%LSPHLC, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2, &
 & YSD_VF=>YDFIELDS%YRSURF%YSD_VF, YSD_VFD=>YDFIELDS%YRSURF%YSD_VFD, YSD_VV=>YDFIELDS%YRSURF%YSD_VV, &
 & YSD_VVD=>YDFIELDS%YRSURF%YSD_VVD, YSP_RR=>YDFIELDS%YRSURF%YSP_RR, YSP_RRD=>YDFIELDS%YRSURF%YSP_RRD, &
 & YSP_SB=>YDFIELDS%YRSURF%YSP_SB, YSP_SBD=>YDFIELDS%YRSURF%YSP_SBD, YSP_SG=>YDFIELDS%YRSURF%YSP_SG, &
 & YSP_SGD=>YDFIELDS%YRSURF%YSP_SGD, &
 & LSIMPH=>YDSIMPHL%LSIMPH, LMSE=>YDARPHY%LMSE, LMPA=>YDARPHY%LMPA, &
 & NGPAR=>YDPARAR%NGPAR)

CALL SC2PRG(1,YFORC(:)%MP,PGFL     ,ZP1FORC)
CALL SC2PRG(YQ%MP        ,PGFL     ,ZPQ)
CALL SC2PRG(YSD_VF%YALBF%MP,PSD_VF ,ZVFALBF)
CALL SC2PRG(YSD_VF%YEMISF%MP,PSD_VF,ZVFEMISF)
CALL SC2PRG(YSD_VF%YLSM%MP,PSD_VF  ,ZVFLSM)
CALL SC2PRG(YSD_VF%YZ0F%MP,PSD_VF  ,ZVFZ0F)
CALL SC2PRG(YSD_VV%YD2%MP,PSD_VV   ,ZVVD2)
CALL SC2PRG(YSD_VV%YIVEG%MP,PSD_VV ,ZVVIVEG)
CALL SC2PRG(YSD_VV%YZ0H%MP,PSD_VV  ,ZVVZ0H)
CALL SC2PRG(YSP_RR%YIC%MP0,PSP_RR  ,ZRRIC0)
CALL SC2PRG(YSP_RR%YT%MP0,PSP_RR   ,ZRRT0)
CALL SC2PRG(YSP_RR%YW%MP0,PSP_RR   ,ZRRW0)
CALL SC2PRG(YSP_SB%YQ%MP0,PSP_SB   ,ZSBQ0)
CALL SC2PRG(YSP_SB%YTL%MP0,PSP_SB  ,ZSBTL0)
CALL SC2PRG(YSP_SB%YT%MP0,PSP_SB   ,ZSBT0)
CALL SC2PRG(YSP_SG%YA%MP0,PSP_SG   ,ZSGA0)
CALL SC2PRG(YSP_SG%YF%MP0,PSP_SG   ,ZSGF0)
CALL SC2PRG(YSP_SG%YR%MP0,PSP_SG   ,ZSGR0)
CALL SC2PRG(YT0%MSP   ,PGMVS   ,ZT0SP)
CALL SC2PRG(YT0%MSPD  ,PGMV    ,ZT0SPD)
CALL SC2PRG(YT0%MSPL  ,PGMVS   ,ZT0SPL)
CALL SC2PRG(YT0%MSPM  ,PGMVS   ,ZT0SPM)
CALL SC2PRG(YT0%MSVD  ,PGMV    ,ZT0SVD)
CALL SC2PRG(YT0%MT    ,PGMV    ,ZT0T)
CALL SC2PRG(YT0%MU    ,PGMV    ,ZT0U)
CALL SC2PRG(YT0%MV    ,PGMV    ,ZT0V)
CALL SC2PRG(YT1%MNHX  ,PGMVT1  ,ZT1NHX)
CALL SC2PRG(YT1%MSP   ,PGMVT1S ,ZT1SP)
CALL SC2PRG(YT1%MSPD  ,PGMVT1  ,ZT1SPD)
CALL SC2PRG(YT1%MSVD  ,PGMVT1  ,ZT1SVD)
CALL SC2PRG(YT1%MT    ,PGMVT1  ,ZT1T)
CALL SC2PRG(YT1%MU    ,PGMVT1  ,ZT1U)
CALL SC2PRG(YT1%MV    ,PGMVT1  ,ZT1V)
CALL SC2PRG(YT9%MDIV  ,PGMV    ,ZT9DIV)
CALL SC2PRG(YT9%MSP   ,PGMVS   ,ZT9SP)
CALL SC2PRG(YT9%MSPD  ,PGMV    ,ZT9SPD)
CALL SC2PRG(YT9%MSPL  ,PGMVS   ,ZT9SPL)
CALL SC2PRG(YT9%MSPM  ,PGMVS   ,ZT9SPM)
CALL SC2PRG(YT9%MSVD  ,PGMV    ,ZT9SVD)
CALL SC2PRG(YT9%MT    ,PGMV    ,ZT9T)
CALL SC2PRG(YT9%MU    ,PGMV    ,ZT9U)
CALL SC2PRG(YT9%MV    ,PGMV    ,ZT9V)

!     ------------------------------------------------------------------

!*       0.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

ZEPS=100.0_JPRB*TINY(1.0_JPRB)
LLSTR=(ABS(RSTRET-1.0_JPRB)>ZEPS)

!     ------------------------------------------------------------------

!*       1.    INTERFACE TO GLOBAL ARRAYS/WORK FILES.
!              --------------------------------------


!*     1.1  TENDENCY COUPLING FOR SURFACE PRESSURE.

IF (YDMODEL%YRML_LBC%LTENC.AND..NOT.YDDYNA%LTWOTL) THEN
  CALL ABOR1('CPG_GP: ABOR1 CALLED: LTENC.AND..NOT.LTWOTL')
ENDIF
IF (NSTEP > 0.AND.YDMODEL%YRML_LBC%LTENC) THEN
  CALL ETENC(YDGEOMETRY,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_LBC,YDFIELDS%YRELBC_FIELDS,&
   & LD_DFISTEP,KENDC,KSTGLO,ZT0SP,&
   & ZT0SPL,ZT0SPM,&
   & ZT9SP,ZT9SPL,ZT9SPM)
ENDIF                 

!*     1.3   SPONGE AT THE TOP OF THE MODEL (ACADEMIC SIMULATIONS).

! new sponge for 2D and 3D models (grid-point GFL only).
IF (YDMODEL%YRML_DYN%YRSPNG%LNSPONGE) CALL GPNSPNG(YGFL,YDMODEL%YRML_DYN%YRSPNG,NPROMA,NFLEVG,KST,KEND,PGFL)

!*     1.4   PART 2 OF TIME FILTER.

IF(NCURRENT_ITER == 0) THEN
  CALL GPTF2(YDGEOMETRY,YDFIELDS%YRGMV,&
   ! --- INPUT ---------------------------------------------------------
   & YDMODEL%YRML_GCONF,YDDYN,YDDYNA,KST,KEND,LDLFSTEP,&
   ! --- INPUT for P.T0.; INPUT/OUTPUT for P.T9. -----------------------
   & PGMV,PGMVS,PGFL)
ENDIF

!*     1.5   MAP FACTOR AND SURFACE PRESSURE VARIABLES.

POROGL(KST:KEND)=YDOROG%OROGL(KST:KEND)*YDGSGEOM%GM(KST:KEND)
POROGM(KST:KEND)=YDOROG%OROGM(KST:KEND)*YDGSGEOM%GM(KST:KEND)
IF(YDDYNA%LNHDYN) THEN
  ZOROGLL(KST:KEND)=YDOROG%OROGLL(KST:KEND)*YDGSGEOM%GM(KST:KEND)**2
  ZOROGMM(KST:KEND)=YDOROG%OROGMM(KST:KEND)*YDGSGEOM%GM(KST:KEND)**2
  ZOROGLM(KST:KEND)=YDOROG%OROGLM(KST:KEND)*YDGSGEOM%GM(KST:KEND)**2
ENDIF

IF(LLSTR.AND.(.NOT.YDDYNA%LPGFSF)) THEN
  IFLAG=0
  CALL GPMPFC(YDFIELDS%YRGMV,YDMODEL%YRML_GCONF,YDDYN,YDDYNA,NPROMA,NFLEVG,KST,KEND,IFLAG,YDGSGEOM%GM,PGMV,PGMVS,PGFL)
ENDIF

CALL GP_SPV(YDGEOMETRY,YDDYN,YDDYNA,YDSIMPHL,.FALSE.,KST,KEND,&
 & ZT0SP,ZT0SPL,ZT0SPM,&
 & ZT9SP,ZT9SPL,ZT9SPM,&
 & PRE0,PRE0L,PRE0M,&
 & PRE9,ZPRE9L,ZPRE9M)

IF (LSFORC .AND. LSPS_FRC) THEN
  CALL CP_FORCING_PS(YDGEOMETRY,KST,KEND,PRE0)
ENDIF

IF (NCURRENT_ITER == 0) THEN 

!*     1.6   SURFACE VARIABLES.

  IF (LDLDIAB.OR.LSPHLC.OR.LSIMPH) THEN
    PQS(KST:KEND)=0.0_JPRB
    IF(LDLFSTEP .AND. .NOT. YDDYNA%LTWOTL) THEN
      CALL GPPOPER(YDDYN,'SET9TO0',YDFIELDS%YRSURF,KBL=KBL)
    ENDIF
  ENDIF

  ! for AROME   
  IF (LMSE.AND.NGPAR/=0) THEN
    PGPAR(1:NPROMA,1:NGPAR) = GPARBUF (:, :, 1+(KSTGLO-1)/YDGEOMETRY%YRDIM%NPROMA)
  ENDIF

!*     1.8   OZON PHYSICO-CHEMICAL PROPERTIES AND SUBGRID SURFACE TEMPERATURE.

  IF (LDLDIAB) THEN
    IF(NTOZ3D > 0.OR.NTOZ2D > 0.OR.NTOZ1D > 0.AND.(.NOT.LEO3CH)) THEN
      CALL GPINOZST(YDGEOMETRY,YDMODEL%YRML_CHEM%YROZO,YDDPHY,KST,KEND,KSTGLO,PKOZO)
    ENDIF
  ENDIF

!*     1.9  NUDGING.

 IF(LNUDG.AND.(.NOT.LMSE)) THEN
    CALL UPDSST(YDMODEL%YRML_AOC%YRMCC,YDMODEL%YRML_PHY_MF%YRPHY1,NPROMA,KST,KEND,NFNUDG,YSP_SBD%NLEVS,KBL,&
     & XPNUDG,&
     & ZRRT0,&
     & ZSBT0,ZRRW0,&
     & ZSBQ0,ZSGF0,&
     & ZSGA0,ZSGR0,&
     & ZRRIC0,ZSBTL0,&
     & ZVFLSM,ZVFALBF,&
     & ZVFEMISF,ZVFZ0F,&
     & ZVVZ0H,ZVVD2,&
     & ZVVIVEG)
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       2.    INTERFACE TO GLOBAL ARRAYS/WORK FILES.
!              --------------------------------------

!     ------------------------------------------------------------------

!*       3.    INITIALISE AUXILIARY VARIABLES.
!              -------------------------------

IF (YDDYNA%LNHQE) THEN
  ! * NHQE model:
  LLGPXX=.NOT.(YDDYNA%LSLAG.AND.(YDDYNA%NVDVAR == 4).AND.(YDDYNA%ND4SYS==2))
  CALL CPG_GP_NHQE(YDGEOMETRY,YDMODEL,YDFIELDS%YRGMV,&
   !---------------------------------------------------------------------
   ! - INPUT .
   & KST,KEND,KSTGLO,LLGPXX,LDLDIAB,LMPA,&
   & KIBL,YDOROG%OROG,POROGL,POROGM,ZOROGLL,ZOROGMM,ZOROGLM,PRE0L,PRE0M,ZPRE9L,ZPRE9M,&
   !---------------------------------------------------------------------
   ! - INPUT/OUTPUT .
   & PGFL,PGMV,PRE0,PRE9,&
   !---------------------------------------------------------------------
   ! - OUTPUT .
   & PRE0F,PNHPRE0F,PNHPRE0H,PXYB0,PUVH0,&
   & PPGTERM%PHI0H,PHIF0,PPGTERM%PHI0FL,PPGTERM%PHI0FM, &
   & PRCP0,PCTY0,PGWFT0,ZGWHT0,PGWT0,PGWS,PGWFL,PGWFM,ZGDW0,PKENE0,&
   & PPGTERM%RT0L,PPGTERM%RT0M,&
   & PRE9F,PNHPRE9F,PNHPRE9H,PXYB9,PHI9,PHIF9,PRCP9,PGWFT9,ZGWHT9,PGWT9,ZGDW9,&
   & PATND,PDBBC,PRDPHI,PNHXT0,PNHXT9,PQCHA0L,PQCHA0M)
ELSEIF (YDDYNA%LNHEE) THEN
  ! * NHEE model:
  LLGPXX=.NOT.(YDDYNA%LSLAG.AND.(YDDYNA%NVDVAR == 4).AND.(YDDYNA%ND4SYS==2))
  CALL CPG_GP_NHEE(YDGEOMETRY,YDMODEL,YDFIELDS%YRGMV,&
   !---------------------------------------------------------------------
   ! - INPUT .
   & KST,KEND,KSTGLO,LLGPXX,LDLDIAB,LMPA,&
   & KIBL,YDOROG%OROG,POROGL,POROGM,ZOROGLL,ZOROGMM,ZOROGLM,PRE0L,PRE0M,&
   !---------------------------------------------------------------------
   ! - INPUT/OUTPUT .
   & PGFL,PGMV,PRE0,PRE9,&
   !---------------------------------------------------------------------
   ! - OUTPUT .
   & PRE0F,PNHPRE0F,PNHPRE0H,PXYB0,PUVH0,&
   & PPGTERM%PHI0H,PHIF0,PPGTERM%PHI0FL,PPGTERM%PHI0FM, &
   & PRCP0,PCTY0,PGWFT0,ZGWHT0,PGWT0,PGWS,PGWFL,PGWFM,ZGDW0,PKENE0,&
   & PPGTERM%RT0L,PPGTERM%RT0M,&
   & PRE9F,PNHPRE9F,PNHPRE9H,PXYB9,PHI9,PHIF9,PRCP9,PGWFT9,ZGWHT9,PGWT9,ZGDW9,&
   & PATND,PDBBC,PRDPHI,PNHXT0,PNHXT9,PQCHA0L,PQCHA0M)

ELSEIF (YDDYNA%LSACC) THEN
  ! * Shallow-atmosphere model with complete Coriolis of Tort & Dubos (2013)
  ! * LLPGX=t is needed because we need [gw] term 
  LLGPXX=.TRUE.
  CALL CPG_GP_SACC(YDGEOMETRY,YDMODEL,YDFIELDS%YRGMV,&
   !---------------------------------------------------------------------
   ! - INPUT .
   & KST,KEND,KSTGLO,LLGPXX,LDLDIAB,LMPA,&
   & KIBL,YDOROG%OROG,POROGL,POROGM,PRE0L,PRE0M,&
   !---------------------------------------------------------------------
   ! - INPUT/OUTPUT .
   & PGFL,PGMV,PRE0,PRE9,&
   !---------------------------------------------------------------------
   ! - OUTPUT .
   & PRE0F,PXYB0,PUVH0,&
   & PPGTERM%PHI0H,PHIF0,PPGTERM%PHI0FL,PPGTERM%PHI0FM,&
   & PRCP0,PCTY0,PGWFT0,ZGWHT0,PKENE0,&
   & PPGTERM%RT0L,PPGTERM%RT0M,&
   & PRE9F,PXYB9,PHI9,PHIF9,PRCP9,PGWFT9,ZGWHT9,&
   & PATND,PNHXT0,PNHXT9)
ELSE
  ! * Hydrostatic model:
  IF(YDCVER%LVERTFE.AND..NOT.YDCVER%LVFE_ECMWF)THEN
    LLGPXX=.FALSE.
  ELSE
    LLGPXX=.TRUE.
  ENDIF
  CALL CPG_GP_HYD(YDGEOMETRY,YDMODEL,YDFIELDS%YRGMV,&
   !--------------------------------------------------------------------
   ! - INPUT .
   & KST,KEND,KSTGLO,LLGPXX,LDLDIAB,LMPA,&
   & KIBL,YDOROG%OROG,POROGL,POROGM,PRE0L,PRE0M,&
   !---------------------------------------------------------------------
   ! - INPUT/OUTPUT .
   & PPGTERM,PGFL,PGMV,PRE0,PRE9,&
   !---------------------------------------------------------------------
   ! - OUTPUT .
   & PRE0F,PXYB0,PUVH0,&
   & PHIF0,PRCP0,PCTY0,PGWFT0,ZGWHT0,PKENE0,&
   & PRE9F,PXYB9,PHI9,PHIF9,PRCP9,PGWFT9,ZGWHT9,&
   & PATND,PNHXT0,PNHXT9) 
ENDIF

!     ------------------------------------------------------------------

!*       4.    TIME DIMENSION.
!              ---------------

CALL GPINISLB(&
 & YDGEOMETRY,YGFL,YDDYN,YDDYNA,YDPTRSLB2,KST,KEND,PTE,.FALSE.,&
 & PGFL,&
 & ZT9U,ZT9V,ZT9T,&
 & ZT9SPD,ZT9SVD,PNHXT9,&
 & ZT9SP,&
 & PCTY0(1,1,YYTCTY0%M_VVEL),PRE0F,&
 & ZT1U,ZT1V,ZT1T,PGFLT1,&
 & ZT1SPD,ZT1SVD,ZT1NHX,&
 & ZT1SP,&
 & PB2,PQICE,PQLI,PQRAIN,PQSNOW,&
 & PGWFT0=PGWFT0,PGDW0=ZGDW0,PGWS0=ZGWHT0(1,NFLEVG),&
 & PGWFT9=PGWFT0,PGDW9=ZGDW9,PGWS9=ZGWHT9(1,NFLEVG))

IF (NCURRENT_ITER == 0) THEN
  IF (LDLDIAB.OR.LSIMPH) THEN
    CALL GPPOPER(YDDYN,'SET1TO9',YDFIELDS%YRSURF,KBL=KBL)
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       5.    TENDENCIES FOR 1D MODEL.
!              ------------------------

! Add large scale forcing for SCUM (Single Column Unified Model)
! - incrementation of PATND for GMV variables.
! - fills ZATND_Q for humidity.

IF (LSFORC) THEN
  IF (NCURRENT_ITER > 0 .AND. YDDYNA%LPC_FULL) THEN
    ! this call to CP_FORCING has not been well checked for the corrector step,
    !  but the corrector step cannot work for GFL variables.
    CALL ABOR1('CPG_GP: part 5')
  ENDIF

  IF (YQ%MP /= YQ%MP1) THEN
    CALL ABOR1('CPG_GP: CP_FORCING NEEDS DOUBLE CHECK ON Q POINTER. REK')
  ENDIF
  CALL CP_FORCING(YDGEOMETRY,YDMODEL,KST,KEND,PDT,PRE0,PRE0F,ZP1FORC,&
   & ZT0U,ZT0V,&
   & ZT0T,ZPQ,&
   ! most likely not correct : & PATND,ZATND_Q,PGFLT1(1,1,YQ%MP),YDDDH)
   & PATND,ZATND_Q,PGFLT1(:,:,YQ%MP1),YDDDH)
ENDIF

!     ------------------------------------------------------------------

!*       6.    STORE T9 ARRAYS FOR NH TRAJECTORY.
!              ----------------------------------

!!! Obsolete way to manage NH trajectory, to recode properly if developpement of some TL/AD NH code.
!!! IF (LNHDYN.AND.LTRAJNH) THEN
!!!   INHFIELDS=NG3NH95*NFLEVG+1
!!!   CALL WRNHTRAJM(YDGEOMETRY,YDMODEL%YRML_DYN%YRTNH,NSTEP,KST,KEND,KSTGLO,INHFIELDS,&
!!!    & ZT0SPD,ZT9SVD,&
!!!    & ZT9DIV,ZT9U,ZT9V,ZT9T,&
!!!    & ZT9SP)
!!! ENDIF

!     ------------------------------------------------------------------

!*       7.    INITIALZE SL DDH TENDENCY ARRAYS.
!              ----------------------------------

IF (LRSLDDH.AND.YDDYNA%LSLAG) THEN

  IF (YDDYNA%LTWOTL) THEN
    PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_U)=-ZT0U(KST:KEND,1:NFLEVG)
    PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_V)=-ZT0V(KST:KEND,1:NFLEVG)
    IF (NFTHER > 0) PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_T)=-ZT0T(KST:KEND,1:NFLEVG)
    IF (YDDYNA%LNHDYN)THEN !!! ?????? ky: LNHDYN or LNHEE there?
      PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_PD)=-ZT0SPD(KST:KEND,1:NFLEVG)
    ENDIF
    IF (YDDYNA%LNHDYN)THEN
      PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_VD)=-ZT0SVD(KST:KEND,1:NFLEVG)
      IF (YDDYNA%LNHX) PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_NHX)=-PNHXT0(KST:KEND,1:NFLEVG)
    ENDIF
    DO JGFL=1,NUMFLDS
      IF (YCOMP(JGFL)%LADV) THEN
        PGFLTNDSL(KST:KEND,1:NFLEVG,JGFL)=-PGFL(KST:KEND,1:NFLEVG,YCOMP(JGFL)%MP)
      ENDIF    
    ENDDO
  ELSE
    PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_U)=-ZT9U(KST:KEND,1:NFLEVG)
    PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_V)=-ZT9V(KST:KEND,1:NFLEVG)
    IF (NFTHER > 0) PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_T)=-ZT9T(KST:KEND,1:NFLEVG)
    IF (YDDYNA%LNHDYN)THEN !!! ?????? ky: LNHDYN or LNHEE there?
      PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_PD)=-ZT9SPD(KST:KEND,1:NFLEVG)
    ENDIF
    IF (YDDYNA%LNHDYN)THEN
      PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_VD)=-ZT9SVD(KST:KEND,1:NFLEVG)
      IF (YDDYNA%LNHX) PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_NHX)=-PNHXT9(KST:KEND,1:NFLEVG)
    ENDIF
    DO JGFL=1,NUMFLDS
      IF (YCOMP(JGFL)%LADV) THEN
        PGFLTNDSL(KST:KEND,1:NFLEVG,JGFL)=-PGFL(KST:KEND,1:NFLEVG,YCOMP(JGFL)%MP9)
      ENDIF
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_GP',1,ZHOOK_HANDLE)
END SUBROUTINE CPG_GP
