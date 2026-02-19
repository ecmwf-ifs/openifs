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

MODULE YOMDYNA

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE INTDYN_MOD,ONLY : TTND,TGMVT,TGFLT

IMPLICIT NONE

SAVE

TYPE TDYNA
! ----------------------------------------------------------------------

!=========== VARIABLES FOR DYNAMICS: FIRST PART ===============================
! We put there not geometry-dependent variables for dynamics.
!==============================================================================
! High level switches (from YOMCT0)
! ----- type of equations:
! LNHQE     : T if NHQE model.
! LNHEE     : T if NHEE model.
! LNHDYN    : .T. if 3-D non-hydrostatic dynamics is active (NHEE or NHQE)
LOGICAL :: LNHEE
LOGICAL :: LNHQE
LOGICAL :: LNHDYN

! ----- type of equations:
! LSACC     : T if hydrostatic shallow-atmosphere complete-Coriolis model: see Tort & Dubos (2013)
LOGICAL :: LSACC
! ----- control of variables which are transformed in spectral space:
! LSPRT     : .T.: if R*T/Rd "virtual temperature" as spectral variable
LOGICAL :: LSPRT
! LTWOTL    : .TRUE. if two-time-level semi-Lagrangian scheme will be used
LOGICAL :: LTWOTL
! ----- advection scheme, semi-Lagrangian scheme:
! LSLAG     : .TRUE. = semi-lagrangian on
LOGICAL :: LSLAG

! ------ NH model ------------------------------------------------------

! NPDVAR : switch for type of variable used for pressure departure
!         2: q_hat = ln(p/pi)

! NVDVAR : switch for type of variable used for pseudo vertical divergence
!         3: d3 = -g (p/(m.Rd.T)) d_parc w/d_parc eta
!         4: d4 = d3 + X = d3 + (p/m.R.T) nabla_phi d_parc V/d_parc eta

! ND4SYS : switch for the way of treatment of term NHX in the NH d4 equation.
!          ND4SYS=1: all contributions of D(NHX)/Dt are treated at the level
!                    of CPG_DYN+CALL_SL, and X is updated via GPXX in CPG_GP.
!          ND4SYS=2 (SL only): only the advective terms are treated at the
!                    level of CPG_DYN+CALL_SL; additional contributions
!                    are done in CPGLAG, and X(t+dt) is updated in CPGLAG.

! LNHX   : variable NHX (grid-point and spectral) is needed at t.
! LNHXDER: horizontal derivatives of variable NHX are needed at t.

! LGWADV: .T. => "vwv" prognostic variable in NH model is "g*w"
!                (but the variable which is transformed into spectral
!                space is "d_nvdvar").
!         .F. => "vwv" prognostic variable in NH model is "d_nvdvar".

! NGWADVSI: when LGWADV, alternate treatments for linear terms.
!           NGWADVSI=1: linear terms set to 0 in LATTEX, LATTES for all eqns;
!            linear terms evaluated at O and linear terms evaluated at F
!            added to the RHS of equations in GPENDTR.
!           NGWADVSI=2 (relevant for LSETTLS but not for LNESC):
!            linear terms set to 0 in LATTEX, LATTES for all eqns;
!            linear terms evaluated at O added to the RHS of equations in
!            LAPINEB (after the call to GNHGW2SVD); linear terms evaluated
!            at F added to the RHS of equations in GPENDTR.

! LRDBBC: .T. if S.-L. diagnostic BBC active for NH (if LSLAG=.T. only)
!         .F. if eulerian development of NH  BBC in LSLAG.

! LSI_NHEE : Helmholtz eqn with DIV as unknown in the NHEE model.

! LNHQE_C2: T/F: constraint C2 is ensured/not ensured in NHQE model.
! LNHQE_SIHYD: T: use HYD SI system in NHQE model for validation, Qcha is kept unchanged from one timestep to the following one.
! LNHEE_SVDLAPL_FIRST: controls order of calculations in RHS of (gw) and vertical divergence equation in NHEE model.
!  T: first compute "LAPL [(pre-prehyd)/prehyd]", then compute RHS of (gw) equation.
!  F: first compute RHS of (gw) equation, then compute "LAPL [(pre-prehyd)/prehyd]".

INTEGER(KIND=JPIM) :: NPDVAR
INTEGER(KIND=JPIM) :: NVDVAR
INTEGER(KIND=JPIM) :: ND4SYS
LOGICAL :: LNHX
LOGICAL :: LNHXDER
LOGICAL :: LGWADV
INTEGER(KIND=JPIM) :: NGWADVSI
LOGICAL :: LRDBBC
LOGICAL :: LSI_NHEE
LOGICAL :: LNHQE_C2
LOGICAL :: LNHQE_SIHYD
LOGICAL :: LNHEE_SVDLAPL_FIRST


! ------ SLHD diffusion ------------------------------------------------

!          SLHD (= horizontal diffusion using damping properties of some
!          semi-Lagrangian interpolatores) keys
! LSLHD     : internal model switch for semi-Lagrangian diffusion computation
! LSLHD_W   : switch for SLHD of horizontal flow.
! LSLHD_T   : switch for SLHD of temperature.
! LSLHD_SPD : switch for SLHD of (NH) pressure departure vs geopotential var.
! LSLHD_SVD : switch for SLHD of (NH) vertical divergence vs 'gw' var.
! LSLHD_GFL : switch for SLHD of GFL var (at least one GFL var has SLHD).
! LSLHD_OLD : use old SLHD interpolator (cubic Lagrange mixed with linear)
! LSLHD_STATIC : do not diagnose kappa from horizontal flow deformation,
!                use static value instead:
!                  SLHDKMIN for undiffused fields
!                  SLHDKMAX for diffused fields
! LSLHDQUAD : internal model switch indicating need to precompute quadratic weights
! SLHDKMIN  : minimum value for the kappa function
! SLHDKMAX  : maximum value for the kappa function
! SLHDKREF  : Reference SLHDKMIN to which SLHD interpolation is relaxed
!               in areas with increased diffusivity (0 stands for Lagrangian cubic)
! SLHDEPSH  : dimensionless strength of horizontal Laplacian smoothing
! SLHDEPSV  : dimensionless strength of vertical Laplacian smoothing

LOGICAL :: LSLHD
LOGICAL :: LSLHD_W
LOGICAL :: LSLHD_T
LOGICAL :: LSLHD_SPD
LOGICAL :: LSLHD_SVD
LOGICAL :: LSLHD_GFL
LOGICAL :: LSLHD_OLD
LOGICAL :: LSLHD_STATIC
LOGICAL :: LSLHDQUAD
REAL(KIND=JPRB) :: SLHDKMIN
REAL(KIND=JPRB) :: SLHDKMAX
REAL(KIND=JPRB) :: SLHDKREF
REAL(KIND=JPRB) :: SLHDEPSH
REAL(KIND=JPRB) :: SLHDEPSV

! ------- SLVF scheme -------------------------------------------------------
!   SLVF = semi-Lagrangian vertical filter 
!          filtering of spurious gravity waves generated due to insufficient vertical resolution

! LSLVF     : main key to activate SLVF on temperature
! SLVFEPS   : dimensionless strength of vertical Laplacian smoothing
! P_SLVF_BOT : pressure level above which SLVF becomes active (half level)
!             (0. means SLVF is set to be active in pure pressure levels)
! SLVFMIN   : value*SLVFEPS to which the SLVF transitios to in the upper stratosphere
! SLVFTR    : the steepness of the transition from SLVFMIN to 1
! P_SLVF_TR : pressure at which SLVF transitions from 1*SLVFEPS to SLVFMIN*SLVFEPS 

LOGICAL :: LSLVF
REAL(KIND=JPRB) :: SLVFEPS
REAL(KIND=JPRB) :: P_SLVF_BOT
REAL(KIND=JPRB) :: SLVFMIN
REAL(KIND=JPRB) :: SLVFTR
REAL(KIND=JPRB) :: P_SLVF_TR

! ------ Other diffusive processes ------------------------------------------

! LGRADSP    : special switch for de-aliasing the pressure gradient term
! LPGFSF     : use de-aliasing filter on the pressure gradient term by spectral fit
! LPGREUSE   : re-use precomputed PG quantities when LPGFSF=true

LOGICAL :: LGRADSP
LOGICAL :: LPGFSF
LOGICAL :: LPGREUSE

! ------ 3D turbulence ------------------------------------------------------

! L3DTURB   : main key to activate 3D turbulence

LOGICAL :: L3DTURB

! ------ Dynamics diagnostics -----------------------------------------------

! LSLDIA    : switch on semi-lagrangian dynamics diagnostics
! LRPRSLTRJ : more detailed printings in LARMES/ELARMES (SL displacement).

LOGICAL :: LSLDIA
LOGICAL :: LRPRSLTRJ

! ------ Semi-Lagrangian scheme ---------------------------------------------

! LRALTVDISP : alternate way to compute vertical displacement (semi-Lagrangian advection),
!              in order to avoid trajectories going out of the atmosphere.
! LVSPLIP    : .T. if vertical spline cubic SL interpolations.
! LRHSVWENO  : .T. if vertical quintic interpolation is used in SL scheme.
! LCOMAD     : swith for COMAD (will be T if any variable needs COMAD interpolation).
! LCOMADH/LCOMADV : COMAD in the horizontal/vertical directions of interpolation
!                   same choice common for all COMAD interpolation
! LCOMAD_W   : COMAD for horizontal flow.
! LCOMAD_T   : COMAD for temperature.
! LCOMAD_SPD : COMAD for (NH) pressure departure vs geopotential var.
! LCOMAD_SVD : COMAD for (NH) vertical divergence vs 'gw' var.
! LCOMAD_SP  : COMAD for surface pressure.
! LCOMAD_GFL : COMAD for GFL (must be T if one of GFL attribut LCOMAD=T).
! LNESCT     : non-extrapolating horizontal displacement in SL2TL.
! LNESCV     : non-extrapolating vertical displacement in SL2TL.
! LNESC      : non-extrapolating RHS in SL2TL.
! LSETTLST   : stable extrapolating horizontal displacement in SL2TL.
! LSETTLSV   : stable extrapolating vertical displacement in SL2TL.
! LSETTLS    : stable extrapolating RHS in SL2TL.
! LELTRA     : alternate ``elegant" calculation of horizontal and vertical displacement in SL2TL.
! LSLINLC2   : separate linear terms from non-linear terms in continuity equation (case (LGWADV,LSETTLS)=(T,T)).
! LSLINL     : separate linear terms from non-linear terms in other equations (case (LGWADV,LSETTLS)=(T,T)).
!
! ------ Higher order interpolation for the SL trajectory research ---------
! LHOISLT    : activates High Order Interpolation of the SL Trajectory
! HOISLTV    : diffusiviry parameter for High Order Interpolation of the SL Trajectory along Vertical
! HOISLTH    : diffusiviry parameter for High Order Interpolation of the SL Trajectory along Horizontal
! LSLTVWENO  : WENO (weighted essentially non-oscillatory) interpolation used along Vertical for the SL 
!              Trajectory research

LOGICAL :: LRALTVDISP
LOGICAL :: LVSPLIP
LOGICAL :: LRHSVWENO
LOGICAL :: LCOMAD
LOGICAL :: LCOMADH
LOGICAL :: LCOMADV
LOGICAL :: LCOMAD_W
LOGICAL :: LCOMAD_T
LOGICAL :: LCOMAD_SPD
LOGICAL :: LCOMAD_SVD
LOGICAL :: LCOMAD_SP
LOGICAL :: LCOMAD_GFL
LOGICAL :: LNESCT
LOGICAL :: LNESCV
LOGICAL :: LNESC
LOGICAL :: LSETTLST
LOGICAL :: LSETTLSV
LOGICAL :: LSETTLS
LOGICAL :: LELTRA
LOGICAL :: LSLINLC2
LOGICAL :: LSLINL
LOGICAL :: LHOISLT
REAL(KIND=JPRB) :: HOISLTV
REAL(KIND=JPRB) :: HOISLTH
LOGICAL :: LSLTVWENO

! ------ Vertical discretisation --------------------------------------------

! LRUBC   : .T. if radiational upper boundary condition

LOGICAL :: LRUBC

! ------ PC (ICI) schemes ---------------------------------------------------

! LPC_FULL  : full PC scheme switch (with reiterations of trajectories)
! LPC_CHEAP : 'cheap PC scheme': when LPC_FULL=T and semi-Lagrangian advection,
!             the PC update is not done on the calculation of the SL trajectory.
! LPC_CHEAP2: alternate way to do "cheaper" FULL PC with semi-Lagragian advection.
!             - trajectory is recomputed at each iteration
!             - for GMV and 3D interpolations (LAITRE_GMV), linear interpolations if ncurrent_iter < nsiter
!             - for GMV and 3D interpolations (LAITRE_GMV), high-order interpolations only for ncurrent_iter = nsiter
!             - high order interpolations are kept for any "ncurrent_iter" for GFL and calls to LAIDDI.
! LMIXETTLS : MIXed Extrapolation for Two Time Level Scheme,
!             for the mixed NESC/SETTLS scheme depending on the weights calc.in LATTE_NL.
! LMIXETTLS_PRINT : print mixed NESC/SETTLS statistics in all levels to listing
! RMIXNL_TRH: treshold for mixed NESC/SETTLS scheme

LOGICAL :: LPC_FULL
LOGICAL :: LPC_CHEAP
LOGICAL :: LPC_CHEAP2
LOGICAL :: LMIXETTLS
LOGICAL :: LMIXETTLS_PRINT
REAL(KIND=JPRB) :: RMIXNL_TRH

!     ------------------------------------------------------------------
! L_RDRY_VD      : .T.: define vertical divergence with R_dry
!                  .F.: define vertical divergence with R_moist
LOGICAL :: L_RDRY_VD

TYPE(TTND) :: YYTTND
TYPE(TGMVT) :: YYTGMVT95
TYPE(TGFLT) :: YYTGFLT95
END TYPE TDYNA
! ----------------------------------------------------------------------
END MODULE YOMDYNA
