! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMDYN

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

! JPSLDIMK  : maximum possible value for NSLDIMK
INTEGER(KIND=JPIM),PARAMETER :: JPSLDIMK=4

TYPE :: TDYN

!     -------------------------------------------------------------------------

!*    Control variables for the DYNAMICS
!     We put there geometry-dependent variables used in dynamics.
!     Values may be different for the different models run under the OOPS layer.

!====== ASSELIN TIMEFILTERING =================================================

! REPS1   : timefiltering constant applied to t-1
! REPS2   : timefiltering constant applied to t+1
! REPSM1  : timefiltering constant applied to t-1 (moisture vars.)
! REPSM2  : timefiltering constant applied to t+1 (moisture vars.)
! REPSP1  : timefiltering constant applied to t-1 for all surface fields

REAL(KIND=JPRB) :: REPS1
REAL(KIND=JPRB) :: REPS2
REAL(KIND=JPRB) :: REPSM1
REAL(KIND=JPRB) :: REPSM2
REAL(KIND=JPRB) :: REPSP1

!====== MAIN HORIZONTAL DIFFUSION SCHEME ======================================

! * CHARACTERISTIC TIMES:
! HDIRVOR  : for diffusion of vorticity.
! HDIRDIV  : for diffusion of divergence.
! HDIRT    : for diffusion of temperature.
! HDIRQ    : for diffusion of humidity.
! HDIRO3   : for diffusion of ozone.
! HDIRPD   : for diffusion of pressure departure (non hydrostatic).
! HDIRVD   : for diffusion of vertical divergence (non hydrostatic).
! HDIRSP   : for diffusion of surface pressure.

! * REVERSE OF CHARACTERISTIC TIMES:
! HRDIRVOR  : for diffusion of vorticity.
! HRDIRDIV  : for diffusion of divergence.
! HRDIRT    : for diffusion of temperature.
! HRDIRQ    : for diffusion of humidity.
! HRDIRO3   : for diffusion of ozone.
! HRDIRPD   : for diffusion of pressure departure (non hydrostatic).
! HRDIRVD   : for diffusion of vertical divergence (non hydrostatic).
! HRDIRSP   : for diffusion of surface pressure.

! RRDXTAU  : overall intensity of HD 
! RDAMPVOR : local enhancing coefficient for diffusion of vorticity.
! RDAMPDIV : local enhancing coefficient for diffusion of divergence.
! RDAMPT   : local enhancing coefficient for diffusion of temperature.
! RDAMPQ   : local enhancing coefficient for diffusion of humidity.
! RDAMPO3  : local enhancing coefficient for diffusion of ozone.
! RDAMPPD  : local enhancing coefficient for diffusion of pressure departure.
! RDAMPVD  : local enhancing coefficient for diffusion of vertical divergence.
! RDAMPSP  : local enhancing coefficient for diffusion of surface pressure.
! LNEWHD   : only for ECMWF: "new" or "historical" values of HD set-up

! Coefficients RDI[X] generally write HRDIR[X]*g(l)*f(n,N(l),n0(X),x0,r) where n0
!  is generally 0 excepted for vorticity where n0 may be 2 in some cases.
! REXPDH   : order "r" of the diffusion (exponent for the wavenumber dependency).
! FRANDH   : threshold "x0" for the wavenumber dependency.
! SLEVDH   : first threshold for the pressure dependency scaled by VP00 used in function g(l).
! SLEVDH1  : first threshold for the pressure dependency scaled by VP00 used in function N(l).
! SLEVDH2  : second threshold for the pressure dependency scaled by VP00 used in function N(l).
! SLEVDH3  : third threshold for the pressure dependency scaled by VP00 used in function g(l)
!            (used to bound the vertical increase of diffusion in the upper stratosphere).
! NSREFDH  : threshold for the truncation dependency used in function N(l).
! RATIO_HDI_TOP: coefficient of additional diffusion near the top, used in function N(l) (used instead of NSREFDH)
!            When both NSREFDH and RATIO_HDI_TOP are given in NAMDYN, RATIO_HDI_TOP has more priority

! NPROFILEHD : type of vertical profile (function g(l)) used in the horizontal diffusion.
!              1: ECMWF-type profile.
!              2: flexible profile using RPROFHDBT,RPROFHDTP,RPROFHDMX,RPROFHDEX.
!              3: 1/prehyd profile (cleanly rewritten version of 4).
!              4: old (and less flexible) version of 3.

! RPROFHDBT to RPROFHDEX are used only if NPROFILEHD=2:
!  For NPROFILEHD=2:
!  * if "standard pressure > RPROFHDBT" same function as for NPROFILEHD=3
!  * if "standard pressure in [RPROFHDTP,RPROFHDBT]" a specific function is applied,
!    the inflexion point of which is controlled by the exponent RPROFHDEX.
!  * if "standard pressure < RPROFHDTP" function is equal to its maximum equal to RPROFHDMX.
!  a specific profile is used between pressure levels RPROFHDBT and RPROFHDTP

! LRDISPE_EC : horizontal function used in the horizontal diffusion.
!  LRDISPE_EC=T: ECMWF type function using PDISPEL, PDISPEE, PDISPE and PDISPEX (not implemented in LAM models).
!  LRDISPE_EC=F: MF type function using PDISPE and PDISPVOR (global model) or BDISPE (LAM model).

! * LEVEL AND WAVENUMBER DEPENDENT INVERSE CHARACTERISTIC TIMES:
! RDIVOR   : for diffusion of vorticity.
! RDIDIV   : for diffusion of divergence.
! RDITG    : for diffusion of temperature.
! RDIGFL   : for diffusion of GFL vars.
! RDIPD    : for diffusion of pressure departure (NH).
! RDIVD    : for diffusion of vertical divergence (NH).
! RDISP    : for diffusion of surface pressure.

! RDHI     : main horizontal diffusion operator used for stretched ARPEGE.

! LSTRHD   : .T.: main horizontal diffusion operator adapted to stretched ARP.
! HDTIME_STRHD: TDT (if not, the main horizontal diffusion operator
!            used for stretched ARPEGE is recomputed).

! LTOP_VOR : optional sponge for vorticity at model top 
!            introduced to remove "solitons" where LLMESO kills the DIV
! NTOP_VOR_TRUNC : cut-off wavenumber for 4th order VOR sponge 
                ! (KMAX in formulae)
! NTOP_VOR_BOT : bottom level of extra sponge for VOR 
                ! (extra sponge linear between NTOP_VOR_BOT and 1) 
! KLSPONGE     : vertical level from which to begin strong DIV damping
! RBEGSPONGE : from which Z_PDILEV to begins sponge, the higher the number the shallower the sponge
! LM0DAMP    : if true, damps the m=0 mode (zonal-mean)


REAL(KIND=JPRB) :: HDIRVOR
REAL(KIND=JPRB) :: HDIRDIV
REAL(KIND=JPRB) :: HDIRT
REAL(KIND=JPRB) :: HDIRQ
REAL(KIND=JPRB) :: HDIRO3
REAL(KIND=JPRB) :: HDIRPD
REAL(KIND=JPRB) :: HDIRVD
REAL(KIND=JPRB) :: HDIRSP
REAL(KIND=JPRB) :: HRDIRVOR
REAL(KIND=JPRB) :: HRDIRDIV
REAL(KIND=JPRB) :: HRDIRT
REAL(KIND=JPRB) :: HRDIRQ
REAL(KIND=JPRB) :: HRDIRO3
REAL(KIND=JPRB) :: HRDIRPD
REAL(KIND=JPRB) :: HRDIRVD
REAL(KIND=JPRB) :: HRDIRSP
REAL(KIND=JPRB) :: RRDXTAU
REAL(KIND=JPRB) :: RDAMPVOR
REAL(KIND=JPRB) :: RDAMPDIV
REAL(KIND=JPRB) :: RDAMPT
REAL(KIND=JPRB) :: RDAMPQ
REAL(KIND=JPRB) :: RDAMPO3
REAL(KIND=JPRB) :: RDAMPPD
REAL(KIND=JPRB) :: RDAMPVD
REAL(KIND=JPRB) :: RDAMPSP
LOGICAL :: LNEWHD
REAL(KIND=JPRB) :: REXPDH
REAL(KIND=JPRB) :: FRANDH
REAL(KIND=JPRB) :: SLEVDH
REAL(KIND=JPRB) :: SLEVDH1
REAL(KIND=JPRB) :: SLEVDH2
REAL(KIND=JPRB) :: SLEVDH3
INTEGER(KIND=JPIM) :: NSREFDH
REAL(KIND=JPRB) :: RATIO_HDI_TOP
INTEGER(KIND=JPIM) :: NPROFILEHD
REAL(KIND=JPRB) :: RPROFHDBT
REAL(KIND=JPRB) :: RPROFHDTP
REAL(KIND=JPRB) :: RPROFHDMX
REAL(KIND=JPRB) :: RPROFHDEX
LOGICAL :: LRDISPE_EC
REAL(KIND=JPRB),ALLOCATABLE:: RDIVOR(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDIDIV(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDITG(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDIGFL(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDIPD(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDIVD(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDISP(:)
REAL(KIND=JPRB),ALLOCATABLE:: RDHI(:,:,:)
LOGICAL :: LSTRHD
REAL(KIND=JPRB) :: HDTIME_STRHD
REAL(KIND=JPRB) :: KLSPONGE
LOGICAL :: LM0DAMP
REAL(KIND=JPRB) :: RBEGSPONGE


LOGICAL :: LTOP_VOR
INTEGER(KIND=JPIM) :: NTOP_VOR_TRUNC
INTEGER(KIND=JPIM) :: NTOP_VOR_BOT

!====== SEMI-LAGRANGIAN SCHEME ====================

! * FOR DP CALCULATION:
! RSAVEDP: SAVE LAT, LON, ETA DP COORDINATES
REAL(KIND=JPRB),ALLOCATABLE:: RSAVEDP(:,:,:,:)

!====== SEMI-LAGRANGIAN HORIZONTAL DIFFUSION SCHEME (SLHD) ====================

! * FOR SLHD INTERPOLATIONS:
! SLHDA   :    Scaling factor of the deformation in f(d) function
!              (including the model resolution correction)
! SLHDA0  :    Namelist variable allowing to compute SLHDA
!              (scaling factor of the deformation in f(d) function
!              without the model resolution correction)
! SLHDA0T :    Namelist variable to influence T
! SLHDB   :    Exponent of the deformation in f(d) function
! SLHDBT  :    Exponent of the deformation in f(d) function for T
! SLHDD0  :    Treshold for deformation tensor enhancement
! SLHDD00 :    Namelist value for deformation tensor enhancement treshold 
! SLHDD00T:    Namelist value for deformation tensor enhancement treshold for T
! SLHDDIV :    Weight for including horizontal divergence into deformation
! SLHDRATDDIV: Nondimensional enhancer of divergence based diffusion
! SLHDHOR :    Switch for computing flow deformation:
!                0 - along eta levels (sloped)
!                1 - along pressure levels (quasi horizontal)
! LSLHDHEAT:   If true, the triggering function for heat variables differs from 
!              the one for momentum variables.
! LSLHDSPONGE  Use SLHD in upper 20 hPa to be preventing grid-point storms there
! SLHD_P_LOW : lowermost pressure level of SLHD activity (in Pa)
! SLHD_P_HIGH: lowermost pressure level of full SLHD activity (in Pa)

! * THE "HDS" CHARACTERISTIC TIMES (obsolete):
! HDSRVOR : for diffusion of vorticity.
! HDSRDIV : for diffusion of divergence.
! HDSRVD  : for diffusion of vertical divergence (NH).

! * REVERSE OF THE "HDS" CHARACTERISTIC TIMES:
! HRDSRVOR : for diffusion of vorticity.
! HRDSRDIV : for diffusion of divergence.
! HRDSRVD  : for diffusion of vertical divergence (NH).

! RDAMPVORS: local enhancing coefficient for HDS diffusion of vorticity
! RDAMPDIVS: local enhancing coefficient for HDS diffusion of divergence
! RDAMPVDS : local enhancing coefficient for HDS diffusion of vert. divergence
! RDAMPHDS : ratio HRDSRDIV/HRDIRDIV.

! Coefficients RDS[X] generally write HRDSR[X]*gs(l)*f(n,N(l),n0(X),x0,r) where n0
!  is generally 0 excepted for vorticity where n0 may be 2 in some cases.
! REXPDHS  : order "r" of the diffusion (exponent for the wavenumber dependency).
! SLEVDHS  : first threshold for the pressure dependency scaled by VP00 used in function gs(l).
! SLEVDHS1 : first threshold for the pressure dependency scaled by VP00 used in function N(l).
! SLEVDHS2 : second threshold for the pressure dependency scaled by VP00 used in function N(l).
! SDRED    : variable modifying the vertical profile based on SLEVDH
!            ( g(l) becomes g(l)-SDRED in the "main" diffusion).

! * "HDS" LEVEL AND WAVENUMBER DEPENDENT INVERSE CHARACTERISTIC TIMES:
! RDSVOR   : for diffusion of vorticity.
! RDSDIV   : for diffusion of divergence.
! RDSVD    : for diffusion of NH vertical divergence variable.
! RDHS     : SLHD additional horizontal diffusion operator used for stretched ARPEGE.

! * MASKING FUNCTION FOR SLHD
! SLHD_MASK_U : controls default variables 
! SLHD_MASK_T : controls heat variables when LSLHDHEAT=true
! NLEV_SPONGE : lowermost level to compute SLHD

REAL(KIND=JPRB),ALLOCATABLE :: SLHDA(:)
REAL(KIND=JPRB) :: SLHDA0
REAL(KIND=JPRB) :: SLHDA0T
REAL(KIND=JPRB) :: SLHDB
REAL(KIND=JPRB) :: SLHDBT
REAL(KIND=JPRB),ALLOCATABLE :: SLHDD0(:)
REAL(KIND=JPRB) :: SLHDD00
REAL(KIND=JPRB) :: SLHDD00T
REAL(KIND=JPRB) :: SLHDDIV
REAL(KIND=JPRB) :: SLHDRATDDIV
REAL(KIND=JPRB) :: SLHDHOR
REAL(KIND=JPRB) :: SLHD_P_LOW
REAL(KIND=JPRB) :: SLHD_P_HIGH
LOGICAL :: LSLHDHEAT
LOGICAL :: LSLHDSPONGE
REAL(KIND=JPRB) :: HDSRVOR
REAL(KIND=JPRB) :: HDSRDIV
REAL(KIND=JPRB) :: HDSRVD
REAL(KIND=JPRB) :: HRDSRVOR
REAL(KIND=JPRB) :: HRDSRDIV
REAL(KIND=JPRB) :: HRDSRVD
REAL(KIND=JPRB) :: RDAMPVORS
REAL(KIND=JPRB) :: RDAMPDIVS
REAL(KIND=JPRB) :: RDAMPVDS
REAL(KIND=JPRB) :: RDAMPHDS
REAL(KIND=JPRB) :: REXPDHS
REAL(KIND=JPRB) :: SLEVDHS
REAL(KIND=JPRB) :: SLEVDHS1
REAL(KIND=JPRB) :: SLEVDHS2
REAL(KIND=JPRB) :: SDRED
REAL(KIND=JPRB),ALLOCATABLE:: RDSVOR(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDSDIV(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDSVD(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: RDHS(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SLHD_MASK_U(:)
REAL(KIND=JPRB),ALLOCATABLE:: SLHD_MASK_T(:)
INTEGER(KIND=JPIM) :: NLEV_SPONGE


!====== OTHER DIFFUSIVE PROCESSES (moved from YOMDYNA at CY45) ====================

! LRFRIC     : .T. = Rayleigh friction in horizontal (zonal) wind.
! LRFRICISOTR: isotropic Rayleigh friction (acts on zonal and meridian horizontal wind).
LOGICAL :: LRFRIC
LOGICAL :: LRFRICISOTR

!======  QUANTITIES TO CHANGE THE VARIABLE IN THE T-EQN =======================

! RCORDIT(NFLEVG)    : correction term at full-levels for diffusion of T.
! RCORDIH(0:NFLEVG)  : correction term at half-levels for SL T-eqn if RCMSMP0/=0
! RCORDIF(NFLEVG)    : correction term at full-levels for SL T-eqn if RCMSMP0/=0

REAL(KIND=JPRB),ALLOCATABLE:: RCORDIT(:)
REAL(KIND=JPRB),ALLOCATABLE:: RCORDIH(:)
REAL(KIND=JPRB),ALLOCATABLE:: RCORDIF(:)

!==== MAXIMUM V-WINDS ALLOWED IN THE SEMI-LAGRANGIAN MODEL ====================

! VMAX1   : if V>VMAX1 (SM) or SQRT(U**2+V**2)>VMAX1 (DM),
!           warning in the SL scheme.
! VMAX2   : if V>VMAX2 (SM) or SQRT(U**2+V**2)>VMAX2 (DM),
!           abort in the SL scheme.

REAL(KIND=JPRB) :: VMAX1
REAL(KIND=JPRB) :: VMAX2

!==== MAXIMUM D3 ALLOWED ======================================================

! LBOUND_D3  : bound absolute value of D3
! RMAX_D3    : maximum value allowed for absolute value of D3 when LBOUND_D3=.T.

REAL(KIND=JPRB) :: RMAX_D3
LOGICAL :: LBOUND_D3

!==== RAYLEIGH FRICTION =======================================================

! RKRF(NFLEVG) : coefficient of Rayleigh friction
! NMAXLEVRF    : maximum level for which Rayleigh friction is applied. If no
!                Rayleigh friction is applied, we set NMAXLEVRF=0
! RRFZ1        : reference value for height of profile
! RRFPLM       : pressure limit - no Rayleigh friction for p > RRFPLM
! RRFTAU       : e-folding time for Rayleigh friction.

REAL(KIND=JPRB),ALLOCATABLE:: RKRF(:) 
INTEGER(KIND=JPIM) :: NMAXLEVRF
REAL(KIND=JPRB) :: RRFZ1
REAL(KIND=JPRB) :: RRFPLM
REAL(KIND=JPRB) :: RRFTAU

!==== UPPER RADIATIVE BOUNDARY CONDITION ======================================

! RTEMRB - tuning temperature for upper radiative b. c. (LRUBC)
! NRUBC   : control of radiative upper boundary condition :
!           =0 <=> non computation
!           =1 <=> computation on the forecast field
!           =2 <=> computation on the departure of the forecast from the coupling field

REAL(KIND=JPRB) :: RTEMRB
INTEGER(KIND=JPIM) :: NRUBC

!==== SEMI-IMPLICIT SCHEME, VERTICAL EIGENMODES, PC SCHEMES ===================

! LSIDG   : .F.: Semi-implicit-scheme with reduced divergence.
!           .T.: Semi-implicit scheme with not reduced divergence.

! BETADT  : coefficient for the semi-implicit treatment of divergence,
!           temperature, continuity (and NH if required) equations.
! RBT     : BETADT multiplied by coefficients depending on VESL, XIDT
! RBTS2   : 0.5*RBT
! REFGEO  : reference geopotentiel for shallow-water model.
! SIPR    : reference surface pressure.
! SITR    : reference temperature.
! SITRA   : namelist acoustic reference temperature.
! SITRAM  : model acoustic reference temperature.
! NOPT_SITRA: option for SITRAM (0 for vertically constant SITRAM, 1 for vertically-dependent SITRAM)
! SITRUB : ref. temper. for SI corr. of temper.(for LRUBC=.T.)
! SIPRUB : coef. for SI corr. of surf. press.  (for LRUBC=.T.)
! SITIME  : =TDT (if not, Helmholtz matrices are recomputed in CNT4).
! SIRPRG  : auxiliary variable for SIGAM,SIGAMA.
! SIRPRN  : auxiliary variable for SITNU,SITNUA
! NSITER  : number of iterations to treat the non linear semi-implicit terms
!           in the non-hydrostatic scheme.
! NCURRENT_ITER : for LNHDYN with PC scheme - current iteration: 
!                   0                 - predictor
!                   1, 2, ..., NSITER - correctors
! LRHDI_LASTITERPC: T (resp. F): when a PC scheme is activated (for example
!  LPC_FULL=.T.), the horizontal diffusion is done at the last iteration
!  of the corrector step (resp. all iterations of the predictor-corrector
!  scheme).
! NITERHELM : in the NH model, when the C1 constraint is not matched,
!  NITERHELM is the number of corrector iterations required to solve
!  the implicit part of the semi-implicit scheme (including the Helmoltz eqn).
!  This variable is not used in the hydrostatic model, and in the NH model
!  when the constraint C1 is matched: in this case there is no corrector
!  iteration in the SI scheme.

! * PRESSURES LINKED TO A REFERENCE PRESSURE = SIPR
! SIALPH(NFLEVG)  : coefficients "alpha" of hydrostatics.
! SILNPR(NFLEVG)  : Log of ratio of pressures between levels.
! SIDELP(NFLEVG)  : pressure differences across layers.
! SIRDEL(NFLEVG)  : their inverse.
! SITLAH(0:NFLEVG): half-level pressures.
! SITLAF(NFLEVG)  : full-level pressures.
! SIDPHI(NFLEVG)  : geopotential differences across layers.

! SIB(NFLEVG,NFLEVG)   : operator "B" of the SI scheme (DIV ===> DP/DT=B.DIV).
! SIMO(NFLEVG,NFLEVG)  : eigenvectors of "B".
! SIMI(NFLEVG,NFLEVG)  : SIMO**-1
! SIVP(NFLEVG)         : eigenvalues of "B".
! SIHEG(NFLEVG,(NSMAX+1)*(NSMAX+2)/2,3), SIHEG2(NFLEVG,NSMAX+1,2:3):
!  Helmholtz operator in case of SI computations with not reduced divergence. 
! SIHEGB(NFLEVG,(NSMAX+1)*(NSMAX+2)/2,3), SIHEGB2(NFLEVG,NSMAX+1,2:3):
!  Additional operators in case of LSIDG=T SI computations in the NH model.

! SIFAC : Used in SI scheme (NHEE model).
!         For example, contains:
!         [ 1 - beta**2 (Delta t)**2 C**2 (SITR/SITRA) (LLstar/H**2) ]
! SIFACI: Inverse of SIFAC. 

! SI_ILAPKSSI: Used in SI scheme (NHQE model).
!  Contains LLstarstar_kappa**-1 and (LLstarstar_kappa**-1 * SSstar_kappa)

! VNORM : constant for new scaling.

LOGICAL :: LSIDG
REAL(KIND=JPRB) :: BETADT
REAL(KIND=JPRB) :: RBT
REAL(KIND=JPRB) :: RBTS2
REAL(KIND=JPRB) :: REFGEO
REAL(KIND=JPRB) :: SIPR
REAL(KIND=JPRB) :: SITR
REAL(KIND=JPRB) :: SITRA
REAL(KIND=JPRB) :: SITRUB
REAL(KIND=JPRB) :: SIPRUB
REAL(KIND=JPRB) :: SITIME
REAL(KIND=JPRB) :: SIRPRG
REAL(KIND=JPRB) :: SIRPRN
INTEGER(KIND=JPIM) :: NSITER
INTEGER(KIND=JPIM) :: NCURRENT_ITER
LOGICAL :: LRHDI_LASTITERPC
INTEGER(KIND=JPIM) :: NITERHELM
INTEGER(KIND=JPIM) :: NOPT_SITRA

REAL(KIND=JPRB),ALLOCATABLE:: SIALPH(:)
REAL(KIND=JPRB),ALLOCATABLE:: SILNPR(:)
REAL(KIND=JPRB),ALLOCATABLE:: SIDELP(:)
REAL(KIND=JPRB),ALLOCATABLE:: SIRDEL(:)
REAL(KIND=JPRB),ALLOCATABLE:: SITLAH(:)
REAL(KIND=JPRB),ALLOCATABLE:: SITLAF(:)
REAL(KIND=JPRB),ALLOCATABLE:: SIDPHI(:)
REAL(KIND=JPRB),ALLOCATABLE:: SIB(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIMO(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIMI(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIVP(:)
REAL(KIND=JPRB),ALLOCATABLE:: SIHEG(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIHEG2(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIHEGB(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIHEGB2(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIFAC(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SIFACI(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SI_ILAPKSSI(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: SITRAM(:)
REAL(KIND=JPRB) :: VNORM

!=========== SEMI-LAGRANGIAN SWITCHES AND WEIGHTS =============================
!=========== + ADDITIONAL "ADVECTION" SWITCHES ALSO USED IN EULERIAN ========== 

! * Switches NxLAG:
! NVLAG   :  switch for formulation or discretisation of continuity equation.
! NWLAG   :  switch for formulation or discretisation of momentum equations.
! NTLAG   :  switch for formulation or discretisation of temperature equation.
! NSPDLAG :  switch for formulation or discretisation of P-hat equation.
! NSVDLAG :  switch for formulation or discretisation of d-hat equation.
! Remarks about NxLAG:
! a) possible value for NxLAG:
!    NxLAG=2 -> averaging of R.H.S. of the corresponding eq.
!               along the trajectory with the part corresponding
!               to the departure point added to the t-dt term
!    NxLAG=3 -> averaging of R.H.S. of the corresponding eq.
!               along the trajectory with the part corresponding
!               to the departure point interpolated linearly
!    NxLAG=4 -> same as NxLAG=4 with the tendency of physics
!               interpolated also linearly (applicable only for x=W,T,SVD)

! b) For NVLAG and 2D model: 
!    NVLAG>0 stands for the conventional formulation of continuity equation.
!    NVLAG<0 stands for the Lagrangian formulation of continuity equation:
!     in this case the remark a) is valid for ABS(NVLAG).

! NSPLTHOI : key to SPLiT High Order Interpolation
!            if NSPLTHOI /= 0 the quantity to be interpolated by high order
!            interpolation is split into:
!            i/  the field itself (note: for d4 it may contain d3 only) and
!            ii/ the remaining part (mostly tendency from physics)
!            if NSPLTHOI = 1 : i/ can use the diffusive interpolation while
!                              the ii/ part uses the non-diffusive int. only
!            if NSPLTHOI = -1: i/ and ii/ are interpolated by the same
!                              interpolation operator.
! LSPLTHOIGFL : model key to interpolate separately GFL fields
!               and their phys. tendencies (used only in SUSLB).
! NSLDIMK   : Number (dimension) of used horizontal non-linear weights in the model:
!            * default is 1
!            * L3DTURB=.T. adds +2 (first for KM, second for KH)
!            * NSPLITHOI=1 adds +1 (always the last one) for weights
!                         applied to physical tendencies

! * Research of semi-Lagrangian trajectory:
! NITMP   : Number of iterations for computing the medium point of the
!           semi-lagrangian trajectory.
! VETAON  : VETAON*eta(layer nr 1)+(1.-VETAON)*eta(top) is the lower
!           value allowed for ETA of the origin/anterior point in
!           the 3D model.
! VETAOX  : VETAOX*eta(bottom layer)+(1.-VETAOX)*eta(ground) is the
!           upper value allowed for ETA of the origin/anterior point
!           in the 3D model.
! RW2TLFF : when computing the refined position of the origin point for
!           Coriolis term, the new wind used is:
!           0.5*RW2TLFF*(V(F)+V(O)) + (1-RW2TLFF)*V(M)

! * WENO vertical interpolation for SL trajectory research
! RALPHA  : power for the linear weights denominator. Higher values
!           represents less overshooting but more diffusive results.
! RALPHA_TOP : The same as above for uppermost NLEV_ZALPHA levels
! NLEV_ZALPHA : Number of uppermost levels with special treatment
! NEDER   : defines the method to obtain smoothness indicator (see report of A. Craciun)
!           1 - L2-norm of high order variations of the polynomials reconstruction (eq.3)
!           2-4 - same as 1 with considering also derivatives (eq. 4-6)
!           5 - undivided differences methods  (eq. 7) 
! LWENOBC : .true. - no special boundary treatment for weno (results in more weight to 
!                    the outstanding result, typically the one from the high order interpolation)
!           .false. - boundaries are interpolated in the standard way (smoother in TL/AD)

! * Quasi-monotonicity for high order wind interpolation of SL trajectory
! NQMHOISLT : 0 no QM treatment
!             1 quasi-monotone along horizontal
!             2 quasi-monotone in 3D
!             3 quasi-monotone along vertical

! * SL trajectory researched by 4th order Runge-Kutta implicit method (Lobatto IIIA)
! LSLDP_RK  : when .TRUE. the iterative algorithm is replaced by 4th order Runge Kutta method

! LSLDP_CURV : on/off - horizontal wind components from lat,lon are transformed to cartesian coords
! LRHS_CURV  : on/off - for the right hand side interpolation
! LSLDP_XYZ  : on/off - when .TRUE. the SL trajectory iterations are performed entirely in 
!              cartesian x, y, z coordinates - needs LSLDP_CURV
! LSLDP_SAVE : on/off  the departure point at each timestep are saved and used as starting values
!              for next step iterations reducing their cost - work with LSLDP_XYZ only
!     in order to mitigate the pole singularity during the SL trajectory research.

! * Uncentering factor in the semi-Lagrangian scheme:
! VESL    : first order uncentering factor applied to non linear and linear
!           terms.
! XIDT    : pseudo-second order uncentering factor applied to linear terms,
!           when an alternative second-order averaging is required in the 
!           2TL SL scheme.

! * Switches for use of quasi-monotone interpolations:
! LQMW    :  Use quasi-monotone three-dimensional interpolations for wind
! LQMHW   :  Use quasi-monotone interpolations in the horizontal for wind
! LQMT    :  Use quasi-monotone three-dimensional interpolations for temperature
! LQMHT   :  Use quasi-monotone interpolations in the horizontal for temperature
! LQMP    :  Use quasi-monotone three-dimensional interpolations for cont. eq
! LQMHP   :  Use quasi-monotone interpolations in the horizontal for cont. eq
! LQMPD   :  Use quasi-monotone three-dimensional interpolations for P-hat eqn.
! LQMHPD  :  Use quasi-monotone interpolations in the horizontal for P-hat eqn.
! LQMVD   :  Use quasi-monotone three-dimensional interpolations for d-hat eqn.
! LQMHVD  :  Use quasi-monotone interpolations in the horizontal for d-hat eqn.

! * Switches for vertical quintic interpolation
! LVWENO_W   : Use vertical quintic interpolation for wind
! WENO_ALPHA_W   : WENO weight for wind
! LVWENO_T   : Use vertical quintic interpolation for temperature
! WENO_ALPHA_T   : WENO weight for temperature
! LVWENO_SP  : Use vertical quintic interpolation for cont. eq
! WENO_ALPHA_SP  : WENO weight for cont. eq
! LVWENO_SPD  : Use vertical quintic interpolation for P-hat eqn.
! WENO_ALPHA_SPD : WENO weight for P-hat eqn.
! LVWENO_SVD  : Use vertical quintic interpolation for d-hat eqn.
! WENO_ALPHA_SVD : WENO weight for d-hat eqn.

! * Treatment of Coriolis term:
! LADVF   : if TRUE then use "advective" treatment of Coriolis terms (SL);
!           in this case 2*Omega*Vec*r is computed analytically.
! LIMPF   : if TRUE then use implicit treatment of Coriolis terms (EUL and SL)
! L2TLFF  : if TRUE then use refined treatment of Coriolis term in 2TLSL scheme
!           (can be currently used also with the 3TL SL vertical interpolating
!           scheme).

! * Change variable with an Eulerian treatment of orography:
! RCMSLP0 : Real for tuning of the Tanguay/Ritchie correction in SL continuity
!           and temperature equations for 3D model.

! * Switch for computation of Moisture Convergence for French deep convection scheme

! NCOMP_CVGQ   :  0 ==> Compute the CVGQ in an Eulerian manner, using spectral
!                       moisture stored in the YQ GFL variable.
!                       In this case YQ must be spectral and
!                       horizontal derivatives are used.
!                 1 ==> Compute the CVGQ in an Eulerian manner, using spectral
!                       moisture stored in the YCVGQ GFL spectral variable and
!                       its horizontal derivatives.
!                       This case is well designed for the case where YQ is
!                       a purely grid-point GFL.
!                 2 ==> Compute the CVGQ in a semi-Lagrangian manner
!                       (Lagrangian tendency - Eulerian tendency), using data
!                       stored in the YCVGQ grid-point variable.
!                       This case is well designed for the case where YQ is
!                       a purely grid-point GFL, and where LSLAG=T.
! remark ky: better to move this variable in SUDYNA/NAMDYNA/YOMDYNA in the
!  future to make it available in SUDIM1 when reading NAMGFL.

! * Stratospheric vertical trajectory smoothing in the SL scheme
! LSVTSM : Stratospheric vertical trajectory smoothed
! RPRES_SVTSM : smoothing done for standard pressure < RPRES_SVTSM

! * SETTLST filter for vertical trajectories in SL scheme 
! LSETTLSVF : filter applied - vertical SL traj with stable scheme
! LSETFSTAT : enable printing of output diagnostics when LSETTLSVF=true
!             specifying the % of pts that SETTLST extrapol switched off 
! RPRES_SETTLSVF : filter applied for standard pressure < RPRES_SETTLSVF
! NFLEVSF : number of levels to apply SETTLST filter
! RSCALE :  scaling factor for transition between 1st and 2nd order zones
!           higher value results in a steeper function more realistically
!           mimicking the original jump but resulting in a less smooth
!           transition between the two extremes
! RSCALEOFF : Offset to extend second order accuracy regime by areas with
!             nearly no vertical velocity
! LSETTLSVF_DIFF : enhance SLHD diffusion for 1st order areas

! * Optimizing dynamics range to help model precision (maaainly in adjoint)
! RINTOPT   : Parameter to multiply all advected GFL fields before interpolation
! RRINTOPT  : 1./RINTOPT

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: NVLAG
INTEGER(KIND=JPIM) :: NWLAG
INTEGER(KIND=JPIM) :: NTLAG
INTEGER(KIND=JPIM) :: NSPDLAG
INTEGER(KIND=JPIM) :: NSVDLAG
INTEGER(KIND=JPIM) :: NSPLTHOI
LOGICAL :: LSPLTHOIGFL
INTEGER(KIND=JPIM) :: NSLDIMK
INTEGER(KIND=JPIM) :: NITMP
REAL(KIND=JPRB) :: VETAON
REAL(KIND=JPRB) :: VETAOX
LOGICAL :: LSETTLSVF
LOGICAL :: LSETFSTAT
LOGICAL :: LSETTLSVF_DIFF
REAL(KIND=JPRB) :: RW2TLFF
REAL(KIND=JPRB) :: RALPHA
REAL(KIND=JPRB) :: RALPHA_TOP
INTEGER(KIND=JPIM) :: NLEV_ZALPHA
INTEGER(KIND=JPIM) :: NEDER
LOGICAL :: LWENOBC
INTEGER(KIND=JPIM) :: NQMHOISLT
LOGICAL :: LSLDP_RK
LOGICAL :: LSLDP_CURV
LOGICAL :: LRHS_CURV
LOGICAL :: LSLDP_XYZ
LOGICAL :: LSLDP_SAVE
REAL(KIND=JPRB) :: VESL
REAL(KIND=JPRB) :: XIDT
LOGICAL :: LQMW
LOGICAL :: LQMHW
LOGICAL :: LQMT
LOGICAL :: LQMHT
LOGICAL :: LQMP
LOGICAL :: LQMHP
LOGICAL :: LQMPD
LOGICAL :: LQMHPD
LOGICAL :: LQMVD
LOGICAL :: LQMHVD
LOGICAL :: LVWENO_W, LVWENO_T, LVWENO_SP, LVWENO_SPD, LVWENO_SVD
REAL(KIND=JPRB) :: WENO_ALPHA_W, WENO_ALPHA_T, WENO_ALPHA_SP, WENO_ALPHA_SPD, WENO_ALPHA_SVD
LOGICAL :: LADVF
LOGICAL :: LIMPF
LOGICAL :: L2TLFF
REAL(KIND=JPRB) :: RCMSLP0
INTEGER(KIND=JPIM) :: NCOMP_CVGQ
LOGICAL :: LSVTSM
REAL(KIND=JPRB) :: RPRES_SVTSM
REAL(KIND=JPRB) :: RPRES_SETTLSVF
INTEGER(KIND=JPIM) :: NFLEVSF
REAL(KIND=JPRB) :: RSCALE
REAL(KIND=JPRB) :: RSCALEOFF
REAL(KIND=JPRB) :: RINTOPT
REAL(KIND=JPRB) :: RRINTOPT

! -----------------------------------------------------------------------------

!=========== VARIABLES FOR ECMWF DIFFUSION SETUP =========================!

! horizontal diffusion serves 4 purposes, 
! a) to compensate for actual viscosity (physics)
! b) to stop the build-up of energy/enstrophy at the smallest scales
! c) as a sponge at the top of the model
! d) to keep the TL and NL evolution similar 
! With LGRADSP=T (b) does not happen anymore, but there is still "physics" missing 
! (i.e. averaged fluctuations on the resolved grid) due to the closure problem. 
! We also still need the sponge layer to prevent unphysical reflections of vertically 
! propagating gravity waves from the model top.
! Given that (b) is solved by the LGRADSP=T filter, less strong (horizontal) diffusion may be applied 
! in the model integration.
! The switches LHDIFFM and NDIFFACT express the strength of horizontal diffusion as a multiple of the time-step.
! The switch LGPSTRESS=T (+tuning parameters RCLSTRESS, RCLPOLE) provide
! an alternative non-linear diffusion (together with spectral viscosity) to add physically based averaged
! fluctuations. For simplicity the sponge is left as tuned in the standard diffusion case.

! LGRADSP   : special switch for de-aliasing the pressure gradient term ( now in YOMDYNA)
! LHDIFFM   : if true, apply horizontal diffusion as a multiple of the time-step (outside sponge)
! NDIFFACT  : specifies the multiple of the time-step for horizontal diffusion (outside sponge)
! LGPSTRESS : switch for computing 2D stress tensor with dynamic similarity model (default==FALSE),
!             in this case spectral viscosity is used instead of horizontal diffusion (outside sponge)
! LSPECVIS  : use spectral viscosity, this is the default for the cubic grid (see Gelb+Gleeson, MWR 129, 2001)

! 1) default for model run: LGRADSP=T, LGPSTRESS=F, LHDIFFM=T, NDIFFACT=6-8
! 2) default for data assimilation: LGRADSP=T, LGPSTRESS=F, LHDIFFM=F
! 3) default for quadratic and cubic grid: LSPECVIS=T, LHDIFFM=F

LOGICAL :: LHDIFFM
LOGICAL :: LSPECVIS
INTEGER(KIND=JPIM) :: NDIFFACT
LOGICAL :: LGPSTRESS
REAL(KIND=JPRB) :: RCLSTRESS
REAL(KIND=JPRB) :: RCLPOLE

!     ------------------------------------------------------------------

!     USED TO IMPROVE VECTORIZATION OF SEMI-LAGRANGIAN ADJOINT

!     (=> safe vertical separation for updates)

! NVSEPC    : number of vertical blocks of layers for vertical loops
!             in the cubic interpolation routines.
! NVSEPL    : number of vertical blocks of layers for vertical loops
!             in the linear interpolation routines.
! LFINDVSEP : computation of NVSEPC and NVSEPL is done if .TRUE.
INTEGER(KIND=JPIM) :: NVSEPC, NVSEPL
LOGICAL ::      LFINDVSEP

!-----------------------------------------------------------------------

!*    Global mass variables required for mass correction,
!     to prevent mass-drift in extended integrations.

!     Global average values are calculated in spnorm
!     Nils Wedi, 2008-02-08

! GMASSI    : mass at start of integration.
! GMASS0    : mass at current timestep.
! GMASSINC  : mass increment to be applied at current time step

! LMASCOR   : .T. apply mass correction
! LMASDRY   : .F. by default, if .T. only correct the dry mass by subtracting
!             the total mass of water vapour (see cormass2 used at Meteo-France)


LOGICAL         :: LMASCOR
LOGICAL         :: LMASDRY

LOGICAL         :: LGPMASCOR

INTEGER(KIND=JPIM) :: NGPMASCOR

REAL(KIND=JPRB) :: GPMASSI

REAL(KIND=JPRB) :: GMASSI
REAL(KIND=JPRB) :: GMASS0
REAL(KIND=JPRB) :: GMASSINC

! Modify CTY equation to allow dry mass conservation
! The total water tendencies from physics are taken into account in GPCTY
! If LNODRYFLX, the global mass fixer must be LMASDRY
LOGICAL         :: LNODRYFLX

!-------------------------------------------------------------------------------
!     SIBI    : INVERSE OF VERTICAL STRUCTURE MATRIX (B) - NEEDED FOR CVAR2...
!
! It is in here, rather than in YOMJG, because it is set up in SUDYN!

REAL(KIND=JPRB), ALLOCATABLE :: SIBI(:,:)

!---------------------------------------------------------------------

CONTAINS

PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

END TYPE TDYN

!!TYPE(TDYN), POINTER :: YRDYN => NULL()

 !---------------------------------------------------------------------

CONTAINS 

  SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)

    USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

    IMPLICIT NONE
    CLASS(TDYN), INTENT(IN) :: SELF
    INTEGER(KIND=JPIM)    , INTENT(IN) :: KDEPTH
    INTEGER(KIND=JPIM)    , INTENT(IN) :: KOUTNO

    INTEGER(KIND=JPIM) :: IDEPTHLOC
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('YOMDYN:PRINT_CONFIGURATION',0,ZHOOK_HANDLE)
    IDEPTHLOC = KDEPTH+2

    WRITE(KOUTNO,*) REPEAT(' ',KDEPTH)    // 'model%yrml_dyn%yrdyn : '
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REPS1 = ', SELF%REPS1
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REPS2 = ', SELF%REPS2
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REPSM1 = ', SELF%REPSM1
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REPSM2 = ', SELF%REPSM2
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REPSP1 = ', SELF%REPSP1
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HDIRVOR = ', SELF%HDIRVOR
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HDIRDIV = ', SELF%HDIRDIV
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HDIRT = ', SELF%HDIRT
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HDIRQ = ', SELF%HDIRQ
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HDIRO3 = ', SELF%HDIRO3
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HDIRPD = ', SELF%HDIRPD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HDIRVD = ', SELF%HDIRVD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HDIRSP = ', SELF%HDIRSP
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HDIRVOR = ', SELF%HRDIRVOR
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HRDIRDIV = ', SELF%HRDIRDIV
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HRDIRT = ', SELF%HRDIRT
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HRDIRQ = ', SELF%HRDIRQ
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HRDIRO3 = ', SELF%HRDIRO3
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HRDIRPD = ', SELF%HRDIRPD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HRDIRVD = ', SELF%HRDIRVD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HRDIRSP = ', SELF%HRDIRSP
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RRDXTAU = ', SELF%RRDXTAU
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDAMPVOR = ', SELF%RDAMPVOR
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDAMPDIV = ', SELF%RDAMPDIV
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDAMPT = ', SELF%RDAMPT
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDAMPQ = ', SELF%RDAMPQ
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDAMPO3 = ', SELF%RDAMPO3
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDAMPPD = ', SELF%RDAMPPD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDAMPVD = ', SELF%RDAMPVD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDAMPSP = ', SELF%RDAMPSP
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LNEWHD = ', SELF%LNEWHD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REXPDH = ', SELF%REXPDH
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'FRANDH = ', SELF%FRANDH
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLEVDH = ', SELF%SLEVDH
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLEVDH1 = ', SELF%SLEVDH1
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLEVDH2 = ', SELF%SLEVDH2
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLEVDH3 = ', SELF%SLEVDH3
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSREFDH = ', SELF%NSREFDH
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NPROFILEHD = ', SELF%NPROFILEHD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPROFHDBT = ', SELF%RPROFHDBT
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPROFHDTP = ', SELF%RPROFHDTP
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPROFHDMX = ', SELF%RPROFHDMX
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPROFHDEX = ', SELF%RPROFHDEX
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LRDISPE_EC = ', SELF%LRDISPE_EC
    IF (ALLOCATED(SELF%RDIVOR)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDIVOR allocated of shape ', SHAPE( SELF%RDIVOR )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDIVOR not allocated '
    ENDIF
    IF (ALLOCATED(SELF%RDIDIV)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDIDIV allocated of shape ', SHAPE( SELF%RDIDIV )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDIDIV not allocated '
    ENDIF
    IF (ALLOCATED(SELF%RDITG)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDITG allocated of shape ', SHAPE( SELF%RDITG )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDITG not allocated '
    ENDIF
    IF (ALLOCATED(SELF%RDIGFL)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDIGFL allocated of shape ', SHAPE( SELF%RDIGFL )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDIGFL not allocated '
    ENDIF
    IF (ALLOCATED(SELF%RDIPD)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDIPD allocated of shape ', SHAPE( SELF%RDIPD )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDIPD not allocated '
    ENDIF
    IF (ALLOCATED(SELF%RDIVD)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDIVD allocated of shape ', SHAPE( SELF%RDIVD )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDIVD not allocated '
    ENDIF
    IF (ALLOCATED(SELF%RDISP)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDISP allocated of shape ', SHAPE( SELF%RDISP )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDISP not allocated '
    ENDIF
    IF (ALLOCATED(SELF%RDHI)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDHI allocated of shape ', SHAPE( SELF%RDHI )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDHI not allocated '
    ENDIF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTRHD = ', SELF%LSTRHD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HDTIME_STRHD = ', SELF%HDTIME_STRHD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LTOP_VOR = ', SELF%LTOP_VOR
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTOP_VOR_TRUNC = ', SELF%NTOP_VOR_TRUNC
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTOP_VOR_BOT = ', SELF%NTOP_VOR_BOT
    IF (ALLOCATED(SELF%RSAVEDP)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSAVEDP allocated of shape ', SHAPE( SELF%RSAVEDP )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSAVEDP not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SLHDA)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHDA allocated of shape ', SHAPE( SELF%SLHDA )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHDA not allocated '
    ENDIF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHDA0 = ', SELF%SLHDA0
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHDA0T = ', SELF%SLHDA0T
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHDB = ', SELF%SLHDB
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHDBT = ', SELF%SLHDBT
    IF (ALLOCATED(SELF%SLHDD0)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHDD0 allocated of shape ', SHAPE( SELF%SLHDD0 )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHDD0 not allocated '
    ENDIF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHDD00 = ', SELF%SLHDD00
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHDD00T = ', SELF%SLHDD00T
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHDDIV = ', SELF%SLHDDIV
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHDRATDDIV = ', SELF%SLHDRATDDIV
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHDHOR = ', SELF%SLHDHOR
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHD_P_LOW = ', SELF%SLHD_P_LOW
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLHD_P_HIGH = ', SELF%SLHD_P_HIGH
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSLHDHEAT = ', SELF%LSLHDHEAT
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HDSRVOR = ', SELF%HDSRVOR
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HDSRDIV = ', SELF%HDSRDIV
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HDSRVD = ', SELF%HDSRVD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HRDSRVOR = ', SELF%HRDSRVOR
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HRDSRDIV = ', SELF%HRDSRDIV
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'HRDSRVD = ', SELF%HRDSRVD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDAMPVORS = ', SELF%RDAMPVORS
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDAMPDIVS = ', SELF%RDAMPDIVS
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDAMPVDS = ', SELF%RDAMPVDS
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDAMPHDS = ', SELF%RDAMPHDS
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REXPDHS = ', SELF%REXPDHS
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLEVDHS = ', SELF%SLEVDHS
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLEVDHS1 = ', SELF%SLEVDHS1
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLEVDHS2 = ', SELF%SLEVDHS2
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SDRED = ', SELF%SDRED
    IF (ALLOCATED(SELF%RDSVOR)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDSVOR allocated of shape ', SHAPE( SELF%RDSVOR )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDSVOR not allocated '
    ENDIF
    IF (ALLOCATED(SELF%RDSDIV)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDSDIV allocated of shape ', SHAPE( SELF%RDSDIV )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDSDIV not allocated '
    ENDIF
    IF (ALLOCATED(SELF%RDSVD)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDSVD allocated of shape ', SHAPE( SELF%RDSVD )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDSVD not allocated '
    ENDIF
    IF (ALLOCATED(SELF%RDHS)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDHS allocated of shape ', SHAPE( SELF%RDHS )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDHS not allocated '
    ENDIF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LRFRIC = ', SELF%LRFRIC
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LRFRICISOTR = ', SELF%LRFRICISOTR

    IF (ALLOCATED(SELF%RCORDIT)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCORDIT allocated of shape ', SHAPE( SELF%RCORDIT )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCORDIT not allocated '
    ENDIF
    IF (ALLOCATED(SELF%RCORDIH)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCORDIH allocated of shape ', SHAPE( SELF%RCORDIH )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCORDIH not allocated '
    ENDIF
    IF (ALLOCATED(SELF%RCORDIF)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCORDIF allocated of shape ', SHAPE( SELF%RCORDIF )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCORDIF not allocated '
    ENDIF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'VMAX1 = ', SELF%VMAX1
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'VMAX2 = ', SELF%VMAX2
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RMAX_D3 = ', SELF%RMAX_D3
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LBOUND_D3 = ', SELF%LBOUND_D3
    IF (ALLOCATED(SELF%RKRF)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RKRF allocated of shape ', SHAPE( SELF%RKRF )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RKRF not allocated '
    ENDIF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NMAXLEVRF = ', SELF%NMAXLEVRF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RRFZ1 = ', SELF%RRFZ1
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RRFPLM = ', SELF%RRFPLM
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RRFTAU = ', SELF%RRFTAU
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RTEMRB = ', SELF%RTEMRB
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRUBC = ', SELF%NRUBC
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSIDG = ', SELF%LSIDG
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'BETADT = ', SELF%BETADT
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RBT = ', SELF%RBT
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RBTS2 = ', SELF%RBTS2
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REFGEO = ', SELF%REFGEO
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIPR = ', SELF%SIPR
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SITR = ', SELF%SITR
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SITRA = ', SELF%SITRA
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SITRUB = ', SELF%SITRUB
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIPRUB = ', SELF%SIPRUB
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SITIME = ', SELF%SITIME
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIRPRG = ', SELF%SIRPRG
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIRPRN = ', SELF%SIRPRN
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSITER = ', SELF%NSITER
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NCURRENT_ITER = ', SELF%NCURRENT_ITER
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LRHDI_LASTITERPC = ', SELF%LRHDI_LASTITERPC
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NITERHELM = ', SELF%NITERHELM

    IF (ALLOCATED(SELF%SIALPH)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIALPH allocated of shape ', SHAPE( SELF%SIALPH )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIALPH not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SILNPR)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SILNPR allocated of shape ', SHAPE( SELF%SILNPR )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SILNPR not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SIDELP)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIDELP allocated of shape ', SHAPE( SELF%SIDELP )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIDELP not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SIRDEL)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIRDEL allocated of shape ', SHAPE( SELF%SIRDEL )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIRDEL not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SITLAH)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SITLAH allocated of shape ', SHAPE( SELF%SITLAH )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SITLAH not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SITLAF)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SITLAF allocated of shape ', SHAPE( SELF%SITLAF )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SITLAF not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SIDPHI)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIDPHI allocated of shape ', SHAPE( SELF%SIDPHI )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIDPHI not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SIB)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIB allocated of shape ', SHAPE( SELF%SIB )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIB not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SIMO)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIMO allocated of shape ', SHAPE( SELF%SIMO )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIMO not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SIMI)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIMI allocated of shape ', SHAPE( SELF%SIMI )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIMI not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SIVP)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIVP allocated of shape ', SHAPE( SELF%SIVP )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIVP not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SIHEG)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIHEG allocated of shape ', SHAPE( SELF%SIHEG )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIHEG not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SIHEG2)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIHEG2 allocated of shape ', SHAPE( SELF%SIHEG2 )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIHEG2 not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SIHEGB)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIHEGB allocated of shape ', SHAPE( SELF%SIHEGB )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIHEGB not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SIHEGB2)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIHEGB2 allocated of shape ', SHAPE( SELF%SIHEGB2 )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIHEGB2 not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SIFAC)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIFAC allocated of shape ', SHAPE( SELF%SIFAC )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIFAC not allocated '
    ENDIF
    IF (ALLOCATED(SELF%SIFACI)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIFACI allocated of shape ', SHAPE( SELF%SIFACI )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIFACI not allocated '
    ENDIF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'VNORM = ', SELF%VNORM
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVLAG = ', SELF%NVLAG
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NWLAG = ', SELF%NWLAG
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTLAG = ', SELF%NTLAG
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSPDLAG = ', SELF%NSPDLAG
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSVDLAG = ', SELF%NSVDLAG
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSPLTHOI = ', SELF%NSPLTHOI
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSPLTHOIGFL = ', SELF%LSPLTHOIGFL
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLDIMK = ', SELF%NSLDIMK
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NITMP = ', SELF%NITMP
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'VETAON = ', SELF%VETAON
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'VETAOX = ', SELF%VETAOX
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSETTLSVF = ', SELF%LSETTLSVF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSETFSTAT = ', SELF%LSETFSTAT
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RW2TLFF = ', SELF%RW2TLFF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'VESL = ', SELF%VESL
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'XIDT = ', SELF%XIDT
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LQMW = ', SELF%LQMW
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LQMHW = ', SELF%LQMHW
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LQMT = ', SELF%LQMT
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LQMHT = ', SELF%LQMHT
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LQMP = ', SELF%LQMP
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LQMHP = ', SELF%LQMHP
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LQMPD = ', SELF%LQMPD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LQMHPD = ', SELF%LQMHPD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LQMVD = ', SELF%LQMVD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LQMHVD = ', SELF%LQMHVD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LVWENO_W = ', SELF%LVWENO_W
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LVWENO_T = ', SELF%LVWENO_T
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LVWENO_SP = ', SELF%LVWENO_SP
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LVWENO_SPD = ', SELF%LVWENO_SPD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LVWENO_SVD = ', SELF%LVWENO_SVD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'WENO_ALPHA_W = ', SELF%WENO_ALPHA_W
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'WENO_ALPHA_T = ', SELF%WENO_ALPHA_T
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'WENO_ALPHA_SP = ', SELF%WENO_ALPHA_SP
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'WENO_ALPHA_SPD = ', SELF%WENO_ALPHA_SPD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'WENO_ALPHA_SVD = ', SELF%WENO_ALPHA_SVD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LADVF = ', SELF%LADVF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LIMPF = ', SELF%LIMPF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'L2TLFF = ', SELF%L2TLFF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCMSLP0 = ', SELF%RCMSLP0
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NCOMP_CVGQ = ', SELF%NCOMP_CVGQ
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSVTSM = ', SELF%LSVTSM
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPRES_SVTSM = ', SELF%RPRES_SVTSM
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPRES_SETTLSVF = ', SELF%RPRES_SETTLSVF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSCALE = ', SELF%RSCALE
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSCALEOFF = ', SELF%RSCALEOFF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RINTOPT = ', SELF%RINTOPT
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RRINTOPT = ', SELF%RRINTOPT
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFLEVSF = ', SELF%NFLEVSF
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LHDIFFM = ', SELF%LHDIFFM
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSPECVIS = ', SELF%LSPECVIS
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDIFFACT = ', SELF%NDIFFACT
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LGPSTRESS = ', SELF%LGPSTRESS
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCLSTRESS = ', SELF%RCLSTRESS
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCLPOLE = ', SELF%RCLPOLE
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVSEPC and NVSEPL = ', SELF%NVSEPC, SELF%NVSEPL
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LFINDVSEP = ', SELF%LFINDVSEP
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LMASCOR = ', SELF%LMASCOR
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LMASDRY = ', SELF%LMASDRY

    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LGPMASCOR = ', SELF%LGPMASCOR
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NGPMASCOR = ', SELF%NGPMASCOR
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GPMASSI = ', SELF%GPMASSI
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GMASSI = ', SELF%GMASSI
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GMASS0 = ', SELF%GMASS0
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GMASSINC = ', SELF%GMASSINC
    IF (ALLOCATED(SELF%SIBI)) THEN
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIBI allocated of shape ', SHAPE( SELF%SIBI )
    ELSE
      WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SIBI not allocated'
    ENDIF
IF (LHOOK) CALL DR_HOOK('YOMDYN:PRINT_CONFIGURATION',1,ZHOOK_HANDLE)
    
  END SUBROUTINE PRINT_CONFIGURATION

!     ------------------------------------------------------------------
END MODULE YOMDYN
