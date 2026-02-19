! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE CALLPAR(YDGEOMETRY,YDVARS,YDSURF,YDMODEL,KDIM,&
 !-----------------------------------------------------------------------
 & PAUX, PRAD, FLUX, PDIAG, PSURF, PCGPP, PCREC, PAG, PRECO, PDDHS, AUXL, SURFL, LLKEYS, PERTL, &
 ! - Model variables (t)
 & PGFL,&
 & PPERT, PSLPHY9,&
 ! - UPDATED TENDENCY
 & PTENGFL,&
 & PHYS_MWAVE,&
 ! Stored quantities
 & PSAVTEND, PGFLSLP,&
 & STATE_T0, TENDENCY_DYN, TENDENCY_CML, STATE_TMP, TENDENCY_TMP,&
 & TENDENCY_VDF, TENDENCY_SATADJ, TENDENCY_LOC, TENDENCY_PHY&
 & )

!**** *CALLPAR * - CALL ECMWF PHYSICS

!     PURPOSE.
!     --------
!     - CALL THE SUBROUTINES OF THE E.C.M.W.F. PHYSICS PACKAGE.

!     ******************************************************************
!     ****** IDIOSYNCRASIES *** IDIOSYNCRASIES *** IDIOSYNCRASIES ******
!     ******************************************************************
!     ***  HEALTH WARNING:                                           ***
!     ***  ===============                                           ***
!     ***  NOTE THAT WITHIN THE E.C.M.W.F. PHYSICS HALF-LEVELS       ***
!     ***  ARE INDEXED FROM 1 TO NFLEVG+1 WHILE THEY ARE BETWEEN     ***
!     ***  0 AND NFLEVG IN THE REST OF THE MODEL. THE CHANGE IS TAKEN***
!     ***  CARE OF IN THE CALL TO THE VARIOUS SUBROUTINES OF THE     ***
!     ***  PHYSICS PACKAGE                                           ***
!     ***                                                            ***
!     ***    THIS IS SUPPOSED TO BE A "TEMPORARY" FEATURE TO BE      ***
!     ***    STRAIGHTENED OUT IN THE "NEAR" FUTURE                   ***
!     ******************************************************************
!     ******************************************************************
!     ***  NOTE THAT WITHIN THE E.C.M.W.F. PHYSICS PACKAGE, SNOW     ***
!     ***  AND MOISTURE CONTENT ARE IN METERS OF WATER, THUS THE     ***
!     ***  CONVERSION FACTOR (1.E-3) FOR THE FIRST LAYER             ***
!     ******************************************************************
!     ******************************************************************
!     ***  MOREOVER, WATER IN DEEPER LAYERS HAS TO BE NORMALIZED TO  ***
!     ***  THE FIRST LAYER DEPTH                                     ***
!     ***  THIS IS EQUIVALENT TO A CHANGE OF UNITS FROM INTEGRATED   ***
!     ***  WATER UNITS (KG/M**2, USED EVERYWEHRE ELSE IN THE CODE)   ***
!     ***  TO VOLUMETRIC WATER UNITS (M/0.07M, USED IN THE PHYSICS)  ***
!     ******************************************************************
!     ******************************************************************

!**   Interface.
!     ----------
!        *CALL* *CALLPAR*

!-----------------------------------------------------------------------

! -   ARGUMENTS.
!  --------------

! KDIM    : derived variable for dimensions
! PAUX    : derived variable for auxiliary quantity
! PRAD    : d.v. for quantities used in radiation scheme
! FLUX    : d.v. for fluxes
! PDIAG   : d.v. for diagnostics quantities
! PSURF   : d.v. for surface and other quantities (sharing the same structure)
! PDDHS   : d.v. for surface DDH (diagnostics) quantities

! PCGPP   : CO2 GPP flux adjustment coefficient
! PCREC   : CO2 REC flux adjustment coefficient
! PAG     : CO2 GPP flux
! PRECO   : CO2 REC flux

! PU      : X-COMPONENT OF WIND.
! PV      : Y-COMPONENT OF WIND.
! PT      : TEMPERATURE.
! PGFL    : GFL FIELDS

! PTDL    : zonal gradient of T
! PTDM    : merid gradient of T
! PTDU    : zonal gradient of U
! PTDV    : zonal gradient of V
! PVOR    : vorticity
! PIV     : divergence

! PPERT   : d.v. for perturbations etc...
! PSLPHY9 : d.v. for contributions from previous timesteps interpolated to O point of SL traj.

! PTENGFL    : TENDENCY OF ALL GFL FIELDS

! PSAVTEND   : ARRAY OF GMV TENDENCIES TO BE SAVED FOR NEXT TIME STEP
! PGFLSLP    : ARRAY OF GFL TENDENCIES TO BE SAVED FOR NEXT TIME STEP

!-----------------------------------------------------------------------

!     Externals.  to be updated...
!     ---------

!     Method. See documentation.
!     -------

!     AUTHOR.
!     -------
!      ORIGINAL 93-10-04 M.HAMRUD/P.VITERBO (FROM APLPAR)

!     MODIFICATIONS.
!     --------------
!     Modified 01-11-27 S. ABDALLA: Filling wind arrays for wave coupling
!                       moved to Sub. CPGLAG ..&.. PZIDLWV (Zi/L) added.
!     Modified 2001-10-10 JJMorcrette CCNs
!     Modified 15-10-01 D.Salmond FULLIMP mods
!     Modified 15-05-02 D.Dent    Code for extra fields
!     Modified 20-11-02 M. Janiskova: Call for new linearized cloud sch.
!     Modified 02-09-30 JJMorcrette PAR, UV, CAPE
!     Modified 03-03-03 M.Ko"hler advection-diffusion PBL
!      Modified 03-12-15 P. Lopez: separate call to simplified convection scheme (CUCALLN2).
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hortal      01-Dec-2003 Intruduce extra-fields coming from the dynamics
!        Y.Tremolet    02-Mar-2004 Check for un-initialised arrays
!        A. Untch  03-2004  EXTRA GFL fields in physics
!        P.Bechtold    11-02-2004 Use simple predictor as SL Physics
!                        and additional first-guess call of cloud scheme,
!                        allow for tracer transport by convection
!        P. Viterbo  24-05-2004  Change surface units
!        M. Ko"hler   3-12-2004  Moist advection-diffusion PBL
!        P. Lopez    28-02-2005  Call VDFMAINS instead of VDFMAIN if LPHYLIN=.T.
!        A. Untch    11-03-2005  Aerosols as named GFL fields
!        K. Yessad and D. Salmond (Feb 2006)  Adapt to LPC_FULL
!        J. Flemming 11--4-2005  Aerosol replaced with reactive gases
!        JJ.Morcrette 2006-02-17 Prognostic aerosols (preliminary v2)
!        P. Bechtold 11-12-2005  Assure positive humidity and reorganize Tracers
!        T.Stockdale / A. Beljaars 31-03-2005  Ocean current b.c.'s
!        M.Janiskova 21-12-2005  Modified conditions for using cloud tendencies
!        D.Salmond   22-Nov-2005 Mods for coarser/finer physics
!        JJMorcrette 20060525    MODIS albedo
!        N. Wedi       06-05-01 phys-dyn coupling
!        M. Ko"hler  6-6-2006 Single Column Model option (LSCMEC)
!        J. Berner   15-Aug-2006 Compute dissipation fields for stochastic physics
!        JJMorcrette 20060721 PP of clear-sky PAR and TOA incident solar radiation
!        S. Serrar   7-9-2006 few tracers added for diagnostics
!        JJMorcrette 20060625 MODIS albedo
!        JJMorcrette 20060925 DU, BC, OM, SO2, SOA climatological fields
!        JJMorcrette 20061002 DDH for aerosol physics
!        G. Balsamo  20070115 Soil type
!        S.Serrar    20070322 physics tendencies for ERA40 are no long
!                      post-processable as extra-fields but as proper GFL fields
!        S.Serrar    17-07-2007 methane introduced
!        A. Tompkins 20070523 Delete conv-cloud iteration on first timestep
!        M. Drusch   21-05-2007  Local arrays for the SEKF / surface analysis
!        A. Geer     25-04-2008  Rainy 4D-Var for multiple sensors
!        S. Serrar   02-05-2008  test on LERA40 removed
!        N. Wedi     08-01-2008  add idealized planet simulations
!        A. Beljaars 27-02-2009  PZIDLWV deleted
!        P. Bechtold 04-10-2008  add call to non-orographic GWD parametrisation
!        M. Leutbecher 09-10-2008  clipping of humidity at initial time (LECLIP*T0)
!        G. Balsamo  08-10-2008  add water holding capacity in the snow-pack
!        A. Geer     01-Oct-2008 rainy 4D-Var now works with physics off; tidied
!        P. Lopez     26-06-2009  Added diagnostic 100m wind components
!        M.Janiskova 12-Mar-2009 not calling SURFTSTP if LPHYLIN=.T.
!        M.Leutbecher 27-02-2009 revised stochastic physics (LSPSDT)
!        P. Lopez     22-10-2008  Ducting diagnostics
!        GMozdzynski/JJMorcrette 20090128 bugfix dynamical extra fields
!        JJMorcrette 20090217 PP of prognostic aerosol diagnsotics and UV processor output fields
!        Y. Takaya   01-Feb-2009 ocean mixed layer model
!        P. de Rosnay 13-02-2009 preliminary passive version of offline Jacobians in surface analysis SEKF
!        JJMorcrette 20091201 Total and clear-sky direct SW radiation flux at surface
!        P. de Rosnay 15-12-2009 Offline Jacobians in surface analysis SEKF
!        J. Munoz Sabater Oct.09 Introduced SMOS data
!        P. Bechtold 07-Oct-2009 enable convective scavenging of tracers
!                                call to diagnostics of diurnal cycle
!        H. Hersbach  04-Dec-2009 10-m neutral wind and friction velocity
!        S. Boussetta/G.Balsamo  05-2009 Add variable LAI fields
!        P.Bechtold/JJMorcrette 11-02-2010 PP of CBASE, 0DEGL and VISIH
!        R. Forbes    07-Apr-2010 cloud / precip fraction for all-sky
!        L. Magnusson 15-feb 2010 Sea-ice (LIM)
!        A. Geer      21-Sep_2010 All-sky AMSU-A
!        P. Lopez     07-Oct-2010 Added obs operator for ground-based radar composites (GBRAD)
!        M. Steinheimer 09-Aug-2010 moved SPBS calculations to subroutine spbsgpupd
!        R. Forbes    01-Mar-2011 Removed code relating to LL3DPRECIPDIAG
!        P. Bechtold  11-Mar-2011 Code cleaning
!        G.Balsamo/S.Boussetta 17-Apr-2011 Added land carbon dioxide
!        J. Hague     21-Mar-2011 YGOM Derived type added
!        M. Ahlgrimm  31-Oct-2011 Add rain,snow and PEXTRA to DDH output
!        M. Ahlgrimm  31-Oct-2011 Clear-sky downward radiation at surface
!        J. Flemming  20-Oct-2011 Call to lightning parameterisation culight
!        J. Flemming  21-Oct-2011 Call to chemical mechanism interface chem_main
!        L. Jones     26-Oct-2011 MACC fluxes no longer passed individually
!        A. Geer      11-Jan-2012 Removed all-sky operator (MWAVE) - it's in HOP now.
!        F. Vana      18-May-2012 Cleaning + few fixes.
!        JJMorcrette  20101125   diagnostic of visibility
!        F. Vana      25-Jul-2012 fix for NANS option
!        JJMorcrette  20120801    Visibility for OPE and MACC
!        N.Semane+P.Bechtold   replace 3600s by RHOUR for small planet
!        F. Vana      05-Dec-2012 rewritten to a bit shorter form
!        JJMorcrette  20130213 PP optical depths GEMS/MACC aerosols
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!        K. Yessad (July 2014): Move some variables.
!        PdeRosnay    201404  SEKF cleaning
!        M. janiskova 10-Jul-2012 Call for simplified surface scheme
!        F. Vana      October-2013 Optimization + fix for simpl. convection
!        M. Ahlgrimm  Apr 2014: write precip fraction for DDH output
!        F. Vana & M. Kharoutdinov 06-Feb-2015: Super-parametrization scheme
!        R. Hogan     14-Nov-2014 Update Tskin each timestep if LApproxLwRadiation
!        P. Lopez     Sept 2015: Modified lightning parameterization.
!        S. Lang      Jan 2016 PPERT2 is passed to SPPT
!        M. Leutbecher & S.-J. Lock (Jan 2016) Introduced SPP scheme (LSPP)
!        SJ Lock      Jan-2016: Enabled independent perturbation patterns (iSPPT)
!        F. Vana      08-Apr-2016  Small fix
!        F. Vana      28-Jun-2016  More consistend iSPPT tendencies
!        R. Hogan     03-Oct-2016  Aerosol diagnostics only every hour
!        SJ Lock      Oct-2016: SPPT option = leave clear-skies radiation UNperturbed
!        F. Vana      19-Nov-2016  Removed useless zeroing of negative q pseudo-flux
!        M Ahlgrimm   2017-11-11 add cloud heterogeneity FSD
!        F. Vana      Jan-2018: 1D model calls radiation in the same way with the 3D case
!        E. Dutra/G.Arduini Jan 2018: snow multi-layer, PSP_SG with 4 dimensions+warm start
!        K. Lonitz    13-Apr-2018  Store PBL and CONV type
!        P. Bechtold  6-Jan-2019 Add clear air turbulence diagnostics
!        B. Ingleby   2019-01-17   Replace PQCFL with Q2M
!        M. Lange     10-Jan-2020 : Adding VARIABLE and FIELD types in EC-physics and surface (IFS-1175)
!        F. Vana      14-Sep-2020 : Updated radiation and input tendency for VDIF.
!        R. Forbes    15-Nov-2020 : Remove first guess cloud satadj (LLCLOUDITER) and TENDENCY_TMP
!        R. Forbes    15-Nov-2020 : Call cloud_satadj just before cloudsc, apply local L/I tend from SLTEND
!        M. Leutbecher   Oct-2020 SPP abstraction
!        S. Massart   March-2021 : Adding BFAS parameters 
!        R. Forbes    15-Oct-2021 : Set up TENDENCY_SATADJ for SLTEND
!     
!-----------------------------------------------------------------------

USE ALGORITHM_STATE_MOD, ONLY : L_OBS_IN_FC
USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0             , ONLY : LIFSMIN, NUNDEFLD, L_OOPS
USE YOMCT3             , ONLY : NSTEP
USE YOMCST             , ONLY : RG       ,RLVTT    ,RLSTT    ,RTT      ,&
 &                              RCPD, RHOUR
USE YOETHF             , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 &                              R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 &                              RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
 &                              RTWAT_RTICE_R      ,RTWAT_RTICECU_R    ,RVTMP2
USE YOECLDP            , ONLY : NCLDQR,NCLDQS,NCLDQI,NCLDQL
USE YOMSEKF            , ONLY : N_SEKF_PT, LUSEKF_REF,&
 &                              FKF_SURF_SO, FKF_SURF_TH, FKF_SURF_CR, FKF_SURF_LR,&
 &                              FKF_TENT, FKF_TENQ, FKF_TENU, FKF_TENV
USE SPP_GEN_MOD        , ONLY : SPP_PERT
USE YOMPHYDER          , ONLY : STATE_TYPE, MASK_GFL_TYPE, DIMENSION_TYPE, AUX_TYPE, SURF_AND_MORE_TYPE,&
 &                              PERTURB_TYPE, MODEL_STATE_TYPE, AUX_RAD_TYPE, FLUX_TYPE, AUX_DIAG_TYPE,&
 &                              DDH_SURF_TYPE, GEMS_LOCAL_TYPE,PERTURB_LOCAL_TYPE,SURF_AND_MORE_LOCAL_TYPE,AUX_DIAG_LOCAL_TYPE,&
 &                              KEYS_LOCAL_TYPE
USE COUPLING
USE YOE_PHYS_MWAVE , ONLY : N_PHYS_MWAVE
USE TM5_CHEM_MODULE    , ONLY : NCHEM2AER
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)           ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELD_VARIABLES)    ,INTENT(INOUT) :: YDVARS
TYPE(TSURF)              ,INTENT(INOUT) :: YDSURF
TYPE(MODEL)              ,INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)    ,INTENT(IN)    :: KDIM
TYPE (AUX_TYPE)          ,INTENT(IN)    :: PAUX
TYPE (AUX_RAD_TYPE)      ,INTENT(INOUT) :: PRAD
TYPE (FLUX_TYPE)         ,INTENT(INOUT) :: FLUX
TYPE (AUX_DIAG_TYPE)     ,INTENT(INOUT) :: PDIAG
TYPE (SURF_AND_MORE_TYPE),INTENT(INOUT) :: PSURF
REAL(KIND=JPRB)          ,INTENT(INOUT) :: PCGPP(KDIM%KLON)
REAL(KIND=JPRB)          ,INTENT(INOUT) :: PCREC(KDIM%KLON)
REAL(KIND=JPRB)          ,INTENT(INOUT) :: PAG(KDIM%KLON)
REAL(KIND=JPRB)          ,INTENT(INOUT) :: PRECO(KDIM%KLON)
TYPE (DDH_SURF_TYPE)     ,INTENT(INOUT) :: PDDHS
TYPE (AUX_DIAG_LOCAL_TYPE),INTENT(INOUT):: AUXL
TYPE (SURF_AND_MORE_LOCAL_TYPE),INTENT(INOUT) :: SURFL
TYPE (KEYS_LOCAL_TYPE)   ,INTENT(INOUT) :: LLKEYS
TYPE (PERTURB_LOCAL_TYPE),INTENT(INOUT) :: PERTL
REAL(KIND=JPRB)          ,INTENT(INOUT) :: PGFL(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM)
TYPE (PERTURB_TYPE)      ,INTENT(INOUT) :: PPERT
TYPE (MODEL_STATE_TYPE)  ,INTENT (IN)   :: PSLPHY9
REAL(KIND=JPRB)          ,INTENT(INOUT) :: PTENGFL(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM1)
REAL(KIND=JPRB)          ,INTENT(INOUT) :: PHYS_MWAVE(KDIM%KLON,KDIM%KLEV,N_PHYS_MWAVE)
REAL(KIND=JPRB)          ,INTENT(INOUT) :: PSAVTEND(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_PHY_G%YRSLPHY%NVTEND)
REAL(KIND=JPRB)          ,INTENT(OUT)   :: PGFLSLP(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_GCONF%YGFL%NDIMSLP)
TYPE (STATE_TYPE)        ,INTENT(INOUT) :: STATE_T0
TYPE (STATE_TYPE)        ,INTENT(INOUT) :: TENDENCY_DYN
TYPE (STATE_TYPE)        ,INTENT(INOUT) :: TENDENCY_CML
TYPE (STATE_TYPE)        ,INTENT(INOUT) :: STATE_TMP
TYPE (STATE_TYPE)        ,INTENT(INOUT) :: TENDENCY_TMP
TYPE (STATE_TYPE)        ,INTENT(INOUT) :: TENDENCY_VDF
TYPE (STATE_TYPE)        ,INTENT(INOUT) :: TENDENCY_SATADJ
TYPE (STATE_TYPE)        ,INTENT(INOUT) :: TENDENCY_LOC
TYPE (STATE_TYPE)        ,INTENT(INOUT) :: TENDENCY_PHY(KDIM%K2DSDT)
!     ------------------------------------------------------------------
LOGICAL :: LLSLPHY, LLRAIN1D, LL_DIAGCLOUDAER, LLISPPT, LL_HRES, LLBUD23
LOGICAL :: LLCLDDIAG=.true.

INTEGER(KIND=JPIM) :: IFLAG, JK, JL, JSW, JEXT, ITRC, J2D, JRF

REAL(KIND=JPRB) :: ZRG, ZRCPD, ZCONS, ZCONS1

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! --------------------------------------
!     Local variables for SPP
! A bunch of SPP variables 
LOGICAL            :: LLPERT_HSDT ! SPP perturbation HSDT on?
INTEGER(KIND=JPIM) :: IPN  ! SPP perturbation pointer
INTEGER(KIND=JPIM) :: IMP  ! SPP random field pointer
TYPE(SPP_PERT)     :: PN1  ! SPP pertn. configs. for HSDT

!     Local arrays for stochastic backscatter
TYPE (GEMS_LOCAL_TYPE) :: GEMSL

! local arrays for the surface analysis (SEKF)
INTEGER(KIND=JPIM)  :: IBLK
REAL(KIND=JPRB)     :: ZU_SAVE_KF(KDIM%KLON,KDIM%KLEV)  ! buffer to save u tendency before moist physics
REAL(KIND=JPRB)     :: ZV_SAVE_KF(KDIM%KLON,KDIM%KLEV)  ! buffer to save v tendency before moist physics
REAL(KIND=JPRB)     :: ZT_SAVE_KF(KDIM%KLON,KDIM%KLEV)  ! buffer to save T tendency before moist physics
REAL(KIND=JPRB)     :: ZQ_SAVE_KF(KDIM%KLON,KDIM%KLEV)  ! buffer to save q tendency before moist physics
!local array for turbulence diagnostics
REAL(KIND=JPRB)     :: ZTENDT_CONV(KDIM%KLON,KDIM%KLEV) ! conv tend for turbulence diag
REAL(KIND=JPRB)     :: ZEDRP(KDIM%KLON,KDIM%KLEV,2)     ! Eddy dissipation rates CAT and MWT
!for SPP
REAL(KIND=JPRB)     :: ZGP2DSPP(KDIM%KLON, YDMODEL%YRML_GCONF%YRSPP_CONFIG%SM%NRFTOTAL)  !SPP pattern

! array for Chemistry to  Aerosol exchange (contains tendencies for gas/wet
! phase reactions
REAL(KIND=JPRB), ALLOCATABLE :: ZCHEM2AER(:,:,:)

!-------------------------------------------------------------

#include "abor1.intfb.h"
#include "local_arrays_ini.intfb.h"
#include "local_arrays_fin.intfb.h"
#include "aerini_layer.intfb.h"
#include "chemini_layer.intfb.h"
#include "compo_apply_emissions_layer.intfb.h"
#include "aer_phy3_layer.intfb.h"
#include "aer_glomapdiag_layer.intfb.h"
#include "cldprg_layer.intfb.h"
#include "clddia_layer.intfb.h"
#include "cloud_s_layer.intfb.h"
#include "cloud_layer.intfb.h"
#include "nocloud.intfb.h"
#include "aer_cloud_layer.intfb.h"
#include "cond_layer.intfb.h"
#include "convection_layer.intfb.h"
#include "convection_s_layer.intfb.h"
#include "noconvection.intfb.h"
#include "noradiation.intfb.h"
#include "nogwdrag.intfb.h"
#include "noturbulence.intfb.h"
#include "state_update.intfb.h"
#include "state_copy.intfb.h"
#include "state_increment.intfb.h"
#include "surfbc_layer.intfb.h"
#include "surftstp_layer.intfb.h"
#include "surftstp_s_layer.intfb.h"
#include "nosurftstp.intfb.h"
#include "gems_tend.intfb.h"
#include "gwdrag_layer.intfb.h"
#include "gwdragwms_layer.intfb.h"
#include "gwdragwms_s_layer.intfb.h"
#include "methox.intfb.h"
#include "o3chem.intfb.h"
#include "qnegat.intfb.h"
#include "qsupersatclip.intfb.h"
#include "climaer_layer.intfb.h"
#include "radflux_layer.intfb.h"
#include "radiation_layer.intfb.h"
#include "radvis_layer.intfb.h"
#include "uvradi_layer.intfb.h"
#include "aerdiag_layer.intfb.h"
#include "satur.intfb.h"
#include "sltend_layer.intfb.h"
#include "stochpert_layer.intfb.h"
#include "surfrad_layer.intfb.h"
#include "turbulence_s_layer.intfb.h"
#include "turbulence_layer.intfb.h"
#include "cuancape2.intfb.h"
#include "vdfvint.intfb.h"
#include "diag_clouds.intfb.h"
#include "diag_dcycle.intfb.h"
#include "ductdia_layer.intfb.h"
#include "convection_ca_layer.intfb.h"
#include "backscatter_layer.intfb.h"
#include "chem_main_layer.intfb.h"
#include "lightning_layer.intfb.h"
#include "update_fields.intfb.h"
#include "gpino3ch.intfb.h"
#include "nemoaddflds_layer.intfb.h"
#include "icestatenemo.intfb.h"
#include "cloud_satadj.intfb.h"
#include "set_ocean_fluxes.intfb.h"
#include "surfws_layer.intfb.h"
#include "diag_turb.intfb.h"

!     ------------------------------------------------------------------

#include "fcttre.func.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CALLPAR',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDSTA=>YDGEOMETRY%YRSTA, YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, &
 & YDLAP=>YDGEOMETRY%YRLAP, YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB, YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB,  YDSPGEOM=>YDGEOMETRY%YSPGEOM, YDEPHLI=>YDMODEL%YRML_PHY_SLIN%YREPHLI, &
 & YDMCC=>YDMODEL%YRML_AOC%YRMCC,YDECUMF=>YDMODEL%YRML_PHY_EC%YRECUMF, YDVDF=>YDMODEL%YRML_PHY_G%YRVDF, &
 & YDERDI=>YDMODEL%YRML_PHY_RAD%YRERDI,YDSLPHY=>YDMODEL%YRML_PHY_G%YRSLPHY,YDEGWD=>YDMODEL%YRML_PHY_EC%YREGWD, &
 & YDSTOPH=>YDMODEL%YRML_PHY_STOCH%YRSTOPH,YDDYNA=>YDMODEL%YRML_DYN%YRDYNA, &
 & YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2,YGFL=>YDMODEL%YRML_GCONF%YGFL, &
 & YDECLDP=>YDMODEL%YRML_PHY_EC%YRECLDP,  &
 & YDECUCONVCA=>YDMODEL%YRML_PHY_EC%YRECUCONVCA,YDEUVRAD=>YDMODEL%YRML_PHY_RAD%YREUVRAD, &
 & YDCUMFS=>YDMODEL%YRML_PHY_SLIN%YRCUMFS,  &
 & YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO, &
 & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM,YDPHNC=>YDMODEL%YRML_PHY_SLIN%YRPHNC,  &
 & YDEGWWMS=>YDMODEL%YRML_PHY_EC%YREGWWMS, YDOZO=>YDMODEL%YRML_CHEM%YROZO, &
 & YDCHEM=>YDMODEL%YRML_CHEM%YRCHEM,YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY,YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD,  &
 & YDNCL=>YDMODEL%YRML_PHY_SLIN%YRNCL, YDSPPT=>YDMODEL%YRML_SPPT, YDSPPT_CONFIG=>YDMODEL%YRML_GCONF%YRSPPT_CONFIG, &
 & YDSPP_CONFIG=>YDMODEL%YRML_GCONF%YRSPP_CONFIG) ! YDSPP=>YDMODEL%YRML_SPP, needed?

ASSOCIATE(NACTAERO=>YGFL%NACTAERO, NAERO=>YGFL%NAERO, NCHEM=>YGFL%NCHEM, &
 & LAERCHEM=>YGFL%LAERCHEM, LAERNITRATE=>YDMODEL%YRML_CHEM%YRCOMPO%LAERNITRATE, &
 & NCHEM_DV=>YGFL%NCHEM_DV, NDIM=>YGFL%NDIM, NDIM1=>YGFL%NDIM1, &
 & NDIMSLP=>YGFL%NDIMSLP, NGEMS=>YGFL%NGEMS, NGFL_EXT=>YGFL%NGFL_EXT, &
 & NGHG=>YGFL%NGHG, YA=>YGFL%YA, YTKE=>YGFL%YTKE,&
 & YAERO=>YGFL%YAERO, YEXT=>YGFL%YEXT, YGHG=>YGFL%YGHG, YI=>YGFL%YI, &
 & YL=>YGFL%YL, YNOGW=>YGFL%YNOGW, NNOGW=>YGFL%NNOGW, YO3=>YGFL%YO3, YEDRP=>YGFL%YEDRP, NEDRP=>YGFL%NEDRP,&
 & YQ=>YGFL%YQ, YR=>YGFL%YR, YS=>YGFL%YS, YUVP=>YGFL%YUVP, &
 & LCHEM_LIGHT=>YDCHEM%LCHEM_LIGHT, &
 & LECUMFS=>YDCUMFS%LECUMFS, &
 & NVEXTR=>YDDPHY%NVEXTR, NVEXTRDYN=>YDDPHY%NVEXTRDYN, &
 & LAERDIAG1=>YDEAERATM%LAERDIAG1, LAERINIT=>YDEAERATM%LAERINIT, &
 & NAERCLD=>YDECLDP%NAERCLD, &
 & LCUCONV_CA=>YDECUCONVCA%LCUCONV_CA, NLIVES=>YDECUCONVCA%NLIVES, &
 & LMFCUCA=>YDECUMF%LMFCUCA, NJKT3=>YDECUMF%NJKT3, NJKT4=>YDECUMF%NJKT4, &
 & LDIAG_STRATO=>YDEGWD%LDIAG_STRATO, &
 & GTPHYGWWMS=>YDEGWWMS%GTPHYGWWMS, &
 & LPHYLIN=>YDEPHLI%LPHYLIN, LPHYSFCLIN=>YDEPHLI%LPHYSFCLIN, LERADIMPL=>YDEPHY%LERADIMPL, &
 & LBUD23=>YDEPHY%LBUD23, LBUDCYCLE=>YDEPHY%LBUDCYCLE, LDIAGTURB_EC=>YDEPHY%LDIAGTURB_EC,&
 & LDUCTDIA=>YDEPHY%LDUCTDIA, LECLIPCLDT0=>YDEPHY%LECLIPCLDT0, &
 & LVEXTRDYNACC=>YDEPHY%LVEXTRDYNACC, LEXTRATEND=>YDEPHY%LEXTRATEND, &
 & LECLIPQT0=>YDEPHY%LECLIPQT0, LECOND=>YDEPHY%LECOND, LECUMF=>YDEPHY%LECUMF, &
 & LEDCLD=>YDEPHY%LEDCLD, LEGWDG=>YDEPHY%LEGWDG, LEGWWMS=>YDEPHY%LEGWWMS, &
 & LEMETHOX=>YDEPHY%LEMETHOX, LEO3CH=>YDEPHY%LEO3CH, LEPCLD=>YDEPHY%LEPCLD, &
 & LEQNGT=>YDEPHY%LEQNGT, LERADI=>YDEPHY%LERADI, LERADS=>YDEPHY%LERADS, LSLPHY=>YDEPHY%LSLPHY,&
 & LERAIN=>YDEPHY%LERAIN, LESURF=>YDEPHY%LESURF, LEVDIF=>YDEPHY%LEVDIF, &
 & LEVDIFSL=>YDEPHY%LEVDIFSL, LEMWAVE=>YDEPHY%LEMWAVE, &
 & LAERVISI=>YDERAD%LAERVISI, LAPPROXLWUPDATE=>YDERAD%LAPPROXLWUPDATE, &
 & LECSRAD=>YDERAD%LECSRAD, NRADFR=>YDERAD%NRADFR, NTSW=>YDERAD%NTSW, &
 & RCCNLND=>YDERAD%RCCNLND, RCCNSEA=>YDERAD%RCCNSEA, &
 & REPCLC=>YDERDI%REPCLC, &
 & LUVPROC=>YDEUVRAD%LUVPROC, NRADUV=>YDEUVRAD%NRADUV, &
 & LNEMOLIMPUT=>YDMCC%LNEMOLIMPUT, LNEMOLIMTHK=>YDMCC%LNEMOLIMTHK, &
 & LNEMOATMFLDS=>YDMCC%LNEMOATMFLDS, &
 & LNCLIN=>YDNCL%LNCLIN, &
 & LENCLD2=>YDPHNC%LENCLD2, LEPCLD2=>YDPHNC%LEPCLD2, &
 & NSTART=>YDRIP%NSTART, RSLWX=>YDSLPHY%RSLWX, &
 & NVTEND=>YDSLPHY%NVTEND, &
 & LEXTRAFIELDS=>YDSTOPH%LEXTRAFIELDS, LFORCENL=>YDSTOPH%LFORCENL, &
 & LSTOPH_CASBS=>YDSTOPH%LSTOPH_CASBS, &
 & LSTOPH_SPBS=>YDSTOPH%LSTOPH_SPBS, LVORTCON=>YDSTOPH%LVORTCON, &
 & NFORCEEND=>YDSTOPH%NFORCEEND, NFORCESTART=>YDSTOPH%NFORCESTART, &
 & YSD_VD=>YDSURF%YSD_VD, YSD_VF=>YDSURF%YSD_VF, YSD_VN=>YDSURF%YSD_VN, &
 & YSP_RR=>YDSURF%YSP_RR, YSP_SG=>YDSURF%YSP_SG, &
 & AERO_SCHEME=>YDCOMPO%AERO_SCHEME, LCLDBUD_TIMEINT=>YDECLDP%LCLDBUD_TIMEINT, &
 & LAERSOA=>YDCOMPO%LAERSOA, &
 & TSPHY=>YDPHY2%TSPHY, LELIGHT=>YDEPHY%LELIGHT, LESNML=>YDEPHY%LESNML,&
 & NSNMLWS=>YDEPHY%NSNMLWS)
!     ------------------------------------------------------------------

!*         0.     INITIALIZATION
!
!          0.1    Initialization of constants

IF(L_OOPS) THEN
  LL_HRES =  L_OBS_IN_FC() .OR. YDEPHY%LFPOS_EC_PHYS
ELSE
  LL_HRES = .NOT. LIFSMIN
ENDIF

LLSLPHY = LSLPHY.AND.YDDYNA%LSLAG
LLRAIN1D = .FALSE. ! key reserved for 1D var rain

IF (YDSPPT_CONFIG%LSPSDT.AND.(.NOT. YDSPPT_CONFIG%LSPPT1)) THEN
  LLISPPT=.TRUE. !  NOT "Standard SPPT", (iSPPT => MPSDT.NE.111111)
ELSE
  LLISPPT=.FALSE.
ENDIF

!LBUD23: storing tendencies
LLBUD23 = LBUD23.OR.LEXTRATEND

! constants
ZRG=1.0_JPRB/RG
ZRCPD=1.0_JPRB/RCPD

! Some sort of security (I wonder why this is not in the setup)
IF (LEPCLD .AND. .NOT.(YA%LACTIVE.AND.YL%LACTIVE.AND.YI%LACTIVE)) THEN
  CALL ABOR1('CALLPAR: How did we get here?')
ENDIF

!          0.2.  Initialization of arrays

! Create the space for local structures, associate the pointers and initialize them
CALL LOCAL_ARRAYS_INI(YDGEOMETRY,YDSURF,YDMODEL,KDIM, LLKEYS, PAUX, AUXL, SURFL, PERTL, GEMSL,&
 & PRAD,PSURF,PGFL,PTENGFL)

IF (.NOT.LEO3CH) THEN
  TENDENCY_CML%O3=0.0_JPRB
ENDIF

! If time integrate is false, then reset PEXTRA so instantaneous values are stored
IF (.NOT. LCLDBUD_TIMEINT) PSURF%PSD_XA(:,:,:)=0.0_JPRB

! Initialize tendency_cml (being simple copy of tendency_dyn at the moment)
CALL STATE_COPY(KDIM,TENDENCY_DYN,TENDENCY_CML)

!     ------------------------------------------------------------------

!*         1.     INITIAL COMPUTATIONS.
!                 ---------------------

!*       1.1   FIRST TIME-STEP
IF (NSTEP == NSTART) THEN
  !*        Clip specific humidity and cloud variables at initial time
  !            (this is specially important for radiation scheme)
  IF (LECLIPCLDT0) THEN
    DO JK=1,KDIM%KLEV
      DO JL=KDIM%KIDIA,KDIM%KFDIA
        STATE_T0%A(JL,JK)=MIN(MAX(STATE_T0%A(JL,JK),REPCLC),1.0_JPRB-REPCLC)
        STATE_T0%CLD(JL,JK,NCLDQL)=MAX( STATE_T0%CLD(JL,JK,NCLDQL), 0.0_JPRB )
        STATE_T0%CLD(JL,JK,NCLDQI)=MAX( STATE_T0%CLD(JL,JK,NCLDQI), 0.0_JPRB )
        STATE_T0%CLD(JL,JK,NCLDQR)=MAX( STATE_T0%CLD(JL,JK,NCLDQR), 0.0_JPRB )
        STATE_T0%CLD(JL,JK,NCLDQS)=MAX( STATE_T0%CLD(JL,JK,NCLDQS), 0.0_JPRB )
      ENDDO
    ENDDO
  ENDIF

  !   clip supersaturation and negative humidity
  IF (LECLIPQT0) THEN
    !   clip supersaturation
    CALL QSUPERSATCLIP(YDECLDP,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON, KDIM%KLEV, STATE_T0%T, STATE_T0%A, PAUX%PAPRSF, STATE_T0%Q)
    !   clip negative humidity
    CALL QNEGAT (&
       & YDMODEL%YRML_PHY_EC%YRECND, &
       & KDIM%KIDIA , KDIM%KFDIA , KDIM%KLON , KDIM%KLEV,&
       & TSPHY, STATE_T0%Q , TENDENCY_LOC%Q, PAUX%PAPRS,&
       ! FLUX OUTPUTS, N.B. not used presently
       & FLUX%PFCQNG )

    CALL UPDATE_FIELDS(YDPHY2,1,KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV,&
       & PTA1=TENDENCY_LOC%Q, PO1=STATE_T0%Q)
  ENDIF

  ! set reasonable defaults for pseudohistoric variables
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    PSURF%PSD_VN(JL,YSD_VN%YTOP%MP)=1.0_JPRB
    PSURF%PSD_VN(JL,YSD_VN%YBAS%MP)=1.0_JPRB
    PSURF%PSD_VN(JL,YSD_VN%YACPR%MP)=0.0_JPRB
    PSURF%PSD_VN(JL,YSD_VN%YACCPR%MP)=0.0_JPRB
    PSURF%PSD_VD(JL,YSD_VD%Y10FGCV%MP)=0.0_JPRB
    PSURF%PSD_VD(JL,YSD_VD%YZ0F%MP)=PSURF%PSD_VF(JL,YSD_VF%YZ0F%MP)
    PSURF%PSD_VD(JL,YSD_VD%YLZ0H%MP)=PSURF%PSD_VF(JL,YSD_VF%YLZ0H%MP)
  ENDDO
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    PSURF%PSD_VN(JL,YSD_VN%YACPR%MP)=0.0_JPRB
    AUXL%IBASC(JL)=1
    AUXL%ITOPC(JL)=1
  ENDDO

ENDIF

!   Various  initializations

!   Initialization of updated state state_tmp
IF (LPHYLIN) THEN
  ! Parallel physics - doing just coppy of state_t0 to state_tmp
  CALL STATE_COPY(KDIM,STATE_T0,STATE_TMP)
ELSE
  ! Sequential physics, at the moment state_tmp updated by dynamics only.
  CALL STATE_UPDATE(YDPHY2,KDIM,STATE_T0,TENDENCY_CML,STATE_TMP)
ENDIF


! liquid water loading in environment neglected in convection
IF ((.NOT.LLSLPHY).AND.(YL%LACTIVE).AND.(YI%LACTIVE)) THEN
  AUXL%ZLISUM(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)=STATE_T0%CLD(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,NCLDQL)&
                                              & +STATE_T0%CLD(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,NCLDQI)
ELSE
  AUXL%ZLISUM(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)=0.0_JPRB
ENDIF

!  Security for GEMS
IF ((NGEMS > 0).AND.(NAERO > 0)) THEN
  DO JEXT=1,NAERO
    DO JK=1,KDIM%KLEV
      DO JL=KDIM%KIDIA,KDIM%KFDIA
        GEMSL%ZAEROP(JL,JK,JEXT) = MAX(0._JPRB, YDVARS%AERO(JEXT)%PH9(JL,JK))
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF(( NVEXTRDYN>0 .AND. NUNDEFLD /=1 ).OR.( NVEXTRDYN>1 .AND. NUNDEFLD==1 )) THEN
  IF (LVEXTRDYNACC) THEN
    PSURF%PSD_XA(KDIM%KIDIA:KDIM%KFDIA,:,NVEXTR-NVEXTRDYN+1:NVEXTR)= &
      & PSURF%PSD_XA(KDIM%KIDIA:KDIM%KFDIA,:,NVEXTR-NVEXTRDYN+1:NVEXTR) + &
      &    PSURF%PEXTRD(KDIM%KIDIA:KDIM%KFDIA,:,1:NVEXTRDYN)
  ELSE
    PSURF%PSD_XA(KDIM%KIDIA:KDIM%KFDIA,:,NVEXTR-NVEXTRDYN+1:NVEXTR)=&
      & PSURF%PEXTRD(KDIM%KIDIA:KDIM%KFDIA,:,1:NVEXTRDYN)
  ENDIF
ENDIF

IF ( (NCHEM > 0 .OR. NGEMS > 0 .OR. NACTAERO > 0) .AND. LL_HRES ) THEN
  ! ZCHEM2AER will be passed as an argument (even if unused)
  IF (NACTAERO > 0 .AND. NCHEM > 0 .AND. (LAERCHEM .OR. LAERNITRATE .OR. LAERSOA) ) THEN
    ! ZCHEM2AER will actually be used
    ! 1 SO4 tendency in kg/(kg*s)
    ! 2 SO2 tendency due to OH  kg/(kg*s)
    ! 3 SO2 tendency total chemistry  kg/(kg*s)
    ! 4 NH4 tendency total chemistry  kg/(kg*s)
    ! 5 SOA 1 tendency from SOG (kg/kgs)
    ! 6 SOA 2 tendency from SOG (kg/kgs)
    ALLOCATE( ZCHEM2AER(KDIM%KLON,KDIM%KLEV,NCHEM2AER) )
    ZCHEM2AER(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,1:NCHEM2AER)=0._JPRB
  ELSE
    ALLOCATE( ZCHEM2AER(1,1,1) )
  ENDIF
ENDIF


!*       1.2   CHANGE UNITS/SCALINGS
!*                and stochastic perturbation
!  Prepare SPP
IF (YDSPP_CONFIG%LSPP) THEN
  IPN = YDSPP_CONFIG%PPTR%HSDT
  LLPERT_HSDT = IPN > 0
  IF (LLPERT_HSDT) THEN
    PN1     = YDSPP_CONFIG%SM%PN(IPN)
    IMP     = PN1%MP
  ENDIF
ELSE
  LLPERT_HSDT  =.FALSE.
ENDIF
SURFL%ZSNM1M(KDIM%KIDIA:KDIM%KFDIA,:) = PSURF%PSP_SG(KDIM%KIDIA:KDIM%KFDIA,:,YSP_SG%YF%MP9)
SURFL%ZRSNM1M(KDIM%KIDIA:KDIM%KFDIA,:)= PSURF%PSP_SG(KDIM%KIDIA:KDIM%KFDIA,:,YSP_SG%YR%MP9)
DO JL=KDIM%KIDIA,KDIM%KFDIA
  GEMSL%ZAZ0M(JL)  = PSURF%PSD_VD(JL,YSD_VD%YZ0F%MP) /RG
  GEMSL%ZAZ0H(JL)  = EXP(PSURF%PSD_VD(JL,YSD_VD%YLZ0H%MP))
  SURFL%ZWLM1M(JL) = PSURF%PSP_RR(JL,YSP_RR%YW%MP9)
  IF (PSURF%PSD_VF(JL,YSD_VF%YLSM%MP) > 0.5_JPRB) THEN
    IF (LLPERT_HSDT) THEN
      !   
      !   perturbation of standard deviation of subgrid orography
      !
      PSURF%PHSTD (JL) = PSURF%PSD_VF(JL,YSD_VF%YGETRL%MP) &
           &*EXP( PN1%MU(1) +PN1%XMAG(1)*PPERT%PGP2DSPP(JL,1,IMP)) /RG
    ELSE
      !
      !   unperturbed standard deviation of subgrid orography
      !
      PSURF%PHSTD (JL) = PSURF%PSD_VF(JL,YSD_VF%YGETRL%MP) /RG
    ENDIF
    SURFL%ZHSDFOR(JL) = PSURF%PSD_VF(JL,YSD_VF%YSDFOR%MP)
  ELSE
    PSURF%PHSTD (JL) = 0._JPRB
    PSURF%PSD_VF(JL,YSD_VF%YSIG%MP) = 0._JPRB
    SURFL%ZHSDFOR(JL) = 0._JPRB
  ENDIF
ENDDO
!*      1.3   AUXILLIARY VARIABLES FOR VDF, SRF AND SURFRAD
CALL SURFBC_LAYER(YDSURF,YDEPHY,KDIM,PAUX,LLKEYS,PSURF,SURFL)

!*      1.3b   INITIALIZE MULTI-LAYER SNOW FIELDS WITH PARAMETRIZED PROFILES
IF ( (NSNMLWS > 0) .AND. (NSTEP==NSTART) .AND. (LESNML) .AND. (KDIM%KLEVSN > 1) ) THEN
  CALL SURFWS_LAYER(YDSURF,YDEPHY,KDIM,PSURF,SURFL,PAUX)
ENDIF

!*         1.6    ATMOSPHERIC COMPOSITION - PRESCRIBED EMISSIONS
!                 ----------------------------------------------

IF (NACTAERO /= 0 .OR. NCHEM /= 0 .OR. NGHG /= 0) THEN
  CALL COMPO_APPLY_EMISSIONS_LAYER(YDSURF,YDMODEL,KDIM,PAUX,PGFL,PSURF,GEMSL)
ENDIF

!*         1.7    PROGNOSTIC CHEMISTRY - INITIAL COMPUTATIONS (CHANGE FLUXES)
!                 ------------------------------------------

IF (NCHEM /= 0) THEN
  CALL CHEMINI_LAYER(YDSURF,YDMODEL,KDIM,PAUX,STATE_TMP,PSURF,SURFL,GEMSL)
ENDIF

!*         1.8    RADIATION TRANSFER - INITIAL COMPUTATIONS
!                 -----------------------------------------

IF (NACTAERO /= 0 .AND. LL_HRES) THEN
  IF (.NOT.LAERINIT) THEN
    !CALL SURFRAD_LAYER(KDIM, PAUX, state_t0, LLKEYS, AUXL, PSURF, SURFL) ! consistent with oper. radiation
    CALL SURFRAD_LAYER(YDSURF,YDMCC,YDERAD,YDEPHY,YDRIP,KDIM,PAUX,STATE_TMP,LLKEYS,AUXL,PSURF,SURFL)

    !*         1.9    PROGNOSTIC AEROSOLS - INITIAL COMPUTATIONS
    CALL AERINI_LAYER(YDGEOMETRY,YDSURF,YDMODEL,KDIM,PAUX,STATE_TMP,YDVARS%R, YDVARS%S, PSURF,SURFL,GEMSL)
  ELSE
    GEMSL%ZCFLX(KDIM%KIDIA:KDIM%KFDIA, GEMSL%IAERO(1):GEMSL%IAERO(NACTAERO))=0._JPRB
    GEMSL%ZTENC(KDIM%KIDIA:KDIM%KFDIA, 1:KDIM%KLEV, GEMSL%IAERO(1):GEMSL%IAERO(NACTAERO))=0._JPRB
  ENDIF
ENDIF

!       Full-Budget-Run: adiabatic tendencies

! Here tendency_cml = tendency_dyn but the former are secured also for YQ%LACTIVE=.F.
IF (LLBUD23) CALL UPDATE_FIELDS(YDPHY2,1,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,&
  & PTA1=TENDENCY_CML%U, PO1=PSURF%PSD_XA(:,:,1),&
  & PTA2=TENDENCY_CML%V, PO2=PSURF%PSD_XA(:,:,2),&
  & PTA3=TENDENCY_CML%T, PO3=PSURF%PSD_XA(:,:,3),&
  & PTA4=TENDENCY_CML%Q, PO4=PSURF%PSD_XA(:,:,4) )

!     ------------------------------------------------------------------

!*         2.     RADIATION TRANSFER
!                 ------------------

! Cloud and aerosol diagnostics typically only needed every hour
IF(MOD(NSTEP*TSPHY,RHOUR)==0.0_JPRB.OR.NSTEP==NSTART)  THEN
  LL_DIAGCLOUDAER=.TRUE.
ELSE
  LL_DIAGCLOUDAER=.FALSE.
ENDIF


IF ( LERADI ) THEN

!*         2.1  FULL RADIATION COMPUTATIONS

!*         2.2 SURFACE BOUNDARY CONDITIONS
  IF (LERADS) THEN
    !CALL SURFRAD_LAYER(KDIM, PAUX, state_t0,LLKEYS, AUXL, PSURF, SURFL) ! this one is consistent with oper. radiation
    CALL SURFRAD_LAYER(YDSURF,YDMCC,YDERAD,YDEPHY,YDRIP,KDIM,PAUX,STATE_TMP,LLKEYS,AUXL,PSURF,SURFL)
  ELSE
    DO JSW=1,NTSW
      SURFL%ZALBD(KDIM%KIDIA:KDIM%KFDIA,JSW)=PSURF%PSD_VF(KDIM%KIDIA:KDIM%KFDIA,YSD_VF%YALBF%MP)
      SURFL%ZALBP(KDIM%KIDIA:KDIM%KFDIA,JSW)=PSURF%PSD_VF(KDIM%KIDIA:KDIM%KFDIA,YSD_VF%YALBF%MP)
    ENDDO
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      PSURF%PSD_VD(JL,YSD_VD%YALB%MP)=PSURF%PSD_VF(JL,YSD_VF%YALBF%MP)
      PSURF%PEMIS(JL)=PSURF%PSD_VF(1,YSD_VF%YEMISF%MP)
      SURFL%ZEMIR(JL)=PSURF%PSD_VF(1,YSD_VF%YEMISF%MP)
      SURFL%ZEMIW(JL)=PSURF%PSD_VF(1,YSD_VF%YEMISF%MP)
      AUXL%ZCCNL(JL)=RCCNLND
      AUXL%ZCCNO(JL)=RCCNSEA
    ENDDO
  ENDIF

  IFLAG=2
  CALL SATUR (KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KTDIA, KDIM%KLEV, YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN, &
   & PAUX%PAPRSF, STATE_T0%T, PDIAG%PQSAT, IFLAG)

  !*         2.3 CALCULATE CLOUD COVER DIAGNOSTICS
  !              (every hour or radiation timestep)
  IF ( (.NOT.YDDYNA%LTWOTL.AND.(NSTEP == 1))&
     & .OR. MOD(NSTEP*TSPHY,RHOUR)==0.0_JPRB .OR. MOD(NSTEP,NRADFR) == 0 .OR. LPHYLIN )    THEN

    ! Calculate cloud cover (total/high/medium/low) diagnostics
    ! and perform numerical security checks for radiation scheme
    IF (LEDCLD) THEN
      ! preliminar modification before prognostic scheme is again called
      ! during NL computation of 4D-Var
      !IF((LEPCLD) .OR. (LENCLD2)) THEN
      IF(LEPCLD .OR. ((LENCLD2.OR.LEPCLD2) .AND. (LNCLIN .OR. (NSTEP /= NSTART)))) THEN
        IF(LERADIMPL) THEN
          CALL CLDPRG_LAYER(YDSURF,YDMODEL,KDIM,PAUX,STATE_TMP,PSURF,PRAD)
        ELSE
          CALL CLDPRG_LAYER(YDSURF,YDMODEL,KDIM,PAUX,STATE_T0,PSURF,PRAD)
        ENDIF
      ELSE
        IF(LERADIMPL) THEN
          CALL CLDDIA_LAYER(YDSURF,YDEPHLI,YDMODEL%YRML_PHY_EC%YRECLD,KDIM,PAUX,STATE_TMP,AUXL,PDIAG,PSURF,PRAD)
        ELSE
          CALL CLDDIA_LAYER(YDSURF,YDEPHLI,YDMODEL%YRML_PHY_EC%YRECLD,KDIM,PAUX,STATE_T0,AUXL,PDIAG,PSURF,PRAD)
        ENDIF
      ENDIF
    ELSE
      DO JK=1,KDIM%KLEV
        DO JL=KDIM%KIDIA,KDIM%KFDIA
          PRAD%PNEB(JL,JK)=0.0_JPRB
          PRAD%PQLI(JL,JK)=0.0_JPRB
          PRAD%PQICE(JL,JK)=0.0_JPRB
        ENDDO
      ENDDO
    ENDIF

  ENDIF

  !*        2.4 CALL RADIATION AT EVERY GRID POINT
  !             (every radiation timestep)
  IF ( (.NOT.YDDYNA%LTWOTL.AND.(NSTEP == 1))&
     & .OR. MOD(NSTEP,NRADFR) == 0 .OR. LPHYLIN )    THEN

     ! Call radiation
     IF (LPHYLIN) THEN
      CALL RADIATION_LAYER(YDGEOMETRY%YRDIMV,YDSURF,YDMODEL,KDIM,STATE_T0,PDIAG,PRAD,PAUX,AUXL,PSURF,SURFL)
    ENDIF

    ! Store the skin T seen by the full radiation timestep
    IF (.NOT. LPHYLIN) PRAD%PEDRO(KDIM%KIDIA:KDIM%KFDIA)=PSURF%PSP_RR(KDIM%KIDIA:KDIM%KFDIA,YSP_RR%YT%MP9)

  ELSE
    DO JK=1,KDIM%KLEV
      DO JL=KDIM%KIDIA,KDIM%KFDIA
        PRAD%PNEB(JL,JK)=0.0_JPRB
        PRAD%PQLI(JL,JK)=0.0_JPRB
        PRAD%PQICE(JL,JK)=0.0_JPRB
      ENDDO
    ENDDO
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      FLUX%PFRSOC(JL,0)=0.0_JPRB
      FLUX%PFRTHC(JL,0)=0.0_JPRB
      FLUX%PFRSOC(JL,1)=0.0_JPRB
      FLUX%PFRTHC(JL,1)=0.0_JPRB
    ENDDO
  ENDIF

  IF (LPHYLIN .OR. LAPPROXLWUPDATE) THEN
    ! Use the most recent skin temperature.  Note that if we are
    ! updating the longwave fluxes every timestep (LApproxLwUpdate)
    ! then it is necessary to update the skin temperature here,
    ! otherwise the surface scheme will linearise around an old skin
    ! temperature that is not consistent with the updated surface
    ! longwave net fluxes
    AUXL%ZTSKRAD(KDIM%KIDIA:KDIM%KFDIA)=PSURF%PSP_RR(KDIM%KIDIA:KDIM%KFDIA,YSP_RR%YT%MP9)
  ELSE
    ! Use the skin T from the last full radiation timestep
    AUXL%ZTSKRAD(KDIM%KIDIA:KDIM%KFDIA)=PRAD%PEDRO(KDIM%KIDIA:KDIM%KFDIA)
  ENDIF

  !*         2.5  RADIATIVE FLUXES AT EVERY TIME-STEP
  CALL RADFLUX_LAYER(YDSURF,YDMODEL,YDSPP_CONFIG,KDIM,PAUX, &
       &             STATE_T0,STATE_TMP,PSURF,SURFL,PRAD, AUXL,FLUX,TENDENCY_LOC,PPERT)

  !* AB: extra variables (4, 3D) to check instantaneous heating rates
  !* 4 new 3D extra variables on model levels need to be set in prepifs
  IF (LECSRAD) THEN
    IF (LEXTRATEND) THEN
      !Store accumulated heating rates in extra fields (30:33)
      CALL UPDATE_FIELDS(YDPHY2,1,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,&
      & PTA1=PRAD%PHRSW, PO1=PSURF%PSD_XA(:,:,30),&
      & PTA2=PRAD%PHRSC, PO2=PSURF%PSD_XA(:,:,31),&
      & PTA3=PRAD%PHRLW, PO3=PSURF%PSD_XA(:,:,32),&
      & PTA4=PRAD%PHRLC, PO4=PSURF%PSD_XA(:,:,33) )

    ELSE
      DO JK=1,KDIM%KLEV
        DO JL=KDIM%KIDIA,KDIM%KFDIA
          PSURF%PSD_XA(JL,JK,1)=PRAD%PHRSW(JL,JK)*86400.0_JPRB
          PSURF%PSD_XA(JL,JK,2)=PRAD%PHRSC(JL,JK)*86400.0_JPRB
          PSURF%PSD_XA(JL,JK,3)=PRAD%PHRLW(JL,JK)*86400.0_JPRB
          PSURF%PSD_XA(JL,JK,4)=0.0_JPRB !PRAD%PHRLC(JL,JK)*86400.0_JPRB
        ENDDO
      ENDDO

      DO JL=KDIM%KIDIA,KDIM%KFDIA
        PSURF%PSD_XA(JL,1,4)=FLUX%PFRSOD(JL) !solar down global
        PSURF%PSD_XA(JL,2,4)=PRAD%PFDIR(JL) !solar direct horiz.
        PSURF%PSD_XA(JL,3,4)=FLUX%PFRSODC(JL)
        PSURF%PSD_XA(JL,4,4)=PAUX%PMU0(JL) !solar zenith angle
        PSURF%PSD_XA(JL,5,4)=PRAD%PISUND(JL) !sunshine duration
        PSURF%PSD_XA(JL,6,4)=PRAD%PDSRP(JL)  !direct orthogonal sw
      ENDDO
    ENDIF  !LEXTRATEND
  ENDIF  !LECSRAD



  IF (LAERVISI .AND. LL_DIAGCLOUDAER) THEN
    CALL CLIMAER_LAYER(YDSURF,YDMODEL,KDIM,PAUX,STATE_T0,STATE_TMP,PSURF,GEMSL)
  ENDIF

  ! -- UV spectral fluxes at the surface at full radiation steps
  IF (MOD(NSTEP,NRADUV) == 0 .AND. LUVPROC) THEN
    CALL UVRADI_LAYER(YDSURF,YDMODEL,KDIM, PAUX, AUXL, GEMSL, SURFL, &
      & STATE_T0, STATE_TMP, PDIAG, PSURF, YDVARS%UVP(1), YDVARS%CHEM)
  ENDIF

  !       Full-Budget-Run: radiation
  IF (LLBUD23) CALL UPDATE_FIELDS(YDPHY2,1,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,&
      & PTA1=TENDENCY_LOC%T,PO1=PSURF%PSD_XA(:,:,5))

  ! Update of tendency_cml after radiation
  CALL STATE_INCREMENT(KDIM,TENDENCY_CML, PA_T=TENDENCY_LOC%T)

  !iSPPT: Store radiation tendencies for stochastic perturbations
  IF (LLISPPT) THEN
    J2D = YDSPPT%MPSDT(1)   !pattern ID for RADiation tendency perturbations
    ! .EQ.0 => do NOT perturb tendency
    IF (J2D > 0) CALL STATE_INCREMENT(KDIM, TENDENCY_PHY(J2D), PA_T=TENDENCY_LOC%T)
  ENDIF

ELSE
  ! By-passed radiation settings
  CALL NORADIATION(YDSURF,KDIM, PRAD, FLUX, AUXL, PSURF, SURFL)
  IF (LLSLPHY.AND.LEVDIFSL) TENDENCY_LOC%T(:,:)=0._JPRB
ENDIF

!-- additional aerosol diagnostics during MACC forecasts
IF (LAERDIAG1) THEN
  CALL AERDIAG_LAYER(YDSURF,YDECLDP,YDERAD,YDMODEL%YRML_GCONF,YDPHY2,KDIM,GEMSL,PAUX,STATE_T0,PDIAG,AUXL,PSURF)
ENDIF


!*         2.8 SAVE NET RADIATION AT THE SURFACE FOR THE SURFACE ANALYSIS (SEKF)

IF (LUSEKF_REF) THEN
! save net radiation at the surface from the unperturbed first sEKF model run
  IF (N_SEKF_PT==0) THEN
    IBLK=(KDIM%KSTGLO-1)/KDIM%KLON + 1
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      FKF_SURF_SO(JL,IBLK,NSTEP) = FLUX%PFRSO(JL,KDIM%KLEV)
      FKF_SURF_TH(JL,IBLK,NSTEP) = FLUX%PFRTH(JL,KDIM%KLEV)
    ENDDO
  ENDIF
! replace radiation with the unperturbed first sEKF values
  IF (N_SEKF_PT>0) THEN
    IBLK=(KDIM%KSTGLO-1)/KDIM%KLON + 1
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      FLUX%PFRSO(JL,KDIM%KLEV) = FKF_SURF_SO(JL,IBLK,NSTEP)
      FLUX%PFRTH(JL,KDIM%KLEV) = FKF_SURF_TH(JL,IBLK,NSTEP)
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*         3.     VERTICAL EXCHANGE OF U,V,T,Q BY TURBULENCE AND
!*                 ----------------------------------------------
!*                     GRAVITY WAVE DRAG PARAMETERISATION

! Create tendecy to feed VD (containing the 0.5 of SLPHYS from previous timestep)
CALL STATE_COPY(KDIM,TENDENCY_CML,TENDENCY_TMP)
IF (LLSLPHY .AND. LEVDIFSL .AND.  NSTEP /= NSTART) THEN
  ZCONS=1.0_JPRB/TSPHY
  DO JK=1,KDIM%KLEV
    ZCONS1=(1._JPRB-RSLWX(JK))*ZCONS
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      ! Add the previous timestep tendency physics
      TENDENCY_TMP%U (JL,JK) = TENDENCY_TMP%U (JL,JK) + ZCONS1*PSLPHY9%U(JL,JK)
      TENDENCY_TMP%V (JL,JK) = TENDENCY_TMP%V (JL,JK) + ZCONS1*PSLPHY9%V(JL,JK)
      TENDENCY_TMP%T (JL,JK) = TENDENCY_TMP%T (JL,JK) + ZCONS1*PSLPHY9%T(JL,JK) &
        &  - (1._JPRB-RSLWX(JK))*TENDENCY_LOC%T (JL,JK)  ! Subtract part of the actual radiation
      IF (YQ%LPHY) THEN
        TENDENCY_TMP%Q (JL,JK) = TENDENCY_TMP%Q (JL,JK) + ZCONS1*PSLPHY9%GFL(JL,JK,YQ%MPSLP)
      ENDIF
      IF (YA%LPHY) THEN
        TENDENCY_TMP%A (JL,JK) = TENDENCY_TMP%A (JL,JK) + ZCONS1*PSLPHY9%GFL(JL,JK,YA%MPSLP)
      ENDIF
      IF (YL%LPHY) THEN
        TENDENCY_TMP%CLD (JL,JK,NCLDQL) = TENDENCY_TMP%CLD (JL,JK,NCLDQL)  &
         &  + ZCONS1*PSLPHY9%GFL(JL,JK,YL%MPSLP)
      ENDIF
      IF (YI%LPHY) THEN
        TENDENCY_TMP%CLD (JL,JK,NCLDQI) = TENDENCY_TMP%CLD (JL,JK,NCLDQI)  &
         &  + ZCONS1*PSLPHY9%GFL(JL,JK,YI%MPSLP)
      ENDIF
      IF (YR%LPHY) THEN
        TENDENCY_TMP%CLD (JL,JK,NCLDQR) = TENDENCY_TMP%CLD (JL,JK,NCLDQR)  &
         &  + ZCONS1*PSLPHY9%GFL(JL,JK,YR%MPSLP)
      ENDIF
      IF (YS%LPHY) THEN
        TENDENCY_TMP%CLD (JL,JK,NCLDQS) = TENDENCY_TMP%CLD (JL,JK,NCLDQS)  &
         &  + ZCONS1*PSLPHY9%GFL(JL,JK,YS%MPSLP)
      ENDIF
      IF (YO3%LPHY) THEN
        TENDENCY_TMP%O3 (JL,JK) = TENDENCY_TMP%O3 (JL,JK) + ZCONS1*PSLPHY9%GFL(JL,JK,YO3%MPSLP)
      ENDIF
    ENDDO
  ENDDO
ENDIF

IF ( LEGWDG ) THEN
  !*         3.1     CALL GWD TO EVALUATE TENDENCY COEFFICIENTS
  CALL GWDRAG_LAYER(YDSURF,YDEPHLI,YDEGWD,KDIM,PAUX,STATE_T0,PSURF,AUXL)
ELSE
  CALL NOGWDRAG(KDIM, FLUX, PDIAG, AUXL)
ENDIF

!*         3.2 Evaluate GWD and VDF tendencies jointly
IF ( LEVDIF ) THEN

  IF (LNEMOLIMTHK) THEN
    CALL ICESTATENEMO(YDMCC,KDIM%KSTGLO,KDIM%KIDIA,KDIM%KFDIA,&
      & PTHKICE=SURFL%ZTHKICE(KDIM%KIDIA:KDIM%KFDIA),PSNTICE=SURFL%ZSNTICE(KDIM%KIDIA:KDIM%KFDIA))
  ENDIF

  IF (LPHYLIN) THEN
    ! Simplified scheme
    CALL TURBULENCE_S_LAYER(YDSURF,YDMODEL,KDIM,PAUX,STATE_T0,TENDENCY_TMP,AUXL,GEMSL,PSURF,&
      & PCGPP, PCREC, PAG, PRECO, SURFL,FLUX,PDIAG,TENDENCY_LOC)
  ELSE
    ! Full scheme
    CALL TURBULENCE_LAYER(YDSURF,YDMODEL,KDIM, PAUX, STATE_T0, TENDENCY_TMP, PRAD, PPERT,&
       & AUXL, GEMSL, PSURF, SURFL, FLUX, PDIAG, PERTL, LLKEYS, PDDHS, TENDENCY_LOC)
  ENDIF

  ! Increment the tendency_cml by local tendencies from GWD+VDF
  CALL STATE_INCREMENT(KDIM,TENDENCY_CML, PA_U=TENDENCY_LOC%U, PA_V=TENDENCY_LOC%V, PA_T=TENDENCY_LOC%T,&
    & PA_Q=TENDENCY_LOC%Q, PA_A=TENDENCY_LOC%A, PA_QL=TENDENCY_LOC%CLD(:,:,NCLDQL),&
    & PA_QI=TENDENCY_LOC%CLD(:,:,NCLDQI),PA_TKE=TENDENCY_LOC%TKE)

  ! Add cloud contribution to T and q when there is no prognostic cloud quantity
  IF (.NOT. LEPCLD) THEN
    DO JK=1,KDIM%KLEV
      DO JL=KDIM%KIDIA,KDIM%KFDIA
        TENDENCY_CML%Q(JL,JK) = TENDENCY_CML%Q(JL,JK)   + TENDENCY_LOC%CLD(JL,JK,NCLDQL)&
                                                      & + TENDENCY_LOC%CLD(JL,JK,NCLDQI)
        TENDENCY_CML%T(JL,JK) = TENDENCY_CML%T(JL,JK)   -  (&
             & RLVTT*TENDENCY_LOC%CLD(JL,JK,NCLDQL) + RLSTT*TENDENCY_LOC%CLD(JL,JK,NCLDQI) ) /(&
             & RCPD * ( 1.0_JPRB + RVTMP2 * STATE_T0%Q(JL,JK) ) )
      ENDDO
    ENDDO

    !iSPPT: Store cloud tendencies for stochastic perturbations
    IF (LLISPPT) THEN
      J2D = YDSPPT%MPSDT(4)   !pattern ID for CLouD tendency perturbations
      IF (J2D > 0) THEN     ! .EQ.0 => do NOT perturb tendency
        DO JK=1,KDIM%KLEV
          DO JL=KDIM%KIDIA,KDIM%KFDIA
            TENDENCY_PHY(J2D)%Q(JL,JK) = TENDENCY_PHY(J2D)%Q(JL,JK) + &
             &  TENDENCY_LOC%CLD(JL,JK,NCLDQL) + TENDENCY_LOC%CLD(JL,JK,NCLDQI)
            TENDENCY_PHY(J2D)%T(JL,JK) = TENDENCY_PHY(J2D)%T(JL,JK) - &
             &  ( RLVTT*TENDENCY_LOC%CLD(JL,JK,NCLDQL) +   &
             &    RLSTT*TENDENCY_LOC%CLD(JL,JK,NCLDQI) ) / &
             &  ( RCPD * ( 1.0_JPRB + RVTMP2 * STATE_T0%Q(JL,JK) ) )
          ENDDO
        ENDDO
      ENDIF
    ENDIF ! LLISPPT
  ENDIF

  ! storing contribution of VDIF (used in diagnostics and SLPHY)
  CALL STATE_COPY(KDIM,TENDENCY_LOC,TENDENCY_VDF)

ELSE

!*         3.1     NECESSARY COMPUTATIONS IF SUBROUTINE IS BY-PASSED.
!*                 note, computations if gwd by-passed are handled above
  CALL NOTURBULENCE(YDSURF,KDIM, PDIAG,FLUX,PSURF,SURFL,AUXL,PDDHS,TENDENCY_VDF)

ENDIF

!       Full-Budget-Run: vertical diffusion + gravity wave drag
IF(LEVDIF.OR.LEGWDG) THEN

  IF(LLBUD23) CALL UPDATE_FIELDS(YDPHY2,1,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,&
   & PTA1=TENDENCY_LOC%U, PO1=PSURF%PSD_XA(:,:,6),&
   & PTA2=TENDENCY_LOC%V, PO2=PSURF%PSD_XA(:,:,7),&
   & PTA3=TENDENCY_LOC%T, PO3=PSURF%PSD_XA(:,:,8),&
   & PTA4=TENDENCY_LOC%Q, PO4=PSURF%PSD_XA(:,:,9),&
   & PTA5=AUXL%ZSOTEU, PO5=PSURF%PSD_XA(:,:,10),&
   & PTA6=AUXL%ZSOTEV, PO6=PSURF%PSD_XA(:,:,11),&
   & LDV7=YL%LT1, PTA7=TENDENCY_LOC%CLD(:,:,NCLDQL), PO7=PSURF%PSD_XA(:,:,21),&
   & LDV8=YI%LT1, PTA8=TENDENCY_LOC%CLD(:,:,NCLDQI), PO8=PSURF%PSD_XA(:,:,22))

  !iSPPT: Store VDF/GWD tendencies for stochastic perturbations
  IF (LLISPPT) THEN
    J2D = YDSPPT%MPSDT(2)   !pattern ID for VDF/GWD tendency perturbations
    ! .EQ.0 => do NOT perturb tendency
    IF (J2D > 0) CALL STATE_INCREMENT(KDIM, TENDENCY_PHY(J2D), &
     & PA_U=TENDENCY_LOC%U, PA_V=TENDENCY_LOC%V, PA_T=TENDENCY_LOC%T, PA_Q=TENDENCY_LOC%Q)
  ENDIF

ENDIF

! Stored quantities for 1D var
IF (LERAIN) CALL UPDATE_FIELDS(YDPHY2, 2,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,&
 & PI1=TENDENCY_CML%T,PO1=PDIAG%PZTENT, PI2=TENDENCY_CML%Q,PO2=PDIAG%PZTENQ)

!     ------------------------------------------------------------------

IF ( LLSLPHY ) THEN
  IF (NGEMS > 0 .OR. NCHEM > 0) THEN
    DO ITRC=1,GEMSL%ITRAC
      DO JK=1,KDIM%KLEV
        DO JL=KDIM%KIDIA,KDIM%KFDIA
          GEMSL%ZCEN(JL,JK,ITRC)=GEMSL%ZCEN(JL,JK,ITRC)+GEMSL%ZTENC(JL,JK,ITRC) *TSPHY
        ENDDO
      ENDDO
    ENDDO
  ENDIF
ENDIF

!          4.4  SAVE TENDENCIES IF THE SDAS sEKF ANALYSIS IS PERFORMED

IF (LUSEKF_REF) CALL UPDATE_FIELDS(YDPHY2, 2,KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV,&
  & PI1=TENDENCY_CML%T, PO1=ZT_SAVE_KF,   PI2=TENDENCY_CML%Q, PO2=ZQ_SAVE_KF,&
  & PI3=TENDENCY_CML%U, PO3=ZU_SAVE_KF,   PI4=TENDENCY_CML%V, PO4=ZV_SAVE_KF)

!     ------------------------------------------------------------------

!*         5a.     get the aerosols to pass down to the cloud schemes
!                  ------------------------------------------------------

IF (NAERCLD > 0) THEN
  CALL AER_CLOUD_LAYER(YDMODEL,KDIM,PAUX,STATE_T0,PDIAG,GEMSL,AUXL)
ELSE
  ! routine is by-passed
  AUXL%ZLCRIT_AER(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)=0.0_JPRB
  AUXL%ZICRIT_AER(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)=0.0_JPRB
  AUXL%ZRE_LIQ(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)=0.0_JPRB
  AUXL%ZRE_ICE(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)=0.0_JPRB
  AUXL%ZCCN(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)=0.0_JPRB
  AUXL%ZNICE(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)=0.0_JPRB
ENDIF

!     ------------------------------------------------------------------

!*         5.     CONVECTION PARAMETRIZATION
!                 --------------------------

! Define state_tmp based on latest tendency_cml (containing dyn + VDF+GWD + radiation)
IF (.NOT. LPHYLIN) CALL STATE_UPDATE(YDPHY2,KDIM,STATE_T0,TENDENCY_CML,STATE_TMP)

!       WRITE CUMULUS CONVECTION

IF (LECUMF) THEN

  IF (LCUCONV_CA.AND.LMFCUCA) THEN
    !Perturb T and Q input profile to convection scheme if CA is active
    CALL CONVECTION_CA_LAYER(YDECUCONVCA,YDECUMF,YDPHY2,KDIM,PAUX,STATE_TMP,TENDENCY_CML,PPERT,PERTL)
  ENDIF

  !*    5.1    CALL CONVECTION SCHEME
  IF (.NOT.LECUMFS) THEN
    ! TIEDTKE'S MASS FLUX CONVECTION SCHEME
    CALL CONVECTION_LAYER(YDSURF,YDMODEL,KDIM, LLSLPHY, STATE_T0, TENDENCY_CML, TENDENCY_DYN, PAUX, PPERT,&
      & LLKEYS, PDIAG, AUXL, PERTL, FLUX, PSURF, GEMSL, TENDENCY_LOC)
  ELSE
    ! SIMPLIFIED MASS FLUX CONVECTION SCHEME (P. Lopez)
    CALL CONVECTION_S_LAYER(YDSURF,YDERAD,YDMODEL%YRML_PHY_SLIN,YDMODEL%YRML_PHY_EC,YDPHY2,KDIM, LLSLPHY, LLRAIN1D, STATE_T0,  &
     &                      TENDENCY_CML, PAUX, YDVDF%RVDIFTS,&
     &                      LLKEYS, PDIAG, AUXL, FLUX, PSURF, GEMSL, TENDENCY_LOC)
  ENDIF

  ! CALCULATE TOTAL AND CLOUD-TO-GROUND LIGHTNING FREQUENCIES AND ASSOCIATED IMPACT ON CHEMISTRY (NOx).
  IF (LCHEM_LIGHT .OR. LELIGHT) THEN
    ITRC=1_JPIM
    DO JL=1,NGFL_EXT
      IF ( TRIM(YEXT(JL)%CNAME) == 'EMILI')    ITRC   = JL
    ENDDO
    CALL LIGHTNING_LAYER(YDMODEL,KDIM,PAUX,LLKEYS,STATE_T0,PDIAG,FLUX,YDVARS%EXT(ITRC))
  ENDIF

  ! UPDATE BACKSCATTER RELATED FIELDS
  IF (LSTOPH_SPBS.OR.LSTOPH_CASBS.OR.LVORTCON) THEN
    ! NOTE: 1/ Input state_tmp not yet updated by convection
    !       2/ Output applies directly to tendency_loc which are then INOUT here
    CALL BACKSCATTER_LAYER(YDGEOMETRY%YRDIM,YDMODEL%YRML_PHY_STOCH,YDMODEL%YRML_DYN%YRDYN,YDPHY2,KDIM,PAUX,&
       & STATE_T0,STATE_TMP,PPERT,PERTL,PDIAG,PSURF, TENDENCY_LOC)
  ENDIF

  ! Add convective contribution to wind gusts : gust factor now 0.3 was 0.6 prior to Cy48r1
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    PSURF%PSD_VD(JL,YSD_VD%Y10FGCV%MP)=0.0_JPRB
    IF(PDIAG%PCAPE(JL)>25._JPRB.AND.PDIAG%ITYPE(JL)==1) THEN
      PSURF%PSD_VD(JL,YSD_VD%Y10FGCV%MP) = 0.3_JPRB* &
                      & MAX(0.0_JPRB,(SQRT(STATE_TMP%U(JL,NJKT4)**2+STATE_TMP%V(JL,NJKT4)**2) &
                      & -SQRT(STATE_TMP%U(JL,NJKT3)**2+STATE_TMP%V(JL,NJKT3)**2)))
      PDIAG%PI10FG(JL)=PDIAG%PI10FG(JL)+ PSURF%PSD_VD(JL,YSD_VD%Y10FGCV%MP)
    ENDIF
  ENDDO
ENDIF

!   Recompute CAPE diagnostic every hour for storage on MARS, erase PCAPE
!   computed previously. Diagnostic always active and based on base state only
IF(MOD(NSTEP*TSPHY,RHOUR)==0.0_JPRB.OR.NSTEP==NSTART)  THEN
   CALL CUANCAPE2(YDECUMF, YDEPHLI, KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV,&
      & PAUX%PRSF1, PAUX%PRS1,  STATE_T0%T, STATE_T0%Q,  PDIAG%ZCAPE, PDIAG%PCIN, PDIAG%PPDEPL)
ELSE
  ! In any case PCIN must be initialized - else Ops fail in DDH.
  PDIAG%PCIN(:,:) = 0.0_JPRB
  PDIAG%ZCAPE(:,:) = 0.0_JPRB
  PDIAG%PPDEPL(:) = 0.0_JPRB
ENDIF

IF (LECUMF) THEN

  !  Initialisations/seeding for convective Cellular Automaton and extra-field output
  IF (LCUCONV_CA.AND.LMFCUCA) THEN
    IF(MOD(NSTEP*TSPHY,3600._JPRB)==0.0_JPRB.OR.NSTEP==NSTART)  THEN
      !calculate number of lives for newborn cells if CA is active
      WHERE (PDIAG%PCIN(KDIM%KIDIA:KDIM%KFDIA,1)<100.0_JPRB)
        PPERT%PCAPECONVCA(KDIM%KIDIA:KDIM%KFDIA)=REAL(NLIVES,JPRB)
      ELSEWHERE
        PPERT%PCAPECONVCA(KDIM%KIDIA:KDIM%KFDIA)=0.0_JPRB
      ENDWHERE
    ENDIF
    IF (LEXTRAFIELDS)  PSURF%PSD_X2(KDIM%KIDIA:KDIM%KFDIA,1)=PPERT%PCUCONVCA(KDIM%KIDIA:KDIM%KFDIA)
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      PPERT%PNLCONVCA(JL)=PPERT%PCAPECONVCA(JL)
      IF(PDIAG%ITYPE(JL)==1 .AND. PPERT%PNLCONVCA(JL) > 1 ) THEN
        PPERT%PCUCONVCA(JL)=1 !nlives
      ELSE
        PPERT%PCUCONVCA(JL)=0
      ENDIF
    ENDDO
  ENDIF

  !       Full-Budget-Run: cumulus convection
  IF (LLBUD23) THEN
    CALL UPDATE_FIELDS(YDPHY2,1,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,&
      & PTA1=TENDENCY_LOC%U, PO1=PSURF%PSD_XA(:,:,13),&
      & PTA2=TENDENCY_LOC%V, PO2=PSURF%PSD_XA(:,:,14),&
      & PTA3=TENDENCY_LOC%T, PO3=PSURF%PSD_XA(:,:,15),&
      & PTA4=TENDENCY_LOC%Q, PO4=PSURF%PSD_XA(:,:,16),&
      & PTA5=FLUX%PFPLCL(:,1:KDIM%KLEV), PO5=PSURF%PSD_XA(:,:,17),&
      & PTA6=FLUX%PFPLCN(:,1:KDIM%KLEV), PO6=PSURF%PSD_XA(:,:,18) )
    CALL UPDATE_FIELDS (YDPHY2,2,KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, 1,&
      & PI1=REAL(PDIAG%ICTOP,JPRB), PO1=PSURF%PSD_XA(:,1,25),&
      & PI2=REAL(PDIAG%ICBOT,JPRB), PO2=PSURF%PSD_XA(:,2,25),&
      & PI3=REAL(PDIAG%ITYPE,JPRB), PO3=PSURF%PSD_XA(:,3,25))
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      SELECT CASE (PDIAG%ITYPE(JL))
      CASE (1)
        PSURF%PSD_XA(JL,4,25) = PSURF%PSD_XA(JL,4,25) + 1.0_JPRB
      CASE (2)
        PSURF%PSD_XA(JL,5,25) = PSURF%PSD_XA(JL,5,25) + 1.0_JPRB
      CASE (3)
        PSURF%PSD_XA(JL,6,25) = PSURF%PSD_XA(JL,6,25) + 1.0_JPRB
      END SELECT
    ENDDO
  ENDIF

  !*    5.2    SET CONVECTIVE CLOUD PARAMETERS FOR NEXT TIME-STEP
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    PSURF%PSD_VN(JL,YSD_VN%YTOP%MP)=REAL ( AUXL%ITOPC(JL),JPRB)  + 0.01_JPRB
    PSURF%PSD_VN(JL,YSD_VN%YBAS%MP)=REAL ( AUXL%IBASC(JL),JPRB)  + 0.01_JPRB
  ENDDO

  !     ------------------------------------------------------------------

  ! Increment tendencies after convection
  CALL STATE_INCREMENT(KDIM,TENDENCY_CML, PA_T=TENDENCY_LOC%T, PA_Q=TENDENCY_LOC%Q,&
     & PA_U=TENDENCY_LOC%U, PA_V=TENDENCY_LOC%V)

  IF(LDIAGTURB_EC) ZTENDT_CONV(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)=TENDENCY_LOC%T(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)

  !iSPPT Store convection tendencies for stochastic perturbations
  IF (LLISPPT) THEN
    J2D = YDSPPT%MPSDT(3)   !pattern ID for CONvection tendency perturbations
    ! .EQ.0 => do NOT perturb tendency
    IF (J2D > 0) CALL STATE_INCREMENT(KDIM, TENDENCY_PHY(J2D), &
     & PA_U=TENDENCY_LOC%U, PA_V=TENDENCY_LOC%V, PA_T=TENDENCY_LOC%T, PA_Q=TENDENCY_LOC%Q)
  ENDIF

ELSE

  !*         5.3     NECESSARY COMPUTATIONS IF SUBROUTINE IS BY-PASSED.
  CALL NOCONVECTION(YDEPHY,KDIM,FLUX,LLKEYS,PDIAG)
  IF(LDIAGTURB_EC) ZTENDT_CONV(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)=0.0_JPRB

ENDIF

!     ------------------------------------------------------------------

!*         6.     LARGE SCALE WATER PHASE CHANGES
!                 -------------------------------

!       WRITE LARGE SCALE CONDENSATION INCREMENT

IF ( LECOND .AND. (.NOT. (LENCLD2.OR.LEPCLD2)) ) THEN

  !*         6.1    CALL COND
  CALL COND_LAYER(YDECLDP,YDMODEL%YRML_PHY_EC%YRECND,YDEPHLI,YDPHY2,KDIM,PAUX,STATE_T0,TENDENCY_CML,FLUX,PDIAG, &
 &                TENDENCY_LOC)
  LLCLDDIAG=.false.

ELSEIF( (LENCLD2 .OR. LEPCLD2) .AND. (.NOT. LEPCLD) ) THEN

  !*           6.2 CALL CLOUDST or CLOUDSC2
  CALL CLOUD_S_LAYER(YDMODEL, KDIM, LLRAIN1D, PAUX, STATE_T0, TENDENCY_CML,&
     &  AUXL, FLUX, PDIAG, TENDENCY_LOC)

ELSEIF( LEPCLD ) THEN

  DO JRF=1, YDSPP_CONFIG%SM%NRFTOTAL
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      ZGP2DSPP(JL,JRF)=PPERT%PGP2DSPP(JL,1,JRF)
    ENDDO
  ENDDO

  !*          6.2a CALL cloud scheme saturation adjustment
  CALL CLOUD_SATADJ(YDECLDP, YDEPHLI, YDECUMF, YDEPHY, YDSPP_CONFIG,&
      ! Array dimensions
    & KDIM%KIDIA,    KDIM%KFDIA,    KDIM%KLON,    KDIM%KLEV, PDIAG%ITYPE, &
      ! Timestep and pressure (full and half-level)
    & YDPHY2%TSPHY, PAUX%PRSF1, PAUX%PRS1, &
      ! State at start of timestep
    & STATE_T0%T, STATE_T0%Q, STATE_T0%A, &
    & STATE_T0%CLD(:,:,NCLDQL),STATE_T0%CLD(:,:,NCLDQI), &
      ! Tendencies so far in timestep
    & TENDENCY_CML%T, TENDENCY_CML%Q, TENDENCY_CML%A, &
    & TENDENCY_CML%CLD(:,:,NCLDQL), TENDENCY_CML%CLD(:,:,NCLDQI), &
      ! Radiation tendencies (SW, LW heating)
    & PRAD%PHRSW, PRAD%PHRLW, &
      ! Vertical diffusion tendencies
    & TENDENCY_VDF%T, TENDENCY_VDF%Q, &
      ! Convection tendencies (mass flux)
    & PDIAG%PMFU, PDIAG%PMFD, PDIAG%ZLUDELI, &
      ! Dynamics tendencies (vertical velocity)
    & PAUX%PVERVEL, ZGP2DSPP, &
      ! Output tendencies
    & TENDENCY_LOC%T, TENDENCY_LOC%Q, TENDENCY_LOC%A, &
    & TENDENCY_LOC%CLD(:,:,NCLDQL), TENDENCY_LOC%CLD(:,:,NCLDQI), &
    & PSURF%PSD_XA, KDIM%KFLDX)

  ! Increment cumulated tendencies after cloud saturation adjustment
  CALL STATE_INCREMENT(KDIM,TENDENCY_CML,PA_T=TENDENCY_LOC%T,PA_Q=TENDENCY_LOC%Q,PA_A=TENDENCY_LOC%A,&
  & PA_QL=TENDENCY_LOC%CLD(:,:,NCLDQL), PA_QI=TENDENCY_LOC%CLD(:,:,NCLDQI) )

  ! Add saturation adjustment tendencies to cloud scheme (LBUD23)
  IF (LBUD23) CALL UPDATE_FIELDS(YDPHY2,1,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,&
      & PTA1=TENDENCY_LOC%T, PO1=PSURF%PSD_XA(:,:,19),&
      & PTA2=TENDENCY_LOC%Q, PO2=PSURF%PSD_XA(:,:,20),&
      & LDV3=YL%LT1, PTA3=TENDENCY_LOC%CLD(:,:,NCLDQL), PO3=PSURF%PSD_XA(:,:,21),&
      & LDV4=YI%LT1, PTA4=TENDENCY_LOC%CLD(:,:,NCLDQI), PO4=PSURF%PSD_XA(:,:,22))

  ! Store tendency from saturation adjustment (used in SLTEND if LSLPHY=T)
  CALL STATE_COPY(KDIM,TENDENCY_LOC,TENDENCY_SATADJ)

  !*           6.2 CALL CLOUDSC
  CALL CLOUD_LAYER(YDSURF,YDECLDP,YDECUMF,YDEPHLI,YDPHY2,YDERAD,YDEPHY,YDVDF,YDSPP_CONFIG,YGFL, &
   &               KDIM,LLSLPHY, PAUX, PPERT, STATE_T0, &
   &               TENDENCY_CML, TENDENCY_DYN, TENDENCY_VDF, &
   &               PRAD, PSURF, LLKEYS, AUXL, FLUX, PDIAG, &
   &               YDVARS%FSD, TENDENCY_LOC)

ELSE

  !*         6.3     NECESSARY COMPUTATIONS IF SUBROUTINE IS BY-PASSED.
  CALL NOCLOUD(YDECLDP,KDIM,FLUX,PDIAG,TENDENCY_LOC,TENDENCY_SATADJ)
  LLCLDDIAG=.FALSE.

ENDIF

! Increment cumulated tendencies after cloud scheme
CALL STATE_INCREMENT(KDIM,TENDENCY_CML,PA_T=TENDENCY_LOC%T,PA_Q=TENDENCY_LOC%Q,PA_A=TENDENCY_LOC%A,&
  & PA_QL=TENDENCY_LOC%CLD(:,:,NCLDQL), PA_QI=TENDENCY_LOC%CLD(:,:,NCLDQI),&
  & PA_QR=TENDENCY_LOC%CLD(:,:,NCLDQR), PA_QS=TENDENCY_LOC%CLD(:,:,NCLDQS) )

!iSPPT: Store cloud tendencies for stochastic perturbations
IF (LLISPPT) THEN
  J2D = YDSPPT%MPSDT(4)   !pattern ID for CLouD tendency perturbations
  ! .EQ.0 => do NOT perturb tendency
  IF (J2D > 0) CALL STATE_INCREMENT(KDIM, TENDENCY_PHY(J2D),&
   & PA_T=TENDENCY_LOC%T, PA_Q=TENDENCY_LOC%Q)
ENDIF

!       Full-Budget-Run: cloud

IF (LLBUD23) CALL UPDATE_FIELDS(YDPHY2,1,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,&
              & PTA1=TENDENCY_LOC%T,               PO1=PSURF%PSD_XA(:,:,19),&
              & PTA2=TENDENCY_LOC%Q,               PO2=PSURF%PSD_XA(:,:,20),&
              & LDV3=YL%LT1, PTA3=TENDENCY_LOC%CLD(:,:,NCLDQL), PO3=PSURF%PSD_XA(:,:,21),&
              & LDV4=YI%LT1, PTA4=TENDENCY_LOC%CLD(:,:,NCLDQI), PO4=PSURF%PSD_XA(:,:,22),&
              & PTA5=FLUX%PFPLSL(:,1:KDIM%KLEV)  , PO5=PSURF%PSD_XA(:,:,23),&
              & PTA6=FLUX%PFPLSN(:,1:KDIM%KLEV)  , PO6=PSURF%PSD_XA(:,:,24))

!     ------------------------------------------------------------------

!          7.0 WARNER-McINTYRE-SCINOCCA NON-OROGRAPHIC GRAVITY WAVE SCHEME

!  Note call could possibly already be placed after diffusion if coupling with
!  with precip or CAPE not necessary, but calling at end of processes preferable

IF ( LEGWWMS ) THEN

  ! Computation acctivated upon an intermittency frequency
  IF(MOD(NSTEP*TSPHY,GTPHYGWWMS)==0.0_JPRB.OR.NSTEP==NSTART) THEN

    ! Update state
    IF (.NOT. LPHYLIN) CALL STATE_UPDATE(YDPHY2,KDIM,STATE_T0,TENDENCY_CML,STATE_TMP)

    ! call the  Warner-McIntyrE-Scinocca non-orographic gravity wave scheme
    IF (LPHYLIN) THEN
      CALL GWDRAGWMS_S_LAYER(YDEGWD,YDEGWWMS,YGFL,KDIM,STATE_TMP,PAUX,FLUX,YDVARS%NOGW)
    ELSE
      CALL GWDRAGWMS_LAYER(YDSTA,YDEGWD,YDEGWWMS,YGFL,KDIM,STATE_TMP,PAUX,FLUX,YDVARS%NOGW)
    ENDIF

  ENDIF

  ! Compute the process tendency, either from stored or from freshly updated quantities
  DO JK=1,KDIM%KLEV
     DO JL=KDIM%KIDIA,KDIM%KFDIA
        TENDENCY_LOC%T(JL,JK)=-( STATE_TMP%U(JL,JK)*YDVARS%NOGW(1)%P(JL,JK)&
                              & +STATE_TMP%V(JL,JK)*YDVARS%NOGW(2)%P(JL,JK) )*ZRCPD
        TENDENCY_LOC%U(JL,JK)=YDVARS%NOGW(1)%P(JL,JK)
        TENDENCY_LOC%V(JL,JK)=YDVARS%NOGW(2)%P(JL,JK)
     ENDDO
  ENDDO

  ! Increment tendency_cml by non orographic GW
  CALL STATE_INCREMENT(KDIM,TENDENCY_CML, PA_T=TENDENCY_LOC%T, PA_U=TENDENCY_LOC%U, PA_V=TENDENCY_LOC%V )

  ! add to Diff and oro GWD
  IF (LLBUD23) CALL UPDATE_FIELDS(YDPHY2, 1, KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV,&
   & PTA1=TENDENCY_LOC%U,  PO1=PSURF%PSD_XA(:,:,10),&
   & PTA2=TENDENCY_LOC%V,  PO2=PSURF%PSD_XA(:,:,11),&
   & PTA3=TENDENCY_LOC%T,  PO3=PSURF%PSD_XA(:,:,12))
  IF (LDIAG_STRATO) CALL UPDATE_FIELDS(YDPHY2, 1, KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV,&
   & PTA1=FLUX%PSTRDU(:,1:KDIM%KLEV),  PO1=PSURF%PSD_XA(:,:,1),&
   & PTA2=FLUX%PSTRDV(:,1:KDIM%KLEV),  PO2=PSURF%PSD_XA(:,:,2))

  !iSPPT: Store non-orographic GWD tendencies for stochastic perturbations
  IF (LLISPPT) THEN
    J2D = YDSPPT%MPSDT(5)   !pattern ID for Non-Orographic GWdrag tendency perturbations
    ! .EQ.0 => do NOT perturb tendency
    IF (J2D > 0) CALL STATE_INCREMENT(KDIM, TENDENCY_PHY(J2D), &
      & PA_U=TENDENCY_LOC%U, PA_V=TENDENCY_LOC%V, PA_T=TENDENCY_LOC%T)
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!      Cloud Diagnostics (Base, top, zero degree level, CIndices, Precip types)

IF (LLCLDDIAG) CALL DIAG_CLOUDS (&
   & YDECUMF, KDIM%KIDIA,  KDIM%KFDIA,  KDIM%KLON,  KDIM%KLEV, LL_DIAGCLOUDAER,&
   & LLKEYS%LLCUM,  PDIAG%ICBOT,  PDIAG%ICTOP,   PAUX%PRSF1,  PAUX%PRS1,&
   & PAUX%PGEOM1, PAUX%PGEOMH,&
   & STATE_T0%T,    STATE_T0%Q,   STATE_T0%CLD(:,:,NCLDQL),&
   & STATE_T0%CLD(:,:,NCLDQI),    STATE_T0%A,&
   & FLUX%PFPLCL,   FLUX%PFPLCN,  FLUX%PFPLSL,  FLUX%PFPLSN,&
   & AUXL%ZRAINFRAC_TOPRFZ,&
   & PSURF%PSD_VD(:,YSD_VD%Y2T%MP),PSURF%PSD_VD(:,YSD_VD%Y2D%MP),&
   & PDIAG%PCBASE,  PDIAG%PCBASEA,PDIAG%PCTOPC, PDIAG%P0DEGL, PDIAG%PM10DEGL, PDIAG%PCONVIND,&
   & PDIAG%PPRECTYPE, PDIAG%PFZRA, PDIAG%PZTWETB, PDIAG%PTROPOTP )

! 	   Turbulence diagnostics (CAT and Mountain Wave Turb)

IF(LDIAGTURB_EC) THEN
 !IF(MOD(NSTEP*TSPHY,RHOUR)==0.0_JPRB.OR.NSTEP==NSTART) THEN

    CALL DIAG_TURB (&
      & YDEPHY, YDECUMF, KDIM%KIDIA,  KDIM%KFDIA,  KDIM%KLON,  KDIM%KLEV, &
      & LEGWWMS, YDEGWWMS%NLAUNCH, PDIAG%ITYPE, PDIAG%ICTOP, &
      & PAUX%POROG,  PSURF%PHSTD, PAUX%PGAW,&
      & PAUX%PRSF1,  PAUX%PRS1,   PAUX%PGEOM1, PAUX%PGEOMH,&
      & STATE_T0%T,  STATE_T0%Q,  STATE_T0%U,    STATE_T0%V,  PAUX%PVERVEL,&
      & YDVARS%T%DL, YDVARS%T%DM, YDVARS%U%DL, YDVARS%V%DL,  YDVARS%VOR%T0,YDVARS%DIV%T0, PERTL%ZDISSCU,&
      & TENDENCY_VDF%U, TENDENCY_VDF%V, ZTENDT_CONV, PDIAG%ZRI,   YDVARS%NOGW(1)%P, YDVARS%NOGW(2)%P, ZEDRP )

     DO JEXT=1,NEDRP
      DO JK=1,KDIM%KLEV
        DO JL=KDIM%KIDIA,KDIM%KFDIA
          YDVARS%EDRP(JEXT)%P(JL,JK)=ZEDRP(JL,JK,JEXT)
        ENDDO
      ENDDO
     ENDDO

 !ENDIF
ENDIF

!      Diagnostics for diurnal cycle

IF(LBUDCYCLE) THEN
  CALL DIAG_DCYCLE (&
     & KDIM%KIDIA,  KDIM%KFDIA,  KDIM%KLON,   KDIM%KLEV,&
     & KDIM%KLEVX,  KDIM%KFLDX,  TSPHY,&
     & PAUX%PRS1,   FLUX%PFPLCL, FLUX%PFPLCN, FLUX%PFPLSL, FLUX%PFPLSN,&
     & FLUX%PDIFTS, FLUX%PDIFTQ, FLUX%PFRTH,  FLUX%PFRTHC,&
     & STATE_T0%CLD(:,:,NCLDQL),     STATE_T0%CLD(:,:,NCLDQI),     PSURF%PSP_RR(:,YSP_RR%YT%MP9),&
     & PSURF%PSD_VD(:,YSD_VD%Y2T%MP),  PSURF%PSD_VD(:,YSD_VD%Y2SH%MP),&
     & PSURF%PSD_XA )
ENDIF

!     ------------------------------------------------------------------

IF (LUSEKF_REF) THEN
! compute and save tendency increments from the unperturbed first sEKF model run
  IF (N_SEKF_PT==0) THEN
    IBLK=(KDIM%KSTGLO-1)/KDIM%KLON + 1
    CALL UPDATE_FIELDS(YDPHY2, 2,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,&
      & PI1=TENDENCY_CML%T, PS1=ZT_SAVE_KF, PO1=FKF_TENT(:,:,IBLK,NSTEP),&
      & PI2=TENDENCY_CML%Q, PS2=ZQ_SAVE_KF, PO2=FKF_TENQ(:,:,IBLK,NSTEP),&
      & PI3=TENDENCY_CML%U, PS3=ZU_SAVE_KF, PO3=FKF_TENU(:,:,IBLK,NSTEP),&
      & PI4=TENDENCY_CML%V, PS4=ZV_SAVE_KF, PO4=FKF_TENV(:,:,IBLK,NSTEP))
  ENDIF
! replace tendency increments with the unperturbed first sEKF increments
  IF (N_SEKF_PT>0) THEN
    IBLK=(KDIM%KSTGLO-1)/KDIM%KLON + 1
    CALL UPDATE_FIELDS(YDPHY2,2,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,&
      & PI1=FKF_TENT(:,:,IBLK,NSTEP), PA1=ZT_SAVE_KF, PO1=TENDENCY_CML%T,&
      & PI2=FKF_TENQ(:,:,IBLK,NSTEP), PA2=ZQ_SAVE_KF, PO2=TENDENCY_CML%Q,&
      & PI3=FKF_TENU(:,:,IBLK,NSTEP), PA3=ZU_SAVE_KF, PO3=TENDENCY_CML%U,&
      & PI4=FKF_TENV(:,:,IBLK,NSTEP), PA4=ZV_SAVE_KF, PO4=TENDENCY_CML%V)
  ENDIF
ENDIF
!     ------------------------------------------------------------------

!*         6.5     DUCTING DIAGNOSTICS
!                 -------------------

IF (LDUCTDIA) THEN
  CALL DUCTDIA_LAYER(YDSURF,KDIM, PAUX, STATE_T0, PSURF)
ENDIF

!*         6.6     100M WIND DIAGNOSTICS
!                 ----------------------

CALL VDFVINT(KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,&
 & STATE_T0%U  ,STATE_T0%V  ,PAUX%PGEOM1, 100._JPRB,&
 ! OUTPUTS
 & PSURF%PSD_VD(:,YSD_VD%Y100U%MP), PSURF%PSD_VD(:,YSD_VD%Y100V%MP) )


!*         6.7     200M WIND DIAGNOSTICS
!                 ----------------------

CALL VDFVINT(KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,&
 & STATE_T0%U  ,STATE_T0%V  ,PAUX%PGEOM1, 200._JPRB,&
 ! OUTPUTS
 & PSURF%PSD_VD(:,YSD_VD%Y200U%MP), PSURF%PSD_VD(:,YSD_VD%Y200V%MP) )


!     ------------------------------------------------------------------

!*         8.     GLOMAP AEROSOLS DIAGNOSTICS
!                 --------------------------------

IF (NACTAERO /=0 ) THEN
  IF (.NOT.LAERINIT .AND. TRIM(AERO_SCHEME) == "glomap")  THEN
    CALL AER_GLOMAPDIAG_LAYER(KDIM, TSPHY, PAUX, STATE_TMP, PDIAG, GEMSL)
  ENDIF
ENDIF



!     ------------------------------------------------------------------

!*         9.     METHANE OXIDATION AND PHOTOLYSIS AND MORE CHEMISTRY
!                 --------------------------------

! Update od prognostic vatiables used in the subsequent computation
! The actual state is consistent with tendencies, i.e. doesn't contain the information from the
!   first cloud pass...
IF (.NOT. LPHYLIN) CALL STATE_UPDATE(YDPHY2,KDIM,STATE_T0,TENDENCY_CML,STATE_TMP)

IF(LEMETHOX) THEN

  TENDENCY_LOC%Q(:,:)= 0._JPRB

  CALL METHOX(KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV, STATE_TMP%Q, TENDENCY_LOC%Q, PAUX%PRSF1 )

  ! Increment tendency
  CALL STATE_INCREMENT(KDIM,TENDENCY_CML, PA_Q=TENDENCY_LOC%Q)
  ! Update state
  IF (.NOT. LPHYLIN) CALL STATE_UPDATE(YDPHY2,KDIM,STATE_T0,TENDENCY_CML,STATE_TMP)

  !iSPPT Store METHOX tendencies for stochastic perturbations
  IF (LLISPPT) THEN
    J2D = YDSPPT%MPSDT(6)   !pattern ID for MethOX tendency perturbations
    ! .EQ.0 => do NOT perturb tendency
    IF (J2D > 0) CALL STATE_INCREMENT(KDIM, TENDENCY_PHY(J2D), PA_Q=TENDENCY_LOC%Q)
  ENDIF

  !       Full-Budget-Run:
  !       ADD HUMIDITY SOURCE FROM METHAN OXIDATION
  IF (LLBUD23) CALL UPDATE_FIELDS(YDPHY2, 1,KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV, &
   & PTA1=TENDENCY_LOC%Q, PO1=PSURF%PSD_XA(:,:,20))

ENDIF

! since this point state_tmp is no longer updated  (if required, do the update at the
!    appropriate position)

!     ------------------------------------------------------------------

! Preparation for chemistry

IF ((NCHEM > 0) .OR. LEO3CH ) THEN
  CALL GPINO3CH(YDGEOMETRY%YRDIMV,YDMODEL%YRML_CHEM%YROZO,YDDPHY,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KGPLAT,GEMSL%ZKOZO)
ENDIF

!     ------------------------------------------------------------------

!*         9.1   Full Chemistry
!                 --------------------------------

! This must run BEFORE prognostic aerosol sources (AER_PHY3) if using LAERCHEM coupling, so that
! the wet oxidation from the chemistry scheme has been calculated.

! also call for GHG or aerosol only for diagnostics
IF ((NCHEM > 0 .OR. NGEMS > 0) .AND. LL_HRES ) THEN

  ! This is not very consistent, having updated state but pressure only appropriate to t0
  CALL CHEM_MAIN_LAYER(YDVAB,YDGEOMETRY%YRDIMV,YDSURF,YDMODEL,KDIM,PAUX,STATE_TMP,PDIAG,SURFL,FLUX,PGFL,PTENGFL, &
   & GEMSL,PSURF,ZCHEM2AER)

  ! This would be then much more consistent alternative...
  !CALL CHEM_MAIN_LAYER(KDIM, PAUX, state_t0, PDIAG, SURFL, FLUX, PGFL, PTENGFL, GEMSL, PSURF)

ENDIF

!     ------------------------------------------------------------------

!*         8.     PROGNOSTIC AEROSOL SOURCES - END
!                 --------------------------------

DO JL=KDIM%KIDIA,KDIM%KFDIA
  GEMSL%ZPRAERS(JL)=0._JPRB
ENDDO

IF (NACTAERO /=0 .AND. LL_HRES) THEN

  IF (.NOT.LAERINIT) THEN
    CALL AER_PHY3_LAYER(YDSURF,YDMODEL,KDIM, PAUX, STATE_TMP, SURFL, AUXL, PDIAG, ZCHEM2AER, PGFL, PSURF, FLUX, GEMSL, PRAD)
  ELSE
    GEMSL%ZTENC(KDIM%KIDIA:KDIM%KFDIA, 1:KDIM%KLEV, GEMSL%IAERO(1):GEMSL%IAERO(NACTAERO))=0._JPRB
  ENDIF

ENDIF

!*         8.1    DIAGNOSTIC OF VISIBILITY
!                 ------------------------

IF (LAERVISI) THEN
  CALL RADVIS_LAYER(YDECLDP,YDERAD,YGFL,YDPHY2,KDIM,PAUX,STATE_TMP,GEMSL,PDIAG)
ENDIF

!     ------------------------------------------------------------------

!*        10.0   Use preciptation from the refernce run (for the sEKF)

IF (LUSEKF_REF) THEN
! compute and save tendency increments from the unperturbed first sEKF model run
  IF (N_SEKF_PT==0) THEN
    IBLK=(KDIM%KSTGLO-1)/KDIM%KLON + 1
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      FKF_SURF_CR(JL,IBLK,NSTEP) = FLUX%PFPLCL(JL,KDIM%KLEV)
      FKF_SURF_LR(JL,IBLK,NSTEP) = FLUX%PFPLSL(JL,KDIM%KLEV)
    ENDDO
  ENDIF
  IF (N_SEKF_PT>0) THEN
    IBLK=(KDIM%KSTGLO-1)/KDIM%KLON + 1
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      FLUX%PFPLCL(JL,KDIM%KLEV) = FKF_SURF_CR(JL,IBLK,NSTEP)
      FLUX%PFPLSL(JL,KDIM%KLEV) = FKF_SURF_LR(JL,IBLK,NSTEP)
    ENDDO
  ENDIF
ENDIF


!*         11.      COMPUTATION OF NEW SURFACE VALUES

IF (LESURF) THEN
  IF (LPHYSFCLIN) THEN
    ! Simplified scheme
    CALL SURFTSTP_S_LAYER(YDSURF,YDEPHY,YDPHY2,PAUX,KDIM,PSURF,SURFL,LLKEYS,FLUX)
    ! Unused fluxes and DDH have to be secured in this case
    CALL NOSURFTSTP(KDIM,FLUX,PDDHS)
  ELSE
    !Coupled surface - full scheme
    CALL SURFTSTP_LAYER(YDSURF,YDEPHY,YDRIP,YDPHY2,KDIM,STATE_T0,PAUX,PSURF,SURFL,LLKEYS,FLUX,PDDHS)
  ENDIF
ELSE
  ! Just securing some fluxes and DDH
  CALL NOSURFTSTP(KDIM,FLUX,PDDHS)
ENDIF

!     ------------------------------------------------------------------

!     12a.  ADD STOCHASTIC PERTURBATION TO PHYSICS TENDENCY IF REQUIRED

!          ----------------------------------------------------------

IF (YDSPPT_CONFIG%LSPSDT.AND.NSTEP >= 0) THEN
  ! This is an exception from the norm: every other process provides only
  !  output tendencies for its own, here the existing tendencies are modified

  IF (KDIM%K2DSDT==1) THEN
   J2D=1
   !Prepare total physics tendencies for SPPT: net minus dynamics
   CALL UPDATE_FIELDS(YDPHY2, 2,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV, &
     & PI1=TENDENCY_CML%T, PS1=TENDENCY_DYN%T, PO1=TENDENCY_PHY(J2D)%T, &
     & PI2=TENDENCY_CML%Q, PS2=TENDENCY_DYN%Q, PO2=TENDENCY_PHY(J2D)%Q, &
     & PI3=TENDENCY_CML%U, PS3=TENDENCY_DYN%U, PO3=TENDENCY_PHY(J2D)%U, &
     & PI4=TENDENCY_CML%V, PS4=TENDENCY_DYN%V, PO4=TENDENCY_PHY(J2D)%V)

    !Default (from 45r1): SPPT does *not* perturb clear-skies radiation tendencies
    IF (.NOT.(YDSPPT_CONFIG%LRADCLR_SDT)) THEN
      !Remove clear-skies heating rates
      ! SW contribution:
      CALL UPDATE_FIELDS(YDPHY2,2,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV, &
        & PS1=PRAD%PHRSC, PO1=TENDENCY_PHY(J2D)%T)
      ! LW contribution:
      CALL UPDATE_FIELDS(YDPHY2,2,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV, &
        & PS1=PRAD%PHRLC, PO1=TENDENCY_PHY(J2D)%T)
    ENDIF
    IF (.NOT.(YDSPPT_CONFIG%LSATADJ_SDT)) THEN
      !Remove saturation adjustment
      CALL UPDATE_FIELDS(YDPHY2,2,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV, &
           & PS1=TENDENCY_SATADJ%T, PO1=TENDENCY_PHY(J2D)%T, &
           & PS2=TENDENCY_SATADJ%Q, PO2=TENDENCY_PHY(J2D)%Q)
    ENDIF
  ENDIF

  CALL STOCHPERT_LAYER(YDMODEL,YGFL,YDPHY2, KDIM, PAUX, STATE_TMP, TENDENCY_DYN, TENDENCY_PHY, &
    & PTENGFL, PPERT, TENDENCY_CML, GEMSL)
ENDIF

!     ------------------------------------------------------------------

!*         13.     TRACER TRANSPORT TENDENCIES BACK INTO EXTRA FIELDS
!                 --------------------------------------------------

IF (NGEMS > 0 .OR. NCHEM > 0) THEN
  CALL GEMS_TEND(YGFL,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLEV,KDIM%KLON,PTENGFL,GEMSL%ZTENC)
ENDIF

!     ------------------------------------------------------------------

!*         14.      OZONE CHEMISTRY
!                   ---------------

IF (LEO3CH .AND. YO3%LGP) THEN
  CALL O3CHEM(YDRIP,YDEPHY,YDOZO,YDPHY2, &
   & KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, 1, KDIM%KLEV, KDIM%KVCLIS,PAUX%PGEMU, PAUX%PMU0,&
   & PAUX%PRS1, PAUX%PRSF1, GEMSL%ZKOZO, PAUX%PDELP, STATE_TMP%T, STATE_TMP%O3, TENDENCY_LOC%O3)

  ! Increment tendency
  CALL STATE_INCREMENT(KDIM,TENDENCY_CML, PA_O3=TENDENCY_LOC%O3)

ENDIF

!     ------------------------------------------------------------------

!*         15. Add non-linear stochastic forcing terms
!*
!*

IF(LFORCENL.AND.(NSTEP*(TSPHY/RHOUR)>=NFORCESTART).AND.(&
   & NSTEP*(TSPHY/RHOUR)<=NFORCEEND)) THEN

  CALL STATE_INCREMENT(KDIM,TENDENCY_CML,PA_U=PPERT%PFORCEU,PA_V=PPERT%PFORCEV,&
                     & PA_T=PPERT%PFORCET,PA_Q=PPERT%PFORCEQ)

ENDIF

!     ------------------------------------------------------------------

!*              XX.    COMPUTE AND ACCUMULATE FLUXES FOR COUPLING TO OPA/LIM
!                      -----------------------------------------------------

IF (LNEMOATMFLDS) CALL NEMOADDFLDS_LAYER(YDSURF,YDMCC,KDIM, SURFL, PSURF)
IF (LNEMOLIMPUT.OR.CPL_NEMO_LIM) CALL SET_OCEAN_FLUXES(YDSURF,YDMCC,KDIM,SURFL,PSURF,FLUX)

!     ------------------------------------------------------------------

!*         16.   SLPHY,    build final profiles,
!*                         create tendencies to save for the next time-step
!*                         remove supersaturation if anything there

IF( LLSLPHY ) THEN

  CALL SLTEND_LAYER(YDMODEL, KDIM, PAUX, STATE_T0, &
     & TENDENCY_DYN, TENDENCY_VDF, TENDENCY_SATADJ, TENDENCY_CML,&
     & PSLPHY9, PSAVTEND, PGFLSLP, PSURF, TENDENCY_LOC)

  CALL STATE_INCREMENT(KDIM,TENDENCY_CML, PA_U=TENDENCY_LOC%U, PA_V=TENDENCY_LOC%V,&
   & PA_T=TENDENCY_LOC%T, PA_Q=TENDENCY_LOC%Q, PA_A=TENDENCY_LOC%A, PA_O3=TENDENCY_LOC%O3,&
   & PA_QL=TENDENCY_LOC%CLD(:,:,NCLDQL), PA_QI=TENDENCY_LOC%CLD(:,:,NCLDQI), &
   & PA_QR=TENDENCY_LOC%CLD(:,:,NCLDQR), PA_QS=TENDENCY_LOC%CLD(:,:,NCLDQS) )
  
  !       Full-Budget-Run:
  !       ADD SATURATION ADJUSTMENT INCREMENTS TO CLOUD BUDGETS
  IF (LLBUD23) CALL UPDATE_FIELDS(YDPHY2,1,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,&
   & PTA1=TENDENCY_LOC%T, PO1=PSURF%PSD_XA(:,:,19),&
   & PTA2=TENDENCY_LOC%Q, PO2=PSURF%PSD_XA(:,:,20) )

ENDIF

!*        17.     ELIMINATION OF NEGATIVE SPECIFIC HUMIDITIES
!                 -------------------------------------------

IF( LEQNGT ) THEN

  !*         17.1   CALL QNEGAT
  CALL QNEGAT (&
   & YDMODEL%YRML_PHY_EC%YRECND, KDIM%KIDIA , KDIM%KFDIA , KDIM%KLON , KDIM%KLEV,&
   & TSPHY, STATE_T0%Q , TENDENCY_LOC%Q, PAUX%PRS1,&
   ! FLUX OUTPUTS
   & FLUX%PFCQNG, PITEND=TENDENCY_CML%Q )

  IF (LLBUD23) CALL UPDATE_FIELDS(YDPHY2, 1,KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV,&
   & PTA1=TENDENCY_LOC%Q, PO1=PSURF%PSD_XA(:,:,20))

  CALL STATE_INCREMENT(KDIM,TENDENCY_CML, PA_Q=TENDENCY_LOC%Q)

ENDIF

! Zeroing tendencies for non prognostic & not advected quantities 
IF (.NOT.LEPCLD .AND. (LENCLD2.OR.LEPCLD2)) THEN
  IF (YR%LACTIVE) PTENGFL(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,YR%MP1)=0.0_JPRB
  IF (YS%LACTIVE) PTENGFL(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,YS%MP1)=0.0_JPRB
ENDIF

!     ------------------------------------------------------------------

!*              18.     END OF E.C.M.W.F. PHYSICS
!                      -------------------------

!*              18.1    UNIT CHANGES/SCALINGS
DO JL=KDIM%KIDIA,KDIM%KFDIA
  PSURF%PSD_VD(JL,YSD_VD%YZ0F%MP) =GEMSL%ZAZ0M(JL)*RG
  PSURF%PSD_VD(JL,YSD_VD%YLZ0H%MP)=LOG(GEMSL%ZAZ0H(JL))
ENDDO


!*      19. DIAGNOSTIC FIELDS USED BY OBSERVATION OPERATOR FOR MW-RADIANCE ASSIMILATION
IF(YDEPHY%LEMWAVE) THEN
  DO JK=1,KDIM%KLEV
   DO JL = KDIM%KIDIA, KDIM%KFDIA
      PHYS_MWAVE(JL,JK,1)=STATE_T0%CLD(JL,JK,NCLDQL)
      PHYS_MWAVE(JL,JK,2)=STATE_T0%CLD(JL,JK,NCLDQI)
      PHYS_MWAVE(JL,JK,3)=STATE_T0%A (JL,JK)
      PHYS_MWAVE(JL,JK,4)=0.5_JPRB * (FLUX%PFPLSL (JL,JK) + FLUX%PFPLSL (JL,JK-1))
      PHYS_MWAVE(JL,JK,5)=0.5_JPRB * (FLUX%PFPLSN (JL,JK) + FLUX%PFPLSN (JL,JK-1))
      PHYS_MWAVE(JL,JK,6)=0.5_JPRB * (FLUX%PFPLCL (JL,JK) + FLUX%PFPLCL (JL,JK-1))
      PHYS_MWAVE(JL,JK,7)=0.5_JPRB * (FLUX%PFPLCN (JL,JK) + FLUX%PFPLCN (JL,JK-1))
      PHYS_MWAVE(JL,JK,8)=PDIAG%PCOVPTOT(JL,JK)
    ENDDO
  ENDDO
  DO JL = KDIM%KIDIA, KDIM%KFDIA
     PHYS_MWAVE(JL,1,9)=PDIAG%IPBLTYPE(JL)
     PHYS_MWAVE(JL,2,9)=PDIAG%ITYPE(JL)
     PHYS_MWAVE(JL,3:KDIM%KLEV,9)=0.0_JPRB
  ENDDO
ENDIF

! Write precip fraction into array to pass up to DDH
DO JK=1,KDIM%KLEV
   DO JL = KDIM%KIDIA, KDIM%KFDIA
      PSURF%PCOVPTOT(JL,JK)=PDIAG%PCOVPTOT(JL,JK)
   ENDDO
ENDDO

! ---------------------------------------------------------------------------------------
 
! Output total physics tendencies (after stochastic physics)
IF (LEXTRATEND) THEN
  CALL UPDATE_FIELDS(YDPHY2,1,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,  &
    & PTA1=TENDENCY_CML%U, PTS1=TENDENCY_DYN%U, PO1=PSURF%PSD_XA(:,:,26),               &
    & PTA2=TENDENCY_CML%V, PTS2=TENDENCY_DYN%V, PO2=PSURF%PSD_XA(:,:,27),               &
    & PTA3=TENDENCY_CML%T, PTS3=TENDENCY_DYN%T, PO3=PSURF%PSD_XA(:,:,28),               &
    & PTA4=TENDENCY_CML%Q, PTS4=TENDENCY_DYN%Q, PO4=PSURF%PSD_XA(:,:,29))
ENDIF


! ---------------------------------------------------------------------------------------

!     -----------------------------------------------------
!      Copy the updated tendencies back into global arrays
!     -----------------------------------------------------

! Note: this couln't be just tendency_cml to tendency_dyn copy.  
!  The cloud part of tendency_dyn is not a pointer to the appropriate GFL structure.
CALL UPDATE_FIELDS(YDPHY2, 2,KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV,&
 & PI1=TENDENCY_CML%U,               PO1=TENDENCY_DYN%U,&
 & PI2=TENDENCY_CML%V,               PO2=TENDENCY_DYN%V,&
 & PI3=TENDENCY_CML%T,               PO3=TENDENCY_DYN%T,&
 & PI4=TENDENCY_CML%Q,               PO4=TENDENCY_DYN%Q,&
 & LDV5=YA%LACTIVE,  PI5=TENDENCY_CML%A,               PO5=TENDENCY_DYN%A,&
 & LDV6=YO3%LACTIVE, PI6=TENDENCY_CML%O3,              PO6=TENDENCY_DYN%O3,&
 & LDV7=YL%LACTIVE,  PI7=TENDENCY_CML%CLD(:,:,NCLDQL), PO7=PTENGFL(:,:,YL%MP1),&
 & LDV8=YI%LACTIVE,  PI8=TENDENCY_CML%CLD(:,:,NCLDQI), PO8=PTENGFL(:,:,YI%MP1),&
 & LDV9=YR%LACTIVE,  PI9=TENDENCY_CML%CLD(:,:,NCLDQR), PO9=PTENGFL(:,:,YR%MP1),&
 & LDV10=YS%LACTIVE, PI10=TENDENCY_CML%CLD(:,:,NCLDQS),PO10=PTENGFL(:,:,YS%MP1),&
 & LDV11=YTKE%LACTIVE,PI11=TENDENCY_CML%TKE,           PO11=TENDENCY_DYN%TKE)


!  ---------------------------------
!  Releasing the allocated space for local structures
IF ( (NCHEM > 0 .OR. NGEMS > 0 .OR. NACTAERO > 0) .AND. LL_HRES ) THEN
  DEALLOCATE(ZCHEM2AER)
ENDIF

CALL LOCAL_ARRAYS_FIN(LLKEYS,AUXL,SURFL,PERTL,GEMSL)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CALLPAR',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE CALLPAR
