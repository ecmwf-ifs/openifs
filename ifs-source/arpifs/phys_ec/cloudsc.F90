! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!OPTIONS XOPT(HSFUN)
#ifdef RS6K
@PROCESS HOT(NOVECTOR) NOSTRICT
#endif
SUBROUTINE CLOUDSC &
 !---input
 & (YDECLDP,YDECUMF,YDEPHLI,YDERAD, YDEPHY, YDVDF, YDSPP_CONFIG,&
 & KIDIA,    KFDIA,    KLON,    KLEV, &
 & PTSPHY, &
 & PT, PQ, &
 & PTENDENCY_CML_T, PTENDENCY_CML_Q, PTENDENCY_CML_A, PTENDENCY_CML_CLD, &
 & PTENDENCY_LOC_T, PTENDENCY_LOC_Q, PTENDENCY_LOC_A, PTENDENCY_LOC_CLD, &
 & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, PDYNR, PDYNS,&
 & PHRSW,    PHRLW, &
 & PVERVEL,  PAP,      PAPH, &
 & PLSM,     PGAW,     LDCUM,    KCTOP,    KTYPE,    KPBLTYPE, PEIS,&
 & PLU,      PLUDE,    PLUDELI,  PSNDE,    PMFU,     PMFD, PGP2DSPP, &
 & LDSLPHY, &
 !---prognostic fields
 & PA, &
 & PCLV,  &
!-- arrays for aerosol-cloud interactions
!!! & PQAER,    KAER, &
 & PLCRIT_AER,PICRIT_AER, &
 & PRE_ICE, &
 & PCCN,     PNICE, &
 !---diagnostic output
 & PCOVPTOT, PFSD, PRAINFRAC_TOPRFZ, &
 !---resulting fluxes
 & PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG, &
 & PFSQRF,   PFSQSF ,  PFCQRNG,  PFCQSNG, &
 & PFSQLTUR, PFSQITUR , &
 & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN, &
 & PEXTRA,   KFLDX)  

!===============================================================================
!**** *CLOUDSC* -  ROUTINE FOR PARAMETRIZATION OF CLOUD PROCESSES
!                  FOR PROGNOSTIC CLOUD SCHEME
!!
!     M.Tiedtke, C.Jakob, A.Tompkins, R.Forbes     (E.C.M.W.F.)
!!
!     PURPOSE
!     -------
!          THIS ROUTINE UPDATES THE CLOUD AND PRECIP FIELDS.
!
!        1. Initial set up
!        2. Tidy up of input values
!        3. Subgrid sources/sinks (convection
!           3.1 Calculate saturation quantities
!           3.2 Detrainment of cloud water from convective updrafts
!           3.3 Vertical advection due to convective subsidence, and 
!               subsequent evaporation due to adiabatic warming 
!           3.4 Erosion at cloud edges by turbulent mixing of cloud air
!               with unsaturated environmental air
!           3.5 Evaporation/condensation of cloud water in connection
!               with heating/cooling such as by subsidence/ascent
!        4. Microphysical processes
!           4.1 Sedimentation of rain, snow and ice
!           4.2 Define subgrid precipiation fractions
!           4.3 Autoconversion of cloud water into rain (collision-coalescence)
!           4.4 Evaporation of rain
!           4.5 Deposition onto ice when liquid water present (Bergeron-Findeison) 
!           4.6 Sublimation of ice 
!           4.7 Deposition onto snow when supersaturated 
!           4.8 Autoconversion of cloud ice to snow (aggregation)
!           4.9 Riming of snow
!           4.10 Snow accretes rain
!           4.11 Rain accretes snow
!           4.12 Melting of snow and ice
!           4.13 Freezing of rain
!           4.14 Freezing of cloud
!           4.15 Sublimation of snow
!
!        5. Implicit solver for all processes
!
!        Note: Turbulent transports of s,q,u,v at cloud tops due to
!           buoyancy fluxes and lw radiative cooling are treated in 
!           the VDF scheme
!!
!     INTERFACE.
!     ----------
!     *CLOUDSC* is called from *CLOUD_LAYER* called from *CALLPAR*
!     The routine takes its input from the prognostic variables:
!     t,q,l,i,a and detrainment of cloud water from the
!     convective clouds (massflux convection scheme)
!     it returns its output to:
!      1.modified tendencies of model variables t and q
!        as well as cloud liquid, ice, rain, snow and cloud fraction
!      2.generates precipitation fluxes from grid-scale clouds
!!
!     EXTERNALS.
!     ----------
!          NONE
!!
!     MODIFICATIONS.
!     -------------
!      M. TIEDTKE    E.C.M.W.F.     8/1988, 2/1990
!     CH. JAKOB      E.C.M.W.F.     2/1994 IMPLEMENTATION INTO IFS
!     A.TOMPKINS     E.C.M.W.F.     2002   NEW NUMERICS
!        01-05-22 : D.Salmond   Safety modifications
!        02-05-29 : D.Salmond   Optimisation
!        03-01-13 : J.Hague     MASS Vector Functions  J.Hague
!        03-10-01 : M.Hamrud    Cleaning
!        04-12-14 : A.Tompkins  New implicit solver and physics changes
!        04-12-03 : A.Tompkins & M.Ko"hler  moist PBL
!     G.Mozdzynski  09-Jan-2006  EXP security fix
!        19-01-09 : P.Bechtold  Changed increased RCLDIFF value for KTYPE=2
!        07-07-10 : A.Tompkins/R.Forbes  4-Phase flexible microphysics
!        01-03-11 : R.Forbes    Mixed phase changes and tidy up
!        01-10-11 : R.Forbes    Melt ice to rain, allow rain to freeze
!        01-10-11 : R.Forbes    Limit supersat to avoid excessive values
!        31-10-11 : M.Ahlgrimm  Add rain, snow and PEXTRA to DDH output
!        17-02-12 : F.Vana      Simplified/optimized LU factorization
!        18-05-12 : F.Vana      Cleaning + better support of sequential physics
!        N.Semane+P.Bechtold     04-10-2012 Add RVRFACTOR factor for small planet
!        01-02-13 : R.Forbes    New params of autoconv/acc,rain evap,snow riming
!        15-03-13 : F. Vana     New dataflow + more tendencies from the first call
!        K. Yessad (July 2014): Move some variables.
!        F. Vana  05-Mar-2015  Support for single precision
!        15-01-15 : R.Forbes    Added new options for snow evap & ice deposition
!        10-01-15 : R.Forbes    New physics for rain freezing
!        23-10-14 : P. Bechtold remove zeroing of convection arrays
!        15-12-2015 : R. Forbes Added inhomog option and variable rain fallspeed
!        27-01-2016 : M. Leutbecher & S.-J. Lock  Introduced SPP scheme (LSPP)
!        01-10-2016 : R. Forbes Tidy up routine
!        01-04-2017 : R. Forbes Modified numerics for rain/snow, 
!                               removed threshold for autoconv/accretion, 
!                               new turbulent erosion, ice evap, snow deposition
!        Oct-2017   : S.-J. Lock  Enabled options for new SPP microphysics perturbations
!        2017-11-11 M Ahlgrimm extent cloud heterogeneity for ice FSD
!        2019-01 : R.Forbes Passed in separate T/Q/A/CLD arrays for improved portability
!        2019-01 : R.Forbes Added additional diagnostics for cloud budget 
!        2019-01 : R.Forbes Added new parameters that can be set in the namelist
!        2019-01 : R.Forbes Tidy up and corrections to cloud budget diagnostics
!        2020-01 : R.Forbes Remove saturation adjustment, now done in CLOUD_SATADJ
!        2020-01 : R.Forbes New snow deposition and ice sublimation processes
!        2020-01 : R.Forbes Move ice deposition to after sedimentation
!        2020-01 : R.Forbes Change turbulent erosion to mixed phase scheme
!        2020-01 : R.Forbes Call CLOUD_SUPERSATCHECK to do final removal of supersaturation
!        2020-01 : R.Forbes Revised cloud budget diagnostics
!        2020-12 : R.Forbes Added PSACR and PRACS processes, modified fallspeed code
!        2020-10 : M.Leutbecher & S. Lang SPP abstraction and revision
!        2021-10 : R.Forbes Removed warm autoconv to snow, turned on PSACR to allow frz drz
!
!     REFERENCES.
!     ----------
!     Tietdke MWR 1993 - original description of the cloud parametrization
!     Jakob PhD 2000
!     Gregory et al. (2000) QJRMS
!     Tompkins el al. (2007) QJRMS - ice supersaturation parametrization
!     Forbes and Tompkins (2011) ECMWF Newsletter 129 - new prognostic liq/ice/rain/snow
!     Forbes et al. (2011) ECMWF Tech Memo 649 - new prognostic liq/ice/rain/snow
!     Forbes et al. (2014) ECMWF Newsletter 141 - freezing rain
!     Ahlgrimm and Forbes (2014) MWR - warm rain processes
!     Forbes and Ahlgrimm (2014) MWR - mixed-phase cloud processes
!     Ahlgrimm and Forbes (2015) MWR - subgrid heterogeneity of cloud condensate
!!
!===============================================================================

USE YOECLDP  , ONLY : TECLDP, NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
USE YOEPHLI  , ONLY : TEPHLI
USE YOERAD   , ONLY : TERAD
USE YOEPHY   , ONLY : TEPHY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : LSCMEC
USE YOMCST   , ONLY : RG, RD, RCPD, RETV, RLVTT, RLSTT, RTT, RV, RA, RPI  
USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 &                    R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE, RTICECU, &
 &                    RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2
USE YOECUMF  , ONLY : TECUMF
USE YOEVDF   , ONLY : TVDF
USE SPP_MOD     , ONLY : TSPP_CONFIG
USE SPP_GEN_MOD , ONLY : SPP_PERT

IMPLICIT NONE

!-------------------------------------------------------------------------------
!                 Declare input/output arguments
!-------------------------------------------------------------------------------
 
TYPE(TECLDP)      ,INTENT(INOUT) :: YDECLDP
TYPE(TECUMF)      ,INTENT(INOUT) :: YDECUMF
TYPE(TEPHLI)      ,INTENT(INOUT) :: YDEPHLI
TYPE(TERAD)       ,INTENT(INOUT) :: YDERAD
TYPE(TEPHY)       ,INTENT(INOUT) :: YDEPHY
TYPE(TVDF)        ,INTENT(IN)    :: YDVDF
TYPE(TSPP_CONFIG) ,INTENT(IN)    :: YDSPP_CONFIG
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON             ! Number of grid points
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV             ! Number of levels
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY           ! Physics timestep
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)    ! T at start of callpar
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV)    ! Q at start of callpar
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDENCY_CML_T(KLON,KLEV)   ! T cumulative tendency
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDENCY_CML_Q(KLON,KLEV)   ! Q cumulative tendency
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDENCY_CML_A(KLON,KLEV)   ! A cumulative tendency
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDENCY_CML_CLD(KLON,KLEV,NCLV) ! CLD cumulative tendency
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDENCY_LOC_T(KLON,KLEV)   ! T local output tendency
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDENCY_LOC_Q(KLON,KLEV)   ! Q local output tendency
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDENCY_LOC_A(KLON,KLEV)   ! A local output tendency
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDENCY_LOC_CLD(KLON,KLEV,NCLV) ! CLD local output tendency
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFA(KLON,KLEV)  ! CC from VDF scheme
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFL(KLON,KLEV)  ! Liq from VDF scheme
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFI(KLON,KLEV)  ! Ice from VDF scheme
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNA(KLON,KLEV) ! CC from Dynamics
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNL(KLON,KLEV) ! Liq tendency from Dynamics
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNI(KLON,KLEV) ! Ice tendency from Dynamics
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNR(KLON,KLEV) ! Rain tendency from Dynamics
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNS(KLON,KLEV) ! Snow tendency from Dynamics
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRSW(KLON,KLEV) ! Short-wave heating rate
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRLW(KLON,KLEV) ! Long-wave heating rate
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) ! Vertical velocity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)   ! Pressure on full levels
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1)! Pressure on half levels
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON)       ! Land fraction (0-1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAW(KLON)       ! Grid area=PGAW*4*RPI*RA**2
! From convection parametrization 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON)      ! Convection active
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON)      ! Convection level top
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE(KLON)      ! Convection type 0,1,2,3
INTEGER(KIND=JPIM),INTENT(IN)    :: KPBLTYPE(KLON)   ! BL type 0,1,2,3
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEIS(KLON)       ! PBL Inversion strength
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLU(KLON,KLEV)   ! Conv. condensate
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLUDE(KLON,KLEV) ! Conv. detrained water 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLUDELI(KLON,KLEV,4)! Conv. detrained liq/ice/vapor/T
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNDE(KLON,KLEV,2)! Conv. detrained snow/rain
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV)  ! Conv. mass flux up
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFD(KLON,KLEV)  ! Conv. mass flux down
! Options and cloud prognostic variables
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP2DSPP(KLON,YDSPP_CONFIG%SM%NRFTOTAL) ! perturbation pattern
LOGICAL           ,INTENT(IN)    :: LDSLPHY          ! True if semi-lag physics
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV)    ! Original Cloud fraction (t)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLV(KLON,KLEV,NCLV) ! Cloud/precip prognostics
! Prognostic aerosol
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLCRIT_AER(KLON,KLEV) ! critical liquid mmr for rain autoconversion process
REAL(KIND=JPRB)   ,INTENT(IN)    :: PICRIT_AER(KLON,KLEV) ! critical liquid mmr for snow autoconversion process
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE_ICE(KLON,KLEV)    ! ice effective radius
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCCN(KLON,KLEV)       ! liquid cloud condensation nuclei
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNICE(KLON,KLEV)      ! ice number concentration (cf. CCN) 
! Precipitation related
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOVPTOT(KLON,KLEV)   ! Precip fraction
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAINFRAC_TOPRFZ(KLON)! Rain/snow fraction at top of refreezing layer 
! Flux diagnostics for DDH budget
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQLF(KLON,KLEV+1)  ! Flux of liquid
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQIF(KLON,KLEV+1)  ! Flux of ice
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQLNG(KLON,KLEV+1) ! -ve corr for liq
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQNNG(KLON,KLEV+1) ! -ve corr for ice
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQRF(KLON,KLEV+1)  ! Flux diagnostics
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQSF(KLON,KLEV+1)  !    for DDH, generic
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQRNG(KLON,KLEV+1) ! rain
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQSNG(KLON,KLEV+1) ! snow
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQLTUR(KLON,KLEV+1) ! liquid flux due to VDF
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQITUR(KLON,KLEV+1) ! ice flux due to VDF
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLSL(KLON,KLEV+1) ! liq+rain sedim flux
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLSN(KLON,KLEV+1) ! ice+snow sedim flux
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPSL(KLON,KLEV+1) ! Enthalpy flux for liq
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPSN(KLON,KLEV+1) ! Enthalp flux for ice
! Extra fields for diagnostics
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(KLON,KLEV,KFLDX) ! extra fields
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX ! Number of extra fields
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSD(KLON,KLEV) ! cloud condensate fractional standard deviation

!-------------------------------------------------------------------------------
!                       Declare local variables
!-------------------------------------------------------------------------------

!  condensation and evaporation terms
REAL(KIND=JPRB) :: ZLEROS
! autoconversion terms
REAL(KIND=JPRB) :: ZRAINAUT(KLON), ZSNOWAUT(KLON)
REAL(KIND=JPRB) :: ZLIQCLD(KLON),  ZICECLD(KLON)

REAL(KIND=JPRB) :: ZFOEALFA(KLON,KLEV+1)
REAL(KIND=JPRB) :: ZICENUCLEI(KLON) ! number concentration of ice nuclei
REAL(KIND=JPRB) :: ZLICLD(KLON)
REAL(KIND=JPRB) :: ZAEROS
REAL(KIND=JPRB) :: ZLFINALSUM(KLON)
REAL(KIND=JPRB) :: ZDQS(KLON)
REAL(KIND=JPRB) :: ZDTGDP(KLON) 
REAL(KIND=JPRB) :: ZRDTGDP(KLON)  
REAL(KIND=JPRB) :: ZTRPAUS(KLON)
REAL(KIND=JPRB) :: ZCOVPCLR(KLON)   
REAL(KIND=JPRB) :: ZPRECLR
REAL(KIND=JPRB) :: ZCOVPTOT(KLON)    
REAL(KIND=JPRB) :: ZCOVPMAX(KLON)
REAL(KIND=JPRB) :: ZQPRETOT(KLON)
REAL(KIND=JPRB) :: ZDPEVAP
REAL(KIND=JPRB) :: ZDTFORC
REAL(KIND=JPRB) :: ZTP1(KLON,KLEV)   
REAL(KIND=JPRB) :: ZLDEFR(KLON)
REAL(KIND=JPRB) :: ZLDIFDT(KLON)
REAL(KIND=JPRB) :: ZDTGDPF(KLON)
REAL(KIND=JPRB) :: ZLCUST(KLON,NCLV)
REAL(KIND=JPRB) :: ZACUST(KLON)
REAL(KIND=JPRB) :: ZMF(KLON) 

REAL(KIND=JPRB) :: ZRHO(KLON)
REAL(KIND=JPRB) :: ZGDP(KLON)

! Accumulators of A,B,and C factors for cloud equations
REAL(KIND=JPRB) :: ZSOLAB(KLON) ! -ve implicit CC
REAL(KIND=JPRB) :: ZSOLAC(KLON) ! linear CC
REAL(KIND=JPRB) :: ZANEW
REAL(KIND=JPRB) :: ZANEWM1(KLON) 

REAL(KIND=JPRB) :: ZDA(KLON)
REAL(KIND=JPRB) :: ZLI(KLON,KLEV)
REAL(KIND=JPRB) :: ZA(KLON,KLEV)
REAL(KIND=JPRB) :: ZAORIG(KLON,KLEV) ! start of scheme value for CC

LOGICAL :: LLO1

INTEGER(KIND=JPIM) :: IK, JK, JL, JM, JN, JO, IS

REAL(KIND=JPRB) :: ZDP(KLON), ZPAPHD(KLON)

REAL(KIND=JPRB) :: ZALFA
REAL(KIND=JPRB) :: ZALFAW
REAL(KIND=JPRB) :: ZBETA,ZBETA1
!REAL(KIND=JPRB) :: ZBOTT
REAL(KIND=JPRB) :: ZCFPR
REAL(KIND=JPRB) :: ZCOR
REAL(KIND=JPRB) :: ZDENOM
REAL(KIND=JPRB) :: ZDPR
REAL(KIND=JPRB) :: ZDTDP
REAL(KIND=JPRB) :: ZE
REAL(KIND=JPRB) :: ZEPSEC
REAL(KIND=JPRB) :: ZFAC, ZFACI, ZFACW
REAL(KIND=JPRB) :: ZGDCP
REAL(KIND=JPRB) :: ZINEW
REAL(KIND=JPRB) :: ZLCRIT
REAL(KIND=JPRB) :: ZMFDN
REAL(KIND=JPRB) :: ZPRECIP
REAL(KIND=JPRB) :: ZQE
REAL(KIND=JPRB) :: ZQTMST, ZRDCP
REAL(KIND=JPRB) :: ZRHC, ZSIG, ZSIGK
REAL(KIND=JPRB) :: ZWTOT
REAL(KIND=JPRB) :: ZZCO, ZZRH, ZQADJ
REAL(KIND=JPRB) :: ZTNEW
REAL(KIND=JPRB) :: ZRG_R,ZGDPH_R,ZCONS1
REAL(KIND=JPRB) :: ZLFINAL
REAL(KIND=JPRB) :: ZMELT
REAL(KIND=JPRB) :: ZEVAP
REAL(KIND=JPRB) :: ZFRZ
REAL(KIND=JPRB) :: ZVPLIQ, ZVPICE
REAL(KIND=JPRB) :: ZADD, ZBDD, ZCVDS, ZICE0, ZDEPOS
REAL(KIND=JPRB) :: ZRE_ICE
REAL(KIND=JPRB) :: ZRLDCP
REAL(KIND=JPRB) :: ZDZ
REAL(KIND=JPRB) :: ZXRAMID
REAL(KIND=JPRB) :: ZQ_UPD(KLON)
REAL(KIND=JPRB) :: ZT_UPD(KLON)
REAL(KIND=JPRB) :: ZA_UPD(KLON)

! Controls T-dependence for liquid/ice production (1=mixedphase, 2=homogfrz)
INTEGER(KIND=JPIM) :: IFTLIQICE 

! A bunch of SPP variables 
LOGICAL            :: LLPERT_RAMID,  LLPERT_RCLDIFF,  LLPERT_RCLCRIT,  LLPERT_RLCRITSNOW
LOGICAL            :: LLPERT_RAINEVAP,  LLPERT_SNOWSUBLIM,  LLPERT_CLOUDINHOM
INTEGER(KIND=JPIM) :: IPRAMID,  IPRCLDIFF,  IPRCLCRIT,  IPRLCRITSNOW              ! SPP random field pointer
INTEGER(KIND=JPIM) :: IPRAINEVAP,  IPSNOWSUBLIM, IPCLOUDINHOMAUT, IPCLOUDINHOMACC ! SPP random field pointer
INTEGER(KIND=JPIM) :: IPN  ! SPP perturbation pointer
TYPE(SPP_PERT)     :: PN1RAMID, PN1RCLDIFF, PN1RCLCRIT, PN1RLCRITSNOW  ! SPP pertn. configs.
TYPE(SPP_PERT)     :: PN1RAINEVAP, PN1SNOWSUBLIM, PN1CLOUDINHOM        ! SPP pertn. configs.

INTEGER(KIND=JPIM) :: IPHASE(NCLV) ! marker for water phase of each species
                                   ! 0=vapour, 1=liquid, 2=ice

INTEGER(KIND=JPIM) :: IMELT(NCLV)  ! marks melting linkage for ice categories
                                   ! ice->liquid, snow->rain

LOGICAL :: LLFALL(NCLV)      ! marks falling species
                             ! LLFALL=0, cloud cover must > 0 for zqx > 0
                             ! LLFALL=1, no cloud needed, zqx can evaporate

REAL(KIND=JPRB) :: ZLIQFRAC(KLON,KLEV)  ! cloud liquid water fraction: ql/(ql+qi)
REAL(KIND=JPRB) :: ZICEFRAC(KLON,KLEV)  ! cloud ice water fraction: qi/(ql+qi)
REAL(KIND=JPRB) :: ZQX(KLON,KLEV,NCLV)  ! water variables
REAL(KIND=JPRB) :: ZQX0(KLON,KLEV,NCLV) ! water variables at start of scheme
REAL(KIND=JPRB) :: ZQXN(KLON,NCLV)      ! new values for zqx at time+1
REAL(KIND=JPRB) :: ZQXFG(KLON,NCLV)     ! first guess values including precip
REAL(KIND=JPRB) :: ZQXNM1(KLON,NCLV)    ! new values for zqx at time+1 at level above
REAL(KIND=JPRB) :: ZFLUXQ(KLON,NCLV)    ! fluxes convergence of species (needed?)

REAL(KIND=JPRB) :: ZPFPLSX(KLON,KLEV+1,NCLV) ! generalized precipitation flux
REAL(KIND=JPRB) :: ZLNEG(KLON,KLEV,NCLV)     ! for negative correction diagnostics
REAL(KIND=JPRB) :: ZMELTMAX(KLON)
REAL(KIND=JPRB) :: ZFRZMAX(KLON)
REAL(KIND=JPRB) :: ZICETOT(KLON)

REAL(KIND=JPRB) :: ZQXN2D(KLON,KLEV,NCLV)   ! water variables store

REAL(KIND=JPRB) :: ZQSMIX(KLON,KLEV) ! diagnostic mixed phase saturation 
REAL(KIND=JPRB) :: ZQSLIQ(KLON,KLEV) ! liquid water saturation
REAL(KIND=JPRB) :: ZQSICE(KLON,KLEV) ! ice water saturation

!REAL(KIND=JPRB) :: ZRHM(KLON,KLEV) ! diagnostic mixed phase RH
!REAL(KIND=JPRB) :: ZRHL(KLON,KLEV) ! RH wrt liq
!REAL(KIND=JPRB) :: ZRHI(KLON,KLEV) ! RH wrt ice

REAL(KIND=JPRB) :: ZFOEEWM(KLON,KLEV)
REAL(KIND=JPRB) :: ZFOEEW(KLON,KLEV)
REAL(KIND=JPRB) :: ZFOEELIQT(KLON,KLEV)

REAL(KIND=JPRB) :: ZDQSLIQDT(KLON), ZDQSICEDT(KLON), ZDQSMIXDT(KLON)
REAL(KIND=JPRB) :: ZCORQSLIQ(KLON)
REAL(KIND=JPRB) :: ZCORQSICE(KLON) 
REAL(KIND=JPRB) :: ZCORQSMIX(KLON)
REAL(KIND=JPRB) :: ZEVAPLIMLIQ(KLON), ZEVAPLIMICE(KLON), ZEVAPLIMMIX(KLON)
REAL(KIND=JPRB) :: ZT_ADJ(KLON) ! Supersat check temperature change
REAL(KIND=JPRB) :: ZQ_ADJ(KLON) ! Supersat check humidity change
REAL(KIND=JPRB) :: ZA_ADJ(KLON) ! Supersat check cloud fraction change
REAL(KIND=JPRB) :: ZL_ADJ(KLON) ! Supersat check cloud water change
REAL(KIND=JPRB) :: ZI_ADJ(KLON) ! Supersat check cloud ice change

!-------------------------------------------------------
! SOURCE/SINK array for implicit and explicit terms
!-------------------------------------------------------
! a POSITIVE value entered into the arrays is a...
!            Source of this variable
!            |
!            |   Sink of this variable
!            |   |
!            V   V
! ZSOLQA(JL,IQa,IQb)  = explicit terms
! ZSOLQB(JL,IQa,IQb)  = implicit terms
! Thus if ZSOLAB(JL,NCLDQL,IQV)=K where K>0 then this is 
! a source of NCLDQL and a sink of IQV
! put 'external' source terms such as PLUDE from 
! detrainment into explicit source/sink array diagnognal
! ZSOLQA(NCLDQL,NCLDQL)= -PLUDE
! i.e. a positive value is a sink! 
!-------------------------------------------------------

REAL(KIND=JPRB) :: ZSOLQA(KLON,NCLV,NCLV) ! explicit sources and sinks
REAL(KIND=JPRB) :: ZSOLQB(KLON,NCLV,NCLV) ! implicit sources and sinks
                        ! e.g. microphysical pathways between ice variables.
REAL(KIND=JPRB) :: ZQLHS(KLON,NCLV,NCLV)  ! n x n matrix storing the LHS of implicit solver
REAL(KIND=JPRB) :: ZVQX(NCLV)        ! fall speeds of three categories
REAL(KIND=JPRB) :: ZEXPLICIT

! REAL(KIND=JPRB) :: ZSINKSUM(KLON,NCLV)

! for sedimentation source/sink terms
REAL(KIND=JPRB) :: ZFALLSINK(KLON,NCLV)
REAL(KIND=JPRB) :: ZFALLSRCE(KLON,NCLV)
REAL(KIND=JPRB) :: ZFALLSPEED(KLON,NCLV)    ! Terminal fall velocity

! for convection detrainment source and subsidence source/sink terms
REAL(KIND=JPRB) :: ZCONVSRCE(KLON,NCLV)
REAL(KIND=JPRB) :: ZCONVSINK(KLON,NCLV)
REAL(KIND=JPRB) :: ZADVW(KLON), ZADVWD(KLON)

! Numerical fit to wet bulb temperature
REAL(KIND=JPRB),PARAMETER :: ZTW1 = 1329.31_JPRB
REAL(KIND=JPRB),PARAMETER :: ZTW2 = 0.0074615_JPRB
REAL(KIND=JPRB),PARAMETER :: ZTW3 = 0.85E5_JPRB
REAL(KIND=JPRB),PARAMETER :: ZTW4 = 40.637_JPRB
REAL(KIND=JPRB),PARAMETER :: ZTW5 = 275.0_JPRB

REAL(KIND=JPRB) :: ZSUBSAT  ! Subsaturation for snow melting term         
REAL(KIND=JPRB) :: ZTDMTW0  ! Diff between dry-bulb temperature and 
                            ! temperature when wet-bulb = 0degC 

! Variables for deposition term
REAL(KIND=JPRB) :: ZTCG ! Temperature dependent function for ice PSD
REAL(KIND=JPRB) :: ZFACX1I, ZFACX1S! PSD correction factor
REAL(KIND=JPRB) :: ZAPLUSB,ZCORRFAC,ZCORRFAC2,ZPR02,ZTERM1,ZTERM2 ! for ice dep
REAL(KIND=JPRB) :: ZCLDTOPDIST(KLON) ! Distance from cloud top
REAL(KIND=JPRB) :: ZINFACTOR         ! No. of ice nuclei factor for deposition
REAL(KIND=JPRB) :: ZOVERLAP_LIQICE   ! Overlap fraction between SLW and ice
REAL(KIND=JPRB) :: ZSUPERSATICE

! Option control variables
INTEGER(KIND=JPIM) :: IWARMRAIN
INTEGER(KIND=JPIM) :: IRAINACC
INTEGER(KIND=JPIM) :: ISUBLICE
INTEGER(KIND=JPIM) :: IDEPICE
INTEGER(KIND=JPIM) :: IDEPSNOW
INTEGER(KIND=JPIM) :: ISUBLSNOW
INTEGER(KIND=JPIM) :: IP_SNOW_ACCRETES_RAIN   ! Snow accretion of rain
INTEGER(KIND=JPIM) :: IP_RAIN_ACCRETES_SNOW   ! Rain accretion of snow
INTEGER(KIND=JPIM) :: IP_ICE_ACCRETES_RAIN    ! Ice accretion of rain
INTEGER(KIND=JPIM) :: IP_RAIN_ACCRETES_ICE    ! Rain accretion of ice
INTEGER(KIND=JPIM) :: IVARFALL
INTEGER(KIND=JPIM) :: ITURBEROSION

! Autoconversion/accretion/riming/evaporation
REAL(KIND=JPRB) :: ZRAINACC(KLON)
REAL(KIND=JPRB) :: ZRAINCLD(KLON)
REAL(KIND=JPRB) :: ZRAINCLDM1(KLON)
REAL(KIND=JPRB) :: ZSNOWRIME(KLON)
REAL(KIND=JPRB) :: ZSNOWCLD(KLON)
REAL(KIND=JPRB) :: ZSNOWCLDM1(KLON)
REAL(KIND=JPRB) :: ZESATLIQ
REAL(KIND=JPRB) :: ZFALLCORR
REAL(KIND=JPRB) :: ZLAMBDA
REAL(KIND=JPRB) :: ZRLAMBDA
REAL(KIND=JPRB) :: ZEVAP_DENOM
REAL(KIND=JPRB) :: ZCORR2
REAL(KIND=JPRB) :: ZKA
REAL(KIND=JPRB) :: ZCONST
REAL(KIND=JPRB) :: ZTEMP
REAL(KIND=JPRB) :: ZN0R(KLON)    ! N0 size distribution intercept for rain
!REAL(KIND=JPRB) :: ZN0R_MP(KLON) ! N0 size distribution intercept for rain(MP)
REAL(KIND=JPRB) :: ZN0S(KLON)    ! N0 size distribution intercept for snow
REAL(KIND=JPRB) :: ZN0I(KLON)    ! N0 size distribution intercept for ice
REAL(KIND=JPRB) :: ZLAMR(KLON)   ! Slope of rain particle size distribution
!REAL(KIND=JPRB) :: ZLAMR_MP(KLON)! Slope of rain particle size distribution (MP)
REAL(KIND=JPRB) :: ZLAMS(KLON)   ! Slope of snow particle size distribution
REAL(KIND=JPRB) :: ZLAMI(KLON)   ! Slope of ice particle size distribution
REAL(KIND=JPRB) :: ZPRACS(KLON)  ! Rain accretion of snow
REAL(KIND=JPRB) :: ZPSACR(KLON)  ! Snow accretion of rain
REAL(KIND=JPRB) :: ZPRACI(KLON)  ! Rain accretion of ice
REAL(KIND=JPRB) :: ZPIACR(KLON)  ! Ice accretion of rain
REAL(KIND=JPRB) :: ZEFF_PSACR    ! PSACR collection efficiency
REAL(KIND=JPRB) :: ZEFF_PRACS    ! PRACS collection efficiency
REAL(KIND=JPRB) :: ZEFF_PIACR    ! PIACR collection efficiency
REAL(KIND=JPRB) :: ZEFF_PRACI    ! PRACI collection efficiency
REAL(KIND=JPRB) :: ZFALLDIFF     ! Differential fallspeed
REAL(KIND=JPRB) :: ZDENSNOW      ! Snow density
REAL(KIND=JPRB) :: ZDENICE       ! Ice density

! Cloud and precipitation inhomogeneity
REAL(KIND=JPRB) :: ZGRIDLEN, ZCLDRAINCORR
REAL(KIND=JPRB) :: ZPHIC, ZFRACSDC, ZEAUT
REAL(KIND=JPRB) :: ZPHIR, ZFRACSDR, ZEACC
REAL(KIND=JPRB) :: ZPHIP1, ZPHIP2, ZPHIP3, ZQTOT

! Rain freezing
LOGICAL :: LLRAINLIQ(KLON)  ! True if majority of raindrops are liquid (no ice core)

! SCM budget statistics 
REAL(KIND=JPRB), ALLOCATABLE :: ZSUMQ0(:,:),  ZSUMQ1(:,:) , ZERRORQ(:,:), &
                               &ZSUMH0(:,:),  ZSUMH1(:,:) , ZERRORH(:,:)
REAL(KIND=JPRB) :: ZRAIN

! Cloud budget
REAL(KIND=JPRB) :: ZBUDCC(KLON,12) ! cloud fraction budget array
REAL(KIND=JPRB) :: ZBUDL(KLON,22)  ! cloud liquid budget array
REAL(KIND=JPRB) :: ZBUDI(KLON,18)  ! cloud ice budget array

REAL(KIND=JPRB) :: ZPSDEP(KLON)                         ! + Deposition of vapour to snow                                    
REAL(KIND=JPRB) :: ZPIEVAP(KLON)                        ! - Evaporation of precipitating ice

! Miscellaneous
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZTMPL,ZTMPI,ZTMPA
REAL(KIND=JPRB) :: ZMM,ZRR

! Cloud heterogeneity variables
! 'FSD' is the fractional standard deviation= standard deviation/mean
! of cloud condensate
REAL(KIND=JPRB) :: ZDELZ              ! Layer thickness in km
REAL(KIND=JPRB) :: ZPHI               ! Factor to account for layer thickness and total water in ice FSD calculation
REAL(KIND=JPRB) :: ZRATFSD(KLON,KLEV) ! detrainment ratio: detrained mass from convection scheme
                                      ! divided by total cloud mass liq/ice in the grid box
REAL(KIND=JPRB) :: ZLIQF(KLON,KLEV)   ! Final cloud liquid at end of cloudsc
REAL(KIND=JPRB) :: ZICEF(KLON,KLEV)   ! Final cloud ice at end of cloudsc
REAL(KIND=JPRB) :: ZLFSD(KLON,KLEV)   ! Liquid condensate FSD
REAL(KIND=JPRB) :: ZIFSD(KLON,KLEV)   ! Ice condensate FSD
REAL(KIND=JPRB) :: ZFSD(KLON,KLEV)    ! Merged FSD for liquid and ice, passed to radiation scheme
REAL(KIND=JPRB) :: ZANEWP(KLON,KLEV)  ! Modified cloud fraction at end of time step, with precip fraction, for ice FSD
REAL(KIND=JPRB) :: ZANEW2(KLON,KLEV)  ! Modified cloud fraction for liquid FSD (avoiding very small cloud fractions)
REAL(KIND=JPRB) :: ZIFSDBACK          ! Background ice FSD dependent on model grid scale

! Factor accounting for parameterized along-track variability being an underestimate of the true area variability.
! 1.3 taken from Hill et al. 2015
REAL(KIND=JPRB), PARAMETER :: ZR12=1.3_JPRB

REAL(KIND=JPRB), PARAMETER :: ZU1=0.446585_JPRB !Parameters used in background ice FSD calculation
REAL(KIND=JPRB), PARAMETER :: ZU2=0.308061_JPRB 
REAL(KIND=JPRB), PARAMETER :: ZU3=0.395736_JPRB
REAL(KIND=JPRB), PARAMETER :: ZU4=-0.744527_JPRB

LOGICAL            :: LLINDEX3(KLON,NCLV) ! index variable
INTEGER(KIND=JPIM) :: IORDER(KLON,NCLV), IORDV(NCLV) ! arrays for sorting explicit terms
REAL(KIND=JPRB) :: ZSINKSUM(KLON) 
REAL(KIND=JPRB) :: ZRATIO(KLON,NCLV), ZZRATIO, ZRAT, ZMAX
REAL(KIND=JPRB) :: ZEPSILON

#include "cloud_supersatcheck.intfb.h"

#include "abor1.intfb.h"

#include "fcttre.func.h"
#include "fccld.func.h"

!===============================================================================
IF (LHOOK) CALL DR_HOOK('CLOUDSC',0,ZHOOK_HANDLE)

ASSOCIATE(LAERICEAUTO=>YDECLDP%LAERICEAUTO, LAERICESED=>YDECLDP%LAERICESED, &
 & LAERLIQAUTOLSP=>YDECLDP%LAERLIQAUTOLSP, LAERLIQCOLL=>YDECLDP%LAERLIQCOLL, &
 & LBUD23=>YDEPHY%LBUD23, &
 & LCLDBUDC=>YDECLDP%LCLDBUDC, LCLDBUDL=>YDECLDP%LCLDBUDL, &
 & LCLDBUDI=>YDECLDP%LCLDBUDI, LCLDBUDT=>YDECLDP%LCLDBUDT, &
 & LCLDBUD_VERTINT=>YDECLDP%LCLDBUD_VERTINT, &
 & LCLDBUD_TIMEINT=>YDECLDP%LCLDBUD_TIMEINT, &
 & LCLDBUDGET=>YDECLDP%LCLDBUDGET, NCLDTOP=>YDECLDP%NCLDTOP, &
 & NSSOPT=>YDECLDP%NSSOPT, RAMID=>YDECLDP%RAMID, RAMIN=>YDECLDP%RAMIN, &
 & RCCN=>YDECLDP%RCCN, RCLCRIT_LAND=>YDECLDP%RCLCRIT_LAND, &
 & RCLCRIT_SEA=>YDECLDP%RCLCRIT_SEA, RCLDIFF=>YDECLDP%RCLDIFF, &
 & RCLDIFF_CONVI=>YDECLDP%RCLDIFF_CONVI, RCLDTOPCF=>YDECLDP%RCLDTOPCF, &
 & RCL_APB1=>YDECLDP%RCL_APB1, RCL_APB2=>YDECLDP%RCL_APB2, &
 & RCL_APB3=>YDECLDP%RCL_APB3, RCL_CDENOM1=>YDECLDP%RCL_CDENOM1, &
 & RCL_CDENOM2=>YDECLDP%RCL_CDENOM2, RCL_CDENOM3=>YDECLDP%RCL_CDENOM3, &
 & RCL_CONST1I=>YDECLDP%RCL_CONST1I, RCL_CONST1R=>YDECLDP%RCL_CONST1R, &
 & RCL_CONST1S=>YDECLDP%RCL_CONST1S, RCL_CONST2I=>YDECLDP%RCL_CONST2I, &
 & RCL_CONST2R=>YDECLDP%RCL_CONST2R, RCL_CONST2S=>YDECLDP%RCL_CONST2S, &
 & RCL_CONST3I=>YDECLDP%RCL_CONST3I, RCL_CONST3R=>YDECLDP%RCL_CONST3R, &
 & RCL_CONST3S=>YDECLDP%RCL_CONST3S, RCL_CONST4I=>YDECLDP%RCL_CONST4I, &
 & RCL_CONST4R=>YDECLDP%RCL_CONST4R, RCL_CONST4S=>YDECLDP%RCL_CONST4S, &
 & RCL_CONST5I=>YDECLDP%RCL_CONST5I, RCL_CONST5R=>YDECLDP%RCL_CONST5R, &
 & RCL_CONST5S=>YDECLDP%RCL_CONST5S, RCL_CONST6I=>YDECLDP%RCL_CONST6I, &
 & RCL_CONST6R=>YDECLDP%RCL_CONST6R, RCL_CONST6S=>YDECLDP%RCL_CONST6S, &
 & RCL_CONST7R=>YDECLDP%RCL_CONST7R, RCL_CONST7S=>YDECLDP%RCL_CONST7S, &
 & RCL_CONST8R=>YDECLDP%RCL_CONST8R, RCL_CONST8S=>YDECLDP%RCL_CONST8S, &
 & RCL_CONST9R=>YDECLDP%RCL_CONST9R, RCL_CONST10R=>YDECLDP%RCL_CONST10R, &
 & RCL_CONST7I=>YDECLDP%RCL_CONST7I, &
 & RCL_EFF_RACW=>YDECLDP%RCL_EFF_RACW, &
 & RCL_DR=>YDECLDP%RCL_DR, &
 & RCL_DI=>YDECLDP%RCL_DI, &
 & RCL_LAMBDA1I=>YDECLDP%RCL_LAMBDA1I, RCL_LAMBDA2I=>YDECLDP%RCL_LAMBDA2I, &
 & RCL_LAM1R=>YDECLDP%RCL_LAM1R, RCL_LAM2R=>YDECLDP%RCL_LAM2R, &
 & RCL_LAM1R_MP=>YDECLDP%RCL_LAM1R_MP, RCL_LAM2R_MP=>YDECLDP%RCL_LAM2R_MP, &
 & RCL_LAM1S=>YDECLDP%RCL_LAM1S, RCL_LAM2S=>YDECLDP%RCL_LAM2S, &
 & RCL_FZRAB=>YDECLDP%RCL_FZRAB, RCL_KA273=>YDECLDP%RCL_KA273, &
 & RCL_KKAAC=>YDECLDP%RCL_KKAAC, RCL_KKAAU=>YDECLDP%RCL_KKAAU, &
 & RCL_KKBAC=>YDECLDP%RCL_KKBAC, RCL_KKBAUN=>YDECLDP%RCL_KKBAUN, &
 & RCL_KKBAUQ=>YDECLDP%RCL_KKBAUQ, &
 & RCL_KK_CLOUD_NUM_LAND=>YDECLDP%RCL_KK_CLOUD_NUM_LAND, &
 & RCL_KK_CLOUD_NUM_SEA=>YDECLDP%RCL_KK_CLOUD_NUM_SEA, &
 & RCL_X3I=>YDECLDP%RCL_X3I, &
 & RCL_X1R=>YDECLDP%RCL_X1R, RCL_X2R=>YDECLDP%RCL_X2R, &
 & RCL_X1R_MP=>YDECLDP%RCL_X1R_MP, RCL_X2R_MP=>YDECLDP%RCL_X2R_MP, &
 & RCL_X1S=>YDECLDP%RCL_X1S, RCL_X2S=>YDECLDP%RCL_X2S, &
 & RCOVPMIN=>YDECLDP%RCOVPMIN, RDENSREF=>YDECLDP%RDENSREF,  &
 & RDENSWAT=>YDECLDP%RDENSWAT, RDEPLIQREFDEPTH=>YDECLDP%RDEPLIQREFDEPTH, &
 & RDEPLIQREFRATE=>YDECLDP%RDEPLIQREFRATE, RICEHI1=>YDECLDP%RICEHI1, &
 & RICEHI2=>YDECLDP%RICEHI2, RICEINIT=>YDECLDP%RICEINIT, RKCONV=>YDECLDP%RKCONV, &
 & RKOOPTAU=>YDECLDP%RKOOPTAU, RLCRITSNOW=>YDECLDP%RLCRITSNOW, &
 & RLMIN=>YDECLDP%RLMIN, RNICE=>YDECLDP%RNICE, RPECONS=>YDECLDP%RPECONS, &
 & RPRC1=>YDECLDP%RPRC1, RPRECRHMAX=>YDECLDP%RPRECRHMAX, &
 & RSNOWLIN1=>YDECLDP%RSNOWLIN1, RSNOWLIN2=>YDECLDP%RSNOWLIN2, &
 & RTAUMEL=>YDECLDP%RTAUMEL, RTHOMO=>YDECLDP%RTHOMO, RVICE=>YDECLDP%RVICE, &
 & RVRAIN=>YDECLDP%RVRAIN, RVRFACTOR=>YDECLDP%RVRFACTOR, &
 & RVSNOW=>YDECLDP%RVSNOW, LMFDSNOW=>YDECUMF%LMFDSNOW, &
 & RMFADVW=>YDECUMF%RMFADVW, RMFADVWDD=>YDECUMF%RMFADVWDD, REISTHSC=>YDVDF%REISTHSC, &
 & LCLOUD_INHOMOG=>YDECLDP%LCLOUD_INHOMOG, &
 & RCL_INHOMOGAUT    => YDECLDP%RCL_INHOMOGAUT, &
 & RCL_INHOMOGACC    => YDECLDP%RCL_INHOMOGACC, &
 & RCL_OVERLAPLIQICE => YDECLDP%RCL_OVERLAPLIQICE, &
 & RCL_EFFRIME       => YDECLDP%RCL_EFFRIME )
!===============================================================================


!######################################################################
!
!             0.  *** SET UP CONSTANTS ***
!
!######################################################################

! Define a small number
ZEPSILON=100._JPRB*EPSILON(ZEPSILON)

! Parameters that control amount of convective subsidence in dynamics
DO JL=KIDIA, KFDIA
  ZADVW(JL)=1.0_JPRB
  ZADVWD(JL)=1.0_JPRB
  IF(KTYPE(JL)==1.AND.RMFADVW>0) THEN
    ZADVW(JL) =1.0_JPRB-RMFADVW
    ZADVWD(JL)=1.0_JPRB-RMFADVW*RMFADVWDD
  ENDIF
ENDDO

! ---------------------------------------------------------------------
! LCLDBUD logicals store enthalpy and cloud water budgets 
! Default to .false. and read in from namelist NAMCLDP in sucldp.F90
! LCLDBUDC        - True = Turn on 3D cloud fraction process budget
! LCLDBUDL        - True = Turn on 3D cloud liquid process budget    
! LCLDBUDI        - True = Turn on 3D cloud ice process budget   
! LCLDBUDT        - True = Turn on 3D cloud process temperature budget   
! LCLDBUD_VERTINT - True = Turn on vertical integrated budget for all terms
! LCLDBUD_TIMEINT - True = Accumulate budget rather than instantaneous. 
!                          Applies to all terms above.
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Hardwired options for microphysical processes
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Set version of warm-rain autoconversion/accretion
! IWARMRAIN = 1 ! Sundquist
! IWARMRAIN = 2 ! Khairoutdinov and Kogan (2000) explicit
! IWARMRAIN = 3 ! Khairoutdinov and Kogan (2000) implicit
! ---------------------------------------------------------------------
IWARMRAIN = 3
! ---------------------------------------------------------------------
! Set version of warm-rain accretion
! Only active for IWARMRAIN = 3
! IRAINACC = 1 ! Khairoutdinov and Kogan (2000) implicit
! IRAINACC = 2 ! Collection equation
! ---------------------------------------------------------------------
IRAINACC = 1
! ---------------------------------------------------------------------
! Version of inhomogeneity parametrization, now switched by namelist logical:
! LCLOUD_INHOMOG = false ! Fixed values for inhomogeneity enhancement
! LCLOUD_INHOMOG = true  ! Parametrization based on Ahlgrimm and Forbes (2016)
! ---------------------------------------------------------------------
!
! ---------------------------------------------------------------------
! Set version of sublimation of snow
! ISUBLSNOW = 1 ! Sundquist
! ISUBLSNOW = 2 ! New
! ---------------------------------------------------------------------
ISUBLSNOW = 1
! ---------------------------------------------------------------------
! Set version of depositional growth of snow
! IDEPSNOW = 0 ! Process turned off
! IDEPSNOW = 1 ! On
! ---------------------------------------------------------------------
IDEPSNOW = 1
! ---------------------------------------------------------------------
! Set version of depositional growth of ice
! IDEPICE = 1 ! Rotstayn (2001)
! IDEPICE = 2 ! New
! ---------------------------------------------------------------------
IDEPICE = 1
! ---------------------------------------------------------------------
! Set version of sublimation of ice
! ISUBLICE = 0 ! None
! ISUBLICE = 1 ! New
! ---------------------------------------------------------------------
ISUBLICE = 0
! ---------------------------------------------------------------------
! Set version of rain and ice/snow collection
! ISACR = 0 ! None - no warm-rain freezing
! ISACR = 1 ! Parametrized scheme based on differential fall velocity
! ISACR = 2 ! Instantaneous freezing of warm-rain production to snow
! ---------------------------------------------------------------------
IP_SNOW_ACCRETES_RAIN = 1
IP_ICE_ACCRETES_RAIN  = 1
! ---------------------------------------------------------------------
! Set version of rain and ice/snow collection
! IRACS = 0 ! None
! IRACS = 1 ! Parametrized scheme based on differential fall velocity
! ---------------------------------------------------------------------
IP_RAIN_ACCRETES_SNOW = 0
IP_RAIN_ACCRETES_ICE  = 0
! ---------------------------------------------------------------------
! Set version of terminal fall speed
! IVARFALL = 0 ! Fixed fall speeds
! IVARFALL = 1 ! Variable fall speeds based on particle size distribution
! ---------------------------------------------------------------------
IVARFALL = 1
! ---------------------------------------------------------------------
! Set version of cloud edge erosion
! ITURBEROSION = 1 ! Original formulation
! ITURBEROSION = 2 ! Morcrette formulation
! ITURBEROSION = 3 ! Morcrette formulation with mixed phase assumption
! ---------------------------------------------------------------------
ITURBEROSION = 3

! ---------------------
! Some simple constants
! ---------------------
ZQTMST  = 1.0_JPRB/PTSPHY
ZGDCP   = RG/RCPD
ZRDCP   = RD/RCPD
ZEPSEC  = 1.E-14_JPRB
ZRG_R   = 1.0_JPRB/RG
ZRLDCP  = 1.0_JPRB/(RALSDCP-RALVDCP)

! Note: Defined in module/yoecldp.F90
! NCLDQL=1    ! liquid cloud water
! NCLDQI=2    ! ice cloud water
! NCLDQR=3    ! rain water
! NCLDQS=4    ! snow
! NCLDQV=5    ! vapour

!-----------------------------------
! Initialize value use in Ice FSD calculation
!-----------------------------------
ZIFSDBACK=0._JPRB

! -----------------------------------------------
! Define species phase, 0=vapour, 1=liquid, 2=ice
! -----------------------------------------------
IPHASE(NCLDQV)=0
IPHASE(NCLDQL)=1
IPHASE(NCLDQR)=1
IPHASE(NCLDQI)=2
IPHASE(NCLDQS)=2

! ---------------------------------------------------
! Set up melting/freezing index, 
! if an ice category melts/freezes, where does it go?
! ---------------------------------------------------
IMELT(NCLDQV)=-99
IMELT(NCLDQL)=NCLDQI
IMELT(NCLDQR)=NCLDQS
IMELT(NCLDQI)=NCLDQR
IMELT(NCLDQS)=NCLDQR

! -------------------------
! Set up fall speeds in m/s
! -------------------------
ZVQX(NCLDQV)=0.0_JPRB 
ZVQX(NCLDQL)=0.0_JPRB 
ZVQX(NCLDQI)=RVICE 
ZVQX(NCLDQR)=RVRAIN
ZVQX(NCLDQS)=RVSNOW
LLFALL(:)=.FALSE.
DO JM=1,NCLV
  IF (ZVQX(JM)>0.0_JPRB) LLFALL(JM)=.TRUE. ! falling species
ENDDO

! ------------------------------------------------
! Prepare parameter perturbations for SPP scheme
! ------------------------------------------------
!

!  Prepare SPP
IF (YDSPP_CONFIG%LSPP) THEN

  ! Critical relative humidity
  !
  IPN = YDSPP_CONFIG%PPTR%RAMID
  LLPERT_RAMID= IPN > 0
  IF (LLPERT_RAMID) THEN
    PN1RAMID  = YDSPP_CONFIG%SM%PN(IPN)
    IPRAMID   = PN1RAMID%MP
  ENDIF

  ! Turbulent erosion at cloud edges
  !
  IPN = YDSPP_CONFIG%PPTR%RCLDIFF
  LLPERT_RCLDIFF= IPN > 0
  IF (LLPERT_RCLDIFF) THEN
    PN1RCLDIFF  = YDSPP_CONFIG%SM%PN(IPN)
    IPRCLDIFF   = PN1RCLDIFF%MP
  ENDIF

  ! Liquid to rain autoconversion threshold
  !
  IPN = YDSPP_CONFIG%PPTR%RCLCRIT
  LLPERT_RCLCRIT= IPN > 0
  IF (LLPERT_RCLCRIT) THEN
    PN1RCLCRIT  = YDSPP_CONFIG%SM%PN(IPN)
    IPRCLCRIT   = PN1RCLCRIT%MP
  ENDIF

  ! Ice to snow autoconversion threshold
  !
  IPN = YDSPP_CONFIG%PPTR%RLCRITSNOW
  LLPERT_RLCRITSNOW= IPN > 0
  IF (LLPERT_RLCRITSNOW) THEN
    PN1RLCRITSNOW  = YDSPP_CONFIG%SM%PN(IPN)
    IPRLCRITSNOW   = PN1RLCRITSNOW%MP
  ENDIF

  ! Rain evaporation
  !
  IPN = YDSPP_CONFIG%PPTR%RAINEVAP
  LLPERT_RAINEVAP= IPN > 0
  IF (LLPERT_RAINEVAP) THEN
    PN1RAINEVAP  = YDSPP_CONFIG%SM%PN(IPN)
    IPRAINEVAP   = PN1RAINEVAP%MP
  ENDIF

  ! Snow sublimation
  !
  IPN = YDSPP_CONFIG%PPTR%SNOWSUBLIM
  LLPERT_SNOWSUBLIM= IPN > 0
  IF (LLPERT_SNOWSUBLIM) THEN
    PN1SNOWSUBLIM  = YDSPP_CONFIG%SM%PN(IPN)
    IPSNOWSUBLIM   = PN1SNOWSUBLIM%MP
  ENDIF

  ! Calculate inhomogeneity
  !
  IPN = YDSPP_CONFIG%PPTR%CLOUDINHOM
  LLPERT_CLOUDINHOM= IPN > 0
  IF (LLPERT_CLOUDINHOM) THEN
    PN1CLOUDINHOM  = YDSPP_CONFIG%SM%PN(IPN)
    IPCLOUDINHOMAUT= PN1CLOUDINHOM%MP
    IPCLOUDINHOMACC= PN1CLOUDINHOM%MP+1
  ENDIF
ELSE
  LLPERT_RAMID      =.FALSE.
  LLPERT_RCLDIFF    =.FALSE.
  LLPERT_RCLCRIT    =.FALSE.
  LLPERT_RLCRITSNOW =.FALSE.
  LLPERT_RAINEVAP   =.FALSE.
  LLPERT_SNOWSUBLIM =.FALSE.
  LLPERT_CLOUDINHOM =.FALSE.
ENDIF

!initialize local FSD arrays
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
     ZIFSD(JL,JK)=YDERAD%RCLOUD_FRAC_STD
     ZLFSD(JL,JK)=YDERAD%RCLOUD_FRAC_STD
     ZFSD(JL,JK)=YDERAD%RCLOUD_FRAC_STD
     ZLIQF(JL,JK)=0.0_JPRB
     ZICEF(JL,JK)=0.0_JPRB
     ZRATFSD(JL,JK)=0.0_JPRB
     ZANEW2(JL,JK)=0.0_JPRB
     ZANEWP(JL,JK)=0.0_JPRB
  ENDDO
ENDDO
!######################################################################
!
!             1.  *** INITIAL VALUES FOR VARIABLES ***
!
!######################################################################

! -----------------------------------------------
! Initialization of output tendencies
! -----------------------------------------------
DO JK=1,KLEV
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
    PTENDENCY_LOC_T(JL,JK)=0.0_JPRB
    PTENDENCY_LOC_Q(JL,JK)=0.0_JPRB
    PTENDENCY_LOC_A(JL,JK)=0.0_JPRB
  ENDDO
ENDDO
DO JM=1,NCLV-1
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PTENDENCY_LOC_CLD(JL,JK,JM)=0.0_JPRB
    ENDDO
  ENDDO
ENDDO

! ----------------------
! non CLV initialization 
! ----------------------
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZTP1(JL,JK)        = PT(JL,JK)+PTSPHY*PTENDENCY_CML_T(JL,JK)
    ZQX(JL,JK,NCLDQV)  = PQ(JL,JK)+PTSPHY*PTENDENCY_CML_Q(JL,JK) 
    ZQX0(JL,JK,NCLDQV) = PQ(JL,JK)+PTSPHY*PTENDENCY_CML_Q(JL,JK)
    ZA(JL,JK)          = PA(JL,JK)+PTSPHY*PTENDENCY_CML_A(JL,JK)
    ZAORIG(JL,JK)      = PA(JL,JK)+PTSPHY*PTENDENCY_CML_A(JL,JK)
  ENDDO
ENDDO

! -------------------------------------
! initialization for CLV family
! -------------------------------------
DO JM=1,NCLV-1
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZQX(JL,JK,JM)  = PCLV(JL,JK,JM)+PTSPHY*PTENDENCY_CML_CLD(JL,JK,JM)
      ZQX0(JL,JK,JM) = PCLV(JL,JK,JM)+PTSPHY*PTENDENCY_CML_CLD(JL,JK,JM)
    ENDDO
  ENDDO
ENDDO

!-------------
! zero arrays
!-------------
ZPFPLSX(:,:,:) = 0.0_JPRB ! precip fluxes
ZQXN2D(:,:,:)  = 0.0_JPRB ! end of timestep values in 2D
ZLNEG(:,:,:)   = 0.0_JPRB ! negative input check
PRAINFRAC_TOPRFZ(:) =0.0_JPRB ! rain fraction at top of refreezing layer
LLRAINLIQ(:) = .TRUE.  ! Assume all raindrops are liquid initially

!---------------------------------
! Find tropopause level (ZTRPAUS)
!---------------------------------
DO JL=KIDIA,KFDIA
  ZTRPAUS(JL)=0.1_JPRB
  ZPAPHD(JL)=1.0_JPRB/PAPH(JL,KLEV+1)
ENDDO
DO JK=1,KLEV-1
  DO JL=KIDIA,KFDIA
    ZSIG=PAP(JL,JK)*ZPAPHD(JL)
    IF (ZSIG>0.1_JPRB.AND.ZSIG<0.4_JPRB.AND.ZTP1(JL,JK)>ZTP1(JL,JK+1)) THEN
      ZTRPAUS(JL)=ZSIG
    ENDIF
  ENDDO
ENDDO

! -------------------------------------------
! Total water and enthalpy budget diagnostics
! -------------------------------------------
IF (LSCMEC.OR.LCLDBUDGET) THEN

  IF (.NOT. ALLOCATED(ZSUMQ0))   ALLOCATE(ZSUMQ0(KLON,KLEV))
  IF (.NOT. ALLOCATED(ZSUMQ1))   ALLOCATE(ZSUMQ1(KLON,KLEV))
  IF (.NOT. ALLOCATED(ZSUMH0))   ALLOCATE(ZSUMH0(KLON,KLEV))
  IF (.NOT. ALLOCATED(ZSUMH1))   ALLOCATE(ZSUMH1(KLON,KLEV))
  IF (.NOT. ALLOCATED(ZERRORQ))  ALLOCATE(ZERRORQ(KLON,KLEV))
  IF (.NOT. ALLOCATED(ZERRORH))  ALLOCATE(ZERRORH(KLON,KLEV))

  ! initialize the flux arrays
  DO JK=1,KLEV
!DIR$ IVDEP
    DO JL=KIDIA,KFDIA
      ZTNEW=PT(JL,JK)+PTSPHY*(PTENDENCY_LOC_T(JL,JK)+PTENDENCY_CML_T(JL,JK))
      IF (JK==1) THEN
        ZSUMQ0(JL,JK)=0.0_JPRB ! total water
        ZSUMH0(JL,JK)=0.0_JPRB ! liquid water temperature
      ELSE
        ZSUMQ0(JL,JK)=ZSUMQ0(JL,JK-1)
        ZSUMH0(JL,JK)=ZSUMH0(JL,JK-1)
      ENDIF

      ! Total for liquid
      ZTMPL = (PCLV(JL,JK,NCLDQL)+PCLV(JL,JK,NCLDQR) &
            &  +(PTENDENCY_LOC_CLD(JL,JK,NCLDQL)+ PTENDENCY_CML_CLD(JL,JK,NCLDQL) &
            &  + PTENDENCY_LOC_CLD(JL,JK,NCLDQR)+ PTENDENCY_CML_CLD(JL,JK,NCLDQR))*PTSPHY)
      ! Total for frozen
      ZTMPI = (PCLV(JL,JK,NCLDQI)+PCLV(JL,JK,NCLDQS) &
            &  +(PTENDENCY_LOC_CLD(JL,JK,NCLDQI)+ PTENDENCY_CML_CLD(JL,JK,NCLDQI) &
            &  + PTENDENCY_LOC_CLD(JL,JK,NCLDQS)+ PTENDENCY_CML_CLD(JL,JK,NCLDQS))*PTSPHY)
      ZTNEW = ZTNEW - RALVDCP*ZTMPL - RALSDCP*ZTMPI
      ZSUMQ0(JL,JK)=ZSUMQ0(JL,JK)*(ZTMPL+ZTMPI)*(PAPH(JL,JK+1)-PAPH(JL,JK))*ZRG_R

      ! detrained water treated here
      ZQE=PLUDE(JL,JK)*PTSPHY*RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
      IF (ZQE>RLMIN) THEN
        ZSUMQ0(JL,JK)=ZSUMQ0(JL,JK)+PLUDE(JL,JK)*PTSPHY
        ZTNEW=ZTNEW-(RALVDCP*PLUDELI(JL,JK,1)+RALSDCP*PLUDELI(JL,JK,2))*ZQE/PLUDE(JL,JK)
      ENDIF

      ZSUMH0(JL,JK)=ZSUMH0(JL,JK)+(PAPH(JL,JK+1)-PAPH(JL,JK))*ZTNEW 
      ZSUMQ0(JL,JK)=ZSUMQ0(JL,JK)+(PQ(JL,JK)+(PTENDENCY_LOC_Q(JL,JK)+PTENDENCY_CML_Q(JL,JK))* &
                    & PTSPHY)*(PAPH(JL,JK+1)-PAPH(JL,JK))*ZRG_R
    ENDDO
  ENDDO
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZSUMH0(JL,JK)=ZSUMH0(JL,JK)/PAPH(JL,JK+1)
    ENDDO
  ENDDO
ENDIF

! ------------------------------
! Define saturation values
! ------------------------------
DO JK=1,KLEV
   DO JL=KIDIA,KFDIA
      !----------------------------------------
      ! old *diagnostic* mixed phase saturation
      !---------------------------------------- 
      ZFOEALFA(JL,JK) = FOEALFA(ZTP1(JL,JK))
      ZFOEEWM(JL,JK)  = MIN(FOEEWM(ZTP1(JL,JK))/PAP(JL,JK),0.5_JPRB)
      ZQSMIX(JL,JK)   = ZFOEEWM(JL,JK)
      ZQSMIX(JL,JK)   = ZQSMIX(JL,JK)/(1.0_JPRB-RETV*ZQSMIX(JL,JK))

      !---------------------------------------------
      ! ice saturation T<273K
      ! liquid water saturation for T>273K 
      !---------------------------------------------
      ZALFA           = FOEDELTA(ZTP1(JL,JK))
      ZFOEEW(JL,JK)   = MIN((ZALFA*FOEELIQ(ZTP1(JL,JK))+ &
                        & (1.0_JPRB-ZALFA)*FOEEICE(ZTP1(JL,JK))) &
                        & /PAP(JL,JK),0.5_JPRB)
      ZFOEEW(JL,JK)   = MIN(0.5_JPRB,ZFOEEW(JL,JK))
      ZQSICE(JL,JK)   = ZFOEEW(JL,JK)/(1.0_JPRB-RETV*ZFOEEW(JL,JK))

      !----------------------------------
      ! liquid water saturation
      !---------------------------------- 
      ZFOEELIQT(JL,JK)= MIN(FOEELIQ(ZTP1(JL,JK))/PAP(JL,JK),0.5_JPRB)
      ZQSLIQ(JL,JK)   = ZFOEELIQT(JL,JK)
      ZQSLIQ(JL,JK)   = ZQSLIQ(JL,JK)/(1.0_JPRB-RETV*ZQSLIQ(JL,JK))
    ENDDO
ENDDO ! on JK


!######################################################################
!
!        2.       *** TIDY UP INPUT VALUES ***
!
!######################################################################

! ----------------------------------------------------
! Tidy up very small cloud cover or total cloud water
! ----------------------------------------------------
DO JK=1,KLEV
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
    IF (ZQX(JL,JK,NCLDQL)+ZQX(JL,JK,NCLDQI)<RLMIN.OR.ZA(JL,JK)<RAMIN) THEN

      ! Evaporate small cloud liquid water amounts
      ZLNEG(JL,JK,NCLDQL)  = ZLNEG(JL,JK,NCLDQL)+ZQX(JL,JK,NCLDQL)
      ZQADJ                = ZQX(JL,JK,NCLDQL)*ZQTMST
      PTENDENCY_LOC_Q(JL,JK)= PTENDENCY_LOC_Q(JL,JK)+ZQADJ
      PTENDENCY_LOC_T(JL,JK)= PTENDENCY_LOC_T(JL,JK)-RALVDCP*ZQADJ
      ZQX(JL,JK,NCLDQV)    = ZQX(JL,JK,NCLDQV)+ZQX(JL,JK,NCLDQL)
      ZQX(JL,JK,NCLDQL)    = 0.0_JPRB

      ! Evaporate small cloud ice water amounts
      ZLNEG(JL,JK,NCLDQI)  = ZLNEG(JL,JK,NCLDQI)+ZQX(JL,JK,NCLDQI)
      ZQADJ                = ZQX(JL,JK,NCLDQI)*ZQTMST
      PTENDENCY_LOC_Q(JL,JK)= PTENDENCY_LOC_Q(JL,JK)+ZQADJ
      PTENDENCY_LOC_T(JL,JK)= PTENDENCY_LOC_T(JL,JK)-RALSDCP*ZQADJ
      ZQX(JL,JK,NCLDQV)    = ZQX(JL,JK,NCLDQV)+ZQX(JL,JK,NCLDQI)
      ZQX(JL,JK,NCLDQI)    = 0.0_JPRB

      ! Set cloud cover to zero
      ZA(JL,JK)            = 0.0_JPRB

    ENDIF
  ENDDO
ENDDO

! ---------------------------------
! Tidy up small CLV variables
! ---------------------------------
!DIR$ IVDEP
DO JM=1,NCLV-1
!DIR$ IVDEP
  DO JK=1,KLEV
!DIR$ IVDEP
    DO JL=KIDIA,KFDIA
      IF (ZQX(JL,JK,JM)<RLMIN) THEN
        ZLNEG(JL,JK,JM)      = ZLNEG(JL,JK,JM)+ZQX(JL,JK,JM)
        ZQADJ                = ZQX(JL,JK,JM)*ZQTMST
        PTENDENCY_LOC_Q(JL,JK)= PTENDENCY_LOC_Q(JL,JK)+ZQADJ
        IF (IPHASE(JM)==1) PTENDENCY_LOC_T(JL,JK) = PTENDENCY_LOC_T(JL,JK)-RALVDCP*ZQADJ
        IF (IPHASE(JM)==2) PTENDENCY_LOC_T(JL,JK) = PTENDENCY_LOC_T(JL,JK)-RALSDCP*ZQADJ
        ZQX(JL,JK,NCLDQV)    = ZQX(JL,JK,NCLDQV)+ZQX(JL,JK,JM)
        ZQX(JL,JK,JM)        = 0.0_JPRB
      ENDIF
    ENDDO
  ENDDO
ENDDO

DO JK=1,KLEV
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA

    !------------------------------------------
    ! Ensure cloud fraction is between 0 and 1
    !------------------------------------------
    ZA(JL,JK)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZA(JL,JK)))

    !-------------------------------------------------------------------
    ! Calculate liq/ice fractions (no longer a diagnostic relationship)
    !-------------------------------------------------------------------
    ZLI(JL,JK)=ZQX(JL,JK,NCLDQL)+ZQX(JL,JK,NCLDQI)
    IF (ZLI(JL,JK)>RLMIN) THEN
      ZLIQFRAC(JL,JK)=ZQX(JL,JK,NCLDQL)/ZLI(JL,JK)
      ZICEFRAC(JL,JK)=1.0_JPRB-ZLIQFRAC(JL,JK)
    ELSE
      ZLIQFRAC(JL,JK)=0.0_JPRB
      ZICEFRAC(JL,JK)=0.0_JPRB
    ENDIF

  ENDDO
ENDDO

!-----------------------------
! Reset single level variables
!-----------------------------

ZQXNM1(:,:) = 0.0_JPRB
ZANEWM1(:)  = 0.0_JPRB
ZDA(:)      = 0.0_JPRB
ZCOVPCLR(:) = 0.0_JPRB
ZCOVPMAX(:) = 0.0_JPRB  
ZCOVPTOT(:) = 0.0_JPRB
ZCLDTOPDIST(:) = 0.0_JPRB

!----------------------------------------------------------------------
!
!                   START OF VERTICAL LOOP OVER LEVELS
!
!----------------------------------------------------------------------

DO JK=NCLDTOP,KLEV

!----------------------------------------------------------------------
! INITIALIZE VARIABLES
!----------------------------------------------------------------------

  !---------------------------------
  ! First guess microphysics
  !---------------------------------
  DO JM=1,NCLV
    DO JL=KIDIA,KFDIA
      ZQXFG(JL,JM)=ZQX(JL,JK,JM)
    ENDDO
  ENDDO

  !---------------------------------
  ! Set KLON arrays to zero
  !---------------------------------

  ZLICLD(:)   = 0.0_JPRB                                
  ZRAINAUT(:) = 0.0_JPRB  ! currently needed for diags  
  ZRAINACC(:) = 0.0_JPRB  ! currently needed for diags  
  ZSNOWAUT(:) = 0.0_JPRB  ! needed                      
  ZSNOWRIME(:) = 0.0_JPRB
  ZLDEFR(:)   = 0.0_JPRB                                
  ZACUST(:)   = 0.0_JPRB  ! set later when needed       
  ZQPRETOT(:) = 0.0_JPRB                                
  ZLFINALSUM(:)= 0.0_JPRB                               

  !-------------------------------------                
  ! solvers for cloud fraction                          
  !-------------------------------------                
  ZSOLAB(:) = 0.0_JPRB
  ZSOLAC(:) = 0.0_JPRB

  !------------------------------------------           
  ! reset matrix so missing pathways are set            
  !------------------------------------------           
  ZSOLQB(:,:,:) = 0.0_JPRB
  ZSOLQA(:,:,:) = 0.0_JPRB

  !----------------------------------                   
  ! reset new microphysics variables                    
  !----------------------------------                   
  ZFALLSPEED(:,:) = 0.0_JPRB
  ZFALLSRCE(:,:) = 0.0_JPRB
  ZFALLSINK(:,:) = 0.0_JPRB
  ZCONVSRCE(:,:) = 0.0_JPRB
  ZCONVSINK(:,:) = 0.0_JPRB
  ZICETOT(:)     = 0.0_JPRB                            

  ! Cloud budget arrays                                 
  ZBUDCC(:,:) = 0.0_JPRB                
  ZBUDL(:,:)  = 0.0_JPRB                 
  ZBUDI(:,:)  = 0.0_JPRB                 
  
  DO JL=KIDIA,KFDIA

    !-------------------------
    ! derived variables needed
    !-------------------------

    ZDP(JL)     = PAPH(JL,JK+1)-PAPH(JL,JK)     ! dp
    ZGDP(JL)    = RG/ZDP(JL)                    ! g/dp
    ZRHO(JL)    = PAP(JL,JK)/(RD*ZTP1(JL,JK))   ! p/RT air density
    ZDTGDP(JL)  = PTSPHY*ZGDP(JL)               ! dt g/dp
    ZRDTGDP(JL) = ZDP(JL)*(1.0_JPRB/(PTSPHY*RG))  ! 1/(dt g/dp)

    IF (JK>1) ZDTGDPF(JL) = PTSPHY*RG/(PAP(JL,JK)-PAP(JL,JK-1))

    !------------------------------------
    ! Calculate dqs/dT correction factor
    !------------------------------------
    ! Reminder: RETV=RV/RD-1
    
    ! liquid
    ZFACW         = R5LES/((ZTP1(JL,JK)-R4LES)**2)
    ZCOR          = 1.0_JPRB/(1.0_JPRB-RETV*ZFOEELIQT(JL,JK))
    ZDQSLIQDT(JL) = ZFACW*ZCOR*ZQSLIQ(JL,JK)
    ZCORQSLIQ(JL) = 1.0_JPRB+RALVDCP*ZDQSLIQDT(JL)

    ! ice
    ZFACI         = R5IES/((ZTP1(JL,JK)-R4IES)**2)
    ZCOR          = 1.0_JPRB/(1.0_JPRB-RETV*ZFOEEW(JL,JK))
    ZDQSICEDT(JL) = ZFACI*ZCOR*ZQSICE(JL,JK)
    ZCORQSICE(JL) = 1.0_JPRB+RALSDCP*ZDQSICEDT(JL)

    ! diagnostic mixed
    ZALFAW        = ZFOEALFA(JL,JK)
    ZFAC          = ZALFAW*ZFACW+(1.0_JPRB-ZALFAW)*ZFACI
    ZCOR          = 1.0_JPRB/(1.0_JPRB-RETV*ZFOEEWM(JL,JK))
    ZDQSMIXDT(JL) = ZFAC*ZCOR*ZQSMIX(JL,JK)
    ZCORQSMIX(JL) = 1.0_JPRB+FOELDCPM(ZTP1(JL,JK))*ZDQSMIXDT(JL)

    ! evaporation/sublimation limits
    ZEVAPLIMMIX(JL) = MAX((ZQSMIX(JL,JK)-ZQX(JL,JK,NCLDQV))/ZCORQSMIX(JL),0.0_JPRB)
    ZEVAPLIMLIQ(JL) = MAX((ZQSLIQ(JL,JK)-ZQX(JL,JK,NCLDQV))/ZCORQSLIQ(JL),0.0_JPRB)
    ZEVAPLIMICE(JL) = MAX((ZQSICE(JL,JK)-ZQX(JL,JK,NCLDQV))/ZCORQSICE(JL),0.0_JPRB)

    !--------------------------------
    ! in-cloud consensate amount
    !--------------------------------
    ZTMPA = 1.0_JPRB/MAX(ZA(JL,JK),0.01_JPRB)
    ZLIQCLD(JL) = ZQX(JL,JK,NCLDQL)*ZTMPA
    ZICECLD(JL) = ZQX(JL,JK,NCLDQI)*ZTMPA
    ZLICLD(JL)  = ZLIQCLD(JL)+ZICECLD(JL)

  ENDDO
  
  !------------------------------------------------
  ! Evaporate very small amounts of liquid and ice
  !------------------------------------------------
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA

    IF (ZQX(JL,JK,NCLDQL) < RLMIN) THEN
      ZSOLQA(JL,NCLDQV,NCLDQL) = ZQX(JL,JK,NCLDQL)
      ZSOLQA(JL,NCLDQL,NCLDQV) = -ZQX(JL,JK,NCLDQL)
    ENDIF

    IF (ZQX(JL,JK,NCLDQI) < RLMIN) THEN
      ZSOLQA(JL,NCLDQV,NCLDQI) = ZQX(JL,JK,NCLDQI)
      ZSOLQA(JL,NCLDQI,NCLDQV) = -ZQX(JL,JK,NCLDQI)
    ENDIF

  ENDDO


  !======================================================================
  !
  !
  !  3.2  DETRAINMENT FROM CONVECTION
  !
  !
  !======================================================================
  ! * Diagnostic T-ice/liq split retained for convection
  !    Note: This link is now flexible and a future convection 
  !    scheme can detrain explicit seperate budgets of:
  !    cloud water, ice, rain and snow
  ! * There is no (1-ZA) multiplier term on the cloud detrainment 
  !    term, since is now written in mass-flux terms  
  !---------------------------------------------------------------------

  IF (JK < KLEV .AND. JK>=NCLDTOP) THEN

!DEC$ IVDEP
    DO JL=KIDIA,KFDIA
    
      PLUDE(JL,JK)=PLUDE(JL,JK)*ZDTGDP(JL)

      IF(LDCUM(JL) .AND. PLUDE(JL,JK)>ZEPSEC .AND. PLU(JL,JK+1)>ZEPSEC) THEN
    
        ZSOLAC(JL)=ZSOLAC(JL)+PLUDE(JL,JK)/PLU(JL,JK+1)
        ZCONVSRCE(JL,NCLDQL) = PLUDELI(JL,JK,1)*ZDTGDP(JL)
        ZCONVSRCE(JL,NCLDQI) = PLUDELI(JL,JK,2)*ZDTGDP(JL)
        ZSOLQA(JL,NCLDQL,NCLDQL) = ZSOLQA(JL,NCLDQL,NCLDQL)+ZCONVSRCE(JL,NCLDQL)
        ZSOLQA(JL,NCLDQI,NCLDQI) = ZSOLQA(JL,NCLDQI,NCLDQI)+ZCONVSRCE(JL,NCLDQI)
        
        ! Store cloud budget diagnostics
        ZBUDL(JL,3)  = ZCONVSRCE(JL,NCLDQL)*ZQTMST
        ZBUDI(JL,3)  = ZCONVSRCE(JL,NCLDQI)*ZQTMST
        ZBUDCC(JL,3) = ZQTMST*PLUDE(JL,JK)/PLU(JL,JK+1)

      ELSE

        PLUDE(JL,JK)=0.0_JPRB
    
      ENDIF

        ! *convective snow/rain detrainment source
      IF(LMFDSNOW) THEN
        ZCONVSRCE(JL,NCLDQR) = PSNDE(JL,JK,1)*ZDTGDP(JL)
        ZCONVSRCE(JL,NCLDQS) = PSNDE(JL,JK,2)*ZDTGDP(JL)
        ZSOLQA(JL,NCLDQR,NCLDQR) = ZSOLQA(JL,NCLDQR,NCLDQR)+ZCONVSRCE(JL,NCLDQR)
        ZSOLQA(JL,NCLDQS,NCLDQS) = ZSOLQA(JL,NCLDQS,NCLDQS)+ZCONVSRCE(JL,NCLDQS)
      ENDIF
    
    ENDDO

  ENDIF ! JK<KLEV


  !======================================================================
  !
  !
  !  3.3  SUBSIDENCE COMPENSATING CONVECTIVE UPDRAUGHTS
  !
  !
  !======================================================================
  ! Three terms:
  ! * Convective subsidence source of cloud from layer above
  ! * Evaporation of cloud within the layer
  ! * Subsidence sink of cloud to the layer below (Implicit solution)
  !---------------------------------------------------------------------

  !-----------------------------------------------
  ! Subsidence source from layer above
  !               and 
  ! Evaporation of cloud within the layer
  !-----------------------------------------------
  IF (JK > NCLDTOP) THEN

    DO JL=KIDIA,KFDIA
      ZMF(JL)=MAX(0.0_JPRB,(PMFU(JL,JK)*ZADVW(JL)+PMFD(JL,JK)*ZADVWD(JL))*ZDTGDP(JL) )
      ZACUST(JL)=ZMF(JL)*ZANEWM1(JL)
    ENDDO

    DO JL=KIDIA,KFDIA
      ZLCUST(JL,NCLDQL) = ZMF(JL)*ZQXNM1(JL,NCLDQL)
      ZLCUST(JL,NCLDQI) = ZMF(JL)*ZQXNM1(JL,NCLDQI)
      ! record total flux for enthalpy budget:
      ZCONVSRCE(JL,NCLDQL) = ZCONVSRCE(JL,NCLDQL)+ZLCUST(JL,NCLDQL)
      ZCONVSRCE(JL,NCLDQI) = ZCONVSRCE(JL,NCLDQI)+ZLCUST(JL,NCLDQI)
    ENDDO

    ! Now have to work out how much liquid evaporates at arrival point 
    ! since there is no prognostic memory for in-cloud humidity, i.e. 
    ! we always assume cloud is saturated. 

    DO JL=KIDIA,KFDIA
      ZDTDP=ZRDCP*0.5_JPRB*(ZTP1(JL,JK-1)+ZTP1(JL,JK))/PAPH(JL,JK)
      ! Limit subsidence warming to layer_dp/dt (CFL=1) - is this necessary?
      ZWTOT   = MIN(RG*(PMFU(JL,JK)*ZADVW(JL)+PMFD(JL,JK)*ZADVWD(JL)),(PAP(JL,JK)-PAP(JL,JK-1))/PTSPHY)
      ZDTFORC = ZDTDP*PTSPHY*ZWTOT
      ![#Note: Diagnostic mixed phase should be replaced below]
      ZDQS(JL)=ZANEWM1(JL)*ZDTFORC*ZDQSMIXDT(JL)
    ENDDO

    ! Cloud liquid (NCLDQL)
    DO JL=KIDIA,KFDIA
      ZLFINAL=MAX(0.0_JPRB,ZLCUST(JL,NCLDQL)-ZDQS(JL)) !lim to zero
      ! no supersaturation allowed incloud ---V
      ZEVAP=MIN((ZLCUST(JL,NCLDQL)-ZLFINAL),ZEVAPLIMMIX(JL)) 
      ZLFINAL=ZLCUST(JL,NCLDQL)-ZEVAP 
      ZLFINALSUM(JL)=ZLFINALSUM(JL)+ZLFINAL ! sum 

      ZSOLQA(JL,NCLDQL,NCLDQL) = ZSOLQA(JL,NCLDQL,NCLDQL)+ZLCUST(JL,NCLDQL)
      ZSOLQA(JL,NCLDQV,NCLDQL) = ZSOLQA(JL,NCLDQV,NCLDQL)+ZEVAP
      ZSOLQA(JL,NCLDQL,NCLDQV) = ZSOLQA(JL,NCLDQL,NCLDQV)-ZEVAP
      ! Store cloud budget diagnostic
      ZBUDL(JL,4) = ZLCUST(JL,NCLDQL)*ZQTMST
      ZBUDL(JL,5) = -ZEVAP*ZQTMST
    ENDDO

    ! Cloud ice (NCLDQI)
    DO JL=KIDIA,KFDIA
      ZLFINAL=MAX(0.0_JPRB,ZLCUST(JL,NCLDQI)-ZDQS(JL)) !lim to zero
      ! no supersaturation allowed incloud ---V
      ZEVAP=MIN((ZLCUST(JL,NCLDQI)-ZLFINAL),ZEVAPLIMMIX(JL)) 
      ZLFINAL=ZLCUST(JL,NCLDQI)-ZEVAP 
      ZLFINALSUM(JL)=ZLFINALSUM(JL)+ZLFINAL ! sum 

      ZSOLQA(JL,NCLDQI,NCLDQI) = ZSOLQA(JL,NCLDQI,NCLDQI)+ZLCUST(JL,NCLDQI)
      ZSOLQA(JL,NCLDQV,NCLDQI) = ZSOLQA(JL,NCLDQV,NCLDQI)+ZEVAP
      ZSOLQA(JL,NCLDQI,NCLDQV) = ZSOLQA(JL,NCLDQI,NCLDQV)-ZEVAP
      ! Store cloud budget diagnostic
      ZBUDI(JL,4) = ZLCUST(JL,NCLDQI)*ZQTMST
      ZBUDI(JL,5) = -ZEVAP*ZQTMST 
    ENDDO
    
    !  Reset the cloud contribution if no cloud water survives to this level:
    DO JL=KIDIA,KFDIA
      ! Update cloud fraction
      ZSOLAC(JL) = ZSOLAC(JL)+ZACUST(JL)
      ! Store cloud fraction diagnostic if required
      ZBUDCC(JL,4) = ZACUST(JL)*ZQTMST
      !  Reset the cloud contribution if no cloud water survives to this level:
      IF (ZLFINALSUM(JL)<ZEPSEC) THEN
        ZSOLAC(JL) = ZSOLAC(JL)-ZACUST(JL)
        ZBUDCC(JL,5) = -ZACUST(JL)*ZQTMST
      ENDIF
    ENDDO

  ENDIF ! on  JK>NCLDTOP

  !---------------------------------------------------------------------
  ! Subsidence sink of cloud to the layer below 
  ! (Implicit - re. CFL limit on convective mass flux)
  !---------------------------------------------------------------------

  DO JL=KIDIA,KFDIA

    IF(JK<KLEV) THEN

      ZMFDN=MAX(0.0_JPRB,(PMFU(JL,JK+1)*ZADVW(JL)+PMFD(JL,JK+1)*ZADVWD(JL))*ZDTGDP(JL) )
      ZSOLAB(JL)=ZSOLAB(JL)+ZMFDN
      ZSOLQB(JL,NCLDQL,NCLDQL)=ZSOLQB(JL,NCLDQL,NCLDQL)+ZMFDN
      ZSOLQB(JL,NCLDQI,NCLDQI)=ZSOLQB(JL,NCLDQI,NCLDQI)+ZMFDN

      ! Record sink for cloud budget and enthalpy budget diagnostics
      ZCONVSINK(JL,NCLDQL) = ZMFDN
      ZCONVSINK(JL,NCLDQI) = ZMFDN

    ENDIF

  ENDDO


  !======================================================================
  !
  !
  ! 3.4  EROSION OF CLOUDS BY TURBULENT MIXING
  !
  !
  !======================================================================
  ! NOTE: This process decreases the cloud area 
  !       but leaves the specific cloud water content
  !       *within clouds* unchanged
  !----------------------------------------------------------------------

 IF (ITURBEROSION == 1) THEN

  ! ------------------------------
  ! Define turbulent erosion rate
  ! ------------------------------
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
    !original version (possibly perturbed by SPP)
    IF (LLPERT_RCLDIFF) THEN !Apply SPP perturbations
      ZLDIFDT(JL)=RCLDIFF*PTSPHY*EXP(PN1RCLDIFF%MU(1)+PN1RCLDIFF%XMAG(1)*PGP2DSPP(JL, IPRCLDIFF))
    ELSE
      ZLDIFDT(JL)=RCLDIFF*PTSPHY ! (unperturbed)
    ENDIF
    !Increase by factor of 5 for convective points
    IF(KTYPE(JL) > 0 .AND. PLUDE(JL,JK) > ZEPSEC) THEN
      IF(.NOT.(KTYPE(JL) >= 2)) &
       & ZLDIFDT(JL)=RCLDIFF_CONVI*ZLDIFDT(JL)  
    ENDIF
  ENDDO

  ! At the moment, works on mixed RH profile and partitioned ice/liq fraction
  ! so that it is similar to previous scheme
  ! Should apply RHw for liquid cloud and RHi for ice cloud separately 
  DO JL=KIDIA,KFDIA
    IF(ZLI(JL,JK) > ZEPSEC) THEN
      ! Calculate environmental humidity
!      ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSMIX(JL,JK))/&
!    &      MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))  
!      ZE=ZLDIFDT(JL)*MAX(ZQSMIX(JL,JK)-ZQE,0.0_JPRB)
      ZE=ZLDIFDT(JL)*MAX(ZQSMIX(JL,JK)-ZQX(JL,JK,NCLDQV),0.0_JPRB)
      ZLEROS=ZA(JL,JK)*ZE
      ZLEROS=MIN(ZLEROS,ZEVAPLIMMIX(JL))
      ZLEROS=MIN(ZLEROS,ZLI(JL,JK))
      ZAEROS=ZLEROS/ZLICLD(JL)  !if linear term

      ! Erosion is -ve LINEAR in L,A
      ZSOLAC(JL)=ZSOLAC(JL)-ZAEROS !linear

      ZSOLQA(JL,NCLDQV,NCLDQL) = ZSOLQA(JL,NCLDQV,NCLDQL)+ZLIQFRAC(JL,JK)*ZLEROS
      ZSOLQA(JL,NCLDQL,NCLDQV) = ZSOLQA(JL,NCLDQL,NCLDQV)-ZLIQFRAC(JL,JK)*ZLEROS
      ZSOLQA(JL,NCLDQV,NCLDQI) = ZSOLQA(JL,NCLDQV,NCLDQI)+ZICEFRAC(JL,JK)*ZLEROS
      ZSOLQA(JL,NCLDQI,NCLDQV) = ZSOLQA(JL,NCLDQI,NCLDQV)-ZICEFRAC(JL,JK)*ZLEROS

      ! Store cloud budget diagnostics if required
      ZBUDL(JL,7)  = -ZLIQFRAC(JL,JK)*ZLEROS*ZQTMST
      ZBUDI(JL,7)  = -ZICEFRAC(JL,JK)*ZLEROS*ZQTMST
      ZBUDCC(JL,7) = -ZAEROS*ZQTMST

    ENDIF
  ENDDO

 ELSEIF (ITURBEROSION == 2) THEN

  ! ------------------------------
  ! Define turbulent erosion rate
  ! Based on Morcrette
  ! Todo:
  ! To reduce calculations could remove duplication of some code for liq/ice
  ! Could also remove if tests as there is a limiter
  ! ZLICLD(JL) is LIQ+ICE/ZA
  ! Should this use ZCORQSLIQ(JL)?
  ! Have to use right subsaturation calculation otherwise duplicate evap
  ! ------------------------------
  DO JL=KIDIA,KFDIA

    ZLDIFDT(JL) = RCLDIFF*PTSPHY

    ! Increase erosion rate if convective
    !  KTYPE=1  Penetrative (deep) convection
    !  KTYPE=2  Shallow convection
    !  KTYPE=3  Mid-level convection
    ! Increase erosion for shallow/mid-level convection only
    IF(KTYPE(JL) >= 2 .AND. PLUDE(JL,JK) > ZEPSEC) THEN
       ZLDIFDT(JL) = RCLDIFF_CONVI*ZLDIFDT(JL)  
    ENDIF

  ENDDO

  ! Turbulent erosion of liquid
  DO JL=KIDIA,KFDIA
    ! If cloud liquid water is present
    IF (ZQX(JL,JK,NCLDQL) > ZEPSEC) THEN

      ! Calculate saturation deficit (wrt water)
      ZE = MAX(ZQSLIQ(JL,JK)-ZQX(JL,JK,NCLDQV),0.0_JPRB)

      ! Following Morcrette (2012)
      ZLEROS = 0.333_JPRB * ZA(JL,JK)*(1._JPRB-ZA(JL,JK))*ZLDIFDT(JL)*ZE
      
      ! Limiter taking account of evaporative cooling reducing saturation
      ZLEROS = MIN(ZLEROS,ZEVAPLIMLIQ(JL))
      
      ! Limiter to not remove more liquid than is present
      ZLEROS = MIN(ZLEROS,ZQX(JL,JK,NCLDQL))
      
      ! Cloud fraction decrease linearly proportional to liq+ice water decrease
      ZAEROS = ZLEROS/ZLICLD(JL)  !if linear term

      ! Erosion is -ve LINEAR in L,A
      ZSOLAC(JL) = ZSOLAC(JL)-ZAEROS !linear

      ! Update source/sink terms
      ZSOLQA(JL,NCLDQV,NCLDQL) = ZSOLQA(JL,NCLDQV,NCLDQL) + ZLEROS
      ZSOLQA(JL,NCLDQL,NCLDQV) = ZSOLQA(JL,NCLDQL,NCLDQV) - ZLEROS

      ! Store cloud budget diagnostics
      ZBUDL(JL,7)  = -ZLEROS*ZQTMST
      ZBUDCC(JL,7) = -ZAEROS*ZQTMST

    ENDIF
  ENDDO
  
  ! Turbulent erosion of ice
  DO JL=KIDIA,KFDIA
    ! If cloud ice is present
    IF (ZQX(JL,JK,NCLDQI) > ZEPSEC) THEN

      ! Calculate saturation deficit (wrt water)
      ZE = MAX(ZQSICE(JL,JK)-ZQX(JL,JK,NCLDQV),0.0_JPRB)

      ! Following Morcrette (2012)
      ZLEROS = 0.333_JPRB * ZA(JL,JK)*(1._JPRB-ZA(JL,JK))*ZLDIFDT(JL)*ZE
      
      ! Limiter taking account of evaporative cooling reducing saturation
      ZLEROS = MIN(ZLEROS,ZEVAPLIMICE(JL))
      
      ! Limiter to not remove more liquid than is present
      ZLEROS = MIN(ZLEROS,ZQX(JL,JK,NCLDQI))
      
      ! Cloud fraction decrease linearly proportional to liq+ice water decrease
      ZAEROS = ZLEROS/ZLICLD(JL)  !if linear term

      ! Erosion is -ve LINEAR in L,A
      ZSOLAC(JL) = ZSOLAC(JL)-ZAEROS !linear

      ! Update source/sink terms
      ZSOLQA(JL,NCLDQV,NCLDQI) = ZSOLQA(JL,NCLDQV,NCLDQI) + ZLEROS
      ZSOLQA(JL,NCLDQI,NCLDQV) = ZSOLQA(JL,NCLDQI,NCLDQV) - ZLEROS

      ! Store cloud budget diagnostics
      ZBUDI(JL,7)  = -ZLEROS*ZQTMST
      ZBUDCC(JL,7) = ZBUDCC(JL,7)-ZAEROS*ZQTMST

    ENDIF
  ENDDO
  
 ELSEIF (ITURBEROSION == 3) THEN

  ! ------------------------------
  ! Define turbulent erosion rate
  ! Based on Morcrette
  ! Todo:
  ! To reduce calculations could remove duplication of some code for liq/ice
  ! Could also remove if tests as there is a limiter
  ! ZLICLD(JL) is LIQ+ICE/ZA
  ! Should this use ZCORQSLIQ(JL)? No because it is the saturation deficit
  ! but do want to limit to max evaporation possible. Hence ZEVAPLIMLIQ
  ! Have to use right subsaturation calculation otherwise duplicate evap
  ! This version works on mixed phase
  ! ------------------------------
  DO JL=KIDIA,KFDIA

  !original version (possibly perturbed by SPP)
    IF (LLPERT_RCLDIFF) THEN
      ZLDIFDT(JL)=RCLDIFF*PTSPHY*EXP(PN1RCLDIFF%MU(1)+PN1RCLDIFF%XMAG(1)*PGP2DSPP(JL, IPRCLDIFF))
    ELSE
      ZLDIFDT(JL)=RCLDIFF*PTSPHY !original version (unperturbed)
    ENDIF

    ! Increase erosion rate if convective
    !  KTYPE=1  Penetrative (deep) convection
    !  KTYPE=2  Shallow convection
    !  KTYPE=3  Mid-level convection
    ! Increase erosion for shallow/mid-level convection only
    ! Alternative 1: IF(KTYPE(JL) >= 2 .AND. PLUDE(JL,JK) > ZEPSEC .AND. PEIS(JL) < REISTHSC) THEN
    ! Alternative 2: IF(KTYPE(JL) >= 2 .AND. JK >= (KCTOP(JL)-1.0_JPRB) .AND. PEIS(JL) < 10.0_JPRB) THEN
    IF(KTYPE(JL) >= 2 .AND. PLUDE(JL,JK) > ZEPSEC .AND. PEIS(JL) < 10.0_JPRB) THEN
       ZLDIFDT(JL) = 20.0_JPRB*RCLDIFF_CONVI*ZLDIFDT(JL)  
    ELSEIF(KTYPE(JL) == 1 .AND. PLUDE(JL,JK) > ZEPSEC) THEN
       ZLDIFDT(JL) = 2.0_JPRB*RCLDIFF_CONVI*ZLDIFDT(JL)  
    ENDIF

  ENDDO

  ! Turbulent erosion of liquid
  DO JL=KIDIA,KFDIA

    ! If cloud liquid water is present
    IF(ZLI(JL,JK) > ZEPSEC) THEN

      ! Calculate environmental humidity
      ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSMIX(JL,JK))/ &
    &      MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))

      ! Calculate saturation deficit (wrt mixed phase)
      ZE = MAX(ZQSMIX(JL,JK)-ZQE,0.0_JPRB)
            
      ! Following Morcrette (2012)
      ZLEROS = 0.333_JPRB * ZA(JL,JK)*(1._JPRB-ZA(JL,JK))*ZLDIFDT(JL)*ZE
      
      ! Limiter taking account of evaporative cooling reducing saturation
      ZLEROS = MIN(ZLEROS,ZEVAPLIMMIX(JL))
      
      ! Limiter to not remove more condensate than is present
      ZLEROS = MIN(ZLEROS,ZLI(JL,JK))
         
      ! Cloud fraction decrease linearly proportional to liq+ice water decrease
      ! but reduces cloud cover less quickly than water content
      ZAEROS = 0.5_JPRB*ZLEROS/ZLICLD(JL)
       
      ! Erosion is -ve LINEAR in L,A
      ZSOLAC(JL) = ZSOLAC(JL)-ZAEROS !linear

      ! Update source/sink terms
      ZSOLQA(JL,NCLDQV,NCLDQL) = ZSOLQA(JL,NCLDQV,NCLDQL)+ZLIQFRAC(JL,JK)*ZLEROS
      ZSOLQA(JL,NCLDQL,NCLDQV) = ZSOLQA(JL,NCLDQL,NCLDQV)-ZLIQFRAC(JL,JK)*ZLEROS
      ZSOLQA(JL,NCLDQV,NCLDQI) = ZSOLQA(JL,NCLDQV,NCLDQI)+ZICEFRAC(JL,JK)*ZLEROS
      ZSOLQA(JL,NCLDQI,NCLDQV) = ZSOLQA(JL,NCLDQI,NCLDQV)-ZICEFRAC(JL,JK)*ZLEROS

      ! Store cloud budget diagnostics if required
      ZBUDL(JL,7)  = -ZLIQFRAC(JL,JK)*ZLEROS*ZQTMST
      ZBUDI(JL,7)  = -ZICEFRAC(JL,JK)*ZLEROS*ZQTMST
      ZBUDCC(JL,7) = -ZAEROS*ZQTMST

    ENDIF
  ENDDO

 ENDIF ! on ITURBEROSION


  !############################################################################
  !#                                                                          #
  !#                                                                          #
  !#                     4.  MICROPHYSICAL PROCESSES                          #
  !#                                                                          #
  !#                                                                          #
  !############################################################################



  !======================================================================
  !
  !
  ! 4.1 SEDIMENTATION/FALLING OF HYDROMETEORS
  !
  !
  !======================================================================
  ! All hydrometeors fall downwards (can be carried upwards by dynamics)
  !----------------------------------------------------------------------
  
  !----------------------------------------------------
  ! Calculate hydrometeor terminal fall velocities
  !----------------------------------------------------
 
  DO JL=KIDIA,KFDIA

    ! Default to constant fall speed
    ZFALLSPEED(JL,NCLDQL) = ZVQX(NCLDQL)
    ZFALLSPEED(JL,NCLDQI) = ZVQX(NCLDQI)
    ZFALLSPEED(JL,NCLDQS) = ZVQX(NCLDQS)

    !-----------
    ! Cloud ice 
    !-----------
    ! if aerosol effect then override 
    !  note that for T>233K this is the same as above.
    IF (LAERICESED) THEN
      ZRE_ICE=PRE_ICE(JL,JK) 
      ! The exponent value is from 
      ! Morrison et al. JAS 2005 Appendix
      ZFALLSPEED(JL,NCLDQI) = 0.002_JPRB*ZRE_ICE**1.0_JPRB
    ENDIF

    ! Option to modify as fn(p,T) Heymsfield and Iaquinta JAS 2000
    !  ZFALL(JL,JM) = ZFALL(JL,JM)*((PAP(JL,JK)*RICEHI1)**(-0.178_JPRB)) &
    !             &*((ZTP1(JL,JK)*RICEHI2)**(-0.394_JPRB))
       
    !---------------------------------
    ! Rain
    !---------------------------------
    IF (IVARFALL == 1) THEN  ! if variable fallspeed on

      ! Fallspeed air density correction 
      ZFALLCORR = (RDENSREF/ZRHO(JL))**0.4_JPRB
     
      ! Rain water content (kg/kg) simple average between this layer and layer above
      ZTEMP = 0.5_JPRB*(ZQX(JL,JK,NCLDQR)+ZQXNM1(JL,NCLDQR))          

      ! Reciprocal of slope of particle size distribution
      ZRLAMBDA = (ZRHO(JL)*ZTEMP/RCL_LAM1R)**RCL_LAM2R 
  
      ! Calculate fallspeed
      ZFALLSPEED(JL,NCLDQR) = ZFALLCORR*RCL_CONST7R*ZRLAMBDA**RCL_DR
    
    ELSE
 
      ZFALLSPEED(JL,NCLDQR) = ZVQX(NCLDQR)
    
    ENDIF
 
  ENDDO


  !-------------------------------------------------------------------
  ! Calculate amount falling in (explicit) and falling out (implicit)
  ! of grid box. Loop over sedimenting hydrometeors
  !-------------------------------------------------------------------

  DO JM = 1,NCLV
    IF (LLFALL(JM)) THEN
      DO JL=KIDIA,KFDIA
        !------------------------
        ! source from layer above 
        !------------------------
        IF (JK > NCLDTOP) THEN
          ZFALLSRCE(JL,JM) = ZPFPLSX(JL,JK,JM)*ZDTGDP(JL) 
          ZSOLQA(JL,JM,JM) = ZSOLQA(JL,JM,JM)+ZFALLSRCE(JL,JM)
          ZQXFG(JL,JM)     = ZQXFG(JL,JM)+ZFALLSRCE(JL,JM)
          ! use first guess precip----------V
          ZQPRETOT(JL)     = ZQPRETOT(JL)+ZQXFG(JL,JM) 
          IF (JM == NCLDQI) THEN
            ZBUDI(JL,12)=ZFALLSRCE(JL,JM)*ZQTMST
          ENDIF
        ENDIF

        !-------------------------------------------------
        ! Calculate sink to next layer (implicit)
        !-------------------------------------------------
        ZFALLSINK(JL,JM)=ZDTGDP(JL)*ZRHO(JL)*ZFALLSPEED(JL,JM)
        ! Cloud budget diagnostic stored at end as implicit
      
      ENDDO ! jl  
    ENDIF ! LLFALL
  ENDDO ! End loop over hydrometeor type


  !======================================================================
  !
  !
  ! 4.2 DEFINE PRECIPITATION GRIDBOX FRACTION
  !
  !
  !======================================================================
  ! Although precipitation (rain/snow) are prognostic variables
  ! precipitation fraction is diagnostic, so needs to be calculated 
  ! using the prognostic cloud cover working down through each grid 
  ! column every timestep. Maximum-random overlap is assumed
  ! Since precipitation may be advected into a column with no cloud above
  ! an arbitrary minimum coverage if precip>0 is assumed (RCOVPMIN).
  ! Since there is no memory of the clear sky precip fraction, 
  ! the precipitation cover ZCOVPTOT, which has the memory in the column,
  ! is reduced proportionally with the precip evaporation rate.
  !---------------------------------------------------------------

  DO JL=KIDIA,KFDIA
    IF (ZQPRETOT(JL) > ZEPSEC) THEN
      ZCOVPTOT(JL)   = 1.0_JPRB - ((1.0_JPRB-ZCOVPTOT(JL))* &
       &              (1.0_JPRB - MAX(ZA(JL,JK),ZA(JL,JK-1)))/ &
       &              (1.0_JPRB - MIN(ZA(JL,JK-1),1.0_JPRB-1.E-06_JPRB)) )  
      ZCOVPTOT(JL)   = MAX(ZCOVPTOT(JL),RCOVPMIN)
      ZCOVPCLR(JL)   = MAX(0.0_JPRB,ZCOVPTOT(JL)-ZA(JL,JK)) ! clear sky proportion
      ZRAINCLD(JL)   = 0.5*(ZQX(JL,JK,NCLDQR)+ZQX(JL,JK-1,NCLDQR))/ZCOVPTOT(JL)
      ZRAINCLDM1(JL) = 0.5*(ZQX(JL,JK,NCLDQR)+ZQXNM1(JL,NCLDQR))/ZCOVPTOT(JL)
      ZSNOWCLD(JL)   = 0.5*(ZQX(JL,JK,NCLDQS)+ZQX(JL,JK-1,NCLDQS))/ZCOVPTOT(JL)
      ZSNOWCLDM1(JL) = 0.5*(ZQX(JL,JK,NCLDQS)+ZQXNM1(JL,NCLDQS))/ZCOVPTOT(JL)
      ZCOVPMAX(JL)   = MAX(ZCOVPTOT(JL),ZCOVPMAX(JL))
    ELSE
      ZRAINCLD(JL)   = 0.0_JPRB ! no precip
      ZRAINCLDM1(JL) = 0.0_JPRB ! no precip
      ZSNOWCLD(JL)   = 0.0_JPRB ! no precip
      ZSNOWCLDM1(JL) = 0.0_JPRB ! no precip
      ZCOVPTOT(JL)   = 0.0_JPRB ! no flux - reset cover
      ZCOVPCLR(JL)   = 0.0_JPRB ! reset clear sky proportion 
      ZCOVPMAX(JL)   = 0.0_JPRB ! reset max cover for ZZRH calc 
    ENDIF
  ENDDO


  !======================================================================
  !
  !
  ! 4.3  AUTOCONVERSION OF LIQUID TO RAIN
  !
  !
  !======================================================================

  DO JL=KIDIA,KFDIA

   IF (ZLIQCLD(JL) > ZEPSEC) THEN

      !----------------------------------------
      ! Calculate inhomogeneity 
      !----------------------------------------
      IF (LCLOUD_INHOMOG .OR. YDEPHY%LRAD_CLOUD_INHOMOG) THEN

        ! Total water (vapour+liquid) in g/kg
        ZQTOT = MIN((ZQX(JL,JK,NCLDQV)+ZQX(JL,JK,NCLDQL))*1000._JPRB,30._JPRB)

        ! Representative grid box length (km)
        ! PGAW = normalised gaussian quadrature weight / no. longitude pts
        ZGRIDLEN = 2*RA*SQRT(RPI*PGAW(JL))*0.001_JPRB

        ! Correlation between cloud and rain
        ZCLDRAINCORR = 1._JPRB-0.8_JPRB*ZA(JL,JK)

        !Maike's new parameters for modified Boutle parameterization
        ZPHIP1=.2_JPRB+.01_JPRB*ZQTOT+.0027_JPRB*ZQTOT**2-.00008_JPRB*ZQTOT**3
        ZPHIP2 = ZPHIP1-0.2_JPRB
        ZPHIP3 = 0.123_JPRB*EXP(-ZPHIP1*0.55_JPRB)
        !cut off very low cloud fraction at 0.1  to avoid very low FSD values
        ZANEW2(JL,JK)=MAX(ZA(JL,JK),0.1_JPRB) ! cut off very low cloud fraction to avoid

        ! Autoconversion enhancement factor from Boutle et al. 2013
        ZPHIC = (ZGRIDLEN*ZANEW2(JL,JK))**(1._JPRB/3._JPRB)* &
     &      ((ZPHIP3*ZGRIDLEN*ZANEW2(JL,JK))**1.5_JPRB + 3._JPRB*ZPHIP3)**(-0.17_JPRB)

        ! Fractional standard deviation for cloud condensate
        ZFRACSDC = (ZPHIP1 - ZPHIP2*ZANEW2(JL,JK))*ZPHIC

        ! Unique value when essentially no cloud edges
        IF (ZA(JL,JK) > 0.95_JPRB) ZFRACSDC = 0.17_JPRB*ZPHIC

        ! multiply by 1D-to-2D variability enhancement factor
        ZLFSD(JL,JK)=ZR12*ZFRACSDC
       ELSE
          ! Default FSD value
          ZLFSD(JL,JK)=YDERAD%RCLOUD_FRAC_STD
       ENDIF 


    !--------------------------------------------------------
    !-
    !- Warm-rain process follow Sundqvist (1989)
    !-
    !--------------------------------------------------------
    ! Implicit in liquid water content 
    
    IF (IWARMRAIN == 1) THEN

      ZZCO=RKCONV*PTSPHY

      IF (LAERLIQAUTOLSP) THEN
        ZLCRIT=PLCRIT_AER(JL,JK)
        ! 0.3 = N**0.333 with N=125 cm-3 
        ZZCO=ZZCO*(RCCN/PCCN(JL,JK))**0.333_JPRB
      ELSE
        ! Modify autoconversion threshold dependent on: 
        !  land (polluted, high CCN, smaller droplets, higher threshold)
        !  sea  (clean, low CCN, larger droplets, lower threshold)
        IF (LLPERT_RCLCRIT) THEN  !Apply SPP perturbations
          IF (PLSM(JL) > 0.5_JPRB) THEN
             ! perturbed land value of RCLCRIT
            ZLCRIT = RCLCRIT_LAND*EXP(PN1RCLCRIT%MU(1)+PN1RCLCRIT%XMAG(1)*PGP2DSPP(JL, IPRCLCRIT))
          ELSE
            ! perturbed ocean value of RCLCRIT
            ZLCRIT = RCLCRIT_SEA *EXP(PN1RCLCRIT%MU(2)+PN1RCLCRIT%XMAG(2)*PGP2DSPP(JL, IPRCLCRIT))
          ENDIF
        ELSE
          IF (PLSM(JL) > 0.5_JPRB) THEN
            ZLCRIT = RCLCRIT_LAND ! land  (unperturbed)
          ELSE
            ZLCRIT = RCLCRIT_SEA  ! ocean (unperturbed)
          ENDIF
        ENDIF
      ENDIF 

      !------------------------------------------------------------------
      ! Parameters for cloud collection by rain and snow.
      ! Note that with new prognostic variable it is now possible 
      ! to REPLACE this with an explicit collection parametrization
      !------------------------------------------------------------------   
      ZPRECIP=(ZPFPLSX(JL,JK,NCLDQS)+ZPFPLSX(JL,JK,NCLDQR))/MAX(ZEPSEC,ZCOVPTOT(JL))
      ZCFPR=1.0_JPRB + RPRC1*SQRT(MAX(ZPRECIP,0.0_JPRB))
!      ZCFPR=1.0_JPRB + RPRC1*SQRT(MAX(ZPRECIP,0.0_JPRB))* &
!       &ZCOVPTOT(JL)/(MAX(ZA(JL,JK),ZEPSEC))

      IF (LAERLIQCOLL) THEN 
        ! 5.0 = N**0.333 with N=125 cm-3 
        ZCFPR=ZCFPR*(RCCN/PCCN(JL,JK))**0.333_JPRB
      ENDIF

      ZZCO=ZZCO*ZCFPR
      ZLCRIT=ZLCRIT/MAX(ZCFPR,ZEPSEC)
  
      IF(ZLIQCLD(JL)/ZLCRIT < 20.0_JPRB )THEN ! Security for exp for some compilers
        ZRAINAUT(JL)=ZZCO*(1.0_JPRB-EXP(-(ZLIQCLD(JL)/ZLCRIT)**2))
      ELSE
        ZRAINAUT(JL)=ZZCO
      ENDIF

      ! rain freezes instantly if T<0C
      IF (IP_SNOW_ACCRETES_RAIN == 2) THEN
        IF(ZTP1(JL,JK) <= RTT) THEN
          ZSOLQB(JL,NCLDQS,NCLDQL)=ZSOLQB(JL,NCLDQS,NCLDQL)+ZRAINAUT(JL)
        ELSE
          ZSOLQB(JL,NCLDQR,NCLDQL)=ZSOLQB(JL,NCLDQR,NCLDQL)+ZRAINAUT(JL)
        ENDIF
      ELSE ! Don't freeze supercooled rain here
        ZSOLQB(JL,NCLDQR,NCLDQL)=ZSOLQB(JL,NCLDQR,NCLDQL)+ZRAINAUT(JL)
      ENDIF

    !--------------------------------------------------------
    !-
    !- Warm-rain process follow Khairoutdinov and Kogan (2000)
    !-
    !--------------------------------------------------------
    ELSEIF (IWARMRAIN == 2 .OR. IWARMRAIN == 3) THEN

      IF (LLPERT_RCLCRIT) THEN  !Apply SPP perturbations
        IF (PLSM(JL) > 0.5_JPRB) THEN
          ZCONST = RCL_KK_CLOUD_NUM_LAND
          ! perturbed land value of RCLCRIT
          ZLCRIT = RCLCRIT_LAND*EXP(PN1RCLCRIT%MU(1)+PN1RCLCRIT%XMAG(1)*PGP2DSPP(JL, IPRCLCRIT))
        ELSE
          ZCONST = RCL_KK_CLOUD_NUM_SEA
          ! perturbed ocean value of RCLCRIT
          ZLCRIT = RCLCRIT_SEA *EXP(PN1RCLCRIT%MU(2)+PN1RCLCRIT%XMAG(2)*PGP2DSPP(JL, IPRCLCRIT))
        ENDIF
      ELSE
        IF (PLSM(JL) > 0.5_JPRB) THEN ! land  (unperturbed)
          ZCONST = RCL_KK_CLOUD_NUM_LAND
          ZLCRIT = RCLCRIT_LAND
        ELSE                          ! ocean (unperturbed)
          ZCONST = RCL_KK_CLOUD_NUM_SEA
          ZLCRIT = RCLCRIT_SEA
        ENDIF
      ENDIF
 
      !----------------------------------------
      ! Calculate inhomogeneity 
      !----------------------------------------
      IF (LCLOUD_INHOMOG) THEN

         ZFRACSDC=ZLFSD(JL,JK)

         ZEAUT = (1._JPRB+ZFRACSDC**2._JPRB)**((RCL_KKBAUQ-1._JPRB) &
              &       *RCL_KKBAUQ/2._JPRB)

        ! Accretion enhancement factor from Boutle et al. 2013
        ZPHIR = (ZGRIDLEN*ZA(JL,JK))**(1._JPRB/3._JPRB)* &
              &       ((0.11_JPRB*ZGRIDLEN*ZA(JL,JK))**1.14_JPRB + 1._JPRB)**(-0.22_JPRB)
         
        ! Fractional standard deviation for rain condensate
         ZFRACSDR = (1.1_JPRB - 0.8_JPRB*ZA(JL,JK))*ZPHIR

        ! Unique value when essentially no cloud edges
         IF (ZA(JL,JK) > 0.95_JPRB) ZFRACSDR = 0.3_JPRB*ZPHIR

         ZEACC = &
     &    (1._JPRB+ZFRACSDC**2._JPRB)**((RCL_KKBAC-1._JPRB)*RCL_KKBAC/2._JPRB) &
     &   *(1._JPRB+ZFRACSDR**2._JPRB)**((RCL_KKBAC-1._JPRB)*RCL_KKBAC/2._JPRB) &
     &   * EXP(ZCLDRAINCORR*RCL_KKBAC*RCL_KKBAC* &
     &   SQRT(LOG(1._JPRB+ZFRACSDC**2._JPRB)*LOG(1._JPRB+ZFRACSDR**2._JPRB)))

       ELSE

        ! Simple constant multiplier to take account of inhomogeneity
        ZEAUT = RCL_INHOMOGAUT
        ZEACC = RCL_INHOMOGACC

       ENDIF

       IF (LLPERT_CLOUDINHOM) THEN
         ZEAUT = ZEAUT*EXP(PN1CLOUDINHOM%MU(1)+PN1CLOUDINHOM%XMAG(1)*PGP2DSPP(JL, IPCLOUDINHOMAUT))
         ZEACC = ZEACC*EXP(PN1CLOUDINHOM%MU(2)+PN1CLOUDINHOM%XMAG(2)*PGP2DSPP(JL, IPCLOUDINHOMACC))
       ENDIF

      !-------------------------------------------------------------------------
      ! Calculate autoconversion of cloud liquid droplets to rain
      ! Calculate accretion of cloud liquid droplets by rain
      !-------------------------------------------------------------------------

      !----------------------
      ! Explicit formulation
      !----------------------
      IF (IWARMRAIN == 2) THEN       
        ! Explicit formulation

        ! Autoconversion
        ZRAINAUT(JL)  = ZEAUT*ZA(JL,JK)*PTSPHY* &
     &                  RCL_KKAAU * ZLIQCLD(JL)**RCL_KKBAUQ * ZCONST**RCL_KKBAUN

        ZRAINAUT(JL) = MIN(ZRAINAUT(JL),ZQXFG(JL,NCLDQL))
        IF (ZRAINAUT(JL) < ZEPSEC) ZRAINAUT(JL) = 0.0_JPRB


        ! Accretion
        ZRAINACC(JL) = ZEACC*ZA(JL,JK)*PTSPHY* &
     &                 RCL_KKAAC * (ZLIQCLD(JL)*ZRAINCLD(JL))**RCL_KKBAC

        ZRAINACC(JL) = MIN(ZRAINACC(JL),ZQXFG(JL,NCLDQL))
        IF (ZRAINACC(JL) < ZEPSEC) ZRAINACC(JL) = 0.0_JPRB

        ! rain freezes instantly if T<0C
        IF (IP_SNOW_ACCRETES_RAIN == 2) THEN
 
          ! If temperature < 0, then autoconversion produces snow rather than rain
          ! Explicit
          IF(ZTP1(JL,JK) <= RTT) THEN
            ZSOLQA(JL,NCLDQS,NCLDQL)=ZSOLQA(JL,NCLDQS,NCLDQL)+ZRAINAUT(JL)
            ZSOLQA(JL,NCLDQS,NCLDQL)=ZSOLQA(JL,NCLDQS,NCLDQL)+ZRAINACC(JL)
            ZSOLQA(JL,NCLDQL,NCLDQS)=ZSOLQA(JL,NCLDQL,NCLDQS)-ZRAINAUT(JL)
            ZSOLQA(JL,NCLDQL,NCLDQS)=ZSOLQA(JL,NCLDQL,NCLDQS)-ZRAINACC(JL)
            ! Store cloud budget diagnostics
            ZBUDL(JL,12) = -ZRAINAUT(JL)*ZQTMST
            ZBUDL(JL,13) = -ZRAINACC(JL)*ZQTMST
          ELSE
            ZSOLQA(JL,NCLDQR,NCLDQL)=ZSOLQA(JL,NCLDQR,NCLDQL)+ZRAINAUT(JL)
            ZSOLQA(JL,NCLDQR,NCLDQL)=ZSOLQA(JL,NCLDQR,NCLDQL)+ZRAINACC(JL)
            ZSOLQA(JL,NCLDQL,NCLDQR)=ZSOLQA(JL,NCLDQL,NCLDQR)-ZRAINAUT(JL)
            ZSOLQA(JL,NCLDQL,NCLDQR)=ZSOLQA(JL,NCLDQL,NCLDQR)-ZRAINACC(JL)
            ! Store cloud budget diagnostics
            ZBUDL(JL,14) = -ZRAINAUT(JL)*ZQTMST
            ZBUDL(JL,15) = -ZRAINACC(JL)*ZQTMST
          ENDIF

        ELSE ! Don't freeze supercooled rain here

            ZSOLQA(JL,NCLDQR,NCLDQL)=ZSOLQA(JL,NCLDQR,NCLDQL)+ZRAINAUT(JL)
            ZSOLQA(JL,NCLDQR,NCLDQL)=ZSOLQA(JL,NCLDQR,NCLDQL)+ZRAINACC(JL)
            ZSOLQA(JL,NCLDQL,NCLDQR)=ZSOLQA(JL,NCLDQL,NCLDQR)-ZRAINAUT(JL)
            ZSOLQA(JL,NCLDQL,NCLDQR)=ZSOLQA(JL,NCLDQL,NCLDQR)-ZRAINACC(JL)
            ! Store cloud budget diagnostics
            ZBUDL(JL,14) = -ZRAINAUT(JL)*ZQTMST
            ZBUDL(JL,15) = -ZRAINACC(JL)*ZQTMST

        ENDIF

      !----------------------
      ! Implicit formulation
      !----------------------
      ELSEIF (IWARMRAIN == 3) THEN
        
        ! (zqxfg taken out and zliqcld=zqxfg/za, so multiply by 1/za
        !  and za cancels out) 
        ZRAINAUT(JL) = ZEAUT*PTSPHY* &
     &                 RCL_KKAAU * ZLIQCLD(JL)**(RCL_KKBAUQ-1.0_JPRB) &
     &                 * ZCONST**RCL_KKBAUN

        ! rain freezes instantly if T<0C
        IF (IP_SNOW_ACCRETES_RAIN == 2) THEN

          ! If temperature < 0, then autoconversion produces snow rather than rain
          ! Implicit
          IF(ZTP1(JL,JK) <= RTT) THEN
            ZSOLQB(JL,NCLDQS,NCLDQL)=ZSOLQB(JL,NCLDQS,NCLDQL)+ZRAINAUT(JL)
          ELSE
            ZSOLQB(JL,NCLDQR,NCLDQL)=ZSOLQB(JL,NCLDQR,NCLDQL)+ZRAINAUT(JL)
          ENDIF

        ELSE ! Don't freeze supercooled rain here

          ZSOLQB(JL,NCLDQR,NCLDQL)=ZSOLQB(JL,NCLDQR,NCLDQL)+ZRAINAUT(JL)

        ENDIF

        IF (ZRAINCLD(JL) > ZEPSEC) THEN

          IF (IRAINACC == 1) THEN

            ! Khairoutdinov and Kogan
            ! (zqxfg taken out and zliqcld=zqxfg/za, so multiply by 1/za
            !  and za cancels out) 
            ZRAINACC(JL) = ZEACC*PTSPHY* &
       &                   RCL_KKAAC * ZLIQCLD(JL)**(RCL_KKBAC-1.0_JPRB) &
       &                   * ZRAINCLD(JL)**RCL_KKBAC
          ELSE

            ! Sweep out
            ! (zqxfg taken out and zliqcld=zqxfg/za, so multiply by 1/za
            !  and za cancels out) - implicit in lwc
          
            ! Fallspeed air density correction 
            ZFALLCORR = (RDENSREF/ZRHO(JL))**0.4

            ! Slope of particle size distribution
            ZLAMBDA = (RCL_LAM1R/(ZRHO(JL)*ZRAINCLD(JL)))**RCL_LAM2R 
                    
            ! Calculate accretion term
            ! Factor of liq water taken out because implicit
            ZRAINACC(JL) = ZEACC*PTSPHY*RCL_EFF_RACW &
                         & *RCL_CONST9R*ZFALLCORR/(ZLAMBDA**RCL_CONST10R)

            ! Limit rain accretion term - needed?
            ZRAINACC(JL)=MIN(ZRAINACC(JL),1.0_JPRB)
          
          ENDIF
          
        ENDIF

        ! rain freezes instantly if T<0C
        IF (IP_SNOW_ACCRETES_RAIN == 2) THEN

          ! If temperature < 0, then autoconversion produces snow rather than rain
          ! Implicit
          IF(ZTP1(JL,JK) <= RTT) THEN
            ZSOLQB(JL,NCLDQS,NCLDQL)=ZSOLQB(JL,NCLDQS,NCLDQL)+ZRAINACC(JL)
            ZBUDL(JL,13) = -ZRAINACC(JL)*ZQTMST
          ELSE
            ZSOLQB(JL,NCLDQR,NCLDQL)=ZSOLQB(JL,NCLDQR,NCLDQL)+ZRAINACC(JL)
            ZBUDL(JL,15) = -ZRAINACC(JL)*ZQTMST
          ENDIF
 
        ELSE ! Don't freeze supercooled rain here

          ZSOLQB(JL,NCLDQR,NCLDQL)=ZSOLQB(JL,NCLDQR,NCLDQL)+ZRAINACC(JL)
          ZBUDL(JL,15) = -ZRAINACC(JL)*ZQTMST

        ENDIF

      ENDIF ! on IWARMRAIN = 2 or 3
    
    ENDIF ! on IWARMRAIN

   ENDIF ! on ZLIQCLD > ZEPSEC
  ENDDO

  !======================================================================
  !
  !
  ! 4.4  EVAPORATION OF RAIN
  !
  !
  !======================================================================
  ! Rain -> Vapour
  !----------------------------------------------------------------------

  DO JL=KIDIA,KFDIA

    !-----------------------------------------------------------------------
    ! Calculate relative humidity limit for rain evaporation 
    ! to avoid cloud formation and saturation of the grid box
    !-----------------------------------------------------------------------
    ! Limit RH for rain evaporation dependent on precipitation fraction 
    ZZRH=RPRECRHMAX+(1.0_JPRB-RPRECRHMAX)*ZCOVPMAX(JL)/MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))
    ZZRH=MIN(MAX(ZZRH,RPRECRHMAX),1.0_JPRB)

    ! Critical relative humidity
    IF (LLPERT_RAMID) THEN !Apply SPP perturbations
      ZXRAMID= RAMID*EXP(PN1RAMID%MU(1)+PN1RAMID%XMAG(1)*PGP2DSPP(JL, IPRAMID))
      IF (YDSPP_CONFIG%LRAMIDLIMIT1) THEN
        ZXRAMID= MIN(1.0_JPRB, ZXRAMID)
      ENDIF
    ELSE
      ZXRAMID= RAMID ! (unperturbed)
    ENDIF
    ZRHC=ZXRAMID
    ZSIGK=PAP(JL,JK)/PAPH(JL,KLEV+1)
    ! Increase RHcrit to 1.0 towards the surface (eta>0.8)
    IF(ZSIGK > 0.8_JPRB) THEN
      ZRHC=ZXRAMID+(1.0_JPRB-ZXRAMID)*((ZSIGK-0.8_JPRB)/0.2_JPRB)**2
    ENDIF
    
    ! Limit evaporation to RHcrit threshold
    ZZRH = MIN(ZRHC,ZZRH)
    ! Limit further to not go above 90%
    ZZRH = MIN(0.9_JPRB,ZZRH)

  
    ZQE=MAX(0.0_JPRB,MIN(ZQX(JL,JK,NCLDQV),ZQSLIQ(JL,JK)))

    ! If there is precipitation in clear sky, and there is rain present 
    ! and the humidity is less than the threshold defined above then... 
    LLO1=ZCOVPCLR(JL)>ZEPSEC .AND. &
       & ZRAINCLDM1(JL)>ZEPSEC .AND. & 
       & ZQE<ZZRH*ZQSLIQ(JL,JK)

    IF(LLO1) THEN

      !-------------------------------------------
      ! Evaporation
      !-------------------------------------------
      ! Calculate local precipitation (kg/kg)
      ZPRECLR = ZRAINCLDM1(JL)

      ! Fallspeed air density correction 
      ZFALLCORR = (RDENSREF/ZRHO(JL))**0.4

      ! Saturation vapour pressure with respect to liquid phase
      ZESATLIQ = RV/RD*FOEELIQ(ZTP1(JL,JK))

      ! Slope of particle size distribution
      ZLAMBDA = (RCL_LAM1R/(ZRHO(JL)*ZPRECLR))**RCL_LAM2R ! ZPRECLR=kg/kg

      ZEVAP_DENOM = RCL_CDENOM1*ZESATLIQ - RCL_CDENOM2*ZTP1(JL,JK)*ZESATLIQ &
              & + RCL_CDENOM3*ZTP1(JL,JK)**3._JPRB*PAP(JL,JK)

      ! Temperature dependent conductivity
      ZCORR2= (ZTP1(JL,JK)/273._JPRB)**1.5_JPRB*393._JPRB/(ZTP1(JL,JK)+120._JPRB)
      ZKA = RCL_KA273*ZCORR2

      ZSUBSAT = MAX(ZZRH*ZQSLIQ(JL,JK)-ZQE,0.0_JPRB)

      ZBETA = (0.15_JPRB/ZQSLIQ(JL,JK))*ZTP1(JL,JK)**2._JPRB*ZESATLIQ* &
     & RCL_CONST1R*(ZCORR2/ZEVAP_DENOM)*(0.78_JPRB/(ZLAMBDA**RCL_CONST4R)+ &
     & RCL_CONST2R*(ZRHO(JL)*ZFALLCORR)**0.5_JPRB/ &
     & (ZCORR2**0.5_JPRB*ZLAMBDA**RCL_CONST3R))
     
      ZDENOM  = 1.0_JPRB+ZBETA*PTSPHY*ZCORQSLIQ(JL)
      ZDPEVAP = ZCOVPCLR(JL)*ZBETA*PTSPHY*ZSUBSAT/ZDENOM

      !Apply SPP perturbations
      IF (LLPERT_RAINEVAP) THEN
        ZDPEVAP = ZDPEVAP*EXP(PN1RAINEVAP%MU(1)+PN1RAINEVAP%XMAG(1)*PGP2DSPP(JL, IPRAINEVAP))
      ENDIF

      !---------------------------------------------------------
      ! Add evaporation term to explicit sink.
      ! this has to be explicit since if treated in the implicit
      ! term evaporation can not reduce rain to zero and model
      ! produces small amounts of rainfall everywhere. 
      !---------------------------------------------------------
      
      ! Limit rain evaporation
      ZEVAP = MIN(ZDPEVAP,ZQXFG(JL,NCLDQR))

      ZSOLQA(JL,NCLDQV,NCLDQR) = ZSOLQA(JL,NCLDQV,NCLDQR)+ZEVAP
      ZSOLQA(JL,NCLDQR,NCLDQV) = ZSOLQA(JL,NCLDQR,NCLDQV)-ZEVAP

      ZBUDL(JL,20) = -ZEVAP*ZQTMST

      !-------------------------------------------------------------
      ! Reduce the total precip coverage proportional to evaporation
      ! to mimic the previous scheme which had a diagnostic
      ! 2-flux treatment, abandoned due to the new prognostic precip
      !-------------------------------------------------------------
      ! Comment out reduction of precipitation fraction with evaporation 
      !ZCOVPTOT(JL) = MAX(RCOVPMIN,ZCOVPTOT(JL)-MAX(0.0_JPRB, &
      ! &            (ZCOVPTOT(JL)-ZA(JL,JK))*ZEVAP/ZQXFG(JL,NCLDQR)))

      ! Update fg field 
      ZQXFG(JL,NCLDQR) = ZQXFG(JL,NCLDQR)-ZEVAP
    
    ENDIF
  ENDDO


  !======================================================================
  !
  !
  ! 4.5 GROWTH OF ICE BY VAPOUR DEPOSITION 
  !
  !
  !======================================================================
  ! does not use the ice nuclei number from cloudaer.F90
  ! but rather a simple Meyers et al. 1992 form based on the 
  ! supersaturation and assuming clouds are saturated with 
  ! respect to liquid water (well mixed), (or Koop adjustment)
  ! Growth considered as sink of liquid water if present
  !----------------------------------------------------------------------

  !--------------------------------------------------------
  !-
  !- Ice deposition following Rotstayn et al. (2001)
  !-  (monodisperse ice particle size distribution)
  !-
  !--------------------------------------------------------
  IF (IDEPICE == 1) THEN

!DIR$ IVDEP  
  DO JL=KIDIA,KFDIA

    !--------------------------------------------------------------
    ! Calculate distance from cloud top 
    ! defined by cloudy layer below a layer with cloud frac <0.01
    ! ZDZ = ZDP(JL)/(ZRHO(JL)*RG)
    !--------------------------------------------------------------
      
    IF (ZA(JL,JK-1) < RCLDTOPCF .AND. ZA(JL,JK) >= RCLDTOPCF) THEN
      ZCLDTOPDIST(JL) = 0.0_JPRB
    ELSE
      ZCLDTOPDIST(JL) = ZCLDTOPDIST(JL) + ZDP(JL)/(ZRHO(JL)*RG)
    ENDIF

    !--------------------------------------------------------------
    ! Set subgrid overlap fraction of supercooled liquid and ice
    ! Reduce in shallow convection because assume SLW in active 
    ! updraught is less overlapped with ice in less active part
    !--------------------------------------------------------------
    ZOVERLAP_LIQICE = RCL_OVERLAPLIQICE
    
    !IF (KTYPE(JL) > 0 .AND. PLUDE(JL,JK) > ZEPSEC) THEN
    !  ZOVERLAP_LIQICE = 0.1_JPRB
    !ENDIF

    !--------------------------------------------------------------
    ! only treat depositional growth if liquid present. due to fact 
    ! that can not model ice growth from vapour without additional 
    ! in-cloud water vapour variable
    !--------------------------------------------------------------
    ZSUPERSATICE = (ZQSLIQ(JL,JK)-ZQSICE(JL,JK))/ZQSICE(JL,JK)

    IF (ZTP1(JL,JK)<(RTT-5._JPRB) .AND. ZQXFG(JL,NCLDQL)>RLMIN) THEN  ! T<273K

      ZVPICE=FOEEICE(ZTP1(JL,JK))*RV/RD
      ZVPLIQ=ZVPICE*FOKOOP(ZTP1(JL,JK))
      ZICENUCLEI(JL)=1000.0_JPRB*EXP(12.96_JPRB*ZSUPERSATICE-0.639_JPRB)

      !------------------------------------------------
      !   2.4e-2 is conductivity of air
      !   8.8 = 700**1/3 = density of ice to the third
      !------------------------------------------------
      ZADD=RLSTT*(RLSTT/(RV*ZTP1(JL,JK))-1.0_JPRB)/(2.4E-2_JPRB*ZTP1(JL,JK))
      ZBDD=RV*ZTP1(JL,JK)*PAP(JL,JK)/(2.21_JPRB*ZVPICE)
      ZCVDS=7.8_JPRB*(ZICENUCLEI(JL)/ZRHO(JL))**0.666_JPRB*ZSUPERSATICE / &
         & (8.87_JPRB*(ZADD+ZBDD))

      !-----------------------------------------------------
      ! RICEINIT=1.E-12_JPRB is initial mass of ice particle
      !-----------------------------------------------------
      ZICE0=MAX(ZICECLD(JL), ZICENUCLEI(JL)*RICEINIT/ZRHO(JL))

      !------------------
      ! new value of ice
      !------------------
      ZINEW=(0.666_JPRB*ZCVDS*PTSPHY+ZICE0**0.666_JPRB)**1.5_JPRB

      !---------------------------
      ! grid-mean deposition rate:
      !--------------------------- 
      ZDEPOS=MAX(ZOVERLAP_LIQICE*ZA(JL,JK)*(ZINEW-ZICE0),0.0_JPRB)

      !--------------------------------------------------------------------
      ! Limit deposition to liquid water amount
      ! If liquid is all frozen, ice would use up reservoir of water 
      ! vapour in excess of ice saturation mixing ratio - However this 
      ! can not be represented without a in-cloud humidity variable. Using 
      ! the grid-mean humidity would imply a large artificial horizontal 
      ! flux from the clear sky to the cloudy area. We thus rely on the 
      ! supersaturation check to clean up any remaining supersaturation
      !--------------------------------------------------------------------
      ZDEPOS=MIN(ZDEPOS,ZQXFG(JL,NCLDQL)) ! limit to liquid water amount
      
      !--------------------------------------------------------------------
      ! At top of cloud, reduce deposition rate near cloud top to account for
      ! small scale turbulent processes, limited ice nucleation and ice fallout 
      !--------------------------------------------------------------------
      ! Include dependence on ice nuclei concentration
      ! to increase deposition rate with decreasing temperatures 
      ZINFACTOR = MIN(ZICENUCLEI(JL)/15000._JPRB, 1.0_JPRB)
      ZDEPOS = ZDEPOS*MIN(ZINFACTOR + (1.0_JPRB-ZINFACTOR)* &
                  & (RDEPLIQREFRATE+ZCLDTOPDIST(JL)/RDEPLIQREFDEPTH),1.0_JPRB)

      !--------------
      ! add to matrix 
      !--------------
      ZSOLQA(JL,NCLDQI,NCLDQL)=ZSOLQA(JL,NCLDQI,NCLDQL)+ZDEPOS
      ZSOLQA(JL,NCLDQL,NCLDQI)=ZSOLQA(JL,NCLDQL,NCLDQI)-ZDEPOS
      ZQXFG(JL,NCLDQI)=ZQXFG(JL,NCLDQI)+ZDEPOS
      ZQXFG(JL,NCLDQL)=ZQXFG(JL,NCLDQL)-ZDEPOS
      ! Store cloud budget diagnostics if required
      ZBUDL(JL,11) = -ZDEPOS*ZQTMST
      ZBUDI(JL,11) = ZDEPOS*ZQTMST

    ENDIF
  ENDDO

  !--------------------------------------------------------
  !-
  !- Ice deposition assuming ice PSD
  !-
  !--------------------------------------------------------
  ELSEIF (IDEPICE == 2) THEN

    DO JL=KIDIA,KFDIA

      !--------------------------------------------------------------
      ! Calculate distance from cloud top 
      ! defined by cloudy layer below a layer with cloud frac <0.01
      ! ZDZ = ZDP(JL)/(ZRHO(JL)*RG)
      !--------------------------------------------------------------

      IF (ZA(JL,JK-1) < RCLDTOPCF .AND. ZA(JL,JK) >= RCLDTOPCF) THEN
        ZCLDTOPDIST(JL) = 0.0_JPRB
      ELSE
        ZCLDTOPDIST(JL) = ZCLDTOPDIST(JL) + ZDP(JL)/(ZRHO(JL)*RG)
      ENDIF

      !--------------------------------------------------------------
      ! Set subgrid overlap fraction of supercooled liquid and ice
      ! Reduce in shallow convection because assume SLW in active 
      ! updraught is less overlapped with ice in less active part
      !--------------------------------------------------------------
      ZOVERLAP_LIQICE = RCL_OVERLAPLIQICE
    
      !IF (KTYPE(JL) > 0 .AND. PLUDE(JL,JK) > ZEPSEC) THEN
      !  ZOVERLAP_LIQICE = 0.1_JPRB
      !ENDIF

      !--------------------------------------------------------------
      ! only treat depositional growth if liquid present. due to fact 
      ! that can not model ice growth from vapour without additional 
      ! in-cloud water vapour variable
      !--------------------------------------------------------------
      IF (ZTP1(JL,JK)<(RTT-5._JPRB) .AND. ZQXFG(JL,NCLDQL)>RLMIN) THEN  ! T<273K
      
        ZVPICE = FOEEICE(ZTP1(JL,JK))*RV/RD
        ZVPLIQ = ZVPICE*FOKOOP(ZTP1(JL,JK))
        ZICENUCLEI(JL)=1000.0_JPRB*EXP(12.96_JPRB*(ZVPLIQ-ZVPICE)/ZVPLIQ-0.639_JPRB)

        !-----------------------------------------------------
        ! RICEINIT=1.E-12_JPRB is initial mass of ice particle
        !-----------------------------------------------------
        ZICE0=MAX(ZICECLD(JL), ZICENUCLEI(JL)*RICEINIT/ZRHO(JL))
        
        ! Particle size distribution
        ZTCG    = 1.0_JPRB
        ZFACX1I = 1.0_JPRB

        ZAPLUSB   = RCL_APB1*ZVPICE-RCL_APB2*ZVPICE*ZTP1(JL,JK)+ &
       &             PAP(JL,JK)*RCL_APB3*ZTP1(JL,JK)**3._JPRB
        ZCORRFAC  = (1.0_JPRB/ZRHO(JL))**0.5_JPRB
        ZCORRFAC2 = ((ZTP1(JL,JK)/273.0_JPRB)**1.5_JPRB) &
       &             *(393.0_JPRB/(ZTP1(JL,JK)+120.0_JPRB))

        ZPR02  = ZRHO(JL)*ZICE0*RCL_CONST1I/(ZTCG*ZFACX1I)

        ZTERM1 = (ZVPLIQ-ZVPICE)*ZTP1(JL,JK)**2.0_JPRB*ZVPICE*ZCORRFAC2*ZTCG* &
       &          RCL_CONST2I*ZFACX1I/(ZRHO(JL)*ZAPLUSB*ZVPICE)
        ZTERM2 = 0.65_JPRB*RCL_CONST6I*ZPR02**RCL_CONST4I+RCL_CONST3I &
       &          *ZCORRFAC**0.5_JPRB*ZRHO(JL)**0.5_JPRB &
       &          *ZPR02**RCL_CONST5I/ZCORRFAC2**0.5_JPRB

        ZDEPOS = MAX(ZOVERLAP_LIQICE*ZA(JL,JK)*ZTERM1*ZTERM2*PTSPHY,0.0_JPRB)

        !--------------------------------------------------------------------
        ! Limit deposition to liquid water amount
        ! If liquid is all frozen, ice would use up reservoir of water 
        ! vapour in excess of ice saturation mixing ratio - However this 
        ! can not be represented without a in-cloud humidity variable. Using 
        ! the grid-mean humidity would imply a large artificial horizontal 
        ! flux from the clear sky to the cloudy area. We thus rely on the 
        ! supersaturation check to clean up any remaining supersaturation
        !--------------------------------------------------------------------
        ZDEPOS=MIN(ZDEPOS,ZQXFG(JL,NCLDQL)) ! limit to liquid water amount

        !--------------------------------------------------------------------
        ! At top of cloud, reduce deposition rate near cloud top to account for
        ! small scale turbulent processes, limited ice nucleation and ice fallout 
        !--------------------------------------------------------------------
        ! Change to include dependence on ice nuclei concentration
        ! to increase deposition rate with decreasing temperatures 
        ZINFACTOR = MIN(ZICENUCLEI(JL)/15000._JPRB, 1.0_JPRB)
        ZDEPOS = ZDEPOS*MIN(ZINFACTOR + (1.0_JPRB-ZINFACTOR)* &
                    & (RDEPLIQREFRATE+ZCLDTOPDIST(JL)/RDEPLIQREFDEPTH),1.0_JPRB)

        !--------------
        ! add to matrix 
        !--------------
        ZSOLQA(JL,NCLDQI,NCLDQL) = ZSOLQA(JL,NCLDQI,NCLDQL)+ZDEPOS
        ZSOLQA(JL,NCLDQL,NCLDQI) = ZSOLQA(JL,NCLDQL,NCLDQI)-ZDEPOS
        ZQXFG(JL,NCLDQI) = ZQXFG(JL,NCLDQI)+ZDEPOS
        ZQXFG(JL,NCLDQL) = ZQXFG(JL,NCLDQL)-ZDEPOS
        ! Store cloud budget diagnostics if required
        ZBUDL(JL,11) = -ZDEPOS*ZQTMST
        ZBUDI(JL,11) = ZDEPOS*ZQTMST

      ENDIF
    ENDDO

  ENDIF ! on IDEPICE

  !======================================================================
  !
  ! 4.6 PIEVAP: ICE SUBLIMATION ASSUMING ICE PSD
  !
  !======================================================================

  ! Initialise process rate
  ZPIEVAP(KIDIA:KFDIA) = 0.0_JPRB

  IF (ISUBLICE == 1) THEN

    DO JL=KIDIA,KFDIA

      !-----------------------------------------------------------------------
      ! Environmental humidity for ice evaporation 
      !-----------------------------------------------------------------------
!      ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSICE(JL,JK))/ &
!      & MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))  
      ZQE=ZQX(JL,JK,NCLDQV)  

      !---------------------------------------------
      ! humidity in moistest ZCOVPCLR part of domain
      !---------------------------------------------
      IF (ZFALLSRCE(JL,NCLDQI) > 0.0_JPRB  .AND. ZCOVPTOT(JL) > ZEPSEC) THEN 

        ! Calculate local in-cloud ice that has fallen into this gridbox (kg/kg)
        ZICE0 = ZFALLSRCE(JL,NCLDQI)/ZCOVPTOT(JL)

        ! Calculate subsaturation
        ZSUPERSATICE = (ZQSICE(JL,JK)-ZQE)/ZQSICE(JL,JK)

        ! Calculate saturation vapour pressure wrt ice
        ZVPICE = FOEEICE(ZTP1(JL,JK))*RV/RD
        
        ! Particle size distribution
        ZTCG    = 1.0_JPRB
        ZFACX1I = 1.0_JPRB

        ZAPLUSB   = RCL_APB1*ZVPICE-RCL_APB2*ZVPICE*ZTP1(JL,JK)+ &
       &             PAP(JL,JK)*RCL_APB3*ZTP1(JL,JK)**3._JPRB
        ZCORRFAC  = (1.0_JPRB/ZRHO(JL))**0.5_JPRB
        ZCORRFAC2 = ((ZTP1(JL,JK)/273.0_JPRB)**1.5_JPRB) &
       &             *(393.0_JPRB/(ZTP1(JL,JK)+120.0_JPRB))

        ZPR02  = ZRHO(JL)*ZICE0*RCL_CONST1I/(ZTCG*ZFACX1I)

        ZTERM1 = ZSUPERSATICE*ZTP1(JL,JK)**2.0_JPRB*ZVPICE*ZCORRFAC2*ZTCG* &
       &          RCL_CONST2I*ZFACX1I/(ZRHO(JL)*ZAPLUSB)
        ZTERM2 = 0.65_JPRB*RCL_CONST6I*ZPR02**RCL_CONST4I+RCL_CONST3I &
       &          *ZCORRFAC**0.5_JPRB*ZRHO(JL)**0.5_JPRB &
       &          *ZPR02**RCL_CONST5I/ZCORRFAC2**0.5_JPRB

        ZPIEVAP(JL) = ZCOVPCLR(JL)*ZTERM1*ZTERM2*PTSPHY
        ! Proportion that is falling in to clear sky is ZCOVPCLR(JL)/ZCOVPTOT(JL) 
        ! and*ZCOVPOT(JL) to go from in-cloud back to gridbox mean gives *ZCOVPCLR(JL)

        !--------------------------------------------------------------------
        ! Limit sublimation to ice water amount or deposition to available supersaturation
        !--------------------------------------------------------------------
        IF (ZPIEVAP(JL) > 0.0_JPRB) THEN
          ! Sublimation
          ZPIEVAP(JL) = MIN(ZPIEVAP(JL),ZEVAPLIMMIX(JL))
          ZPIEVAP(JL) = MIN(ZPIEVAP(JL),ZICE0)
        ELSE
          ! Deposition
!          ZPIEVAP(JL) = MAX(ZPIEVAP(JL),-(ZQE-ZQSMIX(JL,JK))) !/ZCORQSMIX(JL,JK))        
          ZPIEVAP(JL) = 0.0_JPRB
        ENDIF

        !--------------
        ! add to matrix 
        !--------------
        ZSOLQA(JL,NCLDQV,NCLDQI) = ZSOLQA(JL,NCLDQV,NCLDQI)+ZPIEVAP(JL)
        ZSOLQA(JL,NCLDQI,NCLDQV) = ZSOLQA(JL,NCLDQI,NCLDQV)-ZPIEVAP(JL)
        ZQXFG(JL,NCLDQI) = ZQXFG(JL,NCLDQI)-ZPIEVAP(JL)

        ! Decrease cloud amount using RKOOPTAU timescale
        !ZFACI = PTSPHY/RKOOPTAU
        !ZSOLAC(JL) = ZSOLAC(JL)-ZA(JL,JK)*ZFACI

      ENDIF
    ENDDO

  ENDIF ! on ISUBLICE

  !======================================================================
  !
  !
  ! 4.7 PSDEP: DEPOSITION ONTO SNOW
  !
  !
  !======================================================================
  ZPSDEP(KIDIA:KFDIA) = 0.0_JPRB
   
  IF (IDEPSNOW == 1) THEN

    DO JL=KIDIA,KFDIA

      ZSUPERSATICE = MAX((ZQX(JL,JK,NCLDQV)-ZQSMIX(JL,JK))/ZQSMIX(JL,JK),0.0_JPRB)

      IF (ZTP1(JL,JK)<RTT .AND. ZSNOWCLD(JL) > ZEPSEC .AND. ZSUPERSATICE > 0.0_JPRB) THEN  ! T<273K
      
        ZVPICE = FOEEICE(ZTP1(JL,JK))*RV/RD
       
        ! Particle size distribution
        ZTCG    = 1.0_JPRB
        ZFACX1S = 1.0_JPRB

        ZAPLUSB   = RCL_APB1*ZVPICE-RCL_APB2*ZVPICE*ZTP1(JL,JK)+ &
       &             PAP(JL,JK)*RCL_APB3*ZTP1(JL,JK)**3
        ZCORRFAC  = (1.0_JPRB/ZRHO(JL))**0.5_JPRB
        ZCORRFAC2 = ((ZTP1(JL,JK)/273.0_JPRB)**1.5_JPRB) &
       &             *(393.0_JPRB/(ZTP1(JL,JK)+120.0_JPRB))

        ZPR02  = ZRHO(JL)*ZSNOWCLD(JL)*RCL_CONST1S/(ZTCG*ZFACX1S)

        ZTERM1 = ZSUPERSATICE*ZTP1(JL,JK)**2*ZVPICE*ZCORRFAC2*ZTCG* &
       &          RCL_CONST2S*ZFACX1S/(ZRHO(JL)*ZAPLUSB)
        ZTERM2 = 0.65_JPRB*RCL_CONST6S*ZPR02**RCL_CONST4S+RCL_CONST3S &
       &          *ZCORRFAC**0.5_JPRB*ZRHO(JL)**0.5_JPRB &
       &          *ZPR02**RCL_CONST5S/ZCORRFAC2**0.5_JPRB

        ZPSDEP(JL) = MAX(ZCOVPTOT(JL)*ZTERM1*ZTERM2*PTSPHY,0.0_JPRB)

        !--------------------------------------------------------------------
        ! Limit to available supersaturation
        !--------------------------------------------------------------------
        ZPSDEP(JL) = MIN(ZPSDEP(JL),ZQX(JL,JK,NCLDQV)-ZQSMIX(JL,JK))/ZCORQSMIX(JL)
        
        !--------------
        ! add to matrix 
        !--------------
        ZSOLQA(JL,NCLDQS,NCLDQV) = ZSOLQA(JL,NCLDQS,NCLDQV)+ZPSDEP(JL)
        ZSOLQA(JL,NCLDQV,NCLDQS) = ZSOLQA(JL,NCLDQV,NCLDQS)-ZPSDEP(JL)

      ENDIF

    ENDDO

  ENDIF ! on IDEPSNOW
  
  !----------------------------------
  ! revise in-cloud consensate amount
  !----------------------------------
  DO JL=KIDIA,KFDIA
    ZTMPA = 1.0_JPRB/MAX(ZA(JL,JK),0.01_JPRB)
    ZLIQCLD(JL) = ZQXFG(JL,NCLDQL)*ZTMPA
    ZICECLD(JL) = ZQXFG(JL,NCLDQI)*ZTMPA
  ENDDO


  !======================================================================
  !
  !
  ! 4.8  AUTOCONVERSION OF ICE TO SNOW
  !
  !
  !======================================================================
  ! Formulation follows Lin et al. (1983)
  ! Sink of ice, source of snow
  ! Implicit in ice water content 
  !----------------------------------------------------------------------
  
  DO JL=KIDIA,KFDIA
 
    IF(ZTP1(JL,JK) <= RTT) THEN
      IF (ZICECLD(JL)>ZEPSEC) THEN

        ZZCO=PTSPHY*RSNOWLIN1*EXP(RSNOWLIN2*(ZTP1(JL,JK)-RTT))

        IF (LAERICEAUTO) THEN
          ZLCRIT=PICRIT_AER(JL,JK)
          ! 0.3 = N**0.333 with N=0.027 
          ZZCO=ZZCO*(RNICE/PNICE(JL,JK))**0.333_JPRB
        ELSE
          IF (LLPERT_RLCRITSNOW) THEN !Apply SPP perturbations
            ZLCRIT= RLCRITSNOW*EXP(PN1RLCRITSNOW%MU(1)+PN1RLCRITSNOW%XMAG(1)*PGP2DSPP(JL, IPRLCRITSNOW))
          ELSE
            ZLCRIT=RLCRITSNOW    ! (unperturbed)
          ENDIF
        ENDIF

        ZTEMP = -1._JPRB*(ZICECLD(JL)/ZLCRIT)**2
        IF (ZTEMP < 50.0_JPRB ) THEN ! Security for exp
          ZSNOWAUT(JL)=ZZCO*(1.0_JPRB-EXP(ZTEMP))
        ELSE
          ZSNOWAUT(JL)=ZZCO
        ENDIF
 
        ZSOLQB(JL,NCLDQS,NCLDQI)=ZSOLQB(JL,NCLDQS,NCLDQI)+ZSNOWAUT(JL)

      ENDIF
    ENDIF 

  ENDDO

  !======================================================================
  !
  !
  ! 4.12  MELTING OF SNOW AND ICE
  !
  !
  !======================================================================
  ! With implicit solver this also has to treat snow or ice
  ! precipitating from the level above... i.e. local ice AND flux.
  ! in situ ice and snow: could arise from LS advection or warming
  ! falling ice and snow: arrives by precipitation process
  !----------------------------------------------------------------------
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA

    ZICETOT(JL)=ZQXFG(JL,NCLDQI)+ZQXFG(JL,NCLDQS)
    ZMELTMAX(JL) = 0.0_JPRB

    ! If there are frozen hydrometeors present and dry-bulb temperature > 0degC
    IF(ZICETOT(JL) > ZEPSEC .AND. ZTP1(JL,JK) > RTT) THEN

      ! Calculate subsaturation
      ZSUBSAT = MAX(ZQSICE(JL,JK)-ZQX(JL,JK,NCLDQV),0.0_JPRB)
      
      ! Calculate difference between dry-bulb (ZTP1) and the temperature 
      ! at which the wet-bulb=0degC (RTT-ZSUBSAT*....) using an approx.
      ! Melting only occurs if the wet-bulb temperature >0
      ! i.e. warming of ice particle due to melting > cooling 
      ! due to evaporation.
      ZTDMTW0 = ZTP1(JL,JK)-RTT-ZSUBSAT* &
                & (ZTW1+ZTW2*(PAP(JL,JK)-ZTW3)-ZTW4*(ZTP1(JL,JK)-ZTW5))
      ! Not implicit yet... 
      ! Ensure ZCONS1 is positive so that ZMELTMAX=0 if ZTDMTW0<0
      ZCONS1 = ABS(PTSPHY*(1.0_JPRB+0.5_JPRB*ZTDMTW0)/RTAUMEL)
      ZMELTMAX(JL) = MAX(ZTDMTW0*ZCONS1*ZRLDCP,0.0_JPRB)
    ENDIF
  ENDDO

  ! Loop over frozen hydrometeors (ice, snow)
  DO JM=1,NCLV
   IF (IPHASE(JM) == 2) THEN
    JN = IMELT(JM)
    DO JL=KIDIA,KFDIA
      IF(ZMELTMAX(JL)>ZEPSEC .AND. ZICETOT(JL)>ZEPSEC) THEN
        ! Apply melting in same proportion as frozen hydrometeor fractions 
        ZALFA = ZQXFG(JL,JM)/ZICETOT(JL)
        ZMELT = MIN(ZQXFG(JL,JM),ZALFA*ZMELTMAX(JL))
        ! needed in first guess
        ! This implies that zqpretot has to be recalculated below
        ! since is not conserved here if ice falls and liquid doesn't
        ZQXFG(JL,JM)     = ZQXFG(JL,JM)-ZMELT
        ZQXFG(JL,JN)     = ZQXFG(JL,JN)+ZMELT
        ZSOLQA(JL,JN,JM) = ZSOLQA(JL,JN,JM)+ZMELT
        ZSOLQA(JL,JM,JN) = ZSOLQA(JL,JM,JN)-ZMELT
        IF (JM==NCLDQI) ZBUDI(JL,15) = -ZMELT*ZQTMST
        IF (JM==NCLDQI) ZBUDL(JL,17) =  ZMELT*ZQTMST
        IF (JM==NCLDQS) ZBUDI(JL,16) = -ZMELT*ZQTMST
      ENDIF
    ENDDO
   ENDIF
  ENDDO

  
  !======================================================================
  !
  !
  ! 4.13  FREEZING OF RAIN
  !
  !
  !======================================================================
  ! Rain drop freezing rate based on Bigg(1953) and Wisner(1972)
  ! Rain -> Snow
  !----------------------------------------------------------------------
  
!DEC$ IVDEP
  DO JL=KIDIA,KFDIA 

    ! If rain present
    IF (ZQX(JL,JK,NCLDQR) > ZEPSEC) THEN

      IF (ZTP1(JL,JK) <= RTT .AND. ZTP1(JL,JK-1) > RTT) THEN
        ! Base of melting layer/top of refreezing layer so
        ! store rain/snow fraction for precip type diagnosis
        ! If mostly rain, then supercooled rain slow to freeze
        ! otherwise faster to freeze (snow or ice pellets)
        ZQPRETOT(JL) = MAX(ZQX(JL,JK,NCLDQS)+ZQX(JL,JK,NCLDQR),ZEPSEC)
        PRAINFRAC_TOPRFZ(JL) = ZQX(JL,JK,NCLDQR)/ZQPRETOT(JL)
        IF (PRAINFRAC_TOPRFZ(JL) > 0.8) THEN 
          LLRAINLIQ(JL) = .TRUE.
        ELSE
          LLRAINLIQ(JL) = .FALSE.
        ENDIF
      ENDIF
    
      ! If temperature less than zero
      IF (ZTP1(JL,JK) < RTT) THEN

        IF (LLRAINLIQ(JL)) THEN 

          ! Majority of raindrops completely melted
          ! Refreezing is by slow heterogeneous freezing
          
          ! Slope of rain particle size distribution
          ZLAMBDA = (RCL_LAM1R/(ZRHO(JL)*ZQX(JL,JK,NCLDQR)))**RCL_LAM2R

          ! Calculate freezing rate based on Bigg(1953) and Wisner(1972)
          ZTEMP = MIN(RCL_FZRAB*(ZTP1(JL,JK)-RTT), 50._JPRB) ! for EXP security
          ZFRZ  = PTSPHY * (RCL_CONST5R/ZRHO(JL)) * (EXP(ZTEMP)-1._JPRB) &
                  & * ZLAMBDA**RCL_CONST6R
          ZFRZMAX(JL) = MAX(ZFRZ,0.0_JPRB)

        ELSE

          ! Majority of raindrops only partially melted 
          ! Refreeze with a shorter timescale (reverse of melting...for now)
          
          ZCONS1 = ABS(PTSPHY*(1.0_JPRB+0.5_JPRB*(RTT-ZTP1(JL,JK)))/RTAUMEL)
          ZFRZMAX(JL) = MAX((RTT-ZTP1(JL,JK))*ZCONS1*ZRLDCP,0.0_JPRB)

        ENDIF

        IF(ZFRZMAX(JL)>ZEPSEC) THEN
          ZFRZ = MIN(ZQX(JL,JK,NCLDQR),ZFRZMAX(JL))
          ZSOLQA(JL,NCLDQS,NCLDQR) = ZSOLQA(JL,NCLDQS,NCLDQR)+ZFRZ
          ZSOLQA(JL,NCLDQR,NCLDQS) = ZSOLQA(JL,NCLDQR,NCLDQS)-ZFRZ
          ZBUDL(JL,18) = -ZFRZ*ZQTMST
        ENDIF
      ENDIF

    ENDIF

  ENDDO


  !======================================================================
  !
  !
  ! 4.14   FREEZING OF CLOUD LIQUID 
  !
  !
  !======================================================================
  ! All liquid cloud drops assumed to freeze instantaneously to ice crystals
  ! below the homogeneous freezing temperature (-38degC)
  ! Liquid -> Ice
  !----------------------------------------------------------------------
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA 
    ! not implicit yet... 
    ZFRZMAX(JL)=MAX((RTHOMO-ZTP1(JL,JK))*ZRLDCP,0.0_JPRB)
  ENDDO

  DO JL=KIDIA,KFDIA
    IF(ZFRZMAX(JL)>ZEPSEC .AND. ZQXFG(JL,NCLDQL)>ZEPSEC) THEN
      ZFRZ = MIN(ZQXFG(JL,NCLDQL),ZFRZMAX(JL))
      ZSOLQA(JL,NCLDQI,NCLDQL) = ZSOLQA(JL,NCLDQI,NCLDQL)+ZFRZ
      ZSOLQA(JL,NCLDQL,NCLDQI) = ZSOLQA(JL,NCLDQL,NCLDQI)-ZFRZ
      ZBUDL(JL,19) = -ZFRZ*ZQTMST
    ENDIF
  ENDDO

  !======================================================================
  !
  !
  ! 4.9 RIMING - COLLECTION OF CLOUD LIQUID DROPS BY SNOW AND ICE
  !
  !
  !======================================================================
  ! Only active if T<0degC and supercooled liquid water is present
  ! AND if not Sundquist autoconversion (as this includes riming)
  !----------------------------------------------------------------------
  
  IF (IWARMRAIN > 1) THEN

  DO JL=KIDIA,KFDIA
    
    IF(ZTP1(JL,JK) <= RTT .AND. ZLIQCLD(JL)>ZEPSEC) THEN

      ! Fallspeed air density correction 
      ZFALLCORR = (RDENSREF/ZRHO(JL))**0.4

      !------------------------------------------------------------------
      ! Riming of snow by cloud water - implicit in lwc
      !------------------------------------------------------------------
      IF (ZSNOWCLD(JL)>ZEPSEC .AND. ZCOVPTOT(JL)>0.01_JPRB) THEN

        ! Calculate riming term
        ! Factor of liq water taken out because implicit
        ZSNOWRIME(JL) = RCL_EFFRIME*ZCOVPTOT(JL)*PTSPHY*RCL_CONST7S*ZFALLCORR &
     &                  *(ZRHO(JL)*ZSNOWCLD(JL)*RCL_CONST1S)**RCL_CONST8S

        ! Limit snow riming term
        ZSNOWRIME(JL)=MIN(ZSNOWRIME(JL),1.0_JPRB)

        ZSOLQB(JL,NCLDQS,NCLDQL) = ZSOLQB(JL,NCLDQS,NCLDQL) + ZSNOWRIME(JL)

      ENDIF

      !------------------------------------------------------------------
      ! Riming of ice by cloud water - implicit in lwc
      ! NOT YET ACTIVE
      !------------------------------------------------------------------
!      IF (ZICECLD(JL)>ZEPSEC .AND. ZA(JL,JK)>0.01_JPRB) THEN
!
!        ! Calculate riming term
!        ! Factor of liq water taken out because implicit
!        ZSNOWRIME(JL) = ZA(JL,JK)*PTSPHY*RCL_CONST7S*ZFALLCORR &
!     &                  *(ZRHO(JL)*ZICECLD(JL)*RCL_CONST1S)**RCL_CONST8S
!
!        ! Limit ice riming term
!        ZSNOWRIME(JL)=MIN(ZSNOWRIME(JL),1.0_JPRB)
!
!        ZSOLQB(JL,NCLDQI,NCLDQL) = ZSOLQB(JL,NCLDQI,NCLDQL) + ZSNOWRIME(JL)
!
!      ENDIF
    ENDIF
  ENDDO
  
  !======================================================================
  !
  !
  ! Calculate slope of particle size distributions LAMBDA
  !
  !
  !======================================================================
  IF (IP_SNOW_ACCRETES_RAIN == 1 .OR. IP_RAIN_ACCRETES_SNOW == 1 .OR. &
  &   IP_ICE_ACCRETES_RAIN  == 1 .OR. IP_RAIN_ACCRETES_ICE  == 1) THEN

    DO JL=KIDIA,KFDIA

      ! Slope of rain particle size distribution
      IF(ZRAINCLD(JL) > ZEPSEC) THEN
        ZLAMR(JL) = (RCL_LAM1R/(ZRHO(JL)*ZRAINCLD(JL)))**RCL_LAM2R
        ZN0R(JL)  = RCL_X1R*ZLAMR(JL)**RCL_X2R
      ELSE
        ZLAMR(JL) = 0.0_JPRB
        ZN0R(JL)  = 0.0_JPRB 
      ENDIF    

      ! Slope of rain particle size distribution (Marshall-Palmer)
      !IF(ZRAINCLD(JL) > ZEPSEC) THEN
      !  ZLAMR_MP(JL) = (RCL_LAM1R_MP/(ZRHO(JL)*ZRAINCLD(JL)))**RCL_LAM2R_MP
      !  ZN0R_MP(JL)  = RCL_X1R_MP*ZLAMR(JL)**RCL_X2R_MP
      !ELSE
      !  ZLAMR_MP(JL) = 0.0_JPRB
      !  ZN0R_MP(JL)  = 0.0_JPRB 
      !ENDIF    

      ! Slope of snow particle size distribution
      IF(ZSNOWCLD(JL) > ZEPSEC) THEN
        ZLAMS(JL) = (RCL_LAM1S/(ZRHO(JL)*ZSNOWCLD(JL)))**RCL_LAM2S
        ZN0S(JL)  = RCL_X1S ! *ZLAMS(JL)**RCL_X2S ! comment as RCL_X2S=0.0
      ELSE
        ZLAMS(JL) = 0.0_JPRB
        ZN0S(JL)  = 0.0_JPRB 
      ENDIF    

      ! Slope of ice particle size distribution
      IF(ZICECLD(JL) > ZEPSEC) THEN
        ZLAMI(JL) = (RCL_LAM1S/(ZRHO(JL)*ZICECLD(JL)))**RCL_LAM2S
        ZN0I(JL)  = RCL_X1S ! *ZLAMS(JL)**RCL_X2S ! comment as RCL_X2S=0.0
      ELSE
        ZLAMI(JL) = 0.0_JPRB
        ZN0I(JL)  = 0.0_JPRB 
      ENDIF    

    ENDDO  
  
  ENDIF
  
  !======================================================================
  !
  !
  ! PSACR: SNOW ACCRETES RAIN -> SNOW
  !
  !
  !======================================================================
  ! Only active if T<0degC, rain and snow present and no warm layer above.
  ! Therefore only active if there is a mix of rain and snow particles
  ! Don't want this process when rain/snow mix is representing the water
  ! content of individual particles (in the melting layer, or for 
  ! supercooled rain formed from a T>0 warm layer, rather than a mix of rain 
  !----------------------------------------------------------------------
  ZPSACR(:)  = 0.0_JPRB
  ZEFF_PSACR = 1.0_JPRB ! Collection efficiency of snow accretes rain
  
  IF (IP_SNOW_ACCRETES_RAIN == 1) THEN

    DO JL=KIDIA,KFDIA
      IF ( ZTP1(JL,JK) <= RTT .AND. PRAINFRAC_TOPRFZ(JL) < ZEPSEC &
        & .AND. ZRAINCLD(JL) > ZEPSEC .AND. ZSNOWCLD(JL) > ZEPSEC ) THEN

        !------------------------------------------------------------------
        ! PSACR - accretion of rain by snow - explicit
        !------------------------------------------------------------------

        ! Calculate differential fallspeed
        ZFALLDIFF = MAX(ABS(ZFALLSPEED(JL,NCLDQR)-ZFALLSPEED(JL,NCLDQS)), &
                    & (ZFALLSPEED(JL,NCLDQR)+ZFALLSPEED(JL,NCLDQS))/8._JPRB)

        ZFAC = RPI*RPI*RDENSWAT*ZEFF_PSACR*ZRHO(JL)*ZN0R(JL)*ZN0S(JL)/(ZLAMS(JL)*ZLAMR(JL)**3)
        ZRAT = ZLAMR(JL)/ZLAMS(JL)
        ZPSACR(JL) = ZCOVPTOT(JL)*PTSPHY*ZFAC*ZFALLDIFF &
                   & *((0.5_JPRB*ZRAT+2._JPRB)*ZRAT+5._JPRB)/ZLAMR(JL)**3 
          
        ! Limit to amount of rain
        !ZPSACR(JL) = MIN(ZPSACR(JL),ZRAINCLD(JL))
          
        ZSOLQA(JL,NCLDQS,NCLDQR) = ZSOLQA(JL,NCLDQS,NCLDQR) + ZPSACR(JL)
        ZSOLQA(JL,NCLDQR,NCLDQS) = ZSOLQA(JL,NCLDQR,NCLDQS) - ZPSACR(JL)
       
      ENDIF
    ENDDO

  ENDIF ! on IP_SNOW_ACCRETES_RAIN=1
  
  !======================================================================
  !
  !
  ! PRACS: RAIN ACCRETES SNOW -> SNOW
  !
  !
  !======================================================================
  ! Only active if T<0degC, rain and snow present and no warm layer above.
  ! Therefore only active if there is a mix of rain and snow particles
  ! Don't want this process when rain/snow mix is representing the water
  ! content of individual particles (in the melting layer, or for 
  ! supercooled rain formed from a T>0 warm layer, rather than a mix of rain 
  !----------------------------------------------------------------------
  ZPRACS(:)  = 0.0_JPRB
  ZEFF_PRACS = 1.0_JPRB  ! Collection efficiency of rain accretes snow
  ZDENSNOW   = 250._JPRB ! Density of ice particles (currently set to fixed kg/m3)
    
  IF (IP_RAIN_ACCRETES_SNOW == 1) THEN
 
    DO JL=KIDIA,KFDIA
 
      IF ( ZTP1(JL,JK) <= RTT .AND. PRAINFRAC_TOPRFZ(JL) < ZEPSEC &
        & .AND. ZRAINCLD(JL) > ZEPSEC .AND. ZSNOWCLD(JL) > ZEPSEC ) THEN
 
        ! Calculate differential fallspeed
        ZFALLDIFF = MAX(ABS(ZFALLSPEED(JL,NCLDQR)-ZFALLSPEED(JL,NCLDQS)), &
                  & (ZFALLSPEED(JL,NCLDQR)+ZFALLSPEED(JL,NCLDQS))/8._JPRB)

        ZFAC = RPI*RPI*ZDENSNOW*ZEFF_PRACS*ZRHO(JL)*ZN0R(JL)*ZN0S(JL)/(ZLAMR(JL)*ZLAMS(JL)**3)
        ZRAT = ZLAMS(JL)/ZLAMR(JL)
        ZPRACS(JL) = ZCOVPTOT(JL)*PTSPHY*ZFAC*ZFALLDIFF &
                   & *((0.5_JPRB*ZRAT+2._JPRB)*ZRAT+5._JPRB)/ZLAMS(JL)**3 
          
        ! Limit to amount of rain
        !ZPRACS(JL) = MIN(ZPRACS(JL),ZRAINCLD(JL))
          
        ZSOLQA(JL,NCLDQS,NCLDQR) = ZSOLQA(JL,NCLDQS,NCLDQR) + ZPRACS(JL)
        ZSOLQA(JL,NCLDQR,NCLDQS) = ZSOLQA(JL,NCLDQR,NCLDQS) - ZPRACS(JL)
         
      ENDIF
            
    ENDDO

  ENDIF ! on IP_RAIN_ACCRETES_SNOW=1
  
  !======================================================================
  !
  !
  ! PIACR: ICE ACCRETES RAIN -> ICE
  !
  !
  !======================================================================
  ! Only active if T<0degC and rain, ice are present
  !----------------------------------------------------------------------
  ZPIACR(:)  = 0.0_JPRB
  ZEFF_PIACR = 0.1_JPRB ! Collection efficiency of ice accretes rain

  IF (IP_ICE_ACCRETES_RAIN == 1) THEN

    DO JL=KIDIA,KFDIA
      IF ( ZTP1(JL,JK) <= RTT .AND. PRAINFRAC_TOPRFZ(JL) < ZEPSEC &
        & .AND. ZRAINCLD(JL) > ZEPSEC .AND. ZICECLD(JL) > ZEPSEC ) THEN

         ! Calculate differential fallspeed
        ZFALLDIFF = MAX(ABS(ZFALLSPEED(JL,NCLDQR)-ZFALLSPEED(JL,NCLDQI)), &
                    & (ZFALLSPEED(JL,NCLDQR)+ZFALLSPEED(JL,NCLDQI))/8._JPRB)
 
        ZFAC = RPI*RPI*RDENSWAT*ZEFF_PIACR*ZRHO(JL)*ZN0R(JL)*ZN0I(JL)/(ZLAMI(JL)*ZLAMR(JL)**3)
        ZRAT = ZLAMR(JL)/ZLAMI(JL)
        ZPIACR(JL) = ZCOVPTOT(JL)*PTSPHY*ZFAC*ZFALLDIFF &
                   & *((0.5_JPRB*ZRAT+2._JPRB)*ZRAT+5._JPRB)/ZLAMR(JL)**3 
          
        ! Limit to amount of rain
        !ZPIACR(JL) = MIN(ZPIACR(JL),ZRAINCLD(JL))
          
        ZSOLQA(JL,NCLDQI,NCLDQR) = ZSOLQA(JL,NCLDQI,NCLDQR) + ZPIACR(JL)
        ZSOLQA(JL,NCLDQR,NCLDQI) = ZSOLQA(JL,NCLDQR,NCLDQI) - ZPIACR(JL)
       
      ENDIF
    ENDDO

  ENDIF ! on IP_ICE_ACCRETES_RAIN=1

  !======================================================================
  !
  !
  ! PRACI: RAIN ACCRETES ICE -> SNOW
  !
  !
  !======================================================================
  ZPRACI(:)  = 0.0_JPRB
  ZEFF_PRACI = 0.1_JPRB  ! Collection efficiency of snow accretes rain
  ZDENICE    = 250._JPRB ! Density of ice particles (currently set to fixed kg/m3)
    
  IF (IP_RAIN_ACCRETES_ICE == 1) THEN
 
    DO JL=KIDIA,KFDIA
 
      IF(ZTP1(JL,JK) <= RTT .AND. ZRAINCLD(JL)>ZEPSEC &
                          & .AND. ZICECLD(JL)>ZEPSEC) THEN
 
        ! Calculate differential fallspeed
        ZFALLDIFF = MAX(ABS(ZFALLSPEED(JL,NCLDQR)-ZFALLSPEED(JL,NCLDQI)), &
                  & (ZFALLSPEED(JL,NCLDQR)+ZFALLSPEED(JL,NCLDQI))/8._JPRB)

        ZFAC = RPI*RPI*ZDENICE*ZEFF_PRACI*ZRHO(JL)*ZN0R(JL)*ZN0I(JL)/(ZLAMR(JL)*ZLAMI(JL)**3)
        ZRAT = ZLAMI(JL)/ZLAMR(JL)
        ZPRACI(JL) = ZCOVPTOT(JL)*PTSPHY*ZFAC*ZFALLDIFF &
                   & *((0.5_JPRB*ZRAT+2._JPRB)*ZRAT+5._JPRB)/ZLAMI(JL)**3 
          
        ! Limit to amount of rain
        !ZPRACI(JL) = MIN(ZPRACI(JL),ZICECLD(JL))
          
        ZSOLQA(JL,NCLDQI,NCLDQR) = ZSOLQA(JL,NCLDQI,NCLDQR) + ZPRACI(JL)
        ZSOLQA(JL,NCLDQR,NCLDQI) = ZSOLQA(JL,NCLDQR,NCLDQI) - ZPRACI(JL)
         
      ENDIF
            
    ENDDO

  ENDIF ! on IP_RAIN_ACCRETES_ICE=1
  
  ENDIF ! on IWARMRAIN > 1
  

  !======================================================================
  !
  !
  ! 4.15  SUBLIMATION OF SNOW
  !
  !
  !======================================================================
  ! Snow -> Vapour
  !----------------------------------------------------------------------

 IF (ISUBLSNOW == 1) THEN
  
  DO JL=KIDIA,KFDIA
    ZZRH=RPRECRHMAX+(1.0_JPRB-RPRECRHMAX)*ZCOVPMAX(JL)/MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))
    ZZRH=MIN(MAX(ZZRH,RPRECRHMAX),1.0_JPRB)
    ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSICE(JL,JK))/ &
    & MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))  

    !---------------------------------------------
    ! humidity in moistest ZCOVPCLR part of domain
    !---------------------------------------------
    ZQE=MAX(0.0_JPRB,MIN(ZQE,ZQSICE(JL,JK)))
    LLO1=ZCOVPCLR(JL)>ZEPSEC .AND. &
       & ZQXFG(JL,NCLDQS)>ZEPSEC .AND. &
       & ZQE<ZZRH*ZQSICE(JL,JK)

    IF(LLO1) THEN
      ! note: zpreclr is a rain flux a
      ZPRECLR=ZQXFG(JL,NCLDQS)*ZCOVPCLR(JL)/ &
       & SIGN(MAX(ABS(ZCOVPTOT(JL)*ZDTGDP(JL)),ZEPSILON),ZCOVPTOT(JL)*ZDTGDP(JL))

      !--------------------------------------
      ! actual microphysics formula in zbeta
      !--------------------------------------

      ZBETA1=SQRT(PAP(JL,JK)/ &
       & PAPH(JL,KLEV+1))/RVRFACTOR*ZPRECLR/ &
       & MAX(ZCOVPCLR(JL),ZEPSEC)

      ZBETA=RG*RPECONS*(ZBETA1)**0.5777_JPRB  

      ZDENOM=1.0_JPRB+ZBETA*PTSPHY*ZCORQSICE(JL)
      ZDPR = ZCOVPCLR(JL)*ZBETA*(ZQSICE(JL,JK)-ZQE)/ZDENOM*ZDP(JL)*ZRG_R
      ZDPEVAP=ZDPR*ZDTGDP(JL)

      !Apply SPP perturbations
      IF (LLPERT_SNOWSUBLIM) THEN
        ZDPEVAP = ZDPEVAP*EXP(PN1SNOWSUBLIM%MU(1)+PN1SNOWSUBLIM%XMAG(1)*PGP2DSPP(JL, IPSNOWSUBLIM))
      ENDIF

      !---------------------------------------------------------
      ! add evaporation term to explicit sink.
      ! this has to be explicit since if treated in the implicit
      ! term evaporation can not reduce snow to zero and model
      ! produces small amounts of snowfall everywhere. 
      !---------------------------------------------------------
      
      ! Evaporate snow
      ZEVAP = MIN(ZDPEVAP,ZQXFG(JL,NCLDQS))

      ZSOLQA(JL,NCLDQV,NCLDQS) = ZSOLQA(JL,NCLDQV,NCLDQS)+ZEVAP
      ZSOLQA(JL,NCLDQS,NCLDQV) = ZSOLQA(JL,NCLDQS,NCLDQV)-ZEVAP
      ZBUDI(JL,17) = -ZEVAP*ZQTMST
      
      !-------------------------------------------------------------
      ! Reduce the total precip coverage proportional to evaporation
      ! to mimic the previous scheme which had a diagnostic
      ! 2-flux treatment, abandoned due to the new prognostic precip
      !-------------------------------------------------------------
      ZCOVPTOT(JL) = MAX(RCOVPMIN,ZCOVPTOT(JL)-MAX(0.0_JPRB, &
     &              (ZCOVPTOT(JL)-ZA(JL,JK))*ZEVAP/ZQXFG(JL,NCLDQS)))
      
      !Update first guess field
      ZQXFG(JL,NCLDQS) = ZQXFG(JL,NCLDQS)-ZEVAP

    ENDIF
  ENDDO

  !---------------------------------------------------------
  ELSEIF (ISUBLSNOW == 2) THEN

 
   DO JL=KIDIA,KFDIA

    !-----------------------------------------------------------------------
    ! Calculate relative humidity limit for snow evaporation 
    !-----------------------------------------------------------------------
    ZZRH=RPRECRHMAX+(1.0_JPRB-RPRECRHMAX)*ZCOVPMAX(JL)/MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))
    ZZRH=MIN(MAX(ZZRH,RPRECRHMAX),1.0_JPRB)

    ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSICE(JL,JK))/ &
    & MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))  
     
    !---------------------------------------------
    ! humidity in moistest ZCOVPCLR part of domain
    !---------------------------------------------
    ZQE=MAX(0.0_JPRB,MIN(ZQE,ZQSICE(JL,JK)))

    LLO1=ZCOVPCLR(JL)>ZEPSEC .AND. &
       & ZSNOWCLDM1(JL)>ZEPSEC .AND. & 
       & ZQE<ZZRH*ZQSICE(JL,JK)

    IF(LLO1) THEN
      
      ! Calculate local precipitation (kg/kg)
      ZPRECLR = ZSNOWCLDM1(JL)
     
      ! Saturation vapour pressure with respect to ice phase
      ZVPICE = RV/RD*FOEEICE(ZTP1(JL,JK))

      ! Particle size distribution
      ! ZTCG increases Ni with colder temperatures - essentially a 
      ! Fletcher or Meyers scheme? 
      ZTCG=1.0_JPRB !v1 EXP(RCL_X3I*(273.15_JPRB-ZTP1(JL,JK))/8.18_JPRB)
      ! ZFACX1I modification is based on Andrew Barrett's results
      ZFACX1S = 1.0_JPRB !v1 (ZICE0/1.E-5_JPRB)**0.627_JPRB

      ZAPLUSB   = RCL_APB1*ZVPICE-RCL_APB2*ZVPICE*ZTP1(JL,JK)+ &
     &             PAP(JL,JK)*RCL_APB3*ZTP1(JL,JK)**3
      ZCORRFAC  = (1.0/ZRHO(JL))**0.5
      ZCORRFAC2 = ((ZTP1(JL,JK)/273.0)**1.5)*(393.0/(ZTP1(JL,JK)+120.0))

      ZPR02 = ZRHO(JL)*ZPRECLR*RCL_CONST1S/(ZTCG*ZFACX1S)

      ZTERM1 = (ZQSICE(JL,JK)-ZQE)*ZTP1(JL,JK)**2*ZVPICE*ZCORRFAC2*ZTCG* &
     &          RCL_CONST2S*ZFACX1S/(ZRHO(JL)*ZAPLUSB*ZQSICE(JL,JK))
      ZTERM2 = 0.65*RCL_CONST6S*ZPR02**RCL_CONST4S+RCL_CONST3S*ZCORRFAC**0.5 &
     &          *ZRHO(JL)**0.5*ZPR02**RCL_CONST5S/ZCORRFAC2**0.5

      ZDPEVAP = MAX(ZCOVPCLR(JL)*ZTERM1*ZTERM2*PTSPHY,0.0_JPRB)
 
      !--------------------------------------------------------------------
      ! Limit evaporation to snow amount
      !--------------------------------------------------------------------
      ZEVAP = MIN(ZDPEVAP,ZEVAPLIMICE(JL))
      ZEVAP = MIN(ZEVAP,ZQXFG(JL,NCLDQS))
            
      ZSOLQA(JL,NCLDQV,NCLDQS) = ZSOLQA(JL,NCLDQV,NCLDQS)+ZEVAP
      ZSOLQA(JL,NCLDQS,NCLDQV) = ZSOLQA(JL,NCLDQS,NCLDQV)-ZEVAP
      ZBUDI(JL,17) = -ZEVAP*ZQTMST
      
      !-------------------------------------------------------------
      ! Reduce the total precip coverage proportional to evaporation
      ! to mimic the previous scheme which had a diagnostic
      ! 2-flux treatment, abandoned due to the new prognostic precip
      !-------------------------------------------------------------
      ZCOVPTOT(JL) = MAX(RCOVPMIN,ZCOVPTOT(JL)-MAX(0.0_JPRB, &
     &              (ZCOVPTOT(JL)-ZA(JL,JK))*ZEVAP/ZQXFG(JL,NCLDQS)))
      
      !Update first guess field
      ZQXFG(JL,NCLDQS) = ZQXFG(JL,NCLDQS)-ZEVAP

    ENDIF    
  ENDDO
     
ENDIF ! on ISUBLSNOW

  !--------------------------------------
  ! Evaporate small precipitation amounts
  !--------------------------------------
  DO JM=1,NCLV
   IF (LLFALL(JM)) THEN 
!DIR$ IVDEP
    DO JL=KIDIA,KFDIA
      IF (ZQXFG(JL,JM)<RLMIN) THEN
        ZSOLQA(JL,NCLDQV,JM) = ZSOLQA(JL,NCLDQV,JM)+ZQXFG(JL,JM)
        ZSOLQA(JL,JM,NCLDQV) = ZSOLQA(JL,JM,NCLDQV)-ZQXFG(JL,JM)
      ENDIF
    ENDDO
   ENDIF
  ENDDO

  
  !######################################################################
  !
  !            5.  *** SOLVERS FOR A AND L ***
  !
  ! Use an implicit solution rather than exact solution.
  ! Solver is forward in time, upstream difference for advection.
  !######################################################################

  !======================================================================
  !
  ! 5.1 Solver for cloud cover
  !
  !======================================================================
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
    ZANEW=(ZA(JL,JK)+ZSOLAC(JL))/(1.0_JPRB+ZSOLAB(JL))
    ZANEW=MIN(ZANEW,1.0_JPRB)
    IF (ZANEW<RAMIN) ZANEW=0.0_JPRB
    ZDA(JL)=ZANEW-ZAORIG(JL,JK)
    !---------------------------------
    ! variables needed for next level
    !---------------------------------
    ZANEWM1(JL)=ZANEW
  ENDDO

  !======================================================================
  !
  ! 5.2 Solver for the microphysics
  !
  !======================================================================

  !--------------------------------------------------------------
  ! 5.2.1 Truncate explicit sinks to avoid negatives 
  ! Note: Species are treated in the order in which they run out
  ! since the clipping will alter the balance for the other vars
  !--------------------------------------------------------------

  !----------------------------
  ! compute sink terms
  !----------------------------
  DO JM=1,NCLV
    DO JL=KIDIA,KFDIA
      ZSINKSUM(JL)=0.0_JPRB
    ENDDO
    DO JN=1,NCLV
      DO JL=KIDIA,KFDIA
        ZSINKSUM(JL)=ZSINKSUM(JL)-ZSOLQA(JL,JM,JN) ! +ve total is bad
      ENDDO
    ENDDO

    !----------------------------------------------
    ! calculate overshoot and scaling factor
    ! if ZSINKSUM is -ve, no overshoot and ZRATIO=1
    !----------------------------------------------
    DO JL=KIDIA,KFDIA
      ZMAX=MAX(ZQX(JL,JK,JM),ZEPSEC)
      ZRAT=MAX(ZSINKSUM(JL),ZMAX)
      ZRATIO(JL,JM)=ZMAX/ZRAT
    ENDDO
  ENDDO

  !--------------------------------------------------------
  ! sort zratio to find out which species run out first
  !--------------------------------------------------------
  DO JL=KIDIA,KFDIA
    IORDV(1)=1

    DO JM=2,NCLV
      ! Make room to move ZRV(JM) to its final place
      ! in the sorted sequence 
      ! ZRATIO(JL,IORDV(1)) ... ZRATIO(JL,IORDV(JM-1))
      DO JN=JM-1,1,-1
        IF (ZRATIO(JL,IORDV(JN))<=ZRATIO(JL,JM)) EXIT
        IORDV(JN+1)=IORDV(JN)
      ENDDO
        
      IORDV(JN+1)=JM
    ENDDO

    IORDER(JL,1:NCLV)=IORDV(1:NCLV)
  ENDDO

  !----------------
  ! recalculate sum
  !----------------
  DO JM=1,NCLV
    DO JL=KIDIA,KFDIA
      ZSINKSUM(JL)=0.0_JPRB
    ENDDO

    DO JN=1,NCLV
      DO JL=KIDIA,KFDIA
        JO=IORDER(JL,JM)
        LLINDEX3(JL,JN)=ZSOLQA(JL,JO,JN)<0.0_JPRB
        ZSINKSUM(JL)=ZSINKSUM(JL)-ZSOLQA(JL,JO,JN) ! +ve total is bad
      ENDDO
    ENDDO
    !---------------------------
    ! recalculate scaling factor
    !---------------------------
    DO JL=KIDIA,KFDIA
      JO=IORDER(JL,JM)
      ZMM=MAX(ZQX(JL,JK,JO),ZEPSEC)
      ZRR=MAX(ZSINKSUM(JL),ZMM)
      ZRATIO(JL,1)=ZMM/ZRR
    ENDDO
    !------
    ! scale
    !------
    DO JL=KIDIA,KFDIA
      JO=IORDER(JL,JM)
      ZZRATIO=ZRATIO(JL,1)
      !DIR$ IVDEP
      !DIR$ PREFERVECTOR
      DO JN=1,NCLV
        IF (LLINDEX3(JL,JN)) THEN
          ZSOLQA(JL,JO,JN)=ZSOLQA(JL,JO,JN)*ZZRATIO
          ZSOLQA(JL,JN,JO)=ZSOLQA(JL,JN,JO)*ZZRATIO
        ENDIF
      ENDDO
    ENDDO
  ENDDO


  !--------------------------------------------------------------
  ! 5.2.2 Solver
  !------------------------

  !------------------------
  ! set the LHS of equation  
  !------------------------
  DO JM=1,NCLV
     DO JN=1,NCLV
        !----------------------------------------------
        ! diagonals: microphysical sink terms+transport
        !----------------------------------------------
        IF (JN==JM) THEN
           DO JL=KIDIA,KFDIA
              ZQLHS(JL,JN,JM)=1.0_JPRB + ZFALLSINK(JL,JM)
           ENDDO
           !$OMP SIMD PRIVATE(JO)
           DO JL=KIDIA,KFDIA
!DIR$ UNROLL
              DO JO=1,NCLV
                 ZQLHS(JL,JN,JM)=ZQLHS(JL,JN,JM) + ZSOLQB(JL,JO,JN)
              ENDDO
           ENDDO
           !------------------------------------------
           ! non-diagonals: microphysical source terms
           !------------------------------------------
        ELSE
           DO JL=KIDIA,KFDIA
              ZQLHS(JL,JN,JM)= -ZSOLQB(JL,JN,JM) ! here is the delta T - missing from doc.
           ENDDO
        ENDIF
     ENDDO
  ENDDO

  !------------------------
  ! set the RHS of equation  
  !------------------------
  DO JM=1,NCLV
    DO JL=KIDIA,KFDIA
      !---------------------------------
      ! sum the explicit source and sink
      !---------------------------------
      ZEXPLICIT=0.0_JPRB
      DO JN=1,NCLV
        ZEXPLICIT=ZEXPLICIT+ZSOLQA(JL,JM,JN) ! sum over middle index
      ENDDO
      ZQXN(JL,JM)=ZQX(JL,JK,JM)+ZEXPLICIT
    ENDDO
  ENDDO

  !-----------------------------------
  ! *** solve by LU decomposition: ***
  !-----------------------------------
  ! Note: This fast way of solving NCLVxNCLV system
  !       assumes a good behaviour (i.e. non-zero diagonal
  !       terms with comparable orders) of the matrix stored
  !       in ZQLHS. For the moment this is the case but
  !       be aware to preserve it when doing eventual 
  !       modifications.

  ! Non pivoting recursive factorization 
  DO JN = 1, NCLV-1  ! number of steps
    DO JM = JN+1,NCLV ! row index
      ZQLHS(KIDIA:KFDIA,JM,JN)=ZQLHS(KIDIA:KFDIA,JM,JN) &
       &                     / ZQLHS(KIDIA:KFDIA,JN,JN)
      DO IK=JN+1,NCLV ! column index
        DO JL=KIDIA,KFDIA
          ZQLHS(JL,JM,IK)=ZQLHS(JL,JM,IK)-ZQLHS(JL,JM,JN)*ZQLHS(JL,JN,IK)
        ENDDO
      ENDDO
    ENDDO
  ENDDO        

  ! Backsubstitution 
  !  step 1 
  DO JN=2,NCLV
    DO JM = 1,JN-1
      ZQXN(KIDIA:KFDIA,JN)=ZQXN(KIDIA:KFDIA,JN)-ZQLHS(KIDIA:KFDIA,JN,JM) &
       &  *ZQXN(KIDIA:KFDIA,JM)
    ENDDO
  ENDDO
  !  step 2
  ZQXN(KIDIA:KFDIA,NCLV)=ZQXN(KIDIA:KFDIA,NCLV)/ZQLHS(KIDIA:KFDIA,NCLV,NCLV)
  DO JN=NCLV-1,1,-1
    DO JM = JN+1,NCLV
      ZQXN(KIDIA:KFDIA,JN)=ZQXN(KIDIA:KFDIA,JN)-ZQLHS(KIDIA:KFDIA,JN,JM) &
       &  *ZQXN(KIDIA:KFDIA,JM)
    ENDDO
    ZQXN(KIDIA:KFDIA,JN)=ZQXN(KIDIA:KFDIA,JN)/ZQLHS(KIDIA:KFDIA,JN,JN)
  ENDDO

  ! Ensure no small values (including negatives) remain in cloud variables nor
  ! precipitation rates.
  ! Evaporate l,i,r,s to water vapour. Latent heating taken into account below
  DO JN=1,NCLV-1
    DO JL=KIDIA,KFDIA
      IF (ZQXN(JL,JN) < ZEPSEC) THEN
        ZQXN(JL,NCLDQV) = ZQXN(JL,NCLDQV)+ZQXN(JL,JN)
        ZQXN(JL,JN)     = 0.0_JPRB
      ENDIF
    ENDDO
  ENDDO

  !--------------------------------
  ! variables needed for next level
  !--------------------------------
  DO JM=1,NCLV
    DO JL=KIDIA,KFDIA
      ZQXNM1(JL,JM)    = ZQXN(JL,JM)
      ZQXN2D(JL,JK,JM) = ZQXN(JL,JM)
    ENDDO
  ENDDO

  !------------------------------------------------------------------------
  ! 5.3 Precipitation/sedimentation fluxes to next level
  !     diagnostic precipitation fluxes
  !     It is this scaled flux that must be used for source to next layer
  !------------------------------------------------------------------------

  DO JM=1,NCLV
    DO JL=KIDIA,KFDIA
      ZPFPLSX(JL,JK+1,JM) = ZFALLSINK(JL,JM)*ZQXN(JL,JM)*ZRDTGDP(JL)
    ENDDO
  ENDDO

  ! Ensure precipitation fraction is zero if no precipitation
  DO JL=KIDIA,KFDIA
    ZQPRETOT(JL) =ZPFPLSX(JL,JK+1,NCLDQS)+ZPFPLSX(JL,JK+1,NCLDQR)
  ENDDO
  DO JL=KIDIA,KFDIA
    IF (ZQPRETOT(JL)<ZEPSEC) THEN
      ZCOVPTOT(JL)=0.0_JPRB
    ENDIF
  ENDDO
  
  ! Calculate diagnosed process rates for implicit quantities
  DO JL=KIDIA,KFDIA

    IF (IWARMRAIN == 3) THEN
      ! rain freezes instantly if T<0C
      IF (IP_SNOW_ACCRETES_RAIN == 2) THEN
        IF(ZTP1(JL,JK) <= RTT) THEN
          ZBUDL(JL,12) = -ZRAINAUT(JL)*ZQXN(JL,NCLDQL)*ZQTMST
          ZBUDL(JL,13) = -ZRAINACC(JL)*ZQXN(JL,NCLDQL)*ZQTMST
        ELSE
          ZBUDL(JL,14) = -ZRAINAUT(JL)*ZQXN(JL,NCLDQL)*ZQTMST
          ZBUDL(JL,15) = -ZRAINACC(JL)*ZQXN(JL,NCLDQL)*ZQTMST
        ENDIF
      ELSE ! no rain freezing from this process
        ZBUDL(JL,14) = -ZRAINAUT(JL)*ZQXN(JL,NCLDQL)*ZQTMST
        ZBUDL(JL,15) = -ZRAINACC(JL)*ZQXN(JL,NCLDQL)*ZQTMST    
      ENDIF
    ENDIF
    
    ZBUDCC(JL,6) = -ZCONVSINK(JL,NCLDQL)*ZANEWM1(JL)*ZQTMST
    ZBUDL(JL,6) = -ZCONVSINK(JL,NCLDQL)*ZQXN(JL,NCLDQL)*ZQTMST
    ZBUDI(JL,6) = -ZCONVSINK(JL,NCLDQI)*ZQXN(JL,NCLDQI)*ZQTMST
    ZBUDI(JL,13)= -ZFALLSINK(JL,NCLDQI)*ZQXN(JL,NCLDQI)*ZQTMST
    ZBUDL(JL,21)= -ZSNOWRIME(JL)*ZQXN(JL,NCLDQL)*ZQTMST
    ZBUDI(JL,14)= -ZSNOWAUT(JL)*ZQXN(JL,NCLDQI)*ZQTMST
  ENDDO


  !######################################################################
  !
  !              6.  *** UPDATE TENDANCIES ***
  !
  !######################################################################

  !----------------------------------------------
  ! 6.1 Temperature and cloud condensate budgets 
  !----------------------------------------------

  ! Loop over cloud hydrometeors
  DO JM=1,NCLV-1
  
    DO JL=KIDIA,KFDIA
      ! calculate fluxes in and out of box for conservation of TL
      ZFLUXQ(JL,JM) = ZCONVSRCE(JL,JM)+ZFALLSRCE(JL,JM)-&
                   & (ZFALLSINK(JL,JM)+ZCONVSINK(JL,JM))*ZQXN(JL,JM)
    ENDDO

    IF (IPHASE(JM)==1) THEN
      DO JL=KIDIA,KFDIA
        PTENDENCY_LOC_T(JL,JK)=PTENDENCY_LOC_T(JL,JK)+ &
          & RALVDCP*(ZQXN(JL,JM)-ZQX(JL,JK,JM)-ZFLUXQ(JL,JM))*ZQTMST
      ENDDO
    ENDIF

    IF (IPHASE(JM)==2) THEN
      DO JL=KIDIA,KFDIA
        PTENDENCY_LOC_T(JL,JK)=PTENDENCY_LOC_T(JL,JK)+ &
          & RALSDCP*(ZQXN(JL,JM)-ZQX(JL,JK,JM)-ZFLUXQ(JL,JM))*ZQTMST
      ENDDO
    ENDIF

      !----------------------------------------------------------------------
      ! New prognostic tendencies - ice,liquid rain,snow 
      ! Note: CLV arrays use PCLV in calculation of tendency while humidity
      !       uses ZQX. This is due to clipping at start of cloudsc which
      !       include the tendency already in PTENDENCY_LOC_T and PTENDENCY_LOC_q. ZQX was reset
      !----------------------------------------------------------------------
    DO JL=KIDIA,KFDIA
      PTENDENCY_LOC_CLD(JL,JK,JM)=PTENDENCY_LOC_CLD(JL,JK,JM)+(ZQXN(JL,JM)-ZQX0(JL,JK,JM))*ZQTMST
    ENDDO

  ENDDO

!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
  
    !------------------------------
    ! 6.2 Update humidity tendency
    !------------------------------
    PTENDENCY_LOC_Q(JL,JK)=PTENDENCY_LOC_Q(JL,JK)+(ZQXN(JL,NCLDQV)-ZQX(JL,JK,NCLDQV))*ZQTMST

    !-------------------
    ! 6.3  Update cloud cover tendency 
    !-----------------------
    PTENDENCY_LOC_A(JL,JK)=PTENDENCY_LOC_A(JL,JK)+ZDA(JL)*ZQTMST

    !------------------------------------------------------------------
    ! 6.4 Final check for supersaturation based on updated state
    !------------------------------------------------------------------

    ! Calculate updated state
    ZQ_UPD(JL) = PQ(JL,JK) + PTSPHY*(PTENDENCY_LOC_Q(JL,JK)+PTENDENCY_CML_Q(JL,JK))
    ZT_UPD(JL) = PT(JL,JK) + PTSPHY*(PTENDENCY_LOC_T(JL,JK)+PTENDENCY_CML_T(JL,JK))
    ZA_UPD(JL) = PA(JL,JK) + PTSPHY*(PTENDENCY_LOC_A(JL,JK)+PTENDENCY_CML_A(JL,JK))

  ENDDO

  ! Saturation adjustment step if there is any supersaturation above the defined limit
  IFTLIQICE = 1  ! Distribute liquid and ice according to mixed phase function
  !IFTLIQICE = 2  ! Distribute liquid and ice, all liquid warmer than homog freezing temperature
  CALL CLOUD_SUPERSATCHECK(YDECLDP, KIDIA, KFDIA, KLON, KLEV, IFTLIQICE, &
                         & ZT_UPD, ZQ_UPD, ZA_UPD, &
                         & PAP(KIDIA:KFDIA,JK), PAPH(KIDIA:KFDIA,KLEV+1), &
                         & ZT_ADJ, ZQ_ADJ, ZA_ADJ, ZL_ADJ, ZI_ADJ)
                             
  DO JL=KIDIA,KFDIA
      
    PTENDENCY_LOC_Q(JL,JK) = PTENDENCY_LOC_Q(JL,JK) + ZQ_ADJ(JL)*ZQTMST
    PTENDENCY_LOC_T(JL,JK) = PTENDENCY_LOC_T(JL,JK) + ZT_ADJ(JL)*ZQTMST
    PTENDENCY_LOC_A(JL,JK) = PTENDENCY_LOC_A(JL,JK) + ZA_ADJ(JL)*ZQTMST
    PTENDENCY_LOC_CLD(JL,JK,NCLDQL) = PTENDENCY_LOC_CLD(JL,JK,NCLDQL) + ZL_ADJ(JL)*ZQTMST
    PTENDENCY_LOC_CLD(JL,JK,NCLDQI) = PTENDENCY_LOC_CLD(JL,JK,NCLDQI) + ZI_ADJ(JL)*ZQTMST

    ZBUDCC(JL,2) = ZA_ADJ(JL)*ZQTMST
    ZBUDL(JL,2)  = ZL_ADJ(JL)*ZQTMST
    ZBUDI(JL,2)  = ZI_ADJ(JL)*ZQTMST

    ! Tidy up any small negative values due to numerical truncation
    !Cloud liquid
    ZL_ADJ(JL) = ZQX0(JL,JK,NCLDQL)+PTENDENCY_LOC_CLD(JL,JK,NCLDQL)*PTSPHY
    IF (ZL_ADJ(JL) < -RLMIN) WRITE(NULOUT,*) 'CLOUDSC: WARNING Negative cloud liquid! ',JK,ZL_ADJ(JL)
    IF (ZL_ADJ(JL) < 0.0_JPRB) THEN
      PTENDENCY_LOC_CLD(JL,JK,NCLDQL) = PTENDENCY_LOC_CLD(JL,JK,NCLDQL) - ZL_ADJ(JL)*ZQTMST
      PTENDENCY_LOC_Q(JL,JK) = PTENDENCY_LOC_Q(JL,JK) + ZL_ADJ(JL)*ZQTMST
      PTENDENCY_LOC_T(JL,JK) = PTENDENCY_LOC_T(JL,JK) - RALVDCP*ZL_ADJ(JL)*ZQTMST
    ENDIF
    ! Cloud ice
    ZI_ADJ(JL) = ZQX0(JL,JK,NCLDQI)+PTENDENCY_LOC_CLD(JL,JK,NCLDQI)*PTSPHY
    IF (ZI_ADJ(JL) < -RLMIN) WRITE(NULOUT,*) 'CLOUDSC: WARNING Negative cloud ice! ',JK,ZI_ADJ(JL)
    IF (ZI_ADJ(JL) < 0.0_JPRB) THEN
      PTENDENCY_LOC_CLD(JL,JK,NCLDQI) = PTENDENCY_LOC_CLD(JL,JK,NCLDQI) - ZI_ADJ(JL)*ZQTMST
      PTENDENCY_LOC_Q(JL,JK) = PTENDENCY_LOC_Q(JL,JK) + ZI_ADJ(JL)*ZQTMST
      PTENDENCY_LOC_T(JL,JK) = PTENDENCY_LOC_T(JL,JK) - RALSDCP*ZI_ADJ(JL)*ZQTMST
    ENDIF

  ENDDO

  !--------------------------------------------------
  ! Copy precipitation fraction into output variable
  !-------------------------------------------------
  DO JL=KIDIA,KFDIA
    PCOVPTOT(JL,JK) = ZCOVPTOT(JL)
  ENDDO

  ! cloud fraction/precip fraction at the end of the routine
  DO JL=KIDIA,KFDIA
    IF (PA(JL,JK)+ZDA(JL) > 0.001_JPRB) THEN !same cloud fraction threshold as in radiation
      ZANEWP(JL,JK) = MAX(MAX(PA(JL,JK)+ZDA(JL),PCOVPTOT(JL,JK)),0.1_JPRB)
    ENDIF
  ENDDO

!######################################################################
!
!              7.  *** CLOUD BUDGET DIAGNOSTICS ***
!
!######################################################################

  ! set pointer in extra variables array
  IF (LBUD23) THEN
    IS = 26  ! Take account of LBUD23 array diagnostics if turned on
  ELSE
    IS = 0
  ENDIF

  !-----------------------------------------------------------------
  ! Vertical integral of all cloud process terms in one 3D field 
  ! Requires certain number of levels.
  ! At some point need to move this to individual 2D diagnostic fields
  !-----------------------------------------------------------------
  IF (LCLDBUD_VERTINT) THEN

   IF (KLEV < 60) CALL ABOR1('CLOUDSC ERROR: Not enough levels for cloud vertical integral budget.')
   
   DO JL=KIDIA,KFDIA

      ! Layer depth (m)
      ZDZ = ZDP(JL)/(ZRHO(JL)*RG)

      IK = 0
      !PEXTRA(JL,IK+1,1)  = PEXTRA(JL,IK+1,1)  + ZBUDCC(JL,10)*ZDZ  ! + Condensation of new cloud
      !PEXTRA(JL,IK+2,1)  = PEXTRA(JL,IK+2,1)  + ZBUDCC(JL,10)*ZDZ  ! + Evaporation of cloud
      !PEXTRA(JL,IK+3,1)  = PEXTRA(JL,IK+3,1)  + ZBUDCC(JL,2)*ZDZ   ! + Supersat clipping after cloud_satadj
      PEXTRA(JL,IK+4,1)  = PEXTRA(JL,IK+4,1)  + ZBUDCC(JL,2)*ZDZ  ! + Supersat clipping after cloudsc
      !PEXTRA(JL,IK+5,1)  = PEXTRA(JL,IK+5,1)  + ZBUDCC(JL,2)*ZDZ   ! + Supersat clipping after sltend
      PEXTRA(JL,IK+6,1)  = PEXTRA(JL,IK+6,1)  + ZBUDCC(JL,3)*ZDZ   ! + Convective detrainment
      PEXTRA(JL,IK+7,1)  = PEXTRA(JL,IK+7,1)  + (ZBUDCC(JL,4)+ZBUDCC(JL,6))*ZDZ ! +- Convective subsidence source and sink
      PEXTRA(JL,IK+8,1)  = PEXTRA(JL,IK+8,1)  + ZBUDCC(JL,5)*ZDZ   ! - Convective subsidence evaporation
      PEXTRA(JL,IK+9,1)  = PEXTRA(JL,IK+9,1)  + ZBUDCC(JL,7)*ZDZ   ! - Turbulent erosion
      PEXTRA(JL,IK+10,1) = PEXTRA(JL,IK+10,1) + ZBUDCC(JL,12)*ZDZ  ! Tidy up
      PEXTRA(JL,IK+11,1) = PEXTRA(JL,IK+11,1) + PVFA(JL,JK)*ZDZ    ! Vertical diffusion
      PEXTRA(JL,IK+12,1) = PEXTRA(JL,IK+12,1) + PDYNA(JL,JK)*ZDZ   ! Advection from dynamics
      
      IK = 14
      !PEXTRA(JL,IK+1,1)  = PEXTRA(JL,IK+1,1) + ZBUDL(JL,10)*ZDZ ! + Condensation of new cloud (dqs decreasing = supersat)
      !PEXTRA(JL,IK+2,1)  = PEXTRA(JL,IK+2,1) + ZBUDL(JL,9)*ZDZ  ! + Condensation of existing cloud (dqs decreasing = supersat)
      !PEXTRA(JL,IK+3,1)  = PEXTRA(JL,IK+3,1) + ZBUDL(JL,8)*ZDZ  ! - Evaporation of existing cloud (dqs increasing = subsat)
      !PEXTRA(JL,IK+4,1)  = PEXTRA(JL,IK+4,1) + ZBUDL(JL,1)*ZDZ  ! + Supersat clipping after cloud_satadj
      PEXTRA(JL,IK+5,1)  = PEXTRA(JL,IK+5,1) + ZBUDL(JL,2)*ZDZ  ! + Supersat clipping after cloudsc
      !PEXTRA(JL,IK+6,1)  = PEXTRA(JL,IK+6,1) + ZBUDL(JL,2)*ZDZ  ! + Supersat clipping after sltend
      PEXTRA(JL,IK+7,1)  = PEXTRA(JL,IK+7,1) + ZBUDL(JL,3)*ZDZ ! + Convective detrainment
      PEXTRA(JL,IK+8,1)  = PEXTRA(JL,IK+8,1) + (ZBUDL(JL,4)+ZBUDL(JL,6))*ZDZ ! +- Convective subsidence source and sink
      PEXTRA(JL,IK+9,1)  = PEXTRA(JL,IK+9,1) + ZBUDL(JL,5)*ZDZ ! +- Evaporation due to convective subsidence
      PEXTRA(JL,IK+10,1) = PEXTRA(JL,IK+10,1) + ZBUDL(JL,7)*ZDZ ! - Turbulent erosion
      PEXTRA(JL,IK+11,1) = PEXTRA(JL,IK+11,1) + ZBUDL(JL,11)*ZDZ ! - Deposition of liquid to ice
      PEXTRA(JL,IK+12,1) = PEXTRA(JL,IK+12,1) + ZBUDL(JL,14)*ZDZ ! - Autoconversion to rain (IMPLICIT)(ZRAINAUT)
      PEXTRA(JL,IK+13,1) = PEXTRA(JL,IK+13,1) + ZBUDL(JL,15)*ZDZ ! - Accretion of cloud to rain (IMPLICIT)(ZRAINACC)
      PEXTRA(JL,IK+14,1) = PEXTRA(JL,IK+14,1) + ZBUDL(JL,12)*ZDZ ! - Autoconversion to rain+freezing->snow (ZRAINAUT)
      PEXTRA(JL,IK+15,1) = PEXTRA(JL,IK+15,1) + ZBUDL(JL,13)*ZDZ ! - Accretion of cloud to rain+freezing->snow (ZRAINACC)
      PEXTRA(JL,IK+16,1) = PEXTRA(JL,IK+16,1) + ZBUDL(JL,17)*ZDZ ! + Melting of ice to liquid
      PEXTRA(JL,IK+17,1) = PEXTRA(JL,IK+17,1) + (ZBUDL(JL,18)+ZBUDL(JL,19))*ZDZ ! - Freezing of rain-to-snow, liq-to-ice
      PEXTRA(JL,IK+18,1) = PEXTRA(JL,IK+18,1) + ZBUDL(JL,20)*ZDZ ! - Evaporation of rain
      PEXTRA(JL,IK+19,1) = PEXTRA(JL,IK+19,1) + ZBUDL(JL,21)*ZDZ ! - Riming of cloud liquid to snow
      PEXTRA(JL,IK+20,1) = PEXTRA(JL,IK+20,1) + ZPRACS(JL)*ZDZ*ZQTMST ! Rain accretes snow -> snow
      PEXTRA(JL,IK+21,1) = PEXTRA(JL,IK+21,1) + ZPSACR(JL)*ZDZ*ZQTMST ! Snow accretes rain -> snow
      PEXTRA(JL,IK+22,1) = PEXTRA(JL,IK+22,1) + ZPRACI(JL)*ZDZ*ZQTMST ! Rain accretes ice -> ice
      PEXTRA(JL,IK+23,1) = PEXTRA(JL,IK+23,1) + ZPIACR(JL)*ZDZ*ZQTMST ! Ice accretes rain -> ice
      PEXTRA(JL,IK+24,1) = PEXTRA(JL,IK+24,1) + PVFL(JL,JK)*ZDZ  ! +- Vertical diffusion
      PEXTRA(JL,IK+25,1) = PEXTRA(JL,IK+25,1) + PDYNL(JL,JK)*ZDZ ! +- Advection of liquid from dynamics
      !PEXTRA(JL,IK+26,1) = PEXTRA(JL,IK+26,1) + PDYNR(JL,JK)*ZDZ ! +- Advection of rain from dyn

      IK = 39
      !PEXTRA(JL,IK+1,1)  = PEXTRA(JL,IK+1,1) + ZBUDI(JL,10)*ZDZ ! + Condensation of new cloud (dqs decreasing = supersat)
      !PEXTRA(JL,IK+2,1)  = PEXTRA(JL,IK+2,1) + ZBUDI(JL,9)*ZDZ  ! + Condensation of existing cloud (dqs decreasing = supersat)
      !PEXTRA(JL,IK+3,1)  = PEXTRA(JL,IK+3,1) + ZBUDI(JL,8)*ZDZ  ! - Evaporation of existing cloud (dqs increasing = subsat)
      !PEXTRA(JL,IK+4,1)  = PEXTRA(JL,IK+4,1) + ZBUDI(JL,1)*ZDZ  ! + Supersat clipping after cloud_satadj
      PEXTRA(JL,IK+5,1)  = PEXTRA(JL,IK+5,1) + ZBUDI(JL,2)*ZDZ  ! + Supersat clipping after cloudsc
      !PEXTRA(JL,IK+6,1)  = PEXTRA(JL,IK+6,1) + ZBUDI(JL,2)*ZDZ  ! + Supersat clipping after sltend
      PEXTRA(JL,IK+7,1)  = PEXTRA(JL,IK+7,1) + ZBUDI(JL,3)*ZDZ  ! + Convective detrainment
      PEXTRA(JL,IK+8,1)  = PEXTRA(JL,IK+8,1) + (ZBUDI(JL,4)+ZBUDI(JL,6))*ZDZ ! +- Convective subsidence source and sink
      PEXTRA(JL,IK+9,1)  = PEXTRA(JL,IK+9,1) + ZBUDI(JL,5)*ZDZ  ! +- Evaporation due to convective subsidence
      PEXTRA(JL,IK+10,1) = PEXTRA(JL,IK+10,1) + ZBUDI(JL,7)*ZDZ  ! - Turbulent erosion
      PEXTRA(JL,IK+11,1) = PEXTRA(JL,IK+11,1) + ZBUDI(JL,11)*ZDZ ! + Deposition of liquid to ice
      PEXTRA(JL,IK+12,1) = PEXTRA(JL,IK+12,1) + (ZBUDI(JL,12)+ZBUDI(JL,13))*ZDZ  ! Ice sedimentation
      PEXTRA(JL,IK+13,1) = PEXTRA(JL,IK+13,1) + ZBUDI(JL,14)*ZDZ ! - Autoconversion to snow (IMPLICIT)(ZSNOWAUT)
      PEXTRA(JL,IK+14,1) = PEXTRA(JL,IK+14,1) + ZBUDI(JL,16)*ZDZ ! - Melting of snow to rain
      PEXTRA(JL,IK+15,1) = PEXTRA(JL,IK+15,1) + ZBUDI(JL,17)*ZDZ ! - Evaporation of rain
      PEXTRA(JL,IK+16,1) = PEXTRA(JL,IK+16,1) + ZPSDEP(JL)*ZQTMST*ZDZ !+ Deposition of vapour on snow
      PEXTRA(JL,IK+17,1) = PEXTRA(JL,IK+17,1) + PVFI(JL,JK)*ZDZ  ! +- Vertical diffusion
      PEXTRA(JL,IK+18,1) = PEXTRA(JL,IK+18,1) + PDYNI(JL,JK)*ZDZ ! +- Advection of ice from dynamics
      PEXTRA(JL,IK+19,1) = PEXTRA(JL,IK+19,1) + PDYNS(JL,JK)*ZDZ ! +- Advection of snow from dynamics
      
      IK = 60
      PEXTRA(JL,IK,1)   = PRAINFRAC_TOPRFZ(JL)
      
      IK = 70
      IF (KTYPE(JL) == 0)        PEXTRA(JL,IK,1)   = PEXTRA(JL,IK,1) + 1.0_JPRB
      IF (KTYPE(JL) == 1)        PEXTRA(JL,IK+1,1) = PEXTRA(JL,IK+1,1) + 1.0_JPRB
      IF (KTYPE(JL) == 2)        PEXTRA(JL,IK+2,1) = PEXTRA(JL,IK+2,1) + 1.0_JPRB
      IF (KTYPE(JL) == 3)        PEXTRA(JL,IK+3,1) = PEXTRA(JL,IK+3,1) + 1.0_JPRB
      IF (KPBLTYPE(JL) == 0)     PEXTRA(JL,IK+4,1) = PEXTRA(JL,IK+4,1) + 1.0_JPRB
      IF (KPBLTYPE(JL) == 1)     PEXTRA(JL,IK+5,1) = PEXTRA(JL,IK+5,1) + 1.0_JPRB
      IF (KPBLTYPE(JL) == 2)     PEXTRA(JL,IK+6,1) = PEXTRA(JL,IK+6,1) + 1.0_JPRB
      IF (KPBLTYPE(JL) == 3)     PEXTRA(JL,IK+7,1) = PEXTRA(JL,IK+7,1) + 1.0_JPRB
      PEXTRA(JL,IK+8,1) = PEXTRA(JL,IK+8,1) + PEIS(JL)
      IF (PEIS(JL) > 6.0_JPRB)  PEXTRA(JL,IK+9,1) = PEXTRA(JL,IK+9,1) + 1.0_JPRB
      IF (PEIS(JL) > 8.0_JPRB)  PEXTRA(JL,IK+10,1) = PEXTRA(JL,IK+10,1) + 1.0_JPRB
      IF (PEIS(JL) > 10.0_JPRB) PEXTRA(JL,IK+11,1) = PEXTRA(JL,IK+11,1) + 1.0_JPRB
      IF (KTYPE(JL) >= 2 .AND. PEIS(JL)<REISTHSC) PEXTRA(JL,IK+12,1) = PEXTRA(JL,IK+12,1) + 1.0_JPRB
      IF (KTYPE(JL) > 0) PEXTRA(JL,IK+13,1) = KCTOP(JL)   

    ENDDO
    IS = IS + 1
    IF (IS > KFLDX) CALL ABOR1('CLOUDSC ERROR: Not enough PEXTRA variables for cloud vertical integral budget.')
  ENDIF

  !-----------------------------------------------------------------
  ! Cloud fraction budget 
  !-----------------------------------------------------------------
  IF (LCLDBUDC) THEN
    
    DO JL=KIDIA,KFDIA
!      CVEXTRA(IS+1) = 'QSUPSATADJCLDSC'
!      PEXTRA(JL,JK,IS+1)  = PEXTRA(JL,JK,IS+1) + ZBUDCC(JL,10) ! Condensation of new cloud
!      PEXTRA(JL,JK,IS+2)  = PEXTRA(JL,JK,IS+2) + ZBUDCC(JL,11) ! Evaporation of cloud 
!      PEXTRA(JL,JK,IS+3)  = PEXTRA(JL,JK,IS+3) + ZBUDCC(JL,1)  ! Supersat clipping after cloud_satadj
      PEXTRA(JL,JK,IS+4)  = PEXTRA(JL,JK,IS+4)  + ZBUDCC(JL,2)   ! Supersat clipping after cloudsc
!      PEXTRA(JL,JK,IS+5)  = PEXTRA(JL,JK,IS+5) + ZBUDCC(JL,2)  ! Supersat clipping from t-1 sltend
      PEXTRA(JL,JK,IS+6)  = PEXTRA(JL,JK,IS+6)  + ZBUDCC(JL,3)  ! Convective detrainment
      PEXTRA(JL,JK,IS+7)  = PEXTRA(JL,JK,IS+7)  + ZBUDCC(JL,4)+ZBUDCC(JL,6)  ! Convective subsidence
      PEXTRA(JL,JK,IS+8)  = PEXTRA(JL,JK,IS+8)  + ZBUDCC(JL,5)  ! Convective subsidence evaporation
      PEXTRA(JL,JK,IS+9)  = PEXTRA(JL,JK,IS+9)  + ZBUDCC(JL,7)  ! Turbulent erosion
!      PEXTRA(JL,JK,IS+10)  = PEXTRA(JL,JK,IS+10)  + ZBUDCC(JL,?)  ! Tidyup
      PEXTRA(JL,JK,IS+11) = PEXTRA(JL,JK,IS+11) + PVFA(JL,JK)   ! Vertical diffusion
      PEXTRA(JL,JK,IS+12) = PEXTRA(JL,JK,IS+12) + PDYNA(JL,JK)  ! Advection from dynamics
    ENDDO
    IS = IS + 12
    IF (IS > KFLDX) CALL ABOR1('CLOUDSC ERROR: Not enough PEXTRA variables for cloud fraction budget.')
  ENDIF

  !-----------------------------------------------------------------
  ! Cloud liquid condensate budget 
  !-----------------------------------------------------------------
  IF (LCLDBUDL) THEN
    DO JL=KIDIA,KFDIA
!      PEXTRA(JL,JK,IS+1) = PEXTRA(JL,JK,IS+1) + ZBUDL(JL,10) ! + Condensation of new cloud (dqs decreasing = supersat)
!      PEXTRA(JL,JK,IS+2)  = PEXTRA(JL,JK,IS+2) + ZBUDL(JL,9)  ! + Condensation of existing cloud (dqs decreasing = supersat)
!      PEXTRA(JL,JK,IS+3)  = PEXTRA(JL,JK,IS+3) + ZBUDL(JL,8) ! - Evaporation of existing cloud (dqs increasing = subsat)
!      PEXTRA(JL,JK,IS+4)  = PEXTRA(JL,JK,IS+4) + ZBUDL(JL,2) ! + Supersat clipping after cloud_satadj
       PEXTRA(JL,JK,IS+5)  = PEXTRA(JL,JK,IS+5) + ZBUDL(JL,2) ! + Supersat clipping after cloudsc
!      PEXTRA(JL,JK,IS+6)  = PEXTRA(JL,JK,IS+6) + ZBUDL(JL,2) ! + Supersat clipping after sltend (PSUPSAT)
      PEXTRA(JL,JK,IS+7)  = PEXTRA(JL,JK,IS+7) + ZBUDL(JL,3) ! + Convective detrainment
      PEXTRA(JL,JK,IS+8)  = PEXTRA(JL,JK,IS+8) + ZBUDL(JL,4)+ZBUDL(JL,6)! +- Convective subsidence
      PEXTRA(JL,JK,IS+9)  = PEXTRA(JL,JK,IS+9) + ZBUDL(JL,5)! - Convective subsidence evaporation
       ! ZBUDL(JL,4) + Convective subsidence source from layer above
       ! ZBUDL(JL,5) - Convective subsidence source evaporation in layer
       ! ZBUDL(JL,6) - Convective subsidence sink to layer below (IMPLICIT) (ZCONVSINK)
      PEXTRA(JL,JK,IS+10)  = PEXTRA(JL,JK,IS+10) + ZBUDL(JL,7) ! - Turbulent erosion
      PEXTRA(JL,JK,IS+11) = PEXTRA(JL,JK,IS+11) + ZBUDL(JL,11)! - Deposition of liquid to ice
      PEXTRA(JL,JK,IS+12) = PEXTRA(JL,JK,IS+12) + ZBUDL(JL,12)! - Autoconversion to rain+freezing->snow (IMPLICIT)(ZRAINAUT) 
      PEXTRA(JL,JK,IS+13) = PEXTRA(JL,JK,IS+13) + ZBUDL(JL,13)! - Accretion of cloud to rain+freezing->snow (IMPLICIT)(ZRAINACC)
      PEXTRA(JL,JK,IS+14) = PEXTRA(JL,JK,IS+14) + ZBUDL(JL,14)! - Autoconversion to rain (IMPLICIT)(ZRAINAUT)
      PEXTRA(JL,JK,IS+15) = PEXTRA(JL,JK,IS+15) + ZBUDL(JL,15)! - Accretion of cloud to rain (IMPLICIT)(ZRAINACC)
      PEXTRA(JL,JK,IS+16) = PEXTRA(JL,JK,IS+16) + ZBUDL(JL,17)! + Melting of ice to liquid
      PEXTRA(JL,JK,IS+17) = PEXTRA(JL,JK,IS+17) + ZBUDL(JL,18)+ZBUDL(JL,19) ! - Freezing of rain/liq
       ! ZBUDL(JL,18) - Freezing of rain to snow
       ! ZBUDL(JL,19) - Freezing of liquid to ice
      PEXTRA(JL,JK,IS+18) = PEXTRA(JL,JK,IS+18) + ZBUDL(JL,20) ! - Evaporation of rain
      PEXTRA(JL,JK,IS+19) = PEXTRA(JL,JK,IS+19) + ZBUDL(JL,21) ! - Riming of cloud liquid to snow
      PEXTRA(JL,JK,IS+20) = PEXTRA(JL,JK,IS+20) + PVFL(JL,JK)  ! +- Vertical diffusion
      PEXTRA(JL,JK,IS+21) = PEXTRA(JL,JK,IS+21) + PDYNL(JL,JK) ! +- Advection from dynamics
      PEXTRA(JL,JK,IS+22) = PEXTRA(JL,JK,IS+22) + ZQXN(JL,NCLDQL)  ! Final condensate
    ENDDO
    IS = IS + 22
    IF (IS > KFLDX) CALL ABOR1('CLOUDSC ERROR: Not enough PEXTRA variables for cloud liquid budget.')
  ENDIF

  !-----------------------------------------------------------------
  ! Cloud ice condensate budget 
  !-----------------------------------------------------------------
  IF (LCLDBUDI) THEN
    DO JL=KIDIA,KFDIA
!      PEXTRA(JL,JK,IS+1) = PEXTRA(JL,JK,IS+1) + ZBUDI(JL,10) ! + Condensation of new cloud (dqs decreasing = supersat)
!      PEXTRA(JL,JK,IS+2)  = PEXTRA(JL,JK,IS+2) + ZBUDI(JL,9)  ! + Condensation of existing cloud (dqs decreasing = supersat)
!      PEXTRA(JL,JK,IS+3)  = PEXTRA(JL,JK,IS+3) + ZBUDI(JL,8)  ! - Evaporation of existing cloud (dqs increasing = subsat)
!      PEXTRA(JL,JK,IS+4)  = PEXTRA(JL,JK,IS+4) + ZBUDI(JL,2)  ! + Supersat clipping after cloud_satadj
      PEXTRA(JL,JK,IS+5)  = PEXTRA(JL,JK,IS+5) + ZBUDI(JL,2)   ! + Supersat clipping after cloudsc
!      PEXTRA(JL,JK,IS+6)  = PEXTRA(JL,JK,IS+6) + ZBUDI(JL,2)  ! + Supersat clipping after sltend
      PEXTRA(JL,JK,IS+7)  = PEXTRA(JL,JK,IS+7) + ZBUDI(JL,3)   ! + Convective detrainment
      PEXTRA(JL,JK,IS+8)  = PEXTRA(JL,JK,IS+8) + ZBUDI(JL,4)+ZBUDI(JL,6)! +- Convective subsidence
      PEXTRA(JL,JK,IS+9)  = PEXTRA(JL,JK,IS+9) + ZBUDI(JL,5)! - Convective subsidence evaporation
       ! ZBUDI(JL,4) + Convective subsidence source from layer above
       ! ZBUDI(JL,5) - Convective subsidence source evaporation in layer
       ! ZBUDI(JL,6) - Convective subsidence sink to layer below (IMPLICIT) (ZCONVSINK)
      PEXTRA(JL,JK,IS+10)  = PEXTRA(JL,JK,IS+10) + ZBUDI(JL,7)  ! - Turbulent erosion
      ! Microphysics
      PEXTRA(JL,JK,IS+11) = PEXTRA(JL,JK,IS+11) + ZBUDI(JL,11)! + Deposition of liquid to ice
      PEXTRA(JL,JK,IS+12) = PEXTRA(JL,JK,IS+12) + ZBUDI(JL,12)+ZBUDI(JL,13)! + Ice sedimentation
       ! ZBUDI(JL,12) + Ice sedimentation source from above
       ! ZBUDI(JL,13) - Ice sedimentation sink to below (IMPLICIT)(ZFALLSINK)
      PEXTRA(JL,JK,IS+13) = PEXTRA(JL,JK,IS+13) + ZBUDI(JL,14) ! - Autoconversion to snow (IMPLICIT) (ZSNOWAUT)
      PEXTRA(JL,JK,IS+14) = PEXTRA(JL,JK,IS+14) + ZBUDI(JL,15)+ZBUDI(JL,16) ! - Melting of ice/snow to rain
       ! ZBUDI(JL,15)! - Melting of ice to rain
       ! ZBUDI(JL,16)! - Melting of snow to rain
      PEXTRA(JL,JK,IS+15) = PEXTRA(JL,JK,IS+15) + ZPSDEP(JL)*ZQTMST ! Deposition of vapour to snow
      PEXTRA(JL,JK,IS+16) = PEXTRA(JL,JK,IS+16) + ZBUDI(JL,17)! - Evaporation of snow
      PEXTRA(JL,JK,IS+17) = PEXTRA(JL,JK,IS+17) + PVFI(JL,JK) ! Vertical diffusion
      PEXTRA(JL,JK,IS+18) = PEXTRA(JL,JK,IS+18) + ZPIEVAP(JL)*ZQTMST ! Evaporation of precipitating ice
      PEXTRA(JL,JK,IS+19) = PEXTRA(JL,JK,IS+19) + PDYNI(JL,JK)! Advection from dynamics
      PEXTRA(JL,JK,IS+20) = PEXTRA(JL,JK,IS+20) + ZQXN(JL,NCLDQI)  ! Final condensate
    ENDDO
    IS = IS + 20
    IF (IS > KFLDX) CALL ABOR1('CLOUDSC ERROR: Not enough PEXTRA variables for cloud ice budget.')
  ENDIF

  !-----------------------------------------------------------------
  ! Cloud processes temperature budget 
  !-----------------------------------------------------------------
  ! RALVDCP latent heat of condensation (vapour to liquid)
  ! RALSDCP latent heat of sublimation (vapour to solid)
  ! RALFDCP latent heat of melting/freezing (liquid to solid)
  !-----------------------------------------------------------------
  IF (LCLDBUDT) THEN
    DO JL=KIDIA,KFDIA
      ! Note, PEXTRA(:,:,1) is set to convective T tendency in callpar, so start at 2 here

      ! Radiative heating rates (shortwave and longwave)
      PEXTRA(JL,JK,IS+1) = PEXTRA(JL,JK,IS+1) + PHRSW(JL,JK)
      PEXTRA(JL,JK,IS+2) = PEXTRA(JL,JK,IS+2) + PHRLW(JL,JK)

      ! Condensation (vapour to liquid) (heating +ve RALVDCP) 
      !    = "Supersat adjust cloudsc" + "Supersat adjust sltend t-1" 
      !    + "Existing cloud" + "New cloud"

      PEXTRA(JL,JK,IS+3)  = PEXTRA(JL,JK,IS+3) + (ZBUDL(JL,1) + ZBUDL(JL,2) &
     &                     + ZBUDL(JL,9) + ZBUDL(JL,10))*RALVDCP&
     &                     -ZBUDI(JL,11)*RALVDCP

      ! Condensation (vapour to ice) (heating +ve RALSDCP) 
      !    = "Supersat adjust cloudsc" + "Supersat adjust sltend t-1" 
      !    + "Existing cloud" + "New cloud"
      PEXTRA(JL,JK,IS+4)  = PEXTRA(JL,JK,IS+4) + (ZBUDI(JL,1) + ZBUDI(JL,2) &
     &                     + ZBUDI(JL,9) + ZBUDI(JL,10)+ZBUDI(JL,11))*RALSDCP

      ! Deposition of liquid to ice  (Bergeron Findeisen liquid-vapour-ice) (heating +ve RALFDCP)        
      PEXTRA(JL,JK,IS+5) = PEXTRA(JL,JK,IS+5) + ZBUDI(JL,11) * RALFDCP                                      
      	
      ! Evaporation (liquid to vapour) (cooling -ve RALVDCP)
      !    = "Turbulent erosion" + "Evaporation of existing cloud" 
      !    + "Convective subsidence evap" 
      PEXTRA(JL,JK,IS+6)  = PEXTRA(JL,JK,IS+6) + (ZBUDL(JL,7) + ZBUDL(JL,8) + ZBUDL(JL,5))*RALVDCP

      ! Evaporation (ice to vapour) (cooling -ve RALSDCP)
      !    = "Turbulent erosion" + "Evaporation of existing cloud" 
      !    + "Convective subsidence evap" 
      PEXTRA(JL,JK,IS+7)  = PEXTRA(JL,JK,IS+7) + (ZBUDI(JL,7) + ZBUDI(JL,8) + ZBUDI(JL,5))*RALSDCP

      ! Evaporation of rain to vapour (cooling -ve RALVDCP)
      PEXTRA(JL,JK,IS+8) = PEXTRA(JL,JK,IS+8) + ZBUDL(JL,20) * RALVDCP

      ! Evaporation of snow to vapour (cooling -ve RALSDCP)
      PEXTRA(JL,JK,IS+9) = PEXTRA(JL,JK,IS+9) + ZBUDI(JL,17) * RALSDCP 

      ! Melting of ice to liquid (cooling -ve RALFDCP)
      PEXTRA(JL,JK,IS+10) = PEXTRA(JL,JK,IS+10) + ZBUDI(JL,15) * RALFDCP                                       

      ! Melting of snow to rain (cooling -ve RALFDCP)
      PEXTRA(JL,JK,IS+11) = PEXTRA(JL,JK,IS+11) + ZBUDI(JL,16) * RALFDCP                                       
             
      ! Freezing rate of liquid to ice (heating +ve RALFDCP)
      !    "Autoconversion & Accretion to rain+freezing->snow" 
      !    + "Freezing of rain to snow" + "Freezing of liquid to ice"
      PEXTRA(JL,JK,IS+12) = PEXTRA(JL,JK,IS+12) - (ZBUDL(JL,12) + ZBUDL(JL,13) + ZBUDL(JL,18) + &
     &                       ZBUDL(JL,19)) * RALFDCP 

      ! Riming of liquid to snow
      PEXTRA(JL,JK,IS+13) = PEXTRA(JL,JK,IS+13) + ZBUDL(JL,21) * RALFDCP                                       
      
    ENDDO
    IS = IS + 13
    IF (IS > KFLDX) CALL ABOR1('CLOUDSC ERROR: Not enough PEXTRA variables for cloud temperature budget.')
  ENDIF

  IF (YDEPHY%LRAD_CLOUD_INHOMOG) THEN
    DO JL=KIDIA,KFDIA
      ! for ice FSD
      !liquid, ice and detrained condensate at the beginning of timestep
      ! calculate ratio of detrained condensate to all condensate
      IF (ZQX0(JL,JK,NCLDQL)+ZQX0(JL,JK,NCLDQI)+PLUDE(JL,JK) > RLMIN) THEN
        ZRATFSD(JL,JK)=PLUDE(JL,JK)/(PLUDE(JL,JK)+ZQX0(JL,JK,NCLDQL)+ZQX0(JL,JK,NCLDQI))
        !safety - ratio between 0 and 1
        ZRATFSD(JL,JK)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZRATFSD(JL,JK)))
      ENDIF
      ! final liquid and ice condensate
      ZICEF(JL,JK)=ZQXN(JL,NCLDQI)
      ZLIQF(JL,JK)=ZQXN(JL,NCLDQL)
    ENDDO
  ENDIF

ENDDO ! on vertical level JK
!----------------------------------------------------------------------
!
!                  END OF VERTICAL LOOP OVER LEVELS
!
!----------------------------------------------------------------------



!######################################################################
!
!              8.  *** FLUX/DIAGNOSTICS COMPUTATIONS ***
!
!######################################################################

!-------------------------------------
! Enthalpy and total water diagnostics (LSCMEC and LCLDBUDGET normally false)
!-------------------------------------
IF (LSCMEC.OR.LCLDBUDGET) THEN

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZTNEW=PT(JL,JK)+PTSPHY*(PTENDENCY_LOC_T(JL,JK)+PTENDENCY_CML_T(JL,JK))
      IF (JK==1) THEN
        ZSUMQ1(JL,JK)=0.0_JPRB
        ZSUMH1(JL,JK)=0.0_JPRB
      ELSE
        ZSUMQ1(JL,JK)=ZSUMQ1(JL,JK-1)
        ZSUMH1(JL,JK)=ZSUMH1(JL,JK-1)
      ENDIF

      ! cld vars
      DO JM=1,NCLV-1
        IF (IPHASE(JM)==1) ZTNEW=ZTNEW-RALVDCP*(PCLV(JL,JK,JM)+ &
          & (PTENDENCY_LOC_CLD(JL,JK,JM)+PTENDENCY_CML_CLD(JL,JK,JM))*PTSPHY)
        IF (IPHASE(JM)==2) ZTNEW=ZTNEW-RALSDCP*(PCLV(JL,JK,JM)+ &
          & (PTENDENCY_LOC_CLD(JL,JK,JM)+PTENDENCY_CML_CLD(JL,JK,JM))*PTSPHY)
        ZSUMQ1(JL,JK)=ZSUMQ1(JL,JK)+ &
        & (PCLV(JL,JK,JM)+(PTENDENCY_LOC_CLD(JL,JK,JM)+PTENDENCY_CML_CLD(JL,JK,JM))*PTSPHY)* &
        & (PAPH(JL,JK+1)-PAPH(JL,JK))*ZRG_R
      ENDDO
      ZSUMH1(JL,JK)=ZSUMH1(JL,JK)+(PAPH(JL,JK+1)-PAPH(JL,JK))*ZTNEW 

      ! humidity
      ZSUMQ1(JL,JK)=ZSUMQ1(JL,JK)+ &
        &(PQ(JL,JK)+(PTENDENCY_LOC_Q(JL,JK)+PTENDENCY_CML_Q(JL,JK))*PTSPHY)*(PAPH(JL,JK+1)-PAPH(JL,JK))*ZRG_R

      ZRAIN=0.0_JPRB
      DO JM=1,NCLV
        ZRAIN=ZRAIN+PTSPHY*ZPFPLSX(JL,JK+1,JM)
      ENDDO
      ZERRORQ(JL,JK)=ZSUMQ1(JL,JK)+ZRAIN-ZSUMQ0(JL,JK)
    ENDDO
  ENDDO

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZDTGDP(JL)=PTSPHY*RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
      ZRAIN=0.0_JPRB
      DO JM=1,NCLV
        IF (IPHASE(JM)==1) ZRAIN=ZRAIN+RALVDCP*ZDTGDP(JL)*ZPFPLSX(JL,JK+1,JM)* &
                           & (PAPH(JL,JK+1)-PAPH(JL,JK))
        IF (IPHASE(JM)==2) ZRAIN=ZRAIN+RALSDCP*ZDTGDP(JL)*ZPFPLSX(JL,JK+1,JM)* &
                           & (PAPH(JL,JK+1)-PAPH(JL,JK))
      ENDDO
      ZSUMH1(JL,JK)=(ZSUMH1(JL,JK)-ZRAIN)/PAPH(JL,JK+1)
      ZERRORH(JL,JK)=ZSUMH1(JL,JK)-ZSUMH0(JL,JK)
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
    IF (ABS(ZERRORQ(JL,KLEV))>1.E-13_JPRB.OR.ABS(ZERRORH(JL,KLEV))>1.E-13_JPRB) THEN
      ZQADJ=0.0_JPRB ! dummy statement
                     ! place totalview break here to catch non-conservation
    ENDIF
  ENDDO

  IF (ALLOCATED(ZSUMQ0))  DEALLOCATE(ZSUMQ0)
  IF (ALLOCATED(ZSUMQ1))  DEALLOCATE(ZSUMQ1)
  IF (ALLOCATED(ZSUMH0))  DEALLOCATE(ZSUMH0)
  IF (ALLOCATED(ZSUMH1))  DEALLOCATE(ZSUMH1)
  IF (ALLOCATED(ZERRORQ)) DEALLOCATE(ZERRORQ)
  IF (ALLOCATED(ZERRORH)) DEALLOCATE(ZERRORH)

ENDIF

!--------------------------------------------------------------------
! Copy general precip arrays back into PFP arrays for GRIB archiving
! Add rain and liquid fluxes, ice and snow fluxes
!--------------------------------------------------------------------
DO JK=1,KLEV+1
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
    PFPLSL(JL,JK) = ZPFPLSX(JL,JK,NCLDQR)+ZPFPLSX(JL,JK,NCLDQL)
    PFPLSN(JL,JK) = ZPFPLSX(JL,JK,NCLDQS)+ZPFPLSX(JL,JK,NCLDQI)
  ENDDO
ENDDO

!--------
! Fluxes:
!--------
!DIR$ IVDEP
DO JL=KIDIA,KFDIA
  PFSQLF(JL,1)  = 0.0_JPRB
  PFSQIF(JL,1)  = 0.0_JPRB
  PFSQRF(JL,1)  = 0.0_JPRB
  PFSQSF(JL,1)  = 0.0_JPRB
  PFCQLNG(JL,1) = 0.0_JPRB
  PFCQNNG(JL,1) = 0.0_JPRB
  PFCQRNG(JL,1) = 0.0_JPRB !rain
  PFCQSNG(JL,1) = 0.0_JPRB !snow
! fluxes due to turbulence
  PFSQLTUR(JL,1) = 0.0_JPRB
  PFSQITUR(JL,1) = 0.0_JPRB
ENDDO

DO JK=1,KLEV
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA

    ZGDPH_R = -ZRG_R*(PAPH(JL,JK+1)-PAPH(JL,JK))*ZQTMST
    PFSQLF(JL,JK+1)  = PFSQLF(JL,JK)
    PFSQIF(JL,JK+1)  = PFSQIF(JL,JK)
    PFSQRF(JL,JK+1)  = PFSQLF(JL,JK)
    PFSQSF(JL,JK+1)  = PFSQIF(JL,JK)
    PFCQLNG(JL,JK+1) = PFCQLNG(JL,JK)
    PFCQNNG(JL,JK+1) = PFCQNNG(JL,JK)
    PFCQRNG(JL,JK+1) = PFCQLNG(JL,JK)
    PFCQSNG(JL,JK+1) = PFCQNNG(JL,JK)
    PFSQLTUR(JL,JK+1) = PFSQLTUR(JL,JK)
    PFSQITUR(JL,JK+1) = PFSQITUR(JL,JK)

    ZALFAW=ZFOEALFA(JL,JK)

    ! Liquid , LS scheme minus detrainment
    PFSQLF(JL,JK+1)=PFSQLF(JL,JK+1)+ &
   ! &(ZQXN2D(JL,JK,NCLDQL)-ZQX0(JL,JK,NCLDQL)+PVFL(JL,JK)*PTSPHY-PLUDELI(JL,JK,1))*ZGDPH_R
     &(ZQXN2D(JL,JK,NCLDQL)-ZQX0(JL,JK,NCLDQL)+PVFL(JL,JK)*PTSPHY)*ZGDPH_R+PLUDELI(JL,JK,1)
    ! liquid, negative numbers 
    PFCQLNG(JL,JK+1)=PFCQLNG(JL,JK+1)+ZLNEG(JL,JK,NCLDQL)*ZGDPH_R

    ! liquid, vertical diffusion
    PFSQLTUR(JL,JK+1)=PFSQLTUR(JL,JK+1)+PVFL(JL,JK)*PTSPHY*ZGDPH_R

    ! Rain, LS scheme 
    PFSQRF(JL,JK+1)=PFSQRF(JL,JK+1)+(ZQXN2D(JL,JK,NCLDQR)-ZQX0(JL,JK,NCLDQR))*ZGDPH_R 
    ! rain, negative numbers
    PFCQRNG(JL,JK+1)=PFCQRNG(JL,JK+1)+ZLNEG(JL,JK,NCLDQR)*ZGDPH_R

    ! Ice , LS scheme minus detrainment
    PFSQIF(JL,JK+1)=PFSQIF(JL,JK+1)+ &
   ! & (ZQXN2D(JL,JK,NCLDQI)-ZQX0(JL,JK,NCLDQI)+PVFI(JL,JK)*PTSPHY-PLUDELI(JL,JK,2))*ZGDPH_R
     & (ZQXN2D(JL,JK,NCLDQI)-ZQX0(JL,JK,NCLDQI)+PVFI(JL,JK)*PTSPHY)*ZGDPH_R+PLUDELI(JL,JK,2)
     ! ice, negative numbers
    PFCQNNG(JL,JK+1)=PFCQNNG(JL,JK+1)+ZLNEG(JL,JK,NCLDQI)*ZGDPH_R

    ! ice, vertical diffusion
    PFSQITUR(JL,JK+1)=PFSQITUR(JL,JK+1)+PVFI(JL,JK)*PTSPHY*ZGDPH_R

    ! snow, LS scheme
    PFSQSF(JL,JK+1)=PFSQSF(JL,JK+1)+(ZQXN2D(JL,JK,NCLDQS)-ZQX0(JL,JK,NCLDQS))*ZGDPH_R 
    ! snow, negative numbers
    PFCQSNG(JL,JK+1)=PFCQSNG(JL,JK+1)+ZLNEG(JL,JK,NCLDQS)*ZGDPH_R
  ENDDO
ENDDO

!-----------------------------------
! enthalpy flux due to precipitation
!-----------------------------------
DO JK=1,KLEV+1
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
    PFHPSL(JL,JK) = -RLVTT*PFPLSL(JL,JK)
    PFHPSN(JL,JK) = -RLSTT*PFPLSN(JL,JK)
  ENDDO
ENDDO

! Ice FSD calculation
!-----------------------------------
IF (YDEPHY%LRAD_CLOUD_INHOMOG) THEN
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ! Representative grid box length (km)
      ! PGAW = normalised gaussian quadrature weight / no. longitude pts
      ZGRIDLEN = 2*RA*SQRT(RPI*PGAW(JL))*0.001_JPRB
      ZIFSDBACK=ZU1*ZGRIDLEN**(1._JPRB/3._JPRB)*((ZU2*ZGRIDLEN)**ZU3+1._JPRB)**ZU4
       
      !param updated in Sept 2017.
      ZQTOT = max(min((ZQX(JL,JK,NCLDQV)+ZQX(JL,JK,NCLDQI))*1000._JPRB,30._JPRB),0.00001)
      ZDP(JL)     = PAPH(JL,JK+1)-PAPH(JL,JK)     ! dp
      ZRHO(JL)    = PAP(JL,JK)/(RD*ZTP1(JL,JK))   ! p/RT air density
      ZDELZ =1._JPRB/(RG*ZRHO(JL))*ZDP(JL)*0.001_JPRB  !layer thickness in km

      ZDELZ=max(min(ZDELZ,.5),.001)
      ZPHI=(ZDELZ/.24)**.11*(ZQTOT/10.)**.03

      IF (ZANEWP(JL,JK) > 0.95_JPRB) THEN !treat as overcast
        ZIFSD(JL,JK)=ZIFSDBACK*ZPHI
        !add detrainment enhancement, multiply by 1D-to-2D enhancement factor
        ZIFSD(JL,JK)=ZR12*(ZIFSD(JL,JK)+ZRATFSD(JL,JK)*1.5_JPRB)
      ELSEIF (ZANEWP(JL,JK) > 0.001_JPRB .AND. ZANEWP(JL,JK) <= 0.95_JPRB) THEN 
        ZIFSD(JL,JK)=ZIFSDBACK*ZPHI*1.5_JPRB
        !add detrainment enhancement, multiply by 1D-to-2D enhancement factor
        ZIFSD(JL,JK)=ZR12*(ZIFSD(JL,JK)+ZRATFSD(JL,JK)*1.5_JPRB)
      ENDIF 

      !assign liquid or ice fsd based on liquid fraction. 
      !Global default value set to 1. all other cases
      ZFSD(JL,JK)=YDERAD%RCLOUD_FRAC_STD
      IF (ZICEF(JL,JK) > ZEPSEC ) THEN
        ZFSD(JL,JK)=ZIFSD(JL,JK)
      ENDIF
      IF (ZLIQF(JL,JK) > 0._JPRB .AND. ZICEF(JL,JK) <= ZEPSEC) THEN
        ZFSD(JL,JK)=ZLFSD(JL,JK) 
      ENDIF

      ! consistent with limits of possible FSD [0.1,3.575] in 
      ! lookup tables for Gamma/Log-normal functions
      ZFSD(JL,JK)=MIN(MAX(0.1_JPRB,ZFSD(JL,JK)),3.575_JPRB)

    ENDDO
  ENDDO
ENDIF ! LRAD_CLOUD_INHOMOG

! assign FSD to PFSD variable
IF (YDEPHY%LRAD_CLOUD_INHOMOG) THEN
   DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
         PFSD(JL,JK)=ZFSD(JL,JK)
      ENDDO
   ENDDO
ENDIF

!===============================================================================
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CLOUDSC',1,ZHOOK_HANDLE)

END SUBROUTINE CLOUDSC
