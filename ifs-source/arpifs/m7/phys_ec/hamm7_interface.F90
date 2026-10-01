! (C) Copyright 2026- FMI.
! (C) Copyright 2026- KNMI.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

SUBROUTINE HAMM7_INTERFACE( &
 & YDMODEL,   KIDIA,     KFDIA,     KLON,        KTDIA,    KLEV,     KTILES,  &
 & KFLDX,     KLEVX,     KTRAC,     KAERO,       KCHEM,    KSTGLO,   PGEOH,   &
 & PRS1,      PRSF1,     PAEROP,    PCAERO,      PCEN,     PAPHIF,            &
 & PFPLCL,    PFPLCN,    PFPLSL,    PFPLSN,      PGELAT,   PGELAM,            &
 & PAP,       PIP,       PLP,       PRP,         PSP,      PCOVPTOT,          &
 & PLU,       PMFU,      PO3P,      PQP,         PTP,      PTHP,     PTENC,    PCFLX,   &
 & PAERDDP,   PAERSDM,   PAERSRC,   PAERWS,      PAERGUST, PAERUST,  PAERMAP, &
 & PCLAERS,   PPRAERS,   PCHEM2AER,                                           &
 & PFRTI,     PLSM,      PSNS,      PWND,        PWS1,     PAERFLX,  PAERLIF, &
 & PAERODDF,  PTSPHY,    PGFL,                                                &
 & PODTO,     PAERO_WVL_DIAG,                                                 &
 & PAER_TAU,  PAER_SSA,  PAER_ASYM, PAER_TAU_LW,                              &
 & PTAUS_AER, PTAUA_AER, PPMAER,                                              &
 & PEXTRA,    PVERVEL,   PCCNL,     PCCNO,       PAHFSTI,  PCI,      PZ0M,    &
 & PAHFLEV,   PUP,       PVP,       PCVL,        PCVH,     PSO2DD,   PGEMU,   &
 & PBLH,      PSNOWACL,  KCTOP,     KCBOT,                                    &
 & PAEPM1,    PAEPM25,   PAEPM10)

! ╭────────────────────────────────────────────────────────────────────────────╮
! │                                                      (updated 30-APR-2024) │
! │ Purpose :                                                                  │
! │ -------                                                                    │
! │  *hamm7_interface* -                                                       │
! │                                                                            │
! │                                                                            │
! │ Interface :                                                                │
! │ ---------                                                                  │
! │   *HAMM7_INTERFACE* is called from AER_PHY3_LAYER                          │
! │                                                                            │
! │                                                                            │
! │ Input :                                                                    │
! │ -----                                                                      │
! │                                                                            │
! │                                                                            │
! │ Output :                                                                   │
! │ ------                                                                     │
! │                                                                            │
! │                                                                            │
! │ Externals :                                                                │
! │ ---------                                                                  │
! │                                                                            │
! │ Method :                                                                   │
! │ ------                                                                     │
! │                                                                            │
! │ Reference :                                                                │
! │ ---------                                                                  │
! │                                                                            │
! │ Author :                                                                   │
! │ -------                                                                    │
! │     Orginal version:                                                       │
! │     Vicent Huijen (KNMI), Tommi Bergman (FMI), Thomas Kuehn (FMI/UEF)      │
! │                                                                            │
! │ Modifications :                                                            │
! │ -------------                                                              │
! │     May.  2020 - V. Huijnen     : Modifications for TM5M7                  │
! │     Sep.  2020 - T. Bergman     : TM5M7 work                               │
! │     Apr.  2024 - Lianghai Wu    : revision for CY48r1                      │
! │     May.  2024 - R. Checa-Garcia: revision for CY48r1 and refactoring      │
! │                                                                            │
! ╰────────────────────────────────────────────────────────────────────────────╯

! RCHG: (TODO)
!      This subroutine assumes YAEROUT() has a number of elements, we may have 
!      to introduce something that test that the number is consistent with what 
!      is needed. Also the namelist with YAEROUT specification may need to be  
!      consistent with how we fill things here (not sure about it). 
!
! RCHG: (TODO) 
!      For this subroutine the description of the steps is not very uniform 
!      For so long subroutines we may need a table of contents to any reader 
!      can understand what is doing or better per section headers or split 
!      with CONTAINS into different parts. The advantage of the last approach
!      is to indentify and isolate what are the input/outputs of each step.
!      The problem are the associate environments.


!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLEV*         NUMBER OF LEVELS
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEVS*        NUMBER OF SOIL LAYERS
!    *KTILES*       NUMBER OF TILES (I.E. SUBGRID AREAS WITH DIFFERENT 
!                   OF SURFACE BOUNDARY CONDITION)
!    *KVTYPES*      NUMBER OF biomes for land carbon
!    *KTRAC*        Number of tracers
!     KLEVXG  : number of levels to compute
!     KAERO                       : Number of aerosol fields 
!     KCHEM :  Number of chemistry tracers 
!     KFLDX   - number of the last field. of what??? TB
!     PRS1(KLON,0:KLEV)           : HALF-LEVEL PRESSURE           (Pa)
!     PRSF1(KLON,KLEV)             : FULL-LEVEL PRESSURE           (Pa)
!     PAEROP(KLON,KLEV,KAERO)     : Aerosol concentrations  (kg/kg) - Note that fields are only non-zero if NACTAERO > 0
!     PCAERO is the aerosol density                kg cm-3 (KWHAT=0)
!     PCAERO is the optical thickness              ND     (KWHAT=2)
!     PCEN(KLON,KLEV,KCHEM)       :  CONCENTRATION OF TRACERS           (kg/kg) TB:Chemistry
!     PAPHIF         : geopotential height "gz" at full levels.
!     PFPLCL (KLON,KLEV+1)        : CONVECTIVE PRECIPITATION AS RAIN  (kg/m2s)
!     PFPLCN   : convective precipitation as snow.
!    *PFPLSL*       LARGE SCALE RAIN FLUX                        KG/(M2*S)
!    *PFPLSN*       LARGE SCALE SNOW FLUX                        KG/(M2*S)
!     PGELAT(KLON)                : LATITUDE (RADIANS) 
!    *PGELAM*       LONGITUDE                                     RADIANS
!     PAP    : (KLON,KLEV)       ; CLOUD FRACTION
!     PIP     (KLON,KLEV)         :  ICWC                         (kg/kg)
!     PLP     (KLON,KLEV)         :  LWC                         (kg/kg)
!    *PSP*    (KLON,KLEV)         :  Snow water content           (kg/kg) 
!     PRP     (KLON,KLEV)         :  Rain water content           (kg/kg)
!     PCOVPTOT(KLON,KLEV)         :  PRECIP FRACTION               0..1  
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS            KG/KG
!     PQP(KLON,KLEV)      : FULL-LEVEL HUMIDITY (W. DYN.TEND.) (kg kg-1)
!     PTP     (KLON,KLEV)         : TEMPERATURE                   (K)!
!     pthp ?  ZPPTHPW(KPROMA,KOPLEV)      ! theta'w???
!     PTENC  (KLON,KLEV,KTRAC)     : TENDENCY OF CONCENTRATION OF TRACERS including chemistry(kg/kg s-1)
!     PCFLX(KLON,KTRAC)      : SURFACE FLUX OF TRACERS              (xx m-2)
!     PAERDDP(KLON,NACTAERO) : aerosol dry deposition 
!     PAERSDM(KLON,NACTAERO)  : aerosol sedimentation 
!     PAERSRC(KLON,NACTAERO) : SOURCE FLUX                          (xx m-2) 
!     PAERWS                 : Wind speed (average of horizontal wind speed)   (m/s)
!     PAERGUST               : Wind gust (maximum 3 second gust in the hour)   (m/s)
!     PAERMAP(KLON,5)        : DUST MASK-RELATED QUANTITIES
! ??? PCLAERS(KLON)    ! aerosol
!     PPRAERS(KLON) ! radvis.f90
!     PALBD  : (KPROMA,NTSW)        ; DIFFUSE ALBEDO IN THE NSW SW INTERVALS
!    *PFRTI*      TILE FRACTIONS                                   (0-1)
!     PLSM   (KLON) : land-sea mask                       [0-1]
!     PSNS       : MASS OF SNOW PER UNIT SURFACE
!     PWND(KLON)                  : Surface wind
!     PWS1   : REAL     : TOP LAYER SOIL MOISTURE CONTENT
!     PAERFLX(KLON,12,9)     : DIAGNOSTIC DUST SOURCE FLUXES
!     PAERLIF(KLON,9)        : DIAGNOSTIC LIFTING THRESHOLD SPEED
!     PAERODDF    ! Diagnostic array with aerosol fluxes (sources and sinks)
!     PTSPHY                : TIMESTEP                  (s)
!     PGFL         : GFL fields
!     PVERVEL      : Vertical velocity [Pa s-1]
!     PCCNL        : CCN over land
!     PCCNO        : CCN over ocean
!    *PAHFSTI*      SURFACE SENSIBLE HEAT FLUX                    W/M2
!     PCI(KLON)           : FRACTION OF SEA-ICE
!     PZ0M                : roughness length for momentum              (m)
!     PAHFLEV             : latent heat flux                           (W/m2)
!     PUP                 : u-wind (m/s)
!     PVP                 : v-wind (m/s)
!     PCVL                : low vegetation cover
!     PCVH                : high vegetation cover
!     PGEMU               : sine of latitude
!     PGEOH               : geopotential at half levels (M2/S2)
!------------

!-----------------------------------------------------------------------

USE PARKIND1,     ONLY: JPIM, JPRB
USE YOMHOOK,      ONLY: LHOOK, DR_HOOK, JPHOOK
USE TYPE_MODEL,   ONLY: MODEL
USE YOMCST,       ONLY: RD, RG, RPI
USE YOMCT0,       ONLY: LIFSMIN, LIFSTRAJ
USE YOMCT3,       ONLY: NSTEP
USE TM5M7_DATA,   ONLY: MODAL_DATA, MODE_TRACERS, MODE_START,      &
                      & MODE_TRACERS_BY_MODS, MODE_END_SO4, NAERMOD, NMOD, NSOL,&
                      & IISVOC, IELVOC, IACS_N, ISO4, ISO4ACS,ISO4COS
USE YOMCHEM,      ONLY: IEXTR_WD, IEXTR_CH, IEXTR_NG, IEXTR_DD, IEXTR_CHTR

USE YOE_AERODIAG, ONLY: JPAERO_WVL_AOD, JPAERO_WVL_AODABS, JPAERO_WVL_AODFM,   &
                      & JPAERO_WVL_SSA, JPAERO_WVL_ASSIMETRY
USE YOMLUN,       ONLY: NULOUT

! HAM-M7
USE MO_HAM,                  ONLY: nclass, naerocomp, sizeclass, nccndiag, subm_ngasspec
USE OIFS_TO_HAM,             ONLY: ind_oifs_ham
USE MO_HAM_SUBM,             ONLY: HAM_SUBM_INTERFACE                   ! replaced HAM-M7 call with submodel interface
USE MO_ACTIV,                ONLY: activ_updraft,nw, idt_cdnc, idt_icnc ! HAM-M7 activation updraft calculation, effective radii
USE MO_HAM_ACTIV,            ONLY: ham_activ_abdulrazzak_ghan, ham_activ_koehler_ab ! HAM-M7 activation
USE MO_PARAM_SWITCHES,       ONLY: ncd_activ        ! for activation
USE MO_TRACDEF,              ONLY: ntrac, trlist    ! number of tracer for mass/number mixing ratio conversion, trlist for wet deposition flags
USE MO_TRACER_PROCESSES,     ONLY: xt_borrow, xt_conv_massfix ! conserving the negative tracer values from tendency, and convective case
USE MO_TIME_CONTROL,         ONLY: time_step_len    ! time step length for tendency
USE MO_HAMMOZ_WETDEP,        ONLY: wetdep_interface ! wet deposition interface call
USE MO_HAM_WETDEP,           ONLY: ham_conv_lfraq_so2
USE MO_HAMMOZ_SEDIMENTATION, ONLY: sedi_interface   ! sedimentation interface call
USE MO_HAMMOZ_DRYDEP,        ONLY: drydep_interface ! dry deposition interface call
USE MO_HAM_RAD,              ONLY: ham_rad,ham_rad_cache_cleanup,ham_rad_cache

USE YOE_AER_ACTIV,           ONLY: AER_ACTIV ! M&N activation scheme

USE TM5M7_OPTICS_DATA,       ONLY : NWDEP, NASWBAND, ASWBAND !,WDEP, AER_TAU, AER_SSA,AER_ASYM,AER_TAU_LW
USE TM5_PHOTOLYSIS,          ONLY : NBANDS_TROP, WAV_GRID, WAV_GRIDA
USE TM5M7_EMIS_DATA,         ONLY : VKARMAN ! von karman constant for dry deposition
USE TM5_CHEM_MODULE,         ONLY : NCHEM2AER
!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL),        INTENT(IN) :: YDMODEL
INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA, KFDIA, KLON
INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA, KLEV, KFLDX, KLEVX
INTEGER(KIND=JPIM), INTENT(IN) :: KTILES
INTEGER(KIND=JPIM), INTENT(IN) :: KTRAC
INTEGER(KIND=JPIM), INTENT(IN) :: KAERO(YDMODEL%YRML_GCONF%YGFL%NAERO)
INTEGER(KIND=JPIM), INTENT(IN) :: KCHEM(YDMODEL%YRML_GCONF%YGFL%NCHEM)
INTEGER(KIND=JPIM), INTENT(IN) :: KSTGLO
INTEGER(KIND=JPIM), INTENT(IN) :: KCTOP(KLON) ! cloud top level index
INTEGER(KIND=JPIM), INTENT(IN) :: KCBOT(KLON) ! cloud bottom level index

REAL(KIND=JPRB),INTENT(IN)    :: PGEOH(KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PRSF1(KLON,KLEV), PRS1(KLON,0:KLEV), PAPHIF(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PAP(KLON,KLEV), PIP(KLON,KLEV), PLP(KLON,KLEV), PLU(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PRP(KLON,KLEV), PSP(KLON,KLEV), PCOVPTOT(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PAEROP(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO) 
REAL(KIND=JPRB),INTENT(IN)    :: PCAERO(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB),INTENT(IN)    :: PCEN(KLON,KLEV,KTRAC), PCFLX(KLON,KTRAC)
REAL(KIND=JPRB),INTENT(IN)    :: PO3P(KLON,KLEV), PQP(KLON,KLEV), PTP(KLON,KLEV), PTHP(KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PFPLCL(KLON,0:KLEV),PFPLCN(KLON,0:KLEV),PFPLSL(KLON,0:KLEV),PFPLSN(KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PGELAT(KLON)    , PGELAM(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PFRTI(KLON,KTILES)
REAL(KIND=JPRB),INTENT(IN)    :: PAERWS(KLON), PAERGUST(KLON), PAERUST(KLON), PAERMAP(KLON,5)
REAL(KIND=JPRB),INTENT(IN)    :: PAERSRC(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB),INTENT(IN)    :: PCHEM2AER(KLON,KLEV,NCHEM2AER)
REAL(KIND=JPRB),INTENT(IN)    :: PAERFLX(KLON,12,9), PAERLIF(KLON,9), PCLAERS(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PLSM(KLON)  , PSNS(KLON)    , PWND(KLON)   , PWS1(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB),INTENT(IN)    :: PVERVEL(KLON,KLEV) ! vertical velocity (needed in HAM-M7 activation)
REAL(KIND=JPRB),INTENT(IN)    :: PCCNL(KLON) ! CCN over land (needed in liquid effective radius calc.)
REAL(KIND=JPRB),INTENT(IN)    :: PCCNO(KLON) ! CCN over ocean (needed in liquid effective radius calc.)
REAL(KIND=JPRB),INTENT(IN)    :: PAHFSTI(KLON,KTILES) ! surface sensible heat flux for dry deposition
REAL(KIND=JPRB),INTENT(IN)    :: PCI(KLON) ! fraction of sea-ice for dry deposition
REAL(KIND=JPRB),INTENT(IN)    :: PZ0M(KLON) ! roughness length for momentum for dry deposition
REAL(KIND=JPRB),INTENT(IN)    :: PAHFLEV(KLON) ! latent heat flux for dry deposition
REAL(KIND=JPRB),INTENT(IN)    :: PUP(KLON,KLEV) ! u component of wind
REAL(KIND=JPRB),INTENT(IN)    :: PVP(KLON,KLEV) ! v component of wind
REAL(KIND=JPRB),INTENT(IN)    :: PCVL(KLON) ! low vegetation cover
REAL(KIND=JPRB),INTENT(IN)    :: PCVH(KLON) ! high vegetation cover
REAL(KIND=JPRB),INTENT(IN)    :: PGEMU(KLON) ! sine of latitude
REAL(KIND=JPRB),INTENT(IN)    :: PMFU(KLON,KLEV)  ! Conv. mass flux up
REAL(KIND=JPRB),INTENT(IN)    :: PSNOWACL(KLON,KLEV)  ! accretion rate of snow with cloud droplets
REAL(KIND=JPRB),INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)
REAL(KIND=JPRB),INTENT(INOUT) :: PAERDDP(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB),INTENT(INOUT) :: PAERSDM(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
! Total optical depth at various wavelenghts. 
! NOTE!! These wavelength definitions are not necessarily consistent
! with what is used in IFS-AER
REAL(KIND=JPRB),INTENT(OUT)   :: PODTO(KLON)

!REAL(KIND=JPRB),INTENT(OUT)   :: PODTO469(KLON), PODTO670(KLON), PODTO865(KLON), PODTO1240(KLON)
!REAL(KIND=JPRB),INTENT(IN)    :: PALBD(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NTSW), PFRTI(KLON,KTILES)

REAL(KIND=JPRB),INTENT(INOUT)   :: PAERO_WVL_DIAG(KLON,YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG,YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG_TYPES)
REAL(KIND=JPRB),INTENT(OUT)   :: PAERODDF(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO,8)
REAL(KIND=JPRB),INTENT(INOUT) :: PEXTRA(KLON,KLEVX,KFLDX)
REAL(KIND=JPRB),INTENT(INOUT) :: PAER_TAU(KLON,KLEV,14), PAER_SSA(KLON,KLEV,14),PAER_ASYM(KLON,KLEV,14)
REAL(KIND=JPRB),INTENT(INOUT) :: PAER_TAU_LW(KLON,KLEV,16)
REAL(KIND=JPRB),INTENT(OUT)   :: PTAUS_AER(KLON,KLEV,NBANDS_TROP,2),PTAUA_AER(KLON,KLEV,NBANDS_TROP,2)
REAL(KIND=JPRB),INTENT(OUT)   :: PPMAER(KLON,KLEV,NBANDS_TROP,2)
REAL(KIND=JPRB),INTENT(INOUT) :: PGFL(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM), PPRAERS(KLON)
! Simple sulfur scheme variables:
REAL(KIND=JPRB),INTENT(INOUT)   :: PSO2DD(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PBLH(KLON)  ! Boundary layer height
! PM output
REAL(KIND=JPRB),INTENT(OUT)     :: PAEPM1(KLON)
REAL(KIND=JPRB),INTENT(OUT)     :: PAEPM25(KLON)
REAL(KIND=JPRB),INTENT(OUT)     :: PAEPM10(KLON)

!*   0.5    LOCAL VARIABLES
!           ---------------

INTEGER(KIND=JPIM) :: JAER, JK, JL, JWAVL, JT, JB, JN
INTEGER(KIND=JPIM) :: JEXT, ITRC, IKLEVTROP(KLON), IW
INTEGER(KIND=JPIM) :: JO, JH, JY                         ! inside loop index for OIFS contex, HAM context and YAEROUT 
INTEGER(KIND=JPIM) :: JCLASS, JTILE, JMASS, JGAS, JCLOUD ! local loop indice for activation and dry deposition and tracer indexing
INTEGER(KIND=JPIM) :: ISSO2, ISSO4, ISSO4_ACS
INTEGER(KIND=JPIM) :: IMODE 
INTEGER(KIND=JPIM) :: IFLAG
INTEGER(KIND=JPIM) :: KTOP                               ! cloud top index

REAL(KIND=JPRB) :: ZAEROK(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZTAEROK(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZTAERO0(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZFAERO(KLON,ntrac)!YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZAER(KLON,KLEV), ZAERNEG(KLON,KLEV)
REAL(KIND=JPRB) :: ZAP(KLON,KLEV), ZAP_COV(KLON,KLEV)
REAL(KIND=JPRB) :: ZSO2(KLON,KLEV), ZDP(KLON,KLEV), ZDZ(KLON,KLEV) 
REAL(KIND=JPRB) :: ZITSO2(KLON,KLEV)
REAL(KIND=JPRB) :: ZFSO2(KLON)  , ZFSO4(KLON), ZFSO4_AQ(KLON)
REAL(KIND=JPRB) :: ZTSO2(KLON, KLEV)  , ZTSO4(KLON, KLEV,1), ZTSO4_AQ(KLON, KLEV)
REAL(KIND=JPRB) :: ZQSAT(KLON,KLEV), ZRHO(KLON,KLEV)
REAL(KIND=JPRB) :: ZTAER(KLON,KLEV)
REAL(KIND=JPRB) :: ZTAERO(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZRH(KLON,KLEV),ZTENC0(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NCHEM)!,ZTSO4(KLON,KLEV)
REAL(KIND=JPRB) :: ZM6RP(KLON,KLEV,NMOD)! NMOD=7
REAL(KIND=JPRB) :: ZM6DRY(KLON,KLEV,NSOL)
REAL(KIND=JPRB) :: ZWW(KLON,KLEV,NMOD)
REAL(KIND=JPRB) :: ZRHOP(KLON,KLEV,NMOD)
REAL(KIND=JPRB) :: ZSVOC(KLON,KLEV)
REAL(KIND=JPRB) :: ZELVOC(KLON,KLEV)
!THIS-IS-NEVER-USED   REAL(KIND=JPRB) :: ZSO4G(KLON,KLEV)
REAL(KIND=JPRB) :: ZCEN(KLON,KLEV,KTRAC) ! local tracer number and mixing ratios and gas concentrations for not tendency updated values
REAL(KIND=JPRB) :: PODTO469(KLON), PODTO670(KLON), PODTO865(KLON), PODTO1240(KLON)
REAL(KIND=JPRB) :: ZAER_TAU(KLON,KLEV,14,1), ZAER_SSA(KLON,KLEV,14),ZAER_ASYM(KLON,KLEV,14),ZAER_TAU_LW(KLON,KLEV,16)

! Optics output fields (to be used and allocated by methods using the optics)
REAL(KIND=JPRB), DIMENSION(:,:,:),   ALLOCATABLE :: ZTAUS_AER, ZTAUA_AER, ZPMAER ! extinctions
REAL(KIND=JPRB), DIMENSION(:,:,:,:), ALLOCATABLE :: ZAOP_OUT_EXT ! extinctions
REAL(KIND=JPRB), DIMENSION(:,:,:),   ALLOCATABLE :: ZAOP_OUT_A   ! single scattering albedo
REAL(KIND=JPRB), DIMENSION(:,:,:),   ALLOCATABLE :: ZAOP_OUT_G   ! assymetry factor

!Defined here, but should be passed in Call really:
TYPE(MODAL_DATA), DIMENSION(NMOD), TARGET :: RW_MODE
TYPE(MODAL_DATA), DIMENSION(NSOL), TARGET :: RWD_MODE
TYPE(MODAL_DATA), DIMENSION(NSOL), TARGET :: H2O_MODE
TYPE(MODAL_DATA), DIMENSION(NMOD), TARGET :: DENS_MODE

REAL(KIND=JPRB), ALLOCATABLE ::    ZAERNGT(:,:)

REAL(KIND=JPRB) :: ZDEGRAD, ZEPSCOV, ZEPSWAT, ZRWSAT, ZRWPWP
REAL(KIND=JPRB) :: ZQLWP(KLON,KLEV), ZQIWP(KLON,KLEV)
REAL(KIND=JPRB) :: ZTMPA

LOGICAL         :: LLIQCLD(KLON,KLEV) ! logical for liquid cloud
LOGICAL         :: LICECLD(KLON,KLEV) ! logical for ice cloud

REAL(KIND=JPRB), PARAMETER :: ZEPSEC=1e-14_JPRB
REAL(KIND=JPRB), PARAMETER :: ZMIN_CDNC=10.0_JPRB               !eehol: minimum CDNC (can be changed but for now 10 cm-3)
REAL(KIND=JPRB), PARAMETER :: ZDEF_CDNC=125.0_JPRB              !eehol: default CDNC (can be changed but for now 125 cm-3)
REAL(KIND=JPRB), PARAMETER :: ZDEF_RE_LIQ=4.0_JPRB              !eehol: default liq eff radius (can be changed but for now 4 mu m comes from liquid effective radius routine (PP_MIN_RE_UM))
REAL(KIND=JPRB), PARAMETER :: ZDEF_RE_ICE=80._JPRB*0.64952_JPRB !eehol: default ice eff radius (comes from ice effective radius routine (ZDEFAULT_RE_UM))
REAL(KIND=JPRB), PARAMETER :: ZTUNPAR=0.8164965_JPRB            !eehol: tuning parameter for sigma_w derived from TKE (square root of 2/3 (isotropy assumption))

REAL(KIND=JPRB),PARAMETER :: INFINITY=HUGE(1._JPRB)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! variables for the M7 call
REAL(KIND=JPRB) :: ZGRVOL(KLON,KLEV) !grid box volume for diagnostics
REAL(KIND=JPRB) :: ZPBL(KLON) !planetary boundary layer top level (in vdfmain.F90 ITOP=1)
REAL(KIND=JPRB) :: ZTMP(KLON) !temporary array to accumulate diagnostics
REAL(KIND=JPRB) :: ZXTM0(KLON,KLEV,ntrac) !tracer mixing ratios for HAM
REAL(KIND=JPRB) :: ZXTM1(KLON,KLEV,ntrac) !tracer mixing ratios for HAM
REAL(KIND=JPRB) :: ZXTTE(KLON,KLEV,ntrac) !tracer tendency for HAM
REAL(KIND=JPRB) :: ZXTTE_CLD(KLON,KLEV,ntrac) !tracer tendency for cloud vars for HAM
REAL(KIND=JPRB) :: ZXTTEM1(KLON,KLEV,ntrac) !tracer tendency for HAM
! added here variables for HAM-M7 activation
REAL(KIND=JPRB), ALLOCATABLE :: ZW(:,:,:) !mean or bins of updraft velocity [m s-1]
REAL(KIND=JPRB), ALLOCATABLE :: ZWPDF(:,:,:) !updraft velocity PDF over bins
REAL(KIND=JPRB), ALLOCATABLE :: ZRC(:,:,:,:) !critical radius of activation per mode [m]
REAL(KIND=JPRB), ALLOCATABLE :: ZSMAX(:,:,:) !maximum supersaturation
REAL(KIND=JPRB) :: ZTKEM1(KLON,KLEV) !turbulent kinetic energy  as zero as it is not used for now
REAL(KIND=JPRB) :: ZSIGMA_W(KLON,KLEV) !standard deviation of updraft pdf
REAL(KIND=JPRB) :: ZSMAXMN(KLON,KLEV) !maximum supersaturation for M&N scheme
REAL(KIND=JPRB) :: ZWCAPE(KLON) !CAPE  as zero as it is not used
REAL(KIND=JPRB) :: ZRDRY(KLON,KLEV,nclass) !dry radius for each class
REAL(KIND=JPRB) :: ZESW(KLON,KLEV) !saturation water vapor pressure
REAL(KIND=JPRB) :: ZA(KLON,KLEV,nclass) !curvature parameter A of the Koehler equation
REAL(KIND=JPRB) :: ZB(KLON,KLEV,nclass) !hygroscopicity parameter B of the Koehler equation
REAL(KIND=JPRB) :: ZSC(KLON,KLEV,nclass) !critical supersaturation [% 0-1]
REAL(KIND=JPRB) :: ZNACT(KLON,KLEV,nclass) !number of activated particles per mode [m-3]
REAL(KIND=JPRB) :: ZFRACN(KLON,KLEV,nclass) !fraction of activated particles per mode
REAL(KIND=JPRB) :: ZCDNCACT(KLON,KLEV)     !number of activated particles [m-3]
REAL(KIND=JPRB) :: ZCDNC_temp(KLON,KLEV)   !eehol: temporary for tendency of CDNC
REAL(KIND=JPRB) :: ZICNC_temp(KLON,KLEV)   !eehol: temporary for tendency of ICNC
REAL(KIND=JPRB) :: ZRE_LIQ(KLON,KLEV)! liquid effective radius
REAL(KIND=JPRB) :: ZNACT_AS(KLON,KLEV),ZNACT_CS(KLON,KLEV),ZNACT_KS(KLON,KLEV) ! variables for modewise activated fraction calculations
REAL(KIND=JPRB) :: ZFRAC_KS,ZFRAC_AS,ZFRAC_CS  ! variables for modewise activated fraction calculation
REAL(KIND=JPRB) :: ZNACT_TOT  ! variables for modewise activated fraction calculation
! variables for HAM-M7 wet deposition
REAL(KIND=JPRB) :: ZXTP1(KLON,KLEV,ntrac)  !updated tracer mass/number mixing ratio
REAL(KIND=JPRB) :: ZXTP1C(KLON,KLEV,ntrac) !in-cloud tracer mass/number mixing ratio
REAL(KIND=JPRB) :: ZXTP10(KLON,KLEV,ntrac) !ambient tracer mass/number mixing ratio
REAL(KIND=JPRB) :: ZDUMMY(KLON,ntrac) !placeholder for pxtbound which is only necessary in the conv. case (conv. massfix boundary condition)
REAL(KIND=JPRB) :: ZDUM3D(KLON,KLEV,ntrac) !updraft mass flux (for conv case only)

REAL(KIND=JPRB) :: ZWDEP_SCAV_IC(KLON,ntrac) !in-cloud scavenged mr
REAL(KIND=JPRB) :: ZWDEP_SCAV_BC(KLON,ntrac) !below cloud scavenged mr
REAL(KIND=JPRB) :: ZFUXT3D(KLON,KLEV,ntrac) !updraft mass flux (for conv case only)

REAL(KIND=JPRB) :: ZDUM2D(KLON,KLEV) !convective flux needed only for conv. case (see cuflx)
REAL(KIND=JPRB) :: ZLFRAC_SO2(KLON,KLEV) !liquid tracer fraction (SO2) -ham specific-
REAL(KIND=JPRB) :: ZDPG(KLON,KLEV) !dp/g
REAL(KIND=JPRB) :: ZQP(KLON,KLEV) !full level humidity with threshold

LOGICAL :: LSTRAT  !logical switch for stratiform or convective case (TRUE for strat., FALSE for conv.)
REAL(KIND=JPRB) :: ZFEVAPR_cov(KLON,KLEV) !evaporation of rain, convective case [kg/m2/s]
REAL(KIND=JPRB) :: ZFSUBLS_cov(KLON,KLEV) !sublimation of snow, convective case [kg/m2/s]
REAL(KIND=JPRB) :: ZFEVAPR_str(KLON,KLEV) !evaporation of rain, stratiform case [kg/m2/s]
REAL(KIND=JPRB) :: ZFSUBLS_str(KLON,KLEV) !sublimation of snow, stratiform case [kg/m2/s]
REAL(KIND=JPRB) :: ZMRATEPR_cov(KLON,KLEV) !rain formation rate in cloudy part, convective case
REAL(KIND=JPRB) :: ZMRATEPS_cov(KLON,KLEV) ! ice formation rate in cloudy part, convective case
REAL(KIND=JPRB) :: ZMRATEPR_str(KLON,KLEV) !rain formation rate in cloudy part, stratiform case
REAL(KIND=JPRB) :: ZMRATEPS_str(KLON,KLEV) ! ice formation rate in cloudy part, stratiform case

REAL(KIND=JPRB) :: ZMSNOWACL_str(KLON,KLEV) !accretion rate of snow with cloud droplets in cloudy part strat
REAL(KIND=JPRB) :: ZMSNOWACL_cov(KLON,KLEV) !accretion rate of snow with cloud droplets in cloudy part conv
REAL(KIND=JPRB) :: ZFLXR, ZFLXS, ZFLXRB, ZFLXSB !variables to calculate rain and snow evap/formation

REAL(KIND=JPRB) :: ZFPLCL(KLON,KLEV),ZFPLCN(KLON,KLEV),ZFPLSL(KLON,KLEV),ZFPLSN(KLON,KLEV)
REAL(KIND=JPRB) :: ZLP(KLON,KLEV) !temporary variable for cloud water content
REAL(KIND=JPRB) :: ZIP(KLON,KLEV) !temporary variable for cloud ice water content
REAL(KIND=JPRB) :: ZLPU(KLON,KLEV) !temporary variable for cloud water content
REAL(KIND=JPRB) :: ZIPDUM(KLON,KLEV) !temporary variable for cloud ice water content
REAL(KIND=JPRB) :: ZMFU(KLON,KLEV)  ! temporary variable for constraing Conv. mass flux up
! variables for ICNC calculations
REAL(KIND=JPRB) :: ZICNC(KLON,KLEV) ! ice crystal number concentration [#/cm3]
! added here variables for dry deposition and sedimentation
REAL(KIND=JPRB) :: ZTENCIH(KLON,KLEV,ntrac) !for HAM tendencies
REAL(KIND=JPRB) :: ZAHFSM(KLON)
REAL(KIND=JPRB) :: ZWND(KLON)
! variables not needed for HAM aerosol dry deposition (if gas deposition and different surfaces are taken into account these need to be revised!)
REAL(KIND=JPRB) :: ZCFML(KLON), ZCFMW(KLON), ZCFMI(KLON), ZCFNCL(KLON), ZCFNCW(KLON), ZCFNCI(KLON), ZEPDU2, ZKAP
REAL(KIND=JPRB) :: ZGEOM1(KLON,KLEV), ZRIL(KLON), ZRIW(KLON), ZRII(KLON), ZTVIR1(KLON,KLEV), ZTVL(KLON), ZTVW(KLON)
REAL(KIND=JPRB) :: ZTVI(KLON), ZAZ0(KLON), ZFRL(KLON), ZSRFL(KLON), ZFOREST(KLON), ZTSI(KLON), ZAZ0L(KLON)
REAL(KIND=JPRB) :: ZAZ0I(KLON), ZCDNI(KLON)
! variables needed for HAM aerosol dry deposition
LOGICAL         :: ZLOLAND(KLON) !logical land mask
REAL(KIND=JPRB) :: ZXTEMS(KLON,ntrac) !surface emissions modified by dry deposition
REAL(KIND=JPRB) :: ZAZ0W(KLON), ZFRW(KLON), ZCVS(KLON), ZCVW(KLON), ZVGRAT(KLON) !rough. len. wat., wat. frac., snow cov. frac., wet skin frac., veg. ratio
REAL(KIND=JPRB) :: ZCDNL(KLON), ZCDNW(KLON) !ustar (in not used variable), aerodynamic resis. on surface (in not used variable)
REAL(KIND=JPRB) :: ZXTMD1(KLON,KLEV,NTRAC) !tracer mixing ratios for HAM drydep (updated with tend)
! output diagnostics
INTEGER, PARAMETER :: N_NUC_DIAG=5
REAL(KIND=JPRB) :: ZOUT3(KLON,KLEV,2*(NAEROCOMP+NCLASS)), ZOUT_DNUC(KLON,KLEV,N_NUC_DIAG)  
REAL(KIND=JPRB) :: SEDOUT(KLON,KLEV,KTRAC)   ! changed ntrack to ktrac (RCHG)
REAL(KIND=JPRB) :: DDEPOUT(KLON,KLEV,KTRAC)
REAL(KIND=JPRB) :: WDEPOUT(KLON,KLEV,KTRAC)
REAL(KIND=JPRB) :: SEDOUT_2D(KLON,KTRAC)

REAL(KIND=JPRB) :: WDEPOUT_2D(KLON,KTRAC)
REAL(KIND=JPRB) :: WDEPOUT_IC_2D(KLON,KTRAC)
REAL(KIND=JPRB) :: WDEPOUT_BC_2D(KLON,KTRAC)

REAL(KIND=JPRB) :: ZSEDIFLUX(KLON,KLEV,NTRAC)
REAL(KIND=JPRB) :: ZSEDIFLUXSURF(KLON,NTRAC)  
REAL(KIND=JPRB) :: ZDDEPFLUX(KLON,NTRAC)
REAL(KIND=JPRB) :: ZDDEPFLUX_SO2(KLON)
REAL(KIND=JPRB) :: ZVDEP(KLON,NTRAC) !ddep velocity for diagnostics from ham
INTEGER(KIND=JPIM), parameter::ZKROW=1 ! KROW only used in ECHAM but needed inside HAM-codes so set as 1.
INTEGER(KIND=JPIM) :: IBLK
REAL(KIND=JPRB)    :: REFFI(KLON,KLEV,ZKROW), REFFL(KLON,KLEV,ZKROW)
INTEGER(KIND=JPIM) :: LWBANDS ! number of LW bands
REAL(KIND=JPRB)    :: PRS1D(KLON,KLEV)
INTEGER(KIND=JPIM) :: ISO4_C, ISSO4_C ! temporary tracer index of gas-phase SO4 (retrieved from chemistry module)

REAL(KIND=JPRB) :: PAOD(KLON,NASWBAND), PSSA(KLON,NASWBAND), PABS(KLON,NASWBAND), PASY(KLON,NASWBAND),PFAOD(KLON,NASWBAND)
REAL(KIND=JPRB) :: PAOD_LW(KLON,16)

! Boundary layer height index calculation
REAL(KIND=JPRB) :: ZBLHIDX(KLON)   ! index
LOGICAL         :: LBLHFOUND(KLON) ! logical if boundary layer height is found  
REAL(KIND=JPRB) :: ZRG             ! 1/RG

!!! parameters needed for diagnostic aerosol optical properties
LOGICAL, PARAMETER :: LSAFETY=.FALSE. ! If true, impose limits on computed optical properties
LOGICAL         :: LDIAG_AEROPT ! logical for aerosol optics
REAL(KIND=JPRB) :: ZAER_TAU_DIAG(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG)
REAL(KIND=JPRB) :: ZAER_SSA_DIAG(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG)
REAL(KIND=JPRB) :: ZAER_ASYM_DIAG(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG)
REAL(KIND=JPRB) :: ZAER_LidarR_DIAG(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG)
REAL(KIND=JPRB) :: LAMBDA_DIAG(YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG)

REAL(KIND=JPRB) :: ZAOD_DIAG(KLON,YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG), ZSSA_DIAG(KLON,YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG)
REAL(KIND=JPRB) :: ZABS_DIAG(KLON,YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG), ZASY_DIAG(KLON,YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG)
!-----------------------------------------------------------------------

#include "abor1.intfb.h"
#include "surf_inq.h"
#include "satur.intfb.h"
#include "aer_negat.intfb.h"
#include "tm5m7_optics_aop_get.intfb.h"
#include "troplev.intfb.h"
#include "chem_inext.intfb.h"
#include "ice_effective_radius.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('HAMM7_INTERFACE',0,ZHOOK_HANDLE)
ASSOCIATE( &
         & YGFL         => YDMODEL%YRML_GCONF%YGFL,        & 
         & YDPHYRAD     => YDMODEL%YRML_PHY_RAD,           &
         & YDPHYAER     => YDMODEL%YRML_PHY_AER,           &
         & YREAERSRC    => YDMODEL%YRML_PHY_AER%YREAERSRC, &
         & YREAERATM    => YDMODEL%YRML_PHY_RAD%YREAERATM, &
         & YREAERLID    => YDMODEL%YRML_PHY_AER%YREAERLID, &
         & YREAERSNK    => YDMODEL%YRML_PHY_AER%YREAERSNK, &
         & YRERAD       => YDMODEL%YRML_PHY_RAD%YRERAD,    &
         & YDRIP        => YDMODEL%YRML_GCONF%YRRIP,       &
         & YDCHEM       => YDMODEL%YRML_CHEM%YRCHEM,       &
         & YRCOMPO      => YDMODEL%YRML_CHEM%YRCOMPO,      &
         & YREPHY       => YDMODEL%YRML_PHY_EC%YREPHY,     &
         & YRECLDP      => YDMODEL%YRML_PHY_EC%YRECLDP,    &
         & YDSPP_CONFIG => YDMODEL%YRML_GCONF%YRSPP_CONFIG)

 ASSOCIATE( &
         ! --- YGFL ------------------------------------------------
         & NACTAERO  => YGFL%NACTAERO, NAERO    => YGFL%NAERO,     &
         & NGHG      => YGFL%NGHG,     NCHEM    => YGFL%NCHEM,     &
         & NDIM      => YGFL%NDIM,                                 &
         & YAEROUT   => YGFL%YAEROUT,  YR       => YGFL%YR,        &
         & YS        => YGFL%YS,       YCHEM    => YGFL%YCHEM,     &
         & YAERO     => YGFL%YAERO,    YGHG     => YGFL%YGHG,      &
         & YAEROCLIM => YGFL%YAEROCLIM,                            &
         & YCDNC     => YGFL%YCDNC,    YICNC    => YGFL%YICNC,     &
         & YRE_LIQ   => YGFL%YRE_LIQ,  YRE_ICE  => YGFL%YRE_ICE,   &
         !included CDNC, ICNC, liq and ice eff rad
         & NAERO_WVL_DIAG => YGFL%NAERO_WVL_DIAG,                  &
         & LAERCHEM       => YGFL%LAERCHEM,                        &
         ! --- YRCOMPO ---------------------------------------------
         & LCHEM_DIA      => YRCOMPO%LCHEM_DIA,                    &
         & AERO_SCHEME    => YRCOMPO%AERO_SCHEME,                  &
         & CHEM_SCHEME    => YDCHEM%CHEM_SCHEME,                   &
         & LAERNITRATE    => YRCOMPO%LAERNITRATE,                  &
         & NSO4SCHEME     => YREAERSRC%NSO4SCHEME,                 &
         ! --- YREAERATM -------------------------------------------
         & LAERDRYDP      => YREAERATM%LAERDRYDP,                  &
         & LAERSEDIM      => YREAERATM%LAERSEDIM,                  &
         & LAERSURF       => YREAERATM%LAERSURF,                   & !add logicals for dry dep and sedi
         & LAER6SDIA      => YREAERATM%LAER6SDIA,                  &
         & LAERCLIMG      => YREAERATM%LAERCLIMG,                  &
         & LAERCLIMZ      => YREAERATM%LAERCLIMZ,                  &
         & LAERGTOP       => YREAERATM%LAERGTOP,                   &
         & LAERHYGRO      => YREAERATM%LAERHYGRO,                  &
         & LAERLISI       => YREAERATM%LAERLISI,                   &
         & LAERNGAT       => YREAERATM%LAERNGAT,                   &
         & LAERSCAV       => YREAERATM%LAERSCAV,                   &
         & LAERSCAV_CHEM  => YREAERATM%LAERSCAV_CHEM,              &
         & LAERVOL        => YREAERATM%LAERVOL,                    &
         & NXT3DAER       => YREAERATM%NXT3DAER,                   & 
         & LAERRRTM       => YREAERATM%LAERRRTM,                   &       
         & REPSCAER       => YREAERATM%REPSCAER,                   &         
         ! --- YRERAD ----------------------------------------------
         & LAERVISI       => YRERAD%LAERVISI,                      &
         & NTSW           => YRERAD%NTSW,                          &
         & RNS            => YRERAD%RNS,                           &
         & RSIGAIR        => YRERAD%RSIGAIR,                       &
         & NRADFR         => YRERAD%NRADFR,                        & !FREQUENCY OF FULL RADIATION COMPUTATIONS
         & NAEROOPT       => YRERAD%NAEROOPT,                      &
         & NCLOUDACT      => YRERAD%NCLOUDACT,                     & ! integer to switch activation scheme (0=default,1=Morales&Nenes, 2=Abdul-Razzak&Ghan)
         & RCCNSEA        => YRERAD%RCCNSEA,                       & ! default ccn value over sea
         & RCCNLND        => YRERAD%RCCNLND,                       & ! default ccn value over land
         ! --- OTHERS ----------------------------------------------
         & YDAERM7        => YDPHYAER%YREAEROPT,                   & ! use this to transfer AOD, SSA and ASY to rad scheme
         & NWLID          => YREAERLID%NWLID,                      &
         & YSURF          => YREPHY%YSURF,                         &
         & RRHTAB         => YREAERSNK%RRHTAB,                     &
         & RNICE          => YRECLDP%RNICE,                        & !default for ICNC
         & RCLDMAX        => YRECLDP%RCLDMAX,                      & !max cloud value
         & LAERICESED     => YRECLDP%LAERICESED,                   &
         & LAERICEAUTO    => YRECLDP%LAERICEAUTO,                  &
         & NSTART         => YDRIP%NSTART                          ) 

! & NINDSCAV=>YREAERATM%NINDSCAV, NTSCAV=>YREAERATM%NTSCAV, &
! & NDDUST=>YREAERSRC%NDDUST, NTYPAER=>YREAERSRC%NTYPAER, &
!     ------------------------------------------------------------------

!*         0.     PROGNOSTIC AEROSOLS - FINAL COMPUTATIONS
!                 ----------------------------------------


LBLHFOUND(:) = .FALSE.

ZAHFSM  = 0._JPRB
ZEPSCOV = 1.E-03_JPRB
ZEPSWAT = 1.E-18_JPRB
ZDEGRAD = 180._JPRB/RPI

CALL SURF_INQ(YSURF,PRWSAT=ZRWSAT)
CALL SURF_INQ(YSURF,PRWPWP=ZRWPWP)

!*         0.1    SWITCHING ON POSSIBLE DEBUG PRINTS AND ALLOCATING MEMORY
!                 --------------------------------------------------------

ALLOCATE( ZAERNGT(KLON,NACTAERO) )
ZAERNGT(KIDIA:KFDIA,1:NACTAERO) = 0._JPRB

ZOUT3(KIDIA:KFDIA,:,:) = 0._JPRB
ZOUT_dnuc(KIDIA:KFDIA,:,:) = 0._JPRB
! Need to initialize those 3 arrays early in case LAERDRYDP=F (GNU, Lianghai Wu)
ZVDEP(KIDIA:KFDIA,:) = 0._JPRB ! ddep velocity as zero
ZXTEMS(KIDIA:KFDIA,:) = 0._JPRB ! surface emissions as zero for input
ZXTMD1(KIDIA:KFDIA,:,:) = 0._JPRB 

SEDOUT(KIDIA:KFDIA,:,:)      = 0._JPRB
DDEPOUT(KIDIA:KFDIA,:,:)     = 0._JPRB
WDEPOUT(KIDIA:KFDIA,:,:)     = 0._JPRB
ZSEDIFLUX(KIDIA:KFDIA,:,:)   = 0._JPRB
ZSEDIFLUXSURF(KIDIA:KFDIA,:) = 0._JPRB
ZDDEPFLUX(KIDIA:KFDIA,:)     = 0._JPRB
ZDDEPFLUX_SO2(KIDIA:KFDIA)   = 0._JPRB
SEDOUT_2D(KIDIA:KFDIA,:)     = 0._JPRB
WDEPOUT_2D(KIDIA:KFDIA,:)    = 0._JPRB
WDEPOUT_IC_2D(KIDIA:KFDIA,:) = 0._JPRB
WDEPOUT_BC_2D(KIDIA:KFDIA,:) = 0._JPRB
ZMRATEPR_cov(KIDIA:KFDIA,:) = 0._JPRB
ZMRATEPS_cov(KIDIA:KFDIA,:) = 0._JPRB
ZFEVAPR_cov(KIDIA:KFDIA,:) = 0._JPRB ! 
ZFSUBLS_cov(KIDIA:KFDIA,:) = 0._JPRB
ZTAERO(KIDIA:KFDIA,:,:)    = 0._JPRB

ZCDNCACT(KIDIA:KFDIA,:) = 0._JPRB     !number of activated particles [m-3]
ZCEN(KIDIA:KFDIA,:,:)   = 0._JPRB

ZFRACN(KIDIA:KFDIA,:,:) = 0._JPRB !fraction of activated particles per mode

ZRG=1/RG

! computation of tropopause level 
CALL TROPLEV(KLON,KIDIA,KFDIA,KLEV,.FALSE.,PTP,PQP,PRSF1,IKLEVTROP)

! Initializing tracer number and mixing ratios and gas concentrations to not be tendency updated values
ZCEN(KIDIA:KFDIA,1:KLEV,:) = 0._JPRB
ITRC=0
DO JEXT=1,NGHG
  ITRC=ITRC+1
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZCEN(JL,JK,ITRC) = PGFL(JL,JK,YGHG(JEXT)%MP9_PH)
    ENDDO
  ENDDO
ENDDO

DO JEXT=1,NAERO
  ITRC=ITRC+1
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZCEN(JL,JK,ITRC) = PGFL(JL,JK,YAERO(JEXT)%MP9_PH)
    ENDDO
  ENDDO
ENDDO

DO JEXT=1,NCHEM
  ITRC=ITRC+1
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZCEN(JL,JK,ITRC) = PGFL(JL,JK,YCHEM(JEXT)%MP9_PH)
    ENDDO
  ENDDO
ENDDO
  
!DO JAER=1,NACTAERO
!  DO JK=1,KLEV
!    DO JL=KIDIA,KFDIA
!      ! check on max-values here to prevent excessive values in OPTICS_AOP_GET.. probably to be removed again!!
!      ZAEROK(JL,JK,JAER) =MIN(ZCEN (JL,JK,KAERO(JAER)), 1e10)
!      ZTAEROK(JL,JK,JAER)=PTENC(JL,JK,KAERO(JAER))
!    ENDDO
!  ENDDO
!ENDDO

! Initialize output aerosol tendencies
DO JAER=1,NACTAERO
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZAEROK(JL,JK,JAER) =ZCEN (JL,JK,KAERO(JAER))
      ZTAEROK(JL,JK,JAER)=PTENC(JL,JK,KAERO(JAER))
    ENDDO
  ENDDO
ENDDO

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZRHO(JL,JK) = PRSF1(JL,JK)/(RD*PTP(JL,JK))
    ZDP(JL,JK)  = PRS1(JL,JK) - PRS1(JL,JK-1)
    ZDPG(JL,JK) = (PRS1(JL,JK) - PRS1(JL,JK-1))/RG ! calculate dp/g
    ZDZ(JL,JK)  = ZDP(JL,JK) / (ZRHO(JL,JK)*RG)
  ENDDO
ENDDO

ZM6RP(KIDIA:KFDIA,:,:)  = 0.0_JPRB
ZM6DRY(KIDIA:KFDIA,:,:) = 0.0_JPRB
ZRHOP(KIDIA:KFDIA,:,:)  = 0.0_JPRB
ZWW(KIDIA:KFDIA,:,:)    = 0.0_JPRB

DO JB=1, NBANDS_TROP
  PTAUS_AER(KIDIA:KFDIA,1:KLEV,JB,1) = 0.0_JPRB
  PTAUA_AER(KIDIA:KFDIA,1:KLEV,JB,1) = 0.0_JPRB
  PPMAER   (KIDIA:KFDIA,1:KLEV,JB,1) = 0.0_JPRB
  PTAUS_AER(KIDIA:KFDIA,1:KLEV,JB,2) = 0.0_JPRB
  PTAUA_AER(KIDIA:KFDIA,1:KLEV,JB,2) = 0.0_JPRB
  PPMAER   (KIDIA:KFDIA,1:KLEV,JB,2) = 0.0_JPRB
ENDDO

IF (LCHEM_DIA) THEN
  ZTAERO0(KIDIA:KFDIA,1:KLEV,1:NACTAERO) =  ZTAEROK(KIDIA:KFDIA,1:KLEV,1:NACTAERO)
  ZTENC0(KIDIA:KFDIA,1:KLEV, :) = 0._JPRB
ENDIF 

!
!*         1.1    COMPUTE RELATIVE HUMIDITY WITHOUT VERTICAL SMOOTING
!                 ---------------------------------------------------
! Q at saturation for RH calculation
IFLAG=2
CALL SATUR(KIDIA , KFDIA , KLON  , KTDIA , KLEV,  YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN, &
  & PRSF1, PTP    , ZQSAT , IFLAG)  

! RH calculation
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZQP(JL,JK)=MAX(0.0_JPRB,PQP(JL,JK)) ! threshold for full level humid
    ZRH(JL,JK)=ZQP(JL,JK)/(MAX(1.E-30_JPRB,ZQSAT(JL,JK)))
    ZRH(JL,JK)=MIN(1.0_JPRB,MAX(0.0_JPRB,ZRH(JL,JK)))
    ZAP(JL,JK)=MIN(1.0_JPRB,MAX(0.0_JPRB,PAP(JL,JK))) ! threshold for cloud cover
    IF (JK >= KCTOP(JL) .AND. JK <= KCBOT(JL)) THEN ! cloud fraction for convective wet removal (0 outside and 1 inside cloud)
      ZAP_COV(JL,JK) = 1.0_JPRB
    ELSE
      ZAP_COV(JL,JK) = 0.0_JPRB
    END IF
  ENDDO
ENDDO
ZBLHIDX(KIDIA:KFDIA)=1.0_JPRB
! Find top level index of bounrary layer
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ! check the PBL height against the grid point height from surface
    ! PGEOH used here as it is the layer interface geopotential, and when divided by gravitation constat gives height
    ! 
    IF (PBLH(JL)>(((PGEOH(JL,JK)-PGEOH(JL,KLEV))*ZRG)) .and. .not. LBLHFOUND(JL)) THEN
      ZBLHIDX(JL)=REAL(JK, KIND=JPRB)
      LBLHFOUND(JL)=.TRUE.
    END IF
  ENDDO 
ENDDO

!THIS-IS-NEVER-USED   ! TB apparently unnecessary in current implementation, but ISSO4_C still needed for chem_inext in the code.
!THIS-IS-NEVER-USED   ! needs to be reviewed if it can be removed. 
!THIS-IS-NEVER-USED   IF (TRIM(CHEM_SCHEME)=="tm5") THEN
!THIS-IS-NEVER-USED     DO JT=1,NCHEM
!THIS-IS-NEVER-USED       IF(TRIM(YCHEM(JT)%CNAME)== 'SO4' ) THEN
!THIS-IS-NEVER-USED         ISSO4_C=KCHEM(JT)
!THIS-IS-NEVER-USED         ISO4_C=JT
!THIS-IS-NEVER-USED       ENDIF
!THIS-IS-NEVER-USED     ENDDO
!THIS-IS-NEVER-USED     ZSO4G(KIDIA:KFDIA,1:KLEV)=ZCEN(KIDIA:KFDIA,1:KLEV,ISSO4_C)
!THIS-IS-NEVER-USED   ELSE IF (TRIM(CHEM_SCHEME)=="SimChem") THEN
!THIS-IS-NEVER-USED     ZSO4G(KIDIA:KFDIA,1:KLEV)=0._JPRB
!THIS-IS-NEVER-USED   ELSE
!THIS-IS-NEVER-USED     CALL ABOR1(" M7: UNCOUPLED CHEMISTRY SCHEME "//TRIM(CHEM_SCHEME) ) 
!THIS-IS-NEVER-USED   ENDIF


! Convert from kg kg-1 to molec cm-3

SELECT CASE (TRIM(CHEM_SCHEME))
CASE("tm5")
  ZELVOC(KIDIA:KFDIA,1:KLEV)= ZCEN(KIDIA:KFDIA,1:KLEV,KAERO(ielvoc))
  ZSVOC(KIDIA:KFDIA,1:KLEV) = ZCEN(KIDIA:KFDIA,1:KLEV,KAERO(iisvoc))
CASE("SimChem")
  ZELVOC(KIDIA:KFDIA,1:KLEV)= 0.0_JPRB
  ZSVOC(KIDIA:KFDIA,1:KLEV) = 0.0_JPRB
CASE DEFAULT
  ! This should be caught earlier at setup
  CALL ABOR1(" M7: UNCOUPLED CHEMISTRY SCHEME "//TRIM(CHEM_SCHEME) )
END SELECT

!calculate ICNC
ZICNC(KIDIA:KFDIA,1:KLEV) = 0._JPRB
ZICNC(KIDIA:KFDIA,1:KLEV) = RNICE

! Initialize "from-OIFS-to-HAM" tracer mixing ratios and tendencies
ZXTM0(KIDIA:KFDIA,1:KLEV,:) = 0._JPRB
ZXTM1(KIDIA:KFDIA,1:KLEV,:) = 0._JPRB
ZXTTE(KIDIA:KFDIA,1:KLEV,:) = 0._JPRB
ZXTTE_CLD(KIDIA:KFDIA,1:KLEV,:) = 0._JPRB
ZXTTEM1(KIDIA:KFDIA,1:KLEV,:) = 0._JPRB

!number
DO JCLASS=1,nclass
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZXTM1(JL,JK,ind_oifs_ham%ind_class_HAM(JCLASS)) = ZCEN(JL,JK,KAERO(ind_oifs_ham%ind_class_OIFS(JCLASS)))
      ZXTTE(JL,JK,ind_oifs_ham%ind_class_HAM(JCLASS)) = PTENC(JL,JK,KAERO(ind_oifs_ham%ind_class_OIFS(JCLASS)))
    END DO
  END DO
END DO

!mass
DO JMASS=1,naerocomp
  JO=ind_oifs_ham%ind_mass_OIFS(JMASS) ! JO -> index context OIFS
  JH=ind_oifs_ham%ind_mass_HAM(JMASS)  ! JH -> index context HAM
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZXTM1(JL,JK,JH) = ZCEN(JL,JK,KAERO(JO))
      ZXTTE(JL,JK,JH) = PTENC(JL,JK,KAERO(JO))
      ! in case of simple sulfur scheme add SO4_AQ part into SO4_ACS
      ! both original tendency and m7tendency [FIXME: what??]
      
      !ADD SO4 from wet chemistry to tendencies
      if(trim(YAERO(JO)%CNAME)=='SO4_AS') then   
        ZXTTE(JL,JK,JH)=ZXTTE(JL,JK,JH)+PCHEM2AER(JL,JK,2)!!! need to be verrified, Lianghai
      end if
      !if(trim(YAERO(ind_oifs_ham%ind_mass_OIFS(JMASS))%CNAME)=='SO4') then!!! add SO4 into tendency, ugly loop for now,Lianghai
      !  ZXTTE(JL,JK,ind_oifs_ham%ind_mass_HAM(JMASS))=ZXTTE(JL,JK,ind_oifs_ham%ind_mass_HAM(JMASS))+PCHEM2AER(JL,JK,1)
      !end if

    END DO
  END DO
END DO

!gas
DO JGAS=1,subm_ngasspec
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      JO = ind_oifs_ham%ind_gas_OIFS(JGAS) ! JO -> index context OIFS
      JH = ind_oifs_ham%ind_gas_HAM(JGAS)  ! JH -> index context HAM
      
      ZXTM1(JL,JK,JH) = ZCEN(JL,JK,KCHEM(JO))
      ZXTTE(JL,JK,JH) = PTENC(JL,JK,KCHEM(JO))
      
      !IF(TRIM(YCHEM(JO)%CNAME)=='H2SO4')THEN
      !   ZXTM1(JL,JK,JH) = ZCEN(JL,JK,KCHEM(JO))
      !   !ZSO4_PROD(JL,JK) = PCHEM2AER(JL,JK,1)*time_step_len!PTENC(JL,JK,KCHEM(ind_gas_OIFS(JGAS)))*time_step_len
      !   !ZSO4_PROD_TEST(JL,JK) = PTENC(JL,JK,KCHEM(JO))*time_step_len
      !ELSE
      !   ZXTTE(JL,JK,JH) = PTENC(JL,JK,KCHEM(JO))
      !END IF
      
      !add for diagnostics
      !ZSO4_PROD(JL,JK) = PTENC(JL,JK,KCHEM(JO))*time_step_len
      !ZSO4_CONC(JL,JK) = ZCEN(JL,JK,KCHEM(JO)) 
      !         ELSE !simple sulfur scheme
      !            IF(TRIM(YAERO(JO)%CNAME)=='SO2')THEN
      !               ZXTM1(JL,JK,JH) = ZCEN(JL,JK,KAERO(JO))
      !               ! 
      !               ZXTTE(JL,JK,JH) = PTENC(JL,JK,KAERO(JO))
      !               ZXTTEM1(JL,JK,JH) = PTENC(JL,JK,KAERO(JO))
      !
      !            ELSE IF (TRIM(YAERO(JO)%CNAME)=='SO4_gas')THEN
      !               ZXTM1(JL,JK,JH) = ZCEN(JL,JK,KAERO(JO))
      !               ZXTTE(JL,JK,JH) = PTENC(JL,JK,KAERO(JO)) 
      !               ZXTTEM1(JL,JK,JH) = PTENC(JL,JK,KAERO(JO)) 
      !            END IF
    END DO
  END DO
END DO

! RCHG -> This will produce segmentation fault if CDNC are not in the namelist 
!         we need to test these things and do a CALL ABORT1() 
!         for these cases. For this we need to find the best flags or combinations 
!         of flags. FIXME

! IF (LAERRRTM) THEN
!cloud
DO JCLOUD=1,2 !CDNC and ICNC
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZXTM1(JL,JK,ind_oifs_ham%ind_cloud_HAM(JCLOUD)) = ZCEN(JL,JK,KAERO(ind_oifs_ham%ind_cloud_OIFS(JCLOUD)))
      ZXTTE(JL,JK,ind_oifs_ham%ind_cloud_HAM(JCLOUD)) = PTENC(JL,JK,KAERO(ind_oifs_ham%ind_cloud_OIFS(JCLOUD)))
    END DO
  END DO
END DO
!ENDIF

!eehol: init CDNC and ICNC with non updated to temporary HAMM7 variables (unit in #/kg)
ZCDNC_temp(KIDIA:KFDIA,1:KLEV) = (1.0E6_JPRB*PGFL(KIDIA:KFDIA,1:KLEV,YCDNC%MP9_PH))/ZRHO(KIDIA:KFDIA,1:KLEV)
ZICNC_temp(KIDIA:KFDIA,1:KLEV) = (1.0E6_JPRB*PGFL(KIDIA:KFDIA,1:KLEV,YICNC%MP9_PH))/ZRHO(KIDIA:KFDIA,1:KLEV)

! implementation of HAM-M7
ZWND(KIDIA:KFDIA) = 0._JPRB
DO JL=KIDIA,KFDIA
  ZWND(JL)=MAX(1.0E-10_JPRB,PAERWS(JL)) ! make threshold for wind speed
  ! Compute average fluxes over tiles
  DO JTILE=1,KTILES
    ZAHFSM(JL)=ZAHFSM(JL)+PFRTI(JL,JTILE)*PAHFSTI(JL,JTILE)
  ENDDO
ENDDO
 
! --> calling the microphysics scheme
 
    ! Initializations for submodel interface
    ZGRVOL(KIDIA:KFDIA,1:KLEV) = 1.79e12_JPRB ! ZGRVOL is only used for diagnostics (only when HAMMOZ is on)
    ZPBL = 1 ! boundary layer top = 1 (ITOP=1)
    ! Allocate variables for aerosol processes
    ALLOCATE(ZW(KLON,KLEV,nw))
    ALLOCATE(ZWPDF(KLON,KLEV,nw))
    !IF (.NOT. ALLOCATED(w_large)) ALLOCATE(w_large(KLON,KLEV,ZKROW))
    !IF (.NOT. ALLOCATED(w_turb)) ALLOCATE(w_turb(KLON,KLEV,ZKROW))
    ALLOCATE(ZRC(KLON,KLEV,nclass,nw))
    ALLOCATE(ZSMAX(KLON,KLEV,nw))
    !IF (.NOT. ALLOCATED(reffi)) ALLOCATE(reffi(KLON,KLEV,ZKROW))
    !IF (.NOT. ALLOCATED(reffl)) ALLOCATE(reffl(KLON,KLEV,ZKROW))
    DO IMODE=1,NMOD
      ALLOCATE(RW_MODE (IMODE)%d2(KLON,KLEV))
      ALLOCATE(DENS_MODE(IMODE)%d2(KLON,KLEV))
      IF (sizeclass(IMODE)%lsoluble) THEN
        ALLOCATE(RWD_MODE(IMODE)%d2(KLON,KLEV))
        ALLOCATE(H2O_MODE(IMODE)%d2(KLON,KLEV))
      END IF
    ENDDO
    ! End allocate variables for aerosol processes

    ZWPDF(KIDIA:KFDIA,1:KLEV,:) = 0.0_JPRB
    
    !-----------------------------------------------------------------
    ! Submodel interface call (HAM aerosol microphysics)
    
    CALL GSTATS(2501,0)

    CALL HAM_SUBM_INTERFACE(&
         & KFDIA, KLON,  KLEV, ZKROW, & !dimension indices
         & ntrac, PRSF1, PRS1,        & !number of tracers, pressure full levels, pressure half levels
         & PTP,   ZQP,   ZQSAT,       & !temperature, specific humidity, saturation spec. hum.
         & ZXTM1, ZXTTE,              & !tracer mass/number mr, tracer tendencies
         & ZM6RP, ZM6DRY,             & !mean mode actual radius [m], dry radius for soluble modes [m] 
         & ZRHOP, ZWW,                & !mean mode particle density [kg m-3], aerosol water content for each mode [kg(water)-3(air)]
         & ZAP,   ZGRVOL,             & !cloud fraction, grid box volume (only for diagnostics)
         & ZBLHIDX,  ZOUT3,PCVH, zout_dnuc)                !boundary layer top level, outputs, high vegetation
    
    CALL GSTATS(2501,1)

    ! updating wet and dry radii, aerosol density and water content
    DO IMODE=1,NMOD
      RW_MODE(IMODE)%d2(KIDIA:KFDIA,1:KLEV)=ZM6RP(KIDIA:KFDIA,1:KLEV,IMODE) ! m ( , KLEV, NMOD)
      DENS_MODE(IMODE)%d2(KIDIA:KFDIA,1:KLEV)=ZRHOP(KIDIA:KFDIA,1:KLEV,IMODE) ! kg/m3
    ENDDO
    DO IMODE=1,NMOD
      IF (sizeclass(IMODE)%lsoluble) THEN
        RWD_MODE(IMODE)%d2(KIDIA:KFDIA,1:KLEV) = ZM6DRY(KIDIA:KFDIA,1:KLEV,IMODE) ! m
        H2O_MODE(IMODE)%d2(KIDIA:KFDIA,1:KLEV) = ZWW(KIDIA:KFDIA,1:KLEV,IMODE) ! ?
      END IF
    ENDDO
    !-----------------------------------------------------------------

    CALL GSTATS(2502,0)

    ZTKEM1(KIDIA:KFDIA,1:KLEV) = 0._JPRB ! turbulent kinetic energy as zero for now as it is not used (YET!)

    ! LWP
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        IF ( ZAP(JL,JK) >=0.001_JPRB ) THEN
          ZTMPA = 1.0_JPRB/ZAP(JL,JK)
          LLIQCLD(JL,JK) = ( PLP(JL,JK)*ZTMPA  ) > ZEPSEC ! logical for liquid cloud
          LICECLD(JL,JK) = ( PIP(JL,JK)*ZTMPA  ) > ZEPSEC ! logical for ice cloud
          ZQLWP(JL,JK)   = MIN(MAX(0._JPRB, PLP(JL,JK)*ZTMPA), RCLDMAX)   ! lwc
          ZQIWP(JL,JK)   = MIN(MAX(0._JPRB, PIP(JL,JK)*ZTMPA), RCLDMAX)   ! iwc
          ZMSNOWACL_str(JL,JK) = PSNOWACL(JL,JK)*ZTMPA ! accretion rate of snow with cloud drop. in cloudy part for strat. wetdep
        ELSE
          LLIQCLD(JL,JK) = .FALSE.
          LICECLD(JL,JK) = .FALSE.
          ZQLWP(JL,JK)   = 0.0_JPRB
          ZQIWP(JL,JK)   = 0.0_JPRB
          ZMSNOWACL_str(JL,JK) = 0.0_JPRB
        END IF
      END DO
    END DO
    ZMSNOWACL_cov(KIDIA:KFDIA,1:KLEV) = 0.0_JPRB ! accretion rate of snow 0 for conv. wetdep

    !---find highest model level where there is cloud
    KTOP = KLEV !eehol: add as default (basically if no cloud so KTOP is assigned to lowest level)
    DO JK=1,KLEV
       IF ( ANY(ZAP(KIDIA:KFDIA,JK) >=0.001_JPRB) ) THEN ! if cloud fraction is larger than threshold
          KTOP = JK
          EXIT
       ENDIF
    END DO

    ! Default values for effective radii
    REFFL(KIDIA:KFDIA,1:KLEV,ZKROW) = ZDEF_RE_LIQ ! comes from liquid effective radius routine (PP_MIN_RE_UM)
    REFFI(KIDIA:KFDIA,1:KLEV,ZKROW) = ZDEF_RE_ICE ! comes from ice effective radius routine (ZDEFAULT_RE_UM)

    ! Cloud activation scheme
    CLDACT: IF ( NCLOUDACT == 1 ) THEN ! Morales and Nenes
       
       !IF ( LCONSIGW ) THEN !eehol: if using the constant sigma_w it is set to 0.8 otherwise use the TKE to calculate (NOT USED YET!)
       ZSIGMA_W(KIDIA:KFDIA,1:KLEV) = 0.8_JPRB
       !ELSE
       !   ZSIGMA_W(KIDIA:KFDIA,1:KLEV)= MAX(0.1_JPRB, (ZTUNPAR*((ZTKEM1(KIDIA:KFDIA,1:KLEV))**0.5_JPRB))) ! m/s
       !END IF

       ! Update mixing ratios with tendencies
       DO JT =1,NTRAC 
         ZXTP1(KIDIA:KFDIA,1:KLEV,JT) = ZXTM1(KIDIA:KFDIA,1:KLEV,JT) + ZXTTE(KIDIA:KFDIA,1:KLEV,JT) * TIME_STEP_LEN 
       END DO

       CALL AER_ACTIV(KIDIA,   KFDIA,  KTDIA,   KLON,    KLEV,   KSTGLO, &
                    &  PRS1,    PRSF1,    PTP,      ZQP,      ZQSAT,  &
                    &  PVERVEL, ZAP,     PLP,      PIP,              &
                    &  PLSM,    PGELAM,   PGEMU, & !PSLON,   PGEMU,  &
                    &  PGFL, YDMODEL, ZCDNCACT, ZICNC, REFFL(1:KLON,1:KLEV,ZKROW), REFFI(1:KLON,1:KLEV,ZKROW), &
                    &  ZSMAXMN, ZM6DRY, ZXTP1, KTRAC, ZSIGMA_W, ZFRACN, ZMIN_CDNC, ZDEF_CDNC, &
                    &  ZQLWP, LLIQCLD, LICECLD, ZDEF_RE_LIQ, ZDEF_RE_ICE)
       
       ! Store effective radii in PGFL
       PGFL(KIDIA:KFDIA,1:KLEV,YRE_LIQ%MP9_PH) = 1.0E-06_JPRB * REFFL(KIDIA:KFDIA,1:KLEV,ZKROW) ! convert um to meters
       IF (LAERICESED .OR. LAERICEAUTO) THEN
         PGFL(KIDIA:KFDIA,1:KLEV,YRE_ICE%MP9_PH) = 1.0E-06_JPRB * REFFI(KIDIA:KFDIA,1:KLEV,ZKROW) ! convert um to meters
       ENDIF
       
    ELSE IF (NCLOUDACT == 2) THEN ! AR&G scheme

       !IF ( LCONSIGW ) THEN !eehol: if using the constant sigma_w it is set to 0.8 otherwise use the TKE to calculate (NOT USED YET!!)
       ZTKEM1(KIDIA:KFDIA,1:KLEV) = ((1/ZTUNPAR)**2)*((0.8_JPRB)**2) !eehol: this is converted back to sigma_w in mo_activ.F90
       !ELSE
       !   ZTKEM1(KIDIA:KFDIA,1:KLEV) = ZTKEM1(KIDIA:KFDIA,1:KLEV)
       !END IF

       !---calculate updraft velocity
       ZWCAPE(KIDIA:KFDIA) = 0._JPRB !  CAPE as zero as it is not used
       
       CALL ACTIV_UPDRAFT(KFDIA, KLON, KLEV, ZKROW, & ! krow = 1
            ZTKEM1, ZWCAPE, PVERVEL, ZRHO, & ! turbulent kinetic energy, CAPE contr. to conv. vert. veloc. [m s-1], large scale vert. veloc.
            ZW, ZWPDF)
       
       DO JT = 1,NTRAC
          ZXTP1(KIDIA:KFDIA,1:KLEV,JT)  = ZXTM1(KIDIA:KFDIA,1:KLEV,JT) + ZXTTE(KIDIA:KFDIA,1:KLEV,JT) * TIME_STEP_LEN
       END DO
       
       IF (NCD_ACTIV == 2 .OR. NCCNDIAG > 0) THEN
          CALL HAM_ACTIV_KOEHLER_AB(KFDIA, KLON, KLEV, ZKROW, KTDIA, & ! krow=1 ktdia=1
               ZXTP1, PTP, ZA, ZB)
       END IF
       
       DO JCLASS = 1,NCLASS
          IF (SIZECLASS(JCLASS)%LSOLUBLE) THEN 
             ZRDRY(KIDIA:KFDIA,1:KLEV,JCLASS) = ZM6DRY(KIDIA:KFDIA,1:KLEV,JCLASS) !soluble modes rdry from rdry_m7
          ELSE
             ZRDRY(KIDIA:KFDIA,1:KLEV,JCLASS) = ZM6RP(KIDIA:KFDIA,1:KLEV,JCLASS)  !insoluble modes rdry from rwet_m7
          END IF
       END DO
       
       !calculate saturation water vapor pressure from saturation specific humidity
       DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
             ZESW(JL,JK)=(ZQSAT(JL,JK)*PRSF1(JL,JK))/(0.62198_JPRB)
          ENDDO
       ENDDO

       ZWPDF(KIDIA:KFDIA,1:KLEV,:) = MAX(ZWPDF(KIDIA:KFDIA,1:KLEV,:),1.0E-6_JPRB) !eehol: threshold vvel PDF to 1e-6 to get rid of 0 division
       
       CALL HAM_ACTIV_ABDULRAZZAK_GHAN( KFDIA, KLON, KLEV, ZKROW, KTDIA,  & ! in original 1 = ktdia... for diagnostics so krow=1 and ktdia=1
                                      & ZCDNCACT, ZESW, ZRHO,             & ! number of activated particles, saturation vapor pressure, air density
                                      & ZXTP1, PTP, PRSF1, ZQP,           & ! tracer mix rat, temperature, air pressure, spec. humid.
                                      & ZW, ZWPDF, ZA, ZB, ZRDRY,         & ! mean udr veloc, pdf of udr. veloc, Koehler A, Koehler B, dry radius
                                      & ZNACT, ZFRACN, ZSC, ZRC, ZSMAX)     ! num. act. part. per mode, frac ", crit. ssat., crit. radius, max ssat

       ! threshold fraction of activated particles to gridcells with only liquid clouds
       DO JCLASS = 1,NCLASS
         ZFRACN(KIDIA:KFDIA,1:KLEV,JCLASS) = MERGE(ZFRACN(KIDIA:KFDIA,1:KLEV,JCLASS),0._JPRB,LLIQCLD(KIDIA:KFDIA,1:KLEV))
       ENDDO

       ! threshold CDNC and ICNC to grid cells with only liquid or ice clouds
       ZCDNCACT(KIDIA:KFDIA,1:KLEV) = MERGE(ZCDNCACT(KIDIA:KFDIA,1:KLEV),1.0E6_JPRB*ZDEF_CDNC,LLIQCLD(KIDIA:KFDIA,1:KLEV)) !mask only values inside liq cloud
       ZICNC(KIDIA:KFDIA,1:KLEV) = MERGE(ZICNC(KIDIA:KFDIA,1:KLEV), RNICE, LICECLD(KIDIA:KFDIA,1:KLEV)) !mask only values inside ice cloud

    ELSE !Default values if no activation scheme is used
      
       ZFRACN(KIDIA:KFDIA,:,:) = 0._JPRB !init

       DO JL = KIDIA,KFDIA !eehol: add CDNC over land and over ocean as default values
          IF ( PLSM(JL) < 0.5_JPRB ) THEN !over ocean
             ZCDNCACT(JL,1:KLEV) = RCCNSEA*1.0E6_JPRB !from 1/cm3 to 1/m3
          ELSE !over land
             ZCDNCACT(JL,1:KLEV) = RCCNLND*1.0E6_JPRB !from 1/cm3 to 1/m3
          END IF
       END DO
       
       !calculate modewise fraction of activated particles
       !assume only KS, AS, CS modes to be activated
       DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
            ZNACT_TOT = ZCDNCACT(JL,JK) / ZRHO(JL,JK) !calculate total number of activated particles in #/kg
            IF ( ZXTM1(JL,JK,SIZECLASS(4)%IDT_NO) > 1.E-9_JPRB ) THEN
              ZFRAC_CS = MAX(ZNACT_TOT,0._JPRB) / ZXTM1(JL,JK,SIZECLASS(4)%IDT_NO)
            ELSE
              ZFRAC_CS = 0._JPRB
            ENDIF
            ZFRACN(JL,JK,4) = MAX(0._JPRB,MIN(ZFRAC_CS,1._JPRB)) !threshold between 0 and 1
            ZNACT_CS(JL,JK) = ZXTM1(JL,JK,SIZECLASS(4)%IDT_NO) * ZFRACN(JL,JK,4) !calculate activated number for CS mode

            IF ( ZXTM1(JL,JK,SIZECLASS(3)%IDT_NO) > 1.E-9_JPRB ) THEN
              ZFRAC_AS = MAX(ZNACT_TOT - ZNACT_CS(JL,JK),0._JPRB) / ZXTM1(JL,JK,SIZECLASS(3)%IDT_NO)
            ELSE
              ZFRAC_AS = 0._JPRB
            ENDIF
            ZFRACN(JL,JK,3) = MAX(0._JPRB,MIN(ZFRAC_AS, 1._JPRB)) !threshold between 0 and 1
            ZNACT_AS(JL,JK) = ZXTM1(JL,JK,SIZECLASS(3)%IDT_NO) * ZFRACN(JL,JK,3) !calculate activated number for AS mode

            IF ( ZXTM1(JL,JK,SIZECLASS(2)%IDT_NO) > 1.E-9_JPRB ) THEN
              ZFRAC_KS = MAX((ZNACT_TOT - ZNACT_CS(JL,JK) - ZNACT_AS(JL,JK)),0._JPRB) / ZXTM1(JL,JK,SIZECLASS(2)%IDT_NO)
            ELSE
              ZFRAC_KS = 0._JPRB
            ENDIF
            ZFRACN(JL,JK,2) = MAX(0._JPRB,MIN(ZFRAC_KS, 1._JPRB)) !threshold between 0 and 1
            ZNACT_KS(JL,JK) = ZXTM1(JL,JK,SIZECLASS(2)%IDT_NO) * ZFRACN(JL,JK,2) !calculate activated number for KS mode

            ZFRACN(JL,JK,2) = MERGE(ZFRACN(JL,JK,2),0._JPRB,LLIQCLD(JL,JK)) !fraction only where there is cloud
            ZFRACN(JL,JK,3) = MERGE(ZFRACN(JL,JK,3),0._JPRB,LLIQCLD(JL,JK)) !fraction only where there is cloud
            ZFRACN(JL,JK,4) = MERGE(ZFRACN(JL,JK,4),0._JPRB,LLIQCLD(JL,JK)) !fraction only where there is cloud
          END DO
       END DO

       ! threshold CDNC and ICNC to grid cells with only liquid or ice clouds
       ZCDNCACT(KIDIA:KFDIA,1:KLEV) = MERGE(ZCDNCACT(KIDIA:KFDIA,1:KLEV),1.0E6_JPRB*ZDEF_CDNC,LLIQCLD(KIDIA:KFDIA,1:KLEV)) !mask only values inside liq cloud
       ZICNC(KIDIA:KFDIA,1:KLEV) = MERGE(ZICNC(KIDIA:KFDIA,1:KLEV), RNICE, LICECLD(KIDIA:KFDIA,1:KLEV)) !mask only values inside ice cloud

    END IF CLDACT

    CALL GSTATS(2502,1)
    
    ! Store CDNC (number of activated particles) and ICNC as a number mixing ratio to tracer values, and in PGFL
    ZXTM1(KIDIA:KFDIA,1:KLEV,IDT_CDNC) = (MAX(ZCDNCACT(KIDIA:KFDIA,1:KLEV),((1.0E6_JPRB)*ZMIN_CDNC)))/ZRHO(KIDIA:KFDIA,1:KLEV) ! [#/kg] and threshold CDNC
    ZXTM1(KIDIA:KFDIA,1:KLEV,IDT_ICNC) = (1.0E6_JPRB)*ZICNC(KIDIA:KFDIA,1:KLEV)/ZRHO(KIDIA:KFDIA,1:KLEV) !ice crystal number conc = #/cm3 --> number mix rat [#/kg]

    PGFL(KIDIA:KFDIA,1:KLEV,YCDNC%MP9_PH) = MAX((1.0E-6_JPRB)*ZCDNCACT(KIDIA:KFDIA,1:KLEV),ZMIN_CDNC)  ! convert from #/m3 to #/cm3 and threshold minimum value to 1 cm-3
    PGFL(KIDIA:KFDIA,1:KLEV,YICNC%MP9_PH) = MAX( ZICNC(KIDIA:KFDIA,1:KLEV), RNICE) ! no conversion needed: already in #/cm3, just impose minimum value

    !eehol: update tendency of CDNC and ICNC (calculate only the newly formed droplets)
    ZXTTE_CLD(KIDIA:KFDIA,1:KLEV,IDT_CDNC) = (ZXTM1(KIDIA:KFDIA,1:KLEV,IDT_CDNC) - ZCDNC_temp(KIDIA:KFDIA,1:KLEV))/time_step_len
    ZXTTE_CLD(KIDIA:KFDIA,1:KLEV,IDT_ICNC) = (ZXTM1(KIDIA:KFDIA,1:KLEV,IDT_ICNC) - ZICNC_temp(KIDIA:KFDIA,1:KLEV))/time_step_len
    
    !-----------------------------------------------------------------
    !--> Calculation of liquid effective radii (Note: already done if NCLOUDACT=1)
    IF (NCLOUDACT == 2 .OR. NCLOUDACT == 0 ) THEN

      DO JK=1,KLEV
        DO JL=KIDIA,KFDIA
          ! effective radius (in um) calculated similarly as in radlswr.F90 
          ! 2.387e-10 is 3/(4*pi*rho_liq*10^6)  [10^6 for N in right units]
          ZRE_LIQ(JL,JK) = 1.E+06_JPRB*(2.387e-10_JPRB*ZRHO(JL,JK)*ZQLWP(JL,JK)/PGFL(JL,JK,YCDNC%MP9_PH))**0.333_JPRB
        END DO
      END DO

      ! Add liq. eff. rad. to HAM variables (only if there is liquid cloud else minimum value)
      REFFL(KIDIA:KFDIA,1:KLEV,ZKROW) = MERGE(ZRE_LIQ(KIDIA:KFDIA,1:KLEV), ZDEF_RE_LIQ, LLIQCLD(KIDIA:KFDIA,1:KLEV))

      ! add effective radii to PGFL fields
      PGFL(KIDIA:KFDIA,1:KLEV,YRE_LIQ%MP9_PH) = 1.0E-06_JPRB * REFFL(KIDIA:KFDIA,1:KLEV,ZKROW) ! convert um to meters and save to PGFL fields
    ENDIF

    !--> Calculation of ice effective radii (Note: already done if NCLOUDACT=1 *and* (LAERICESED .OR. LAERICEAUTO) )
    IF (.NOT. (NCLOUDACT==1 .AND. (LAERICESED .OR. LAERICEAUTO) )) THEN
      CALL ICE_EFFECTIVE_RADIUS(YRERAD, YDSPP_CONFIG, KIDIA, KFDIA, KLON, KLEV, &
           &  PRSF1, PTP, ZAP, PIP, PSP, PGEMU, & ! pressure, temp, cloud fr., IWC, SWC, sine of latitude
           &  reffi(1:KLON,1:KLEV,ZKROW)) ! ice effective radius (updated to mo_activ variable 'reffi' which used in mo_ham_wetdep)

      ! only if there is ice cloud else minimum value
      REFFI(KIDIA:KFDIA,1:KLEV,ZKROW) = MERGE(REFFI(KIDIA:KFDIA,1:KLEV,ZKROW), ZDEF_RE_ICE, LICECLD(KIDIA:KFDIA,1:KLEV))

      ! add effective radii to PGFL fields
      PGFL(KIDIA:KFDIA,1:KLEV,YRE_ICE%MP9_PH) = 1.0E-06_JPRB * REFFI(KIDIA:KFDIA,1:KLEV,ZKROW) ! convert um to meters and save to PGFL fields
    ENDIF

    !-----------------------------------------------------------------
    !--- Mass conserving correction of negative tracer values:
    CALL XT_BORROW(KFDIA, KLON,  KLEV, KLEV+1, NTRAC, &
         PRSF1, PRS1, &
         ZXTM1, ZXTTE)

    !-----------------------------------------------------------------
    !--> Wet deposition for HAM-M7
    CALL GSTATS(2503,0)

    IF ( LAERSCAV ) THEN

      !--> initialize mixing ratios for wet deposition
      !    Only ZXTP1, in case no tracers subject to wet dep, since it may be used in drydep!
      DO JT = 1,NTRAC
        ZXTP1(KIDIA:KFDIA,1:KLEV,JT)  = ZXTM1(KIDIA:KFDIA,1:KLEV,JT) + ZXTTE(KIDIA:KFDIA,1:KLEV,JT) * TIME_STEP_LEN
      END DO

      !<-- call wetdep interface for wet deposition
      !-- interface to wet deposition routine (also from cuflx_subm)
      IF ( ANY(trlist%ti(:)%nwetdep > 0) ) THEN

        ! for calculating the rain and snow evaporation/formation variables used in wet deposition
        ! convective wet deposition values
        DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
            ZFLXR=PFPLCL(JL,JK-1)
            ZFLXS=PFPLCN(JL,JK-1)
            ZFLXRB=PFPLCL(JL,JK)
            ZFLXSB=PFPLCN(JL,JK)
            ZMRATEPR_cov(JL,JK) =  MAX(ZFLXRB-ZFLXR,1.E-10_JPRB)/ZDPG(JL,JK) ! [kg/kg/s]
            ZMRATEPS_cov(JL,JK) =  MAX(ZFLXSB-ZFLXS,1.E-10_JPRB)/ZDPG(JL,JK) ! [kg/kg/s]
            ZMRATEPR_cov(JL,JK) = ZMRATEPR_cov(JL,JK) * TIME_STEP_LEN ! time integrated
            ZMRATEPS_cov(JL,JK) = ZMRATEPS_cov(JL,JK) * TIME_STEP_LEN ! time integrated
            !same formula negatives/positives for evap or formation
            ZFEVAPR_cov(JL,JK) =  -1._JPRB*MIN(ZFLXRB-ZFLXR,0._JPRB) ! [kg/m2.s]
            ZFSUBLS_cov(JL,JK) =  -1._JPRB*MIN(ZFLXSB-ZFLXS,0._JPRB) ! [kg/m2.s]
            ZFPLCL(JL,JK) = 0.5_JPRB*(ZFLXR+ZFLXRB)
            ZFPLCN(JL,JK) = 0.5_JPRB*(ZFLXS+ZFLXSB)
          END DO
        END DO
        ! stratiform wet deposition values
        DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
            ZFLXR=PFPLSL(JL,JK-1)
            ZFLXS=PFPLSN(JL,JK-1)
            ZFLXRB=PFPLSL(JL,JK)
            ZFLXSB=PFPLSN(JL,JK)
            IF ( ZAP(JL,JK) >=0.001_JPRB ) THEN
              ZTMPA = 1.0_JPRB/ZAP(JL,JK)
              ZMRATEPR_str(JL,JK) = MAX(ZFLXRB-ZFLXR,1.E-10_JPRB)/ZDPG(JL,JK) ! [kg/kg/s]
              ZMRATEPS_str(JL,JK) = MAX(ZFLXSB-ZFLXS,1.E-10_JPRB)/ZDPG(JL,JK) ! [kg/kg/s]              
              ZMRATEPR_str(JL,JK) = ZMRATEPR_str(JL,JK) * TIME_STEP_LEN * ZTMPA ! time integrated
              ZMRATEPS_str(JL,JK) = ZMRATEPS_str(JL,JK) * TIME_STEP_LEN * ZTMPA ! time integrated
            ELSE
              ZMRATEPR_str(JL,JK) = 0.0_JPRB
              ZMRATEPS_str(JL,JK) = 0.0_JPRB
            END IF
            !same formula negatives/positives for evap or formation
            IF ( PCOVPTOT(JL,JK) > ZEPSEC ) THEN
              ZTMPA = 1.0_JPRB/PCOVPTOT(JL,JK)
              ZFEVAPR_str(JL,JK) =  -1._JPRB*MIN(ZFLXRB-ZFLXR,0._JPRB) * ZTMPA ! [kg/m2.s]
              ZFSUBLS_str(JL,JK) =  -1._JPRB*MIN(ZFLXSB-ZFLXS,0._JPRB) * ZTMPA ! [kg/m2.s]
              ZFPLSL(JL,JK) = 0.5_JPRB*(ZFLXR+ZFLXRB) * ZTMPA
              ZFPLSN(JL,JK) = 0.5_JPRB*(ZFLXS+ZFLXSB) * ZTMPA
            ELSE
              ZFEVAPR_str(JL,JK) = 0.0_JPRB
              ZFSUBLS_str(JL,JK) = 0.0_JPRB
              ZFPLSL(JL,JK) = 0.0_JPRB
              ZFPLSN(JL,JK) = 0.0_JPRB
            END IF
          END DO
        END DO

        ZLFRAC_SO2(KIDIA:KFDIA,:) = 0._JPRB ! zlfrac_so2 only needed in gas scavenging and this is off for now (put this zero)

        ! In cloudsc.f90 PLU is defined as convective condensate (understood as total condensate, ie "include liquid and ice phases").
        !  So we use ZIPDUM to set the ice part to zero and use total condensate as liquid for the convective case of wetdep.
        !  In the future, it would be better to calculate the actual fraction based on 7.6. This may be important to describe the
        !  liquid and ice fractions correctly.
        ZLP(KIDIA:KFDIA,1:KLEV) = ZQLWP(KIDIA:KFDIA,1:KLEV)  ! temporary variable for cloud water content     (modified in wetdep - stratif case)
        ZIP(KIDIA:KFDIA,1:KLEV) = ZQIWP(KIDIA:KFDIA,1:KLEV)  ! temporary variable for cloud ice water content (modified in wetdep - stratif case)
        ZIPDUM(KIDIA:KFDIA,1:KLEV) = 0._JPRB               ! temporary variable for cloud ice water content (modified in wetdep - convec case)
        ZLPU(KIDIA:KFDIA,1:KLEV) = PLU(KIDIA:KFDIA,1:KLEV) ! temporary variable for cloud water content     (modified in wetdep - convec case)

        ! Double call to wet deposition. One for convective case and one for stratiform case:
        
        ! WETDEP CONVECTIVE CASE
        
        !-- initialise in-cloud and interstitial mixing ratios
        !   set both equal to tracer mixing ratio as starting point
        !   ham_wet_chemistry will re-compute these values if lham=true
        DO JT = 1,NTRAC
          ZXTP1(KIDIA:KFDIA,1:KLEV,JT)  = ZXTM1(KIDIA:KFDIA,1:KLEV,JT) + ZXTTE(KIDIA:KFDIA,1:KLEV,JT) * TIME_STEP_LEN
          ZXTP1C(KIDIA:KFDIA,1:KLEV,JT) = ZXTP1(KIDIA:KFDIA,1:KLEV,JT)
          ZXTP10(KIDIA:KFDIA,1:KLEV,JT) = ZXTP1(KIDIA:KFDIA,1:KLEV,JT)
        END DO

        ZDUMMY(KIDIA:KFDIA,:) = 0._JPRB         ! output: massfix boundary condition (updated in mo_hammoz_wetdep)
        ZWDEP_SCAV_IC(KIDIA:KFDIA,:) = 0._JPRB  ! output: diagnostic wdep in-cloud
        ZWDEP_SCAV_BC(KIDIA:KFDIA,:) = 0._JPRB  ! output: diagnostic wdep below cloud
        ZDUM2D(KIDIA:KFDIA,1:KLEV) = 0._JPRB    ! dummy "fraction of grid box covered by precip" for conv. case (zero at start as in ECHAM)
        ZDUM3D(KIDIA:KFDIA,1:KLEV,:) = 0._JPRB  ! dummy "previous tracer mixing ratio" for conv. case
        ZFUXT3D(KIDIA:KFDIA,1:KLEV,:) = 0._JPRB ! initialize updr mass flux for conv. case only (updated in wetdep so put 0 at start)

        LSTRAT = .FALSE. !False for convective case
        IF (.NOT. LSTRAT) THEN
          CALL XT_CONV_MASSFIX(KFDIA, KLON, KLEV, KLEV+1, NTRAC, ZKROW, PRSF1, PRS1, ZXTTE, .TRUE., ZDUMMY) ! call convective mass conserving (init zxtte_old)
                
          ZMFU(KIDIA:KFDIA,1:KLEV) = MIN(PMFU(KIDIA:KFDIA,1:KLEV),ZDPG(KIDIA:KFDIA,1:KLEV)/TIME_STEP_LEN)!constraing Conv. mass flux up

          CALL WETDEP_INTERFACE(KFDIA, KLON, KLEV, KTOP, ZKROW, LSTRAT, & ! ktop = cloud top level index, lstrat = FALSE for conv. case
                  ZDPG,  ZMRATEPR_COV, ZMRATEPS_COV, ZMSNOWACL_cov,     & ! dp/g, evap. of rain, subl. of snow, accr. rate of snow with cl. drop in-cl.
                  ZLPU,  ZIPDUM,                                        & ! cloud water content, cloud ice water content
                  ZM6RP,  ZM6DRY,                                    & ! m7 aerosol: to replace rwet_m7, dry radius for soluble modes [cm]
                  REFFI,  REFFL,                                     & ! effective radii
                  ZCDNCACT, ZFRACN,                                  & ! number/fraction of activated particles per mode
                  PTP, ZXTM1, ZLFRAC_SO2,                            & ! temperature, prev. mixing ratio, zlfrac_so2 only needed in gas scavenging (0 for now)
                  ZXTTE, ZXTP10, ZXTP1C,                             & ! tendencies/mixing ratios (in/out)
                  ZFPLCL, ZFPLCN, ZFEVAPR_cov, ZFSUBLS_cov,          & ! rain flux, snow flux, 
                  ZMFU, ZFUXT3D,                                     & ! conv flux, updraft mass flux (updated in wetdep)
                  ZAP_COV,  ZDUM2D, ZRHO, ZDUMMY, ZWDEP_SCAV_IC, ZWDEP_SCAV_BC)  ! cloud frac., precip. frac., air dens., in/output*3

          CALL XT_CONV_MASSFIX(KFDIA, KLON, KLEV, KLEV+1, NTRAC, ZKROW, PRSF1, PRS1, ZXTTE, .FALSE., ZDUMMY) ! call convective mass conserving
        END IF

        !Add convective case wet removal fluxes to diagnostics
        DO JMASS=1,NAEROCOMP
          JY=KAERO(ind_oifs_ham%ind_mass_OIFS(JMASS))
          WDEPOUT_2D   (KIDIA:KFDIA,JY) = WDEPOUT_2D   (KIDIA:KFDIA,JY) + ZDUMMY       (KIDIA:KFDIA, ind_oifs_ham%ind_mass_ham(JMASS))
          WDEPOUT_IC_2D(KIDIA:KFDIA,JY) = WDEPOUT_IC_2D(KIDIA:KFDIA,JY) + ZWDEP_SCAV_IC(KIDIA:KFDIA, ind_oifs_ham%ind_mass_ham(JMASS))
          WDEPOUT_BC_2D(KIDIA:KFDIA,JY) = WDEPOUT_BC_2D(KIDIA:KFDIA,JY) + ZWDEP_SCAV_BC(KIDIA:KFDIA, ind_oifs_ham%ind_mass_ham(JMASS))
        END DO

        DO JCLASS=1,NCLASS
          JY=KAERO(ind_oifs_ham%ind_class_OIFS(JCLASS))
          WDEPOUT_2D   (KIDIA:KFDIA,JY) = WDEPOUT_2D   (KIDIA:KFDIA,JY) + ZDUMMY       (KIDIA:KFDIA, ind_oifs_ham%ind_class_ham(JCLASS))
          WDEPOUT_IC_2D(KIDIA:KFDIA,JY) = WDEPOUT_IC_2D(KIDIA:KFDIA,JY) + ZWDEP_SCAV_IC(KIDIA:KFDIA, ind_oifs_ham%ind_class_ham(JCLASS))
          WDEPOUT_BC_2D(KIDIA:KFDIA,JY) = WDEPOUT_BC_2D(KIDIA:KFDIA,JY) + ZWDEP_SCAV_BC(KIDIA:KFDIA, ind_oifs_ham%ind_class_ham(JCLASS))
        END DO

        ! WETDEP STRATIFORM CASE

        !-- initialise in-cloud and interstitial mixing ratios
        !   set both equal to tracer mixing ratio as starting point
        !   ham_wet_chemistry will re-compute these values if lham=true
        DO JT = 1,NTRAC
          ZXTP1(KIDIA:KFDIA,1:KLEV,JT)  = ZXTM1(KIDIA:KFDIA,1:KLEV,JT) + ZXTTE(KIDIA:KFDIA,1:KLEV,JT) * TIME_STEP_LEN
          ZXTP1C(KIDIA:KFDIA,1:KLEV,JT) = ZXTP1(KIDIA:KFDIA,1:KLEV,JT)
          ZXTP10(KIDIA:KFDIA,1:KLEV,JT) = ZXTP1(KIDIA:KFDIA,1:KLEV,JT)
        END DO

        ZDUMMY(KIDIA:KFDIA,:) = 0._JPRB        ! output: massfix boundary condition (updated in mo_hammoz_wetdep)
        ZWDEP_SCAV_IC(KIDIA:KFDIA,:) = 0._JPRB ! output: diagnostic wdep in-cloud
        ZWDEP_SCAV_BC(KIDIA:KFDIA,:) = 0._JPRB ! output: diagnostic wdep below cloud
        ZDUM2D(KIDIA:KFDIA,1:KLEV) = 0._JPRB   ! dummy conv flux for strat. case
        ZDUM3D(KIDIA:KFDIA,1:KLEV,:) = 0._JPRB ! dummy updraft mass flux for strat. case

        LSTRAT = .TRUE. !True for strat case, large scale
        CALL WETDEP_INTERFACE(KFDIA, KLON, KLEV, KTOP, ZKROW, LSTRAT, & ! ktop = cloud top level index, lstrat = TRUE for strat. case
                ZDPG,  ZMRATEPR_STR, ZMRATEPS_STR, ZMSNOWACL_str,  & ! dp/g, evap. of rain, subl. of snow, accr. rate of snow with cl. drop in-cl.
                ZLP,  ZIP,                                         & ! cloud water content, cloud ice water content
                ZM6RP,  ZM6DRY,                                    & ! m7 aerosol: to replace rwet_m7, dry radius for soluble modes [cm]
                REFFI,  REFFL,                                     & ! effective radii
                ZCDNCACT, ZFRACN,                                  & ! number/fraction of activated particles per mode
                PTP, ZXTM1, ZLFRAC_SO2,                            & ! temperature, prev. mixing ratio, zlfrac_so2 only needed in gas scavenging (0 for now)
                ZXTTE, ZXTP10, ZXTP1C,                             & ! tendencies/mixing ratios (in/out)
                ZFPLSL, ZFPLSN, ZFEVAPR_str, ZFSUBLS_str,          & ! rain flux, snow flux, 
                ZDUM2D, ZDUM3D,                                    & ! zeroes as these are not needed in strat. case
                ZAP,  PCOVPTOT, ZRHO, ZDUMMY, ZWDEP_SCAV_IC, ZWDEP_SCAV_BC)  ! cloud frac., precip. frac., air dens., in/output*3

        ! Add stratiform case wet removal fluxes to diagnostics
        DO JMASS=1,NAEROCOMP
          JY=KAERO(ind_oifs_ham%ind_mass_OIFS(JMASS))
          WDEPOUT_2D   (KIDIA:KFDIA,JY) = WDEPOUT_2D   (KIDIA:KFDIA,JY) + ZDUMMY       (KIDIA:KFDIA, ind_oifs_ham%ind_mass_ham(JMASS))
          WDEPOUT_IC_2D(KIDIA:KFDIA,JY) = WDEPOUT_IC_2D(KIDIA:KFDIA,JY) + ZWDEP_SCAV_IC(KIDIA:KFDIA, ind_oifs_ham%ind_mass_ham(JMASS))
          WDEPOUT_BC_2D(KIDIA:KFDIA,JY) = WDEPOUT_BC_2D(KIDIA:KFDIA,JY) + ZWDEP_SCAV_BC(KIDIA:KFDIA, ind_oifs_ham%ind_mass_ham(JMASS))
        END DO

        DO JCLASS=1,NCLASS
          JY=KAERO(ind_oifs_ham%ind_class_OIFS(JCLASS))
          WDEPOUT_2D   (KIDIA:KFDIA,JY) = WDEPOUT_2D   (KIDIA:KFDIA,JY) + ZDUMMY       (KIDIA:KFDIA, ind_oifs_ham%ind_class_ham(JCLASS))
          WDEPOUT_IC_2D(KIDIA:KFDIA,JY) = WDEPOUT_IC_2D(KIDIA:KFDIA,JY) + ZWDEP_SCAV_IC(KIDIA:KFDIA, ind_oifs_ham%ind_class_ham(JCLASS))
          WDEPOUT_BC_2D(KIDIA:KFDIA,JY) = WDEPOUT_BC_2D(KIDIA:KFDIA,JY) + ZWDEP_SCAV_BC(KIDIA:KFDIA, ind_oifs_ham%ind_class_ham(JCLASS))
        END DO

      END IF
    END IF
    CALL GSTATS(2503,1)

    !<-- End wet deposition for HAM-M7
    !-----------------------------------------------------------------

    !--- Mass conserving correction of negative tracer values:
    CALL XT_BORROW(KFDIA, KLON,  KLEV, KLEV+1, ntrac, &
         PRSF1, PRS1, &
         ZXTM1, ZXTTE)

    !-----------------------------------------------------------------
    !--> Sedimentation for HAM-M7
    CALL GSTATS(2504,0)

    IF (LAERSEDIM) THEN
      IF ( ANY(trlist%ti(:)%nsedi > 0) ) THEN

        ZTENCIH(KIDIA:KFDIA,1:KLEV,1:NTRAC)=ZXTTE(KIDIA:KFDIA,1:KLEV,1:NTRAC)

        CALL SEDI_INTERFACE(KLON, KFDIA, KLEV, ZKROW,   &
             PTP, ZQP, PRSF1, PRS1, & ! temperature, specific humidity, pressure at full level, pressure at half level
             ZM6RP, ZRHOP, & ! mean mode actual radius [m], mean mode particle density [kg m-3]
             ZXTM1, ZXTTE, ZSEDIFLUX, ZSEDIFLUXSURF) ! tracer mixing ratios and tendency (sediflux for diagnostics)

        SEDOUT(KIDIA:KFDIA, 1:KLEV,1:NTRAC) = ZTENCIH(KIDIA:KFDIA, 1:KLEV,1:NTRAC) - ZXTTE(KIDIA:KFDIA, 1:KLEV,1:NTRAC)
        DO JK=1,KLEV
          DO JCLASS=1,NCLASS
            SEDOUT_2D(KIDIA:KFDIA,KAERO(ind_oifs_ham%ind_class_OIFS(JCLASS)))=SEDOUT_2D(KIDIA:KFDIA,KAERO(ind_oifs_ham%ind_class_OIFS(JCLASS))) + ZSEDIFLUX(KIDIA:KFDIA, JK,ind_oifs_ham%ind_class_HAM(JCLASS))
          END DO
          DO JMASS=1,NAEROCOMP
            SEDOUT_2D(KIDIA:KFDIA,KAERO(ind_oifs_ham%ind_mass_OIFS(JMASS)))=SEDOUT_2D(KIDIA:KFDIA,KAERO(ind_oifs_ham%ind_mass_OIFS(JMASS))) + ZSEDIFLUX(KIDIA:KFDIA,JK,ind_oifs_ham%ind_mass_HAM(JMASS))
          END DO
        END DO
      END IF
    ENDIF
    CALL GSTATS(2504,1)

    !<-- End sedimentation for HAM-M7
    !-----------------------------------------------------------------
    
    !--- Mass conserving correction of negative tracer values:
    CALL xt_borrow(KFDIA, KLON,  KLEV, KLEV+1, ntrac, &
         PRSF1, PRS1, &
         ZXTM1, ZXTTE)

    !-----------------------------------------------------------------
    !--> Dry deposition for HAM-M7
    IF (.NOT.LAERSURF) THEN
      PAERDDP(:,:)  =0._JPRB
    ELSEIF (LAERSURF) THEN

      !*             DRY DEPOSITION INCLUDED AS MODIFICATION TO SURFACE FLUXES
      !              ---------------------------------------------------------
      CALL GSTATS(2505,0)

      IF (LAERDRYDP) THEN

        !--> variables not needed for aerosol dry deposition
        ZCFML(:) = 0._JPRB
        ZCFMW(:) = 0._JPRB
        ZCFMI(:) = 0._JPRB
        ZCFNCL(:) = 0._JPRB
        ZCFNCW(:) = 0._JPRB
        ZCFNCI(:) = 0._JPRB
        ZEPDU2 = 0._JPRB
        ZKAP = 0._JPRB
        ZGEOM1(:,:) = 0._JPRB
        ZRIL(:) = 0._JPRB
        ZRIW(:) = 0._JPRB
        ZRII(:) = 0._JPRB
        ZTVIR1(:,:) = 0._JPRB
        ZTVL(:) = 0._JPRB
        ZTVW(:) = 0._JPRB
        ZTVI(:) = 0._JPRB
        ZAZ0(:) = 0._JPRB
        ZFRL(:) = 0._JPRB
        ZSRFL(:) = 0._JPRB
        ZFOREST(:) = 0._JPRB
        ZTSI(:) = 0._JPRB
        ZAZ0L(:) = 0._JPRB
        ZAZ0I(:) = 0._JPRB
        ZCDNI(:) = 0._JPRB

        !--> variables calculated for dry deposition
        DO JL = KIDIA,KFDIA
          IF ( PLSM(JL) < 0.99_JPRB ) THEN
            ZLOLAND(JL) = .FALSE.
            ZAZ0W(JL) = PZ0M(JL)
          ELSE
            ZLOLAND(JL) = .TRUE.
            ZAZ0W(JL) = 0._JPRB
          END IF
          ZAZ0W(JL)  = MAX(1.0E-5_JPRB,ZAZ0W(JL))  ! threshold roughness length to min value
          ZFRW(JL)   = MAX(0.,1.-PLSM(JL)-PCI(JL)) ! water fraction = 1 - land mask - sea ice fraction
          ZCVS(JL)   = PFRTI(JL,5)+PFRTI(JL,7)     ! snow cover fraction = Snow on low-veg + snow on bare-soil + snow under high-veg
          ZCVW(JL)   = PFRTI(JL,3)                 ! wet skin fraction
          ZVGRAT(JL) = PCVL(JL)+PCVH(JL)           ! vegetation ratio = low veg. cover + high veg. cover
          ZCDNL(JL)  = PAERUST(JL)                 ! adding ustar to not used variable
          ZCDNW(JL)  = LOG(ZDZ(JL,KLEV)/PZ0M(JL))/(VKARMAN*PAERUST(JL)) ! calculate aerodyn. resistance on surface to not used variable
        END DO
        
        !--> init values
        ZTENCIH(KIDIA:KFDIA,1:KLEV,:) = ZXTTE(KIDIA:KFDIA,1:KLEV,:) ! init tendency before drydep
        ZXTEMS(KIDIA:KFDIA,:)         = 0._JPRB                     ! surface emissions as zero for input
        ZXTMD1(KIDIA:KFDIA,1:KLEV,:)  = 0._JPRB
        ZXTMD1(KIDIA:KFDIA,1:KLEV,:)  = ZXTM1(KIDIA:KFDIA,1:KLEV,:) + (ZTENCIH(KIDIA:KFDIA,1:KLEV,:) * TIME_STEP_LEN) ! update mixrat with tendency
        ZVDEP(KIDIA:KFDIA,:)          = 0._JPRB                     ! ddep velocity as zero

        ! RCHG: Recommendation, those subroutines specific of m7 should have 
        !       m7 in the name not sure if this is specific or general/common 
        !       but adapted to m7.
        CALL DRYDEP_INTERFACE(KLON, KFDIA,  KLEV, ZKROW,                    &
            & ZQP(:,KLEV), ZQSAT(:,KLEV), PTP(:,KLEV), ZCFML, ZCFMW, ZCFMI, &
            & ZCFNCL, ZCFNCW, ZCFNCI,                                       &
            & ZEPDU2, ZKAP, PUP, PVP, ZGEOM1, ZRIL, ZRIW,                   &
            & ZRII,                                                         &
            & ZTVIR1, ZTVL, ZTVW, ZTVI, ZAZ0,                               &
            & PTP(:,KLEV), ZLOLAND,                                         &
            & ZM6RP, ZRHOP,                                                 & ! M7
            & ZFRL,   ZFRW,  PCI,     ZCVS,   ZCVW,     ZVGRAT,             &
            & ZSRFL,  PUP(:,KLEV),     PVP(:,KLEV),                         & !eehol: FIXME 10m u and v wind from lowest level.. needs to be revised in future!!
            & ZXTEMS, ZXTMD1, ZRHO(:,KLEV), PRS1, ZFOREST, ZTSI,            & !air dens lowest, air press at int.
            & ZAZ0L, ZAZ0W, ZAZ0I, ZCDNL, ZCDNW, ZCDNI, ZDDEPFLUX, ZVDEP)     !ZCDNL and ZCDNW used for ustar and aerodyn. resist.
           
        
        !--> modify tendency at surface according to changes in surface emissions
        DO JT = 1,NTRAC
          DO JL = KIDIA,KFDIA
            ZXTTE(JL,KLEV,JT) = ZTENCIH(JL,KLEV,JT) + ((ZXTEMS(JL,JT)*RG)/(ZDP(JL,KLEV)))
          END DO
        END DO

      ENDIF ! LAERDRYDP
    END IF
    CALL GSTATS(2505,1)

    !<-- End dry deposition for HAM-M7
    !-----------------------------------------------------------------

    !--- Mass conserving correction of negative tracer values:
    CALL xt_borrow(KFDIA, KLON,  KLEV, KLEV+1, ntrac, &
         PRSF1, PRS1, &
         ZXTM1, ZXTTE)

    !-----------------------------------------------------------------
    !--> Add HAM modified tendency back to PTENC (OIFS values)

    !number
    DO JCLASS=1,NCLASS
      PTENC(KIDIA:KFDIA,1:KLEV,KAERO(ind_oifs_ham%ind_class_OIFS(JCLASS))) = ZXTTE(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_class_HAM(JCLASS))
    END DO
    !mass
    DO JMASS=1,NAEROCOMP
      PTENC(KIDIA:KFDIA,1:KLEV,KAERO(ind_oifs_ham%ind_mass_OIFS(JMASS))) = ZXTTE(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_mass_HAM(JMASS))
    END DO
    !gas
    IF(LAERCHEM ) THEN
      DO JGAS=1,SUBM_NGASSPEC
        PTENC(KIDIA:KFDIA,1:KLEV,KCHEM(ind_oifs_ham%ind_gas_OIFS(JGAS))) = ZXTTE(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_gas_HAM(JGAS))
      END DO    
    END IF
    !!ELSE
    !!  DO JGAS=1,SUBM_NGASSPEC
    !!    PTENC(KIDIA:KFDIA,1:KLEV,KAERO(ind_oifs_ham%ind_gas_OIFS(JGAS))) = ZXTTE(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_gas_HAM(JGAS))
    !!  END DO
    !!END IF

    ! RCHG -> not sure best way to solve here. I commented to avoid segmentation fault 
    !         but it may be avoided with other more specific flag. Anyway something was 
    !         needed to avoid core-dump. 
    !         The problem with these arrays in loop below is that the indices:
    !          ind_oifs_ham%ind_cloud_HAM(JCLOUD)) 
    !          KAERO(ind_oifs_ham%ind_cloud_OIFS(JCLOUD)))
    !         are not propoperly set up by hamm7_init.F90
    !         -- this need to be solved probably in hamm7_init.F90 which detect the 
    !         problems with these tracers about CCN. 
    !cloud variables
    DO JCLOUD=1,2 !CDNC and ICNC
       !PTENC(KIDIA:KFDIA,1:KLEV,KAERO(ind_oifs_ham%ind_cloud_OIFS(JCLOUD))) = ZXTTE(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_cloud_HAM(JCLOUD))
       PTENC(KIDIA:KFDIA,1:KLEV,KAERO(ind_oifs_ham%ind_cloud_OIFS(JCLOUD))) = (1.0E-6_JPRB) * ZRHO(KIDIA:KFDIA,1:KLEV) * ZXTTE_CLD(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_cloud_HAM(JCLOUD))
    END DO
    !<-- End adding HAM modified tendency back to PTENC
    !-----------------------------------------------------------------
    

! ZTSO4 is filled only if LAERCHEM=F, which is not possible with M7 - Commented out
!THIS-IS-NEVER-USED   
!THIS-IS-NEVER-USED   ! write flux to extra fields for diagnostic of aerosol 'chemical' conversion  
!THIS-IS-NEVER-USED   IF (.NOT. LAERCHEM .AND. LCHEM_DIA) THEN
!THIS-IS-NEVER-USED     CALL CHEM_INEXT( KIDIA , KFDIA  , KLON , KLEV , 1, 1,  &
!THIS-IS-NEVER-USED        &    ZDP, PTSPHY, ZTSO4, ZTENC0,PEXTRA(:,ISO4_C,IEXTR_CH))
!THIS-IS-NEVER-USED   END IF


!*         4.    ELIMINATION OF NEGATIVE PROGNOSTIC AEROSOL CONCENTRATIONS
!                 ---------------------------------------------------------

IF (LAERNGAT) THEN

  IF (LCHEM_DIA) THEN
    ZTAERO0(KIDIA:KFDIA,1:KLEV,1:NACTAERO) =  ZTAEROK(KIDIA:KFDIA,1:KLEV,1:NACTAERO)
  ENDIF

  DO JAER=1,NACTAERO

    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZAER(JL,JK) = ZCEN(JL,JK,KAERO(JAER))
        ZTAER(JL,JK)= PTENC(JL,JK,KAERO(JAER))
      ENDDO
    ENDDO

    CALL AER_NEGAT &
         & ( YREAERATM, KIDIA   , KFDIA, KLON , KLEV, &
         &   PTSPHY  , &
         &   ZAER    , ZTAER, PRS1, &
         &   ZAERNEG )  

    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZTAERO(JL,JK,JAER)=ZTAER(JL,JK)
      ENDDO
    ENDDO
    DO JL=KIDIA,KFDIA
      ZAERNGT(JL,JAER)=ZAERNEG(JL,KLEV)
    ENDDO

  ENDDO

  ! collect neg fix tendencies
  IF (LCHEM_DIA) THEN
    CALL CHEM_INEXT( KIDIA , KFDIA  , KLON , KLEV , NACTAERO, NACTAERO ,  &
         &    ZDP, PTSPHY, ZTAERO ,ZTAERO0, PEXTRA(:,NCHEM+1:NCHEM+NACTAERO,IEXTR_NG))
  ENDIF

  ! do not fix the tendencies for now, number concentration fixes will break the
  ! correlation between mass and number
  PTENC(KIDIA:KFDIA,1:KLEV,KAERO(1):KAERO(NACTAERO)) = ZTAERO(KIDIA:KFDIA,1:KLEV,1:NACTAERO)

ENDIF

!------------------------------------------------------------------------------
!*         5.       STORE ALL AEROSOL VERTICALLY INTEGRATED FLUXES
!                   ----------------------------------------------

DO JAER=1,NACTAERO
  DO JL=KIDIA,KFDIA
    PAERODDF(JL,JAER,1)=PAERSRC(JL,JAER) ! aerosol so4 source term
    PAERODDF(JL,JAER,2)=PAERDDP(JL,JAER) ! aerosol dry deposition
    PAERODDF(JL,JAER,3)=PAERSDM(JL,JAER) ! aerosol sedimentation 
    PAERODDF(JL,JAER,4)=0.0              ! (todo) so2 sink added to scavenging
    PAERODDF(JL,JAER,5)=0.0              ! (todo) scavenging (in-cloud & below cloud) so wet deposition
    PAERODDF(JL,JAER,6)=ZAERNGT(JL,JAER)
    PAERODDF(JL,JAER,7)=0.0              ! (todo) total AOD?
  ENDDO
ENDDO

!-----------------------------------------------------------------------
!*         6.      OPTICAL PROPERTIES
!                  -------------------------------------------------- 
IBLK=(KSTGLO-1)/KLON + 1

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    DO JAER=1,14
      PAER_TAU(JL,JK,JAER)   = YDAERM7%M7AOD(JL,JK,JAER,IBLK)
      PAER_SSA(JL,JK,JAER)   = YDAERM7%M7SSA(JL,JK,JAER,IBLK)
      PAER_ASYM(JL,JK,JAER ) = YDAERM7%M7ASYM(JL,JK,JAER,IBLK)
    ENDDO
    DO JAER=1,16    
      PAER_TAU_LW(JL,JK,JAER)= YDAERM7%M7AODLW(JL,JK,JAER,IBLK)
    ENDDO
  ENDDO
ENDDO


!*         6.1      Calculate optical properties only when radiation is called
!                   ----------------------------------------------------------
IF(MOD(NSTEP,NRADFR) == 0) THEN
  
  CALL GSTATS(2506,0)
  ZAER_TAU(KIDIA:KFDIA,:,:,:)  = 0.0_JPRB
  ZAER_SSA(KIDIA:KFDIA,:,:)    = 0.0_JPRB
  ZAER_ASYM(KIDIA:KFDIA,:,:)   = 0.0_JPRB
  ZAER_TAU_LW(KIDIA:KFDIA,:,:) = 0.0_JPRB

  ZAOD_DIAG(KIDIA:KFDIA,:)= 0._JPRB
  ZSSA_DIAG(KIDIA:KFDIA,:)= 0._JPRB
  ZASY_DIAG(KIDIA:KFDIA,:)= 0._JPRB
  ZAER_TAU_DIAG(KIDIA:KFDIA,:,:)  = 0.0_JPRB
  ZAER_SSA_DIAG(KIDIA:KFDIA,:,:)  = 0.0_JPRB
  ZAER_ASYM_DIAG(KIDIA:KFDIA,:,:) = 0.0_JPRB
  ZAER_LidarR_DIAG(KIDIA:KFDIA,:,:) = 0.0_JPRB

SELECT CASE (NAEROOPT) 

CASE (0)

   ! Optical properties is not calculated

CASE (1)
   ! Use TM5 codes to calculate optical properties (optical properties for LW = 0)
      !--> Add HAM updated tendency to ZTAERO and use that in optics
   ZTAERO(KIDIA:KFDIA,:,:) = 0._JPRB
   DO JAER=1,NACTAERO
      DO JK=1,KLEV
         DO JL=KIDIA,KFDIA
            ZTAERO(JL,JK,JAER)= PTENC(JL,JK,KAERO(JAER))
         ENDDO
      ENDDO
   ENDDO

   ZAEROK(KIDIA:KFDIA,:,:) = ZAEROK(KIDIA:KFDIA,:,:) + ZTAERO(KIDIA:KFDIA,:,:)*TIME_STEP_LEN

   CALL TM5M7_OPTICS_AOP_GET( YGFL, YREAERSRC, KIDIA,KFDIA, KLON, KLEV,NACTAERO, &
   &                          NASWBAND, ASWBAND, 1, .false., &
   &                          ZRHO, ZAEROK,RW_MODE,RWD_MODE,H2O_MODE,&
   &                          ZAER_TAU, ZAER_SSA, ZAER_ASYM)

   PAOD(KIDIA:KFDIA,:)=0._JPRB
   DO JK = 1, KLEV
     DO JL = KIDIA,KFDIA
       DO IW=1,NASWBAND
         PAER_TAU(JL,JK,IW) = ZAER_TAU(JL,JK,IW,1)*(PGEOH(JL,JK-1) - PGEOH(JL,JK))/RG
         PAER_SSA(JL,JK,IW) = ZAER_SSA(JL,JK,IW)
         PAER_ASYM(JL,JK,IW)= ZAER_ASYM(JL,JK,IW)
         !PAOD(JL,IW)=ZAER_TAU(JL,JK,IW,1)*(PGEOH(JL,JK-1) - PGEOH(JL,JK))+PAOD(JL,IW)
       ENDDO
       DO IW=1,16
         PAER_TAU_LW(JL,JK,IW)=0.0_JPRB
       END DO
     ENDDO
   ENDDO

 CASE (2)
   ! Use HAM codes to calculate optical properties
   LWBANDS=16
   PRS1D(KIDIA:KFDIA,:) = PRS1(KIDIA:KFDIA,1:KLEV)-PRS1(KIDIA:KFDIA,0:KLEV-1)

   LDIAG_AEROPT = (NAERO_WVL_DIAG >0) .AND. (YGFL%NAERO_WVL_DIAG_TYPES>0)

   DO IW = 1, NAERO_WVL_DIAG
      LAMBDA_DIAG(IW) = YGFL%YAERO_WVL_DIAG_NL(IW)%IWVL *1.0E-9_JPRB ! nm to m
   ENDDO
   !CALL ham_rad_cache(KLON,KLEV)
   ZXTM0(KIDIA:KFDIA,1:KLEV,:) = ZXTM1(KIDIA:KFDIA,1:KLEV,:) + ZXTTE(KIDIA:KFDIA,1:KLEV,:)*time_step_len

   CALL HAM_RAD(KFDIA, KLON, KLEV, ZKROW, LWBANDS, NASWBAND, ZXTM0, PRS1D, &
        & ZAER_TAU(:,:,:,1), ZAER_SSA, ZAER_ASYM, ZAER_TAU_LW, ZM6RP, &
        & LDIAG_AEROPT,NAERO_WVL_DIAG,YGFL%NAERO_WVL_DIAG_TYPES, &
        & LAMBDA_DIAG, ZAER_TAU_DIAG, ZAER_SSA_DIAG, ZAER_ASYM_DIAG, ZAER_LidarR_DIAG)
   !CALL ham_rad_cache_cleanup

   DO JK = 1, KLEV
     DO JL = KIDIA,KFDIA
       DO IW=1,NASWBAND
         PAER_TAU(JL,JK,IW) = ZAER_TAU(JL,JK,IW,1)!*(PGEOH(JL,JK-1) - PGEOH(JL,JK))
         PAER_SSA(JL,JK,IW) = ZAER_SSA(JL,JK,IW)
         PAER_ASYM(JL,JK,IW)= ZAER_ASYM(JL,JK,IW)
       ENDDO
       DO IW=1,16
         PAER_TAU_LW(JL,JK,IW)=ZAER_TAU_LW(JL,JK,IW)
       END DO

     ENDDO
   ENDDO

   ! Vertically integrated optical properties at diagnostic wavelenghts
   DO JK = 1, KLEV
     DO JL = KIDIA,KFDIA
       DO IW=1,NAERO_WVL_DIAG
         ZAOD_DIAG(JL,IW)= ZAER_TAU_DIAG(JL,JK,IW)+ZAOD_DIAG(JL,IW)
         ZSSA_DIAG(JL,IW)= ZAER_SSA_DIAG(JL,JK,IW)*ZAER_TAU_DIAG(JL,JK,IW)+ZSSA_DIAG(JL,IW)
         ZASY_DIAG(JL,IW)= ZAER_ASYM_DIAG(JL,JK,IW)*ZAER_SSA_DIAG(JL,JK,IW)*ZAER_TAU_DIAG(JL,JK,IW)+ZASY_DIAG(JL,IW)
       END DO
     END DO
   END DO

   DO JL = KIDIA,KFDIA
     DO IW=1,NAERO_WVL_DIAG
       IF(ZSSA_DIAG(JL,IW)>0._JPRB) THEN
         ZASY_DIAG(JL,IW) = ZASY_DIAG(JL,IW)/ZSSA_DIAG(JL,IW)! AVERAGE
       ENDIF
       IF(ZAOD_DIAG(JL,IW)>0._JPRB) THEN
         ZSSA_DIAG(JL,IW) = ZSSA_DIAG(JL,IW)/ZAOD_DIAG(JL,IW)! AOD AVERAGE
       ENDIF
     END DO
   END DO

 END SELECT

 ! Impose safety limits
 IF ((NAEROOPT .NE. 0).AND. LSAFETY) THEN
   ! REPSCAER was needed to avoid corruption in radiation_scheme at
   ! some point. Such limit could be removed once we are satisfied
   ! with computed optical depths
   DO JK = 1, KLEV
     DO JL = KIDIA,KFDIA
       DO IW=1,NASWBAND
         PAER_TAU(JL,JK,IW) = MAX(REPSCAER,      PAER_TAU(JL,JK,IW) )
         PAER_SSA(JL,JK,IW) = MIN(MAX(0._JPRB,   PAER_SSA(JL,JK,IW) ), 1._JPRB)
         PAER_ASYM(JL,JK,IW)= MIN( MAX(-1._JPRB, PAER_ASYM(JL,JK,IW)), 1._JPRB)
       ENDDO
       DO IW=1,16
         PAER_TAU_LW(JL,JK,IW) = MAX(REPSCAER, PAER_TAU_LW(JL,JK,IW))
       END DO
     ENDDO
   ENDDO
 ENDIF
 
 CALL GSTATS(2506,1)
ENDIF  ! It's a time step when radiation is called

!*         6.2      Vertically integrated optical properties
!                   ----------------------------------------
PAOD (KIDIA:KFDIA,:)=0._JPRB
PABS (KIDIA:KFDIA,:)=0._JPRB
PFAOD(KIDIA:KFDIA,:)=0._JPRB
PSSA (KIDIA:KFDIA,:)=0._JPRB
PASY (KIDIA:KFDIA,:)=0._JPRB
PAOD_LW(KIDIA:KFDIA,:)=0._JPRB

DO JK = 1, KLEV
  DO JL = KIDIA,KFDIA
    DO IW=1,NASWBAND
      PAOD(JL,IW)=PAER_TAU(JL,JK,IW)+PAOD(JL,IW)
      PSSA(JL,IW)=PAER_SSA(JL,JK,IW)*PAER_TAU(JL,JK,IW)+PSSA(JL,IW)
      PASY(JL,IW)=PAER_ASYM(JL,JK,IW)*PAER_TAU(JL,JK,IW)*PAER_SSA(JL,JK,IW)+PASY(JL,IW)
    END DO

    DO IW=1,16
      PAOD_LW(JL,IW)=PAOD_LW(JL,IW)+PAER_TAU_LW(JL,JK,IW)
    END DO

  END DO
END DO

DO JL = KIDIA,KFDIA
  DO IW=1,NASWBAND
    IF(PSSA(JL,IW)>0._JPRB) THEN
      PASY(JL,IW) = PASY(JL,IW)/PSSA(JL,IW) ! AVERAGE
    ENDIF
    IF(PAOD(JL,IW)>0._JPRB) THEN
      PSSA(JL,IW) = PSSA(JL,IW)/PAOD(JL,IW) ! AOD AVERAGE
    ENDIF
  END DO
END DO

!*         6.3     Fill selective aerosol OD fields in structure as available in IFS-AER
!                  ---------------------------------------------------------------------

IF(MOD(NSTEP,NRADFR) == 0) THEN ! Use computed values
  DO JWAVL=1,NAERO_WVL_DIAG
    DO JL=KIDIA,KFDIA
      IF (YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_AOD) THEN
        PAERO_WVL_DIAG(JL,JWAVL,JPAERO_WVL_AOD)    = ZAOD_DIAG(JL,JWAVL)
      ENDIF
      IF (YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_AODABS) THEN
        !PAERO_WVL_DIAG(JL,JWAVL,JPAERO_WVL_AODABS) = ZABS_DIAG(JL,JWAVL)! 0.0
        PAERO_WVL_DIAG(JL,JWAVL,JPAERO_WVL_AODABS) = ZAOD_DIAG(JL,JWAVL)*(1._JPRB-ZSSA_DIAG(JL,JWAVL))! absorption
      ENDIF
      IF (YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_AODFM) THEN! not implemented yet
        PAERO_WVL_DIAG(JL,JWAVL,JPAERO_WVL_AODFM)  = 0._JPRB!PFAOD(JL,JWAVL)! 0.0
      ENDIF
      IF (YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_SSA) THEN
        PAERO_WVL_DIAG(JL,JWAVL,JPAERO_WVL_SSA)    = ZSSA_DIAG(JL,JWAVL)
      ENDIF
      IF (YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_ASSIMETRY) THEN
        PAERO_WVL_DIAG(JL,JWAVL,JPAERO_WVL_ASSIMETRY) = ZASY_DIAG(JL,JWAVL)
      ENDIF
    ENDDO
  ENDDO
ELSE                            ! Use stored values (PAERO_WVL_DIAG get corrupted - see ticket OIFS-668)

  IF (YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_AOD) THEN
    PAERO_WVL_DIAG(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, JPAERO_WVL_AOD) = PGFL(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, YAEROUT(6)%MP)
  ENDIF
  IF (YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_AODABS) THEN
    PAERO_WVL_DIAG(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, JPAERO_WVL_AODABS) = PGFL(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, YAEROUT(7)%MP)
  ENDIF
  IF (YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_AODFM) THEN
    PAERO_WVL_DIAG(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, JPAERO_WVL_AODFM) = PGFL(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, YAEROUT(8)%MP)
  ENDIF
  IF (YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_SSA) THEN
    PAERO_WVL_DIAG(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, JPAERO_WVL_SSA) = PGFL(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, YAEROUT(9)%MP)
  ENDIF
  IF (YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_ASSIMETRY) THEN
    PAERO_WVL_DIAG(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, JPAERO_WVL_ASSIMETRY) = PGFL(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, YAEROUT(10)%MP)
  ENDIF  
ENDIF

!*
!* VH - Requires full checking by expert if this makes sense - FIXME
!*
DO JL=KIDIA,KFDIA
  PODTO(JL)    =PTAUS_AER(JL,KLEV,1,1)
  PODTO469(JL) =PTAUS_AER(JL,KLEV,2,1)
  PODTO670(JL) =PTAUS_AER(JL,KLEV,3,1)
  PODTO865(JL) =PTAUS_AER(JL,KLEV,4,1)
  PODTO1240(JL)=PTAUS_AER(JL,KLEV,5,1)
ENDDO


CALL HAMM7_DIAG_PM &
 &( YDMODEL, KIDIA  , KFDIA  , KLON   , KLEV , NAERO ,&
 &  PAEROP, &
 &  PAEPM1, PAEPM25, PAEPM10, &
 &  ZRHO, ZM6DRY, ZM6RP, ZRHOP) 
!*         7.     STORE IN AEROUTs
!                 ------------------------------

! LIFSMIN (T if running minimisation) and LIFSTRAJ (T if running high
! resolution trajectory integration) are both assimilation flags
IF(.NOT.LIFSMIN  .AND. .NOT.LIFSTRAJ) THEN

  !** YAEROUT(1) : RADIATIVE PROPERTIES
  
  ! AOD/SSA/ASY of one internal (short) wavelengths - 533 nm
  PGFL(KIDIA:KFDIA,1,YAEROUT(1)%MP)=PAOD(KIDIA:KFDIA,10)
  PGFL(KIDIA:KFDIA,2,YAEROUT(1)%MP)=PSSA(KIDIA:KFDIA,10)
  PGFL(KIDIA:KFDIA,3,YAEROUT(1)%MP)=PASY(KIDIA:KFDIA,10)

  ! AOD of 14 internal (short) wavelengths
  PGFL(KIDIA:KFDIA,4:17,YAEROUT(1)%MP)=PAOD(KIDIA:KFDIA,1:14)

  ! AOD of 16 internal (long) wavelengths
  PGFL(KIDIA:KFDIA,18:33,YAEROUT(1)%MP)= PAOD_LW(KIDIA:KFDIA,1:16)

  !** YAEROUT(2) : DRY DEPOSITION

  DO JN=1,NAEROCOMP
    PGFL(KIDIA:KFDIA,ind_oifs_ham%ind_mass_OIFS(JN),YAEROUT(2)%MP)=ZDDEPFLUX(KIDIA:KFDIA,ind_oifs_ham%IND_mass_HAM(JN))
  END DO
  DO JN=1,NCLASS
    PGFL(KIDIA:KFDIA,ind_oifs_ham%ind_class_OIFS(JN),YAEROUT(2)%MP)=ZDDEPFLUX(KIDIA:KFDIA,ind_oifs_ham%IND_class_HAM(JN))
  END DO
  
  !** YAEROUT(3) : COLUMN INTEGRATED WET-DEPOSITION

  DO JN=1,NACTAERO
    PGFL(KIDIA:KFDIA,KAERO(JN),YAEROUT(3)%MP)  = WDEPOUT_2D(KIDIA:KFDIA,KAERO(JN))
  END DO
  
  !** YAEROUT(4) : SEDIMENTATION

  DO JN=1,NAEROCOMP
    PGFL(KIDIA:KFDIA,ind_oifs_ham%ind_mass_OIFS(JN),YAEROUT(4)%MP)  = ZSEDIFLUXSURF(KIDIA:KFDIA,ind_oifs_ham%IND_mass_HAM(JN))
  END DO
  DO JN=1,NCLASS
    PGFL(KIDIA:KFDIA,ind_oifs_ham%ind_class_OIFS(JN),YAEROUT(4)%MP) = ZSEDIFLUXSURF(KIDIA:KFDIA,ind_oifs_ham%IND_class_HAM(JN))
  END DO

  !** YAEROUT(5) : NET AEROSOLS FLUXES (EMISSIONS - ???) ; level index of top of boundary layer ; boundary layer height
  
  DO JN=1,NACTAERO
    PGFL(KIDIA:KFDIA,KAERO(JN),YAEROUT(5)%MP)  = PAERSRC(KIDIA:KFDIA,KAERO(JN)) - PCFLX(KIDIA:KFDIA,KAERO(JN))* ZDPG(KIDIA:KFDIA,KLEV)
  END DO
  PGFL(KIDIA:KFDIA,NACTAERO+2,YAEROUT(5)%MP)  = ZBLHIDX(KIDIA:KFDIA)
  PGFL(KIDIA:KFDIA,NACTAERO+3,YAEROUT(5)%MP)  = PBLH(KIDIA:KFDIA)

  !** YAEROUT(6:10) : Store all requested AOP at selected (diagnostic) wavelengths
  
  IF(MOD(NSTEP,NRADFR) == 0) THEN
    IF (YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_AOD) THEN
      PGFL(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, YAEROUT(6)%MP) = PAERO_WVL_DIAG(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, JPAERO_WVL_AOD)
    ENDIF
    IF (YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_AODABS) THEN
      PGFL(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, YAEROUT(7)%MP) = PAERO_WVL_DIAG(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, JPAERO_WVL_AODABS)
    ENDIF
    IF (YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_AODFM) THEN
      PGFL(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, YAEROUT(8)%MP) = PAERO_WVL_DIAG(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, JPAERO_WVL_AODFM)
    ENDIF
    IF (YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_SSA) THEN
      PGFL(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, YAEROUT(9)%MP) = PAERO_WVL_DIAG(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, JPAERO_WVL_SSA)
    ENDIF
    IF (YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_ASSIMETRY) THEN
      PGFL(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, YAEROUT(10)%MP) = PAERO_WVL_DIAG(KIDIA:KFDIA, 1:NAERO_WVL_DIAG, JPAERO_WVL_ASSIMETRY)
    ENDIF
  ENDIF
  
  !** YAEROUT(11) : Total column mass and number concentration
  
  DO JN=1,NAEROCOMP
    JO=ind_oifs_ham%ind_mass_OIFS(JN)  ! JO -> index context OIFS 
    JH=ind_oifs_ham%IND_mass_HAM(JN)   ! JH -> index context HAM 
    JY=YAEROUT(11)%MP
    ZTMP(KIDIA:KFDIA)=0.0_JPRB
    DO JK=1,KLEV
      ZTMP(KIDIA:KFDIA)= ZTMP(KIDIA:KFDIA) + (ZXTM1(KIDIA:KFDIA,JK,JH)+(ZXTTE(KIDIA:KFDIA,JK,JH)*TIME_STEP_LEN)) * ZDPG(KIDIA:KFDIA,JK)
    END DO
    PGFL(KIDIA:KFDIA,JO,JY) = ZTMP(KIDIA:KFDIA)
  END DO

  DO JN=1,NCLASS
    JO=ind_oifs_ham%ind_class_OIFS(JN)  ! JO -> index context OIFS 
    JH=ind_oifs_ham%IND_class_HAM(JN)   ! JH -> index context HAM 
    JY=YAEROUT(11)%MP
    ZTMP(KIDIA:KFDIA)=0.0_JPRB
    DO JK=1,KLEV
      ZTMP(KIDIA:KFDIA) = ZTMP(KIDIA:KFDIA) + (ZXTM1(KIDIA:KFDIA,JK,JH)+(ZXTTE(KIDIA:KFDIA,JK,JH)*TIME_STEP_LEN)) * ZDPG(KIDIA:KFDIA,JK)
    END DO
    PGFL(KIDIA:KFDIA,JO,JY) = ZTMP(KIDIA:KFDIA)
  END DO

  !** YAEROUT(12) : mass and number tendency
  ! kg/kg -> kg/m2 N/kg-> N/m2
  DO JN=1,NAEROCOMP    !ntrac!NACTAERO   
    JO=ind_oifs_ham%ind_mass_OIFS(JN)  ! JO -> index context OIFS 
    JH=ind_oifs_ham%IND_mass_HAM(JN)   ! JH -> index context HAM 
    JY=YAEROUT(12)%MP
    ZTMP(KIDIA:KFDIA)=0.0_JPRB
    DO JK=1,KLEV
      ZTMP(KIDIA:KFDIA) = ZTMP(KIDIA:KFDIA) + ZXTTE(KIDIA:KFDIA,JK,JH)
    END DO
    PGFL(KIDIA:KFDIA,JO,JY) = ZTMP(KIDIA:KFDIA)
  END DO

  DO JN=1,NCLASS
    JO=ind_oifs_ham%ind_class_OIFS(JN)  ! JO -> index context OIFS 
    JH=ind_oifs_ham%IND_class_HAM(JN)   ! JH -> index context HAM 
    JY=YAEROUT(12)%MP
    ZTMP(KIDIA:KFDIA)=0.0_JPRB
    DO JK=1,KLEV
      ZTMP(KIDIA:KFDIA) = ZTMP(KIDIA:KFDIA) + ZXTTE(KIDIA:KFDIA,JK,JH)
    END DO
    PGFL(KIDIA:KFDIA,JO,JY) = ZTMP(KIDIA:KFDIA)
  END DO

  !** YAEROUT(13) : Surface fluxes of tracers (not emissions)
  
  DO JN=1,NACTAERO
    !ZTMP(KIDIA:KFDIA)=0.0_JPRB
    !DO JK=1,KLEV
    !  ZTMP(KIDIA:KFDIA)=ZTMP(KIDIA:KFDIA)+PCEN(KIDIA:KFDIA,JK,KAERO(JN))
    !END DO
    !PGFL(KIDIA:KFDIA,KAERO(JN),YAEROUT(13)%MP)  = ZTMP(KIDIA:KFDIA)
    PGFL(KIDIA:KFDIA,KAERO(JN),YGFL%YAEROUT(13)%MP)= - PCFLX(KIDIA:KFDIA,KAERO(JN))* ZDPG(KIDIA:KFDIA,KLEV)
  END DO

  !** YAEROUT(14) : Total column tracer/number PREVIOUS (before call to M7) concentration
  
  DO JN=1,NACTAERO
    ZTMP(KIDIA:KFDIA)=0.0_JPRB
    DO JK=1,KLEV
      ZTMP(KIDIA:KFDIA)=ZTMP(KIDIA:KFDIA)+ZCEN(KIDIA:KFDIA,JK,KAERO(JN))
    END DO
    PGFL(KIDIA:KFDIA,KAERO(JN),YAEROUT(14)%MP)  = ZTMP(KIDIA:KFDIA)
  END DO

  !** YAEROUT(15) : Total column UPDATED TENDENCIES tracer
  
  DO JN=1,NACTAERO
    ZTMP(KIDIA:KFDIA)=0.0_JPRB
    DO JK=1,KLEV
      ZTMP(KIDIA:KFDIA)=ZTMP(KIDIA:KFDIA)+PTENC(KIDIA:KFDIA,JK,KAERO(JN))
    END DO
    PGFL(KIDIA:KFDIA,KAERO(JN),YAEROUT(15)%MP)  = ZTMP(KIDIA:KFDIA)
  END DO

  !** YAEROUT(16) : M7 mass and number mixing ratio at surface
  
  DO JN=1,NAEROCOMP    !ntrac!NACTAERO   
    JO=ind_oifs_ham%ind_mass_OIFS(JN)  ! JO -> index context OIFS 
    JH=ind_oifs_ham%IND_mass_HAM(JN)   ! JH -> index context HAM 
    JY=YAEROUT(16)%MP
    ZTMP(KIDIA:KFDIA)= ZXTM1(KIDIA:KFDIA,KLEV,JH)+(ZXTTE(KIDIA:KFDIA,KLEV,JH)*TIME_STEP_LEN)!*ZDPG(KIDIA:KFDIA,KLEV)
    PGFL(KIDIA:KFDIA,JO,JY) = ZTMP(KIDIA:KFDIA)
  END DO

  DO JN=1,NCLASS
    JO=ind_oifs_ham%ind_class_OIFS(JN)  ! JO -> index context OIFS 
    JH=ind_oifs_ham%IND_class_HAM(JN)   ! JH -> index context HAM 
    JY=YAEROUT(16)%MP
    ZTMP(KIDIA:KFDIA) =  ZXTM1(KIDIA:KFDIA,KLEV,JH)+(ZXTTE(KIDIA:KFDIA,KLEV,JH)*TIME_STEP_LEN)!*PRHO(KIDIA:KFDIA,KLEV)
    PGFL(KIDIA:KFDIA,JO,JY) = ZTMP(KIDIA:KFDIA)
  END DO
  
  !** YAEROUT(17-18) : IN-CLOUD & BELOW CLOUD WET DEPOSITION
  
  DO JN=1,NACTAERO
    PGFL(KIDIA:KFDIA,KAERO(JN),YAEROUT(17)%MP) = WDEPOUT_IC_2D(KIDIA:KFDIA,KAERO(JN))
    PGFL(KIDIA:KFDIA,KAERO(JN),YAEROUT(18)%MP) = WDEPOUT_BC_2D(KIDIA:KFDIA,KAERO(JN))
  END DO

  !** YAEROUT(19-20) : 
  
  PGFL(KIDIA:KFDIA,KLEV,YAEROUT(19)%MP)   = ZXTTE(KIDIA:KFDIA,KLEV,3)      ! tendency SS CS ham after update surface
  PGFL(KIDIA:KFDIA,KLEV-1,YAEROUT(20)%MP) = ZTENCIH(KIDIA:KFDIA,KLEV,17)   ! tendency SS CS ham before update surface

  !** YAEROUT(21) :  height of each level top
  
  DO JK=1,KLEV
    PGFL(KIDIA:KFDIA,JK,YAEROUT(21)%MP) = (PGEOH(KIDIA:KFDIA,JK-1)-PGEOH(KIDIA:KFDIA,KLEV))*ZRG
  END DO
  
  !** YAEROUT(22-26) : Nucleation diagnostics
  
  PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(22)%MP) = ZOUT_dnuc(KIDIA:KFDIA,1:KLEV,1)  ! NS-mass production kg/s from nucleatiom (limited by amount of SO4)
  PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(23)%MP) = ZOUT_dnuc(KIDIA:KFDIA,1:KLEV,2)  ! NS-number production rate #/s from nucleatiom (limited by amount of SO4 vs original #/s)
  PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(24)%MP) = ZOUT_dnuc(KIDIA:KFDIA,1:KLEV,3)  ! original #/s  from H2SO4/H2O nucleatiom
  PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(25)%MP) = ZOUT_dnuc(KIDIA:KFDIA,1:KLEV,4)  ! original #/s from organic nucleatiom
  PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(26)%MP) = ZOUT_dnuc(KIDIA:KFDIA,1:KLEV,5)  ! original #/s sum from organic and vehkamäki nucleation schemes

  !** YAEROUT(27) : M7 GAS MIXING RATIOS at SURFACE

  DO JN=1,SUBM_NGASSPEC
    PGFL(KIDIA:KFDIA,JN,YAEROUT(27)%MP)=zxtm1(KIDIA:KFDIA,KLEV,ind_oifs_ham%ind_gas_HAM(JN))
  END DO

  !** YAEROUT(28) :  EMISSIONS
  
  DO JN=1,NACTAERO
    PGFL(KIDIA:KFDIA, JN, YAEROUT(28)%MP) = PAERSRC(KIDIA:KFDIA,KAERO(JN))
  END DO

  !** YAEROUT(29-30) : LIDAR RATIO (sr) AND LAYER OPTICAL DEPTH AT DIAGNOSTIC WAVELENGTH 2
  IF(MOD(NSTEP,NRADFR) == 0 .AND. NAEROOPT == 2) THEN
    IF(NAERO_WVL_DIAG > 1 .AND. YGFL%NAERO_WVL_DIAG_TYPES > 0) THEN
      PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(29)%MP) = ZAER_LidarR_DIAG(KIDIA:KFDIA,1:KLEV,2)
      PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(30)%MP) = ZAER_TAU_DIAG(KIDIA:KFDIA,1:KLEV,2)
    ENDIF
  ENDIF
  
  DO JN=1,NTRAC
    PGFL(KIDIA:KFDIA,JN,YAEROUT(39)%MP)=ZXTEMS(KIDIA:KFDIA,JN)
  END DO

  !** YAEROUT(31) : --
  !** YAEROUT(32) : --
  !** YAEROUT(33) : --
  !** YAEROUT(34) : --
  !** YAEROUT(35) : --
  !** YAEROUT(36) : --
  !** YAEROUT(37) : --
  !** YAEROUT(38) : --
  !** YAEROUT(39) : --
  !** YAEROUT(46) 
  !   !---------------DUST FIELDS DIAGNOSTICS-----------------------
  ! !** YAEROUT(47) : --
  ! !3D Dust DU_AI in KG/M3
  !  JO=  20 ! JO -> index context OIFS``
  !  JH=  20  ! JH -> index context HAM
  !  JY=YAEROUT(47)%MP
  !  DO JK=1,KLEV
  !    PGFL(KIDIA:KFDIA,JK,JY)=  (ZXTM1(KIDIA:KFDIA,JK,JH)+(ZXTTE(KIDIA:KFDIA,JK,JH)*TIME_STEP_LEN)) * ZDPG(KIDIA:KFDIA,JK)
  !  END DO
   

  !  ! ** YAEROUT(48) : -- 
  ! !DU_CI
  !  JO=  26 ! JO -> index context OIFS``
  !  JH=  21  ! JH -> index context HAM
  !  JY=YAEROUT(48)%MP
  !  DO JK=1,KLEV
  !    PGFL(KIDIA:KFDIA,JK,JY)=  (ZXTM1(KIDIA:KFDIA,JK,JH)+(ZXTTE(KIDIA:KFDIA,JK,JH)*TIME_STEP_LEN)) * ZDPG(KIDIA:KFDIA,JK)
  !  END DO

  ! ! **YAEROUT(49)
  !  !DU_AI NUM
  !  JO=  19 ! JO -> index context OIFS``
  !  JH=  27  ! JH -> index context HAM
  !  JY=YAEROUT(49)%MP
  !  DO JK=1,KLEV
  !    PGFL(KIDIA:KFDIA,JK,JY)=  ZXTM1(KIDIA:KFDIA,JK,JH)+(ZXTTE(KIDIA:KFDIA,JK,JH)*TIME_STEP_LEN)
  !  END DO
  
  !  ! **YAEROUT(50)
  !  !DU_CI NUM
  !  JO=  25 ! JO -> index context OIFS``
  !  JH=  28  ! JH -> index context HAM
  !  JY=YAEROUT(50)%MP
  !  DO JK=1,KLEV
  !    PGFL(KIDIA:KFDIA,JK,JY)=  ZXTM1(KIDIA:KFDIA,JK,JH)+(ZXTTE(KIDIA:KFDIA,JK,JH)*TIME_STEP_LEN)
  !  END DO
   
!Commented because of overlap      DO JGAS=1,SUBM_NGASSPEC
!Commented because of overlap    !!! BUG: JL? JK? !!!    PGFL(KIDIA:KFDIA,JGAS,YAEROUT(28+JGAS)%MP)=ZCEN(JL,JK,KCHEM(ind_oifs_ham%ind_gas_OIFS(JGAS)))
!Commented because of overlap        PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(29+JGAS)%MP)= PTENC(KIDIA:KFDIA,1:KLEV,KCHEM(ind_oifs_ham%ind_gas_OIFS(JGAS)))
!Commented because of overlap      END DO
!Commented because of overlap    
!Commented because of overlap    DO IMODE=1,NMOD
!Commented because of overlap      PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(30+IMODE)%MP)=RW_MODE(IMODE)%d2(KIDIA:KFDIA,1:KLEV) ! m
!Commented because of overlap    ENDDO

  !IF (TRIM(CHEM_SCHEME)=="SimChem")THEN
  !  PGFL(KIDIA:KFDIA,KLEV,YAEROUT(40)%MP)   = ZFSO2(KIDIA:KFDIA)           ! tendency SS CS ham after update surface
  !  PGFL(KIDIA:KFDIA,KLEV,YAEROUT(41)%MP)   = ZFSO4(KIDIA:KFDIA)           ! tendency SS CS ham after update surfac
  !  PGFL(KIDIA:KFDIA,KLEV,YAEROUT(42)%MP)   = ZFSO4_AQ(KIDIA:KFDIA)        ! tendency SS CS ham after update surface
  !  PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(43)%MP) = ZTSO4(KIDIA:KFDIA,1:KLEV,1)  ! tendency SS CS ham after update surface
  !  PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(44)%MP) = ZTSO4_AQ(KIDIA:KFDIA,1:KLEV) ! tendency SS CS ham after update surface
  !  PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(45)%MP) = ZTSO2(KIDIA:KFDIA,1:KLEV)    ! tendency SS CS ham after update surface
  !END IF

ENDIF

!------------------------------------------------------------------------------
!*         9.       RELEASE LOCAL MEMORY
!                   --------------------
IF (ALLOCATED(ZAERNGT)) DEALLOCATE( ZAERNGT )
IF (ALLOCATED(ZWPDF ) ) DEALLOCATE( ZWPDF )
IF (ALLOCATED(ZW ) )    DEALLOCATE( ZW )
IF (ALLOCATED(ZRC) )    DEALLOCATE( ZRC )
IF (ALLOCATED(ZSMAX))   DEALLOCATE( ZSMAX )

DO IMODE=1,NMOD
  IF (ASSOCIATED(RW_MODE(IMODE)%D2)) DEALLOCATE(RW_MODE(IMODE)%D2)
  IF (ASSOCIATED(DENS_MODE(IMODE)%D2)) DEALLOCATE(DENS_MODE(IMODE)%D2)
ENDDO

DO IMODE=1,NMOD
  IF (SIZECLASS(IMODE)%LSOLUBLE) THEN
    IF (ASSOCIATED(RWD_MODE(IMODE)%D2)) DEALLOCATE(RWD_MODE(IMODE)%D2)
    IF (ASSOCIATED(H2O_MODE(IMODE)%D2)) DEALLOCATE(H2O_MODE(IMODE)%D2)
  END IF
ENDDO

END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('HAMM7_INTERFACE',1,ZHOOK_HANDLE)

END SUBROUTINE HAMM7_INTERFACE
