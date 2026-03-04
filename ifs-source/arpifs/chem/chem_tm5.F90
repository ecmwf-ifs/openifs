! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

 SUBROUTINE CHEM_TM5&
 &   (YDVAB,YDDIMV,YDMODEL,KSTEP, KIDIA  , KFDIA , KLON, KLEV ,KVCLIS,KAERO,&
 &    PTSTEP ,PDELP, PRS1, PRSF1, PGEOH, PQP, PTP,&
 &    PLP, PIP, PAP,  PPTROPO, PALB, PWND, PLSM, PCSZA, PGELAT,&
 &    PGELAM, PGEMU,  PKOZO, KLEVTROP, PDV, PCEN , PTENC1, PBUDR, PBUDJ, PBUDX, POUT,&
 &    PAEROP,PWETDIAM, PWETVOL,PND,PAERAOT, PAERAAOT, PAERASY,PSOGTOSOA, & 
 &    PCHEM2GHG) 
 

!**   DESCRIPTION 
!     ----------
!
!   TM5 routine for IFS chemistry 
!
!
!
!**   INTERFACE.
!     ----------
!          *CHEM_TM5* IS CALLED FROM *CHEM_MAIN*.


! INPUTS:
! -------
! KSTEP                       : Time step 
! KIDIA                       : Start of Array  
! KFDIA                       : End  of Array 
! KLON                        : Length of Arrays 
! KLEV                        : Number of Levels
! KVCLIS                      : Number Cariolle chemistry coefficinets 
! KAERO                       : Number of aerosol fields 
! PTSTEP                      : Time step in seconds 
! PDELP(KLON,KLEV)            : PRESSURE DELTA in PRESSURE UNITES      (Pa)
! PRS1(KLON,0:KLEV)           : HALF-LEVEL PRESSURE           (Pa)
! PRSF1(KLON,KLEV)            : FULL-LEVEL PRESSURE           (Pa)
! PGEOH(KLON,0:KLEV)          : GEOPOTENTIAL                 (m*m/s*s)
! PQP     (KLON,KLEV)         : SPECIFIC HUMIDITY             (kg/kg)
! PTP     (KLON,KLEV)         : TEMPERATURE                   (K)
! PLP     (KLON,KLEV)         : LCWC                          (kg/kg)
! PIP     (KLON,KLEV)         : ICWC                          (kg/kg)
! PAP     (KLON,KLEV)         : CLOUD FRACTION                0..1  
! PPTROPO  (KLON)              : Tropopause pressure           (Pa)
! PLSM    (KLON)              : land-sea-mask
! PALB(KLON)                  : Surface albedo
! PWND(KLON)                  : Surface wind
! PLSM(KLON)                  : Land Sea Mask albedo
! PCSZA(KLON)                 : COS of Solar Zenit Angle
! PGELAM(KLON)                : LONGITUDE (RADIANS)
! PGELAT(KLON)                : LATITUDE (RADIANS) 
! PGEMU(KLON)                 : SINE OF LATITUDE
! PDV(KLON,NCHEM_DV)          : DEPOSITION VELOCITIES in m/s (positive) 
! PCEN(KLON,KLEV,NCHEM)       : CONCENTRATION OF TRACERS           (kg/kg)
! PKOZO(KLON,KLEV,KVCLIS)     : PHOTOCHEMICAL COEFFICIENTS COMPUTED FROM A 2D PHOTOCHEMICAL MODEL (KVCLIS=8)!
! PAEROP(KLON,KLEV,KAERO)     : Aerosol concentrations  (kg/kg) - Note that fields are only non-zero if NACTAERO > 0
! PWETDIAM(KLON,KLEV,NMODE)   : Glomap geometric mean wet diameter per mode
! PWETVOL(KLON,KLEV,NMODE)    : Glomap avg wet volume of size mode (m3)
! PAERAOT(KLON,KLEV,6)        : Glomap extinction AOD per model level at 6 wavelengths
! PAERAAOT(KLON,KLEV,6)       : Glomap absorption AOD per model levelat 6 wavelengths
! PAERASY(KLON,KLEV,6)        : Glomap asymetry factor
!
! OUTPUTS:
! -------
! PTENC1(KLON,KLEV,NCHEM)     : TENDENCY OF CONCENTRATION OF TRACERS BECAUSE OF CHEMISTRY (kg/kg s-1)
! PBUDR (KLON,KLEV,NCHEM)     : TENDENCIES DUE TO GAS-PHASE REACTIONS WITH OH (kg/kg/s)
! PBUDJ (KLON,KLEV,NPHOTO+NSOA_BUDG)     : TENDENCIES (loss) DUE TO PHOTOLYSIS (kg/kg/s)
! PBUDX(KLON,KLEV,NBUD_EXTRA) : Extra chemical TENDENCIES (kg/kg/s)
! POUT (KLON,KLEV,5)         : additiional output, e.g. UBC contribution , Photolysis rates O3 , NO2, tau for output 
! PSOGTOSOA(KLON,KLEV,2)     : SOG to SOA conversion tendency
! PCHEM2GHG(KLON,KLEV,NCHEM2GHG): Information from chemistry to GHG. 
!                                   1. atmospheric CH4 loss rate              [s-1] 
!                                   2. tropospheric CO2 production tendency due to CO oxidation [kg CO2/kg/s] 
!
! LOCAL:
! -------
!
! ZCVM0(KLON,NCHEM)           : initial volume ratios OF TRACERS           (molec/cm3)
! ZCVM1(KLON,NCHEM)           :         volume ratios OF TRACERS ; corrected for NOx  (molec/cm3)
! ZCVM (KLON,NCHEM+3)         : final   volume ratios OF TRACERS           (molec/cm3)
! ZRJ(KLON,NPHOTO)            : photolysis reaction rates
! ZRR(KLON,NREAC)             : reaction rates
! ZCC (KPROMA,KFLEV)          : OVERHEAD CLOUD COVER 
! ZTCTAUC (KPROMA)            : Total optical depth
! ZPCO3 (KPROMA,KFLEV)        : O3 column        
! ZHPLUS(KPROMA)              : concentration H+; optionally to be outputted...        
!
!  ZTAUS_AER, ZTAUA_AER, ZPMAER : Aerosol optical properties
!  ZTAUA_CLD,ZTAUS_CLD,ZPMCLD   : Cloud optical properties
!  ZSAD_AER,ZSAD_CLD,ZSAD_ICE   : Surface area densities, available for diagnostics purposes 


! ZCR2(KLON,NPHOTO)           : accumulators for budget calculations: photolysis
! ZCR3(KLON,NREAC)            : accumulators for budget calculations: gas phase chemistry
!
!     AUTHOR.
!     -------
!        JOHANNES FLEMMING  *ECMWF*
!        VINCENT HUIJNEN    *KNMI*
!        This source code originates from TM5 (*The TM5-community*)
!        ORIGINAL : 2009-07-22

!     MODIFICATIONS.
!     --------------
!        TM5-implementation of RnPb   scheme               : 2009-09-07
!        TM5-implementation of TM5-chemistry code structure: 2009-09-08
!        TM5-implementation of CBM-IV scheme               : 2009-09-23
!        Further code cleaning                             : 2011-10-25
!        K. Yessad (July 2014): Move some variables.
!-----------------------------------------------------------------------



USE TYPE_MODEL , ONLY : MODEL
USE YOMVERT  , ONLY : TVAB
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN    ,ONLY : NULERR, NULOUT
! NCHEM : number of chemical species
! YCHEM : Data structure with Chemistry meta data 
USE TM5_CHEM_MODULE, ONLY : NREAC, IO3, ICH4, IHNO3, &
 & XMH2O, IACID, IAIR, IH2O, IO3S, INO, INO2, IPSC, &
 & INO3_A, ISO4, INH3 , INH4, IMSA, IHO2, IOH, ICO, NBUD_EXTRA,&
 & ISOG1,ISOG2A,ISOG2B, NSOA_BUDG, &
 & KCOOH, KHO2L, KHO2_AER, &
 & NCHEM2GHG
USE TM5_PHOTOLYSIS , ONLY : NPHOTO,NBANDS_TROP,NGRID, JNO2
USE YOMCST   , ONLY : RD, RMD, RG , RPI, RMCO2
USE YOMRIP0  , ONLY : NINDAT 

! TM5 chemistry ...
USE TM5_KPP_Parameters, ONLY : NREACT, NVAR
USE TM5_KPP_global    , ONLY : RTOL,ATOL, ROUNDOFF_STORE

! General KPP parameters
USE CIFS_KPP_IntParam  , ONLY : HMIN,HSTART,RTOLS_G,IAUTONOM,IROSMETH, VMR_BAD_LARGE

! Glomap specifics
USE UKCA_MODE_SETUP, ONLY: NMODES

!!   use ieee_arithmetic 
!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVAB)        ,INTENT(IN) :: YDVAB
TYPE(TDIMV)       ,INTENT(IN) :: YDDIMV
TYPE(MODEL)       ,INTENT(INOUT):: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN) :: KSTEP, KIDIA , KFDIA , KLON , KLEV, KVCLIS, KAERO
REAL(KIND=JPRB)   ,INTENT(IN) :: PTSTEP
REAL(KIND=JPRB)   ,INTENT(IN) :: PDELP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PRSF1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PRS1(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PGEOH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PQP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PTP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PLP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PIP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PAP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT):: PTENC1(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NCHEM)
REAL(KIND=JPRB)   ,INTENT(IN) :: PCEN(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NCHEM) 
REAL(KIND=JPRB)   ,INTENT(IN) :: PCSZA(KLON)                  
REAL(KIND=JPRB)   ,INTENT(IN) :: PPTROPO(KLON)                   
REAL(KIND=JPRB)   ,INTENT(IN) :: PALB(KLON)                   
REAL(KIND=JPRB)   ,INTENT(IN) :: PWND(KLON)                   
REAL(KIND=JPRB)   ,INTENT(IN) :: PLSM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PGELAT(KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PGELAM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PGEMU(KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PDV(KLON,YDMODEL%YRML_GCONF%YGFL%NCHEM_DV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PKOZO(KLON,KLEV,KVCLIS)
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVTROP(KLON)                  
REAL(KIND=JPRB)   ,INTENT(OUT):: PBUDJ(KLON,KLEV,NPHOTO+NSOA_BUDG) 
REAL(KIND=JPRB)   ,INTENT(OUT):: PBUDR(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NCHEM) 
REAL(KIND=JPRB)   ,INTENT(OUT):: PBUDX(KLON,KLEV,NBUD_EXTRA)  
REAL(KIND=JPRB)   ,INTENT(OUT):: POUT(KLON,KLEV,5)  
REAL(KIND=JPRB)   ,INTENT(IN) :: PAEROP(KLON,KLEV,KAERO)
REAL(KIND=JPRB)   ,INTENT(IN) :: PWETDIAM(KLON,KLEV,NMODES)
REAL(KIND=JPRB)   ,INTENT(IN) :: PWETVOL(KLON,KLEV,NMODES)
REAL(KIND=JPRB)   ,INTENT(IN) :: PND(KLON,KLEV,NMODES)
REAL(KIND=JPRB)   ,INTENT(IN) :: PAERAOT(KLON,KLEV,6)
REAL(KIND=JPRB)   ,INTENT(IN) :: PAERAAOT(KLON,KLEV,6)
REAL(KIND=JPRB)   ,INTENT(IN) :: PAERASY(KLON,KLEV,6)
REAL(KIND=JPRB)   ,INTENT(OUT):: PSOGTOSOA(KLON,KLEV,2)
REAL(KIND=JPRB)   ,INTENT(OUT):: PCHEM2GHG(KLON,KLEV,NCHEM2GHG)

!-----------------------------------------------------------------------

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! * Lat /Lon
REAL(KIND=JPRB) , DIMENSION(KLON)   :: ZLAT
REAL(KIND=JPRB) , DIMENSION(KLON)   :: ZLON

! * counters
INTEGER(KIND=JPIM) :: JK, JL, JT, JLEV
!INTEGER(KIND=JPIM) :: IDRYDEP
INTEGER(KIND=JPIM) :: ISSO4, ISNH4, ISNO3
INTEGER(KIND=JPIM) :: ISSOA1,ISSOA2

! * chemical data 
REAL(KIND=JPRB) , DIMENSION(KLON,YDMODEL%YRML_GCONF%YGFL%NCHEM+3)   :: ZCVM
REAL(KIND=JPRB) , DIMENSION(KLON,YDMODEL%YRML_GCONF%YGFL%NCHEM)     :: ZCVM0, ZCVM1 
REAL(KIND=JPRB) , DIMENSION(KLON)           :: ZAIRDM

!REAL(KIND=JPRB)     :: ZHGT(KLON)       ! geopotential layer bottom  

! * Photolysis data: cloud/ozone info, H+ concentration
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV)   :: ZCC
REAL(KIND=JPRB) , DIMENSION(KLON)        :: ZTCTAUC 
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV)   :: ZPCO3
REAL(KIND=JPRB)                          :: ZCOLO3(KLON,0:KLEV)
REAL(KIND=JPRB) , DIMENSION(KLON)        :: ZHPLUS

REAL(KIND=JPRB) , DIMENSION(KLON,KLEV,NBANDS_TROP,NGRID) ::&
      &                                       ZTAUA_AER, ZTAUS_AER, ZPMAER   
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV)   ::   ZTAUA_CLD, ZTAUS_CLD, ZPMCLD 
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV)   ::   ZCLOUD_REFF
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV)   ::   ZSAD_AER, ZSAD_CLD, ZSAD_ICE

! * reaction rates
REAL(KIND=JPRB) , DIMENSION(KLON,NREAC)       :: ZRR 
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV,NPHOTO) :: ZRJ 
REAL(KIND=JPRB) , DIMENSION(KLON,NPHOTO)      :: ZRJ_IN 

!* Saturation pressure for water vapour
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV) :: ZQSAT, ZRHCL

!* O3S help tracer
REAL(KIND=JPRB)                          :: ZCEN_O3

! * budget accumulators
REAL(KIND=JPRB) , DIMENSION(KLON,NPHOTO) :: ZCR2 
REAL(KIND=JPRB) , DIMENSION(KLON,NREAC)  :: ZCR3 

! * EQSAM input / output parameters
REAL(KIND=JPRB)                          :: ZRH,ZCCS !,ZTR,ZWV,ZRRH
REAL(KIND=JPRB)                          :: ZNH,ZNO3,ZSO4
REAL(KIND=JPRB) , DIMENSION(4)           :: ZYEQ

! * deposition velocities, optionally provided via IFS
REAL(KIND=JPRB) , DIMENSION(KLON,YDMODEL%YRML_GCONF%YGFL%NCHEM)  :: ZVD

! *  EBI-solver variables
INTEGER(KIND=JPIM)                       :: IMAXIT,ITEREBI,IB
REAL(KIND=JPRB)                          :: ZDT_EBI

! * O3 tendency from Cariolle chemistry 
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV)   :: ZTENO3COR
REAL(KIND=JPRB)                          :: ZPMAXO3CAR,ZPMINO3CAR,ZPO3CAR_RATIO
! * tendencies from surface / stratospheric BC
INTEGER(KIND=JPIM), DIMENSION(KLON,4)     :: JLEVBC
REAL(KIND=JPRB)   , DIMENSION(KLON,3)     :: ZTENBC, ZCENBC1, ZCENBC2
INTEGER(KIND=JPIM)                        :: IMONTH, IDAY, IYEAR
INTEGER(KIND=JPIM)                        :: ILEVMIN , ILEV_CH4
INTEGER(KIND=JPIM)                        :: ILEVMIN_GAS(KLON)

!KPP related
INTEGER(KIND=JPIM),DIMENSION(20)          :: ICNTRL, ISTATUS
REAL(KIND=JPRD),   DIMENSION(20)          :: ZCNTRL, ZCNTRL_P, ZSTATE
! Troposphere...
REAL(KIND=JPRD),   DIMENSION(NREACT)      :: ZRCONST
REAL(KIND=JPRD),   DIMENSION(NVAR)        :: ZVAR
INTEGER(KIND=JPIM)                        :: IERR, IFLAG

INTEGER(KIND=JPIM), DIMENSION(3)          :: JAER_TRACER
INTEGER(KIND=JPIM), DIMENSION(2)          :: JRATE_TRACER


REAL(KIND=JPRB)                           :: ZTAU_NUDGE
REAL(KIND=JPRB)                           :: ZRGI, ZDELP, ZT0NO, ZT0NO2, ZCONC_HO2 
LOGICAL                                   :: LLCOD_TM5
! Switch for selection of aerosol parameterization for photolysis
INTEGER(KIND=JPIM), PARAMETER             :: ITAU_MACC=0
! Switch for selection of HNO3 BC parameterization. 1=updated
INTEGER(KIND=JPIM), PARAMETER             :: IMODE_HNO3BC=1

!Help variable
REAL(KIND=JPRB)                           :: ZFAC

! Variables needed for SOA computation
INTEGER(KIND=JPIM), PARAMETER             :: INSOG=3
INTEGER(KIND=JPIM), PARAMETER             :: INSOA=2
REAL(KIND=JPRB)                           :: ZORGAERO, ZRHO
REAL(KIND=JPRB), DIMENSION(INSOG)         :: ZSOG
REAL(KIND=JPRB), DIMENSION(INSOA)         :: ZSOGH, ZSOA
INTEGER(KIND=JPIM)                        :: INBDU,INBOM
REAL(KIND=JPRB)                           :: ZTSCALI,ZXLSOA,ZSOA_TMP,ZJSOA
INTEGER(KIND=JPIM), DIMENSION(INSOG)      :: JSOG_TRACER
INTEGER(KIND=JPIM), DIMENSION(INSOA)      :: JSOA_TRACER
! ------------------------------------------------------------------
#include "fcttim.func.h"
!-------------------------------------------------------------------
#include "satur.intfb.h"
#include "tm5_aerosol_info.intfb.h"
#include "tm5_boundary_ch4.intfb.h"
#include "tm5_boundary_hno3.intfb.h"
#include "tm5_calrates.intfb.h"
#include "tm5_do_ebi.intfb.h"
! #include "tm5_do_ebi_tc02b.intfb.h"
#include "tm5_eqsam.intfb.h"
#include "tm5_glomap_aerosol.intfb.h"
#include "tm5_ibud.intfb.h"
#include "tm5_macc_aerosol.intfb.h"
#include "tm5_o3s.intfb.h"
#include "tm5_photo_flux.intfb.h"
#include "tm5_rbud.intfb.h"
#include "tm5_soa.intfb.h"
#include "tm5_slingo.intfb.h"
#include "tm5_stratbc_ch4.intfb.h"
#include "tm5_sundis.intfb.h"
#include "tm5_wetchem.intfb.h"
#include "cod_op_tm5.intfb.h"
#include "o3chem.intfb.h"
! KPP code - TROP - Code already included in chem_tm5bascoe.F90...
#include "tm5_kpp_rates.intfb.h"
#include "tm5_v0_kpp_initialize.intfb.h"
#include "tm5_kpp_integrator.intfb.h"
#include "tm5_v0_kpp_update_cifs_conc.intfb.h"
#include "cifs_kpp_wlamch.intfb.h"

IF (LHOOK) CALL DR_HOOK('CHEM_TM5',0,ZHOOK_HANDLE )
ASSOCIATE(YDECLD=>YDMODEL%YRML_PHY_EC%YRECLD,YGFL=>YDMODEL%YRML_GCONF%YGFL,&
 & YDCHEM=>YDMODEL%YRML_CHEM%YRCHEM, YDERDI=>YDMODEL%YRML_PHY_RAD%YRERDI, &
 & YDEAERSRC=>YDMODEL%YRML_PHY_AER%YREAERSRC, YDEAERSNK=>YDMODEL%YRML_PHY_AER%YREAERSNK, &
 & YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM, YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO, &
 & LPHYLIN=>YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN)
ASSOCIATE(NACTAERO=>YGFL%NACTAERO, NCHEM=>YGFL%NCHEM, NCHEM_DV=>YGFL%NCHEM_DV, &
 & YCHEM=>YGFL%YCHEM, LAERCHEM=>YGFL%LAERCHEM,LAERSOA=>YDCOMPO%LAERSOA, &
 & LAERSOA_COUPLED=>YDCOMPO%LAERSOA_COUPLED, &
 & NTYPAER=>YDEAERATM%NTYPAER, &
 & LAERNITRATE => YDCOMPO%LAERNITRATE, &
 & KCHEM_SOLVE=>YDCHEM%KCHEM_SOLVE, LCHEM_0NOX=>YDCHEM%LCHEM_0NOX, LCHEM_ANACH4=>YDCHEM%LCHEM_ANACH4, &
 & LCHEM_WEAK_CH4_RELAXATION=>YDCHEM%LCHEM_WEAK_CH4_RELAXATION, &
 & LCHEM_DIAC=>YDCHEM%LCHEM_DIAC, LCHEM_JOUT=>YDCHEM%LCHEM_JOUT,LCHEM_AEROI=>YDCHEM%LCHEM_AEROI, &
 & LCHEM_REVCHEM=>YDCHEM%LCHEM_REVCHEM, REPSEC=>YDECLD%REPSEC, &
 & REPCLC=>YDERDI%REPCLC,AERO_SCHEME=>YDCOMPO%AERO_SCHEME) 
!-----------------------------------------------------------------------
! chemistry scheme name - this will later also come from external input
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Preparation for kpp - solver
!
! Set kpp parameters to default, taken from cifs_kpp_IntParam.F90 module 
! ----------------------------------------------------------------------
!VH RTOLS = RTOLS_G

RTOL(1:NVAR) = 0.5 * RTOLS_G
ATOL(1:NVAR) = 10.  !  was 1.e-16*cfactor before v3s04 or v3d06

ICNTRL(:) = 0_JPIM
ICNTRL(1) = IAUTONOM
! Select Integrator
ICNTRL(3) = IROSMETH
ICNTRL(4) = 200 ! set max. no of steps ?
ICNTRL(7) = 1 ! Currently no adjoint

ZCNTRL(:) = 0._JPRB
ZCNTRL(1) = HMIN
ZCNTRL(2) = PTSTEP
ZCNTRL(3) = HSTART

! Compute here 'roundoff' number - there is a paralellization
! issue when calling WLMACH as part of kpp-code.
IF ( KSTEP == 0_JPIM ) THEN
  CALL CIFS_KPP_WLAMCH(ROUNDOFF_STORE , 'E')
ENDIF

JAER_TRACER(1)=ISO4
JAER_TRACER(2)=IMSA
JAER_TRACER(3)=INO3_A

JRATE_TRACER(1)=IO3
JRATE_TRACER(2)=IHO2


IF (LAERSOA .AND. LAERSOA_COUPLED) THEN

  JSOG_TRACER(1)=ISOG1
  JSOG_TRACER(2)=ISOG2A
  JSOG_TRACER(3)=ISOG2B

  IF (TRIM(AERO_SCHEME)=="aer" )THEN
    ! Specify aerosol types: SOA - double-check indices..
    ISSOA1 =NTYPAER(1)+NTYPAER(2)+NTYPAER(3)+NTYPAER(4)+NTYPAER(5)+NTYPAER(6)+NTYPAER(7)+1
    ISSOA2 =NTYPAER(1)+NTYPAER(2)+NTYPAER(3)+NTYPAER(4)+NTYPAER(5)+NTYPAER(6)+NTYPAER(7)+2

    JSOA_TRACER(1)=ISSOA1
    JSOA_TRACER(2)=ISSOA2
  ELSE
    WRITE(NULERR,*)'No valid Organic Aerosol interaction available:',LAERSOA,TRIM(AERO_SCHEME)     
  ENDIF
ENDIF


POUT(KIDIA:KFDIA,:,:)  = 0.0_JPRB
PSOGTOSOA(KIDIA:KFDIA,:,:) = 0.0_JPRB

! Lat / Lon
DO JL=KIDIA,KFDIA
  ZLAT(JL)=(180.0_JPRB/RPI)*PGELAT(JL)
  ZLON(JL)=(180.0_JPRB/RPI)*PGELAM(JL)
ENDDO


!  1.1  Compute full level O3 TC

ZRGI=1.0_JPRB/RG

DO JL=KIDIA,KFDIA
  ZCOLO3(JL,0)=0.0_JPRB
ENDDO

DO JLEV=1,KLEV
  DO JL=KIDIA,KFDIA
    ZDELP=PDELP(JL,JLEV)
    ZCOLO3(JL,JLEV)=ZCOLO3(JL,JLEV-1)+MAX(0.0_JPRB,PCEN(JL,JLEV,IO3))*ZDELP*ZRGI
    ZPCO3(JL,JLEV)=ZCOLO3(JL,JLEV)
  ENDDO
ENDDO

! 1.2  Calculate integrated cloud cover above level - adapted from cldpp.f90

ZCC=0.0_JPRB

DO JL=KIDIA,KFDIA
  ZCC(JL,1)=1.0_JPRB-MIN(MAX(PAP(JL,1),REPCLC),1.0_JPRB-REPCLC)
ENDDO

DO JLEV=2,KLEV
  DO JL=KIDIA,KFDIA
    IF(PAP(JL,JLEV-1) < 1.0_JPRB-REPSEC) THEN
      ZCC(JL,JLEV)=ZCC(JL,JLEV-1)*&
       & (1.0_JPRB-MAX(PAP(JL,JLEV),PAP(JL,JLEV-1)))&
       & /(1.0_JPRB-MIN(PAP(JL,JLEV-1),1.0_JPRB-REPSEC))
    ELSE
      ZCC(JL,JLEV) = 0.0_JPRB
    ENDIF
  ENDDO
ENDDO

DO JLEV=1,KLEV
  DO JL=KIDIA,KFDIA
    ZCC(JL,JLEV)=MAX(0.0_JPRB,1.0_JPRB-ZCC(JL,JLEV))
  ENDDO
ENDDO


! Compute Qsat
IFLAG=2
CALL SATUR (KIDIA , KFDIA , KLON  , 1 , KLEV , LPHYLIN,&
  & PRSF1, PTP    , ZQSAT , IFLAG)

! Relative humidity
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZRHCL(JL,JK)=PQP(JL,JK)/(MAX(1.E-30_JPRB, ZQSAT(JL,JK)))
  ENDDO
ENDDO

! 1.3 Initialize tendencies and budget accumulation...
PTENC1(KIDIA:KFDIA,:,:)   = 0.0_JPRB
IF (LCHEM_DIAC) THEN
  PBUDJ(KIDIA:KFDIA,:,:) = 0.0_JPRB 
  PBUDR(KIDIA:KFDIA,:,:) = 0.0_JPRB 
  PBUDX(KIDIA:KFDIA,:,:) = 0.0_JPRB
ENDIF

! 1.4 calculate cloud optical depth 
! please note that the arguments to cod_op range from 0:KLEV 

LLCOD_TM5=.TRUE.
IF ( LLCOD_TM5 ) THEN
! * IFS scheme - modified (optimized) for TM5 format...
  CALL COD_OP_TM5(YDDIMV,YDMODEL%YRML_PHY_RAD%YRERAD,KIDIA,KFDIA, KLON,KLEV, 1,PQP, PTP, PAP, PRS1, PRSF1 , PLSM, PWND, PLP, PIP,&
     &     ZCLOUD_REFF, ZTCTAUC,ZTAUS_CLD, ZTAUA_CLD, ZPMCLD)
ELSE
!* Use New Photolysis scheme 
!* new cloud optical depth routine
  CALL TM5_SLINGO(KIDIA,KFDIA,KLON,KLEV,PIP,PLP,PAP,PGEOH,PRS1,PTP,&
     & ZTAUA_CLD,ZTAUS_CLD,ZPMCLD,ZCLOUD_REFF)
ENDIF

! * calculate the aerosol scattering/absorption 

! ITAU_MACC = 0_JPIM ! Allways switch off aerosol optical depth 
! ITAU_MACC = 1_JPIM  ! redundant , use MACC aerosol if LCHEM_AEROI  
! ITAU_MACC = 2_JPIM  ! Use simple (but wrong) climatology for aerosol optical depth if no prognostic MACC aerosol is used 

IF ( LCHEM_AEROI ) THEN 
  IF (TRIM(AERO_SCHEME)=="aer" )THEN
    ! * from MACC fields
    CALL TM5_MACC_AEROSOL(KIDIA,KFDIA,KLON,KLEV, KAERO, &
      &   PRS1   , PAEROP  , ZRHCL   ,   &
      &  ZTAUS_AER,ZTAUA_AER,ZPMAER)
  ELSEIF (TRIM(AERO_SCHEME)=="glomap") THEN
    CALL TM5_GLOMAP_AEROSOL(KIDIA,KFDIA,KLON,KLEV,   &
    &  PAERAOT, PAERAAOT, PAERASY, &
    &  ZTAUS_AER,ZTAUA_AER,ZPMAER)
  ELSE
    WRITE(NULERR,*)'No valid Aerosol interaction available:',LCHEM_AEROI,TRIM(AERO_SCHEME)     
  ENDIF
ELSE 
  IF ( ITAU_MACC==0_JPIM ) THEN
! * Zero aerosol fields...
     ZTAUS_AER=0._JPRB
     ZTAUA_AER=0._JPRB
     ZPMAER=0._JPRB
  ELSEIF ( ITAU_MACC==2_JPIM ) THEN
! * Interpolate/ extrapolate aerosol info
    CALL TM5_AEROSOL_INFO(KIDIA,KFDIA,KLON,KLEV,&
     &   PGEOH,PTP, PLSM,PRS1 ,PRSF1, PAP, PQP, PGELAT,&
     &   ZTAUS_AER,ZTAUA_AER, ZPMAER  )
  ELSE
    WRITE(NULERR,*)'This option for AOD is not available: ITAU_MACC, Aerosol interaction ',ITAU_MACC, LCHEM_AEROI     
  ENDIF
ENDIF 

!
! 1.5 Calculate photolysis rates 
! CHECK negative cloud taus

  IF ( ANY ( ZTAUS_CLD(KIDIA:KFDIA,1:KLEV) < 0.0_JPRB)) THEN 
!      Print*, ' negative cloud taus', MINVAL(ZTAUS_CLD(KIDIA:KFDIA,1:KLEV)), MINLOC(ZTAUS_CLD(KIDIA:KFDIA,1:KLEV))
      WHERE ( ZTAUS_CLD(KIDIA:KFDIA,1:KLEV)  < 0.0_JPRB ) ZTAUS_CLD(KIDIA:KFDIA,1:KLEV) = 0.0_JPRB     
  ENDIF 

!
CALL TM5_PHOTO_FLUX(KIDIA, KFDIA, KLON, KLEV, PTP, PCSZA, PALB,& 
   &          ZPCO3, PRSF1, PRS1,&
   &          ZTAUA_CLD,ZTAUS_CLD,ZPMCLD,&
   &          ZTAUA_AER,ZTAUS_AER,ZPMAER,&
   &          PAP,PGEOH, ZRJ )


! 1.6 adjust photolysis rate to sund earth distance variation
IMONTH=NMM(NINDAT)
IDAY=NDD(NINDAT)
CALL TM5_SUNDIS(KIDIA,KFDIA,KLON,KLEV,IMONTH,IDAY, ZRJ)

!
! 1.7  for PZRJ output: Output photolysis rates
IF (LCHEM_JOUT) THEN
  DO JLEV=1,KLEV
    DO JL=KIDIA,KFDIA
     POUT(JL,JLEV,2) = ZRJ(JL,JLEV,1)
     POUT(JL,JLEV,3) = ZRJ(JL,JLEV,2)
     POUT(JL,JLEV,4) = ZTAUS_AER(JL,JLEV,1,1)
     POUT(JL,JLEV,5) = ZTAUA_AER(JL,JLEV,1,1)
    ENDDO
  ENDDO
!  DO JLEV=1,KLEV
!    DO JL=KIDIA,KFDIA
!     POUT(JL,JLEV,2) = ZRJ(JL,JLEV,JNO2)
!    ENDDO
!  ENDDO
ENDIF


! resrtict active chemistry to below 1 hPa 
ILEVMIN = 1
DO JK=1,KLEV
  IF ( YDVAB%VAH(JK) + YDVAB%VBH(JK) * 101300.0_JPRB <= 100_JPRB ) ILEVMIN =  JK
ENDDO
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF ( YDVAB%VAH(JK) + YDVAB%VBH(JK) * 101300.0_JPRB <= PPTROPO(JL) ) ILEVMIN_GAS(JL) = JK 
  ENDDO
ENDDO


ZSAD_AER(KIDIA:KFDIA,:)=0.0_JPRB
ZSAD_CLD(KIDIA:KFDIA,:)=0.0_JPRB
ZSAD_ICE(KIDIA:KFDIA,:)=0.0_JPRB

!2.0 loop over levels for solving chemistry
DO JK=ILEVMIN,KLEV

!2.2 calculate reaction rates ZRR using pre-calculated, temperature-dependent values 
!    also finalize evaluation of ZRJ
   DO JT=1,NPHOTO
     ZRJ_IN(KIDIA:KFDIA,JT)=ZRJ(KIDIA:KFDIA,JK,JT)
   ENDDO

    ZRR(KIDIA:KFDIA,:) = 0.0_JPRB ! 97-100 are not assigned TM5_CALRATES

   CALL TM5_CALRATES(YDCHEM,YDCOMPO,YDEAERSNK,YDEAERSRC,YDEAERATM,YGFL, &
     & KIDIA, KFDIA, KLON, JK, KLEV, JRATE_TRACER, KAERO, NMODES, PTP, PRSF1 ,&
     & PQP,PAP, PIP,PLP, ZRJ_IN, ZCLOUD_REFF, PCEN(:,JK,IO3), PCEN(:,JK,IHO2),&
     & PCEN(:,JK,INH4),PCEN(:,JK,INO3_A), PCEN(:,JK,ISO4),PWETDIAM,&
     & PWETVOL,PND, &
     & PAEROP,ZSAD_AER,ZSAD_CLD, ZSAD_ICE,&
     & ZRR)     

!2.3 convert concentrations (mass mixing ratios) to #/cm3
   DO JL=KIDIA,KFDIA
!*   ZAIRD(JL) = 7.24291e24_JPR * PRSF1(JL,JK)/PTP(JL,JK) * 1e-3_JPRB * 1.0e-6_JPRB * 10.0_JPRB
!*   multiply with RMD (dry air molar mass) for efficiency  
     ZAIRDM(JL) = (7.24291E16_JPRB * PRSF1(JL,JK) / PTP(JL,JK)) * RMD
!*   fill some extra fields, required for budget calculations:
     ZCVM(JL,IAIR) = ZAIRDM(JL) / RMD
     ZCVM(JL,IH2O) = PQP(JL,JK)*ZAIRDM(JL)/XMH2O
     ZCVM(JL,IACID) = 0._JPRB

     DO JT=1,NCHEM
!*       assure positivity for initial concentrations
      ZCVM0(JL,JT) = MAX(PCEN(JL,JK,JT) / YCHEM(JT)%RMOLMASS * ZAIRDM(JL), 0._JPRB) 
     ENDDO

     IF (LAERCHEM) THEN
!*       use aerosol scheme SO4 for WETCHEM and EQSAM
      ISSO4 = SUM(NTYPAER(1:4)) + 1
      ZCVM0(JL,ISO4) = MAX(PAEROP(JL,JK,ISSO4) / YCHEM(ISO4)%RMOLMASS *ZAIRDM(JL), 0._JPRB)        
     ENDIF

     IF (LAERNITRATE) THEN
!*       use aerosol scheme NO3_A for WETCHEM and hetchem. Add both coarse and fine mode into one
      ISNO3 = SUM(NTYPAER(1:5)) + 1
      ZCVM0(JL,INO3_A) =                    MAX(PAEROP(JL,JK,ISNO3)   / YCHEM(INO3_A)%RMOLMASS * ZAIRDM(JL), 0._JPRB)
      ZCVM0(JL,INO3_A) = ZCVM0(JL,INO3_A) + MAX(PAEROP(JL,JK,ISNO3+1) / YCHEM(INO3_A)%RMOLMASS * ZAIRDM(JL), 0._JPRB)
!*       use aerosol scheme NH4 for WETCHEM and hetchem.
      ISNH4 = SUM(NTYPAER(1:6)) + 1
      ZCVM0(JL,INH4) = MAX(PAEROP(JL,JK,ISNH4) / YCHEM(INH4)%RMOLMASS * ZAIRDM(JL), 0._JPRB)
     END IF

!*       PSC is not part of chemistry 
     ZCVM0(JL,IPSC) = 0.0_JPRB  

!*     initialize final (ZCVM) and intermediate (ZCVM1) concentrations
     ZCVM(JL,1:NCHEM)  = ZCVM0(JL,1:NCHEM)   
     ZCVM1(JL,1:NCHEM) = ZCVM0(JL,1:NCHEM)
   ENDDO

!2.4 Treatment of dry deposition within chemistry:
   ZVD(KIDIA:KFDIA,:) = 0.0_JPRB
! see new routine CHEM_DRYDEP 
!   IF ( .NOT. LCHEM_DDFLX )THEN
!   !   deposition velocities are applied in chemistry!!'
!      IF (JK == KLEV) then 
!         DO JL=KIDIA,KFDIA
!           ZHGT(JL) =  PGEOH(JL,KLEV-1)  * ZRGI  ! height in meter
!         ENDDO
!         IDRYDEP=0 
!         DO  JT=1,NCHEM
!           IF (YCHEM(JT)%IGRIBDV > 0 ) THEN
!             IDRYDEP=IDRYDEP+1
!             DO JL=KIDIA,KFDIA
!               ZVD(JL,JT) =  PDV(JL,IDRYDEP) / ZHGT(JL)  
!             ENDDO 
!           ENDIF 
!         ENDDO 
!      ENDIF 
!   ENDIF 

!2.5 set ZCR2 / ZCR3 budget accumulators for this level to zero.
   IF (LCHEM_DIAC) THEN
      ZCR2(KIDIA:KFDIA,:) = 0._JPRB
      ZCR3(KIDIA:KFDIA,:) = 0._JPRB
   ENDIF

!2.6 Initialize NOy components - should not be used!!!
   ! CALL TM5_NOy_mass_init(KIDIA,KFDIA,KLON,ZCVM1)

!3.0 Select chemistry solver 
IF ( KCHEM_SOLVE == 1_JPIM) THEN
!  3.0 Perform iterations for chemical solver
   ITEREBI=MAX(1_JPIM,NINT(PTSTEP/1350_JPIM))  !needed if PTSTEP > 2400
   ZDT_EBI=PTSTEP/ITEREBI
   DO IB=1_JPIM,ITEREBI  


!3.1  wet sulphur/ammonia chemistry
      CALL TM5_WETCHEM(YGFL,KIDIA,KFDIA,KLON,ZDT_EBI,PTP(:,JK),PAP(:,JK),PRSF1(:,JK),PLP(:,JK),ZCVM1,ZHPLUS,ZCVM)

      IMAXIT=8_JPIM ! maximum EBI iterations
      IF(KLEV+1-JK<=4) IMAXIT = IMAXIT*2   ! surface layers more iterations

!3.2  Call 'euler backward integration' - solver for individual levels
      IF (LCHEM_REVCHEM) THEN
        ! Troposphierc Chemistry version 'tc02b':
        ! Updated isoprene chemistry - depreciated. No longer maintained
        !VH CALL TM5_DO_EBI_TC02B(YGFL,KIDIA,KFDIA,KLON,IMAXIT,ZDT_EBI,ZRR,ZRJ_IN,ZCVM1,ZCVM,ZVD,PRSF1(:,JK))
      ELSE
        ! Troposphierc Chemistry version 'tc02a':
        ! Standard chemistry
        CALL TM5_DO_EBI(YGFL,KIDIA,KFDIA,KLON,IMAXIT,ZDT_EBI,ZRR,ZRJ_IN,ZCVM1,ZCVM,ZVD,PRSF1(:,JK))
      ENDIF

!      ZCVM(KIDIA:KFDIA,1:NREAC) =ZCVM1(KIDIA:KFDIA,1:NREAC)

!3.3   Perform NOy mass correction step - depreciated
       ! CALL TM5_NOy_MASS_END(KIDIA,KFDIA,KLON,ZCVM1,ZCVM,ZVD,ZDT_EBI)


      IF (LCHEM_DIAC) THEN
!3.4     increase budget accumulators ZCR2 and ZCR3 (photolysis / gas-phase chem)
         CALL TM5_IBUD(YGFL,KIDIA,KFDIA,KLON,ZRR,ZRJ_IN,ZCVM,ZCR2,ZCR3)
      ENDIF

!3.5 Update initial concentrations in solver
      DO JL=KIDIA,KFDIA
        DO JT=1,NCHEM 
          ZCVM1(JL,JT) = ZCVM(JL,JT)
        ENDDO
      ENDDO

   ENDDO ! iterebi
ELSEIF( KCHEM_SOLVE == 2_JPIM) THEN
   ! Required for budget accumulation 'RBUD'
   ITEREBI = 1_JPIM
   ! select KPP...
!3.1  wet sulphur/ammonia chemistry
      CALL TM5_WETCHEM(YGFL, KIDIA, KFDIA, KLON,PTSTEP,PTP(:,JK),PAP(:,JK),PRSF1(:,JK),PLP(:,JK),ZCVM1,ZHPLUS,ZCVM)

!3.2  Call selected kpp solver 
DO JL=KIDIA,KFDIA

          ! starting value for integration time step
          ZCNTRL_P=ZCNTRL
  !VH --- to be switched on...    ZCNTRL_P(3) = ZHSAVE_KPP(JL,JLEV)

          ! Try to initialize timestep with something reasonable...
          !ZCNTRL_P(3) = ZDT_CHEM/20._JPRB

          ! Near surface reduce initial time step near surface (at least lowest 2 levels)
          IF ( PRSF1(JL,JK)> 80000_JPRB .OR. JK > KLEV-6 ) THEN
            ZCNTRL_P(3)=ZCNTRL_P(3)/4.
          ENDIF

          ! Special fix for HO2+HO2 -> H2O2 reactions.
          ! In tm5_calrates.F90 and in tm5_do_ebi.F90 these are treated as single-body reactions,
          ! but in KPP these are treated as second-order reaction.
          ! Now scale kho2_aer and kho2_liq with HO2 concentrations, 
          ! to account for it as a self-reaction in KPP solver
          ! To ensure a reasonably small reaction rate assume a minimum HO2 concentration 
          ZCONC_HO2 = MAX(ZCVM(JL,IHO2), 1E5_JPRB)
          ZRR(JL,KHO2L)   =ZRR(JL,KHO2L)    / ZCONC_HO2
          ZRR(JL,KHO2_AER)=ZRR(JL,KHO2_AER) / ZCONC_HO2
  
          ! Update kpp rates... (Assume dry deposition is done outside chemistry!!!
          CALL TM5_KPP_RATES(ZRR(JL,1:NREAC),ZRJ_IN(JL,1:NPHOTO),ZRCONST)
               
          ! Initialize concentrations to KPP (merge with prepare_kpp_conc)
          CALL TM5_V0_KPP_INITIALIZE(YGFL,ZCVM(JL,1:NCHEM),ZVAR(1:NVAR),LAERSOA_COUPLED)

          ! Call kpp integrator...
          CALL TM5_KPP_INTEGRATOR(0._JPRD, PTSTEP, ICNTRL,ZCNTRL_P, ISTATUS,ZSTATE,IERR,ZVAR,ZRCONST,ZCVM(JL,iair) )

          ! update new preferred time step...
          ! ZHSAVE_KPP(JL,JLEV) = ZSTATE(Nhnew)
          
          !- Filter error due to bad concentrations
         IF (IERR>0) THEN

           !No errors: update concentrations...
           CALL TM5_V0_KPP_UPDATE_CIFS_CONC(YGFL,ZVAR(1:NVAR),ZCVM(JL,1:NCHEM),LAERSOA_COUPLED)
     
          ELSEIF( IERR==-9) THEN
            write(NULERR,'(a)') '     ZVAR below are out of range:'
            DO JT= 1, NVAR
               IF( ZVAR(JT)/ZCVM(JL,iair) > VMR_BAD_LARGE) THEN
                 WRITE(NULERR,'(a,2(i5),a,es12.5)') '  vmr-idx ',JT,JK, ' ; reached value: ', ZVAR(JT)/ZCVM(JL,iair)
               ENDIF
            ENDDO
            WRITE(NULERR,*) '  -> chem integrator skipped'
          ENDIF

          ! Check on concentrations
          DO JT=1,NCHEM 
            ZCVM(JL,JT) = MAX(ZCVM(JL,JT),0.0_JPRB)
          ENDDO

        ENDDO

        IF (LCHEM_DIAC) THEN
!3.4           increase budget accumulators ZCR2 and ZCR3 (photolysis / gas-phase chem)
           CALL TM5_IBUD(YGFL, KIDIA, KFDIA, KLON, ZRR, ZRJ_IN, ZCVM, ZCR2, ZCR3 )
        ENDIF

!3.5 Update initial concentrations in solver
        DO JL=KIDIA,KFDIA
          DO JT=1,NCHEM 
            ZCVM1(JL,JT) = ZCVM(JL,JT)
          ENDDO
        ENDDO

 ELSE
   WRITE(NULOUT,*) ' WARNING: No valid option for solver selected!'  
 ENDIF
 
!3.6 Check on final concentrations  - should not be necessary
   !  DO JL=KIDIA,KFDIA
   !    DO JT=1,NCHEM 
   !      ZCVM(JL,JT) = max(ZCVM(JL,JT),0.0_JPRB)
   !    ENDDO
   !  ENDDO

!3.7 add budgets for this timestep
   IF (LCHEM_DIAC) THEN
      CALL TM5_RBUD(YGFL,KIDIA,KFDIA,KLON,KLEV,JK,IOH,ITEREBI,ZAIRDM,ZCR2,ZCR3,PBUDJ,PBUDR,PBUDX)
   ENDIF

!4.0 Eqsam solver for aerosol- gas phase interaction 
! done in aerosols in LAERNITRATE
   IF (.NOT. LAERNITRATE) THEN
     DO JL=KIDIA,KFDIA
       ! Calculate RH - relative humidity
!       ZTR = 1._JPRB - 373.15/PTP(JL,JK)
!       ZWV=EXP((((-.1299_JPRB*ZTR-.6445)*ZTR-1.976_JPRB)*ZTR+13.3185_JPRB)*ZTR)
!       ZRRH = ZCVM(JL,IH2O)*PTP(JL,JK)/(1013.25*ZWV*7.24E16)
!       ZRH = 0.01_JPRB*MAX(0.01_JPRB, MIN(ZRRH, 99.9_JPRB ) )   ! 0-0.999 scale!
       ZRH = MAX(0.01_JPRB,MIN(ZRHCL(JL,JK),0.99_JPRB))
       ! scale relative humidity to cloudfree part
       ! assuming 100% rh in the cloudy part, but never smaller than 0.75!
       IF (ZRH > 0.90) THEN 
         ZCCS = PAP(JL,JK)
         IF((1._JPRB - ZCCS) > TINY(ZCCS)) ZRH = MAX(0.75_JPRB, (ZRH-ZCCS)/(1._JPRB - ZCCS)) 
       ENDIF
      
       ZNH = (ZCVM(JL,INH3) + ZCVM(JL,INH4))/6.02E23_JPRB * 1E6_JPRB ! molec/cm3 -> mol/m3
       ZNO3 = (ZCVM(JL,INO3_A) + ZCVM(JL,IHNO3))/6.02E23_JPRB * 1E6_JPRB ! molec/cm3 -> mol/m3
       ZSO4 = ZCVM(JL,ISO4)/6.02E23_JPRB*1E6_JPRB ! molec/cm3 -> mol/m3 . SO4 is ony used in input
       CALL TM5_EQSAM(PTP(JL,JK),ZRH,ZNH,ZNO3,ZSO4,ZYEQ)
       ZCVM(JL,IHNO3) = ZYEQ(1)*6.02E23*1E-6 !  mol/m3 -> molec/cm3 
       ZCVM(JL,INH3)  = ZYEQ(2)*6.02E23*1E-6
       ZCVM(JL,INH4)  = ZYEQ(3)*6.02E23*1E-6
       ZCVM(JL,INO3_A)= ZYEQ(4)*6.02E23*1E-6
     ENDDO
   ENDIF

   IF (LAERSOA .AND. LAERSOA_COUPLED) THEN
     ! apply gas-aerosol partitioning for Sec Org Aerosol, and gas-phase version
     INBDU=NTYPAER(1)+NTYPAER(2)
     INBOM=NTYPAER(1)+NTYPAER(2)+NTYPAER(3)
     DO JL=KIDIA,KFDIA
       IF (JK > ILEVMIN_GAS(JL)) THEN
         ZRHO=PRSF1(JL,JK)/(RD*PTP(JL,JK))
  
         ! Call to simple SOG-solver, based on Euler Backward Integrator
         ! No longer required, as SOG production is now integral part of KPP-solver..
         !CALL TM5_SOG(NCHEM,1_JPIM, PTSTEP,ZRR(JL,1:NREAC),ZCVM1(JL,1:NCHEM),ZCVM(JL,1:NCHEM))
  
         DO JT = 1,INSOG
           ! SOG: conversion from molec/cm3 -> mol/cm3-> gr/cm3-> ug/m3
           ZSOG(JT)=ZCVM(JL,JSOG_TRACER(JT)) /6.02E23_JPRB * 1E12_JPRB * YCHEM(JSOG_TRACER(JT))%RMOLMASS   ! ug/m3
         ENDDO
         DO JT = 1,INSOA
           ! SOA: conversion from kg/kg to ug/m3
           ZSOA(JT)=PAEROP(JL,JK,JSOA_TRACER(JT))*ZRHO* 1E9_JPRB  
         ENDDO
         ! Input SOG for SOA/SOG equilibrium: Only first (biogenic) and third (anthro, low volatility) arrays
         ZSOGH(1)=ZSOG(1)
         ZSOGH(2)=ZSOG(3)

         ! compute total pre-existing aerosol , so far only organic aerosol
         ZORGAERO=0._JPRB
         DO JT=INBDU+1,INBOM
           ZORGAERO=ZORGAERO+PAEROP(JL,JK,JT)*ZRHO*1E9_JPRB ! kg/kg -> ug/m3
         ENDDO
  
         ! Call to simple SOA-solver for two SOG/SOA tracers
         CALL TM5_SOA(INSOA, PTP(JL,JK),ZSOGH,ZSOA,ZORGAERO)
   
         ! update ZSOG: Only first (biogenic) and third (anthro, low volatility) arrays
         ZSOG(1)=ZSOGH(1)
         ! ZSOG(2) is not changed: All assumed to remain in gas-phase.
         ZSOG(3)=ZSOGH(2)

         DO JT=1,INSOG
           !conversion back to molec/cm3 
           ZCVM(JL,JSOG_TRACER(JT)) =ZSOG(JT)*6.02E23_JPRB / 1E12_JPRB / YCHEM(JSOG_TRACER(JT))%RMOLMASS   
         ENDDO
         DO JT=1,INSOA  
           !conversion back to kg/kg 
           ZSOA(JT)=ZSOA(JT)/(ZRHO*1E9_JPRB) 

           ZSOA_TMP=ZSOA(JT)
           ! Apply photolytic loss, Hodzic et al. (2016)
           ZJSOA=4E-4_JPRB*ZRJ_IN(JL,JNO2)
           ZXLSOA=1._JPRB/(1._JPRB+ZJSOA*PTSTEP)
           ZSOA(JT) = ZSOA(JT) *ZXLSOA
            ! Budget diagnostics for photolysis loss...
           ! For testing purposes only process photolysis for second tracer
           IF (LCHEM_DIAC ) THEN
             ! Put budget in additional fields JT, 1-2.. (SOA loss is positive)
             PBUDJ(JL,JK,NPHOTO+JT)=PBUDJ(JL,JK,NPHOTO+JT)- (ZSOA(JT)-ZSOA_TMP)/PTSTEP  !units kg/kg/s  
           ENDIF

           ! Compute aerosol tendencies (kg/kg/sec), to be applied in aer_phy3
           PSOGTOSOA(JL,JK,JT) =  (ZSOA(JT)-PAEROP(JL,JK,JSOA_TRACER(JT)))/PTSTEP
         ENDDO
       ELSE ! JK > ILEVMIN_GAS(JL)

        !Five-day lifetime in stratosphere, (representing photolysis) to prevent build-up..
        ZTSCALI =  1./(1._JPRB+2.315E-6_JPRB*PTSTEP)
        DO JT=1,INSOG
           ZCVM(JL,JSOG_TRACER(JT)) =  ZCVM(JL,JSOG_TRACER(JT)) *ZTSCALI
        ENDDO
        DO JT=1,INSOA
           !lifetime in stratosphere to prevent build-up..
           ZSOA(JT)=PAEROP(JL,JK,JSOA_TRACER(JT))*ZTSCALI
           ! Compute aerosol tendencies (kg/kg/sec), to be applied in aer_phy3
           PSOGTOSOA(JL,JK,JT) =  (ZSOA(JT)-PAEROP(JL,JK,JSOA_TRACER(JT)))/PTSTEP
        ENDDO


       ENDIF ! JK > ILEVMIN_GAS(JL)
     ENDDO
   ENDIF ! LAERSOA

   DO JL=KIDIA,KFDIA
     ! Fill tropospheric CO2 production tendency computed from CO + OH reaction budget
     ! units: kg CO2 / kg / sec
     ZFAC=RMCO2/(ZAIRDM(JL)*PTSTEP)
     PCHEM2GHG(JL,JK,2) = ZCVM(JL,ICO)*ZCVM(JL,IOH)*ZRR(JL,KCOOH)  * ZFAC
   ENDDO


!4.9: O3S handling

  DO JL=KIDIA,KFDIA
    IF ( JK >= KLEVTROP(JL) ) THEN
      ! Compute loss of O3S in troposphere at location JL
      CALL TM5_O3S(YDCHEM,NCHEM,PTSTEP,ZRR(JL,1:NREAC),ZRJ_IN(JL,1:NPHOTO),ZCVM(JL,1:NCHEM) )
    ENDIF
  ENDDO


!5.0 convert concentration tendencies to mass mixing ratio
   DO JL=KIDIA,KFDIA
     DO JT=1,NCHEM 
       IF (JT /= IPSC) THEN
         PTENC1(JL,JK,JT) =   (ZCVM(JL,JT)-ZCVM0(JL,JT)) * YCHEM(JT)%RMOLMASS /(ZAIRDM(JL)*PTSTEP)
       ELSE
         PTENC1(JL,JK,JT) = 0.0_JPRB        
       ENDIF  
     ENDDO
   ENDDO

! check for large values 
!   IF (LLCHECKMAX) THEN
!     DO JL=KIDIA,KFDIA
!       IF ( PTENC1(JL,JK,IHNO3) * PTSTEP > 1.0e-4 ) THEN 
!          print*, ' HNO3 large after ebi ', ZCVM(JL,IHNO3) , PTENC1(JL,JK,IHNO3) * PTSTEP, ZCVM0(JL,IHNO3) , JL, JK 
!        ENDIF
!        IF ( PTENC1(JL,JK,INO2) * PTSTEP > 1.0e-4 ) THEN 
!           print*, ' NO2 large after ebi ', ZCVM(JL,INO2) , PTENC1(JL,JK,INO2) * PTSTEP, ZCVM0(JL,INO2) , JL, JK 
!        ENDIF  
!     ENDDO
!   ENDIF   

ENDDO ! loop over levels

!6.0  cariolle ozone chemistry  

ZTENO3COR(KIDIA:KFDIA,1:KLEV) = 0.0_JPRB

! for O3S - No longer activated
! CALL O3CHEM (YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_CHEM%YROZO,YDMODEL%YRML_PHY_MF%YRPHY2, &
!    & KIDIA, KFDIA, KLON, 1, KLEV, KVCLIS, PGEMU, PCSZA,&
!    & PRS1, PRSF1, PKOZO, PDELP, PTP, PCEN(:,:,IO3S),&
!    & ZTENO3SCOR)

! for O3
 CALL O3CHEM (YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_CHEM%YROZO,YDMODEL%YRML_PHY_MF%YRPHY2, &
    & KIDIA, KFDIA, KLON, 1, KLEV, KVCLIS, PGEMU, PCSZA,&
    & PRS1, PRSF1, PKOZO, PDELP, PTP, PCEN(:,:,IO3),&
    & ZTENO3COR)

!ENDIF

! 6.1 replace tendecies for ozone above ZPMAXO3CAR with Cariolle tendency.
! The stratosphere is defined at zonal mean ozone levels exceeding 150 ppb
! based on a climatology. This boundary can approximately be described
! by the function P = 230-148*( cos(lat) )^4 (hPa)
! stratospheric Cariolle chemistry at (P-20) hPa is applied.

 DO JK=1,KLEV
   DO JL=KIDIA,KFDIA
     ZPMAXO3CAR=PPTROPO(JL) 
     ZPMINO3CAR=ZPMAXO3CAR-4000_JPRB
     IF (PRSF1(JL,JK)<= ZPMAXO3CAR)  THEN
         ! At least partialy stratospheric tendencies
      IF (PRSF1(JL,JK)> ZPMINO3CAR)  THEN
         ! Transition region ZPMAX > P > ZPMIN
         ZPO3CAR_RATIO=(PRSF1(JL,JK)-ZPMINO3CAR)/(ZPMAXO3CAR-ZPMINO3CAR)
         PTENC1(JL,JK,IO3) =  ZPO3CAR_RATIO* PTENC1(JL,JK,IO3) +(1_JPRB-ZPO3CAR_RATIO)* ZTENO3COR(JL,JK) 
      ELSE
        ! completely stratospheric region: P < ZPMIN
        PTENC1(JL,JK,IO3) =  ZTENO3COR(JL,JK) 
     ENDIF

    ENDIF
   ENDDO
 ENDDO

! Special treatment for O3S in stratosphere
DO JL=KIDIA,KFDIA
  DO JK=1,KLEV
    IF (JK < KLEVTROP(JL)) THEN
      ! Stratosphere: compute tendency to arrive at same concentration as O3 itself
      ZCEN_O3=PCEN(JL,JK,IO3)+PTENC1(JL,JK,IO3)*PTSTEP
      PTENC1(JL,JK,IO3S)=(ZCEN_O3-PCEN(JL,JK,IO3S))/PTSTEP
    ENDIF
  ENDDO
ENDDO



!
! 6.2 replace CH4 at surface and in stratopshere with value from boundary condition. 
IF (.NOT. LCHEM_ANACH4) THEN
  IYEAR=NCCAA(NINDAT)
  CALL TM5_BOUNDARY_CH4(YGFL,KIDIA,KFDIA,KLON,IMONTH,IYEAR,ICH4,PGELAT,PCEN(:,KLEV,ICH4),ZTENBC(1:KLON,1))

  IF (LCHEM_WEAK_CH4_RELAXATION) THEN
    ZTAU_NUDGE=86400._JPRB
  ELSE
    ZTAU_NUDGE=2500._JPRB
  ENDIF
  
  DO JL=KIDIA,KFDIA
    ! Introduce more relaxed constraint on CH4 for background conditions 
    ! This should be better when including CH4 emissions and dry deposition.  

    PTENC1(JL,KLEV,ICH4) =  ((1._JPRB - EXP(-PTSTEP / (ZTAU_NUDGE)))*ZTENBC(JL,1))  / PTSTEP
    ! store this special LBC CH4 budget contribution  POUT(:,1,1) 
    POUT(JL,1,1)=PTSTEP *PTENC1(JL,KLEV,ICH4)*PDELP(JL,KLEV) / RG 

    ! Original, stringent nudging:
    !PTENC1(JL,KLEV,ICH4) =  ZTENBC(JL,1)  / PTSTEP
    ! store this special LBC CH4 budget contribution  POUT(:,1,1) 
    !POUT(JL,1,1)=PTSTEP *PTENC1(JL,KLEV,ICH4)*PDELP(JL,KLEV) / RG 
  ENDDO

! Apply stratospheric boundary conditions for CH4 at 2 altitude levels: 90 & 45 hPa.

  ! initialize JLEVBC to prevent any non-initialized values..
  JLEVBC(KIDIA:KFDIA,1:4)=1

  ILEV_CH4 =1
  !* estimate level at which pressure is 45 hPa
  DO JL=KIDIA, KFDIA   
    DO JK=1,KLEV
      IF (PRSF1(JL,JK) <= 4500_JPRB) THEN 
        JLEVBC(JL,1) = JK
      ENDIF
    ENDDO
    ZCENBC1(JL,1) =  PCEN(JL,JLEVBC(JL,1),ICH4)
  ENDDO
 
  CALL TM5_STRATBC_CH4(YGFL,KIDIA,KFDIA,KLON,PTSTEP,IMONTH,PGELAT,ILEV_CH4,ZCENBC1(1:KLON,1),ZTENBC(1:KLON,1))
   
  DO JL=KIDIA,KFDIA
    PTENC1(JL, JLEVBC(JL,1),ICH4) =  ZTENBC(JL,1)
    ! store this special UBC CH4 ibudget contribution into POUT(:,2,1) 
    POUT(JL,2,1) = PTSTEP * PTENC1(JL,JLEVBC(JL,1),ICH4)*PDELP(JL,JLEVBC(JL,1)) / RG 
  ENDDO

  ILEV_CH4 =2
  !* estimate level at which pressure is 90 hPa
  DO JL=KIDIA, KFDIA   
    JLEVBC(JL,1) = 1 ! for safety, take model top if not found
    DO JK=1,KLEV
      IF (PRSF1(JL,JK) <= 9000_JPRB) THEN 
        JLEVBC(JL,1) = JK
      ENDIF
    ENDDO
    ZCENBC1(JL,1) =  PCEN(JL,JLEVBC(JL,1),ICH4)
  ENDDO

  CALL TM5_STRATBC_CH4(YGFL,KIDIA,KFDIA,KLON,PTSTEP,IMONTH,PGELAT,ILEV_CH4,ZCENBC1(1:KLON,1),ZTENBC(1:KLON,1))

! Only overwrite extra-tropical tendencies (90hPa in tropics is stil troposphere)
  DO JL=KIDIA,KFDIA
    IF (ABS(ZLAT(JL)) > 30._JPRB) THEN
      PTENC1(JL, JLEVBC(JL,1),ICH4) =  ZTENBC(JL,1)
      ! store this special UBC CH4 ibudget contribution into POUT(:,2,1) 
      POUT(JL,2,1) = PTSTEP * PTENC1(JL,JLEVBC(JL,1),ICH4)*PDELP(JL,JLEVBC(JL,1)) / RG 
    ENDIF
  ENDDO
ENDIF


! IF (LGHG_CHEMTEND_CH4) THEN
  ! Fill field with CH4 tendency term for analysis purposes. units are in kg/kg/sec. 
  ! Different to version coded in chem_bascoetm5, this tendency term includes
  ! contribution of surface nudging (if any), as computed above.
  ! Also the stratospheric nudging is included. This tendency output from this chem_tm5 module
  ! should be considered with great care, i.e. it is not well-defined.
  ! Change units from kg/kg/sec to sec-1. Loss term is positive.
  DO JLEV=1,KLEV
    DO JL=KIDIA,KFDIA
      PCHEM2GHG(JL,JLEV,1) = - PTENC1(JL,JLEV,ICH4) / PCEN(JL,JLEV,ICH4)
    ENDDO
  ENDDO
! ENDIF

IMONTH=NMM(NINDAT)

!* estimate level at which pressure is 10, 50, 100 and 1 hPa. 

DO JL=KIDIA, KFDIA   
  JLEVBC(JL,:)=1 ! for safety, take model top if not found
  DO JK=1,KLEV
    IF (PRSF1(JL,JK) <= 1000_JPRB) THEN 
      JLEVBC(JL,1) = JK
    ENDIF
    IF (PRSF1(JL,JK) <= 5000_JPRB) THEN 
      JLEVBC(JL,2) = JK
    ENDIF
    IF (PRSF1(JL,JK) <= 10000_JPRB) THEN 
      JLEVBC(JL,3) = JK
    ENDIF
    IF (PRSF1(JL,JK) <= 100_JPRB) THEN 
      JLEVBC(JL,4) = JK
    ENDIF
  ENDDO
  ZCENBC1(JL,1) =  PCEN(JL,JLEVBC(JL,1),IHNO3)
  ZCENBC2(JL,1) =  PCEN(JL,JLEVBC(JL,1),IO3)
  ZCENBC1(JL,2) =  PCEN(JL,JLEVBC(JL,2),IHNO3)
  ZCENBC2(JL,2) =  PCEN(JL,JLEVBC(JL,2),IO3)
  ZCENBC1(JL,3) =  PCEN(JL,JLEVBC(JL,3),IHNO3)
  ZCENBC2(JL,3) =  PCEN(JL,JLEVBC(JL,3),IO3)
ENDDO

!6.3 replace tm5-HNO3 concentrations at 10 hPa with value from boundary condition. 
CALL TM5_BOUNDARY_HNO3(YGFL,IMODE_HNO3BC,KIDIA,KFDIA,KLON,IMONTH,PGELAT,ZCENBC1,ZCENBC2,ZTENBC)

!*  HNO3 constraint at 10 hpa (nothing at higher altitudes)
DO JL=KIDIA,KFDIA
  PTENC1(JL,JLEVBC(JL,1),IHNO3) =  ZTENBC(JL,1)  / PTSTEP
! store this special tendency into POUT(:,3,1) 
  POUT(JL,3,1)=PTSTEP * PTENC1(JL,JLEVBC(JL,1),IHNO3)*PDELP(JL,JLEVBC(JL,1)) / RG 
ENDDO
IF (IMODE_HNO3BC > 0) THEN
  !* additional constraints at 50  and 100 hpa, with 1-day nudging time step
  ZTAU_NUDGE=86400.
  DO JL=KIDIA,KFDIA
    PTENC1(JL,JLEVBC(JL,2),IHNO3) = ((1._JPRB - EXP(-PTSTEP / (ZTAU_NUDGE)))*ZTENBC(JL,2))  / PTSTEP
    POUT(JL,3,1)=POUT(JL,3,1)+PTSTEP * PTENC1(JL,JLEVBC(JL,2),IHNO3)*PDELP(JL,JLEVBC(JL,2)) / RG 
    ! Only overwrite extra-tropical tendencies (100hPa in tropics is stil troposphere)
    IF (ABS(ZLAT(JL)) > 30._JPRB) THEN
      PTENC1(JL,JLEVBC(JL,3),IHNO3) = ((1._JPRB - EXP(-PTSTEP / (ZTAU_NUDGE)))*ZTENBC(JL,3))  / PTSTEP
      POUT(JL,3,1)=POUT(JL,3,1)+PTSTEP * PTENC1(JL,JLEVBC(JL,3),IHNO3)*PDELP(JL,JLEVBC(JL,3)) / RG 
    ENDIF
  ENDDO

  !* final constraint on NO and NO2 at altitudes above 1 hPa adopting a 10-day nudging time step
  !* on top of existing diurnal cycle.
  ZTAU_NUDGE=864000.
  DO JL=KIDIA,KFDIA
    DO JK = 1,JLEVBC(JL,4) 
      PTENC1(JL,JK,INO2) = PTENC1(JL,JK,INO2) - PCEN(JL,JK,INO2) * (1._JPRB - EXP(-PTSTEP / (ZTAU_NUDGE))) / PTSTEP
      PTENC1(JL,JK,INO)  = PTENC1(JL,JK,INO)   - PCEN(JL,JK,INO) * (1._JPRB - EXP(-PTSTEP / (ZTAU_NUDGE))) / PTSTEP
    ENDDO
  ENDDO
ENDIF


! check for large values 
!   IF (LLCHECKMAX) THEN
!     DO JL=KIDIA,KFDIA
!       IF ( PTENC1(JL,JK,IHNO3) * PTSTEP > 1.0e-4 ) THEN 
!          print*, ' HNO3 large after bc ', ZCVM(JL,IHNO3) , PTENC1(JL,JK,IHNO3) * PTSTEP, ZCVM0(JL,IHNO3) , JL, JK 
!        ENDIF
!        IF ( PTENC1(JL,JK,INO2) * PTSTEP > 1.0e-4 ) THEN 
!           print*, ' NO2 large after bc ', ZCVM(JL,INO2) , PTENC1(JL,JK,INO2) * PTSTEP, ZCVM0(JL,INO2) , JL, JK 
!        ENDIF  
!     ENDDO
!   ENDIF   

!6.4 set NO2 and NO tendencies to ensure zero NOx in Stratopshere
IF (LCHEM_0NOX) THEN
 DO JK=1,KLEV
   DO JL=KIDIA,KFDIA
     ZPMAXO3CAR=PPTROPO(JL) 
     ZPMINO3CAR=ZPMAXO3CAR-4000_JPRB
     ZT0NO= -1.0_JPRB *  ( PCEN(JL,JK,INO) / PTSTEP )
     ZT0NO2= -1.0_JPRB *  ( PCEN(JL,JK,INO2) / PTSTEP )

     IF (PRSF1(JL,JK)<= ZPMAXO3CAR)  THEN
         ! At least partialy stratospheric tendencies
      IF (PRSF1(JL,JK)> ZPMINO3CAR)  THEN
         ! Transition region ZPMAX > P > ZPMIN
         ZPO3CAR_RATIO=(PRSF1(JL,JK)-ZPMINO3CAR)/(ZPMAXO3CAR-ZPMINO3CAR)
         PTENC1(JL,JK,INO) =  ZPO3CAR_RATIO* PTENC1(JL,JK,INO) +(1_JPRB-ZPO3CAR_RATIO)* ZT0NO  
         PTENC1(JL,JK,INO2) =  ZPO3CAR_RATIO* PTENC1(JL,JK,INO2) +(1_JPRB-ZPO3CAR_RATIO)* ZT0NO2  
      ELSE
        ! completely stratospheric region: P < ZPMIN
        PTENC1(JL,JK,INO) =  ZT0NO  
        PTENC1(JL,JK,INO2) =  ZT0NO2  
      ENDIF

    ENDIF
    
!    IF (PRSF1(JL,JK) <= ZPMAXO3CAR) THEN 
!      IF ( ZTENO3COR(JL,JK)*3600.0_JPRB > PCEN(JL,JK,IO3) * 0.01_JPRB )  & 
!   &  WRITE(NULOUT,*) ' O3Chem high tend ', JL, JK,  ZTENO3COR(JL,JK)*3600.0_JPRB, PCEN(JL,JK,IO3)  
!    ENDIF
   ENDDO
 ENDDO
ENDIF

!-----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHEM_TM5',1,ZHOOK_HANDLE )
END SUBROUTINE CHEM_TM5 
