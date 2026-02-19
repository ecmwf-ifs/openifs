! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

 SUBROUTINE CHEM_BASCOETM5 &
 &    (YDDIMV, YDMODEL,KSTEP, KIDIA  , KFDIA , KLON, KLEV , KAERO, &
 &     PTSTEP, PDELP, PRS1, PRSF1, PGEOH, PQP, PTP, &
 &     PLP, PIP, PAP, KLEVTROP, PALB, PWND, PLSM, PCSZA, PGELAT, &
 &     PGELAM, PGEMU,   PCEN , PTENC1,PBUDR, PBUDJ, PBUDX,  POUT, &
 &     PAEROP, PWETDIAM, PWETVOL,PND,PAERAOT, PAERAAOT, PAERASY,PSOGTOSOA, &
 &     PCHEM2GHG)


!**   DESCRIPTION
!     ----------
!
!   routine for C-IFS-BASCOE stratospheric chemistry merged with CB05 trop chem
!
!
!
!**   INTERFACE.
!     ----------
!          *CHEM_BASCOETM5* IS CALLED FROM *CHEM_MAIN*.

! INPUTS:
! -------
! KSTEP : Time step
! KIDIA :  Start of Array
! KFDIA :  End  of Array
! KLON  :  Length of Arrays
! KLEV  :  Number of Levels
! PTSTEP:  Time step in seconds
! PDELP(KLON,KLEV)            : PRESSURE DELTA in PRESSURE UNITES      (Pa)
! PRS1(KLON,0:KLEV)           : HALF-LEVEL PRESSURE           (Pa)
! PRSF1(KLON,KLEV)            : FULL-LEVEL PRESSURE           (Pa)
! PGEOH (KLON,KLEV)           :  Geopotential ??              (gpm ?)
! PQP     (KLON,KLEV)         :  SPECIFIC HUMIDITY            (kg/kg)
! PTP     (KLON,KLEV)         :  TEMPERATURE                  (K)
! PLP     (KLON,KLEV)         :  LCWC                         (kg/kg)
! PIP     (KLON,KLEV)         :  ICWC                         (kg/kg)
! PAP     (KLON,KLEV)         :  CLOD FRACTION                 0..1
! KLEVTROP  (KLON)            : Index for referring to tropopause level (humidity)      (-)
! PLSM    (KLON)              : land-sea-mask
! PALB(KLON)                  : Surface albedo
! PWND(KLON)                  : Surface wind
! PLSM(KLON)                  : Land Sea Mask albedo
! PCSZA(KLON)                 : COS of Solar Zenit Angle
! PGELAM(KLON)                : LONGITUDE (RADIANS)
! PGELAT(KLON)                : LATITUDE (RADIANS)
! PGEMU(KLON)                 : SINE OF LATITUDE
! PCEN(KLON,KLEV,NCHEM)       : CONCENTRATION OF TRACERS           (kg/kg)
! PAEROP(KLON,KLEV,KAERO)     : Aerosol concentrations  (kg/kg) - Note that fields are only non-zero if NACTAERO > 0
! PWETDIAM(KLON,KLEV,NMODE)   : Glomap geometric mean wet diameter per mode (real dims 1,1,1 if GLOMAP not used)
! PWETVOL(KLON,KLEV,NMODE)    : Glomap avg wet volume of size mode (m3) (real dims 1,1,1 if GLOMAP not used)
! PND(KLON,KLEV,NMODE)        : Glomap number concentration (cm-3) (real dims 1,1,1 if GLOMAP not used)
! PAERAOT(KLON,KLEV,6)        : Glomap extinction AOD per model level at 6 wavelengths
! PAERAAOT(KLON,KLEV,6)       : Glomap absorption AOD per model levelat 6 wavelengths
! PAERASY(KLON,KLEV,6)        : Glomap asymetry factor
!
! OUTPUTS:
! -------
! PTENC1  (KLON,KLEV,NCHEM)     : TENDENCY OF CONCENTRATION OF TRACERS BECAUSE OF CHEMISTRY (kg/kg s-1), no update
! PBUDR (KLON,KLEV,NCHEM)       : TENDENCIES DUE TO GAS-PHASE REACTIONS WITH OH (kg/kg/s)
! PBUDJ (KLON,KLEV,NPHOTO)      : TENDENCIES (loss) DUE TO PHOTOLYSIS (kg/kg/s)
! PBUDX(KLON,KLEV,NBUD_EXTRA)   : Extra chemical TENDENCIES (kg/kg/s)
! POUT (KLON,KLEV,5)            : additional output, e.g. UBC contribution , Photolysis rates O3 , NO2, tau for output
! PSOGTOSOA(KLON,KLEV,2)        : SOG to SOA conversion tendency
! PCHEM2GHG(KLON,KLEV,NCHEM2GHG): Information from chemistry to GHG.
!                                   1. atmospheric CH4 loss rate              [s-1]
!                                   2. tropospheric CO2 production tendency due to CO oxidation [kg CO2/kg/s]

!
! LOCAL:
! -------
!
! ZCVM0(KLON,NCHEM)       : initial volume ratios OF TRACERS           (molec/cm3)
! ZCVM (KLON,NCHEM+3)     : final   volume ratios OF TRACERS           (molec/cm3)
! ZCC (KPROMA,KFLEV)      : OVERHEAD CLOUD COVER
! ZTAUC (KPROMA,KFLEV)    : CLOUD OPTICAL THICKNESS
! ZTCTAUC (KPROMA)        : Total optical depth
! ZPCO3 (KPROMA,KFLEV)    : O3 column
!
!
!     AUTHOR.
!     -------
!        JOHANNES FLEMMING  *ECMWF*
!        VINCENT HUIJNEN    *KNMI*
!        Quentin Errera     *BIRA*
!        ORIGINAL : 2014-02-01

!     MODIFICATIONS.
!     --------------
!        TM5 -code structure                               : 2009-10-25
!        BASCOE-implementation of BASCOE scheme            : 2014-03-03
!        BASCOE sb15b J online                             : 2018-03-25
!
!     NOTES
!     -----
!   POUT is used for extra fields to save:
!       consistency of its usage should be carefully checked !
!
!-----------------------------------------------------------------------


USE TYPE_MODEL , ONLY : MODEL
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
! NCHEM : number of chemical species
! YCHEM : Data structure with Chemistry meta data
USE YOMCST   , ONLY : RD, RMD, RG , RPI , RNAVO, RMCO2
USE YOMRIP0  , ONLY : NINDAT
USE YOMLUN   , ONLY : NULERR, NULOUT

USE BASCOETM5_MODULE, ONLY : IO3, IO3S, ICO, IH2O, ICO2, INO2, ICH4, ISTRATAER, &
  &  INH3, INH4, INO3_A, ISO4, IMSA, IHO2, IOH, &
  &  IN2O5, IHCL, IHOCL, ICLONO2, IHOBR, IHBR, IBRONO2, IHNO3, &
  &  IAIR, IACID, NBC, BASCOE_BC, &
  &  ISOG1, ISOG2A, ISOG2B, &
  &  IN, INO

! BASCOE chemistry...
USE BASCOE_MODULE, ONLY : NHET, NBINS, NAER
USE BASCOE_LBC_MODULE   , ONLY : NLATBOUND_LBC, XLATBOUND_LBC,  &
  &                              MONTH_LBC, VALUES_LBC
USE BASCOE_J_MODULE, ONLY: NDISS, J_O3_O1D,J_NO2,J_H2O2_OH,J_HO2NO2_HO2,    &
  &  J_N2O5,J_CH2O_CO,J_CH2O_HCO,J_NO3_O,J_NO3_O2,J_O2_O,J_CH3OOH, J_HNO3
USE BASCOE_TUV_MODULE    , ONLY : mxwvn, NABSPEC, fbeamr, fbeamr2d, fbeam_dates, daily_solflux, &
  & AIR, O2ABS, O3ABS, NOABS, CO2ABS, NO2ABS
USE BASCOE_KPP_PARAMETERS, ONLY : NREACT_BASCOE=>NREACT,NVAR_BASCOE=>NVAR, NFIX_BASCOE=> NFIX
USE BASCOE_KPP_GLOBAL    , ONLY : RTOL_BASCOE=>RTOL,ATOL_BASCOE=>ATOL,ROUNDOFF_BASCOE=> ROUNDOFF_STORE

! General KPP settings
USE CIFS_KPP_INTPARAM    , ONLY : HMIN,HSTART,RTOLS_G,IAUTONOM,IROSMETH, VMR_BAD_LARGE

! TM5 chemistry ...
USE TM5_KPP_PARAMETERS, ONLY : NREACT_TM5=>NREACT, NVAR_TM5=>NVAR
USE TM5_KPP_GLOBAL    , ONLY : RTOL_TM5=>RTOL,ATOL_TM5=>ATOL, ROUNDOFF_TM5=> ROUNDOFF_STORE
USE TM5_CHEM_MODULE   , ONLY : NREAC,NBUD_EXTRA, KHO2L, KHO2_AER, KCOOH, NCHEM2GHG

USE TM5_PHOTOLYSIS    , ONLY : JO3D,JNO2,JH2O2,JHNO3,JHNO4,JN2O5,JACH2O,JBCH2O,JANO3, &
    &  JBNO3,JO2,JMEPE, NPHOTO,NBANDS_TROP,NGRID

! Glomap specifics
USE UKCA_MODE_SETUP, ONLY: NMODES

!-----------------------------------------------------------------------

IMPLICIT NONE

!*       0.1  ARGUMENTS
!             ---------

TYPE(TDIMV)    ,INTENT(IN)    :: YDDIMV
TYPE(MODEL)    ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN) :: KSTEP, KIDIA , KFDIA , KLON , KLEV, KAERO
REAL(KIND=JPRB),INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB),INTENT(IN)    :: PDELP(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PRSF1(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PRS1(KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PGEOH(KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PQP(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PTP(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PLP(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PIP(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PAP(KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT)   :: PTENC1(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NCHEM)
REAL(KIND=JPRB),INTENT(IN)    :: PCEN(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NCHEM)
REAL(KIND=JPRB),INTENT(IN)    :: PCSZA(KLON)
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVTROP(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PALB(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PWND(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PLSM(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PGELAT(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PGELAM(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PGEMU(KLON)
REAL(KIND=JPRB),INTENT(OUT)   :: PBUDJ(KLON,KLEV,NPHOTO)
REAL(KIND=JPRB),INTENT(OUT)   :: PBUDR(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NCHEM)
REAL(KIND=JPRB),INTENT(OUT)   :: PBUDX(KLON,KLEV,NBUD_EXTRA)
REAL(KIND=JPRB),INTENT(OUT)   :: POUT(KLON,KLEV,5)
REAL(KIND=JPRB),INTENT(IN)    :: PAEROP(KLON,KLEV,KAERO)
REAL(KIND=JPRB),INTENT(IN)    :: PWETDIAM(KLON,KLEV,NMODES)
REAL(KIND=JPRB),INTENT(IN)    :: PWETVOL(KLON,KLEV,NMODES)
REAL(KIND=JPRB),INTENT(IN)    :: PND(KLON,KLEV,NMODES)
REAL(KIND=JPRB),INTENT(IN)    :: PAERAOT(KLON,KLEV,6)
REAL(KIND=JPRB),INTENT(IN)    :: PAERAAOT(KLON,KLEV,6)
REAL(KIND=JPRB),INTENT(IN)    :: PAERASY(KLON,KLEV,6)
REAL(KIND=JPRB),INTENT(OUT)   :: PSOGTOSOA(KLON,KLEV,2)
REAL(KIND=JPRB),INTENT(OUT)   :: PCHEM2GHG(KLON,KLEV,NCHEM2GHG)

!*       0.5   LOCAL VARIABLES
!              ---------------

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! * Lat /Lon time
REAL(KIND=JPRB) , DIMENSION(KLON)   :: ZLAT
REAL(KIND=JPRB) , DIMENSION(KLON)   :: ZLON
INTEGER(KIND=JPIM)                  :: IMONTH0, IDAY0, IYEAR0, ILMONTHS(12)
INTEGER(KIND=JPIM)                  :: IMONTH, IDAY, IYEAR, IYYYYMM, IYYYYMMDD, IYEAR_CH4, IJUL, IJULMAX
REAL(KIND=JPRB)                     :: ZJUL1, ZUT

! * counters
INTEGER(KIND=JPIM) :: JK, JL, JT, JLEV, JB
!INTEGER(KIND=JPIM) :: IDRYDEP
INTEGER(KIND=JPIM) :: ISSO4, ISNH4, ISNO3
INTEGER(KIND=JPIM) :: ISSOA1, ISSOA2

! * chemical data
REAL(KIND=JPRB) , DIMENSION(KLON,YDMODEL%YRML_GCONF%YGFL%NCHEM+3)   :: ZCVM
REAL(KIND=JPRB) , DIMENSION(KLON,YDMODEL%YRML_GCONF%YGFL%NCHEM)     :: ZCVM0
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV)                              :: ZDENS
REAL(KIND=JPRB) , DIMENSION(KLON)                :: ZAIRDM
REAL(KIND=JPRB)                                  :: ZAIRDM1
REAL(KIND=JPRD) :: ZDENS_DP

! * METEO-INFO; should in final version come from IFS
!REAL(KIND=JPRB)     :: ZHGT(KLON)       ! geopotential layer bottom

! * Photolysis data: cloud/ozone info, actinic fluxes, photolysis rates
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV)   :: ZPCO3
REAL(KIND=JPRB)                          :: ZCOLO3(KLON,0:KLEV)
REAL(KIND=JPRB)                          :: ZCOLO3_DU(KLON,KLEV)
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV)   :: ZCC
REAL(KIND=JPRB) , DIMENSION(KLON)        :: ZTCTAUC
REAL(KIND=JPRB) , DIMENSION(KLON)        :: ZHPLUS

REAL(KIND=JPRB) , DIMENSION(KLON,KLEV,NBANDS_TROP,NGRID) :: &
      &                                       ZTAUA_AER, ZTAUS_AER, ZPMAER
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV)   ::   ZTAUA_CLD, ZTAUS_CLD, ZPMCLD
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV)   ::   ZCLOUD_REFF
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV)   ::   ZSAD_AER, ZSAD_CLD, ZSAD_ICE
!   ...correction factor for photo rates in function of sun-earth distance
REAL(KIND=JPRB)                          :: ZDISTFAC
!   ... stratosphere: absorbing species, lowest level to compute, solar flux
REAL(KIND=JPRB), DIMENSION(KLEV,NABSPEC) :: ZDENSA
INTEGER(KIND=JPIM)                       :: JBOTJ
REAL(KIND=JPRB), DIMENSION(MXWVN)        :: ZFBEAMR
INTEGER(KIND=JPIM)                       :: IDXDATE(1)
!   ...blending zone for strato-tropo photolysis rates

INTEGER(KIND=JPIM), PARAMETER                 :: JPNLBLEND = 4         !number of levels under tropopause
INTEGER(KIND=JPIM), PARAMETER                 :: INJ_STRAT_TROP = 12   ! number of photo rates
INTEGER(KIND=JPIM), DIMENSION(INJ_STRAT_TROP) :: IBASCOE_DIS, ITM5_DIS    ! corresponding indices

! * reaction rates - troposphere
REAL(KIND=JPRB) , DIMENSION(NREAC)            :: ZRR
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV,NPHOTO) :: ZRJ
REAL(KIND=JPRB) , DIMENSION(NPHOTO)           :: ZRJ_IN

!* Saturation pressure for water vapour
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV) :: ZQSAT, ZRHCL

! * budget accumulators
REAL(KIND=JPRB) , DIMENSION(KLON,NPHOTO) :: ZCR2
REAL(KIND=JPRB) , DIMENSION(KLON,NREAC)  :: ZCR3

! * EQSAM input / output parameters
REAL(KIND=JPRB)                          :: ZRH,ZCCS
REAL(KIND=JPRB)                          :: ZNH,ZNO3,ZSO4
REAL(KIND=JPRB) , DIMENSION(4)           :: ZYEQ


REAL(KIND=JPRB)                          :: ZTAU_NUDGE
REAL(KIND=JPRB)                          :: ZRGI, ZDELP , ZCONC_HO2
INTEGER(KIND=JPIM)                       :: IERR, IFLAG

! BASCOE Conversions  : 1.255e+21 *.767e-17   from kg/m2 --> mol/cm2 --> DU
!REAL(KIND=JPRB), PARAMETER   :: ZTODU = 47296. ! O3 from kg/m2 --> DU
REAL(KIND=JPRB)                           :: ZCST  ! molec/cm2 -> DU
REAL(KIND=JPRB)                           :: ZVMRO3,ZSZA, ZHGT, ZH2O_ECMWF
REAL(KIND=JPRB), DIMENSION(KLON,KLEV,NDISS):: ZAJVAL
REAL(KIND=JPRB), DIMENSION(NHET)          :: ZRHET
REAL(KIND=JPRB)                           :: ZALB,ZTROPO_FAC

!KPP related
INTEGER(KIND=JPIM),DIMENSION(20)            :: ICNTRL, ISTATUS
REAL(KIND=JPRD),   DIMENSION(20)            :: ZCNTRL, ZCNTRL_P, ZSTATE

!Stratosphere
REAL(KIND=JPRD),   DIMENSION(NREACT_BASCOE) :: ZRCONST_BASCOE
REAL(KIND=JPRD),   DIMENSION(NVAR_BASCOE)   :: ZVAR_BASCOE
REAL(KIND=JPRD),   DIMENSION(NFIX_BASCOE)   :: ZFIX_BASCOE

! Troposphere...
REAL(KIND=JPRD),   DIMENSION(NREACT_TM5)     :: ZRCONST_TM5
REAL(KIND=JPRD),   DIMENSION(NVAR_TM5)       :: ZVAR_TM5

! Strat. PSC / aerosol related
LOGICAL,            DIMENSION(KLON)        :: LL_PSC_POSSIBLE
INTEGER(KIND=JPIM), DIMENSION(KLON)        :: JBOT_PSC,JTOP_PSC
REAL(KIND=JPRB),    DIMENSION(KLON,KLEV,NBINS)   :: ZSA_SIZEDIST
REAL(KIND=JPRB),    DIMENSION(KLON,KLEV,NAER)    :: ZAER
REAL(KIND=JPRB),    DIMENSION(KLON,KLEV)  :: ZAER_INFO
INTEGER(KIND=JPIM), DIMENSION(KLON)       :: JTROPOP
INTEGER(KIND=JPIM), DIMENSION(9)          :: JHET_TRACER
INTEGER(KIND=JPIM), DIMENSION(2)          :: JPSC_TRACER
INTEGER(KIND=JPIM), DIMENSION(3)          :: JAER_TRACER
INTEGER(KIND=JPIM), DIMENSION(2)          :: JRATE_TRACER

! Indices to slice GLOMAP arrays with. These must be forced to
! 1 if GLOMAP is not in use because the arrays are then not
! fully allocated. (This is only necessary because we're calling
! TM5_CALRATES individually on each grid point. If we passed the
! whole array like CHEM_TM5 does, it would never be accessed
! unless GLOMAP was active, so not be an issue.)
INTEGER(KIND=JPIM)                        :: IL_GLOMAP, IK_GLOMAP, IMODES
LOGICAL                                   :: LLGLOMAP

! Boundary condition (CH4 @ surface, BASCOE)
REAL(KIND=JPRB) , DIMENSION(KLON)          :: ZTENBC
INTEGER(KIND=JPIM)                         :: JBC
REAL(KIND=JPRB), DIMENSION(NLATBOUND_LBC-1):: ZBCVAL      ! for every latitude band
INTEGER(KIND=JPIM)                         :: IMONTH_LBC(1),IMONTH_LBC_H, JLAT_LBC
! LOGICAL                                    :: LLBC_PREINDUST

! Variable to check NOy mass tendency
!REAL(KIND=JPRB)                           :: ZNOYTEND

! Variables used to compute altitude
REAL(KIND=JPRB)       :: ZPSURF_STD,ZSURF_H, ZTHKNESS,ZHGT_BASCOE(KLON,KLEV)

LOGICAL                                   :: LLCOD_TM5, LLSTRATAIR
! Switch for selection of aerosol parameterization for photolysis
INTEGER(KIND=JPIM), PARAMETER             :: ITAU_MACC=0
! Range for overwriting H2O in tropopause region.
!   - Should be approx. 4 model levels for 60-level IFS version
!   - Should be approx. 8 model levels for 137-level IFS version
INTEGER(KIND=JPIM)                        :: IRANGE_TROPOP

!Help variable
REAL(KIND=JPRB)                           :: ZFAC

! Variables needed for SOA computation
INTEGER(KIND=JPIM), PARAMETER             :: INSOG=3
INTEGER(KIND=JPIM), PARAMETER             :: INSOA=2
REAL(KIND=JPRB)                           :: ZORGAERO, ZRHO
REAL(KIND=JPRB), DIMENSION(INSOG)         :: ZSOG
REAL(KIND=JPRB), DIMENSION(INSOA)         :: ZSOGH, ZSOA
INTEGER(KIND=JPIM)                        :: INBDU,INBOM
REAL(KIND=JPRB)                           :: ZXLSOA,ZSOA_TMP,ZJSOA
INTEGER(KIND=JPIM), DIMENSION(INSOG)      :: JSOG_TRACER
INTEGER(KIND=JPIM), DIMENSION(INSOA)      :: JSOA_TRACER
! ------------------------------------------------------------------
#include "fcttim.func.h"
!-------------------------------------------------------------------
#include "satur.intfb.h"
#include "bascoe_gs_liq.intfb.h"
! #include "bascoe_zenith_fct.intfb.h"
#include "bascoe_hetconst.intfb.h"
#include "bascoe_j_interp.intfb.h"
#include "bascoe_j_calc.intfb.h"
#include "bascoe_psc_param.intfb.h"
#include "bascoe_psc_possible.intfb.h"
#include "bascoe_tropopause.intfb.h"
!#include "bascoetm5_noymass.intfb.h"
! KPP code - STRAT
#include "bascoe_kpp_rates.intfb.h"
#include "bascoe_kpp_initialize.intfb.h"
#include "bascoe_kpp_integrator.intfb.h"
#include "bascoe_kpp_update_cifs_conc.intfb.h"
#include "cifs_kpp_wlamch.intfb.h"
! KPP code - TROP
#include "tm5_kpp_rates.intfb.h"
#include "tm5_kpp_initialize.intfb.h"
#include "tm5_kpp_integrator.intfb.h"
#include "tm5_kpp_update_cifs_conc.intfb.h"
! TM5 code
#include "tm5_aerosol_info.intfb.h"
#include "tm5_boundary_ch4.intfb.h"
#include "tm5_calrates.intfb.h"
#include "tm5_eqsam.intfb.h"
#include "tm5_glomap_aerosol.intfb.h"
#include "tm5_macc_aerosol.intfb.h"
#include "tm5_ibud.intfb.h"
#include "tm5_o3s.intfb.h"
#include "tm5_rbud.intfb.h"
#include "tm5_soa.intfb.h"
#include "tm5_photo_flux.intfb.h"
#include "tm5_slingo.intfb.h"
! #include "tm5_sundis.intfb.h"
#include "sundistcorr.intfb.h"
#include "tm5_wetchem_point.intfb.h"
#include "tm5_stratoloss.intfb.h"
#include "cod_op_tm5.intfb.h"

IF (LHOOK) CALL DR_HOOK('CHEM_BASCOETM5',0,ZHOOK_HANDLE )
ASSOCIATE(YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO,YDERDI=>YDMODEL%YRML_PHY_RAD%YRERDI, &
 & YDCHEM=>YDMODEL%YRML_CHEM%YRCHEM, YDEAERSNK=>YDMODEL%YRML_PHY_AER%YREAERSNK, &
 & YDEAERSRC=>YDMODEL%YRML_PHY_AER%YREAERSRC, YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM, &
 & YDECLD=>YDMODEL%YRML_PHY_EC%YRECLD, &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL, YDRIP=>YDMODEL%YRML_GCONF%YRRIP, &
 & LPHYLIN=>YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN)
ASSOCIATE(NACTAERO=>YGFL%NACTAERO, NCHEM=>YGFL%NCHEM, NCHEM_DV=>YGFL%NCHEM_DV, &
 & YCHEM=>YGFL%YCHEM, LAERCHEM=>YGFL%LAERCHEM,LAERSOA=>YDCOMPO%LAERSOA, &
 & LAERSOA_COUPLED=>YDCOMPO%LAERSOA_COUPLED, &
 & AERO_SCHEME=>YDCOMPO%AERO_SCHEME, &
 & YEMIS2D_DESC=>YDCOMPO%YEMIS2D_DESC, NEMIS2D_DESC=>YDCOMPO%NEMIS2D_DESC, &
 & NTYPAER=>YDEAERATM%NTYPAER, &
 & LAERNITRATE => YDCOMPO%LAERNITRATE, &
 & LCHEM_ANACH4=>YDCHEM%LCHEM_ANACH4, &
 & LCHEM_WEAK_CH4_RELAXATION=>YDCHEM%LCHEM_WEAK_CH4_RELAXATION, &
 & LCHEM_TROPO=>YDCOMPO%LCHEM_TROPO, &
! & LGHG_CHEMTEND_CH4=>YDCOMPO%LGHG_CHEMTEND_CH4, &
 & LCHEM_JOUT=>YDCHEM%LCHEM_JOUT, LCHEM_AEROI=>YDCHEM%LCHEM_AEROI, &
 & LCHEM_DIAC=>YDCHEM%LCHEM_DIAC, KCHEM_YEARPI=>YDCHEM%KCHEM_YEARPI, &
 & REPSEC=>YDECLD%REPSEC, &
 & REPCLC=>YDERDI%REPCLC, &
 & NSTADD=>YDRIP%NSTADD, RHGMT=>YDRIP%RHGMT)
!-----------------------------------------------------------------------
! chemistry scheme name - this will later also come from external input
!-----------------------------------------------------------------------


! Preparation for kpp-solver.
! Set kpp parameters to default, taken from cifs_kpp_IntParam module
! See comments in Integrator module for a list of the defaults.

RTOL_BASCOE(1:NVAR_BASCOE) = 0.05_JPRB * RTOLS_G
ATOL_BASCOE(1:NVAR_BASCOE) = 10._JPRB  !  was 1.e-16*cfactor before v3s04 or v3d06

RTOL_TM5(1:NVAR_TM5) = 0.5_JPRB * RTOLS_G
ATOL_TM5(1:NVAR_TM5) = 10._JPRB  !  was 1.e-16*cfactor before v3s04 or v3d06

ICNTRL(:) = 0_JPIM
ICNTRL(1) = IAUTONOM
! Change some parameters from the default to new values
! Select Integrator
!    ICNTRL(3)  -> selection of a particular method.
!               For Rosenbrock, options are:
!        = 0 :  default method is Rodas3
!        = 1 :  method is  Ros2
!        = 2 :  method is  Ros3
!        = 3 :  method is  Ros4
!        = 4 :  method is  Rodas3
!        = 5 :  method is  Rodas4

! ----------------------------------------------------------------------
!  Set Integrator input parameters. Values set in chem_IntParam.f90
! ----------------------------------------------------------------------
ICNTRL(3) = IROSMETH
!VH ICNTRL(4) = 200 ! set max. no of steps ?
ICNTRL(7) = 1 ! Currently no adjoint

ZCNTRL(:) = 0._JPRB
ZCNTRL(1) = HMIN
ZCNTRL(2) = PTSTEP
ZCNTRL(3) = HSTART
!VH ZCNTRL(3) = PTSTEP

! Preferred time step - should be made flexible (see HSAVE_KPP)
!ZSTATE(Nhexit) = 0._JPRB

! Set range for tropopause 4 for 60 and 91 level version, 8 for 137 level version
IRANGE_TROPOP=4
IF (KLEV > 130) THEN
  IRANGE_TROPOP=8
ENDIF


POUT(KIDIA:KFDIA,:,:)  = 0.0_JPRB
PSOGTOSOA(KIDIA:KFDIA,:,:) = 0.0_JPRB

! Lat / Lon
DO JL=KIDIA,KFDIA
  ZLAT(JL)=(180.0_JPRB/RPI)*PGELAT(JL)
  ZLON(JL)=(180.0_JPRB/RPI)*PGELAM(JL)
ENDDO

! Initialize appropriate tracers for strat. heterogeneous and psc chemistry
JHET_TRACER(1)= IH2O
JHET_TRACER(2)= IN2O5
JHET_TRACER(3)= IHCL
JHET_TRACER(4)= IHOCL
JHET_TRACER(5)= ICLONO2
JHET_TRACER(6)= IHOBR
JHET_TRACER(7)= IHBR
JHET_TRACER(8)= IBRONO2
JHET_TRACER(9)= IHNO3

JPSC_TRACER(1)= IH2O
JPSC_TRACER(2)= IHNO3

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


! Time stuff: Start time...
IYEAR=NCCAA(NINDAT)
IMONTH=NMM(NINDAT)
IDAY=NDD(NINDAT)
! Start Jul day + number of days since start of model
IYEAR0=NCCAA(NINDAT)
IMONTH0=NMM(NINDAT)
IDAY0=NDD(NINDAT)
! ...current date...
CALL UPDCAL (IDAY0,IMONTH0,IYEAR0,YDRIP%NSTADD,IDAY,IMONTH,IYEAR,ILMONTHS,-1)
IYYYYMM=IYEAR*100+IMONTH
IYYYYMMDD=IYYYYMM*100+IDAY
! ...day of year, number of days in year
ZJUL1=RJUDAT(IYEAR,1,1)
IJUL=NINT(RJUDAT(IYEAR,IMONTH,IDAY) - ZJUL1) + 1
IJULMAX=NINT(RJUDAT(IYEAR,12,31) - ZJUL1) + 1

! RHGMT: GMT time of model - between 0 and 86400.
! ZUT: Time of day in hours.
ZUT = RHGMT /3600.

! Compute correction factor for photo rates
CALL SUNDISTCORR(IMONTH,IDAY,ZDISTFAC)

! Compute here 'roundoff' number - there is a paralellization
! issue when calling WLMACH as part of kpp-code.
IF ( KSTEP == 0_JPIM ) THEN
    CALL CIFS_KPP_WLAMCH(ROUNDOFF_BASCOE , 'E')
    ROUNDOFF_TM5=ROUNDOFF_BASCOE
ENDIF

! 1.0 Initialize output tendencies to 0.
PTENC1(KIDIA:KFDIA,:,:)=0.0_JPRB

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

! Set a few matching photolysis rates, for merging
!
IBASCOE_DIS = (/J_O3_O1D,J_NO2,J_H2O2_OH,J_HNO3,J_HO2NO2_HO2,J_N2O5,J_CH2O_CO,J_CH2O_HCO,J_NO3_O,J_NO3_O2,J_O2_O, J_CH3OOH/)
ITM5_DIS    = (/JO3D ,JNO2 ,JH2O2 ,JHNO3 ,JHNO4 ,JN2O5 ,JACH2O ,JBCH2O ,JANO3 ,JBNO3 ,JO2 , JMEPE   /)


IF (.NOT. YDCHEM%LCHEM_BASCOE_JON) THEN
  ! BASCOE variant to compute overhead column [DU]
  !-----------------------------------------------------------------------
  ! Calculate overhead ozone columns *at* levels:  Zcst * SUM( vmr * delta_p )
  ! where delta_p is between levels and vmr are at mid-levels
  !-----------------------------------------------------------------------
  ZCST  = (1.0E-4/RG)*( RNAVO / (1.0E-3* RMD) ) * 1.0E3 / 2.687E19  ! molec/cm2 -> DU
  DO JL=KIDIA,KFDIA
    ! convert to mixing ratio
    ZVMRO3=MAX(PCEN(JL,1,IO3) / YCHEM(IO3)%RMOLMASS *RMD ,0._JPRB)
    ! convert to DU
    ZCOLO3_DU(JL,1)  = ZCST * ZVMRO3 * ( PRSF1(JL,1) - 0._JPRB )
  ENDDO

  DO JLEV=2,KLEV
    DO JL=KIDIA,KFDIA
      ! convert to mixing ratio
      ZVMRO3=MAX(0.5_JPRB*(PCEN(JL,JLEV,IO3)+PCEN(JL,JLEV-1,IO3)) / YCHEM(IO3)%RMOLMASS *RMD ,0._JPRB)
      ! convert to DU
      ZCOLO3_DU(JL,JLEV)=ZCOLO3_DU(JL,JLEV-1)+ZCST * ZVMRO3* ( PRSF1(JL,JLEV) - PRSF1(JL,JLEV-1))
    ENDDO
  ENDDO
ENDIF

! BASCOE way to compute model altitude, first surface model level:
ZPSURF_STD=101325._JPRB ! std p at surf (Pa)
DO JL=KIDIA,KFDIA

  IF( PRS1(JL,KLEV) < ZPSURF_STD ) THEN
    ZSURF_H = 7._JPRB*LOG( ZPSURF_STD / PRS1(JL,KLEV) )
  ELSE
    ZSURF_H=0.0_JPRB
  ENDIF
  ZTHKNESS = PTP(JL,KLEV)*287./9.806*LOG(PRS1(JL,KLEV)/PRSF1(JL,KLEV))
  ZHGT_BASCOE(JL,KLEV) = ZSURF_H + 1.E-3*ZTHKNESS
ENDDO

DO JLEV=KLEV-1,1,-1
  DO JL=KIDIA,KFDIA
     ZTHKNESS=0.5*(PTP(JL,JLEV+1)+PTP(JL,JLEV))*287./9.806*&
     &                 LOG(PRSF1(JL,JLEV+1)/PRSF1(JL,JLEV))
     ZHGT_BASCOE(JL,JLEV)=ZHGT_BASCOE(JL,JLEV+1)+1E-3*ZTHKNESS
  ENDDO
ENDDO

!--------
! Find tropopause level
CALL BASCOE_TROPOPAUSE(KIDIA,KFDIA,KLON,KLEV,PTP,PRSF1,PGEOH,ZLAT,JTROPOP)


! ------ GET SAD climatology -

IF (KSTEP /= 0_JPIM ) THEN
  ! ZAER_INFO will be used in BASCOE_GS_LIQ
  ZAER_INFO(KIDIA:KFDIA,1:KLEV)=PCEN(KIDIA:KFDIA,1:KLEV,ISTRATAER)
ELSE
  ! ZAER_INFO will be updated in BASCOE_GS_LIQ
  ZAER_INFO(KIDIA:KFDIA,1:KLEV)=0._JPRB
ENDIF

CALL BASCOE_GS_LIQ(KSTEP, IMONTH, KIDIA, KFDIA, KLON, KLEV, JTROPOP, PRSF1, ZLAT, PTP, ZAER, ZSA_SIZEDIST, ZAER_INFO)

IF (KSTEP == 0_JPIM) THEN
  ! Create 'tendency' to arrive at aerosol field
  PTENC1(KIDIA:KFDIA,1:KLEV,ISTRATAER)=(ZAER_INFO(KIDIA:KFDIA,1:KLEV) - PCEN(KIDIA:KFDIA,1:KLEV,ISTRATAER))/ PTSTEP
ENDIF

! --------------------
!  Find range where PSC can be present
! --------------------

CALL BASCOE_PSC_POSSIBLE(KIDIA, KFDIA, KLON, KLEV,IJUL,PTP,PRSF1,ZLAT,&
                        & JTROPOP,LL_PSC_POSSIBLE,JTOP_PSC,JBOT_PSC)

! ---------------------
! initialize air densities (molec/cm3)

DO JLEV=1,KLEV
  DO JL=KIDIA,KFDIA
    ZDENS(JL,JLEV) = 7.24291e16_JPRB*PRSF1(JL,JLEV)/PTP(JL,JLEV)
  ENDDO
ENDDO
! ---------------------
!  Compute strato photolysis rates
!   but ONLY where needed, to avoid unnecessary processing:
!   - up to level where strato-tropo blending of photo rates is applied
!   - for solar zenith angle < 96 degrees)

ZAJVAL(KIDIA:KFDIA,:,:) = 0._JPRB
IF (YDCHEM%LCHEM_BASCOE_JON) THEN
  ! Use online calculations...
  !
  !   ...obtain solflux for current date if daily solflux is requested...
  !
  IF (daily_solflux) THEN
    ! check if date in range
    !
    IF ( IYYYYMMDD < fbeam_dates(1) .OR. IYYYYMMDD > fbeam_dates(SIZE(fbeam_dates)) ) THEN
      WRITE(NULOUT,*) 'CHEM_BASCOETM5: error solflux date ',fbeam_dates(1),'-',fbeam_dates(SIZE(fbeam_dates)),&
      &               ' outside range. current date ',IYYYYMMDD
      CALL ABOR1(' error: date outside range of solflux dates read from file')
    ENDIF
    ! get closest date
    !
    IDXDATE = MINLOC( ABS(fbeam_dates - IYYYYMMDD ))
    IF ( fbeam_dates(IDXDATE(1)) /= IYYYYMMDD ) THEN
        WRITE(NULOUT,*) 'CHEM_BASCOETM5: warning solflux date ',fbeam_dates(IDXDATE(1)), &
        &               ' differs from current date ',IYYYYMMDD
    ENDIF
    ZFBEAMR(:) = FBEAMR2D(IDXDATE(1),:)
  ELSE
    ZFBEAMR = FBEAMR
  ENDIF

  DO JL=KIDIA,KFDIA
    ZSZA = ACOS(PCSZA(JL))*180_JPRB/RPI

    !VH propose to take over BASCOE version (double check that it works for all (FC-)times !
    ! CALL BASCOE_ZENITH_FCT(ZLAT(JL),ZLON(JL),IJUL,IJULMAX,ZUT,ZSZA)

    IF( ZSZA < 96._JPRB ) THEN
      JBOTJ = JTROPOP(JL) + JPNLBLEND       ! lowest level to compute photolysis rates

      ! compute number densities profiles for the absorbing species
      DO JK = 1, JBOTJ
        ZDENSA(JK,AIR)   = ZDENS(JL,JK)
        ZDENSA(JK,O2ABS) = ZDENS(JL,JK) * 0.209_JPRB
        ZAIRDM1 = ZDENS(JL,JK) * RMD
        ZDENSA(JK,O3ABS) = PCEN(JL,JK,IO3)  / YCHEM(IO3)%RMOLMASS  *ZAIRDM1
        ZDENSA(JK,NOABS) = PCEN(JL,JK,INO)  / YCHEM(INO)%RMOLMASS  *ZAIRDM1
        ZDENSA(JK,CO2ABS)= PCEN(JL,JK,ICO2) / YCHEM(ICO2)%RMOLMASS *ZAIRDM1
        ZDENSA(JK,NO2ABS)= PCEN(JL,JK,INO2) / YCHEM(INO2)%RMOLMASS *ZAIRDM1
      ENDDO

      ! IF (LPDEFAULTALB) THEN
      !   ZALB = PPALB
      ! ELSE
      ZALB = PALB(JL)
      ! ENDIF

      CALL BASCOE_J_CALC( JBOTJ, ZHGT_BASCOE(JL,1:JBOTJ), ZSZA, PTP(JL,1:JBOTJ), &
                      &   ZDENSA(1:JBOTJ,1:NABSPEC), ZALB, ZFBEAMR, ZAJVAL(JL,1:JBOTJ,1:NDISS) )
    ENDIF
  ENDDO
ELSE
  ! use offline (faster lookup-table) approach
  DO JL=KIDIA,KFDIA
    ZSZA = ACOS(PCSZA(JL))*180_JPRB/RPI

    !VH propose to take over BASCOE version (double check that it works for all (FC-)times !
    ! CALL BASCOE_ZENITH_FCT(ZLAT(JL),ZLON(JL),IJUL,IJULMAX,ZUT,ZSZA)

    IF( ZSZA < 96._JPRB ) THEN
      JBOTJ =  JTROPOP(JL) + JPNLBLEND       ! lowest level to compute photolysis rates

      ! compute number densities profiles for the absorbing species
      DO JK = 1, JBOTJ
        !VH ZHGT = PGEOH(JL,JK-1)  * ZRGI *1e-3 ! height in km
        ! ZHGT = 0.5*(PGEOH(JL,JK-1) + PGEOH(JL,JK))  * ZRGI *1e-3 ! height in km
        !
        ZHGT = ZHGT_BASCOE(JL,JK)

        ! Set maximimum height to ~110 km altitude
        ZHGT = MIN(ZHGT,109.9_JPRB)
        ! ZCOLO3_DU: o3 overhead column in DU, converted from kg/m2
        ! compute photolysis rates (ZAJVAL)
        CALL BASCOE_J_INTERP( ZSZA, ZHGT  , ZCOLO3_DU(JL,JK), ZAJVAL(JL,JK,1:NDISS) )
      ENDDO
    ENDIF
  ENDDO
ENDIF


! adjust photolysis rate to sun earth distance variation ? (as in troposphere)
ZAJVAL(KIDIA:KFDIA,1:KLEV,1:NDISS) = ZAJVAL(KIDIA:KFDIA,1:KLEV,1:NDISS) * ZDISTFAC

! --------------------



! 1.2  Calculate integrated cloud cover above level - adapted from cldpp.f90

ZCC(KIDIA:KFDIA,1:KLEV)=0.0_JPRB

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

! 1.3 Initialize budget accumulation...
PBUDJ(KIDIA:KFDIA,:,:) = 0.0_JPRB
PBUDR(KIDIA:KFDIA,:,:) = 0.0_JPRB
PBUDX(KIDIA:KFDIA,:,:) = 0.0_JPRB

! 1.4 calculate cloud optical depth
! please note that the arguments to cod_op range from 0:KLEV


LLCOD_TM5=.TRUE.
IF ( LLCOD_TM5 ) THEN
! * IFS scheme - modified (optimized) for TM5 format...
  CALL COD_OP_TM5(YDDIMV,YDMODEL%YRML_PHY_RAD%YRERAD,KIDIA,KFDIA,KLON,KLEV,1,PQP,PTP,PAP,PRS1,PRSF1,PLSM, &
    & PWND,PLP,PIP,ZCLOUD_REFF,ZTCTAUC,ZTAUS_CLD,ZTAUA_CLD,ZPMCLD)
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
CALL TM5_PHOTO_FLUX(KIDIA, KFDIA, KLON, KLEV, PTP, PCSZA, PALB,&
   &          ZPCO3, PRSF1, PRS1,&
   &          ZTAUA_CLD,ZTAUS_CLD,ZPMCLD,&
   &          ZTAUA_AER,ZTAUS_AER,ZPMAER,&
   &          PAP,PGEOH, ZRJ )

! 1.6 adjust photolysis rate to sun earth distance variation

! CALL TM5_SUNDIS(KIDIA,KFDIA,KLON,KLEV,IMONTH,IDAY, ZRJ)
ZRJ(KIDIA:KFDIA,1:KLEV,1:NPHOTO) = ZRJ(KIDIA:KFDIA,1:KLEV,1:NPHOTO) * ZDISTFAC

! 1.6.1 special case for extension of JNO2 to stratosphere
!       -> copy tropo NO2 photo rate from TM5 (ZRJ) to BASCOE (ZAJVAL)
IF (YDCHEM%LCHEM_EXTENDJNO2) THEN
  ZAJVAL(KIDIA:KFDIA,1:KLEV,J_NO2)= ZRJ(KIDIA:KFDIA,1:KLEV,JNO2)
ENDIF

ZSAD_AER(KIDIA:KFDIA,1:KLEV)=0.0_JPRB
ZSAD_CLD(KIDIA:KFDIA,1:KLEV)=0.0_JPRB
ZSAD_ICE(KIDIA:KFDIA,1:KLEV)=0.0_JPRB

! Default value used for slicing dummy GLOMAP arrays.
! Will be overridden with the real one if GLOMAP is in use.
IK_GLOMAP=1
IL_GLOMAP=1
IMODES=1

LLGLOMAP = (NACTAERO > 0 .AND. TRIM(AERO_SCHEME)=="glomap")
IF (LLGLOMAP) THEN



  CALL ABOR1("OIFS - GLOMAP should never be called, EXIT")

ENDIF

!2.0 loop over levels for solving chemistry
DO JK=1,KLEV


  !2.01 set ZCR2 / ZCR3 budget accumulators for this level to zero.
  IF (LCHEM_DIAC) THEN
    ZCR2(KIDIA:KFDIA,1:NPHOTO) = 0._JPRB
    ZCR3(KIDIA:KFDIA,1:NREAC) = 0._JPRB
  ENDIF

  IF (LLGLOMAP) IK_GLOMAP=JK

  !2.1 Loop over lon/lat
  DO JL=KIDIA,KFDIA

    IF (LLGLOMAP) IL_GLOMAP=JL

    !VH Computation from PCSZA
    ZSZA = ACOS(PCSZA(JL))*180_JPRB/RPI

    ! Option to use BASCOE version (double check that it works for all (FC-)times !
    ! CALL BASCOE_ZENITH_FCT(ZLAT(JL),ZLON(JL),IJUL,IJULMAX,ZUT,ZSZA)

    ! only in those situations compute strat. photolysis rates (for efficiency)
    IF ( JK < JTROPOP(JL)+JPNLBLEND ) THEN

      ! Blend with tropospheric photolysis rates within range
      ! Only merge rates when SZA<94
      IF ( ZSZA < 94._JPRB) THEN

         IF (JK > JTROPOP(JL) - JPNLBLEND .AND. JK <= JTROPOP(JL) ) THEN

          ! Compute factor: 1 at top, JPNLBLEND=4 levels above tropopause interface (low pressure end)
          ! 0 at bottom of interface (high pressure end, tropopause interface)
          ZTROPO_FAC=FLOAT(JTROPOP(JL)-JK)/FLOAT( JPNLBLEND  )
          ZTROPO_FAC=MAX(0._JPRB,MIN(1._JPRB,ZTROPO_FAC))

          ! Merge BASCOE rates with selection of corresponding TM5 photodissociation rates
          DO JT = 1,INJ_STRAT_TROP
            ZAJVAL(JL,JK,IBASCOE_DIS(JT))=(1._JPRB-ZTROPO_FAC) * ZRJ(JL,JK,ITM5_DIS(JT)) &
            &                            +(        ZTROPO_FAC) * ZAJVAL(JL,JK,IBASCOE_DIS(JT))
          ENDDO
         ELSEIF (JK > JTROPOP(JL) ) THEN
          ! Take selection of corresponding TM5 photodissociation rates
          DO JT = 1,INJ_STRAT_TROP
            ZAJVAL(JL,JK,IBASCOE_DIS(JT))=ZRJ(JL,JK,ITM5_DIS(JT))
          ENDDO

        ENDIF

      ENDIF

    ENDIF

    IF (LCHEM_JOUT) THEN
        POUT(JL,JK,2) =  ZAJVAL(JL,JK,J_O3_O1D)
        POUT(JL,JK,3) =  ZAJVAL(JL,JK,J_NO2)
        POUT(JL,JK,4) =  ZSZA
    ENDIF

!*  Air density mutiplied with RMD (dry air molar mass) for efficiency
    ZAIRDM(JL) = ZDENS(JL,JK) * RMD
    ZCVM(JL,IAIR)=ZDENS(JL,JK)
    ZCVM(JL,IACID)=0._JPRB
    ZCVM(JL,NCHEM+3)=0._JPRB

!*  convert tracer concentrations from kg/kg to molec/cm3
!*  and assure positivity for initial concentrations
    DO JT=1,NCHEM
      ZCVM0(JL,JT) = MAX(PCEN(JL,JK,JT) / YCHEM(JT)%RMOLMASS * ZAIRDM(JL), 0._JPRB)
    ENDDO

!*  special treatment for SO4, NO3_A and NH4
    IF (LAERCHEM) THEN
!*       use aerosol scheme SO4 for WETCHEM and EQSAM
      ISSO4 = SUM(NTYPAER(1:4)) + 1
      ZCVM0(JL,ISO4) = MAX(PAEROP(JL,JK,ISSO4) / YCHEM(ISO4)%RMOLMASS * ZAIRDM(JL), 0._JPRB)
    ENDIF

    IF (LAERNITRATE) THEN
!*       use aerosol scheme NO3_A for WETCHEM and hetchem. Add both coarse and fine mode into one
      ISNO3 = SUM(NTYPAER(1:5)) + 1
      ZCVM0(JL,INO3_A) =                    MAX(PAEROP(JL,JK,ISNO3)   / YCHEM(INO3_A)%RMOLMASS * ZAIRDM(JL), 0._JPRB)
      ZCVM0(JL,INO3_A) = ZCVM0(JL,INO3_A) + MAX(PAEROP(JL,JK,ISNO3+1) / YCHEM(INO3_A)%RMOLMASS * ZAIRDM(JL), 0._JPRB)
!*       use aerosol scheme NH4 for WETCHEM and hetchem.
      ISNH4 = SUM(NTYPAER(1:6)) + 1
      ZCVM0(JL,INH4) = MAX(PAEROP(JL,JK,ISNH4) / YCHEM(INH4)%RMOLMASS * ZAIRDM(JL), 0._JPRB)
    ENDIF

!*     ZCVM not strictly needed for KPP, but necessary to initialize tropospheric concentrations
!*     in case tropospheric chemistry is switched off
    ZCVM(JL,1:NCHEM) = ZCVM0(JL,1:NCHEM)

!*   Overwrite ECMWF H2O in troposphere as well as tropopause region
    IF (JK > JTROPOP(JL)-IRANGE_TROPOP) THEN
!* Overwrite only in tropopause region, where Strat chemistry is solved:
!      IF (JK < JTROPOP(JL) + 4) THEN
      ZCVM0(JL,IH2O)=PQP(JL,JK)/ YCHEM(IH2O)%RMOLMASS * ZAIRDM(JL)
!      ENDIF
      ZCVM(JL,IH2O)=ZCVM0(JL,IH2O)
    ENDIF

    ! Baseline: We are in Stratosphere
    LLSTRATAIR=.TRUE.

    IF ( LCHEM_TROPO ) THEN

      ! Selection criterion based on humidity (LCHEM_TROPO=TRUE)
      IF ( JK >= KLEVTROP(JL) ) THEN
        ! Now we are in troposphere
        LLSTRATAIR=.FALSE.
      ENDIF

    ELSE
      ! Selection criterion based on chemical composition

      IF (PRSF1(JL,JK)> 4000 .AND.&
        &   PCEN(JL,JK,IO3)/ YCHEM(IO3)%RMOLMASS *RMD < 200E-9 .AND.&
        &   PCEN(JL,JK,ICO)/ YCHEM(ICO)%RMOLMASS *RMD >  40E-9) THEN
          ! Such parcels would suggest tropospheric air, Starting from pressure levels higher than 40hPa
          LLSTRATAIR=.FALSE.
      ENDIF
      ! Make sure that no stratospheric air parcels fall below lower boundary
      IF (PRSF1(JL,JK) > PRSF1(JL,JTROPOP(JL))) THEN
        ! Here definitely no stratospheric air should be present...
        LLSTRATAIR=.FALSE.
      ENDIF

    ENDIF


    ! Try Only solve stratospheric chemistry from tropopause onwards?!
    IF (LLSTRATAIR) THEN

      IF (LCHEM_JOUT) THEN
        ! Fill index 1 of POUT field 2 with the lowest pressure level
        ! where stratospheric chem is still active.
        POUT(JL,1,2) = PRSF1(JL,JK)
      ENDIF

      !VH IF (JK < JTROPOP(JL)+4) THEN

      !* Fixed concentrations
      ZFIX_BASCOE(1) = 0.209*ZDENS(JL,JK)    ! O2 number density
      ZFIX_BASCOE(2) = 0.781*ZDENS(JL,JK)    ! N2 number density


      ! ----------------------------------------------------------------------
      !  Compute strato heterogeneous reaction rates (ZRHET)
      ! ----------------------------------------------------------------------
      IF (YDCHEM%LCHEM_BASCOE_HETCHEM) THEN

        CALL BASCOE_HETCONST(YGFL,JHET_TRACER,PTP(JL,JK),PRSF1(JL,JK),ZDENS(JL,JK),&
         & LL_PSC_POSSIBLE(JL),JTOP_PSC(JL),JBOT_PSC(JL),JK,ZCVM0(JL,1:NCHEM), &
         & ZSA_SIZEDIST(JL,JK,1:NBINS),ZAER(JL,JK,1:NAER),PTSTEP,ZRHET)

      ELSE

        ! in case of emergency switch off:
        ZRHET(:) = 0._JPRB

      ENDIF

      ! Initialize concentrations to KPP (fill VAR)
      CALL BASCOE_KPP_INITIALIZE(YGFL,ZCVM0(JL,1:NCHEM),ZVAR_BASCOE(1:NVAR_BASCOE))

      ! Initialize all kpp reaction rates... (fill RCONST)
      CALL BASCOE_KPP_RATES(ZAJVAL(JL,JK,:),ZRHET, PTP(JL,JK), ZDENS(JL,JK), &
      & ZRCONST_BASCOE, ZVAR_BASCOE(1:NVAR_BASCOE))

      ! starting value for integration time step
      ZCNTRL_P=ZCNTRL
      !VH Maybe to be switched on...    ZCNTRL_P(3) = ZHSAVE_KPP(JL,JLEV)

      !VH - Overwrite initial timestep  - not needed here
      !ZCNTRL_P(3) = MIN(PTSTEP, ZHSTART)


      ! ----------------------------------------------------------------------
      !  Now call the chem box solver
      ! ----------------------------------------------------------------------
      ! Call kpp integrator... (provide 'VAR' and 'RCONST' !)
      ZDENS_DP=ZDENS(JL,JK)
      CALL BASCOE_KPP_INTEGRATOR(0._JPRD, PTSTEP, ICNTRL,ZCNTRL_P,&
       & ISTATUS,ZSTATE,IERR, ZVAR_BASCOE, ZFIX_BASCOE, ZRCONST_BASCOE, ZDENS_DP)

      !- Filter error due to bad concentrations
      IF (IERR>0) THEN

        !No errors: update concentrations...
        CALL BASCOE_KPP_UPDATE_CIFS_CONC(YGFL,ZVAR_BASCOE(1:NVAR_BASCOE),ZCVM(JL,1:NCHEM))

      ELSEIF( IERR==-9) THEN
        WRITE(NULERR,'(a)') '     ZVAR_BASCOE below are out of range:'
        DO JT= 1, NVAR_BASCOE
          IF( ZVAR_BASCOE(JT)/ZDENS(JL,JK) > VMR_BAD_LARGE) THEN
            WRITE(NULERR,'(a,2(i5),a,es12.5)') '  vmr-idx ',JT,JK, ' ; reached value: ',  ZVAR_BASCOE(JT)/ZDENS(JL,JK)
          ENDIF
        ENDDO
        WRITE(NULERR,*) '  -> chem integrator skipped'
      ENDIF

    ! Double check convergence of total N-mass in stratosphere
    ! Provided that RTOL is sufficiently small, this check is no longer required.
    !
    ! Compute NOy tendency (molec N cm-3 dt-1)
    ! N2O+O1D-> 2 NO - Double check indices!
    ! ZNOYTEND=2.*ZCVM(JL,IN2O)*ZCVM(JL,IO1D)*ZRCONST_BASCOE(14)*PTSTEP
    ! NO+ N ->N2O - Double check indices!
    ! ZNOYTEND=ZNOYTEND - 2*ZCVM(JL,IN)*ZCVM(JL,INO)*ZRCONST_BASCOE(50)*PTSTEP
    !
    ! CALL BASCOETM5_NOYMASS(YGFL,ZCVM0(JL,1:NCHEM),ZCVM(JL,1:NCHEM),ZNOYTEND)


    !Apply simplified loss term of tropospheric tracers not active in stratosphere
    CALL TM5_STRATOLOSS(YGFL,ZCVM0(JL,1:NCHEM),ZCVM(JL,1:NCHEM),PTSTEP)

    ! Set O3S tracer to be identical to O3 tracer in stratosphere
    ZCVM(JL,IO3S)=ZCVM(JL,IO3)

    ! Posibly to be added in future:
    ! Fill stratospheric CO2 loss rate [s-1] computed from CO2_C tendency
    ! Here ICO2 refers to the chemical CO2 tendency. Loss= negative.
    ! For now only allow CO2 loss, i.e. no source.
    ! units: [s-1]
    ! PCHEM2GHG(JL,JK,3) = MIN( (ZCVM(JL,ICO2) - ZCVM0(JL,ICO2)) / (ZCVM(JL,ICO2) * PTSTEP) , 0._JPRB)

   ELSE
      ! ----------------------------------------------------------------------
      !  Troposphere
      !
      ! First setup photolysis rates
      DO JT=1,NPHOTO
        ZRJ_IN(JT)=ZRJ(JL,JK,JT)
      ENDDO

      IF (JK < JTROPOP(JL)+JPNLBLEND) THEN
        ! Merge TM5 rates with selection of corresponding BASCOE photodissociation rates
        ! Note that these BASCOE rates were already merged with the TM5 photolysis rates!

        ! Only merge rates when SZA<94
        IF ( ZSZA < 94._JPRB) THEN
          DO JT = 1,INJ_STRAT_TROP
            ZRJ_IN(ITM5_DIS(JT))=ZAJVAL(JL,JK,IBASCOE_DIS(JT))
          ENDDO
        ENDIF
      ENDIF

    IF (LCHEM_JOUT) THEN
        POUT(JL,JK,2) =  ZRJ_IN(JO3D)
        POUT(JL,JK,3) =  ZRJ_IN(JNO2)
    ENDIF


      !calculate reaction rates ZRR using pre-calculated, temperature-dependent values
      ZRR(:)=0.0_JPRB
      CALL TM5_CALRATES(YDCHEM,YDCOMPO,YDEAERSNK,YDEAERSRC,YDEAERATM,YGFL, &
        & 1_JPIM, 1_JPIM, 1_JPIM, 1_JPIM, 1_JPIM, JRATE_TRACER, &
        & KAERO, IMODES, PTP(JL,JK), PRSF1(JL,JK) , &
        & PQP(JL,JK),PAP(JL,JK), PIP(JL,JK),PLP(JL,JK), ZRJ_IN(1:NPHOTO), &
        & ZCLOUD_REFF(JL,JK), PCEN(JL,JK,IO3), PCEN(JL,JK,IHO2),&
        & PCEN(JL,JK,INH4),PCEN(JL,JK,INO3_A), PCEN(JL,JK,ISO4),&
        & PWETDIAM(IL_GLOMAP,IK_GLOMAP,1:IMODES), &
        & PWETVOL(IL_GLOMAP,IK_GLOMAP,1:IMODES),PND(IL_GLOMAP,IK_GLOMAP,1:IMODES), &
        & PAEROP(JL,JK,:),ZSAD_AER(JL,JK),ZSAD_CLD(JL,JK), ZSAD_ICE(JL,JK), &
        & ZRR)

        ! Special fix for HO2+HO2 -> H2O2 reactions.
        ! In tm5_calrates.F90 and in tm5_do_ebi.F90 these are treated as single-body reactions,
        ! but in KPP these are treated as second-order reaction.
        ! Now scale kho2_aer and kho2_liq with HO2 concentrations,
        ! to account for it as a self-reaction in KPP solver
        ! To ensure a reasonably small reaction rate assume a minimum HO2 concentration

        ZCONC_HO2 = MAX(ZCVM(JL,IHO2), 1E5_JPRB)
        ZRR(KHO2L)   =ZRR(KHO2L)    / ZCONC_HO2
        ZRR(KHO2_AER)=ZRR(KHO2_AER) / ZCONC_HO2



!3.1  First solve wet sulphur/ammonia chemistry
      CALL TM5_WETCHEM_POINT(YGFL,PTSTEP,PTP(JL,JK),PAP(JL,JK),PRSF1(JL,JK),PLP(JL,JK),&
         & ZHPLUS(JL),ZCVM(JL,1:NCHEM))

      ! starting value for integration time step
      ZCNTRL_P=ZCNTRL
      !VH --- to be switched on...    ZCNTRL_P(3) = ZHSAVE_KPP(JL,JLEV)

      ! Try to initialize timestep with something reasonable...
      !ZCNTRL_P(3) = PTSTEP/20._JPRB

      ! Near surface reduce initial time step near surface (at least lowest 2 levels)
      IF ( JK > KLEV-3 ) THEN
        ZCNTRL_P(3) = ZCNTRL(3)/6._JPRB
      ELSEIF ( PRSF1(JL,JK)> 90000_JPRB .OR. JK > KLEV-6 ) THEN
        ZCNTRL_P(3) = ZCNTRL_P(3)/4._JPRB
      ENDIF

      ! Update kpp rates... (Assume dry deposition is done outside chemistry!!!
      CALL TM5_KPP_RATES(ZRR(1:NREAC),ZRJ_IN(1:NPHOTO),ZRCONST_TM5 )

      ! Initialize concentrations to KPP (merge with prepare_kpp_conc)
      CALL TM5_KPP_INITIALIZE(YGFL,ZCVM(JL,1:NCHEM),ZVAR_TM5(1:NVAR_TM5),LAERSOA_COUPLED)

      ! Call kpp integrator...
      CALL TM5_KPP_INTEGRATOR(0._JPRD, PTSTEP, ICNTRL,ZCNTRL_P, ISTATUS,ZSTATE,IERR,ZVAR_TM5,ZRCONST_TM5,ZDENS(JL,JK))


!- Filter error due to bad concentrations
      IF (IERR>0) THEN

         !No errors: update concentrations...
         CALL TM5_KPP_UPDATE_CIFS_CONC(YGFL,ZVAR_TM5(1:NVAR_TM5),ZCVM(JL,1:NCHEM),LAERSOA_COUPLED)

      ELSEIF( IERR==-9) THEN
        WRITE(NULERR,'(a)') '     ZVAR_TM5 below are out of range:'
        DO JT= 1, NVAR_TM5
          IF( ZVAR_TM5(JT)/ZDENS(JL,JK) > VMR_BAD_LARGE) THEN
            WRITE(NULERR,'(a,2(i5),a,es12.5)') '  vmr-idx ',JT,JK, ' ; reached value: ', ZVAR_TM5(JT)/ZDENS(JL,JK)
          ENDIF
        ENDDO
        DO JT=1,NCHEM
          IF( ZCVM(JL,JT)/ZDENS(JL,JK) > VMR_BAD_LARGE/1E3) THEN
            WRITE(NULERR,'(a,2(i5),a,es12.5)') ' Input Tracer ',JT,JK, ' ; reached value: ', ZCVM(JL,JT)/ZDENS(JL,JK)
          ENDIF
        ENDDO
        WRITE(NULERR,*) '  -> chem integrator skipped'
      ENDIF


      IF (LCHEM_DIAC) THEN
!3.4           increase budget accumulators ZCR2 and ZCR3 (photolysis / gas-phase chem) at indec JL
        CALL TM5_IBUD(YGFL,1,1,1,ZRR,ZRJ_IN,ZCVM(JL,1:NCHEM+3),ZCR2(JL,1:NPHOTO),ZCR3(JL,1:NREAC))
      ENDIF

      ! update new preferred time step...
      ! ZHSAVE_KPP(JL,JLEV) = ZSTATE(NHNEW)


      !4.0 Eqsam solver for aerosol- gas phase interaction
      ! done in aerosols in LAERNITRATE
      IF (.NOT. LAERNITRATE) THEN

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

      ENDIF ! NOT LAERNITRATE

      IF (LAERSOA .AND. LAERSOA_COUPLED) THEN
        ! apply gas-aerosol partitioning for Sec Org Aerosol, and gas-phase version
        INBDU=NTYPAER(1)+NTYPAER(2)
        INBOM=NTYPAER(1)+NTYPAER(2)+NTYPAER(3)

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
           ZJSOA=4E-4_JPRB*ZRJ_IN(JNO2)
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
   ENDIF ! LAERSOA

   ! Compute loss of O3S in troposphere at location JL
   CALL TM5_O3S(YDCHEM,NCHEM,PTSTEP,ZRR,ZRJ_IN,ZCVM(JL,1:NCHEM) )

   ! Fill tropospheric CO2 production tendency computed from CO + OH reaction budget
   ! units: kg CO2 / kg / sec
   ZFAC=RMCO2/(ZAIRDM(JL)*PTSTEP)
   PCHEM2GHG(JL,JK,2) = ZCVM(JL,ICO)*ZCVM(JL,IOH)*ZRR(KCOOH)  * ZFAC

  ENDIF ! Within stratosphere or troposphere

!VH make sure that within troposphere and around tropopause
!VH we get the tendencies to arrive at ECMWF-H2O
!VH so impose initial and final concentration fields.
!    IF (JK > JTROPOP(JL)-8) THEN
!Test for lower altitude wrt H2O fix
    IF (JK > JTROPOP(JL)-IRANGE_TROPOP) THEN
      ZCVM0(JL,IH2O)=PCEN(JL,JK,IH2O) / YCHEM(IH2O)%RMOLMASS *ZAIRDM(JL)

!* Nudge towards ECMWF-H2O, but don't enforce it...
      ZH2O_ECMWF=PQP(JL,JK)/ YCHEM(IH2O)%RMOLMASS * ZAIRDM(JL)
      IF (JK < JTROPOP(JL)+IRANGE_TROPOP) THEN
        ! 3 hour decay time near tropopause
        ZCVM(JL,IH2O) = ZCVM0(JL,IH2O)+(1._JPRB-EXP(-PTSTEP/(10800._JPRB)))*(ZH2O_ECMWF-ZCVM0(JL,IH2O))
      ELSE
        ! 1 day decay time in toposphere
        ZCVM(JL,IH2O) = ZCVM0(JL,IH2O)+(1._JPRB-EXP(-PTSTEP/(86400._JPRB)))*(ZH2O_ECMWF-ZCVM0(JL,IH2O))
      ENDIF
!      IF (JK < JTROPOP(JL)+4 ) THEN
!        ZCVM(JL,IH2O)=PQP(JL,JK)/ YCHEM(IH2O)%RMOLMASS * ZAIRDM
!*Overwrite only in tropopause region, where Strat chemistry is solved:
!        ZCVM(JL,IH2O)=PQP(JL,JK)/ YCHEM(IH2O)%RMOLMASS * ZAIRDM
!      ELSE
!* Elsewhere (in troposphere) let H2O decay with 2-day decay rate...
!        ZCVM(JL,IH2O)=ZCVM0(JL,IH2O)*EXP( -5.787E-6_JPRB*PTSTEP )
!      ENDIF
    ENDIF




  ENDDO ! loop over JL



  ! ----------------------------------------------------------------------
  !4.0 Compute sedimentation in stratosphere... (Here? Or before chemistry?)
  !
  CALL BASCOE_PSC_PARAM(YGFL,KIDIA,KFDIA,KLON,JPSC_TRACER,PTSTEP,JTROPOP,JK,PTP(1:KLON,JK), &
    & PRSF1(1:KLON,JK),ZCVM(1:KLON,1:NCHEM))



  ! ----------------------------------------------------------------------
  !5.0 convert concentration tendencies to mass mixing ratio
  !
  DO JL=KIDIA,KFDIA
!*  Air density mutiplied with RMD (dry air molar mass) for efficiency
    ZAIRDM(JL) = ZDENS(JL,JK) * RMD
    DO JT=1,NCHEM
      IF (JT /= ISTRATAER) THEN
        PTENC1(JL,JK,JT) =   (ZCVM(JL,JT)-ZCVM0(JL,JT)) * YCHEM(JT)%RMOLMASS /(ZAIRDM(JL)*PTSTEP)
      ENDIF
    ENDDO
  ENDDO

  ! ----------------------------------------------------------------------
  !6.0 add budgets for this timestep
  !
  IF (LCHEM_DIAC) THEN
    CALL TM5_RBUD(YGFL,KIDIA,KFDIA,KLON,KLEV,JK,IOH,1,ZAIRDM,ZCR2,ZCR3,PBUDJ,PBUDR,PBUDX)
  ENDIF

ENDDO ! loop over levels



! ----------------------------------------------------------------------
! 7.0 Simple boundary conditions at surface for a list of BASCOE species -
! Check this when introducing emissions!
!
! YC 20180830: allow monthly varying, latitude band dependent LBC
!
! get the index of current month
!   assume current month in LBC coordinates, otherwise nearest neighbour
!   *this could be refined, as trend and seasonality is lost*
IMONTH_LBC = MINLOC( ABS(MONTH_LBC - IYYYYMM) )
IMONTH_LBC_H=IMONTH_LBC(1)

DO JB=1,NBC
  JBC    = BASCOE_BC(JB)
  ! Select appropriate BC and convert to mass ratio [kg/kg]
  !
  ZBCVAL(1:NLATBOUND_LBC-1) = VALUES_LBC(JB,IMONTH_LBC_H,1:NLATBOUND_LBC-1)* YCHEM(JBC)%RMOLMASS / RMD

  ZTAU_NUDGE=2500._JPRB
  IF (ANY(YEMIS2D_DESC(1:NEMIS2D_DESC)%SPECIES == YCHEM(JBC)%CNAME)) THEN
    ! Apply less stringent relaxation in case of using emissions...
    ZTAU_NUDGE=86400._JPRB
  ENDIF

  DO JL=KIDIA,KFDIA
    DO JLAT_LBC = 1, NLATBOUND_LBC-1
      IF( ZLAT(JL) >= XLATBOUND_LBC(JLAT_LBC) .AND. ZLAT(JL) <= XLATBOUND_LBC(JLAT_LBC+1) ) EXIT
    ENDDO
    ZTENBC(JL) = ZBCVAL(JLAT_LBC) - PCEN(JL,KLEV,JBC)
    PTENC1(JL,KLEV,JBC) = ((1._JPRB - EXP(-PTSTEP / (ZTAU_NUDGE)))*ZTENBC(JL))  / PTSTEP
    ! store this special LBC budget contribution in POUT(:,JB+1,1)
    ! Preserve first index for CH4 boundary condition, see below!
    POUT(JL,JB+1,1)=PTSTEP *PTENC1(JL,KLEV,JBC)*PDELP(JL,KLEV) / RG
  ENDDO
ENDDO



! IF (LGHG_CHEMTEND_CH4) THEN
  ! Fill field with CH4 atmospheric loss rate [s-1] , loss term is positive.
  ! Exclude contribution of surface nudging (if any), as computed next
  DO JLEV=1,KLEV
    DO JL=KIDIA,KFDIA
      PCHEM2GHG(JL,JLEV,1) = - PTENC1(JL,JLEV,ICH4) / PCEN(JL,JLEV,ICH4)
    ENDDO
  ENDDO
! ENDIF

!
! 8.0 replace CH4 at surface only with value from boundary condition.
IF (.NOT. LCHEM_ANACH4) THEN
  IYEAR_CH4=IYEAR
  IF (KCHEM_YEARPI < 2000)  IYEAR_CH4=KCHEM_YEARPI

  CALL TM5_BOUNDARY_CH4(YGFL,KIDIA,KFDIA,KLON,IMONTH,IYEAR_CH4,ICH4,PGELAT,PCEN(:,KLEV,ICH4),ZTENBC)

  IF (LCHEM_WEAK_CH4_RELAXATION) THEN
    ZTAU_NUDGE=86400._JPRB
  ELSE
    ZTAU_NUDGE=2500._JPRB
  ENDIF

  DO JL=KIDIA,KFDIA
    ! Introduce more specific constraint on CH4 for background conditions
    ! This should be useable when including CH4 emissions and dry deposition.

    PTENC1(JL,KLEV,ICH4) =  ((1._JPRB - EXP(-PTSTEP / (ZTAU_NUDGE)))*ZTENBC(JL))  / PTSTEP
    ! store this special LBC CH4 budget contribution  POUT(:,1,1)
    POUT(JL,1,1)=PTSTEP *PTENC1(JL,KLEV,ICH4)*PDELP(JL,KLEV) / RG

    ! Simple, direct tendency
    ! PTENC1(JL,KLEV,ICH4) =  ZTENBC(JL)  / PTSTEP
    ! store this special LBC CH4 budget contribution  POUT(:,1,1)
    ! POUT(JL,1,1)=PTSTEP *PTENC1(JL,KLEV,ICH4)*PDELP(JL,KLEV) / RG
  ENDDO


ENDIF

IF (LCHEM_JOUT) THEN
  ! Fill field with CH4 tendency term for analysis purposes
  DO JLEV=1,KLEV
    DO JL=KIDIA,KFDIA
      POUT(JL,JLEV,5) = PTENC1(JL,JLEV,ICH4)
    ENDDO
  ENDDO
ENDIF


END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHEM_BASCOETM5',1,ZHOOK_HANDLE)
END SUBROUTINE CHEM_BASCOETM5
