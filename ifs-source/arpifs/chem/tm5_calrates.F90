! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_CALRATES(YDCHEM,YDCOMPO,YDEAERSNK,YDEAERSRC,YDEAERATM,YGFL,KIDIA, KFDIA, KLON, KL, KLEV, KTRACER, &
     & KAEROC, KMODES, PTP, PRSF1, &
     & PQP, PAP, PIP, PLP, PRJ, PCLOUD_REFF, PO3, PHO2, &
     & PNH4,PNO3_A, PSO4, PWETDIAM, PWETVOL,PND, &
     & PAEROP,PSAD_AER,PSAD_CLD,PSAD_ICE, &
     & PRR)


!**   DESCRIPTION
!     ----------
!
!   Part of TM5 routines for IFS chemistry:
!   evaluation of reaction rates (and photolysis rates)
!
!
!
!**   INTERFACE.
!     ----------
!          *TM5_calrates* IS CALLED FROM *CHEM_tm5*.

! INPUTS:
! -------
! KIDIA :  Start of Array
! KFDIA :  End  of Array
! KLON  :  Length of Arrays
! KL    :  Current level
! KLEV  :  Number of Levels
! KAEROC : dimension of aerosol field
! KMODES : number of GLOMAP modes in arrays
! PQP     (KLON,KLEV)         :  SPECIFIC HUMIDITY            (kg/kg)
! PTP     (KLON,KLEV)         :  TEMPERATURE                  (K)
! PAP     (KLON,KLEV)         :  CLOUD FRACTION               0..1
! PIP     (KLON,KLEV)         :  ICWC                         (kg/kg)
! PLP     (KLON,KLEV)         :  LCWC                         (kg/kg)
! PRSF1(KLON,KLEV)             : FULL-LEVEL PRESSURE          (Pa)
! PO3, PHO2, PNH4,PNO3_A, PSO4: concentrations of ozone, HO2, and ammonium, nitrate and sulphate chem. aerosol

! PRJ  (KLON,NPHOTO)     : photolysis rates
! PAEROP(KLON,KLEV,KAEROC)  : Aerosol concentrations  (kg/kg)
!
!
! OUTPUTS:
! -------
! PRR  (KLON,NREAC)      : reaction rates
! PRJ  (KLON,NPHOTO)     : photolysis rates
! PSAD_AER,PSAD_CLD,PSAD_ICE: surface area density for aerosol, cloud and ice, for diagnostics purposes
!
! LOCAL:
! -------
!
!
!     AUTHOR.
!     -------
!        VINCENT HUIJNEN    *KNMI*
!        TM5-community
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2009-09-08



USE PARKIND1 , ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMCST   , ONLY : RPI  , RMD , RD, RNAVO
USE TM5_CHEM_MODULE , ONLY : &
 & KNOO3,      KHO2NO,      KMO2NO,      KNO2OH,      KOHHNO3,    KNO2O3,     KNONO3,      KNO2NO3,    &
 & KHNO4A,     KHNO4B,      KN2O5A,      KN2O5B,      &
 & KN2O5,      KHNO4OH,     KNO2HO2,     KHNO4M,      KMENO2,     KMENO2M,    KODM,        KH2OOD,     &
 & KO3HO2,     KCOOH,       KO3OH,       KHPOH,       KFRMOH,     KCH4OH,     KOHMPER,     KOHROOH,    &
 & KMO2MO2,    KHO2OH,      KHO2HO2,     KN2O5AQ,     KN2O5L,     KH2OH,      KC41,        KC43,       &
 & KC44,       KC46,        KC47,        KC48,        KC49,       KC50A,      KC50B,       KC52,       &
 & KC53,       KC54,        KC57,        KC58,        KC59,       KC61,       KC62,        KC73,       &
 & KC76,       KC77,        KC78,        KC79,        KC80,       KC81,       KC82,        KC83,       &
 & KC84,       KC85,        KDMSOHA,     KDMSOHB,     KDMSNO3,    KSO2OH,     KNH3SO4,     KNH3OH,     &
 & KNH2NO,     KNH2NO2,     KNH2HO2,     KNH2OH,      KNH2O2,     KRN222,     KO3PO3,      KO3PO2,     &
 & KNO2OHA,    KNO2OHB,     KOHHNO3A,    KOHHNO3B,    KOHHNO3C,   KNO2NO3A,   KNO2NO3B,    KNO2HO2A,   &
 & KNO2HO2B,   KHO2HO2A,    KHO2HO2B,    KHO2HO2C,    KC47A,      KC47B,      KC61A,      &
 & KC61B,      KSO2OHA,     KSO2OHB,     KDMSOHC,     KCOOHA,     KCOOHB,     KCOOHC,      KCOOHD,     &
 & KHONO,      KHONOA,      KHONOB,      KOHHONO,     KMO2HO2A,   KMO2HO2B,   KMENO2A,     KMENO2B,    &
 & KOHMCHO,    KOHMCH2CHO,  KNO3MCHO,    KNO3MCH2CHO, KOHMVK,     KOHOLE,     KOHMACR,     KO3MVK,     &
 & KO3OLE,     KO3MACR,     KNO3OLE,     KOHC3H8,     KO3C3H6,    KNO3C3H6,   KOHETHOH,    KOHTERP,    &
 & KO3TERP,    KNO3TERP,    KOHHCOOH,    KOHCH3OH,    KNO3HO2,    KOHC3H6A,   KOHC3H6B,    KOHC3H6,    &
 & KNO3MO2,    KOHC2H6,     KNO3C2O3,    KNO3XO2,     KOHMCOOH,   KOHISPD,    KO3ISPD,     KNO3ISPD,   &
 & KOHACET,    KOHACETA,    KOHACETB,    KACO2HO2,    KACO2MO2,   KACO2NO,    KACO2XO2,    KXO2XO2N,   &
 & KXO2N,      KOHHCN,      KOHCH3CN,    KOHTOL,      KO3TOL,     KNO3TOL,    KOHXYL,      KO3XYL,     &
 & KNO3XYL,    KAROO2NO,    K2AROO2,     KAROO2HO2,   KAROO2XO2,  KCXYLO3A,   KCXYLO3B,    KNOIC3H7O2, &
 & KHO2IC3H7O2,KNOHYPROPO2, KHO2HYPROPO2,KISOPBO2A,   KISOPBO2B,  KISOPDO2A,  KISOPDO2B,   KISOPBO2HO2,&
 & KISOPBO2NO, KISOPDO2HO2, KISOPDO2NO,  KISOPOOHOH,  KHPALD1OH,  KHPALD2OH,  KGLYOH,      KGLYALDOH,  &
 & KHYACOH,    &
 & KHO2L, KHO2_AER , KN2O5_AER, KNO3_AER, &
 & KSO2O3G, KSO3H2O,KOCSOH, KSOA_AG, &
 & RATES_LUT, NTLOW, NTEMP, &
 & XMH2O, XMAIR, NREAC
USE TM5_PHOTOLYSIS, ONLY : NPHOTO, JO2, JO3D

USE YOEAERATM, ONLY : TEAERATM, YREAERATM
USE YOEAERSRC, ONLY : TEAERSRC
USE YOEAERSNK ,ONLY : TEAERSNK
USE YOMCHEM  , ONLY : TCHEM
USE YOMCOMPO , ONLY : TCOMPO


USE UKCA_MODE_SETUP, ONLY: NMODES, SIGMAG, MODE

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(TCHEM),       INTENT(INOUT) :: YDCHEM
TYPE(TCOMPO),       INTENT(IN)   :: YDCOMPO
TYPE(TEAERSNK),    INTENT(INOUT) :: YDEAERSNK
TYPE(TEAERSRC),    INTENT(INOUT) :: YDEAERSRC
TYPE(TEAERATM),    INTENT(INOUT) :: YDEAERATM
TYPE(TYPE_GFLD),   INTENT(INOUT) :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA, KLON, KL, KLEV, KAEROC, KMODES
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRACER(2)
REAL(KIND=JPRB),   INTENT(IN)    :: PRSF1(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN)    :: PTP(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN)    :: PQP(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN)    :: PAP(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN)    :: PIP(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN)    :: PLP(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN)    :: PCLOUD_REFF(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN)    :: PO3(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PHO2(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PNH4(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PNO3_A(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PSO4(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PWETDIAM(KLON,KLEV,KMODES)
REAL(KIND=JPRB),   INTENT(IN)    :: PWETVOL(KLON,KLEV,KMODES)
REAL(KIND=JPRB),   INTENT(IN)    :: PND(KLON,KLEV,KMODES)

REAL(KIND=JPRB),   INTENT(INOUT) :: PRJ(KLON,NPHOTO)
REAL(KIND=JPRB),   INTENT(IN)    :: PAEROP(KLON,KLEV,KAEROC)
REAL(KIND=JPRB),   INTENT(INOUT) :: PSAD_AER(KLON,KLEV), PSAD_CLD(KLON,KLEV), PSAD_ICE(KLON,KLEV)

REAL(KIND=JPRB),   INTENT(OUT)   :: PRR(KLON,NREAC)

! * LOCAL

INTEGER(KIND=JPIM)            :: JAER,ITYP,IBIN
INTEGER(KIND=JPIM)            :: IO3, IHO2
INTEGER(KIND=JPIM)            :: IRH, JTAB

! * cloud and aerosol chemistry of n2o5, ho2 and no3
REAL(KIND=JPRB)               :: ZKT_LIQ_N2O5, ZKT_ICE_N2O5, ZKT_AER_N2O5, ZM_AER, ZRHO_P, ZVFRAC
REAL(KIND=JPRB)               :: ZKT_LIQ_HO2,ZKT_ICE_HO2, ZKT_AER_HO2, ZRHO_AER
REAL(KIND=JPRB)               :: ZRH,ZWV, ZTR, ZR_EFF, ZKTNO3_AER
REAL(KIND=JPRB), PARAMETER    :: ZRHO_SS=2.2 ! gr/cm^3 ! mean density sea salt. Is this correct?
REAL(KIND=JPRB), PARAMETER    :: ZRHO_DD=2.5 ! gr/cm^3 ! mean density dessert dust. Is this correct?
REAL(KIND=JPRB), PARAMETER    :: ZRHO_OM=1.8  ! gr/cm^3
REAL(KIND=JPRB), PARAMETER    :: ZRHO_BC=1.0  ! gr/cm^3
REAL(KIND=JPRB), PARAMETER    :: ZRHO_SO4=1.7  ! gr/cm^3
REAL(KIND=JPRB), PARAMETER    :: ZRHO_H2O=1.0  ! gr/cm^3 (water)
REAL(KIND=JPRB), PARAMETER    :: ZRHO_SOA=1.8_JPRB  ! gr/m^3

!* tentative sea salt effective radii (cm).
REAL(KIND=JPRB), DIMENSION(3), PARAMETER    :: ZR_SS=(/0.1E-4,0.7E-4,8.E-4/)
!* Taken from J.J.Morcrette
REAL(KIND=JPRB), DIMENSION(3), PARAMETER    :: ZR_DD=(/0.051E-4,0.81E-4,18.E-4/)
REAL(KIND=JPRB), DIMENSION(3),PARAMETER     :: ZR_SOA=(/0.03E-4_JPRB,0.15E-4_JPRB,0.15E-4_JPRB/)
REAL(KIND=JPRB), PARAMETER    :: ZR_OM = 0.13E-4 !cm
REAL(KIND=JPRB), PARAMETER    :: ZR_BC = 0.04E-4 !cm
!REAL(KIND=JPRB), PARAMETER    :: ZR_SO4 = 0.355e-4 !cm SO4 dry particle radius, email Morcrette
REAL(KIND=JPRB), PARAMETER    :: ZR_SO4 = 0.18E-4   !cm SO4 dry particle radius, Martin et al., 2003

REAL(KIND=JPRB), PARAMETER    :: ZR_NO3_A=0.2E-4 ! cm (guess value for nitrate particles. Assume also for ammonia)

!* Growth factors corresponding to RH table as given in RRHTAB, according to J.J. Morcrette
REAL(KIND=JPRB), DIMENSION(12), PARAMETER    :: ZRH_GROWTH_SO4= &
         & (/1.00,1.00,1.00,1.00,1.169,1.220,1.282,1.363,1.485,1.581,1.732,2.085/)
! According to Chin et al., AMS 2002
REAL(KIND=JPRB), DIMENSION(12), PARAMETER    :: ZRH_GROWTH_BC= &
         & (/1.00,1.00,1.00,1.00,1.00 ,1.000,1.000,1.000,1.200,1.300,1.400,1.500/)
REAL(KIND=JPRB), DIMENSION(12), PARAMETER    :: ZRH_GROWTH_OM= &
         & (/1.00,1.00,1.00,1.00,1.169,1.200,1.300,1.400,1.500,1.550,1.600,1.800/)
REAL(KIND=JPRB), DIMENSION(12), PARAMETER    :: ZRH_GROWTH_SS= &
         & (/1.00,1.00,1.00,1.00,1.200,1.600,1.700,1.800,2.000,2.200,2.400,2.900/)
! SOA Martin et al 2003, Table 1
REAL(KIND=JPRB), DIMENSION(12), PARAMETER    :: ZRH_GROWTH_SOA= &
         & (/1.00_JPRB,1.00_JPRB,1.00_JPRB,1.00_JPRB,1.00_JPRB,1.2_JPRB,&
         & 1.3_JPRB,1.4_JPRB,1.5_JPRB,1.6_JPRB,1.7_JPRB,1.9_JPRB/)


! Assumed uptake coefficients.
REAL(KIND=JPRB)               :: ZG_HO2_L, ZG_HO2_AER ! gamma (HO2) on cloud according to Thornton
REAL(KIND=JPRB), PARAMETER    :: ZG_HO2_I=0.025_JPRB ! cooper & abbatt, 1996.
!REAL(KIND=JPRB), PARAMETER    :: ZG_HO2_L=0.06_JPRB  ! Bulk liquid uptake for HO2, Kolb et al., ACP 2010
REAL(KIND=JPRB), PARAMETER    :: ZG_HO2_DD=0.06_JPRB ! Bedjanian et al., ACPD 2013 consider this as an upper limit...
! REAL(KIND=JPRB), PARAMETER    :: ZG_HO2 = 0.2       ! Jacob, 2000 (but value not used)
! REAL(KIND=JPRB), PARAMETER    :: ZG_HO2_HIGH = 0.7   !Liang et al., ACPD 2013; valid for 'anthropogenic' aerosol
! Specifics related to HO2 heterogeneous chemistry (Thornton et al., JGR 2008)
REAL(KIND=JPRB):: ZDG, ZM_CLD,ZM_ICE,ZSAD_CLOUD,ZSAD_ICE,ZR_ICE,ZRHO
REAL(KIND=JPRD):: ZCONC_HO2, ZNUM, ZK_EFF
REAL(KIND=JPRD):: ZDENOM
REAL(KIND=JPRB):: ZHPL_FIX=3.16227E-6_JPRB ! is rain water ph=5.5 H+
REAL(KIND=JPRB), PARAMETER :: ZK_EQ  = 2.1E-5_JPRB
REAL(KIND=JPRB), PARAMETER :: ZHNRYA=4E3_JPRB ! Hanson et al., 1992
REAL(KIND=JPRB), PARAMETER :: ZHNRYB=5900_JPRB ! Hanson et al., 1992
REAL(KIND=JPRB)            :: ZH_HO2 ! Henry constant for  HO2
REAL(KIND=JPRB)            :: ZH_EFF      ! Effective HO2 henry constant
REAL(KIND=JPRB), PARAMETER :: ZK1A=2.4E9, ZK1B=-2.36E3    ! Reaction rate constants
REAL(KIND=JPRB), PARAMETER :: ZK2A=1.6E10, ZK2B=-1.51E3   ! for aq. phase HO2 ractions
REAL(KIND=JPRB)            :: ZK1, ZK2
!REAL(KIND=JPRB), PARAMETER :: ZALPH_HO2I  = 1./ 0.5       ! Inverse of Alpha
REAL(KIND=JPRB), PARAMETER :: ZALPH_HO2I  = 1._JPRB/ 0.05_JPRB   ! Inverse of Alpha. Reduced value according to Whalley

!VH Standard settings for g_n2o5
REAL(KIND=JPRB), PARAMETER    :: ZG_N2O5_I=0.02_JPRB !, g_n2o5_l=0.02_JPRB
REAL(KIND=JPRB), PARAMETER    :: ZG_N2O5_SS=0.02
REAL(KIND=JPRB), PARAMETER    :: ZG_N2O5_DD=0.01 ! Tang et al., ACP 2010
REAL(KIND=JPRB), PARAMETER    :: ZG_N2O5_OM=0.02
REAL(KIND=JPRB), PARAMETER    :: ZG_N2O5_BC=0.01 ! Tang et al., ACP 2010
REAL(KIND=JPRB), PARAMETER    :: ZG_N2O5_SO4=0.02 ! e.g. Macintyre and Evans, ACP 2010
REAL(KIND=JPRB), PARAMETER    :: ZG_N2O5_NO3A=0.002 ! Factor 10 lower than SO4, see Zhang, Jacob et al., ACP 2012
REAL(KIND=JPRB)               :: ZG_N2O5_L, ZG_N2O5
REAL(KIND=JPRB), PARAMETER    :: ZM_HO2 =  (1.+2*16.)/6.02E23_JPRB *1E-3_JPRB ! HO2 mass, kg/molec
REAL(KIND=JPRB), PARAMETER    :: ZM_N2O5 =  (2*14+5*16.)/6.02E23_JPRB *1E-3_JPRB ! N2O5 mass, kg/molec
REAL(KIND=JPRD)               :: ZC_THERMAL_N2O5, ZC_THERMAL_HO2
! NO3 uptake...
REAL(KIND=JPRB)               :: ZG_NO3 ! Type=dependent

! Coefficients when using GLOMAP aerosol, specified per mode
! MODE_NAMES=(/'nuc_sol','ait_sol','acc_sol','cor_sol','ait_ins','acc_ins','cor_ins'/)
REAL(KIND=JPRB),DIMENSION(7), PARAMETER    :: ZG_N2O5_MODE= &
  &         (/0.02     ,  0.02    , 0.02    , 0.02    , 0.01    , 0.01    , 0.01     /)
REAL(KIND=JPRB),DIMENSION(7), PARAMETER    :: ZG_NO3_MODE= &
  &         (/0.01     ,  0.01    , 0.01    , 0.01    , 0.00    , 0.00    , 0.00     /)

REAL(KIND=JPRB)               :: ZR_CLOUD
REAL(KIND=JPRB)               :: ZDIAMW, ZSAD_AER,ZSAD_AERT

! * Reduced efficiency of heterogeneous chemistry on clouds due to subgrid scale effects
REAL(KIND=JPRB)               :: ZCC, ZSGS_MIX_HO2,ZSGS_MIX_N2O5
! Ratio of ambient HO2 concenrtations outside / inside cloud
!REAL(KIND=JPRB), PARAMETER    :: ZHO2_Ratio = 10.0


! * help variables
INTEGER (KIND=JPIM)  :: ITEMP
REAL(KIND=JPRB)      :: ZTEMP, ZRX1, ZRX2, ZRX3, ZNN, ZNN_B
REAL(KIND=JPRB)      :: ZH2OX, ZAIRD, ZO2, ZO3, ZO3_P


REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! * counters
INTEGER(KIND=JPIM) :: JL,JMODE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#include "abor1.intfb.h"


IF (LHOOK) CALL DR_HOOK('TM5_CALRATES',0,ZHOOK_HANDLE )
ASSOCIATE(NACTAERO=>YGFL%NACTAERO, YCHEM=>YGFL%YCHEM, &
 & RRHTAB=>YDEAERSNK%RRHTAB, RSSGROWTH_RHTAB=>YDEAERSNK%RSSGROWTH_RHTAB, &
 & RSSDENS_RHTAB=>YDEAERSNK%RSSDENS_RHTAB, &
 & NTYPAER=>YDEAERATM%NTYPAER, &
 & YAERO_DESC=>YREAERATM%YAERO_DESC, &
 & RSS_DRY_DIAFAC=>YDEAERATM%RSS_DRY_DIAFAC, &
 & RSS_DRY_MASSFAC=>YDEAERATM%RSS_DRY_MASSFAC,AERO_SCHEME=>YDCOMPO%AERO_SCHEME )

IO3 =KTRACER(1)
IHO2=KTRACER(2)

! Introduce check on NMODES
IF (TRIM(AERO_SCHEME)=="glomap") THEN
  CALL ABOR1("OIFS - GLOMAP should never be called, EXIT")

ENDIF

ZRX3=0.6_JPRB
ZNN=0.75_JPRB-(1.27_JPRB*LOG10(0.41_JPRB))
ZNN_B=0.75_JPRB-(1.27_JPRB*LOG10(0.35_JPRB))

DO JL=KIDIA,KFDIA
    ZTEMP=PTP(JL,KL)
    ITEMP=NINT(ZTEMP-FLOAT(NTLOW))
    ITEMP=MIN(MAX(ITEMP,1_JPIM),NTEMP) !limit temperatures in loop-up table range

    !
    ! 1.0 Calculation of relative humidity [%]
    ! Richardson's approximation for water vapor pressure
    ! should be calculated in subroutine rates
    !
    ZAIRD = 7.24291E16_JPRB*PRSF1(JL,KL)/ZTEMP
    ZH2OX = PQP(JL,KL)*ZAIRD*XMAIR/XMH2O
    ! ZTR=1.0_JPRB-373.15_JPRB/ZTEMP
    ! ZWV=exp((((-.1299_JPRB*ZTR-0.6445_JPRB)*ZTR-1.976_JPRB)*ZTR+13.3185_JPRB)*ZTR)
    ! ZRH=ZH2OX*ZTEMP/(1013.25_JPRB*ZWV*7.24291e16_JPRB)
    ! ZRH = max(min(ZRH,100.0_JPRB),0.0_JPRB)   !limit rh between 0-100%

    ZO2=0.209476_JPRB*ZAIRD
    ZO3 = MAX(PO3(JL) / YCHEM(IO3)%RMOLMASS *ZAIRD*RMD ,0._JPRB)
    !
    !2.0 **** calculate three body reaction rates
    !
    ! ZRX1=RATES_LUT(KNO2OHA,ITEMP)*ZAIRD
    ! ZRX2=RATES_LUT(KNO2OHB,ITEMP)
    ! PRR(JL,KNO2OH)=ZRX1/(1._JPRB+ZRX1/ZRX2)*ZRX3**(1._JPRB/(1._JPRB+LOG10(ZRX1/ZRX2)**2_JPIM))
    ZRX1=RATES_LUT(KNO2OHA,ITEMP)*ZAIRD
    ZRX2=RATES_LUT(KNO2OHB,ITEMP)
    PRR(JL,KNO2OH)=(ZRX1*ZRX2)/(ZRX1+ZRX2)*10**(LOG10(0.41_JPRB)/(1_JPRB+(LOG10(ZRX1/ZRX2)/ZNN)**2_JPIM))

    ZRX1=RATES_LUT(KOHHNO3C,ITEMP)
    ZRX2=RATES_LUT(KOHHNO3B,ITEMP)*ZAIRD
    PRR(JL,KOHHNO3)=RATES_LUT(KOHHNO3A,ITEMP)+ZRX1*ZRX2/(ZRX1+ZRX2)
    ZRX1=RATES_LUT(KNO2NO3A,ITEMP)*ZAIRD
    ZRX2=RATES_LUT(KNO2NO3B,ITEMP)
    PRR(JL,KNO2NO3)=(ZRX1*ZRX2)/(ZRX1+ZRX2)*10**(LOG10(0.35)/(1._JPRB+(lOG10(ZRX1/ZRX2)/ZNN_B)**2_JPIM))
    ZRX1=RATES_LUT(KHONOA,ITEMP)*ZAIRD
    ZRX2=RATES_LUT(KHONOB,ITEMP)

    PRR(JL,KHONO)= ZRX1/(1._JPRB+ZRX1/ZRX2)*0.81_JPRB**(1._JPRB/(1._JPRB+LOG10(ZRX1/ZRX2)**2_JPIM))
    ZRX1=RATES_LUT(KNO2HO2A,ITEMP)*ZAIRD
    ZRX2=RATES_LUT(KNO2HO2B,ITEMP)
    PRR(JL,KNO2HO2)=(ZRX1*ZRX2)/(ZRX1+ZRX2)*10**(log10(0.4_JPRB)/(1._JPRB+(log10(ZRX1/ZRX2)/ZNN)**2_JPIM))
    ZRX1=RATES_LUT(KMENO2A,ITEMP)*ZAIRD
    ZRX2=RATES_LUT(KMENO2B,ITEMP)
    PRR(JL,KMENO2)=ZRX1/(1._JPRB+ZRX1/ZRX2)*ZRX3**(1._JPRB/(1._JPRB+LOG10(ZRX1/ZRX2)**2_JPIM))
    ZRX1=RATES_LUT(KCOOHA,ITEMP)*ZAIRD
    ZRX2=RATES_LUT(KCOOHB,ITEMP)
    PRR(JL,KCOOH)=ZRX1/(1._JPRB+ZRX1/ZRX2)*ZRX3**(1._JPRB/(1._JPRB+LOG10(ZRX1/ZRX2)**2_JPIM))
    !JEW: now add the second term for CO + OH
    ZRX1=RATES_LUT(KCOOHC,ITEMP)
    ZRX2=RATES_LUT(KCOOHD,ITEMP)/ZAIRD
    PRR(JL,KCOOH)=PRR(JL,KCOOH)+(ZRX1/(1._JPRB+ZRX1/ZRX2)*ZRX3**(1._JPRB/(1._JPRB+LOG10(ZRX1/ZRX2)**2_JPIM)))
    ZRX1=RATES_LUT(KC61A,ITEMP)*ZAIRD
    ZRX2=RATES_LUT(KC61B,ITEMP)
    PRR(JL,KC61)=ZRX1/(1._JPRB+ZRX1/ZRX2)*ZRX3**(1._JPRB/(1._JPRB+LOG10(ZRX1/ZRX2)**2_JPIM))

    ZRX1=RATES_LUT(KOHC3H6A,ITEMP)*ZAIRD
    ZRX2=RATES_LUT(KOHC3H6B,ITEMP)
    PRR(JL,KOHC3H6)=ZRX1/(1._JPRB+ZRX1/ZRX2)*0.5_JPRB**(1._JPRB/(1._JPRB+LOG10(ZRX1/ZRX2)**2_JPIM))

    !
    ! 3.0 photolysis rates and 2 body reactions
    !
    PRR(JL,KNOO3)=RATES_LUT(KNOO3,ITEMP)
    PRR(JL,KHO2NO)=RATES_LUT(KHO2NO,ITEMP)
    PRR(JL,KMO2NO)=RATES_LUT(KMO2NO,ITEMP)
    PRR(JL,KNO2O3)=RATES_LUT(KNO2O3,ITEMP)
    PRR(JL,KNONO3)=RATES_LUT(KNONO3,ITEMP)

    ZRX1=RATES_LUT(KN2O5A,ITEMP)*ZAIRD*EXP(-11000._JPRB/ZTEMP)
    ZRX2=RATES_LUT(KN2O5B,ITEMP)*exp(-11080_JPRB/ZTEMP)
    PRR(JL,KN2O5)=(ZRX1*ZRX2)/(ZRX1+ZRX2)*10**(LOG10(0.35)/(1+(LOG10(ZRX1/ZRX2)/ZNN_B)**2_JPIM))

    PRR(JL,KOHHONO)=RATES_LUT(KOHHONO,ITEMP)
    PRR(JL,KHNO4OH)=RATES_LUT(KHNO4OH,ITEMP)

    ZRX1=RATES_LUT(KHNO4A,ITEMP)*ZAIRD
    ZRX2=RATES_LUT(KHNO4B,ITEMP)
    !
    !=((4.1e-5*exp(-10650/T))*M*(6.0e15*exp(-11170/T)))/((4.1e-5*exp(-10650/T)) &
    ! *M+(6.0e15*exp(-11170/T)))*10^(log10(0.4)/(1+(log10((4.1e-5*exp(-10650/T))&
    ! *M/(6.0e15*exp(-11170/T)))/(0.75-1.27*log10(0.4)))^2))
    !
    PRR(JL,KHNO4M)=(ZRX1*ZRX2)/(ZRX1+ZRX2)*10**(LOG10(0.4_JPRB)/(1._JPRB+(LOG10(ZRX1/ZRX2)/ZNN)**2_JPIM))

    PRR(JL,KODM)=RATES_LUT(KODM,ITEMP)
    PRR(JL,KH2OOD)=RATES_LUT(KH2OOD,ITEMP)
    PRR(JL,KO3HO2)=RATES_LUT(KO3HO2,ITEMP)
    !old formulation    PRR(JL,kcooh)=1.5E-13+9E-14*airp/101325.
    PRR(JL,KO3OH)=RATES_LUT(KO3OH,ITEMP)
    PRR(JL,KHPOH)=RATES_LUT(KHPOH,ITEMP)
    PRR(JL,KFRMOH)=RATES_LUT(KFRMOH,ITEMP)
    PRR(JL,KCH4OH)=RATES_LUT(KCH4OH,ITEMP)
    PRR(JL,KH2OH)=RATES_LUT(KH2OH,ITEMP)*550.E-9_JPRB*ZAIRD !h2=550ppbv JEW update
    PRR(JL,KOHMPER)=RATES_LUT(KOHMPER,ITEMP)
    PRR(JL,KOHROOH)=RATES_LUT(KOHROOH,ITEMP)
    PRR(JL,KMO2HO2A)=RATES_LUT(KMO2HO2A,ITEMP)
    PRR(JL,KMO2HO2B)=RATES_LUT(KMO2HO2B,ITEMP)
    PRR(JL,KMO2MO2)=RATES_LUT(KMO2MO2,ITEMP)
    PRR(JL,KHO2OH)=RATES_LUT(KHO2OH,ITEMP)
    PRR(JL,KHO2HO2)=(RATES_LUT(KHO2HO2A,ITEMP) +&
      &  RATES_LUT(KHO2HO2B,ITEMP)*ZAIRD)*(1._JPRB+RATES_LUT(KHO2HO2C,ITEMP)*ZH2OX)
    PRR(JL,KC41)=RATES_LUT(KC41,ITEMP)
    ! JEW for ALD take the average of the C2 and C3 rate co-efficients; measurements suggest
    ! that CH3CHO only comprises 50% of [higher aldehydes] - grossman et al, JGR, 2003.
    PRR(JL,KC43)=(RATES_LUT(KOHMCHO,ITEMP)+RATES_LUT(KOHMCH2CHO,ITEMP))/2._JPRB
    PRR(JL,KC44)=(RATES_LUT(KNO3MCHO,ITEMP)+ RATES_LUT(KNO3MCH2CHO,ITEMP))/2._JPRB
    PRR(JL,KC46)=RATES_LUT(KC46,ITEMP)
    ZRX1=RATES_LUT(KC47A,ITEMP)*ZAIRD
    ZRX2=RATES_LUT(KC47B,ITEMP)
    PRR(JL,KC47)=ZRX1/(1._JPRB+ZRX1/ZRX2)*0.3_JPRB**(1._JPRB/(1._JPRB+LOG10(ZRX1/ZRX2)**2_JPIM))
    !
    ! JEW : bring consistency for N2O5, HNO4 and PAN
    !
    !VH ZRX1=RATES_LUT(KC48A,ITEMP)*ZAIRD
    !VH ZRX2=RATES_LUT(KC48B,ITEMP)
    !VH PRR(JL,KC48)=ZRX1/(1._JPRB+ZRX1/ZRX2)*0.3_JPRB**(1._JPRB/(1._JPRB+LOG10(ZRX1/ZRX2)**2_JPIM))
    PRR(JL,KC48)=PRR(JL,KC47)/RATES_LUT(KC48,ITEMP)
    PRR(JL,KMENO2M)=PRR(JL,KMENO2)/RATES_LUT(KMENO2M,ITEMP)
    PRR(JL,KC49)=RATES_LUT(KC49,ITEMP)
    PRR(JL,KC50A)=RATES_LUT(KC50A,ITEMP)
    PRR(JL,KC50B)=RATES_LUT(KC50B,ITEMP)
    PRR(JL,KC52)=RATES_LUT(KC52,ITEMP)
    PRR(JL,KC53)=RATES_LUT(KC53,ITEMP)
    PRR(JL,KC54)=RATES_LUT(KC54,ITEMP)
    ! JEW updated as MVC and MACR are now the species ispd
    PRR(JL,KC57)=RATES_LUT(KOHOLE,ITEMP)
    PRR(JL,KC58)=RATES_LUT(KO3OLE,ITEMP)
    PRR(JL,KC59)=RATES_LUT(KNO3OLE,ITEMP)
    PRR(JL,KC62)=RATES_LUT(KC62,ITEMP)
    PRR(JL,KC73)=RATES_LUT(KC73,ITEMP)
    PRR(JL,KC76)=RATES_LUT(KC76,ITEMP)
    PRR(JL,KC77)=RATES_LUT(KC77,ITEMP)
    PRR(JL,KC78)=RATES_LUT(KC78,ITEMP)
    ! JEW use rates in cb05 for XO2 reactions
    PRR(JL,KC79)=RATES_LUT(KC79,ITEMP)
    PRR(JL,KC80)=RATES_LUT(KC80,ITEMP)
    PRR(JL,KC81)=RATES_LUT(KC81,ITEMP)
    PRR(JL,KC82)=RATES_LUT(KC82,ITEMP)
    PRR(JL,KC83)=RATES_LUT(KC83,ITEMP)
    PRR(JL,KC84)=RATES_LUT(KC84,ITEMP)
    PRR(JL,KC85)=RATES_LUT(KC85,ITEMP)
    ! fraction rjo3d=>OH
    PRJ(JL,JO3D)=PRJ(JL,JO3D)*PRR(JL,KH2OOD)*ZH2OX /&
      &        ( PRR(JL,KODM)*ZAIRD + PRR(JL,KH2OOD)*ZH2OX )
    ZRX1=RATES_LUT(KSO2OHA,ITEMP)*ZAIRD
    ZRX2=RATES_LUT(KSO2OHB,ITEMP)
    PRR(JL,KSO2OH)=ZRX1/(1._JPRB+ZRX1/ZRX2)*0.6_JPRB**(1._JPRB/(1._JPRB+LOG10(ZRX1/ZRX2)**2_JPIM))
    !
    ! dmsoha => so2
    ! dmsohb => 0.75 SO2 + 0.25 MSA
    !
    PRR(JL,KDMSOHA)=RATES_LUT(KDMSOHA,ITEMP)
    PRR(JL,KDMSOHB)=RATES_LUT(KDMSOHB,ITEMP)*ZO2/&
      &    (1.E30_JPRB+RATES_LUT(KDMSOHC,ITEMP)*ZO2)
    PRR(JL,KDMSNO3)=RATES_LUT(KDMSNO3,ITEMP)
    ! Sulphur reactions
    PRR(JL,KSO2O3G)=RATES_LUT(KSO2O3G,ITEMP)
    ! First ensure DP computation.. In troposphere assume kpp.eqn-reaction rate formulation "SO3 = SO4"
    ZRX1=RATES_LUT(KSO3H2O,ITEMP)*ZH2OX*ZH2OX
    PRR(JL,KSO3H2O)=ZRX1
    !PRR(JL,KSO3H2O)=RATES_LUT(KSO3H2O,ITEMP)*ZH2OX*ZH2OX
    PRR(JL,KOCSOH)=RATES_LUT(KOCSOH,ITEMP)
    !3.0 ammonia chemistry
    PRR(JL,KNH3OH)=RATES_LUT(KNH3OH,ITEMP)
    PRR(JL,KNH2NO)=RATES_LUT(KNH2NO,ITEMP)
    PRR(JL,KNH2NO2)=RATES_LUT(KNH2NO2,ITEMP)
    PRR(JL,KNH2HO2)=RATES_LUT(KNH2HO2,ITEMP)
    PRR(JL,KNH2O2)=RATES_LUT(KNH2O2,ITEMP)*ZO2
    PRR(JL,KNH2OH)=RATES_LUT(KNH2OH,ITEMP) ! This is actually zero.
    ! Hauglustaine et al (2014) mention NH2O and not NH2O2
    ! Reactions involving NH2O2 currently switched off:
    ! PRR(JL,KNH2O3)=RATES_LUT(KNH2O3,ITEMP)
    ! PRR(JL,KNH2O2NO)=RATES_LUT(KNH2O2NO,ITEMP)
    ! PRR(JL,KNH2O2O3)=RATES_LUT(KNH2O2O3,ITEMP)
    ! PRR(JL,KNH2O2HO2)=RATES_LUT(KNH2O2HO2,ITEMP)
    !
    ! O3P in molecules cm3
    !
    PRR(JL,KO3PO2)=RATES_LUT(KO3PO2,ITEMP)*ZAIRD
    PRR(JL,KO3PO3)=0._JPRB
    ! JEW: Reformulated June 2012
    !
    ! Runaway O3 above 50 hPa due to missing stratopsheric chemistry
    ! therefore use pressure as an index
    !
    ZO3_P=0._JPRB
    IF(PRSF1(JL,KL) > 5000_JPRB .AND. PRSF1(JL,KL) < 35000_JPRB) ZO3_P=(2*PRJ(JL,JO2)*ZO2)
    !
    !  [O2] not used in the ebi
    !
    PRJ(JL,JO2)=(ZO3_P/(PRR(JL,KO3PO2)*ZO2+PRR(JL,KO3PO3)*ZO3))*ZO2*PRR(JL,KO3PO2)
    PRR(JL,KO3PO3)=PRR(JL,KO3PO3)*ZO3_P
    !VH PRJ(JL,jo2)=0._JPRB
    !VH PRR(JL,ko3po3)=0._JPRB
    !
    !4.0 additional biogenic reactions
    PRR(JL,KOHHCOOH)=4.0E-13_JPRB
    PRR(JL,KOHCH3OH)=RATES_LUT(KOHCH3OH,ITEMP)

    !4.1 additional no3 peroxy reactions
    PRR(JL,KNO3HO2)=4.0E-12 ! IUPAC 2017

    PRR(JL,KNO3MO2)=1.2E-12_JPRB
    PRR(JL,KNO3C2O3)=4.0E-12_JPRB
    PRR(JL,KNO3XO2)=2.5E-12_JPRB ! Zaveri and Peters

    !4.2 additional C2 compounds
    PRR(JL,KOHMCOOH)=RATES_LUT(KOHMCOOH,ITEMP)
    PRR(JL,KOHC2H6)=RATES_LUT(KOHC2H6,ITEMP)
    PRR(JL,KOHETHOH)=RATES_LUT(KOHETHOH,ITEMP)
    PRR(JL,KOHC3H8)=RATES_LUT(KOHC3H8,ITEMP)
    PRR(JL,KO3C3H6)=RATES_LUT(KO3C3H6,ITEMP)
    PRR(JL,KNO3C3H6)=RATES_LUT(KNO3C3H6,ITEMP)
    ! SOA Aging
    PRR(JL,KSOA_AG)=8E-12_JPRB
    !
    !4.5  TERPENE reactions
    !
    PRR(JL,KOHTERP) = RATES_LUT(KOHTERP,ITEMP)
    PRR(JL,KO3TERP) = RATES_LUT(KO3TERP,ITEMP)
    PRR(JL,KNO3TERP)= RATES_LUT(KNO3TERP,ITEMP)

    ! s-1 decay time Rn222 to Pb210
    !PRR(JL,krn222)=2.11e-6_JPRB
    ! polmip value
    PRR(JL,KRN222)=2.10E-6_JPRB
    !
    !4.6 ISPD reactions
    !
    PRR(JL,KOHISPD) = (RATES_LUT(KOHMVK,ITEMP) +  RATES_LUT(KOHMACR,ITEMP))/2.0_JPRB
    PRR(JL,KO3ISPD) = (RATES_LUT(KO3MVK,ITEMP) +  RATES_LUT(KO3MACR,ITEMP))/2.0_JPRB
    PRR(JL,KNO3ISPD)= (6.0E-16_JPRB+3.5E-15_JPRB)/2.0_JPRB

    ! Explicit IC3H7O2 and HYPROPO2

    PRR(JL,KNOHYPROPO2)=RATES_LUT(KNOHYPROPO2,ITEMP)
    PRR(JL,KHO2HYPROPO2)=RATES_LUT(KHO2HYPROPO2,ITEMP)
    PRR(JL,KNOIC3H7O2)=RATES_LUT(KNOIC3H7O2,ITEMP)
    PRR(JL,KHO2IC3H7O2)=RATES_LUT(KHO2IC3H7O2,ITEMP)

    !4.7 acetone
    PRR(JL,KOHACET)=(RATES_LUT(KOHACETA,ITEMP)+RATES_LUT(KOHACETB,ITEMP))

    PRR(JL,KACO2HO2)=1.0E-11_JPRB
    PRR(JL,KACO2MO2)=3.8E-12_JPRB
    PRR(JL,KACO2NO)=8.0E-12_JPRB
    PRR(JL,KACO2XO2)=6.8E-14_JPRB ! CB05  peroxy-loss
    PRR(JL,KXO2XO2N)=6.8E-14_JPRB
    PRR(JL,KXO2N)=6.8E-14_JPRB
    !
    ! BB tracers
    !
    PRR(JL,KOHHCN)=RATES_LUT(KOHHCN,ITEMP)
    PRR(JL,KOHCH3CN)=RATES_LUT(KOHCH3CN,ITEMP)
    !
    !   TOL and XYL
    !
    PRR(JL,KOHTOL)=RATES_LUT(KOHTOL,ITEMP)
    PRR(JL,KO3TOL)=RATES_LUT(KO3TOL,ITEMP)
    PRR(JL,KNO3TOL)=RATES_LUT(KNO3TOL,ITEMP)


    PRR(JL,KOHXYL)=RATES_LUT(KOHXYL,ITEMP)
    PRR(JL,KO3XYL)=(RATES_LUT(KCXYLO3A,ITEMP)+RATES_LUT(KCXYLO3B,ITEMP)+(0.796*RATES_LUT(KCXYLO3A,ITEMP)))/3.0_JPRB
    PRR(JL,KNO3XYL)=RATES_LUT(KNO3XYL,ITEMP)

    PRR(JL,KAROO2NO)=RATES_LUT(KAROO2NO,ITEMP)
    PRR(JL,K2AROO2)=RATES_LUT(K2AROO2,ITEMP)
    PRR(JL,KAROO2HO2)=RATES_LUT(KAROO2HO2,ITEMP)
    PRR(JL,KAROO2XO2)=RATES_LUT(KAROO2XO2,ITEMP)

    PRR(JL,KISOPBO2A)=RATES_LUT(KISOPBO2A,ITEMP)
    PRR(JL,KISOPBO2B)=RATES_LUT(KISOPBO2B,ITEMP)
    PRR(JL,KISOPDO2A)=RATES_LUT(KISOPDO2A,ITEMP)
    PRR(JL,KISOPDO2B)=RATES_LUT(KISOPDO2B,ITEMP)
    PRR(JL,KISOPBO2HO2)=RATES_LUT(KISOPBO2HO2,ITEMP)
    PRR(JL,KISOPBO2NO)=RATES_LUT(KISOPBO2NO,ITEMP)
    PRR(JL,KISOPDO2HO2)=RATES_LUT(KISOPDO2HO2,ITEMP)
    PRR(JL,KISOPDO2NO)=RATES_LUT(KISOPDO2NO,ITEMP)

    PRR(JL,KISOPOOHOH)=RATES_LUT(KISOPOOHOH,ITEMP)
    PRR(JL,KHPALD1OH)=RATES_LUT(KHPALD1OH,ITEMP)
    PRR(JL,KHPALD2OH)=RATES_LUT(KHPALD2OH,ITEMP)
    PRR(JL,KGLYOH)=RATES_LUT(KGLYOH,ITEMP)
    PRR(JL,KGLYALDOH)=RATES_LUT(KGLYALDOH,ITEMP)
    PRR(JL,KHYACOH)=RATES_LUT(KHYACOH,ITEMP)

ENDDO


!5.0 calculate heterogeneous removal constants of n2o5, ho2 and no3
!*   via kt=(r/Dg + 4/vgamma)^{-1} * A [cm2/cm3]
!*   on liquid and ice cloud and aerosol
!*   to be included in gas phase chemistry

DO JL = KIDIA,KFDIA

    ZR_CLOUD = 0._JPRB
    ZR_ICE = 0._JPRB
    !* 5.1 Liquid / Ice Cloud uptake

    IF(PCLOUD_REFF(JL,KL)>=9.0_JPRB) ZR_CLOUD = (PCLOUD_REFF(JL,KL)-2.0_JPRB)/1E4_JPRB
    IF(PCLOUD_REFF(JL,KL)>=6.0_JPRB  .AND. PCLOUD_REFF(JL,KL)<9.0_JPRB) THEN
      ZR_CLOUD=(PCLOUD_REFF(JL,KL)-1.0_JPRB)/1E4_JPRB
    ENDIF
    IF(PCLOUD_REFF(JL,KL) < 6.0_JPRB) ZR_CLOUD = (PCLOUD_REFF(JL,KL)-0.25_JPRB)/1E4_JPRB

    ZTEMP=PTP(JL,KL)
    ! air density (kg/m3)
    ZRHO = PRSF1(JL,KL)/(RD*ZTEMP)
! added - is this correct?
    ZAIRD = 7.24291E16_JPRB*PRSF1(JL,KL)/ZTEMP
    ZSAD_ICE = 0._JPRB

    ! default ice particle radius [cm]
    ZR_ICE=5.E-3_JPRB
    IF(PIP(JL,KL)>1.0E-10_JPRB .AND. PAP(JL,KL) > 1E-2_JPRB ) THEN
      !  [kg/kg] -> [gr/m3]
      ZM_ICE=PIP(JL,KL)*ZRHO * 1E3_JPRB / PAP(JL,KL)

      ZSAD_ICE=1.0E-4_JPRB*ZM_ICE**0.9_JPRB

      ! calculate the ZR_ICE using the relationship of Fu (1996)
      ZR_ICE=(0.8179_JPRB)*((ZM_ICE*1E-6_JPRB)/ZSAD_ICE)

      ! Maximize radius...
      ZR_ICE = MIN(ZR_ICE,1._JPRB)

      ! the value adopted in von Kuhlmann and Lawrence is too low
      ZSAD_ICE=10_JPRB*ZSAD_ICE
    ENDIF


    ! kg/kg -> g/cm^3
    ZM_CLD = 0._JPRB
    IF (PAP(JL,KL) > 1E-2_JPRB)  ZM_CLD=PLP(JL,KL)*ZRHO * 1E-3_JPRB / PAP(JL,KL)

    ! cloud SAD  [cm^2/cm^3] in cloudy part only
    ZSAD_CLOUD=3.0_JPRB * ZM_CLD / ( ZRHO_H2O * ZR_CLOUD )

    ! The surface area density for ice particles in now linked to
    ! the IWC by the parameterization of Heymsfield and McFarquar (1996)
    ! and the effective radii from Fu (1996)
    ZDG=0.1E5_JPRB/PRSF1(JL,KL)     !simple approximation for diffusion coeff. [cm2/s]
    ZG_N2O5_L = 2.7E-5_JPRB*EXP(1800._JPRB/ZTEMP)
    ZC_THERMAL_N2O5 = SQRT (8_JPRB*1.38E-23_JPRB *ZTEMP /(RPI* ZM_N2O5 ) )  * 1E2_JPRB ! m/s -> cm/sec
    ZC_THERMAL_HO2  = SQRT (8_JPRB*1.38E-23_JPRB *ZTEMP /(RPI* ZM_HO2  ) )  * 1E2_JPRB ! m/s -> cm/sec
    ZKT_LIQ_N2O5 =1._JPRB/(ZR_CLOUD/ZDG  + 4._JPRB/(ZC_THERMAL_N2O5*ZG_N2O5_L))
    ZKT_ICE_N2O5 =1._JPRB/(ZR_ICE/ZDG    + 4._JPRB/(ZC_THERMAL_N2O5* ZG_N2O5_I))

    ! Scaling factor accounting for subgrid-scale mixing within time step.
    ! ZSGS_MIX=0 if CC=0, and 1 if CC=1; ZSGS_MIX is lower than CC for 0<CC<1 ,reducing the effective reaction rate

    ZCC = MIN(0.99_JPRB,PAP(JL,KL))

! Sub grid scale effects for HO2 and N2O5 - switched off i, i.e. set to 1.0
  ! ZSGS_MIX for HO2: inspired by Cariolle (JGR paper)
  ! ZSGS_MIX_HO2=ZCC * (1._JPRB - EXP(-0.1_JPRB/(1._JPRB-ZCC)))
    ZSGS_MIX_HO2=0.0_JPRB
    IF (ZCC > 1E-2_JPRB) THEN
!      ZSGS_MIX_HO2 = ZCC / ((1.-ZCC)*ZHO2_Ratio+ZCC)
!     no sub scale effect
!VH       ZSGS_MIX_HO2 = 1.0_JPRB
      ZSGS_MIX_HO2 = ZCC
    ENDIF


!    No sub scale effect
    ZSGS_MIX_N2O5=ZCC
    ! Compute N2O5 reaction on cloud and ice:
    PRR(JL,KN2O5L)=(ZKT_LIQ_N2O5*ZSAD_CLOUD+ZKT_ICE_N2O5*ZSAD_ICE)* ZSGS_MIX_N2O5  ! Apply SGS scaling

    ! Effective Henry's constant for HO2, enhanced by H+
    ZTR=(1.0_JPRB/ZTEMP-1.0_JPRB/298.0_JPRB)
    ZH_HO2=ZHNRYA*EXP(ZHNRYB * ZTR)
    ZH_EFF=ZH_HO2*(1._JPRB+ZK_EQ/ZHPL_FIX)

    ! Effective reaction rate defining partitioning between HO2(aq) and O2-(aq)
    ZK1 = ZK1A* EXP(ZK1B/ZTEMP)
    ZK2 = ZK2A* EXP(ZK2B/ZTEMP)
    ZNUM  = ZK1 + ZK_EQ/ZHPL_FIX * ZK2
    ZDENOM= ( 1_JPRB +  ZK_EQ / ZHPL_FIX)**2_JPIM
    ZK_EFF =ZNUM/ZDENOM

    ! Now compute eq. 7 from Thornton et al.
    ! For stability purposes introduce lower bound of concentration of 1.
    ZCONC_HO2 = MAX(PHO2(JL) / YCHEM(IHO2)%RMOLMASS *ZAIRD*RMD ,1.0_JPRB)
    ZNUM = 3._JPRB * ZC_THERMAL_HO2 * RNAVO
    ZDENOM = 8000._JPRB * (ZH_EFF *8.314E-2_JPRB*ZTEMP)**2_JPIM * ZK_EFF * ZCONC_HO2 *ZR_CLOUD
    ZG_HO2_L=1._JPRB/( ZALPH_HO2I + ZNUM/ZDENOM  )
    ZG_HO2_L=MAX(1E-10_JPRB,MIN(1._JPRB,ZG_HO2_L))


    ! Apply HO2 uptake on cloud & ice?
    !VH IF (KCHEM_HET >= 2_JPIM) THEN
      ZKT_LIQ_HO2=1._JPRB/(ZR_CLOUD/ZDG+4._JPRB/(ZC_THERMAL_HO2*ZG_HO2_L))
      ZKT_ICE_HO2=1._JPRB/((ZR_ICE/ZDG)+(4._JPRB/(ZC_THERMAL_HO2*ZG_HO2_I)))
      PRR(JL,KHO2L)=(ZKT_LIQ_HO2*ZSAD_CLOUD+ZKT_ICE_HO2*ZSAD_ICE) * ZSGS_MIX_HO2
    !VH ELSE
    !VH  PRR(JL,kho2l)=0._JPRB !cloud
    !VH ENDIF
    ! export cloud surface area
    PSAD_CLD(JL,KL)=ZSAD_CLOUD
    ! export ice SAD...
    PSAD_ICE(JL,KL)=ZSAD_ICE


    !* 5.2 uptake on aerosol

    PRR(JL,KN2O5_AER) = 0._JPRB
    PRR(JL,KHO2_AER) = 0._JPRB
    PRR(JL,KNO3_AER) = 0._JPRB
    ZSAD_AERT=0.

    ! relative humidity
    ZTR = 1._JPRB - 373.15_JPRB/ZTEMP
    ZWV = EXP((((-.1299_JPRB*ZTR-.6445_JPRB)*ZTR-1.976_JPRB)*ZTR+13.3185_JPRB)*ZTR)
    ZRH = PQP(JL,KL)*PRSF1(JL,KL)*XMAIR /(XMH2O* 1013.25_JPRB*ZWV)
    ZRH = MAX(MIN(ZRH,99.9_JPRB),0.01_JPRB)   !limit rh between 0-100%

    ! Get table index...
    DO JTAB=1,12
      IF (ZRH > RRHTAB(JTAB)) THEN
        IRH = JTAB
      ENDIF
    ENDDO

    IF ( YDCHEM%LCHEM_AEROI .AND. TRIM(AERO_SCHEME)=="aer" ) THEN

      DO JAER=1,KAEROC

        ! Baseline values
        ZVFRAC=1.0_JPRB
        ZRHO_AER=ZRHO_SO4

        ! standard gamma values
        ZG_NO3=0.01
        ZG_N2O5=ZG_N2O5_SO4

        !kg/kg->g/cm^3 baseline value, if no trace-gas specific modification
        ZM_AER=PAEROP(JL,KL,JAER)*ZRHO * 1E-3_JPRB

        ITYP=YAERO_DESC(JAER)%NTYP
        IBIN=YAERO_DESC(JAER)%NBIN
        IF (ITYP==1 ) THEN
          ! Sea salt...
          ! kg/kg -> g/cm^3 ... specific for sea salt
          ZM_AER=PAEROP(JL,KL,JAER)*RSS_DRY_MASSFAC*ZRHO * 1E-3_JPRB
          ! effective wet particle radius
          ZR_EFF=ZR_SS(IBIN)*RSS_DRY_DIAFAC * ZRH_GROWTH_SS(IRH)
          ! Aerosol density
          ZRHO_AER=ZRHO_SS
          ! Volume fraction dry over full aerosol particle
          ZVFRAC = 1.0_JPRB / ZRH_GROWTH_SS(IRH)**3_JPIM
          ! Gamma values..
          ZG_N2O5=ZG_N2O5_SS

        ELSEIF (ITYP==2) THEN
          ! Dust... (independent to RH)
          ZR_EFF=ZR_DD(IBIN)
          ! Aerosol density
          ZRHO_AER=ZRHO_DD
          ! Volume fraction dry over full aerosol particle
          ZVFRAC = 1.0_JPRB
          ! Gamma values..
          IF (ZRH < 50._JPRB) THEN
            ! For now low reaction for NO3 adopted
            ZG_NO3=1E-4
          ENDIF
          ZG_N2O5=ZG_N2O5_DD

        ELSEIF (ITYP==3) THEN
          ! POM
          ZR_EFF=ZR_OM * ZRH_GROWTH_OM(IRH)
          ! Aerosol density
          ZRHO_AER=ZRHO_OM
          ! Volume fraction dry over full aerosol particle
          ZVFRAC = 1.0_JPRB / ZRH_GROWTH_OM(IRH)**3_JPIM
          ! Gamma values..
          IF (IBIN == 1) THEN
            ! High reaction on hydrophylic
            ZG_N2O5=ZG_N2O5_OM
          ELSE
            ! For now hydrophobic reaction on OM not included
            CYCLE
            ! un future consider some reaction
            ! ZG_NO3 =0.0001
            ! ZG_N2O5=0.0001
          ENDIF

        ELSEIF (ITYP==4) THEN
          ! Black Carbon... (only hydrophylic...)

          ! Volume fraction dry over full aerosol particle
          ZR_EFF=ZR_BC * ZRH_GROWTH_BC(IRH)
          ! Aerosol density
          ZRHO_AER=ZRHO_BC
          ! Volume fraction dry over full aerosol particle
          ZVFRAC = 1.0_JPRB / ZRH_GROWTH_BC(IRH)**3_JPIM
          ! Gamma values..
          IF (IBIN == 1) THEN
            ! reaction on hydrophylic
            ZG_N2O5=ZG_N2O5_BC
          ELSE
            ! So far hydrophobic BC not included
            CYCLE
            ! In future consider reaction on hydrophobic
            ! ZG_NO3 =0.0001
            ! ZG_N2O5=0.0001
          ENDIF

        ELSEIF (ITYP==5) THEN
          ! Sulphate

          ! Volume fraction dry over full aerosol particle
          ZR_EFF=ZR_SO4 * ZRH_GROWTH_SO4(IRH)
          ! Aerosol density
          ZRHO_AER=ZRHO_SO4
          ! Volume fraction dry over full aerosol particle
          ZVFRAC = 1.0_JPRB / ZRH_GROWTH_SO4(IRH)**3_JPIM
          ! Gamma values..
          ZG_N2O5=ZG_N2O5_SO4

        ELSEIF ( ITYP==6 .OR. ITYP==7) THEN
          ! Nitrate or Ammonium

          ! Volume fraction dry over full aerosol particle
          ZR_EFF=ZR_NO3_A * ZRH_GROWTH_SO4(IRH)
          ! Aerosol density
          ZRHO_AER=ZRHO_SO4
          ! Volume fraction dry over full aerosol particle
          ZVFRAC = 1.0_JPRB / ZRH_GROWTH_SO4(IRH)**3_JPIM
          ! Gamma values..
          ZG_N2O5=ZG_N2O5_NO3A

        ELSEIF (ITYP==8) THEN
          ! Secondary organic aerosol

          ! Volume fraction dry over full aerosol particle
          ZR_EFF=ZR_SOA(IBIN) * ZRH_GROWTH_SOA(IRH)
          ! Aerosol density
          ZRHO_AER=ZRHO_SOA

          ! Volume fraction dry over full aerosol particle
          ZVFRAC = 1.0_JPRB / ZRH_GROWTH_SOA(IRH)**3_JPIM
          ! Gamma values..
          ZG_N2O5=ZG_N2O5_OM

        ENDIF

        ! particle density (combination of dry aerosol + water)
        ZRHO_P = ZRHO_H2O*(1.0_JPRB-ZVFRAC) + ZVFRAC*ZRHO_AER
        !aerosol SAD (cm^2/cm^3)
        ZSAD_AER=3.0_JPRB * ZM_AER / ( ZRHO_P * ZR_EFF )
        ! accumulate total SAD
        ZSAD_AERT=ZSAD_AERT+ZSAD_AER

        ! NO3
        ZKTNO3_AER=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(3.0E4_JPRB*ZG_NO3)))
        PRR(JL,KNO3_AER) = PRR(JL,KNO3_AER) + ZKTNO3_AER*ZSAD_AER

        ! N2O5
        ZKT_AER_N2O5=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(ZC_THERMAL_N2O5*ZG_N2O5)))
        PRR(JL,KN2O5_AER) = PRR(JL,KN2O5_AER) + ZKT_AER_N2O5*ZSAD_AER

        ! HO2
        IF (ITYP == 2) THEN
          ! Bedjanian et al., ACPD 2013, suggest an inverse dependency on RH...
          ! with lower uptake with higher RH. This is not accounted for here.
          ZKT_AER_HO2=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(ZC_THERMAL_HO2*ZG_HO2_DD)))
        ELSE
          ! Following Thornton et al.
          ZDENOM = 8000._JPRB * (ZH_EFF *8.314E-2_JPRB*ZTEMP)**2_JPIM * ZK_EFF * ZCONC_HO2 *ZR_EFF
          ZG_HO2_AER=1._JPRB/( ZALPH_HO2I + ZNUM/ZDENOM  )
          ZG_HO2_AER=MAX(1E-10_JPRB,MIN(1._JPRB,ZG_HO2_AER))
          ZKT_AER_HO2=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(ZC_THERMAL_HO2*ZG_HO2_AER)))
        ENDIF
        PRR(JL,KHO2_AER) = PRR(JL,KHO2_AER) + ZKT_AER_HO2*ZSAD_AER
      ENDDO

      ! Output total aerosol SAD...
      PSAD_AER(JL,KL)=ZSAD_AERT




    ELSEIF ( YDCHEM%LCHEM_AEROI .AND. TRIM(AERO_SCHEME)=="glomap" ) THEN
      ! Alternative computation, in case of glomap aerosol scheme

      DO JMODE = 1,KMODES
        IF (MODE(JMODE)) THEN
          IF (PND(JL,KL,JMODE) > 1E-10_JPRB .AND. PWETDIAM(JL,KL,JMODE) > 1E-12_JPRB ) THEN
            ! Median wet diameter (m -> cm)
            ZDIAMW=PWETDIAM(JL,KL,JMODE)*1E2_JPRB ! *EXP((-(SIGMAG(JMODE))**2_JPIM)/2._JPRB)
            ! SAD = pi*N * D^2*exp(2 sigma_g) (cm2/cm3)
            ZSAD_AER=(RPI*PND(JL,KL,JMODE) * ZDIAMW**2_JPIM ) ! EXP(2*(SIGMAG(JMODE))**2)
            !
            ZR_EFF=ZDIAMW*0.5_JPRB

           ! alternative.. leads to too high SAD?!
           ! SAD = pi*N * D^2*exp(2 sigma_g) (cm2/cm3)
           ! i.e. conversion from 'mean' to 'median'
           !ZDIAMW=PWETDIAM(JL,KL,JMODE)*1E2_JPRB  *EXP((-(SIGMAG(JMODE))**2_JPIM)/2._JPRB)
           !ZSAD_AER=(RPI*PND(JL,KL,JMODE) * ZDIAMW**2_JPIM ) * EXP(2*(SIGMAG(JMODE))**2_JPIM)

            ! accumulate total SAD
            ZSAD_AERT=ZSAD_AERT+ZSAD_AER

            !N2O5...
            ZKT_AER_N2O5=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(ZC_THERMAL_N2O5*ZG_N2O5_MODE(JMODE))))
            PRR(JL,KN2O5_AER) = PRR(JL,KN2O5_AER) + ZKT_AER_N2O5*ZSAD_AER

            ! HO2...
            ! Following Thornton et al.
            ! ZCONC_HO2 = MAX(PHO2(JL) / YCHEM(IHO2)%RMOLMASS *ZAIRD*RMD ,0._JPRB)
            ! ZNUM = 3. * ZC_THERMAL_HO2 * RNAVO
            ZDENOM = 8000._JPRB * (ZH_EFF *8.314E-2_JPRB*ZTEMP)**2_JPIM * ZK_EFF * ZCONC_HO2 *ZR_EFF
            ZG_HO2_AER=1._JPRB/( ZALPH_HO2I + ZNUM/ZDENOM  )
            ZG_HO2_AER=MAX(1E-10_JPRB,MIN(1._JPRB,ZG_HO2_AER))
            ZKT_AER_HO2=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(ZC_THERMAL_HO2*ZG_HO2_AER)))
            PRR(JL,KHO2_AER) = PRR(JL,KHO2_AER) + ZKT_AER_HO2*ZSAD_AER

            IF (ZG_NO3_MODE(JMODE) > 0._JPRB) THEN
              ZKTNO3_AER=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(3.0E4_JPRB*ZG_NO3_MODE(JMODE))))
              PRR(JL,KNO3_AER) = PRR(JL,KNO3_AER) + ZKTNO3_AER*ZSAD_AER
            ENDIF
          ENDIF ! PND > 0
        ENDIF
      ENDDO

      !  additionally NO3_A and NH4 aerosol...
      !  Conform that these chemistry tracers are used in GLOMAP.
      ! kg/kg -> g/cm^3
      ZM_AER=(PNO3_A(JL)+PNH4(JL))*ZRHO * 1E-3_JPRB

      ! Volume fraction dry over full aerosol particle. Assume similar to SO4...
      ZR_EFF=ZR_NO3_A * ZRH_GROWTH_SO4(IRH)
      ZVFRAC = 1.0/ ZRH_GROWTH_SO4(IRH)**3_JPIM
      ZRHO_P = ZRHO_H2O*(1.0-ZVFRAC) + ZVFRAC*ZRHO_SO4

      ZSAD_AER=3.0_JPRB * ZM_AER / ( ZRHO_P * ZR_EFF )
      ZSAD_AERT=ZSAD_AERT+ZSAD_AER


      ZKT_AER_N2O5=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(ZC_THERMAL_N2O5*ZG_N2O5_NO3A)))

      ! Apply pozzoli et al. 2008:
      ! ZKT_AER_N2O5=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(2.4e4_JPRB*Zg_n2o5_Liao)))
      PRR(JL,KN2O5_AER) = PRR(JL,KN2O5_AER) + ZKT_AER_N2O5*ZSAD_AER

      ZDENOM = 8000._JPRB * (ZH_EFF *8.314E-2_JPRB*ZTEMP)**2_JPIM * ZK_EFF * ZCONC_HO2 *ZR_EFF
      ZG_HO2_AER=1._JPRB/( ZALPH_HO2I + ZNUM/ZDENOM  )
      ZG_HO2_AER=MAX(1E-10_JPRB,MIN(1._JPRB,ZG_HO2_AER))
      ZKT_AER_HO2=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(ZC_THERMAL_HO2*ZG_HO2_AER)))

      PRR(JL,KHO2_AER) = PRR(JL,KHO2_AER) + ZKT_AER_HO2*ZSAD_AER
      !VH ENDIF

      ZG_NO3=0.01
      ZKTNO3_AER=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(3.0E4_JPRB*ZG_NO3)))
      PRR(JL,KNO3_AER) = PRR(JL,KNO3_AER) + ZKTNO3_AER*ZSAD_AER

      ! Output total aerosol SAD...
      PSAD_AER(JL,KL)=ZSAD_AERT




    ELSE ! LCHEM_AEROI false

      !include chemical SO4 for heterogeneous reactions, as well as NO3 and NH4:
      ! kg/kg -> g/cm^3
       ZM_AER=(PSO4(JL)+PNO3_A(JL)+PNH4(JL))*ZRHO * 1E-3_JPRB

      ! Volume fraction dry over full aerosol particle. Assume similar to SO4...
      ZR_EFF=ZR_NO3_A * ZRH_GROWTH_SO4(IRH)
      ZVFRAC = 1.0/ ZRH_GROWTH_SO4(IRH)**3_JPIM
      ZRHO_P = ZRHO_H2O*(1.0-ZVFRAC) + ZVFRAC*ZRHO_SO4

      ZSAD_AER=3.0_JPRB * ZM_AER / ( ZRHO_P * ZR_EFF )
      ZSAD_AERT=ZSAD_AERT+ZSAD_AER

      ZKT_AER_N2O5=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(ZC_THERMAL_N2O5*ZG_N2O5_NO3A)))

      ! Apply pozzoli et al. 2008:
      ! ZKT_AER_N2O5=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(2.4e4_JPRB*Zg_n2o5_Liao)))
      PRR(JL,KN2O5_AER) = PRR(JL,KN2O5_AER) + ZKT_AER_N2O5*ZSAD_AER

      ZDENOM = 8000._JPRB * (ZH_EFF *8.314E-2_JPRB*ZTEMP)**2_JPIM * ZK_EFF * ZCONC_HO2 *ZR_EFF
      ZG_HO2_AER=1._JPRB/( ZALPH_HO2I + ZNUM/ZDENOM  )
      ZG_HO2_AER=MAX(1E-10_JPRB,MIN(1._JPRB,ZG_HO2_AER))
      ZKT_AER_HO2=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(ZC_THERMAL_HO2*ZG_HO2_AER)))

      PRR(JL,KHO2_AER) = PRR(JL,KHO2_AER) + ZKT_AER_HO2*ZSAD_AER
      !VH ENDIF

      ZG_NO3=0.01
      ZKTNO3_AER=1._JPRB/((ZR_EFF/ZDG)+(4._JPRB/(3.0E4_JPRB*ZG_NO3)))
      PRR(JL,KNO3_AER) = PRR(JL,KNO3_AER) + ZKTNO3_AER*ZSAD_AER

      ! Output total aerosol SAD...
      PSAD_AER(JL,KL)=ZSAD_AERT

    ENDIF


    ! obsolete...
    PRR(JL,KN2O5AQ)=0._JPRB
    PRR(JL,KNH3SO4)=0._JPRB
    !
ENDDO


END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TM5_CALRATES',1,ZHOOK_HANDLE )

END SUBROUTINE TM5_CALRATES
