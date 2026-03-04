! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_SO2SO4_V2 &
  &( YDRIP, KIDIA , KFDIA , KLON  , KLEV  ,         &
  &  PTSPHY, PTP   , PRSF1 , PNEB  , PQLI  , PGELAT, PGELAM, &
  &  PSO2  , PITSO2, POH, PO3, PH2O2    ,         &
  &  PTSO2 , PTSO4, PFSO2, PFSO4, PDP )

!*** *AER_SO2SO4_V2* - GAS-TO-PARTICLE (SULPHATE AEROSOLS)

!**   INTERFACE.
!     ----------
!          *AER_SO2SO4_V2* IS CALLED BY *AER_PHY3*.

!     AUTHOR.
!     -------
!        Josué Bock, Samuel Remy

!     SOURCE.
!     -------

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : February 2017
!        28-Mar-2017 : in this test version, oxidants are not depleted
!                      reason: C-IFS climatologies are already depleted after S(iv) oxidation, thus input oxidants concentration
!                              can be considered as "background" concentrations, which should not be depleted during oxidation

!     FUTURE IMPROVEMENTS / To Do LIST
!     --------------------------------
!        - check that SO2 (g) + OH (g) is the limiting step: see notebook p. 93
!        - add SO2 (g) + O3 (g) if necessary
!        - transition metals?
!        - het. react? see for instance Li et al, ACP 2017, 10.5194/acp-17-3301-2017
!        - notebook p. 24: at low H2O2, HNO4 plays a significant role in SO2 oxidation
!        - notebook p. 24: oxidants limitation
!-----------------------------------------------------------------------

USE PARKIND1    ,ONLY : JPIM     ,JPRB
USE YOMHOOK     ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMRIP   , ONLY : TRIP

USE YOMCST      ,ONLY : &
  & RNAVO, & ! Avogadro number      [mol**-1]
  & R, &     ! gas constant         [J / K / mol]
  & RD,RG         ! dry air gas constant [J / K / kg]


IMPLICIT NONE


!*       0.1   ARGUMENTS
!              ---------

TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV

REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY           ! timestep                               [s]
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTP  (KLON,KLEV) ! temperature                            [K]
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNEB (KLON,KLEV) ! fractional cloudiness                  [-]
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSF1(KLON,KLEV) ! pressure                               [Pa]
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQLI (KLON,KLEV) ! liquid water                           [kg/kg]
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM (KLON)    ! Longitude
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAT (KLON)    ! Latitude

REAL(KIND=JPRB)   ,INTENT(IN)    :: PSO2(KLON,KLEV)   ! SO2 mass mixing ratio                 [kg / kg(air)]
REAL(KIND=JPRB)   ,INTENT(IN)    :: PITSO2(KLON,KLEV) ! previous tendency for SO2             [kg / kg(air) / s]
REAL(KIND=JPRB)   ,INTENT(IN)    :: POH(KLON,KLEV)    ! oxidants [kg / kg]
REAL(KIND=JPRB)   ,INTENT(IN)    :: PH2O2(KLON,KLEV)    ! oxidants [kg / kg]
REAL(KIND=JPRB)   ,INTENT(IN)    :: PO3(KLON,KLEV)    ! oxidants [kg / kg]
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDP(KLON,KLEV) 


REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PTSO2(KLON,KLEV), PTSO4(KLON,KLEV) ! new tendencies       [kg / kg(air) / s]
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSO2(KLON), PFSO4(KLON)


!*       0.3   LOCAL PARAMETERS
!              ----------------

!==================!
! Henry's law data !
!==================!
! reference temperature for Henry's law solubility coefficients                                   [K]
REAL(KIND=JPRB), PARAMETER :: ZTREF_H = 298.15_JPRB
! inverse reference temperature for Henry's law solubility coefficients, and some rate constants  [K**-1]
REAL(KIND=JPRB), PARAMETER :: ZITREF_H = 1._JPRB / ZTREF_H

! Henry's law solubility coefficients at T=Tref, and their temperature dependency (if relevant)
REAL(KIND=JPRB), PARAMETER :: ZHCP_H2O2_REF = 9.1E+2_JPRB ! H^cp(H2O2), Sander, ACP 2015          [mol/m3/Pa]
REAL(KIND=JPRB), PARAMETER :: ZHCP_H2O2_TD = 6600._JPRB   ! temperature dependency for H^cp(H2O2) [K]
REAL(KIND=JPRB), PARAMETER :: ZHCP_O3_REF = 1.0E-4_JPRB   ! H^cp(O3), Sander, ACP 2015            [mol/m3/Pa]
REAL(KIND=JPRB), PARAMETER :: ZHCP_O3_TD = 2800._JPRB     ! temperature dependency for H^cp(O3)   [K]
REAL(KIND=JPRB), PARAMETER :: ZHCP_OH_REF = 3.8E-1_JPRB   ! H^cp(OH), Sander, ACP 2015            [mol/m3/Pa]
!REAL(KIND=JPRB), PARAMETER :: ZHCP_OH_TD =                ! no temperature dependency for H^cp(OH)
REAL(KIND=JPRB), PARAMETER :: ZHCP_SO2_REF = 1.3E-2_JPRB  ! H^cp(SO2), Sander, ACP 2015           [mol/m3/Pa]
REAL(KIND=JPRB), PARAMETER :: ZHCP_SO2_TD = 2100._JPRB    ! temperature dependency for H^cp(SO2)  [K]

!=================================!
! Acid-base equilibrium constants !
!=================================!
! Keq_1(SO2(aq)<=>HSO3-)   Seinfeld & Pandis, 1998, table 6.A.1 p.394
REAL(KIND=JPRB), PARAMETER :: ZKEQ1_SO2_REF = 1.3E-2_JPRB ! Keq_1(SO2(aq)<=>HSO3-) at T=Tref      [mol / L]
REAL(KIND=JPRB), PARAMETER :: ZKEQ1_SO2_TD = 1960._JPRB   ! temperature dependency for Keq_1      [K]
! Keq_2(HSO3-<=>SO3=)      Seinfeld & Pandis, 1998, table 6.A.1 p.394
REAL(KIND=JPRB), PARAMETER :: ZKEQ2_SO2_REF = 6.6E-8_JPRB ! Keq_2(HSO3-<=>SO3=) at T=Tref         [mol / L]
REAL(KIND=JPRB), PARAMETER :: ZKEQ2_SO2_TD = 1500._JPRB   ! temperature dependency for Keq_2      [K]


!=========================!
! Reaction rate constants !
!=========================!
! <jjb 10-02-2017  technical note: RNAVO is declared in YOMCST, but initialised in sucst.
!                  Thus, it is not a parameter, and this hinders using it in parameters initialisations.
!                  As a consequence, the next 3 variables cannot be defined as parameters (while they are parameters...)
!                  nor be initialised in the declaration part of the subroutine.
!  jjb 10-02-2017>

! conversion factor, multiply by ZCONV1 to convert cm**3/molec into m**3/mol
!REAL(KIND=JPRB), PARAMETER :: ZCONV1 = 1.E-6_JPRB * RNAVO ! doesn't work, cf comment above
!REAL(KIND=JPRB)            :: ZCONV1 = 1.E-6_JPRB * RNAVO ! doesn't work, cf comment above
REAL(KIND=JPRB)            :: ZCONV1
! conversion factor, multiply by ZCONV3 to convert L into m3
!   for instance, 2nd order reaction rate in aqueous phase [L / mol / s] into SI units [m3 / mol / s]
REAL(KIND=JPRB)            :: ZCONV3 = 1.E-3_JPRB ! [m**3 / L]


! k_OH reaction rate, low pressure limit, at 300 K, converted in [m**6 / (mol**2 * s)]
!REAL(KIND=JPRB), PARAMETER :: ZKOH_LOW300 = 3.3E-31_JPRB * ZCONV1 * ZCONV1 ! doesn't work, cf comment above
!REAL(KIND=JPRB)            :: ZKOH_LOW300 = 3.3E-31_JPRB * ZCONV1 * ZCONV1 ! doesn't work, cf comment above
REAL(KIND=JPRB)            :: ZKOH_LOW300
! k_OH reaction rate, high pressure limit, no temp. dependency, converted in [m**3 / (mol * s)]
!REAL(KIND=JPRB), PARAMETER :: ZKOH_HIGH   = 1.6E-12_JPRB * ZCONV1 ! doesn't work, cf comment above
!REAL(KIND=JPRB)            :: ZKOH_HIGH   = 1.6E-12_JPRB * ZCONV1 ! doesn't work, cf comment above
REAL(KIND=JPRB)            :: ZKOH_HIGH


! Reaction rate constant k_(S(iv)+H2O2(aq))   Seinfeld & Pandis, 1998, pp.366, 378, and 396
!   second one: see for instance Berglen et al 2004, JGR vol. 109, D19310
!   Notice the different units!
!   note: in S & P, there is a mistake on the unit of k, which is a 3rd order rate constant (not a 2nd order)
INTEGER(KIND=JPIM), PARAMETER :: IKH2O2_select = 1           ! Case selector, choose one of the following kinetics
REAL(KIND=JPRB), PARAMETER :: ZKH2O2_REF_v1 = 7.5E+7_JPRB    ! k_(S(iv)+H2O2) at Tref=298.15 K       [L**2 / mol**2 / s]
REAL(KIND=JPRB), PARAMETER :: ZKH2O2_TD_v1 = -4430._JPRB     ! temperature dependency                [K]
REAL(KIND=JPRB), PARAMETER :: ZKH2O2_REF_v2 = 8.0E+4_JPRB    ! k_(S(iv)+H2O2) at Tref=298.15 K       [L / mol / s]
REAL(KIND=JPRB), PARAMETER :: ZKH2O2_TD_v2 = -3650._JPRB     ! temperature dependency                [K]



! Reaction rate constants O3(aq) + S(iv)(aq) is split in 3 parts (reaction with SO2(aq), HSO3-(aq) and SO3=(aq)
REAL(KIND=JPRB), PARAMETER :: ZKO3_REF1 = 2.4E+4_JPRB    ! k_(SO2(aq))+O3(aq))                    [L / mol / s]
REAL(KIND=JPRB), PARAMETER :: ZKO3_REF2 = 3.7E+5_JPRB    ! k_(HSO3-(aq))+O3(aq)) at Tref=298.15 K [L / mol / s]
REAL(KIND=JPRB), PARAMETER :: ZKO3_TD2 = -5530._JPRB     ! temperature dependency                 [K]
REAL(KIND=JPRB), PARAMETER :: ZKO3_REF3 = 1.5E+9_JPRB    ! k_(SO3=(aq))+O3(aq)) at Tref=298.15 K  [L / mol / s]
REAL(KIND=JPRB), PARAMETER :: ZKO3_TD3 = -5280._JPRB     ! temperature dependency                 [K]


!==============!
! Molar masses !
!==============!
!REAL(KIND=JPRB), PARAMETER :: ZRMD    = 28.9644E-3_JPRB ! air  molar mass      [kg / mol]
REAL(KIND=JPRB), PARAMETER :: ZRMH2O2 = 34.0147E-3_JPRB ! H2O2 molar mass      [kg / mol]
REAL(KIND=JPRB), PARAMETER :: ZRMO3   = 47.9982E-3_JPRB ! O3   molar mass      [kg / mol]
REAL(KIND=JPRB), PARAMETER :: ZRMOH   = 17.008E-3_JPRB  ! OH   molar mass      [kg / mol]
REAL(KIND=JPRB), PARAMETER :: ZRMSO2  = 64.056E-3_JPRB  ! SO2  molar mass      [kg / mol] ! jjb quick fix: change unit in sucst.F90
REAL(KIND=JPRB), PARAMETER :: ZRMSO4  = 96.052E-3_JPRB  ! SO4= molar mass      [kg / mol]

! Other parameters
! ----------------
REAL(KIND=JPRB), PARAMETER :: ZRHOLW = 1000._JPRB ! liquid water density                          [kg / m3]

!REAL(KIND=JPRB), PARAMETER :: ZPH = 5._JPRB         ! pH of cloud liquid water, assumed
REAL(KIND=JPRB), PARAMETER :: ZHP = 1.0E-5_JPRB     ! proton concentration                        [mol / L]

REAL(KIND=JPRB), PARAMETER :: ZPTSCHEM = 2._JPRB   ! timestep for chemistry                              [s]



!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JL, JK    ! running indexes for latitude and level

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! jjb
!
REAL(KIND=JPRB) :: ZTEMPERAT ! grid cell temperature [K]
REAL(KIND=JPRB) :: ZPRESSURE ! grid cell pressure [Pa]

REAL(KIND=JPRB) :: ZAIR_CONC ! air concentration [mol / m3]
REAL(KIND=JPRB) :: ZAIR_DENS ! air density       [kg  / m3]

REAL(KIND=JPRB) :: ZCLW_VFRAC ! volume fraction of cloud liquid water [m3(aq) / m3(g)]

! conversion factor, ZCONV2 = R * T. Multiply by ZCONV2 to convert H^cp (in mol/m3/Pa) into H^cc (in m3(g) / m3(aq))
REAL(KIND=JPRB) :: ZCONV2   !                    [J/mol]

! Oxidants concentrations. Local arrays to split PCHEMSULF1

REAL(KIND=JPRB) :: ZFACT1   ! intermediate factor 1
REAL(KIND=JPRB) :: ZFACT2   ! intermediate factor 2

REAL(KIND=JPRB) :: ZFAQ_H2O2 ! fraction of H2O2 dissolved in aqueous phase  [dimensionless]
REAL(KIND=JPRB) :: ZFAQ_O3   ! fraction of O3   dissolved in aqueous phase  [dimensionless]
REAL(KIND=JPRB) :: ZFAQ_OH   ! fraction of OH   dissolved in aqueous phase  [dimensionless]
REAL(KIND=JPRB) :: ZFAQ_SO2  ! fraction of SO2  dissolved in aqueous phase  [dimensionless]

REAL(KIND=JPRB) :: ZHCC_H2O2    ! dimensionless Henry's law solubility for H2O2 [m3(g) / m3(aq)]
REAL(KIND=JPRB) :: ZHCC_O3      ! dimensionless Henry's law solubility for O3   [m3(g) / m3(aq)]
REAL(KIND=JPRB) :: ZHCC_OH      ! dimensionless Henry's law solubility for OH   [m3(g) / m3(aq)]
REAL(KIND=JPRB) :: ZHCC_SO2     ! dimensionless Henry's law solubility for SO2  [m3(g) / m3(aq)]
REAL(KIND=JPRB) :: ZHCC_SO2_EFF ! dimensionless Henry's law effective solubility for SO2 [m3(g) / m3(aq)]

REAL(KIND=JPRB) :: ZKEQ1_SO2 ! equilibrium constant SO2(aq) <=> HSO3-   [mol / L]
REAL(KIND=JPRB) :: ZKEQ2_SO2 ! equilibrium constant HSO3- <=> SO3=      [mol / L]
REAL(KIND=JPRB) :: ZKEQ1_FACT
REAL(KIND=JPRB) :: ZKEQ2_FACT

REAL(KIND=JPRB) :: ZKOH_LOW  ! low pressure limit = f(T)          [m**6 / (mol**2 * s)]
REAL(KIND=JPRB) :: ZKOH      ! gas phase reaction rate SO2 + OH         [m**3 / (mol * s)]
!REAL(KIND=JPRB) :: ZKPOH     ! modified gas phase reaction rate SO2 + OH  including OH and SO2 gas fractions  [m**3 / (mol * s)]

REAL(KIND=JPRB) :: ZKH2O2    ! Reaction rate constant k_(S(iv)+H2O2(aq))   [mol**2 / L**2 / s]

REAL(KIND=JPRB) :: ZKO3_1    ! Reaction rate constant k_(O3+SO2(aq)) (T)   [mol / m3 / s]
REAL(KIND=JPRB) :: ZKO3_2    ! Reaction rate constant k_(O3+HSO3-(aq)) (T) [mol / m3 / s]
REAL(KIND=JPRB) :: ZKO3_3    ! Reaction rate constant k_(O3+SO3--(aq)) (T) [mol / m3 / s]

REAL(KIND=JPRB) :: ZTFACT ! temperature factor for Henry's law solubility calculation [K**-1]

REAL(KIND=JPRB) :: ZC_H2O2_gas ! concentration of H2O2 in the gas phase  [mol/m3(air)]
REAL(KIND=JPRB) :: ZC_O3_gas   ! concentration of O3   in the gas phase  [mol/m3(air)]
REAL(KIND=JPRB) :: ZC_OH_gas   ! concentration of OH   in the gas phase  [mol/m3(air)]
REAL(KIND=JPRB) :: ZC_SO2_gas  ! concentration of SO2  in the gas phase  [mol/m3(air)]
REAL(KIND=JPRB) :: ZC_SO2_tot  ! total concentration of SO2              [mol/m3(air)]

REAL(KIND=JPRB) :: ZC_H2O2_gas_ini ! initial concentration of H2O2 in the gas phase  [mol/m3(air)]
REAL(KIND=JPRB) :: ZC_O3_gas_ini   ! initial concentration of O3   in the gas phase  [mol/m3(air)]
REAL(KIND=JPRB) :: ZC_OH_gas_ini   ! initial concentration of OH   in the gas phase  [mol/m3(air)]

REAL(KIND=JPRB) :: ZC_H2O2_aqp  ! "potential" concentration of H2O2  in the aqueous phase  [mol/m3(aq)]
REAL(KIND=JPRB) :: ZC_O3_aqp    ! "potential" concentration of O3    in the aqueous phase  [mol/m3(aq)]
REAL(KIND=JPRB) :: ZC_Siv_aqp   ! "potential" concentration of S(iv) in the aqueous phase  [mol/m3(aq)]
REAL(KIND=JPRB) :: ZC_SO2_aqp   ! "potential" concentration of SO2   in the aqueous phase  [mol/m3(aq)]
REAL(KIND=JPRB) :: ZC_HSO3m_aqp ! "potential" concentration of HSO3- in the aqueous phase  [mol/m3(aq)]
REAL(KIND=JPRB) :: ZC_SO3mm_aqp ! "potential" concentration of SO3=  in the aqueous phase  [mol/m3(aq)]

REAL(KIND=JPRB) :: ZTend_OH     ! tendency for S(iv) + OH(g)                  [mol/m3(air)]
REAL(KIND=JPRB) :: ZTend_H2O2   ! tendency for S(iv) + H2O2(aq)               [mol/m3(air)]  <== final unit
REAL(KIND=JPRB) :: ZTend_O3     ! tendency for S(iv) + O3(aq)                 [mol/m3(air)]  <== final unit
REAL(KIND=JPRB) :: ZTend_O3_r1  ! tendency for SO2(aq) + O3(aq)                 [mol/m3(air)]  <== final unit
REAL(KIND=JPRB) :: ZTend_O3_r2  ! tendency for HS03m + O3(aq)                 [mol/m3(air)]  <== final unit
REAL(KIND=JPRB) :: ZTend_O3_r3  ! tendency for SO3mm + O3(aq)                 [mol/m3(air)]  <== final unit
REAL(KIND=JPRB) :: ZTend_O3_ul  ! unlimited
REAL(KIND=JPRB) :: ZTend_O3_li  ! limited
REAL(KIND=JPRB) :: ZLimit_fact_O3
REAL(KIND=JPRB) :: ZSum_HSO3m_ox_ul
REAL(KIND=JPRB) :: ZSum_HSO3m_ox_li
REAL(KIND=JPRB) :: ZLimit_fact_HSO3m
!REAL(KIND=JPRB) :: ZTend_TOT_1  ! tendency for S(iv) --> S(vi) including oxidants and "branch-specific" S(iv) limitation
!REAL(KIND=JPRB) :: ZTend_TOT_2  ! tendency for S(iv) --> S(vi) further including total S(iv) limitation
REAL(KIND=JPRB) :: ZTend_H2O2_Sum  !
REAL(KIND=JPRB) :: ZTend_O3_Sum  !
REAL(KIND=JPRB) :: ZTend_OH_Sum  !
REAL(KIND=JPRB) :: ZTend_Gas  !
REAL(KIND=JPRB) :: ZTend_Aq  !
REAL(KIND=JPRB) :: ZTend_Aq_Sum  !
REAL(KIND=JPRB) :: ZTend_Sum  !

REAL(KIND=JPRB) :: ZPNEB, ZSCALEO3(KLON), ZSCALEOH(KLON)

REAL(KIND=JPRB) :: ZmmrH2O2
REAL(KIND=JPRB) :: ZmmrO3
REAL(KIND=JPRB) :: ZmmrOH
REAL(KIND=JPRB) :: ZmmrSO2

INTEGER(KIND=JPIM) :: JPTS ! Number of sub-timesteps for chemistry
INTEGER(KIND=JPIM) :: JTS  ! Running index for chemistry do-loop

#include "compo_diurnal.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_SO2SO4_V2',0,ZHOOK_HANDLE)

! Initialisation of "parameters" that cannot be initialised in the declaration part (cf. remarks above)
ZCONV1 = 1.E-6_JPRB * RNAVO                  ! [m3 / cm3 / mol]
ZKOH_LOW300 = 3.3E-31_JPRB * ZCONV1 * ZCONV1 ! now in [m**6 / mol**2 / s]
ZKOH_HIGH   = 1.6E-12_JPRB * ZCONV1          ! now in [m**3 / mol    / s]

! Define the number of intermediate timesteps for chemistry
IF(MODULO(PTSPHY,ZPTSCHEM) > 0._JPRB) THEN
   CALL ABOR1('ABORT: IN AER_SO2SO4_V2, ZPTSCHEM must divide PTSTEMP')
ENDIF
JPTS = INT(PTSPHY/ZPTSCHEM)


! Split PCHEMSULF1 array

!! Initialise tendencies (necessary when using sub-timesteps for chemistry,
!!  since the total tendency will be the sum of "sub-tendencies"
PTSO2(:,:) = 0._JPRB
PTSO4(:,:) = 0._JPRB

CALL COMPO_DIURNAL(YDRIP, KIDIA, KFDIA, KLON, 'Sine', PGELAM, PGELAT, ZSCALEOH, PAMPLITUDE=0.7_JPRB, PHOURPEAK=15.0_JPRB)
CALL COMPO_DIURNAL(YDRIP, KIDIA, KFDIA, KLON, 'Sine', PGELAM, PGELAT, ZSCALEO3, PAMPLITUDE=0.7_JPRB, PHOURPEAK=15.0_JPRB)

DO JK=1,KLEV
   DO JL=KIDIA,KFDIA
      ! Initialise all scalar variables to 0._JPRB (for debug only)
      ZmmrH2O2 = 0._JPRB
      ZmmrO3   = 0._JPRB
      ZC_H2O2_gas = 0._JPRB
      ZC_O3_gas   = 0._JPRB
      ZFAQ_H2O2 = 0._JPRB
      ZFAQ_O3 = 0._JPRB
      ZFAQ_OH = 0._JPRB
      ZFAQ_SO2 = 0._JPRB
      ZHCC_H2O2 = 0._JPRB
      ZHCC_OH = 0._JPRB
      ZHCC_O3 = 0._JPRB
      ZHCC_SO2 = 0._JPRB
      ZHCC_SO2_EFF = 0._JPRB
      ZTend_OH_Sum = 0._JPRB
      ZTend_Gas = 0._JPRB
      ZTend_Aq = 0._JPRB
      ZTend_O3_Sum =  0._JPRB
      ZTend_H2O2_Sum =  0._JPRB

! -- Miscellaneous variables
      ! -- pressure and temperature copied in scalars
      ZPRESSURE = PRSF1(JL,JK)
      ZTEMPERAT = PTP  (JL,JK)
      ! -- air concentration and density
      ZAIR_CONC = ZPRESSURE / (R *ZTEMPERAT) ! [mol / m3]
      ZAIR_DENS = ZPRESSURE / (RD*ZTEMPERAT) ! [kg  / m3]
      ! -- fractional cloudiness
      ZPNEB = PNEB(JL,JK)
      ! -- liquid volume fraction in the cloudy part [m3(aq) / m3(g)]

      ZCLW_VFRAC = 0._JPRB
      IF (ZPNEB > 1.0E-12_JPRB) THEN
        ZCLW_VFRAC = PQLI(JL,JK) / ZPNEB * ZAIR_DENS / ZRHOLW
      ENDIF


      !write (*,*) "SO2SO4",JL,JK,ZPNEB,ZCLW_VFRAC,POH(JL,JK),PO3(JL,JK),PH2O2(JL,JK)
!=======================================================================
! Two cases:
!   - gas phase reaction only (if no cloud, or no liquid water)
!   - gas + aqueous phase reactions
!=======================================================================
      IF (.NOT.(ZPNEB > 1.0E-12_JPRB )) THEN  ! 1.0E-12 is the default PNEB value

!=======================================================================
!    GAS PHASE ONLY:
!=======================================================================

! -- Reactants concentration (gas phase only)
         ! -- mmr [kg / kg]
         IF (POH(JL,JK) > 1.0E-14_JPRB ) THEN  
           ZmmrOH   = POH(JL,JK)*ZSCALEOH(JL)
         ELSE
           ZmmrOH   = 0.0_JPRB
         ENDIF
         ZmmrSO2  = PSO2(JL,JK) + PTSPHY * PITSO2(JL,JK)
         ! -- gas phase concentration [mol / m3]
         ZC_OH_gas   = ZmmrOH   * ZAIR_DENS / ZRMOH   ! OH concentration from climatology is a gas concentration
         ZC_SO2_tot  = ZmmrSO2  * ZAIR_DENS / ZRMSO2  ! SO2 concentration : in this case, ZC_SO2_gas = ZC_SO2_tot
         ZC_SO2_gas  = ZC_SO2_tot

         IF (ZC_OH_gas < 1E-18_JPRB) THEN
           ZC_OH_gas=0._JPRB
         ENDIF
         IF (ZC_SO2_gas < 1E-18_JPRB) THEN
           ZC_SO2_gas=0._JPRB
         ENDIF
         IF (ZC_SO2_tot < 1E-18_JPRB) THEN
           ZC_SO2_tot=0._JPRB
         ENDIF


         ZC_OH_gas_ini = ZC_OH_gas
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

! -- reaction rate constants
         ! -- k_OH reaction rate [m3 / mol / s]
         !         reference: JPL Pub. 15-10, Eval Nb 18

         ! low pressure limit depends on temperature
         ZKOH_LOW = ZKOH_LOW300 * (300._JPRB / ZTEMPERAT)**4.3_JPRB

         ! complete expression: see doc
         ZFACT1 = ZKOH_LOW * ZAIR_CONC
         ZFACT2 = ZFACT1 / ZKOH_HIGH

         ZKOH = ZFACT1 / (1._JPRB + ZFACT2)
         ZKOH = ZKOH * .6_JPRB ** (1._JPRB / (1._JPRB + LOG10(ZFACT2)**2) )
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
         IF (ZKOH < 1E-18_JPRB) THEN
           ZKOH=0._JPRB
         ENDIF


         DO JTS=1,JPTS

            ! -- ZTend_OH [mol / m3(air)]
            ZTend_OH = ZKOH * ZC_OH_gas * ZC_SO2_gas * ZPTSCHEM
            ZTend_OH = MIN(ZTend_OH, ZC_OH_gas, ZC_SO2_gas)     ! in this case, ZC_SO2_tot = ZC_SO2_gas

            ! -- ZTend_OH_Sum [mol / m3(air)]
            ZTend_OH_Sum = ZTend_OH_Sum + ZTend_OH

            ! -- Update gas phase concentrations
            ZC_OH_gas  = ZC_OH_gas  - ZTend_OH ! in this test version, no depletion of OH during oxidation
            ZC_SO2_gas = ZC_SO2_gas - ZTend_OH

         ENDDO

         ! -- PT(xxx) in [kg(xxx) / kg(air) / s]
         !    using PTSPHY=sum(ZPTSCHEM)
         PTSO2(JL,JK) = PTSO2(JL,JK) - ZTend_OH_Sum * ZRMSO2 / ZAIR_DENS / PTSPHY
         PTSO4(JL,JK) =  PTSO4(JL,JK)+ ZTend_OH_Sum * ZRMSO4 / ZAIR_DENS / PTSPHY


! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --


      ELSE
!=======================================================================
!    GAS + AQUEOUS PHASE REACTION:
!=======================================================================

! -- Miscellaneous variables only used in aqueous phase
         ! -- temperature factor for Henry's law and reaction rates
         ZTFACT = 1._JPRB / ZTEMPERAT - ZITREF_H
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --


! -- dimensionless Henry's law solubility constants [m3(g) / m3(aq)]
         ! -- Conversion factor, to get dimensionless Henry coefficients (H^cc)
         ZCONV2 = R * ZTEMPERAT

         ! -- of H2O2
         ZHCC_H2O2 = ZHCP_H2O2_REF * EXP(ZHCP_H2O2_TD * ZTFACT) * ZCONV2
         ! -- of O3
         ZHCC_O3 = ZHCP_O3_REF * EXP(ZHCP_O3_TD * ZTFACT) * ZCONV2
         ! -- of OH
         ZHCC_OH = ZHCP_OH_REF * ZCONV2                           ! no temperature dependency for K^cp(OH)

         ! -- of SO2
         ZHCC_SO2 = ZHCP_SO2_REF * EXP(ZHCP_SO2_TD * ZTFACT) * ZCONV2
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

! -- acid-base equilibrium constants for S(iv) aqueous species [mol / L]
         ! NB: non SI units are kept here, since they are immediately converted into a dimensionless variable...
         ZKEQ1_SO2 = ZKEQ1_SO2_REF * EXP(ZKEQ1_SO2_TD * ZTFACT)
         ZKEQ2_SO2 = ZKEQ2_SO2_REF * EXP(ZKEQ2_SO2_TD * ZTFACT)
         !     ... since ZHP is also in [mol / L]  ==> ZKEQ1/2 are dimensionless
         ZKEQ1_FACT = ZKEQ1_SO2 / ZHP
         ZKEQ2_FACT = ZKEQ2_SO2 / ZHP
         ! Effective Henry solubility coefficient for SO2
         ZHCC_SO2_EFF = ZHCC_SO2 * ( 1._JPRB + ZKEQ1_FACT + ZKEQ1_FACT * ZKEQ2_FACT )
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

! -- aqueous fractions in the cloudy part [dimensionless]
         ZFAQ_H2O2 = 1._JPRB / (1._JPRB + 1._JPRB / (ZHCC_H2O2    * ZPNEB ) )
         ZFAQ_O3   = 1._JPRB / (1._JPRB + 1._JPRB / (ZHCC_O3      * ZPNEB ) )
         ZFAQ_OH   = 1._JPRB / (1._JPRB + 1._JPRB / (ZHCC_OH      * ZPNEB ) )
         ZFAQ_SO2  = 1._JPRB / (1._JPRB + 1._JPRB / (ZHCC_SO2_EFF * ZPNEB ) )
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --


! -- Oxidants concentrations / partial pressure / ...
         ! -- mmr [kg / kg]
         IF (PH2O2(JL,JK) > 1.0E-14_JPRB ) THEN  
           ZmmrH2O2   = PH2O2(JL,JK)
         ELSE
           ZmmrH2O2   = 0.0_JPRB
         ENDIF
         IF (PO3(JL,JK) > 1.0E-14_JPRB ) THEN  
           ZmmrO3   = PO3(JL,JK)*ZSCALEO3(JL)
         ELSE
           ZmmrO3   = 0.0_JPRB
         ENDIF
         IF (POH(JL,JK) > 1.0E-14_JPRB ) THEN  
           ZmmrOH   = POH(JL,JK)*ZSCALEOH(JL)
         ELSE
           ZmmrOH   = 0.0_JPRB
         ENDIF
         ZmmrSO2  = PSO2(JL,JK) + PTSPHY * PITSO2(JL,JK)
    
! -- gas phase concentration [mol / m3(g)]
         ! -- oxidants concentrations from climatologies are gas phase concentration.
         !    use aqueous fractions calculated above to get total concentrations in the grid
         ZC_H2O2_gas = ZmmrH2O2 * ZAIR_DENS / ZRMH2O2
         ZC_O3_gas   = ZmmrO3 * ZAIR_DENS / ZRMO3
         ZC_OH_gas   = ZmmrOH * ZAIR_DENS / ZRMOH

         ZC_H2O2_gas_ini = ZC_H2O2_gas
         ZC_O3_gas_ini   = ZC_O3_gas
         ZC_OH_gas_ini   = ZC_OH_gas


       ! -- S(iv) concentrations from climatology is total concentration
         ZC_SO2_tot  = ZmmrSO2  * ZAIR_DENS / ZRMSO2
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --


! -- reaction rate constants
         ! -- k_OH reaction rate [m3 / mol / s]
         !         reference: JPL Pub. 15-10, Eval Nb 18

         ! low pressure limit depends on temperature
         ZKOH_LOW = ZKOH_LOW300 * (300._JPRB / ZTEMPERAT)**4.3_JPRB

         ! complete expression: see doc
         ZFACT1 = ZKOH_LOW * ZAIR_CONC
         ZFACT2 = ZFACT1 / ZKOH_HIGH

         ZKOH = ZFACT1 / (1._JPRB + ZFACT2)
         ZKOH = ZKOH * .6_JPRB ** (1._JPRB / (1._JPRB + LOG10(ZFACT2)**2) )

         ! -- k_H2O2 reaction rate
         !           reference: Seinfeld and Pandis, 2nd Ed., 2006, table 7.6 p.316
         SELECT CASE(IKH2O2_select)
         CASE(1)
            ZKH2O2 = ZKH2O2_REF_v1 * EXP(ZKH2O2_TD_v1 * ZTFACT)      ! unit is [L**2 / mol**2 / s]
            ZKH2O2 = ZKH2O2 * ZHP / (1._JPRB + 13._JPRB * ZHP) ! unit is [L    / mol    / s]
            ZKH2O2 = ZKH2O2 * ZCONV3                           ! unit is [m**3 / mol    / s]
         CASE(2)
            ZKH2O2 = ZKH2O2_REF_v2 * EXP(ZKH2O2_TD_v2 * ZTFACT)      ! unit is [L / mol / s] ???
            ZKH2O2 = ZKH2O2 / (0.1_JPRB + ZHP)                 ! unit is [L / mol / s] ??? according to Tsai et al ACP 2010
            ZKH2O2 = ZKH2O2 * ZCONV3                           ! unit is [m**3 / mol / s]
         CASE DEFAULT
            CALL ABOR1('ABORT: IN AER_SO2SO4_V2, IKH2O2_select has incorrect value')
         END SELECT

         ! -- k_O3 reaction rate constants
         !         reference: Seinfeld and Pandis, 2nd Ed., 2006, table 7.A.7 p.329 (in [L / mol / s])
         !                    converted in [m3 / mol / s] with ZCONV3
         ZKO3_1 = ZKO3_REF1                          * ZCONV3 ! no temperature dependency
         ZKO3_2 = ZKO3_REF2 * EXP(ZKO3_TD2 * ZTFACT) * ZCONV3
         ZKO3_3 = ZKO3_REF3 * EXP(ZKO3_TD3 * ZTFACT) * ZCONV3
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --


         DO JTS=1,JPTS


! -- gas phase concentration [mol / m3(g)]
            ZC_SO2_gas  = ZC_SO2_tot * (1._JPRB - ZFAQ_SO2*ZCLW_VFRAC)

! -- aqueous phase concentration [mol / m3(aq)]    ! WARNING: these are "potential" aqueous concentrations
            ZC_H2O2_aqp = ZC_H2O2_gas * ZHCC_H2O2
            ZC_O3_aqp = ZC_O3_gas * ZHCC_O3

            ZC_Siv_aqp = ZC_SO2_gas * ZHCC_SO2_EFF

            ! Update fractioning between aqueous forms of S(iv)
            ZC_SO2_aqp = ZC_Siv_aqp / (1._JPRB + ZKEQ1_FACT + ZKEQ1_FACT * ZKEQ2_FACT )
            ZC_HSO3m_aqp = ZC_SO2_aqp * ZKEQ1_FACT
            ZC_SO3mm_aqp = ZC_HSO3m_aqp * ZKEQ2_FACT
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

! -- tendencies = reaction rates * timestep
            ! -- ZTend_OH [mol / m3(air)]
            ZTend_OH = ZKOH * ZC_OH_gas * ZC_SO2_gas * ZPTSCHEM
            ZTend_OH = MIN(ZTend_OH, ZC_OH_gas, ZC_SO2_gas)

            ! -- ZTend_OH_Sum [mol / m3(air)]
            ZTend_OH_Sum = ZTend_OH_Sum + ZTend_OH

            ! -- Update gas phase concentrations
            !ZC_OH_gas  = ZC_OH_gas  - ZTend_OH ! in this test version, no depletion of OH during oxidation
            !ZC_SO2_gas = ZC_SO2_gas - ZTend_OH ! Will be updated using ZC_SO2_tot


! WARNING: the aqueous tendencies are calculated using "potential" aqueous concentrations (ZC_xxx_aqp) (in [mol / m3(aq)] )
!      Then, the actual aqueous quantities need to be evaluated so that the actual oxidation is limited by the aqueous quantities,
!      which are calculated using ZC_xxx_tot * ZFAQ_xxx * ZPNEB

            ! -- Tend_H2O2
            ZTend_H2O2 = ZKH2O2 * ZC_H2O2_aqp * ZC_HSO3m_aqp * ZPTSCHEM  ! unit is [mol / m**3(aq)]

            ! note that when ~2 < pH < ~7, HSO3- is the dominant S(iv) aqueous species
            !           when ~3 < pH < ~6, HSO3- is over 90% of S(iv) aqueous species
            !    ref: Seinfeld & Pandis, 1st edition, 1998, Fig. 6.8 p. 352
            ZTend_H2O2 = MIN(ZTend_H2O2, ZC_H2O2_aqp, ZC_HSO3m_aqp) ! limited by reactant availability ("potential" concentrations)
            !ZTend_H2O2 = ZTend_H2O2 * ZCLW_VFRAC * ZPNEB             ! convert to [mol / m**3(g)]

            ! -- Tend_O3
            ZTend_O3_r1 = ZKO3_1 * ZC_SO2_aqp * ZC_O3_aqp * ZPTSCHEM   ! unit is [mol / m**3(aq)]
            ZTend_O3_r2 = ZKO3_2 * ZC_HSO3m_aqp * ZC_O3_aqp * ZPTSCHEM ! unit is [mol / m**3(aq)]
            ZTend_O3_r3 = ZKO3_3 * ZC_SO3mm_aqp * ZC_O3_aqp * ZPTSCHEM ! unit is [mol / m**3(aq)]

            ZTend_O3_r1 = MIN(ZTend_O3_r1, ZC_SO2_aqp, ZC_O3_aqp)   ! limited by reactants availability
            ZTend_O3_r2 = MIN(ZTend_O3_r2, ZC_HSO3m_aqp, ZC_O3_aqp) ! limited by reactants availability
            ZTend_O3_r3 = MIN(ZTend_O3_r3, ZC_SO3mm_aqp, ZC_O3_aqp) ! limited by reactants availability

            ! Limiting factor due to O3 availability (in the aqueous phase)
            !   --> will be applied to the sum of 3 reaction paths with O3
            ZTend_O3_ul = ZTend_O3_r1 + ZTend_O3_r2 + ZTend_O3_r3
            ZTend_O3_li = MIN(ZTend_O3_ul, ZC_O3_aqp)
            IF(ZTend_O3_ul > 0._JPRB) THEN
               ZLimit_fact_O3 = ZTend_O3_li / ZTend_O3_ul
            ELSE
               ZLimit_fact_O3 = 0._JPRB
            ENDIF

            ! Limiting factor due to HSO3- availability
            !   --> will be applied to the 2 reaction with HSO3-
            ZSum_HSO3m_ox_ul = ZTend_H2O2 + ZTend_O3_r2 * ZLimit_fact_O3
            ZSum_HSO3m_ox_li = MIN(ZSum_HSO3m_ox_ul, ZC_HSO3m_aqp)
            IF(ZSum_HSO3m_ox_ul > 0._JPRB) THEN
               ZLimit_fact_HSO3m = ZSum_HSO3m_ox_li / ZSum_HSO3m_ox_ul
            ELSE
               ZLimit_fact_HSO3m = 0._JPRB
            ENDIF

            ! Final tendencies expressed as "potential" concentrations [mol / m3(aq)]
            ZTend_O3 = (ZTend_O3_r1 + ZTend_O3_r2*ZLimit_fact_HSO3m + ZTend_O3_r3) * ZLimit_fact_O3
            ZTend_H2O2 = ZTend_H2O2*ZLimit_fact_HSO3m

            ZTend_O3_Sum = ZTend_O3_Sum + ZTend_O3       ! These variables record the fraction of each oxidation path
            ZTend_H2O2_Sum = ZTend_H2O2_Sum + ZTend_H2O2 !    --   (here expressed in  [mol / m3(aq)] )

            ZTend_Aq = ZTend_O3 + ZTend_H2O2

            ! Update aqueous "potential" concentrations
            !ZC_H2O2_aqp = ZC_H2O2_aqp - ZTend_H2O2 ! in this test version, no depletion of OH during oxidation
            !ZC_O3_aqp = ZC_O3_aqp - ZTend_O3       ! in this test version, no depletion of OH during oxidation

            ZC_Siv_aqp = ZC_Siv_aqp - ZTend_Aq

            ! -- Update concentrations that will be used at the beginning of loop
            ZTend_Aq = ZTend_Aq * ZCLW_VFRAC ! convert to [mol / m**3(g)] for the last calculations

            ! S(iv)tot : initial - gas - aq
            ZC_SO2_tot = ZC_SO2_tot - ZTend_OH - ZTend_Aq
            ! H2O2 and O3 : initial - aq
            ZC_H2O2_gas = ZC_H2O2_gas - ZTend_Aq
            ZC_O3_gas = ZC_O3_gas - ZTend_Aq
            
            
         ENDDO

         ZTend_H2O2_Sum = ZTend_H2O2_Sum * ZCLW_VFRAC ! convert to [mol / m**3(g)]
         ZTend_O3_Sum = ZTend_O3_Sum * ZCLW_VFRAC     ! convert to [mol / m**3(g)]
         ZTend_Aq_Sum = ZTend_H2O2_Sum + ZTend_O3_Sum

         ZTend_Sum = ZTend_Aq_Sum + ZTend_OH_Sum

         PTSO2(JL,JK) = PTSO2(JL,JK) -ZTend_Sum * ZRMSO2 / ZAIR_DENS / PTSPHY
         PTSO4(JL,JK) =  PTSO4(JL,JK) + ZTend_Sum * ZRMSO4 / ZAIR_DENS / PTSPHY

! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

      ENDIF


   ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  PFSO4(JL)=0.0_JPRB
  PFSO2(JL)=0.0_JPRB
ENDDO


DO JK=1,KLEV
   DO JL=KIDIA,KFDIA
      PFSO4(JL) = PFSO4(JL) + PTSO4(JL,JK)*(PDP(JL,JK))/RG
      PFSO2(JL) = PFSO2(JL) + PTSO2(JL,JK)*(PDP(JL,JK))/RG
   ENDDO
ENDDO


!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_SO2SO4_V2',1,ZHOOK_HANDLE)
END SUBROUTINE AER_SO2SO4_V2
