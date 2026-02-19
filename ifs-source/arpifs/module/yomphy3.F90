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

MODULE YOMPHY3

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

! number of spectral bands used in cloud model
INTEGER(KIND=JPIM), PARAMETER :: N_SPBANDL=2

TYPE TPHY3
!*
!     ------------------------------------------------------------------
!     CONSTANTES PHYSIQUES REGLABLES UTILISEES POUR LES CALCULS
!     RADIATIFS :
!       BSFSA   : "BACK-SCATTERED FRACTION" SOLAIRE POUR LES AEROSOLS.
!               : SOLAR "BACK-SCATTERED FRACTION" FOR AEROSOLS.
!       BSFSI   : "BACK-SCATTERED FRACTION" SOLAIRE POUR LA GLACE.
!               : SOLAR "BACK-SCATTERED FRACTION" FOR ICE CLOUDS.
!       BSFSN   : "BACK-SCATTERED FRACTION" SOLAIRE POUR LES NUAGES.
!               : SOLAR "BACK-SCATTERED FRACTION" FOR CLOUDS.
!       BSFTA   : "BACK-SCATTERED FRACTION" THERMIQUE POUR LES AEROSOLS.
!               : THERMAL "BACK-SCATTERED FRACTION" FOR AEROSOLS.
!       BSFTI   : "BACK-SCATTERED FRACTION" THERMIQUE POUR LA GLACE.
!               : THERMAL "BACK-SCATTERED FRACTION" FOR ICE CLOUDS.
!       BSFTN   : "BACK-SCATTERED FRACTION" THERMIQUE POUR LES NUAGES.
!               : THERMAL "BACK-SCATTERED FRACTION" FOR CLOUDS.
!       EARRT   : EPAISSEUR DE L'ATMOSPHERE  /  RAYON DE LA TERRE.
!               : RATIO "DEPTH OF THE ATMOSPHERE / EARTH'S RADIUS".
!       EOASA   : COEFFICIENT D'ABSORPTION SOLAIRE PAR LES AEROSOLS.
!               : SOLAR ABSORPTION COEFFICIENT FOR AEROSOLS.
!       EOASI   : COEFFICIENT D'ABSORPTION SOLAIRE PAR LA GLACE.
!               : SOLAR ABSORPTION COEFFICIENT FOR ICE CLOUDS.
!       EOASN   : COEFFICIENT D'ABSORPTION SOLAIRE PAR LES NUAGES.
!               : SOLAR ABSORPTION COEFFICIENT FOR CLOUDS.
!       EOATA   : COEFFICIENT D'ABSORPTION THERMIQUE PAR LES AEROSOLS.
!               : THERMAL ABSORPTION COEFFICIENT FOR AEROSOLS.
!       EOATI   : COEFFICIENT D'ABSORPTION THERMIQUE PAR LA GLACE.
!               : THERMAL ABSORPTION COEFFICIENT FOR ICE CLOUDS.
!       EOATN   : COEFFICIENT D'ABSORPTION THERMIQUE PAR LES NUAGES.
!               : THERMAL ABSORPTION COEFFICIENT FOR CLOUDS.
!       EODSA   : COEFFICIENT DE DIFFUSION SOLAIRE PAR LES AEROSOLS.
!               : SOLAR SCATTERING COEFFICIENT FOR AEROSOLS.
!       EODSI   : COEFFICIENT DE DIFFUSION SOLAIRE PAR LA GLACE.
!               : SOLAR SCATTERING COEFFICIENT FOR ICE CLOUDS.
!       EODSN   : COEFFICIENT DE DIFFUSION SOLAIRE PAR LES NUAGES.
!               : SOLAR SCATTERING COEFFICIENT FOR CLOUDS.
!       EODTA   : COEFFICIENT DE DIFFUSION THERMIQUE PAR LES AEROSOLS.
!               : THERMAL SCATTERING COEFFICIENT FOR AEROSOLS.
!       EODTI   : COEFFICIENT DE DIFFUSION THERMIQUE PAR LA GLACE.
!               : THERMAL SCATTERING COEFFICIENT FOR ICE CLOUDS.
!       EODTN   : COEFFICIENT DE DIFFUSION THERMIQUE PAR LES NUAGES.
!               : THERMAL SCATTERING COEFFICIENT FOR CLOUDS.
!       EORAY   : COEFFICIENT DE DIFFUSION RAYLEIGH.
!               : RAYLEIGH SCATTERING COEFFICIENT.
!       GCA(6)  : POUR LE CALCUL "WEAK LINE" DE LA LARGEUR EQUIVALENTE.
!               : FOR THE "WEAK LINE" PART OF THE EQUIVALENT WIDTH.
!       GCB(6)  : POUR LE CALCUL "STRONG LINE" DE LA LARGEUR EQUIVAL..
!               : FOR THE "STRONG LINE" PART OF THE EQUIVALENT WIDTH.
!       GCC(6)  : POUR LE CALCUL "CONTINUUM" DE LA LARGEUR EQUIVALENTE.
!               : FOR THE "CONTINUUM" PART OF THE EQUIVALENT WIDTH.
!       GCD4    : POUR LA CONTRIBUTION "E-TYPE" A GCC(4) (H2O THERM.).
!               : FOR THE E-TYPE CONTRIBUTION TO GCC(4) (H2O THERM.).
!       GCE4    : POUR LA DEPENDANCE EN TEMPERATURE DU "E-TYPE" (GCD4).
!               : FOR THE TEMPERATURE DEPENDENCY OF THE E-TYPE (GCD4).
!       GIREC*  : JEU DE COEFFICIENTS MODULANT L'INTERACTION INFRA-ROUGE ENTRE COUCHES.
!               : COEFFICIENTS SET TO TUNE THE INFRA-RED EXCHANGE BETWEEN LAYERS.
!       QCO2    : CONCENTRATION MASSIQUE DU CO2.
!               : SPECIFIC RATIO OF CO2.
!       RII0    : VALEUR INSTANTANNEE DE LA CONST. SOLAIRE (CYCLE ANN.).
!               : INSTANTANEOUS VALUE OF THE SOLAR CONST. (ANN. CYCLE).
!       USAA    : AU NUMERATEUR DE "L'UPSCATTERED FRACTION" CAS AEROS.
!               : AT THE UPPER CASE OF THE UPSCATTERED FRACTION, AEROS.
!       USAI    : AU NUMERATEUR DE "L'UPSCATTERED FRACTION" CAS GLACE.
!               : AT THE UPPER CASE OF THE UPSCATTERED FRACTION, ICE.
!       USAN    : AU NUMERATEUR DE "L'UPSCATTERED FRACTION" CAS NUAGES.
!               : AT THE UPPER CASE OF THE UPSCATTERED FRACTION, CLOUDS.
!       USBA    : AU DENOMINATEUR DE "L'UPSCATTERED FRACT." CAS AEROS.
!               : AT THE LOWER CASE OF THE UPSCATTERED FRACTION, AEROS.
!       USBI    : AU DENOMINATEUR DE "L'UPSCATTERED FRACT." CAS GLACE.
!               : AT THE LOWER CASE OF THE UPSCATTERED FRACTION, ICE.
!       USBN    : AU DENOMINATEUR DE "L'UPSCATTERED FRACT." CAS NUAGES.
!               : AT THE LOWER CASE OF THE UPSCATTERED FRACTION, CLOUDS.
!       VDP(5,6): AU DENOMINATEUR DES FONCTIONS DE PADE POUR LES GAZ.
!               : AT THE LOWER CASE OF PADE FUNCTIONS FOR GASES.
!       VNP(5,6): AU NUMERATEUR DES FONCTIONS DE PADE POUR LES GAZ.
!               : AT THE UPPER CASE OF PADE FUNCTIONS FOR GASES.
!       GOLC*   : PARAMETERS FOR GASEOUS OVERLAPS (OLD ACRANEB).
! Parameters for cloud model:

!   Notations:
!     g      - asymmetry factor            (unscaled)
!     k_abs  - mass absorption coefficient (delta-scaled)
!     k_scat - mass scattering coefficient (delta-scaled)
!     delta0 - unsaturated optical depth
!     c_abs  - saturation factor for k_abs
!     c_scat - saturation factor for k_scat
!     iwc    - ice water content
!     lwc    - liquid water content

!   First index of FCM arrays (FCM = Fitting parameters for Cloud Model)
!   denotes spectral band:
!     1      - solar
!     2      - thermal

! exclusively ACRANEB (NRAY=1):
!   REXP_NEB          : Scaling exponent for cloud fraction in definition
!                       of effective delta0.
!   FCM_DEL_A(2)      : Critical value of delta0 for computation of c_abs.
!   FCM_DEL_D(2)      : Critical value of delta0 for computation of c_scat.
!   FCM_MU_A(2)       : Exponent mu for computation of c_abs. 
!   FCM_MU_D(2)       : Exponent mu for computation of c_scat.
!   FCM_N_I           : Scaling exponent for iwc.
!   FCM_N_L           : Scaling exponent for lwc.

! exclusively ACRANEB2 (NRAY=2):
!   FCM_AI            : Tuning parameter for gas-ice cloud spectral overlap.
!   FCM_AL            : Tuning parameter for gas-liquid cloud spectral overlap.
!   FCM_B_AI          : Tuning parameter for saturation of ice clouds.
!   FCM_B_AL          : Tuning parameter for saturation of liquid clouds.
!   FCM_B_BI          : Tuning parameter for saturation of ice clouds.
!   FCM_B_BL          : Tuning parameter for saturation of liquid clouds.
!   FCM_DEL_AI        : Critical value of delta0 for computation of c_abs_i.
!   FCM_DEL_AL        : Critical value of delta0 for computation of c_abs_l.
!   FCM_DEL_DI        : Critical value of delta0 for computation of c_scat_i.
!   FCM_DEL_DL        : Critical value of delta0 for computation of c_scat_l.
!   FCM_MU_AI         : Exponent mu for computation of c_abs_i.
!   FCM_MU_AL         : Exponent mu for computation of c_abs_l.
!   FCM_MU_DI         : Exponent mu for computation of c_scat_i.
!   FCM_MU_DL         : Exponent mu for computation of c_scat_l.
!   FCM_NU_AI         : Exponent nu for computation of c_abs_i.
!   FCM_NU_AL         : Exponent nu for computation of c_abs_l.
!   FCM_NU_DI         : Exponent nu for computation of c_scat_i.
!   FCM_NU_DL         : Exponent nu for computation of c_scat_l.
!   FCM_IWC2DE(2,-2:3): Parameters for IWC to D_e conversion.
!   FCM_LWC2RE(2,-2:3): Parameters for LWC to R_e conversion.

! both ACRANEB and ACRANEB2:
!   FCM_P_AI(2,0:3)   : Pade coefficients in numerator for k_abs, ice.
!   FCM_P_AL(2,0:3)   : Pade coefficients in numerator for k_abs, liquid.
!   FCM_P_DI(2,0:3)   : Pade coefficients in numerator for k_scat, ice.
!   FCM_P_DL(2,0:3)   : Pade coefficients in numerator for k_scat, liquid.
!   FCM_P_GI(2,0:3)   : Pade coefficients in numerator for g, ice.
!   FCM_P_GL(2,0:3)   : Pade coefficients in numerator for g, liquid.
!   FCM_Q_AI(2,1:3)   : Pade coefficients in denominator for k_abs, ice.
!   FCM_Q_AL(2,1:3)   : Pade coefficients in denominator for k_abs, liquid.
!   FCM_Q_DI(2,1:3)   : Pade coefficients in denominator for k_scat, ice.
!   FCM_Q_DL(2,1:3)   : Pade coefficients in denominator for k_scat, liquid.
!   FCM_Q_GI(2,1:3)   : Pade coefficients in denominator for g, ice.
!   FCM_Q_GL(2,1:3)   : Pade coefficients in denominator for g, liquid.

! Fitting parameters for statistical model (NER):
! 1) new (NRAY=2):
!   FSM_AA   - coefficients for EBL_min
!   FSM_BB   - coefficients for EBL_max
! 2) old (NRAY=1, LNEWSTAT=.T.):
!   FSM_CC   - parameter for tuning profile ZMAN
!   FSM_DD   - parameter for tuning profile ZMAN
!   FSM_EE   - parameter for tuning profile ZMAK
!   FSM_FF   - parameter for tuning profile ZMAK
!   FSM_GG   - parameter for tuning profile ZMAK
!   FSM_HH   - parameter for tuning profile ZMAK
!   FSM_II   - parameter for tuning profile ZMAK

! Fitting parameters for gaseous transmissions (FGT):
!   - symbol [X] denotes spectral band:
!       S - solar
!       C - thermal with B(T0) weights
!       T - thermal with dB/dT(T0) weights
!   - reference values are at T=1K and p=1Pa

!   FGT[X]_A0      - reference value for a (weak   line Malkmus coefficient)
!   FGT[X]_B0      - reference value for b (strong line Malkmus coefficient)
!   FGT[X]_C0      - reference value for c (continuum)
!   FGT[X]_A_ALPHA - temperature exponent for a
!   FGT[X]_B_ALPHA - temperature exponent for b
!   FGT[X]_B_GAMMA - pressure    exponent for b
!   FGT[X]_C_ALPHA - temperature exponent for c
!   FGT[X]_ALPHA   - power alpha           for scaling optical depth
!   FGT[X]_DELTA0  - critical value delta0 for scaling optical depth
!   FGT[X]_D1      - parameter d1 for H4O2
!   FGT[X]_D2      - parameter d2 for H4O2
!   FGT[X]_OA      - H2O-H4O2 overlap parameter a
!   FGT[X]_OB      - H2O-H4O2 overlap parameter b
!   FGT[X]_ODEL    - H2O-H4O2 overlap parameter delta0

! Fitting parameters for Rayleigh scattering (FRS, solar band only):
!   FRS_K_SCAT0 - unsaturated value of scattering coefficient [1/Pa]
!   FRS_P_CRIT  - critical pressure                           [Pa]
!   FRS_BETA    - saturation exponent                         [1]

! Albedo related parameters:
!   RLAMB_SOLID - proportion of Lambertian scattering for solid surfaces
!   RLAMB_WATER - proportion of Lambertian scattering for water

! Bracketing parameters:
!   RMIXD  - critical distance between EBL_min and EBL_max [W/m^2]
!   RMIXP0 - clearsky weight when EBL_min = EBL_max        [1]

! number of spectral bands used in cloud model
INTEGER(KIND=JPIM) :: N_SPBAND=N_SPBANDL

! linearization temperature T0 used in Planck weights [K]
REAL(KIND=JPRB) :: RTL=255.8_JPRB

REAL(KIND=JPRB) :: BSFSA(6)
REAL(KIND=JPRB) :: BSFSI
REAL(KIND=JPRB) :: BSFSN
REAL(KIND=JPRB) :: BSFTA(6)
REAL(KIND=JPRB) :: BSFTI
REAL(KIND=JPRB) :: BSFTN
REAL(KIND=JPRB) :: EARRT
REAL(KIND=JPRB) :: EOASA(6)
REAL(KIND=JPRB) :: EOASI
REAL(KIND=JPRB) :: EOASN
REAL(KIND=JPRB) :: EOATA(6)
REAL(KIND=JPRB) :: EOATI
REAL(KIND=JPRB) :: EOATN
REAL(KIND=JPRB) :: EODSA(6)
REAL(KIND=JPRB) :: EODSI
REAL(KIND=JPRB) :: EODSN
REAL(KIND=JPRB) :: EODTA(6)
REAL(KIND=JPRB) :: EODTI
REAL(KIND=JPRB) :: EODTN
REAL(KIND=JPRB) :: EORAY
REAL(KIND=JPRB) :: GCA(6)
REAL(KIND=JPRB) :: GCB(6)
REAL(KIND=JPRB) :: GCC(6)
REAL(KIND=JPRB) :: GCD4
REAL(KIND=JPRB) :: GCE4
REAL(KIND=JPRB) :: GIREC1
REAL(KIND=JPRB) :: GIREC2
REAL(KIND=JPRB) :: GIREC3
REAL(KIND=JPRB) :: GIREC4
REAL(KIND=JPRB) :: QCO2
REAL(KIND=JPRB) :: RII0
REAL(KIND=JPRB) :: USAA(6)
REAL(KIND=JPRB) :: USAI
REAL(KIND=JPRB) :: USAN
REAL(KIND=JPRB) :: USBA
REAL(KIND=JPRB) :: USBI
REAL(KIND=JPRB) :: USBN
REAL(KIND=JPRB) :: VDP(5,6)
REAL(KIND=JPRB) :: VNP(5,6)
REAL(KIND=JPRB) :: REXP_NEB
REAL(KIND=JPRB) :: FCM_DEL_A(N_SPBANDL)
REAL(KIND=JPRB) :: FCM_DEL_D(N_SPBANDL)
REAL(KIND=JPRB) :: FCM_MU_A(N_SPBANDL)
REAL(KIND=JPRB) :: FCM_MU_D(N_SPBANDL)
REAL(KIND=JPRB) :: FCM_N_I
REAL(KIND=JPRB) :: FCM_N_L
REAL(KIND=JPRB) :: FCM_AI
REAL(KIND=JPRB) :: FCM_AL
REAL(KIND=JPRB) :: FCM_B_AI
REAL(KIND=JPRB) :: FCM_B_AL
REAL(KIND=JPRB) :: FCM_B_BI
REAL(KIND=JPRB) :: FCM_B_BL
REAL(KIND=JPRB) :: FCM_DEL_AI
REAL(KIND=JPRB) :: FCM_DEL_AL
REAL(KIND=JPRB) :: FCM_DEL_DI
REAL(KIND=JPRB) :: FCM_DEL_DL
REAL(KIND=JPRB) :: FCM_MU_AI
REAL(KIND=JPRB) :: FCM_MU_AL
REAL(KIND=JPRB) :: FCM_MU_DI
REAL(KIND=JPRB) :: FCM_MU_DL
REAL(KIND=JPRB) :: FCM_NU_AI
REAL(KIND=JPRB) :: FCM_NU_AL
REAL(KIND=JPRB) :: FCM_NU_DI
REAL(KIND=JPRB) :: FCM_NU_DL
REAL(KIND=JPRB) :: FCM_IWC2DE(N_SPBANDL,-2:3)
REAL(KIND=JPRB) :: FCM_LWC2RE(N_SPBANDL,-2:3)
REAL(KIND=JPRB) :: FCM_P_AI(N_SPBANDL,0:3)
REAL(KIND=JPRB) :: FCM_P_AL(N_SPBANDL,0:3)
REAL(KIND=JPRB) :: FCM_P_DI(N_SPBANDL,0:3)
REAL(KIND=JPRB) :: FCM_P_DL(N_SPBANDL,0:3)
REAL(KIND=JPRB) :: FCM_P_GI(N_SPBANDL,0:3)
REAL(KIND=JPRB) :: FCM_P_GL(N_SPBANDL,0:3)
REAL(KIND=JPRB) :: FCM_Q_AI(N_SPBANDL,1:3)
REAL(KIND=JPRB) :: FCM_Q_AL(N_SPBANDL,1:3)
REAL(KIND=JPRB) :: FCM_Q_DI(N_SPBANDL,1:3)
REAL(KIND=JPRB) :: FCM_Q_DL(N_SPBANDL,1:3)
REAL(KIND=JPRB) :: FCM_Q_GI(N_SPBANDL,1:3)
REAL(KIND=JPRB) :: FCM_Q_GL(N_SPBANDL,1:3)
REAL(KIND=JPRB) :: FSM_AA(0:2)
REAL(KIND=JPRB) :: FSM_BB(0:2)
REAL(KIND=JPRB) :: FSM_CC
REAL(KIND=JPRB) :: FSM_DD
REAL(KIND=JPRB) :: FSM_EE
REAL(KIND=JPRB) :: FSM_FF
REAL(KIND=JPRB) :: FSM_GG
REAL(KIND=JPRB) :: FSM_HH
REAL(KIND=JPRB) :: FSM_II
REAL(KIND=JPRB) :: GOLCA(6)
REAL(KIND=JPRB) :: GOLCB(6)
REAL(KIND=JPRB) :: GOLCC(6)
REAL(KIND=JPRB) :: FGTS_A     (3,0:2)
REAL(KIND=JPRB) :: FGTS_B     (3,0:2)
REAL(KIND=JPRB) :: FGTS_C     (5)
REAL(KIND=JPRB) :: FGTS_ALPHA (3)
REAL(KIND=JPRB) :: FGTS_DELTA0(3)
REAL(KIND=JPRB) :: FGTS_D     (3)
REAL(KIND=JPRB) :: FGTS_P00   (3,0:2)
REAL(KIND=JPRB) :: FGTS_P     (3,0:5,0:2)
REAL(KIND=JPRB) :: FGTS_Q     (3,0:2)
REAL(KIND=JPRB) :: FGTS_OA    (3)
REAL(KIND=JPRB) :: FGTS_OB    (3)
REAL(KIND=JPRB) :: FGTS_OC    (3)
REAL(KIND=JPRB) :: FGTS_OD    (3)
REAL(KIND=JPRB) :: FGTC_A     (3,0:2)
REAL(KIND=JPRB) :: FGTC_B     (3,0:2)
REAL(KIND=JPRB) :: FGTC_C     (5)
REAL(KIND=JPRB) :: FGTC_ALPHA (3)
REAL(KIND=JPRB) :: FGTC_DELTA0(3)
REAL(KIND=JPRB) :: FGTC_D     (3)
REAL(KIND=JPRB) :: FGTC_P00   (3,0:2)
REAL(KIND=JPRB) :: FGTC_P     (3,0:5,0:2)
REAL(KIND=JPRB) :: FGTC_Q     (3,0:2)
REAL(KIND=JPRB) :: FGTC_OA    (3)
REAL(KIND=JPRB) :: FGTC_OB    (3)
REAL(KIND=JPRB) :: FGTC_OC    (3)
REAL(KIND=JPRB) :: FGTC_OD    (3)
REAL(KIND=JPRB) :: FGTT_A     (3,0:2)
REAL(KIND=JPRB) :: FGTT_B     (3,0:2)
REAL(KIND=JPRB) :: FGTT_C     (5)
REAL(KIND=JPRB) :: FGTT_ALPHA (3)
REAL(KIND=JPRB) :: FGTT_DELTA0(3)
REAL(KIND=JPRB) :: FGTT_D     (3)
REAL(KIND=JPRB) :: FGTT_P00   (3,0:2)
REAL(KIND=JPRB) :: FGTT_P     (3,0:5,0:2)
REAL(KIND=JPRB) :: FGTT_Q     (3,0:2)
REAL(KIND=JPRB) :: FGTT_OA    (3)
REAL(KIND=JPRB) :: FGTT_OB    (3)
REAL(KIND=JPRB) :: FGTT_OC    (3)
REAL(KIND=JPRB) :: FGTT_OD    (3)
REAL(KIND=JPRB) :: FRS_K_SCAT0
REAL(KIND=JPRB) :: FRS_P_CRIT
REAL(KIND=JPRB) :: FRS_BETA
REAL(KIND=JPRB) :: RLAMB_SOLID
REAL(KIND=JPRB) :: RLAMB_WATER
REAL(KIND=JPRB) :: RMIXD
REAL(KIND=JPRB) :: RMIXP0

END TYPE TPHY3

!!TYPE(TPHY3), POINTER :: YRPHY3 => NULL()

!     ------------------------------------------------------------------
END MODULE YOMPHY3
