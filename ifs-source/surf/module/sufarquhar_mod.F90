MODULE SUFARQUHAR_MOD
CONTAINS
SUBROUTINE SUFARQUHAR(PRVCMAX25,PRHUMREL,PRA1,PRB1,PRG0,PRGM25,PRE_VCMAX,PRE_JMAX,YDVEG,YDAGS,YDAGF)
!***
! (C) Copyright 2021- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**   *SUFARQUHAR* - DOES THE INITIALISATION OF AGS PARAMETERS
 
!    PURPOSE
!    -------

!     Initialize model to calculate net assimilation of 
!     CO2 and stomatal conductance.
!              
!     METHOD
!    ------
!     Yin et al. (2009) [from model of Farquhar(1980)]

!     EXTERNAL
!     --------
!     none

!     REFERENCE
!     ---------

!     Yin et al. (2009)

!     MODIFICATIONS
!     -------------
!     V.Bastrikov,F.Maignan,P.Peylin,A.Agusti-Panareda/S.Boussetta Feb 2021 Add Farquhar photosynthesis model


!     -------------------------------------------------------------------------
!
USE PARKIND1, ONLY : JPIM, JPRB, JPRM
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOS_AGS         ,ONLY : TAGS
USE YOS_AGF         ,ONLY : TAGF
USE YOS_VEG         ,ONLY : TVEG

IMPLICIT NONE

REAL(KIND=JPRB)   ,INTENT(IN)    :: PRVCMAX25(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHUMREL(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRA1(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRB1(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRG0(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRGM25(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE_VCMAX(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE_JMAX(:,:)
TYPE(TAGS)        ,INTENT(IN)    :: YDAGS
TYPE(TVEG)        ,INTENT(IN)    :: YDVEG
TYPE(TAGF)        ,INTENT(INOUT) :: YDAGF

INTEGER(KIND=JPIM) :: ICO2TYP,JSFC

!REAL(KIND=JPRB), ALLOCATABLE, SAVE :: temp_growth(:)         ! Growth temperature (°C) - Is equal to t2m_month
!LOGICAL, SAVE                      :: FIRST_CALL=.TRUE.
!INTEGER(KIND=JPIM), SAVE           :: KSTEP_SAV=-1_JPIM

INTEGER(KIND=JPIM), PARAMETER :: MAX_TEMP=370_JPIM           ! Maximum temperature for saturated humidity (K) and also used as 
							     ! the size of local array to keep saturated humidity (unitless)
INTEGER(KIND=JPIM), PARAMETER :: NLAI = 20_JPIM              ! Number of LAI levels (unitless)
INTEGER(KIND=JPIM), PARAMETER :: NLAI1 = 21_JPIM              ! Number of LAI levels + 1 (unitless)
INTEGER(KIND=JPIM), PARAMETER :: NPATH = 2_JPIM              ! Number of photosunthetic pathways, C3 and C4 (unitless)


!notused REAL(KIND=JPRB), DIMENSION(max_temp),SAVE :: qsfrict         ! Array to keep water vapor pressure at saturation for each temperature level 
!notused                                                              ! (hPa)
!notused REAL(KIND=JPRB), PARAMETER :: zero = 0._JPRB                 ! Numerical constant set to 0 (unitless)
!notused REAL(KIND=JPRB), PARAMETER :: one = 1._JPRB                  ! Numerical constant set to 1 (unitless)
!notused REAL(KIND=JPRB), PARAMETER :: two = 2._JPRB                  ! Numerical constant set to 2 (unitless)
!notused REAL(KIND=JPRB), PARAMETER :: three = 3._JPRB                ! Numerical constant set to 3 (unitless)
!notused REAL(KIND=JPRB), PARAMETER :: four = 4._JPRB                 ! Numerical constant set to 4 (unitless)
!notused REAL(KIND=JPRB), PARAMETER :: min_sechiba = 1.E-8_JPRB       ! Epsilon to detect a near zero floating point (unitless)

!defined in suscst.F90 REAL(KIND=JPRB), PARAMETER :: RR = 8.314_JPRB                ! Ideal gas constant (J.mol^{-1}.K^{-1})
!notused REAL(KIND=JPRB), PARAMETER :: min_wind = 0.1_JPRB            ! The minimum wind (m.s^{-1})

!defined in yos_ags (RMAIR) (they do not seem to be used in farquhar_mod)
!REAL(KIND=JPRB), PARAMETER :: msmlr_air = 28.964E-03_JPRB    ! Molecular weight of dry air (kg.mol^{-1})
!REAL(KIND=JPRB), PARAMETER :: msmlr_h2o = 18.02E-03_JPRB     ! Molecular weight of water vapor (kg.mol^{-1}) 

REAL(KIND=JPRB),PARAMETER :: RDOWNREGULATION_CO2_BASELEVEL = 380._JPRB ! CO2 base level (ppm)
REAL(KIND=JPRB),PARAMETER :: RDOWNREGULATION_CO2_MINIMUM = 280._JPRB   ! CO2 value above which downregulation is taken into account (ppm)


! delete two lines below (PI is already defined in YOS_CTS)
! Mathematical and numerical constants
!REAL(KIND=JPRB), PARAMETER :: PI = 3.141592653589793238_JPRB! pi souce : http://mathworld.wolfram.com/Pi.html (unitless)

!notused REAL(KIND=JPRB), PARAMETER :: one_point_five = 1.5_JPRB     ! Numerical constant set to 1.5 (unitless)
!notused REAL(KIND=JPRB), PARAMETER :: nine = 9._JPRB                ! Numerical constant set to 9 (unitless)
!notused REAL(KIND=JPRB), PARAMETER :: twenty_seven = 27._JPRB       ! Numerical constant set to 27 (unitless)
!notused REAL(KIND=JPRB), PARAMETER :: fifty_four = 54._JPRB         ! Numerical constant set to 54 (unitless)

REAL(KIND=JPRB), PARAMETER :: RUNDEF = -9999._JPRB           ! Undefined/missing value
REAL(KIND=JPRB), PARAMETER :: NIUNDEF = 9999_JPIM          ! undefined/missing integer value


REAL(KIND=JPRB), PARAMETER :: RONE_MONTH=60._JPRB*60._JPRB*24._JPRB*30._JPRB ! One month (s)
REAL(KIND=JPRB), PARAMETER :: RUNDEF_SECHIBA = 1.E+20_JPRB   ! The undef value used in SECHIBA (unitless)

                                                            
REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RARJV= &     ! a coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio (mu mol e- (mu mol CO2)-1)
 & (/2.59_JPRB, 1.715_JPRB/)                                ! See Table 3 of Kattge & Knorr (2007)
                                                            ! For C4 plants, we assume that there is no
                                                            ! acclimation and that for a temperature of 25°C, ASV is the same for both C4 and C3 plants (no strong jusitification - need further parametrization)
REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RBRJV= &     ! b coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio ((mu mol e- (mu mol CO2)-1) (°C)-1)
 & (/-0.035_JPRB, 0._JPRB/)                                    ! See Table 3 of Kattge & Knorr (2007)
                                                            ! We assume No acclimation term for C4 plants.

REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RASJ= &      ! a coefficient of the linear regression (a+bT) defining the Entropy term for Jmax (J K-1 mol-1)
 & (/659.7_JPRB, 630._JPRB/)                                ! See Table 3 of Kattge & Knorr (2007)
                                                            ! and Table 2 of Yin et al. (2009) for C4 plants
REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RBSJ= &      ! b coefficient of the linear regression (a+bT) defining the Entropy term for Jmax (J K-1 mol-1 °C-1)
 & (/-0.75_JPRB, 0._JPRB/)                                  ! See Table 3 of Kattge & Knorr (2007)
                                                            ! We assume no acclimation term for C4 plants.

REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RASV= &      ! a coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax (J K-1 mol-1)
 & (/668.39_JPRB, 641.64_JPRB/)                             ! See Table 3 of Kattge & Knorr (2007)
                                                            ! For C4 plants, we assume that there is no
                                                            ! acclimation and that at for a temperature of 25°C, ASV is the same for both C4 and C3 plants (no strong jusitification - need further parametrization)
REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RBSV= &      ! b coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax (J K-1 mol-1 °C-1)
 & (/-1.07_JPRB, 0.0_JPRB/)                                     ! See Table 3 of Kattge & Knorr (2007)
                                                            ! We assume No acclimation term for C4 plants.

REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RD_VCMAX= &  ! Energy of deactivation for Vcmax (J mol-1)
 & (/200000._JPRB, 192000._JPRB/)                           ! Medlyn et al. (2002) also uses 200000. for C3 plants (same value than D_JMAX)
                                                            ! 'Consequently', we use the value of D_JMAX for C4 plants
REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RD_JMAX= &   ! Energy of deactivation for Jmax (J mol-1)
 & (/200000._JPRB, 192000._JPRB/)                           ! See Table 2 of Yin et al. (2009)
                                                            ! Medlyn et al. (2002) also uses 200000. for C3 plants

REAL(KIND=JPRB), PARAMETER :: RE_GAMMA_STAR=37830._JPRB      ! Energy of activation for gamma_star (J mol-1)
                                                            ! See Medlyn et al. (2002) from Bernacchi al. (2001) 
                                                            ! for C3 plants - We use the same values for C4 plants.

REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RE_GM= &     ! Energy of activation for gm (J mol-1) 
 & (/49600._JPRB, RUNDEF/)                                   ! See Table 2 of Yin et al. (2009) 

REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RS_GM= &     ! Entropy term for gm (J K-1 mol-1) 
 & (/1400._JPRB, RUNDEF/)                                    ! See Table 2 of Yin et al. (2009) 

REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RS_JMAX= &     ! Entropy term for Jmax (J K-1 mol-1) 
 & (/650._JPRB, 630._JPRB/)                                    ! See Table 2 of Yin et al. (2009) 

REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RD_GM= &     ! Energy of deactivation for gm (J mol-1) 
 & (/437400._JPRB, RUNDEF/)                                  ! See Table 2 of Yin et al. (2009) 
                                                            
!orig REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RE_VCMAX= &  ! Energy of activation for Vcmax (J mol-1)
!orig  & (/71513._JPRB, 67300._JPRB/)                             ! See Table 2 of Yin et al. (2009) for C4 plants
                                                            ! and Kattge & Knorr (2007) for C3 plants (table 3)
!orig REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RE_JMAX= &   ! Energy of activation for Jmax (J mol-1)
!orig  & (/49884._JPRB, 77900._JPRB/)                             ! See Table 2 of Yin et al. (2009) for C4 plants
                                                            ! and Kattge & Knorr (2007) for C3 plants (table 3)
REAL(KIND=JPRB), PARAMETER :: RE_KMC=79430._JPRB             ! Energy of activation for KmC (J mol-1)
                                                            ! See Medlyn et al. (2002) 
                                                            ! from Bernacchi al. (2001)
REAL(KIND=JPRB), PARAMETER :: RE_KMO=36380._JPRB             ! Energy of activation for KmO (J mol-1)
                                                            ! See Medlyn et al. (2002) 
                                                            ! from Bernacchi al. (2001)

REAL(KIND=JPRB), PARAMETER :: RE_RD=46390._JPRB              ! Energy of activation for Rd (J mol-1)
                                                            ! See Table 2 of Yin et al. (2009)
  
REAL(KIND=JPRB), PARAMETER :: RE_SCO=-24460._JPRB            ! Energy of activation for Sco (J mol-1)
                                                            ! See Table 2 of Yin et al. (2009)
                                                            ! Value for C4 plants is not mentioned - We use C3 for all plants.
                                                            
!orig REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RG0= &       ! Residual stomatal conductance when irradiance approaches zero (mol CO2 m−2 s−1 bar−1)
!orig  & (/0.00625_JPRB, 0.01875_JPRB/)                           ! Value from ORCHIDEE - No other reference.
                                                            ! modify to account for the conversion for conductance to H2O to CO2 
REAL(KIND=JPRB), PARAMETER :: RGAMMA_STAR25=42.75_JPRB       ! Ci-based CO2 compensation point in the absence of Rd at 25°C (ubar)
                                                            ! See Medlyn et al. (2002) for C3 plants - For C4 plants, we use the same value (probably uncorrect)
!orig REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RGM25= &     ! Mesophyll diffusion conductance at 25°C (mol m-2 s-1 bar-1) 
!orig  & (/0.4_JPRB, RUNDEF/)                                      ! See legend of Figure 6 of Yin et al. (2009) 
                                                            ! and review by Flexas et al. (2008) - gm is not used for C4 plants 
REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RKMC25= &    ! Michaelis–Menten constant of Rubisco for CO2 at 25°C (ubar)
 & (/404.9_JPRB, 650._JPRB/)                                ! See Table 2 of Yin et al. (2009) for C4
                                                            ! and Medlyn et al (2002) for C3
REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RKMO25= &    ! Michaelis–Menten constant of Rubisco for O2 at 25°C (ubar)
 & (/278400._JPRB, 450000._JPRB/)                           ! See Table 2 of Yin et al. (2009) for C4 plants and Medlyn et al. (2002) for C3


REAL(KIND=JPRB), PARAMETER, DIMENSION(NPATH) :: RSCO25= &    ! Relative CO2 /O2 specificity factor for Rubisco at 25Â°C (bar bar-1)
 & (/2800._JPRB, 2590._JPRB/)                               ! See Table 2 of Yin et al. (2009)

REAL(KIND=JPRB), PARAMETER :: RSTRESS_GM=1._JPRB             ! Water stress on gm
REAL(KIND=JPRB), PARAMETER :: RSTRESS_GS=1._JPRB             ! Water stress on gs
REAL(KIND=JPRB), PARAMETER :: RSTRESS_VCMAX=1._JPRB          ! Water stress on vcmax

REAL(KIND=JPRB), PARAMETER :: RTPHOTO_MAX=55._JPRB           ! maximum photosynthesis temperature (deg C) 
REAL(KIND=JPRB), PARAMETER :: RTPHOTO_MIN=-4._JPRB           ! minimum photosynthesis temperature (deg C) 
 
REAL(KIND=JPRB), PARAMETER :: RLAI_LEVEL_DEPTH = 0.15_JPRB 
REAL(KIND=JPRB), PARAMETER :: RLAIMAX = 12._JPRB             ! Maximal LAI used for splitting LAI into N layers (m^2.m^{-2})
REAL(KIND=JPRB), PARAMETER :: REXT_COEFF = 0.5               ! extinction coefficient of the Monsi&Saeki relationship (1953) (unitless)    !PFT-dependant

REAL(KIND=JPRB), PARAMETER :: RGB_REF = 1./25._JPRB          ! Leaf bulk boundary layer resistance (s m-1)
REAL(KIND=JPRB), PARAMETER :: RTP_00= 273.15_JPRB            ! 0 degree Celsius in degree Kelvin (K)
REAL(KIND=JPRB), PARAMETER :: RPB_STD = 1013._JPRB           ! standard pressure (hPa)
REAL(KIND=JPRB), PARAMETER :: RRG_TO_PAR = 0.5_JPRB
REAL(KIND=JPRB), PARAMETER :: RW_TO_MOL = 4.6_JPRB           ! W_to_mmol * RG_TO_PAR = 2.3

REAL(KIND=JPRB), PARAMETER :: RALPHA_LL = 0.3_JPRB           ! Conversion efficiency of absorbed light into J at strictly limiting light (mol e− (mol photon)−1)
                                                            ! See comment from Yin et al. (2009) after eq. 4
                                                            ! ALPHA value from Medlyn et al. (2002)   
                                                            ! 0.093 mol CO2 fixed per mol absorbed photons
                                                            ! times 4 mol e- per mol CO2 produced
                                                            ! PFT-dependant
REAL(KIND=JPRB), PARAMETER :: RTHETA = 0.7_JPRB              ! Convexity factor for response of J to irradiance (-)        
                                                            ! See Table 2 of Yin et al. (2009) 
                                                            ! PFT-dependant

!for C4 species
REAL(KIND=JPRB), PARAMETER :: RFPSIR = 0.4_JPRB              ! Fraction of PSII e− transport rate partitioned to the C4 cycle (-)     
                                                            ! See Table 2 of Yin et al. (2009) - x parameter
                                                            ! PFT-dependant

REAL(KIND=JPRB), PARAMETER :: RFQ = 1._JPRB                  ! Fraction of electrons at reduced plastoquinone    
                                                            ! that follow the Q-cycle (-) - Values for C3 plants are not used.
                                                            ! See Table 2 of Yin et al. (2009) 
                                                            ! PFT-dependant

REAL(KIND=JPRB), PARAMETER :: RFPSEUDO = 0.1_JPRB            ! Fraction of electrons at PSI that follow 
                                                            ! pseudocyclic transport (-) - Values for C3 plants are not used.
                                                            ! See Table 2 of Yin et al. (2009) 
                                                            ! PFT-dependant
REAL(KIND=JPRB), PARAMETER :: RH_PROTONS = 4._JPRB           ! Number of protons required to produce one ATP (mol mol-1)
                                                            ! See Table 2 of Yin et al. (2009) - h parameter
                                                            ! PFT-dependant

REAL(KIND=JPRB), PARAMETER :: RGBS = 0.003_JPRB              ! Bundle-sheath conductance (mol m−2 s−1 bar−1)
                                                            ! See legend of Figure 6 of Yin et al. (2009)
                                                            ! PFT-dependant
REAL(KIND=JPRB), PARAMETER :: ROI = 210000._JPRB             ! Intercellular oxygen partial pressure (ubar)
REAL(KIND=JPRB), PARAMETER :: RALPHA = 0.1_JPRB              ! Fraction of PSII activity in the bundle sheath (-)
                                                            ! See legend of Figure 6 of Yin et al. (2009)
                                                            ! PFT-dependant

REAL(KIND=JPRB), PARAMETER :: RKP = 0.7_JPRB                 ! Initial carboxylation efficiency of the PEP carboxylase (mol m−2 s−1 bar−1) 
                                                            ! See Table 2 of Yin et al. (2009)
                                                            ! PFT-dependant

REAL(KIND=JPRB), PARAMETER :: RTETENS_1 = 0.622_JPRB         ! Ratio between molecular weight of water vapor and molecular weight  
                                                            ! of dry air (unitless)
REAL(KIND=JPRB), PARAMETER :: RTETENS_2 = 0.378_JPRB         ! 1-TETENS_1
REAL(KIND=JPRB), PARAMETER :: RREF_TEMP = 298._JPRB
REAL(KIND=JPRB), PARAMETER :: RMOL_TO_M_1 = 0.0244           ! ??? conversion from mol/m^2/s to m/s (for stomatal conductance, Pearcy, Schulze and Zimmermann)

REAL(KIND=JPRB), PARAMETER :: RRATIO_H2O_TO_CO2 = 1.6        ! Ratio of water vapor diffusivity to the CO2 diffusivity (unitless)

REAL(KIND=JPRB), PARAMETER :: RN_VERT_ATT = .7_JPRB          ! N vertical attenuation factor within the canopy

REAL(KIND=JPRB), PARAMETER :: RRVEG_PFT = 1._JPRB                ! Potentiometer to set vegetation resistance (unitless)

REAL(KIND=JPRB), PARAMETER :: RTAU_TEMP_AIR_MONTH = 60._JPRB*60._JPRB*24._JPRB*20. ! Relaxation constant for computing monthly temperature average (s)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUFARQUHAR_MOD:SUFARQUHAR',0,ZHOOK_HANDLE)
!*         1.2 INITIALIZE PARAMETERS DEPENDING ON VEGETATION TYPE
!              ---------- ---------- --------- -- ---------- ----


! BATS table matching with ECOCLIMAP table 
! (1)  ! Crops, Mixed Farming			=>! 6  C3 CROPS
! (2)  ! Short Grass				    =>! 4  C3 GRASS
! (3)  ! Evergreen Needleleaf Trees		=>! 2  CONIFEROUS
! (4)  ! Deciduous Needleleaf Trees		=>! 2  CONIFEROUS
! (5)  ! Deciduous Broadleaf Trees		=>! 1  DECIDUOUS
! (6)  ! Evergreen Broadleaf Trees		=>! 3  EVERGREEN
! (7)  ! Tall Grass				        =>! 5  C4 GRASS7
! (8)  ! Desert					        =>! 
! (9)  ! Tundra					        =>! 4  C3 GRASS
! (10) ! Irrigated Crops			    =>! 6  C3 CROPS
! (11) ! Semidesert				        =>! 4  C3 GRASS
! (12) ! Ice Caps and Glaciers			=>!
! (13) ! Bogs and Marshes			    =>! 4  C3 GRASS
! (14) ! Inland Water				    =>!
! (15) ! Ocean					        =>!
! (16) ! Evergreen Shrubs		   	    =>! 5  C4 GRASS
! (17) ! Deciduous Shrubs			    =>! 4  C3 GRASS
! (18) ! Mixed Forest/woodland			=>! 2  CONIFEROUS
! (19) ! Interrupted Forest			    =>! 2  CONIFEROUS
! (20) ! Water and Land Mixtures		=>! 4  C3 GRASS


ASSOCIATE(NVTYPES=>YDVEG%NVTYPES)





!YDAGF%Vcmax25(1:20)=Vcmax25(1:20)
!YDAGF%DOWNREGULATION_CO2_COEF(1:20)=DOWNREGULATION_CO2_COEF(1:20)
!YDAGF%RSTRUCT_CONST(1:20)=RSTRUCT_CONST(1:20)
!YDAGF%NCO2TYP(1:20)=NCO2TYP(1:20)

!IF (.NOT.ALLOCATED(YDAGF%NCO2TYP)) ALLOCATE (YDAGF%NCO2TYP(NVTYPES,2)) 
IF (.NOT.ALLOCATED(YDAGF%RVCMAX25)) ALLOCATE (YDAGF%RVCMAX25(NVTYPES,2)) 
IF (.NOT.ALLOCATED(YDAGF%RHUMREL)) ALLOCATE (YDAGF%RHUMREL(NVTYPES,2)) 
IF (.NOT.ALLOCATED(YDAGF%RA1)) ALLOCATE (YDAGF%RA1(NVTYPES,2)) 
IF (.NOT.ALLOCATED(YDAGF%RB1)) ALLOCATE (YDAGF%RB1(NVTYPES,2)) 
IF (.NOT.ALLOCATED(YDAGF%RG0)) ALLOCATE (YDAGF%RG0(NVTYPES,2)) 
IF (.NOT.ALLOCATED(YDAGF%RGM25)) ALLOCATE (YDAGF%RGM25(NVTYPES,2)) 
IF (.NOT.ALLOCATED(YDAGF%RE_VCMAX)) ALLOCATE (YDAGF%RE_VCMAX(NVTYPES,2)) 
IF (.NOT.ALLOCATED(YDAGF%RE_JMAX)) ALLOCATE (YDAGF%RE_JMAX(NVTYPES,2)) 

IF (.NOT.ALLOCATED(YDAGF%RDOWNREGULATION_CO2_COEF)) ALLOCATE (YDAGF%RDOWNREGULATION_CO2_COEF(NVTYPES,2)) 
IF (.NOT.ALLOCATED(YDAGF%RSTRUCT_CONST)) ALLOCATE (YDAGF%RSTRUCT_CONST(NVTYPES,2))

 
! Replace with optimized values pass through namelist defined in yoephy

DO JSFC=1,NVTYPES ! Loop over PFTs
 DO ICO2TYP=1,2
  YDAGF%RVCMAX25(JSFC,ICO2TYP) = PRVCMAX25(JSFC,ICO2TYP)
  YDAGF%RHUMREL(JSFC,ICO2TYP) = PRHUMREL(JSFC,ICO2TYP)
  YDAGF%RG0(JSFC,ICO2TYP)=PRG0(JSFC,ICO2TYP)
  YDAGF%RGM25(JSFC,ICO2TYP)=PRGM25(JSFC,ICO2TYP)
  YDAGF%RA1(JSFC,ICO2TYP)=PRA1(JSFC,ICO2TYP)
  YDAGF%RB1(JSFC,ICO2TYP)=PRB1(JSFC,ICO2TYP)
  YDAGF%RE_VCMAX(JSFC,ICO2TYP)=PRE_VCMAX(JSFC,ICO2TYP)
  YDAGF%RE_JMAX(JSFC,ICO2TYP)=PRE_JMAX(JSFC,ICO2TYP)
  YDAGF%RDOWNREGULATION_CO2_COEF(JSFC,ICO2TYP)=RUNDEF
  YDAGF%RSTRUCT_CONST(JSFC,ICO2TYP)=RUNDEF

 ENDDO
ENDDO

YDAGF%RDOWNREGULATION_CO2_COEF(1,1)=0.24_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(1,2)=0.03_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(2,1)=0.24_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(3,1)=0.20_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(4,1)=0.20_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(5,1)=0.26_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(6,1)=0.26_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(7,2)=0.03_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(9,1)=0.24_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(10,1)=0.24_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(11,1)=0.03_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(13,1)=0.24_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(16,1)=0.03_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(17,1)=0.24_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(18,1)=0.26_JPRB
YDAGF%RDOWNREGULATION_CO2_COEF(19,1)=0.26_JPRB

! Coefficient for CO2 downregulation if downregulation_co2 (used for CMIP6 6.1.11) (unitless)
!REAL(KIND=JPRB), PARAMETER, DIMENSION(20) :: DOWNREGULATION_CO2_COEF = &  
!& (/0.03_JPRB,   0.24_JPRB,   0.20_JPRB,    0.20_JPRB,   0.26_JPRB, &
!&   0.26_JPRB,   0.03_JPRB, -9999._JPRB,    0.24_JPRB,   0.24_JPRB, &
!&   0.03_JPRB, -9999._JPRB,    0.24_JPRB, -9999._JPRB, -9999._JPRB, &
!&   0.03_JPRB,   0.24_JPRB,    0.26_JPRB,   0.26_JPRB, -9999._JPRB/) 


YDAGF%RSTRUCT_CONST(1,1)=2.0_JPRB
YDAGF%RSTRUCT_CONST(1,2)=2.0_JPRB
YDAGF%RSTRUCT_CONST(2,1)=2.5_JPRB
YDAGF%RSTRUCT_CONST(3,1)=25.0_JPRB
YDAGF%RSTRUCT_CONST(4,1)=25.0_JPRB
YDAGF%RSTRUCT_CONST(5,1)=25.0_JPRB
YDAGF%RSTRUCT_CONST(6,1)=25.0_JPRB
YDAGF%RSTRUCT_CONST(7,2)=2.0_JPRB
YDAGF%RSTRUCT_CONST(8,1)=0.0_JPRB
YDAGF%RSTRUCT_CONST(9,1)=2.5_JPRB
YDAGF%RSTRUCT_CONST(10,1)=2.0_JPRB
YDAGF%RSTRUCT_CONST(11,1)=2.0_JPRB
YDAGF%RSTRUCT_CONST(13,1)= 2.5_JPRB
YDAGF%RSTRUCT_CONST(16,1)=2.0_JPRB
YDAGF%RSTRUCT_CONST(17,1)=2.5_JPRB
YDAGF%RSTRUCT_CONST(18,1)=25.0_JPRB
YDAGF%RSTRUCT_CONST(19,1)=25.0_JPRB

!REAL(KIND=JPRB), PARAMETER, DIMENSION(20) :: RSTRUCT_CONST = & ! Structural resistance (s.m^{-1})
! & (/2.0_JPRB,     2.5_JPRB, 25.0_JPRB,    25.0_JPRB,   25.0_JPRB, &
! &  25.0_JPRB,     2.0_JPRB,  0.0_JPRB,     2.5_JPRB,    2.0_JPRB, &       
! &   2.0_JPRB, -9999._JPRB,   2.5_JPRB, -9999._JPRB,  -9999._JPRB, &
! &   2.0_JPRB,     2.5_JPRB, 25.0_JPRB,    25.0_JPRB, -9999._JPRB/)


YDAGF%NLAI = NLAI
YDAGF%NLAI1 = NLAI1
YDAGF%NPATH = NPATH

!delete (defined in ydcst) YDAGF%RR = RR
YDAGF%RDOWNREGULATION_CO2_BASELEVEL = RDOWNREGULATION_CO2_BASELEVEL
YDAGF%RDOWNREGULATION_CO2_MINIMUM = RDOWNREGULATION_CO2_MINIMUM
YDAGF%RUNDEF_SECHIBA=RUNDEF_SECHIBA
YDAGF%RUNDEF=RUNDEF

YDAGF%RONE_MONTH=RONE_MONTH
YDAGF%RARJV(1:2)=RARJV(1:2)
YDAGF%RBRJV(1:2)=RBRJV(1:2)
YDAGF%RASJ(1:2)=RASJ(1:2)
YDAGF%RBSJ(1:2)=RBSJ(1:2)
YDAGF%RASV(1:2)=RASV(1:2)
YDAGF%RBSV(1:2)=RBSV(1:2)
YDAGF%RD_VCMAX(1:2)=RD_VCMAX(1:2)
YDAGF%RD_JMAX(1:2)=RD_JMAX(1:2)
YDAGF%RE_GAMMA_STAR=RE_GAMMA_STAR
YDAGF%RE_GM(1:2)=RE_GM(1:2)
YDAGF%RS_GM(1:2)=RS_GM(1:2)

YDAGF%RS_JMAX(1:2)=RS_JMAX(1:2)

YDAGF%RD_GM(1:2)=RD_GM(1:2)
YDAGF%RE_KMC=RE_KMC
YDAGF%RE_KMO=RE_KMO
YDAGF%RE_RD=RE_RD
YDAGF%RE_SCO=RE_SCO
YDAGF%RGAMMA_STAR25=RGAMMA_STAR25

YDAGF%RKMC25(1:2)=RKMC25(1:2)
YDAGF%RKMO25(1:2)=RKMO25(1:2)
YDAGF%RSCO25(1:2)=RSCO25(1:2)

YDAGF%RSTRESS_GM=RSTRESS_GM
YDAGF%RSTRESS_GS=RSTRESS_GS
YDAGF%RSTRESS_VCMAX=RSTRESS_VCMAX

YDAGF%RTPHOTO_MAX=RTPHOTO_MAX
YDAGF%RTPHOTO_MIN=RTPHOTO_MIN

YDAGF%RLAI_LEVEL_DEPTH=RLAI_LEVEL_DEPTH
YDAGF%RLAIMAX=RLAIMAX
YDAGF%REXT_COEFF=REXT_COEFF
YDAGF%RGB_REF=RGB_REF
YDAGF%RTP_00=RTP_00
YDAGF%RPB_STD=RPB_STD
YDAGF%RRG_TO_PAR=RRG_TO_PAR
YDAGF%RW_TO_MOL=RW_TO_MOL
YDAGF%RALPHA_LL=RALPHA_LL
YDAGF%RTHETA=RTHETA

YDAGF%RFPSIR=RFPSIR
YDAGF%RFQ=RFQ
YDAGF%RFPSEUDO=RFPSEUDO
YDAGF%RH_PROTONS=RH_PROTONS
YDAGF%RGBS=RGBS
YDAGF%ROI=ROI
YDAGF%RALPHA=RALPHA
YDAGF%RKP=RKP

YDAGF%RTETENS_1=RTETENS_1
YDAGF%RTETENS_2=RTETENS_2
YDAGF%RREF_TEMP=RREF_TEMP
YDAGF%RMOL_TO_M_1=RMOL_TO_M_1
YDAGF%RRATIO_H2O_TO_CO2=RRATIO_H2O_TO_CO2
YDAGF%RN_VERT_ATT=RN_VERT_ATT
YDAGF%RRVEG_PFT=RRVEG_PFT
YDAGF%RTAU_TEMP_AIR_MONTH=RTAU_TEMP_AIR_MONTH

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SUFARQUHAR_MOD:SUFARQUHAR_MOD',1,ZHOOK_HANDLE)

END SUBROUTINE SUFARQUHAR
 
END MODULE SUFARQUHAR_MOD

