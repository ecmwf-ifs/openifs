MODULE YOS_AGF
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE
SAVE
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!     -----------------------------------------------------------------
!*    ** *AGF* - 
!     -----------------------------------------------------------------

TYPE :: TAGF
!REAL(KIND=JPRB), ALLOCATABLE :: temp_growth(:)     ! Growth temperature (°C) - Is equal to t2m_month
!LOGICAL            :: FIRST_CALL
!INTEGER(KIND=JPIM) :: KSTEP_SAV


INTEGER(KIND=JPIM) :: NLAI,NLAI1                         ! Number of LAI levels (unitless)
INTEGER(KIND=JPIM) :: NPATH                        ! Number of photosunthetic pathways, C3 and C4 (unitless)
INTEGER(KIND=JPIM) :: NIUNDEF                        ! Undefined/missing integer value


REAL(KIND=JPRB)    :: RDOWNREGULATION_CO2_BASELEVEL ! CO2 base level (ppm)
REAL(KIND=JPRB)    :: RDOWNREGULATION_CO2_MINIMUM   ! CO2 value above which downregulation is taken into account (ppm)

REAL(KIND=JPRB)    :: RUNDEF_SECHIBA                ! The undef value used in SECHIBA (unitless)
REAL(KIND=JPRB)    :: RUNDEF                        ! Undefined/missing real value

REAL(KIND=JPRB)    :: RONE_MONTH                    ! One month (s)

! Physics
                                                   ! See Table 2 of Yin et al. (2009)
                                                            
REAL(KIND=JPRB), DIMENSION(2) :: RARJV              ! a coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio (mu mol e- (mu mol CO2)-1)
                                                   ! See Table 3 of Kattge & Knorr (2007)
                                                   ! For C4 plants, we assume that there is no
                                                   ! acclimation and that for a temperature of 25°C, aSV is the same for both C4 and C3 plants (no strong jusitification - need further parametrization)
REAL(KIND=JPRB), DIMENSION(2) :: RBRJV              ! b coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio ((mu mol e- (mu mol CO2)-1) (°C)-1)
                                                   ! See Table 3 of Kattge & Knorr (2007)
                                                   ! We assume No acclimation term for C4 plants.

REAL(KIND=JPRB), DIMENSION(2) :: RASJ               ! a coefficient of the linear regression (a+bT) defining the Entropy term for Jmax (J K-1 mol-1)
                                                   ! See Table 3 of Kattge & Knorr (2007)
                                                   ! and Table 2 of Yin et al. (2009) for C4 plants
REAL(KIND=JPRB), DIMENSION(2) :: RBSJ               ! b coefficient of the linear regression (a+bT) defining the Entropy term for Jmax (J K-1 mol-1 °C-1)
                                                   ! See Table 3 of Kattge & Knorr (2007)
                                                   ! We assume no acclimation term for C4 plants.

REAL(KIND=JPRB), DIMENSION(2) :: RASV               ! a coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax (J K-1 mol-1)
                                                   ! See Table 3 of Kattge & Knorr (2007)
                                                   ! For C4 plants, we assume that there is no
                                                   ! acclimation and that at for a temperature of 25°C, aSV is the same for both C4 and C3 plants (no strong jusitification - need further parametrization)
REAL(KIND=JPRB), DIMENSION(2) :: RBSV               ! b coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax (J K-1 mol-1 °C-1)
                                                   ! See Table 3 of Kattge & Knorr (2007)
                                                   ! We assume No acclimation term for C4 plants.

REAL(KIND=JPRB), DIMENSION(2) :: RD_VCMAX          ! Energy of deactivation for Vcmax (J mol-1)
                                                   ! Medlyn et al. (2002) also uses 200000. for C3 plants (same value than D_Jmax)
                                                   ! 'Consequently', we use the value of D_Jmax for C4 plants
REAL(KIND=JPRB), DIMENSION(2) :: RD_JMAX            ! Energy of deactivation for Jmax (J mol-1)
                                                   ! See Table 2 of Yin et al. (2009)
                                                   ! Medlyn et al. (2002) also uses 200000. for C3 plants

REAL(KIND=JPRB) :: RE_GAMMA_STAR                   ! Energy of activation for gamma_star (J mol-1)
                                                   ! See Medlyn et al. (2002) from Bernacchi al. (2001) 
                                                   ! for C3 plants - We use the same values for C4 plants.

REAL(KIND=JPRB), DIMENSION(2) :: RE_GM             ! Energy of activation for gm (J mol-1) 
                                                   ! See Table 2 of Yin et al. (2009) 
REAL(KIND=JPRB), DIMENSION(2) :: RS_GM             ! Entropy term for gm (J K-1 mol-1) 
                                                   ! See Table 2 of Yin et al. (2009) 

REAL(KIND=JPRB), DIMENSION(2) :: RS_JMAX           ! Entropy term for Jmax (J K-1 mol-1) 
                                                   ! See Table 2 of Yin et al. (2009) 

REAL(KIND=JPRB), DIMENSION(2) :: RD_GM             ! Energy of deactivation for gm (J mol-1) 
                                                   ! See Table 2 of Yin et al. (2009) 
                                                            
REAL(KIND=JPRB) :: RE_KMC                          ! Energy of activation for KmC (J mol-1)
                                                   ! See Medlyn et al. (2002) 
                                                   ! from Bernacchi al. (2001)
REAL(KIND=JPRB) :: RE_KMO                          ! Energy of activation for KmO (J mol-1)
                                                   ! See Medlyn et al. (2002) 
                                                   ! from Bernacchi al. (2001)

REAL(KIND=JPRB) :: RE_RD                           ! Energy of activation for Rd (J mol-1)
                                                   ! See Table 2 of Yin et al. (2009)
  
REAL(KIND=JPRB) :: RE_SCO                          ! Energy of activation for Sco (J mol-1)
                                                   ! See Table 2 of Yin et al. (2009)
                                                   ! Value for C4 plants is not mentioned - We use C3 for all plants.
                                                            
                                                   ! Value from ORCHIDEE - No other reference.
                                                   ! modify to account for the conversion for conductance to H2O to CO2 
REAL(KIND=JPRB) :: RGAMMA_STAR25                   ! Ci-based CO2 compensation point in the absence of Rd at 25°C (ubar)
                                                   ! See Medlyn et al. (2002) for C3 plants - For C4 plants, we use the same value (probably uncorrect)
                                                   ! See legend of Figure 6 of Yin et al. (2009) 
                                                   ! and review by Flexas et al. (2008) - gm is not used for C4 plants 
REAL(KIND=JPRB), DIMENSION(2) :: RKMC25            ! Michaelis–Menten constant of Rubisco for CO2 at 25°C (ubar)
                                                   ! See Table 2 of Yin et al. (2009) for C4
                                                   ! and Medlyn et al (2002) for C3
REAL(KIND=JPRB), DIMENSION(2) :: RKMO25            ! Michaelis–Menten constant of Rubisco for O2 at 25°C (ubar)
                                                   ! See Table 2 of Yin et al. (2009) for C4 plants and Medlyn et al. (2002) for C3

REAL(KIND=JPRB), DIMENSION(2) :: RSCO25            ! Relative CO2 /O2 specificity factor for Rubisco at 25Â°C (bar bar-1)
                                                   ! See Table 2 of Yin et al. (2009)

REAL(KIND=JPRB) :: RSTRESS_GM                       ! Water stress on gm
REAL(KIND=JPRB) :: RSTRESS_GS                       ! Water stress on gs
REAL(KIND=JPRB) :: RSTRESS_VCMAX                    ! Water stress on vcmax

REAL(KIND=JPRB) :: RTPHOTO_MAX                      ! maximum photosynthesis temperature (deg C) 
REAL(KIND=JPRB) :: RTPHOTO_MIN                      ! minimum photosynthesis temperature (deg C) 

REAL(KIND=JPRB), ALLOCATABLE :: RVCMAX25(:,:)                  ! Maximum rate of Rubisco activity-limited carboxylation at 25°C (\mu mol.m^{-2}.s^{-1})
REAL(KIND=JPRB), ALLOCATABLE :: RHUMREL(:,:)                  ! scaling factor for soil moisture stress (optimized per PFT with FLUXNET data)
REAL(KIND=JPRB), ALLOCATABLE :: RA1(:,:)                ! Empirical factor involved in the calculation of fvpd (-)
                                                   ! See Table 2 of Yin et al. (2009)
REAL(KIND=JPRB), ALLOCATABLE :: RB1(:,:)                ! Empirical factor involved in the calculation of fvpd (-)

REAL(KIND=JPRB), ALLOCATABLE :: RG0(:,:)               ! Residual stomatal conductance when irradiance approaches zero (mol CO2 m−2 s−1 bar−1)
REAL(KIND=JPRB), ALLOCATABLE :: RGM25(:,:)             ! Mesophyll diffusion conductance at 25Â°C (mol m-2 s-1 bar-1) 
REAL(KIND=JPRB), ALLOCATABLE :: RE_VCMAX(:,:)          ! Energy of activation for Vcmax (J mol-1)
                                                   ! See Table 2 of Yin et al. (2009) for C4 plants
                                                   ! and Kattge & Knorr (2007) for C3 plants (table 3)
REAL(KIND=JPRB), ALLOCATABLE :: RE_JMAX(:,:)           ! Energy of activation for Jmax (J mol-1)
                                                   ! See Table 2 of Yin et al. (2009) for C4 plants
                                                   ! and Kattge & Knorr (2007) for C3 plants (table 3)


REAL(KIND=JPRB), ALLOCATABLE :: RDOWNREGULATION_CO2_COEF(:,:)  ! Coefficient for CO2 downregulation if downregulation_co2 (used for CMIP6 6.1.11) (unitless)
REAL(KIND=JPRB), ALLOCATABLE :: RSTRUCT_CONST(:,:)               ! Structural resistance (s.m^{-1})
 
REAL(KIND=JPRB) :: RLAI_LEVEL_DEPTH 
REAL(KIND=JPRB) :: RLAIMAX                         ! Maximal LAI used for splitting LAI into N layers (m^2.m^{-2})
REAL(KIND=JPRB) :: REXT_COEFF                       ! extinction coefficient of the Monsi&Saeki relationship (1953) (unitless)    !PFT-dependant

REAL(KIND=JPRB) :: RGB_REF                          ! Leaf bulk boundary layer resistance (s m-1)
REAL(KIND=JPRB) :: RTP_00                           ! 0 degree Celsius in degree Kelvin (K)
REAL(KIND=JPRB) :: RPB_STD                          ! standard pressure (hPa)
REAL(KIND=JPRB) :: RRG_TO_PAR  
REAL(KIND=JPRB) :: RW_TO_MOL                        ! W_to_mmol * RG_to_PAR = 2.3

REAL(KIND=JPRB) :: RALPHA_LL                        ! Conversion efficiency of absorbed light into J at strictly limiting light (mol e− (mol photon)−1)
                                                   ! See comment from Yin et al. (2009) after eq. 4
                                                   ! alpha value from Medlyn et al. (2002)   
                                                   ! 0.093 mol CO2 fixed per mol absorbed photons
                                                   ! times 4 mol e- per mol CO2 produced
                                                   ! PFT-dependant
REAL(KIND=JPRB) :: RTHETA                           ! Convexity factor for response of J to irradiance (-)        
                                                   ! See Table 2 of Yin et al. (2009) 
                                                   ! PFT-dependant

!for C4 species
REAL(KIND=JPRB) :: RFPSIR                           ! Fraction of PSII e− transport rate partitioned to the C4 cycle (-)     
                                                   ! See Table 2 of Yin et al. (2009) - x parameter
                                                   ! PFT-dependant

REAL(KIND=JPRB) :: RFQ                              ! Fraction of electrons at reduced plastoquinone    
                                                   ! that follow the Q-cycle (-) - Values for C3 plants are not used.
                                                   ! See Table 2 of Yin et al. (2009) 
                                                   ! PFT-dependant

REAL(KIND=JPRB) :: RFPSEUDO                        ! Fraction of electrons at PSI that follow 
                                                   ! pseudocyclic transport (-) - Values for C3 plants are not used.
                                                   ! See Table 2 of Yin et al. (2009) 
                                                   ! PFT-dependant
REAL(KIND=JPRB) :: RH_PROTONS                      ! Number of protons required to produce one ATP (mol mol-1)
                                                   ! See Table 2 of Yin et al. (2009) - h parameter
                                                   ! PFT-dependant

REAL(KIND=JPRB) :: RGBS                            ! Bundle-sheath conductance (mol m−2 s−1 bar−1)
                                                   ! See legend of Figure 6 of Yin et al. (2009)
                                                   ! PFT-dependant
REAL(KIND=JPRB) :: ROI                             ! Intercellular oxygen partial pressure (ubar)
REAL(KIND=JPRB) :: RALPHA                          ! Fraction of PSII activity in the bundle sheath (-)
                                                   ! See legend of Figure 6 of Yin et al. (2009)
                                                   ! PFT-dependant

REAL(KIND=JPRB) :: RKP                             ! Initial carboxylation efficiency of the PEP carboxylase (mol m−2 s−1 bar−1) 
                                                   ! See Table 2 of Yin et al. (2009)
                                                   ! PFT-dependant

REAL(KIND=JPRB) :: RTETENS_1                       ! Ratio between molecular weight of water vapor and molecular weight  
                                                   ! of dry air (unitless)
REAL(KIND=JPRB) :: RTETENS_2        
REAL(KIND=JPRB) :: RREF_TEMP  
REAL(KIND=JPRB) :: RMOL_TO_M_1        
REAL(KIND=JPRB) :: RRATIO_H2O_TO_CO2                ! Ratio of water vapor diffusivity to the CO2 diffusivity (unitless)
REAL(KIND=JPRB) :: RN_VERT_ATT                      ! N vertical attenuation factor within the canopy
REAL(KIND=JPRB) :: RRVEG_PFT                        ! Potentiometer to set vegetation resistance (unitless)
REAL(KIND=JPRB) :: RTAU_TEMP_AIR_MONTH              ! Relaxation constant for computing monthly temperature average (s)

END TYPE TAGF
END MODULE YOS_AGF
