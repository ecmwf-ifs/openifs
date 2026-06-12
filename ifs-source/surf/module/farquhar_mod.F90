MODULE FARQUHAR_MOD

  USE PARKIND1, ONLY : JPIM, JPRB, JPRM, JPRD
  USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

  USE YOS_AGS,  ONLY : TAGS
  USE YOS_AGF,  ONLY : TAGF
  USE YOS_CST,  ONLY : TCST

  USE SUFARQUHAR_MOD
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

  IMPLICIT NONE

  INTERFACE ARRHENIUS
     MODULE PROCEDURE  ARRHENIUS_FN, ARRHENIUS_MODIFIED_FN
  END INTERFACE


CONTAINS

!     --------------------------------------------------------------------------
! FUNCTION    : Arrhenius
!     --------------------------------------------------------------------------
FUNCTION ARRHENIUS_FN (PTEMP,PREF_TEMP,PENERGY_ACT,YDAGF,YDCST) RESULT(ARRHENIUS_VAL)

IMPLICIT NONE
REAL(KIND=JPRB), INTENT(IN) :: PTEMP             ! temperature (K)
REAL(KIND=JPRB), INTENT(IN) :: PREF_TEMP         ! temperature of reference (K)
REAL(KIND=JPRB), INTENT(IN) :: PENERGY_ACT       ! Activation Energy (J mol-1)
TYPE(TAGF)     , INTENT(IN) :: YDAGF
TYPE(TCST)     , INTENT(IN) :: YDCST

        
!*         0.     RESULT
!                 ------

REAL(KIND=JPRB) :: ARRHENIUS_VAL ! temperature dependance based on an Arrhenius function (-)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

ARRHENIUS_VAL=EXP(((PTEMP-PREF_TEMP)*PENERGY_ACT)/(PREF_TEMP*YDCST%R*(PTEMP)))

END FUNCTION ARRHENIUS_FN


!     --------------------------------------------------------------------------
! FUNCTION    : Arrhenius_modified_1d
!     --------------------------------------------------------------------------
FUNCTION ARRHENIUS_MODIFIED_FN (PTEMP,PREF_TEMP,PENERGY_ACT,PENERGY_DEACT,PENTROPY,YDAGF,YDCST) RESULT (VAL_ARRHENIUS)

REAL(KIND=JPRB), INTENT(IN) :: PTEMP          ! Temperature (K)
REAL(KIND=JPRB), INTENT(IN) :: PREF_TEMP      ! Temperature of reference (K)
REAL(KIND=JPRB), INTENT(IN) :: PENERGY_ACT    ! Activation Energy (J mol-1)
REAL(KIND=JPRB), INTENT(IN) :: PENERGY_DEACT  ! Deactivation Energy (J mol-1)
REAL(KIND=JPRB), INTENT(IN) :: PENTROPY       ! Entropy teZRM (J K-1 mol-1)
TYPE(TAGF)     , INTENT(IN) :: YDAGF
TYPE(TCST)     , INTENT(IN) :: YDCST
       
!*         0.     RESULT
!                 ------

REAL(KIND=JPRB) :: VAL_ARRHENIUS ! Temperature dependance based on an Arrhenius function (-)
    
VAL_ARRHENIUS = EXP(((PTEMP-PREF_TEMP)*PENERGY_ACT)/(PREF_TEMP*YDCST%R*(PTEMP))) &
       &  * (1._JPRB + EXP( (PREF_TEMP * PENTROPY - PENERGY_DEACT) / (PREF_TEMP * YDCST%R ))) &
       &  / (1._JPRB + EXP( (PTEMP * PENTROPY - PENERGY_DEACT) / ( YDCST%R*PTEMP)))
         
END FUNCTION ARRHENIUS_MODIFIED_FN



!! ========================================================================================

SUBROUTINE FARQUHAR (KIDIA, KFDIA, KLON, LDLAND, KVTYPE, KCO2TYP, YDAGS, YDAGF, YDCST, PSRFD, &
          &          PAPHM, PQS, PQSURF, PQM , PTM, PTSKM, PTSOIL, PF2, &
                                PCO2, PLAI, PGSTOT, PGPP, PRD, PAN)

!! This subroutine computes carbon assimilation and stomatal 
!! conductance, following respectively Farqhuar et al. (1980) and Ball et al. (1987).
!!
!! DESCRIPTION:
!! *** General:
!! The equations are different depending on the photosynthesis mode (C3 versus C4).
!! Assimilation and conductance are computed over 20 levels of LAI and then 
!! integrated at the canopy level.
!! This routine also computes partial beta coefficient: transpiration for each 
!! type of vegetation.
!! There is a main loop on the PFTs, then inner loops on the points where 
!! assimilation has to be calculated.
!! This subroutine is called at each sechiba time step by sechiba_main.
!! *** Details:
!! - Integration at the canopy level
!! - Light's extinction 
!! The available light follows a simple Beer extinction law. 
!! The extinction coefficients (ext_coef) are PFT-dependant constants.
!! - Estimation of relative humidity of air (for calculation of the stomatal conductance)
!! - Calculation of the water limitation factor
!! - Calculation of temperature dependent parameters for C4 plants
!! - Calculation of temperature dependent parameters for C3 plants
!! - Vmax scaling
!! Vmax is scaled into the canopy due to reduction of nitrogen 
!! (Johnson and Thornley,1984).
!! - Assimilation for C4 plants (Collatz et al., 1992)
!! - Assimilation for C3 plants (Farqhuar et al., 1980)
!! - Estimation of the stomatal conductance (Ball et al., 1987)
!!
!!
!! REFERENCE(S) :
!! - Ball, J., T. Woodrow, and J. Berry (1987), A model predicting stomatal 
!! conductance and its contribution to the control of photosynthesis under 
!! different environmental conditions, Prog. Photosynthesis, 4, 221– 224.
!! - Collatz, G., M. Ribas-Carbo, and J. Berry (1992), Coupled photosynthesis 
!! stomatal conductance model for leaves of C4 plants, Aust. J. Plant Physiol.,
!! 19, 519–538.
!! - Farquhar, G., S. von Caemmener, and J. Berry (1980), A biochemical model of 
!! photosynthesis CO2 fixation in leaves of C3 species, Planta, 149, 78–90.
!! - Johnson, I. R., and J. Thornley (1984), A model of instantaneous and daily
!! canopy photosynthesis, J Theor. Biol., 107, 531 545
!! - McMurtrie, R.E., Rook, D.A. and Kelliher, F.M., 1990. Modelling the yield of Pinus radiata on a
!! site limited by water and nitrogen. For. Ecol. Manage., 30: 381-413
!! - Bounoua, L., Hall, F. G., Sellers, P. J., Kumar, A., Collatz, G. J., Tucker, C. J., and Imhoff, M. L. (2010), Quantifying the 
!! negative feedback of vegetation to greenhouse warming: A modeling approach, Geophysical Research Letters, 37, Artn L23701, 
!! Doi 10.1029/2010gl045338
!! - Bounoua, L., Collatz, G. J., Sellers, P. J., Randall, D. A., Dazlich, D. A., Los, S. O., Berry, J. A., Fung, I., 
!! Tucker, C. J., Field, C. B., and Jensen, T. G. (1999), Interactions between vegetation and climate: Radiative and physiological 
!! effects of doubled atmospheric co2, Journal of Climate, 12, 309-324, Doi 10.1175/1520-0442(1999)012<0309:Ibvacr>2.0.Co;2
!! - Sellers, P. J., Bounoua, L., Collatz, G. J., Randall, D. A., Dazlich, D. A., Los, S. O., Berry, J. A., Fung, I., 
!! Tucker, C. J., Field, C. B., and Jensen, T. G. (1996), Comparison of radiative and physiological effects of doubled atmospheric
!! co2 on climate, Science, 271, 1402-1406, DOI 10.1126/science.271.5254.1402
!! - Lewis, J. D., Ward, J. K., and Tissue, D. T. (2010), Phosphorus supply drives nonlinear responses of cottonwood 
!! (populus deltoides) to increases in co2 concentration from glacial to future concentrations, New Phytologist, 187, 438-448, 
!! DOI 10.1111/j.1469-8137.2010.03307.x
!! - Kattge, J., Knorr, W., Raddatz, T., and Wirth, C. (2009), Quantifying photosynthetic capacity and its relationship to leaf 
!! nitrogen content for global-scale terrestrial biosphere models, Global Change Biology, 15, 976-991, 
!! DOI 10.1111/j.1365-2486.2008.01744.x
!!
!_ ===================================================================================================

!     INTERFACE
!     ---------
!
!   INPUT PARAMETERS
!   ----------------
!   KIDIA        Begin point in arrays                                   
!   KFDIA        End point in arrays                                     
!   KLON         Length of arrays                                        
!   LDLAND       Indicator of land grid point
!   KVTYPE       VEGETATION TYPE CORRESPONDING TO TILE                   
!   KCO2TYP      TYPE OF PHOTOSYNTHETIC PATHWAY FOR LOW VEGETATION(C3/C4) (INDEX 1/2)
!   
!   PSRFD        DOWNWARD SHORT WAVE RADIATION FLUX AT SURFACE W/M2      	  					 
!   PAPHM        SURFACE PRESSURE 			       PA        
!   PQS          SATURATION Q AT SURFACE			           KG/KG
!   PQSURF       SPECIFIC HUMIDITY NEAR THE SURFACE            KG/KG
!   PQM          SPECIFIC HUMIDITY                             KG/KG     
!   PTM          TEMPERATURE                                   K         
!   PTSKM        SURFACE (SKIN) TEMPERATURE                    K         
!   PTSOIL       SOIL TEMPERATURE LEVEL 3 (28 - 100cm)         K         
!   PF2	         SOIL MOISTURE STRESS FUNCTION 	               -            
!   PCO2         atmospheric CO2 concentration (kgCO2 kgAir-1)    
!   PLAI         LEAF AREA INDEX                               M2/M2     


!
!   OUTPUT PARAMETERS
!   -----------------
!   PGSTOT      CANOPY CONDUCTANCE TO H20 (CUTICULAR AND STOMATAL) M/S
!   PGPP        GROSS PRIMARY PRODUCTION                           KG_CO2 KG_AIR-1 M S-1
!   PRD         DARK AND DAY RESPIRATION OF CO2 BY LEAVES          KG_CO2 KG_AIR-1 M S-1 
!   PAN         NET CO2 ASSIMILATION OVER CANOPY                   KG_CO2/M2/S
!
!
!   MATCHING TABLE WITH ORCHIDEE PFTS
!   ---------------------------------
!   (1)  ! Crops, Mixed Farming          => PFT13 (C4 crops)
!   (2)  ! Short Grass                   => PFT10 (C3 grass)
!   (3)  ! Evergreen Needleleaf Trees    => PFT7 (Boreal Needleleaf Evergreen)
!   (4)  ! Deciduous Needleleaf Trees    => PFT9 (Boreal Needleleaf Deciduous)
!   (5)  ! Deciduous Broadleaf Trees     => PFT6 (Temperate Broadleaf Summergreen)
!   (6)  ! Evergreen Broadleaf Trees     => PFT5 (Temperate Broadleaf Evergreen)
!   (7)  ! Tall Grass                    => PFT11 (C4 grass)
!   (8)  ! Desert                        => none
!   (9)  ! Tundra                        => PFT10 (C3 grass)
!   (10) ! Irrigated Crops               => PFT12 (C3 crops)
!   (11) ! Semidesert                    => PFT11 (C4 grass)
!   (12) ! Ice Caps and Glaciers         => none
!   (13) ! Bogs and Marshes              => PFT10 (C3 grass)
!   (14) ! Inland Water                  => none
!   (15) ! Ocean                         => none
!   (16) ! Evergreen Shrubs              => PFT11 (C4 grass)
!   (17) ! Deciduous Shrubs              => PFT10 (C3 grass)
!   (18) ! Mixed Forest/woodland         => PFT6 (Temperate Broadleaf Summergreen)
!   (19) ! Interrupted Forest            => PFT6 (Temperate Broadleaf Summergreen)
!   (20) ! Water and Land Mixtures       => none
!
!     REFERENCE.
!     ----------
!
!     MODIFICATIONS
!     -------------
!     V.Bastrikov,F.Maignan,P.Peylin,A.Agusti-Panareda/S.Boussetta Feb 2021 Add Farquhar photosynthesis model



  INTEGER(KIND=JPIM), INTENT(IN)  :: KLON 
  INTEGER(KIND=JPIM), INTENT(IN)  :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN)  :: KFDIA
  LOGICAL           , INTENT(IN)  :: LDLAND(:)
  INTEGER(KIND=JPIM), INTENT(IN)  :: KVTYPE(:)
  INTEGER(KIND=JPIM), INTENT(IN)  :: KCO2TYP(:)
  TYPE(TAGS)        , INTENT(IN)  :: YDAGS
  TYPE(TAGF)        , INTENT(IN)  :: YDAGF
  TYPE(TCST)        , INTENT(IN)  :: YDCST
  REAL(KIND=JPRB),    INTENT(IN)  :: PSRFD(:)  
  REAL(KIND=JPRB),    INTENT(IN)  :: PAPHM(:)                  
  REAL(KIND=JPRB),    INTENT(IN)  :: PQS(:)                
  REAL(KIND=JPRB),    INTENT(IN)  :: PQSURF(:)                                                                 
  REAL(KIND=JPRB),    INTENT(IN)  :: PQM(:)                                                                      
  REAL(KIND=JPRB),    INTENT(IN)  :: PTM(:)                
  REAL(KIND=JPRB),    INTENT(IN)  :: PTSKM(:)               
  REAL(KIND=JPRB),    INTENT(IN)  :: PTSOIL(:) 
  REAL(KIND=JPRB),    INTENT(IN)  :: PCO2(:)      ! CO2 concentration inside the canopy [ppm]
  REAL(KIND=JPRB),    INTENT(IN)  :: PF2(:)       ! Soil moisture stress (0-1,unitless) 
  REAL(KIND=JPRB),    INTENT(IN)  :: PLAI(:)      ! Leaf area index (m2/m2)
  REAL(KIND=JPRB),    INTENT(OUT) :: PGSTOT(:)    ! canopy conductance to H2O.Final unit is (m/s)
  REAL(KIND=JPRB),    INTENT(OUT) :: PGPP(:)      ! Gross primary production(kgCO2 m^{-2} s^{-1})
  REAL(KIND=JPRB),    INTENT(OUT) :: PRD(:)       ! Total Day respiration (respiratory CO2 release other than by photorespiration)
                                                  ! (KgCO2 m−2 s−1). We assume the dark respiration at nighttime 
                                                  ! is the same as the day respiration.
  REAL(KIND=JPRB),    INTENT(OUT) :: PAN(:)       ! Net CO2 assimilation over canopy (kgCO2 m^{-2} s^{-1})

   
!*         0.     LOCAL VARIABLES.
!                 ----- ----------

  INTEGER(KIND=JPIM) :: JL                          ! indices (unitless)
  INTEGER(KIND=JPIM) :: LIMIT_PHOTO, JLEVEL         ! indices (unitless)
  INTEGER(KIND=JPIM) :: ILAI(KLON)                  ! counter for loops on LAI levels (unitless)
  INTEGER(KIND=JPIM) :: KVTYPE_TEMPO(KLON)          ! to deal with KVTYPE=0
  INTEGER(KIND=JPIM) :: NIC, INIC, ICINIC           ! counter/indices (unitless)
  INTEGER(KIND=JPIM) :: INDEX_CALC(KLON)            ! index (unitless)
  INTEGER(KIND=JPIM) :: NIA,INIA,NINA               ! counter/indices (unitless)
  INTEGER(KIND=JPIM) :: ININA,IAINIA                ! counter/indices (unitless)
  INTEGER(KIND=JPIM) :: INDEX_ASSI(KLON)            ! indices (unitless)
  INTEGER(KIND=JPIM) :: INDEX_NON_ASSI(KLON)        ! indices (unitless)
  INTEGER(KIND=JPIM) :: ZCTYPE(KLON)                ! indices for C3/C4 photosynthesis (unitless)
  INTEGER(KIND=JPIM), ALLOCATABLE :: INFO_LIMITPHOTO(:,:) ! ???
  
  LOGICAL            :: LACCLIMATION, LOSM_ACCLIMATION ! Logical switch to enable acclimation of photosynthetic traits  
  LOGICAL            :: LDOWNREGULATION_CO2      ! Logical switch to enable downregulation of CO2 (in atm co2 coupled run) necessary for long climate runs
  LOGICAL            :: LSKIN_TEMP  	         ! Logical to select leaf temperature between skin and surface air temperatures
  LOGICAL            :: LASSIMILATE(KLON)        ! where assimilation is to be calculated (unitless) 
  LOGICAL            :: LCALCULATE(KLON)         ! where assimilation is to be calculated for in the PFTs loop (unitless)
  
  !REAL(KIND=JPRB) :: laisum(KLON)           ! when calculating cim over nlai
  !REAL(KIND=JPRB) :: wind(KLON)             ! Wind (m s^{-1})
!notused  REAL(KIND=JPRB) :: vbeta23(KLON)        ! Beta for fraction of wetted foliage that will transpire (unitless) 
!notused  REAL(KIND=JPRB) :: LEAF_CI_LOWEST(KLON) ! intercellular CO2 concentration at the lowest LAI level (mmol/mol)
!notused  REAL(KIND=JPRB) :: speed                   ! wind speed @tex ($m s^{-1}$) @endtex
!notused  REAL(KIND=JPRB) :: templeafci(KLON,20)  

! REAL(KIND=JPRB) :: ZLEAF_CI(KLON,NLAI)          ! intercellular CO2 concentration (ppm)
 ! REAL(KIND=JPRB) :: ZVC2(KLON,NLAI1)            ! rate of carboxylation (at a specific LAI level) (micromol CO2/m2/s)
 ! REAL(KIND=JPRB) :: ZVJ2(KLON,NLAI1)            ! rate of RubiZSCO regeneration (at a specific LAI (micromol e-/m2/s)
 ! REAL(KIND=JPRB) :: ZASSIMI(KLON,NLAI1)         ! Net assimilation (at a specific LAI level) (micromol/m2/s)
 ! REAL(KIND=JPRB) :: ZGS(KLON,NLAI1)             ! stomatal conductance to CO2 (mol/m2/s)
 ! REAL(KIND=JPRB) :: ZLAITAB(NLAI1)              ! tabulated LAI steps (m2/m2)
 ! REAL(KIND=JPRB) :: ZLIGHT(NLAI)                ! fraction of light that gets through upper LAI levels (0-1,unitless)
 ! REAL(KIND=JPRB) :: ZRD(KLON,NLAI1)             ! Day respiration (respiratory CO2 release other than by photorespiration) (mumol CO2 m−2 s−1)
 ! REAL(KIND=JPRB) :: ZJJ(KLON,NLAI1)             ! Rate of e− transport (umol e− m−2 s−1)
 ! REAL(KIND=JPRB) :: ZCC(KLON,NLAI1)             ! Chloroplast CO2 partial pressure (ubar)
 
  REAL(KIND=JPRB) :: ZTEMP_GROWTH(KLON)       ! Growth temperature (°C) - In offline mode equal to t2m_month
  REAL(KIND=JPRB) :: ZPB(KLON)                   ! Surface pressure [hPa]
  REAL(KIND=JPRB) :: ZCO2(KLON)                  ! CO2 concentration inside the canopy [ppm]
  REAL(KIND=JPRB) :: ZRVEGET(KLON)               ! stomatal resistance of vegetation  (s/m)
  REAL(KIND=JPRB) :: ZVCMAX25(KLON)              ! maximum rate of carboxylation at 25oC (micromol CO2 /m2/s)
  REAL(KIND=JPRB) :: ZJMAX25(KLON)               ! maximum rate of electron transport at 25oC (micromol e-/m2/s)
 
  REAL(KIND=JPRB) :: ZGSTOP(KLON)                ! stomatal conductance to H2O at topmost level (m/s)
  REAL(KIND=JPRB) :: ZGAMMA_STAR(KLON)           ! CO2 compensation point (ppm)
  REAL(KIND=JPRB) :: ZVPD(KLON)                  ! Vapor Pressure Deficit (kPa)
  REAL(KIND=JPRB) :: ZWATER_LIM(KLON)            ! water limitation factor (0-1,unitless)
  REAL(KIND=JPRB) :: ZASSIMTOT(KLON)             ! total apparent assimilation GPP=An+Rd (micromol CO2/m2/s)
  REAL(KIND=JPRB) :: ZRDTOT(KLON)                ! total mitocondrial (day and dark) respiration by leaves (micromol CO2/m2/s)
  REAL(KIND=JPRB) :: ZLEAF_GS_TOP(KLON)          ! leaf stomatal conductance to H2O at topmost level (mol H2O/m2/s)
  REAL(KIND=JPRB) :: ZT_VCMAX                    ! Temperature dependance of Vcmax (unitless)
  REAL(KIND=JPRB) :: ZS_VCMAX_ACCLIMTEMP         ! Entropy term for Vcmax 
  REAL(KIND=JPRB) :: ZT_GAMMA_STAR               ! Temperature dependance of gamma_star (unitless)    
                                                 ! accounting for acclimation to temperature (J K-1 mol-1)
  REAL(KIND=JPRB) :: ZT_JMAX                     ! Temperature dependance of Jmax
  REAL(KIND=JPRB) :: ZS_JMAX_ACCLIMTEMP          ! Entropy term for Jmax accounting for acclimation to temperature (J K-1 mol-1)
  REAL(KIND=JPRB) :: ZT_GM                       ! Temperature dependance of gmw
  REAL(KIND=JPRB) :: ZT_RD                       ! Temperature dependance of Rd (unitless)
  REAL(KIND=JPRB) :: ZT_KMC                      ! Temperature dependance of KmC (unitless)
  REAL(KIND=JPRB) :: ZT_KMO                      ! Temperature dependance of KmO (unitless)
  REAL(KIND=JPRB) :: ZT_SCO                      ! Temperature dependance of ZSCO
  REAL(KIND=JPRB) :: ZVC(KLON)                   ! Maximum rate of RubiZSCO activity-limited carboxylation (mumol CO2 m−2 s−1)
  REAL(KIND=JPRB) :: ZVJ(KLON)                   ! Maximum rate of e- transport under saturated light (mumol CO2 m−2 s−1)
  REAL(KIND=JPRB) :: ZGM(KLON)                   ! Mesophyll diffusion conductance (molCO2 m−2 s−1 bar−1)
  REAL(KIND=JPRB) :: ZG0VAR(KLON)                !  ???? Residual stomatal conductance when irradiance approaches zero (mol CO2 m−2 s−1 bar−1)???? Is this the same as cuticle conductance?
  REAL(KIND=JPRB) :: ZKMC(KLON)                  ! Michaelis–Menten constant of RubiZSCO for CO2 (mubar)
  REAL(KIND=JPRB) :: ZKMO(KLON)                  ! Michaelis–Menten constant of RubiZSCO for O2 (mubar)
  REAL(KIND=JPRB) :: ZSCO(KLON)                  ! Relative CO2 /O2 specificity factor for Rubisco (bar bar-1)
  REAL(KIND=JPRB) :: ZGBCO2(KLON)                ! Boundary-layer conductance (molCO2 m−2 s−1 bar−1)
  REAL(KIND=JPRB) :: ZGBH2O(KLON)                ! Boundary-layer conductance (molH2O m−2 s−1 bar−1)
  REAL(KIND=JPRB) :: ZFVPD(KLON)                 ! Factor for describing the effect of leaf-to-air vapour difference on gs (-)
  REAL(KIND=JPRB) :: ZLOW_GAMMA_STAR(KLON)       ! Half of the reciprocal of Sc/o (bar bar-1)
  REAL(KIND=JPRB) :: ZNITRO_VCMAX                ! Nitrogen level dependance of Vcmacx and Jmax 
  REAL(KIND=JPRB) :: ZIABS(KLON)                 ! Photon flux density absorbed by leaf photosynthetic pigments (umol photon m−2 s−1)
  REAL(KIND=JPRB) :: ZJMAX(KLON)                 ! Maximum value of J under saturated light (umol e− m−2 s−1)
  REAL(KIND=JPRB) :: ZLEAF_TEMP(KLON)            ! ZLEAF_TEMPerature (????K)


  REAL(KIND=JPRD) :: ZFCYC                       ! Fraction of electrons at PSI that follow cyclic transport around PSI (-)
  REAL(KIND=JPRD) :: ZZ                          ! A lumped parameter (see Yin et al. 2009) ( mol mol-1)     
  REAL(KIND=JPRD) :: ZRM                         ! Day respiration in the mesophyll (umol CO2 m−2 s−1)
  REAL(KIND=JPRD) :: ZCS_STAR                    ! Cs -based CO2 compensation point in the absence of Rd (ubar)
  REAL(KIND=JPRD) :: ZJ2                         ! Rate of all e− transport through PSII (umol e− m−2 s−1)
  REAL(KIND=JPRD) :: ZVPJ2                       ! e− transport-limited PEP carboxylation rate (umol CO2 m−2 s−1)
  REAL(KIND=JPRD) :: ZA_1, ZA_3                  ! Lowest First and third roots of the analytical solution for a general cubic equation 
                                                 ! (see Appendix A of Yin et al. 2009) (umol CO2 m−2 s−1)
  REAL(KIND=JPRD) :: ZA_1_tmp, ZA_3_tmp          ! Temporary First and third roots of the analytical solution for a general cubic equation 
                                                 ! (see Appendix A of Yin et al. 2009) (umol CO2 m−2 s−1)
  REAL(KIND=JPRD) :: ZLEAF_CI                    ! intercellular CO2 concentration (ppm)
  REAL(KIND=JPRD) :: ZCC                         ! Chloroplast CO2 partial pressure (ubar)
  REAL(KIND=JPRD) :: ZOBS                        ! Bundle-sheath oxygen partial pressure (ubar)
  REAL(KIND=JPRD) :: ZCI_STAR                    ! Ci-based CO2 compensation point in the absence of Rd (ubar)    
  REAL(KIND=JPRD) :: ZA,ZB,ZC,ZD,ZM,ZF,ZJ,ZG,ZH,ZI,ZL,ZP,ZQ,ZR ! Variables used for solving the cubic equation (see Yin et al. (2009))
  REAL(KIND=JPRD) :: ZQQ,ZUU,ZPSI,ZX1,ZX2,ZX3    ! Variables used for solving the cubic equation (see Yin et al. (2009))
!notused  REAL(KIND=JPRB) :: ZCRESIST                    ! coefficient for resistances (????)
  

  REAL(KIND=JPRB), ALLOCATABLE :: ZVC2(:,:)      ! rate of carboxylation (at a specific LAI level) (micromol CO2/m2/s)
  REAL(KIND=JPRB), ALLOCATABLE :: ZVJ2(:,:)      ! rate of RubiZSCO regeneration (at a specific LAI (micromol e-/m2/s)
  REAL(KIND=JPRB), ALLOCATABLE :: ZASSIMI(:,:)   ! net assimilation (at a specific LAI level) (micromol/m2/s)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGS(:,:)       ! stomatal conductance to CO2 (mol/m2/s)
  REAL(KIND=JPRB), ALLOCATABLE :: ZLAITAB(:)     ! tabulated LAI steps (m2/m2)
  REAL(KIND=JPRB), ALLOCATABLE :: ZLIGHT(:)      ! fraction of light that gets through upper LAI levels (0-1,unitless)
  REAL(KIND=JPRB), ALLOCATABLE :: ZRD(:,:)       ! Day respiration (respiratory CO2 release other than by photorespiration) (mumol CO2 m−2 s−1)
  REAL(KIND=JPRB), ALLOCATABLE :: ZJJ(:,:)       ! Rate of e− transport (umol e− m−2 s−1)

!ORIG  REAL(KIND=JPRB), ALLOCATABLE :: ZLEAF_CI(:,:)  ! intercellular CO2 concentration (ppm)
!ORIG  REAL(KIND=JPRB), ALLOCATABLE :: ZCC(:,:)       ! Chloroplast CO2 partial pressure (ubar)

!

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('FARQUHAR_MOD:FARQUHAR',0,ZHOOK_HANDLE)
  
  ASSOCIATE(NLAI=>YDAGF%NLAI, NLAI1=>YDAGF%NLAI1, NPATH=>YDAGF%NPATH, &
   & RDOWNREGULATION_CO2_BASELEVEL=>YDAGF%RDOWNREGULATION_CO2_BASELEVEL, &
   & RDOWNREGULATION_CO2_MINIMUM=>YDAGF%RDOWNREGULATION_CO2_MINIMUM, &
   & RUNDEF_SECHIBA=>YDAGF%RUNDEF_SECHIBA, RUNDEF=>YDAGF%RUNDEF, NIUNDEF=>YDAGF%NIUNDEF,&
   & RONE_MONTH=>YDAGF%RONE_MONTH,  &
   & RA1=>YDAGF%RA1,RB1=>YDAGF%RB1, RARJV=>YDAGF%RARJV, RBRJV=>YDAGF%RBRJV, &
   & RASJ=>YDAGF%RASJ, RBSJ=>YDAGF%RBSJ, RASV=>YDAGF%RASV, RBSV=>YDAGF%RBSV, &
   & RD_VCMAX=>YDAGF%RD_VCMAX, RD_JMAX=>YDAGF%RD_JMAX, RE_GAMMA_STAR=>YDAGF%RE_GAMMA_STAR, &
   & RE_GM=>YDAGF%RE_GM, RS_GM=>YDAGF%RS_GM, RS_JMAX=>YDAGF%RS_JMAX,&
   & RD_GM=>YDAGF%RD_GM, RE_VCMAX=>YDAGF%RE_VCMAX, RE_JMAX=>YDAGF%RE_JMAX, RE_KMC=>YDAGF%RE_KMC, &
   & RE_KMO=>YDAGF%RE_KMO, RE_RD=>YDAGF%RE_RD, &
   & RE_SCO=>YDAGF%RE_SCO, RG0=>YDAGF%RG0, RGAMMA_STAR25=>YDAGF%RGAMMA_STAR25, RGM25=>YDAGF%RGM25,&
   & RKMC25=>YDAGF%RKMC25, RKMO25=>YDAGF%RKMO25, &
   & RSCO25=>YDAGF%RSCO25, RSTRESS_GM=>YDAGF%RSTRESS_GM, RSTRESS_GS=>YDAGF%RSTRESS_GS, &
   & RSTRESS_VCMAX=>YDAGF%RSTRESS_VCMAX, &
   & RTPHOTO_MAX=>YDAGF%RTPHOTO_MAX, RTPHOTO_MIN=>YDAGF%RTPHOTO_MIN, RVCMAX25=>YDAGF%RVCMAX25, &
   & RHUMREL=>YDAGF%RHUMREL, &
   & RDOWNREGULATION_CO2_COEF=>YDAGF%RDOWNREGULATION_CO2_COEF, RSTRUCT_CONST=>YDAGF%RSTRUCT_CONST, &
   & RLAI_LEVEL_DEPTH=>YDAGF%RLAI_LEVEL_DEPTH, RLAIMAX=>YDAGF%RLAIMAX, &
   & REXT_COEFF=>YDAGF%REXT_COEFF, RGB_REF=>YDAGF%RGB_REF, &
   & RTP_00=>YDAGF%RTP_00, RPB_STD=>YDAGF%RPB_STD, RRG_TO_PAR=>YDAGF%RRG_TO_PAR, &
   & RW_TO_MOL=>YDAGF%RW_TO_MOL, RALPHA_LL=>YDAGF%RALPHA_LL, &
   & RTHETA=>YDAGF%RTHETA, RFPSIR=>YDAGF%RFPSIR, RFQ=>YDAGF%RFQ, RFPSEUDO=>YDAGF%RFPSEUDO,&
   & RH_PROTONS=>YDAGF%RH_PROTONS, RGBS=>YDAGF%RGBS, &
   & ROI=>YDAGF%ROI, RALPHA=>YDAGF%RALPHA, RKP=>YDAGF%RKP, RTETENS_1=>YDAGF%RTETENS_1, &
   & RTETENS_2=>YDAGF%RTETENS_2, RREF_TEMP=>YDAGF%RREF_TEMP, &
   & RMOL_TO_M_1=>YDAGF%RMOL_TO_M_1, RRATIO_H2O_TO_CO2=>YDAGF%RRATIO_H2O_TO_CO2, &
   & RN_VERT_ATT=>YDAGF%RN_VERT_ATT, &
   & RTAU_TEMP_AIR_MONTH=>YDAGF%RTAU_TEMP_AIR_MONTH, RPI=>YDCST%RPI, &
   & RMAIR=>YDAGS%RMAIR,RMCO2=>YDAGS%RMCO2,RRDCF=>YDAGS%RRDCF) 
  
 
! set default configuration of Farquhar model
!

  LACCLIMATION = .TRUE.
  LOSM_ACCLIMATION = .FALSE. ! This can only work in the offline model (acclimation as implemented in ORCHIDEE, default is false)
  LDOWNREGULATION_CO2 = .FALSE.    ! Set to .TRUE. if you want CO2 downregulation. Note this is mainly relevant for climate change runs with doubling of CO2
  LSKIN_TEMP = .FALSE.          ! Leaf temperature is equal to skin temperature (else to surface air temperature)


! initialization of global variables

DO JL=KIDIA,KFDIA

  PGPP(JL) = 0._JPRB
  PRD(JL) = 0._JPRB
  PAN(JL) = 0._JPRB
  PGSTOT(JL) = 0._JPRB

ENDDO
  
! initialization of local variables

   IF (.NOT.ALLOCATED(INFO_LIMITPHOTO)) ALLOCATE (INFO_LIMITPHOTO(KLON,NLAI1)) 
   IF (.NOT.ALLOCATED(ZVC2))            ALLOCATE (ZVC2(KLON,NLAI1)) 
   IF (.NOT.ALLOCATED(ZVJ2))            ALLOCATE (ZVJ2(KLON,NLAI1)) 
   IF (.NOT.ALLOCATED(ZASSIMI))         ALLOCATE (ZASSIMI(KLON,NLAI1)) 
   IF (.NOT.ALLOCATED(ZGS))             ALLOCATE (ZGS(KLON,NLAI1)) 
   IF (.NOT.ALLOCATED(ZLAITAB))         ALLOCATE (ZLAITAB(NLAI1)) 
   IF (.NOT.ALLOCATED(ZLIGHT))          ALLOCATE (ZLIGHT(NLAI)) 
   IF (.NOT.ALLOCATED(ZRD))             ALLOCATE (ZRD(KLON,NLAI1)) 
   IF (.NOT.ALLOCATED(ZJJ))             ALLOCATE (ZJJ(KLON,NLAI1)) 
!ORIG   IF (.NOT.ALLOCATED(ZLEAF_CI))        ALLOCATE (ZLEAF_CI(KLON,NLAI)) 
!ORIG   IF (.NOT.ALLOCATED(ZCC))             ALLOCATE (ZCC(KLON,NLAI1)) 
   
 
DO JL=KIDIA,KFDIA

  ZVCMAX25(JL)=0._JPRB
  ZGAMMA_STAR(JL)=0._JPRB
  ZKMO(JL)=0._JPRB
  ZKMC(JL)=0._JPRB
  ZGM(JL)=0._JPRB
  ZG0VAR(JL) =0._JPRB
  ZGSTOP(JL) = 0._JPRB
  ZASSIMTOT(JL) = 0._JPRB
  ZRDTOT(JL) = 0._JPRB
  ZLEAF_GS_TOP(JL) = 0._JPRB
!  PGSMEAN(JL) = 0._JPRB

  ! Change units to hPa
  ZPB(JL) = PAPHM(JL)/100.0
  ! Change units to ppm
  ZCO2(JL) = REAL(REAL(PCO2(JL),KIND=JPRD)*1.E6_JPRD*RMAIR/RMCO2 ,KIND=JPRB)
  IF ((ZCO2(JL) > 600.0_JPRB) .OR. (ZCO2(JL) < 200.0_JPRB)) THEN
    WRITE(*,*) 'farquhar : PCO2(JL), ZCO2 = ', PCO2(JL), ZCO2(JL)
  ENDIF
  ZVC(JL)=0._JPRB
  ZVJ(JL)=0._JPRB
  LCALCULATE(JL) = .FALSE.

  ZCTYPE=1
  IF (KCO2TYP(JL).EQ.4) THEN
      ZCTYPE(JL)=2
  ENDIF

  DO JLEVEL = 1, NLAI+1  
!orig    ZCC(JL,JLEVEL)=0._JPRB
    ZVC2(JL,JLEVEL)=0._JPRB
    ZVJ2(JL,JLEVEL)=0._JPRB
    ZJJ(JL,JLEVEL)=0._JPRB
    INFO_LIMITPHOTO(JL,JLEVEL)=0
    ZGS(JL,JLEVEL)=0._JPRB
!notused    templeafci(JL,JLEVEL)=0._JPRB
    ZASSIMI(JL,JLEVEL)=0._JPRB
    ZRD(JL,JLEVEL)=0._JPRB
  ENDDO

  IF (LSKIN_TEMP) THEN
    ZLEAF_TEMP(JL)=PTSKM(JL)
  ELSE
    ZLEAF_TEMP(JL)=PTM(JL)
  ENDIF

  ! Originally: If acclimation is on (offline mode), ZTEMP_GROWTH is updated only once per time step 
  ! to get a monthly mean. In online mode the acclimation this is not possible, so 
  ! tempgrowth is taken as the temperature from the 3rd soil layer as a proxy from mean 20-30 day temperature
  ! from model level closest to surface

!   IF (LOSM_ACCLIMATION) THEN
!   
!  !   ZTEMP_GROWTH(JL) = ZLEAF_TEMP(JL) - RTP_00
! 
!    ! We use temp_air (Kattge & Knorr, 2007: "the plants’ average ambient temperature").
!      IF (KSTEP .NE. KSTEP_SAV) THEN
!        IF (KSTEP .LT. RONE_MONTH/PTSTEP) THEN
!          !Simple mean while lower than one month
!            ZTEMP_GROWTH(JL) = ( ZTEMP_GROWTH(JL) * KSTEP + PTM(JL) - RTP_00) / (KSTEP+1)
!        ELSE
!         !Weighted average
!            ZTEMP_GROWTH(JL) = ( ZTEMP_GROWTH(JL) * ( TAU_TEMP_AIR_MONTH - PTSTEP ) + &
!              &                  (PTM(JL) - RTP_00) * PTSTEP ) / TAU_TEMP_AIR_MONTH     
!        ENDIF
! 
!         WHERE ( ABS(ZTEMP_GROWTH(JL) + TP_00) .LT. EPSILON(0._JPRM) )
!           ZTEMP_GROWTH(JL) = -TP_00
!         ENDWHERE
!       ENDIF
! 
!  ELSE

    ZTEMP_GROWTH(JL) = PTSOIL(JL) - RTP_00
      
!   ENDIF ! osm_acclimation

ENDDO !JL


    ! 1. Preliminary calculations
! =========================================================================================

  !
  ! 1.1 Calculate LAI steps
  ! The integration at the canopy level is done over nlai fixed levels.
  DO JLEVEL = 1, NLAI+1
    ZLAITAB(JLEVEL) = RLAIMAX*(EXP(RLAI_LEVEL_DEPTH*REAL(JLEVEL-1._JPRB))-1._JPRB)/ &
       &                      (EXP(RLAI_LEVEL_DEPTH*REAL(NLAI,JPRB))-1._JPRB)
  ENDDO
  !
  ! 1.2 Calculate light fraction for each LAI step
  ! The available light follows a simple Beer extinction law. 
  ! The extinction coefficients (ext_coef) are PFT-dependant constants.

  DO JLEVEL = 1, NLAI
    ZLIGHT(JLEVEL) = exp( -REXT_COEFF*ZLAITAB(JLEVEL) )
  ENDDO 
 
  ! Here we set a protection for when KVTYPE(JL) is zero. We set the KVTYPE_TEMPO array to one 
  ! in this case to be able to index the arrays that depend on KVTPE
  ! and start at 1. A protection against KVTYPE being zero is explicitely done is the loops 
  ! below to ensure that there are no caculations done in this case.
  WHERE(KVTYPE(KIDIA:KFDIA) .GT. 0_JPIM ) 
    KVTYPE_TEMPO(KIDIA:KFDIA)=KVTYPE(KIDIA:KFDIA)
  ELSEWHERE
    KVTYPE_TEMPO(KIDIA:KFDIA)=1_JPIM
  ENDWHERE
     

  IF (LDOWNREGULATION_CO2) THEN
    DO JL=KIDIA,KFDIA    
      IF (LDLAND(JL)) THEN
          ZVCMAX25(JL) = RVCMAX25(KVTYPE_TEMPO(JL),ZCTYPE(JL))* &
   &                   (1._JPRB-RDOWNREGULATION_CO2_COEF(KVTYPE_TEMPO(JL),ZCTYPE(JL)) * &
   &	               (MAX(ZCO2(JL),RDOWNREGULATION_CO2_MINIMUM)-RDOWNREGULATION_CO2_BASELEVEL) / &
   &	               (MAX(ZCO2(JL),RDOWNREGULATION_CO2_MINIMUM)+20._JPRB))
      ENDIF
    ENDDO  
  ELSE
    DO JL=KIDIA,KFDIA    
      IF (LDLAND(JL)) THEN
          ZVCMAX25(JL) = RVCMAX25(KVTYPE_TEMPO(JL),ZCTYPE(JL))
      ENDIF
    ENDDO
  ENDIF

  !
  ! 1.3 Estimate vapour pressure deficit of air (for calculation of the stomatal conductance).
  !
  !Note: qsat is passed as input (qsurf or PQSURF)
  !orig    CALL qsatcalc (KIDIA, KFDIA, KLON, YDAGF, ZLEAF_TEMP, ZPB, PQSURF)


  DO JL=KIDIA,KFDIA
 
    ZVPD(JL) = (PQS(JL)*ZPB(JL) / (RTETENS_1+PQS(JL)*RTETENS_2 ) ) &
       - ( PQM(JL) * ZPB(JL) / (RTETENS_1+PQM(JL)* RTETENS_2) )


    ! VPD is needed in kPa
    ZVPD(JL) = ZVPD(JL)/10._JPRB

    !
    ! structural resistance
    !
    !TO BE DEpALT WITH FOR BARE SOIL ??
!!    ZRSTRUCT(:,1) = RSTRUCT_CONST(1)
    !
    ! stomatal resistance
    !
    ZRVEGET(JL) = RUNDEF_SECHIBA
!    ZRSTRUCT(JL) = RSTRUCT_CONST(KVTYPE_TEMPO(JL),ZCTYPE(JL))
    !
    ! 2. beta coefficient for vegetation transpiration
    !
    !
    !vbeta3(:,:) = 0._JPRB
    !vbeta3pot(:,:) = 0._JPRB

  ENDDO !JL

  NIA=0
  NINA=0
  !
  DO JL=KIDIA,KFDIA

!    WRITE(*,*) 'farquhar : LDLAND/PLAI/KVTYPE = ', LDLAND(JL), PLAI(JL), KVTYPE(JL)
    IF ( (LDLAND(JL)) .AND. ( PLAI(JL) .GT. 0.01_JPRB ) .AND. (KVTYPE(JL) .GT. 0_JPIM)) THEN


       !! Calculates the water limitation factor.
!               ZWATER_LIM(JL) = PF2(JL)
        ZWATER_LIM(JL) = MIN(1._JPRB,PF2(JL)*RHUMREL(KVTYPE_TEMPO(JL),ZCTYPE(JL)))


       ! Day respiration (assumed to be equal to dark respiration at nighttime)
       ! Parametrization of Yin et al. (2009) - from Bernacchi et al. (2001)

        ZT_RD = ARRHENIUS_FN(ZLEAF_TEMP(JL),RREF_TEMP,RE_RD,YDAGF,YDCST)
        DO JLEVEL = 1, NLAI
             ! see Comment in legend of Fig. 6 of Yin et al. (2009)
             ! Rd25 is assumed to equal 0.01 Vcmax25 
             ZNITRO_VCMAX = ( 1._JPRB - RN_VERT_ATT * ( 1._JPRB - ZLIGHT(JLEVEL) ) )
             ZRD(JL,JLEVEL) = ZVCMAX25(JL) * ZNITRO_VCMAX * 0.01_JPRB * ZT_RD &
                                &           * MAX(1._JPRB-RSTRESS_VCMAX, ZWATER_LIM(JL))
        ENDDO

            
       IF ( ( PSRFD(JL) .GT. EPSILON(1._JPRM) )   .AND. &
                 ( PF2(JL) .GT. EPSILON(1._JPRM) ) .AND. &
                ( ZTEMP_GROWTH(JL) .GT. RTPHOTO_MIN ) .AND. &
                 ( ZTEMP_GROWTH(JL) .LT. RTPHOTO_MAX ) ) THEN

              LASSIMILATE(JL) = .TRUE.
              NIA=NIA+1
              INDEX_ASSI(NIA)=JL          
              ILAI(JL) = 1

              ! Here is the calculation of ZASSIMIlation and stomatal conductance
              ! based on the work of Farquahr, von Caemmerer and Berry (FvCB model) 
              ! as described in Yin et al. 2009
              ! Yin et al. developed a extended version of the FvCB model for C4 plants
              ! and proposed an analytical solution for both photosynthesis pathways (C3 and C4)
              ! Photosynthetic parameters used are those reported in Yin et al. 
              ! Except For Vcmax25, relationships between Vcmax25 and Jmax25 for which we use 
              ! Medlyn et al. (2002) and Kattge & Knorr (2007)
              ! Because these 2 references do not consider mesophyll conductance, we neglect this term
              ! in the formulations developed by Yin et al. 
              ! Consequently, gm (the mesophyll conductance) tends to the infinite
              ! This is of importance because as stated by Kattge & Knorr and Medlyn et al.,
              ! values of Vcmax and Jmax derived with different model parametrizations are not 
              ! directly comparable and the published values of Vcmax and Jmax had to be standardized
              ! to one consistent formulation and parametrization
        
              ! See eq. 6 of Yin et al. (2009)
              ! Parametrization of Medlyn et al. (2002) - from Bernacchi et al. (2001)
              ZT_KMC        = ARRHENIUS_FN(ZLEAF_TEMP(JL),RREF_TEMP,RE_KMC,YDAGF,YDCST)
              ZT_KMO        = ARRHENIUS_FN(ZLEAF_TEMP(JL),RREF_TEMP,RE_KMO,YDAGF,YDCST)
              ZT_SCO        = ARRHENIUS_FN(ZLEAF_TEMP(JL),RREF_TEMP,RE_SCO,YDAGF,YDCST)
              ZT_GAMMA_STAR = ARRHENIUS_FN(ZLEAF_TEMP(JL),RREF_TEMP,RE_GAMMA_STAR,YDAGF,YDCST)
        

              ! For C3 plants, we assume that the Entropy term for Vcmax and Jmax 
              ! acclimates to temperature as shown by Kattge & Knorr (2007) - Eq. 9 and 10
              ! and that Jmax and Vcmax respond to temperature following a modified Arrhenius function
              ! (with a decrease of these parameters for high temperature) as in Medlyn et al. (2002) 
              ! and Kattge & Knorr (2007).
              ! In Yin et al. (2009), temperature dependance to Vcmax is based only on a Arrhenius function
              ! Concerning this apparent unconsistency, have a look to the section 'Limitation of 
              ! Photosynthesis by gm' of Bernacchi (2002) that may provide an explanation
              
              ! Growth temperature tested by Kattge & Knorr range from 11 to 35°C
              ! So, we limit the relationship between these lower and upper limits


              ! As shown by Kattge & Knorr (2007), we make use
              ! of Jmax25/Vcmax25 ratio (rJV) that acclimates to temperature for C3 plants
              ! rJV is written as a function of the growth temperature
              ! rJV = ARJV + BRJV * T_month 
              ! See eq. 10 of Kattge & Knorr (2007)
              ! and Table 3 for Values of ARJV anf BRJV 
              ! Growth temperature is monthly temperature (expressed in °C) - See first paragraph of
              ! section Methods/Data of Kattge & Knorr


              IF (ZCTYPE(JL) .EQ. 1_JPIM) THEN
              !C3
                IF (LACCLIMATION) THEN      
                   ! acclimation of entropy term for Vcmax25 and Jmax25 modified Arrhenius fn
                   ! and acclimation of ratio Jmax25/Vcmax25

                  ZS_JMAX_ACCLIMTEMP = RASJ(1) + RBSJ(1) * MAX(11._JPRB, MIN(ZTEMP_GROWTH(JL),35._JPRB))      
             	  ZT_JMAX = ARRHENIUS_MODIFIED_FN(ZLEAF_TEMP(JL),RREF_TEMP,&
                                    &              RE_JMAX(KVTYPE_TEMPO(JL),1),RD_JMAX(1), &
                                    &              ZS_JMAX_ACCLIMTEMP,YDAGF,YDCST)

          	  ZS_VCMAX_ACCLIMTEMP = RASV(1) + RBSV(1) * MAX(11._JPRB, MIN(ZTEMP_GROWTH(JL),35._JPRB))   
     	          ZT_VCMAX = ARRHENIUS_MODIFIED_FN(ZLEAF_TEMP(JL),RREF_TEMP,RE_VCMAX(KVTYPE_TEMPO(JL),1), &
                                    &              RD_VCMAX(1),ZS_VCMAX_ACCLIMTEMP,YDAGF,YDCST)


                  ZJMAX25(JL)=(RARJV(1)+ RBRJV(1)*MAX(11._JPRB,MIN(ZTEMP_GROWTH(JL),35._JPRB)))*ZVCMAX25(JL)

                ELSE 
                  ! No acclimation: 

                  ! Use standard entropy term S value independent of temperature is not available 
                  ! from Yin and Struik (2009) table 2
             	  ZT_JMAX = ARRHENIUS_MODIFIED_FN(ZLEAF_TEMP(JL),RREF_TEMP,RE_JMAX(KVTYPE_TEMPO(JL),1),&
                                    &               RD_JMAX(1), RS_JMAX(1),YDAGF,YDCST)

                  ZT_VCMAX = ARRHENIUS_FN(ZLEAF_TEMP(JL),RREF_TEMP,RE_VCMAX(KVTYPE_TEMPO(JL),1),YDAGF,YDCST)

                  ZJMAX25(JL)=(RARJV(1)+ RBRJV(1)*25._JPRB)*ZVCMAX25(JL)

                ENDIF

                ZVJ(JL) = ZJMAX25(JL) * ZT_JMAX
                ZVC(JL) = ZVCMAX25(JL) * ZT_VCMAX


                ZT_GM = ARRHENIUS_MODIFIED_FN(ZLEAF_TEMP(JL),RREF_TEMP,RE_GM(1),RD_GM(1),RS_GM(1),YDAGF,YDCST)
                ZGM(JL) = RGM25(KVTYPE_TEMPO(JL),1) * ZT_GM * MAX(1._JPRB-RSTRESS_GM, ZWATER_LIM(JL))

                ZG0VAR(JL) = RG0(KVTYPE_TEMPO(JL),1)* MAX(1._JPRB-RSTRESS_GS, ZWATER_LIM(JL))
                ZKMC(JL)=RKMC25(1)*ZT_KMC
                ZKMO(JL)=RKMO25(1)*ZT_KMO
                ZSCO(JL)=RSCO25(1)*ZT_SCO
                ! VPD expressed in kPa
                ! Note : MIN(1.-min_sechiba,MAX(min_sechiba,(A1(JVT) - B1(JVT) * VPD(:)))) 
                ! is always between 0-1 not including 0 and 1
                ZFVPD(JL) = 1._JPRM / ( 1._JPRM / MIN(1._JPRM-EPSILON(1._JPRM), &
                         &    MAX(EPSILON(1._JPRM),(RA1(KVTYPE_TEMPO(JL),1) - RB1(KVTYPE_TEMPO(JL),1) * &
                         &      ZVPD(JL)))) - 1._JPRB ) &
                         &     * MAX(1._JPRM-RSTRESS_GS, ZWATER_LIM(JL))

	      ELSE IF (ZCTYPE(JL) .EQ. 2_JPIM) THEN
      	      !C4
                IF (LACCLIMATION) THEN      
     	          ZS_JMAX_ACCLIMTEMP = RASJ(2) + RBSJ(2) * MAX(11., MIN(ZTEMP_GROWTH(JL),35.))   
             	  ZT_JMAX = ARRHENIUS_MODIFIED_FN(ZLEAF_TEMP(JL),  RREF_TEMP,RE_JMAX(KVTYPE_TEMPO(JL),2),&
                          &                    RD_JMAX(2),ZS_JMAX_ACCLIMTEMP,YDAGF,YDCST)
   
        	  ZS_VCMAX_ACCLIMTEMP = RASV(2) + RBSV(2) * MAX(11., MIN(ZTEMP_GROWTH(JL),35.))  
          	  ZT_VCMAX = ARRHENIUS_MODIFIED_FN(ZLEAF_TEMP(JL),RREF_TEMP,&
                          &                      RE_VCMAX(KVTYPE_TEMPO(JL),2),RD_VCMAX(2),&
                          &                      ZS_VCMAX_ACCLIMTEMP,YDAGF,YDCST)

                  ZVJ(JL)=(RARJV(2)+RBRJV(2)* MAX(11.,MIN(ZTEMP_GROWTH(JL),35.)))*ZVCMAX25(JL)*ZT_JMAX
                  ZVC(JL) = ZVCMAX25(JL) * ZT_VCMAX

                ELSE

                  ZJMAX25(JL)=(RARJV(2)+ RBRJV(2)*25._JPRB)*ZVCMAX25(JL)


                  ! Use standard entropy term S value independent of temperature is not available 
                  ! from Yin and Struik (2009) table 2
             	  ZVJ(JL) = ZJMAX25(JL) * ARRHENIUS_MODIFIED_FN(ZLEAF_TEMP(JL),RREF_TEMP, &
                           & RE_JMAX(KVTYPE_TEMPO(JL),2),RD_JMAX(2), RS_JMAX(2),YDAGF,YDCST)

                  ZVC(JL) = ZVCMAX25(JL) * ARRHENIUS_FN(ZLEAF_TEMP(JL),RREF_TEMP,&
                           &                            RE_VCMAX(KVTYPE_TEMPO(JL),2),YDAGF,YDCST)

                ENDIF

               ! ZT_GM = ARRHENIUS_MODIFIED_FN(ZLEAF_TEMP(JL),RREF_TEMP,RE_GM(2),RD_GM(2),RS_GM(2),YDAGF,YDCST)
               !gm not neded for C4 plants
               !ZGM(JL) = RGM25(2) * ZT_GM * MAX(1-RSTRESS_GM, ZWATER_LIM(:))
               ZG0VAR(JL) = RG0(KVTYPE_TEMPO(JL),2)* MAX(1._JPRB-RSTRESS_GS, ZWATER_LIM(JL))
               ZKMC(JL)=RKMC25(2)*ZT_KMC
               ZKMO(JL)=RKMO25(2)*ZT_KMO
               ZSCO(JL)=RSCO25(2)*ZT_SCO
              ! VPD expressed in kPa
              ! Note : MIN(1.-min_sechiba,MAX(min_sechiba,(A1(JVT) - B1(JVT) * VPD(:)))) 
              ! is always between 0-1 not including 0 and 1
               ZFVPD(JL) = 1._JPRM / ( 1._JPRM / MIN(1._JPRM-EPSILON(1._JPRM), &
                         &      MAX(EPSILON(1._JPRM),(RA1(KVTYPE_TEMPO(JL),2) - RB1(KVTYPE_TEMPO(JL),2) &
                         &       * ZVPD(JL)))) - 1._JPRM ) &
                         &      * MAX(1._JPRM-RSTRESS_GS, ZWATER_LIM(JL))
              ELSE

  	       ZSCO(JL)=RSCO25(2)

              ENDIF !C3/C4


              ZGAMMA_STAR(JL) = RGAMMA_STAR25*ZT_GAMMA_STAR

              ! low_gamma_star is defined by Yin et al. (2009)
              ! as the half of the reciprocal of Sco - See Table 2
              ZLOW_GAMMA_STAR(JL) = 0.5_JPRB / ZSCO(JL)

              ! leaf boundary layer conductance 
              ! conversion from a conductance in (m s-1) to (mol H2O m-2 s-1)
              ! from Pearcy et al. (1991, see below)
              ZGBH2O(JL) = RGB_REF * 44.6_JPRB * (RTP_00/PTM(JL)) * (ZPB(JL)/RPB_STD) 

              ! conversion from (mol H2O m-2 s-1) to (mol CO2 m-2 s-1)
              ZGBCO2(JL) = ZGBH2O(JL) / RRATIO_H2O_TO_CO2

              DO JLEVEL = 1, NLAI
                 ! 2.4.1 Vmax is scaled into the canopy due to reduction of nitrogen 
                 !! (Johnson and Thornley,1984).

                  ZNITRO_VCMAX = ( 1._JPRB - RN_VERT_ATT * ( 1._JPRB - ZLIGHT(JLEVEL) ) )
                  ZVC2(JL,JLEVEL) = ZVC(JL) * ZNITRO_VCMAX * MAX(1._JPRB-RSTRESS_VCMAX, ZWATER_LIM(JL))
                  ZVJ2(JL,JLEVEL) = ZVJ(JL) * ZNITRO_VCMAX * MAX(1._JPRB-RSTRESS_VCMAX, ZWATER_LIM(JL))

!These lines have been moved outside the assim check because Rd needs to be computed at nighttime as well
!orig                 ! see Comment in legend of Fig. 6 of Yin et al. (2009)
!orig                 ! Rd25 is assumed to equal 0.01 Vcmax25 
!orig                 ZRD(JL,JLEVEL) = ZVCMAX25(JL) * ZNITRO_VCMAX * 0.01_JPRB * ZT_RD &
!orig                               &           * MAX(1._JPRB-RSTRESS_VCMAX, ZWATER_LIM(JL))

                 ZIABS(JL)=PSRFD(JL)*RW_TO_MOL*RRG_TO_PAR*REXT_COEFF*ZLIGHT(JLEVEL)
         
                 ! eq. 4 of Yin et al (2009)
                 ZJMAX(JL)=ZVJ2(JL,JLEVEL)
                 ZJJ(JL,JLEVEL) = ( RALPHA_LL * ZIABS(JL) + ZJMAX(JL) - &
                               &   SQRT((RALPHA_LL * ZIABS(JL) + ZJMAX(JL) )**2._JPRB &
                               &   - 4._JPRB * RTHETA * ZJMAX(JL) * RALPHA_LL * ZIABS(JL)) ) &
                               &    / ( 2._JPRB * RTHETA)

              ENDDO !NLAI
               !
       ELSE
               !
               LASSIMILATE(JL) = .FALSE.
               NINA=NINA+1
               INDEX_NON_ASSI(NINA)=JL

!orig              ! give a default value of ci for all pixel that do not ZASSIMIlate
!orig               DO JLEVEL=1,NLAI
!orig                  ZLEAF_CI(JL,JLEVEL) = ZCO2(JL)
!orig               ENDDO
 
               !
       ENDIF !radiation check (assimilate true)

    ELSE
            !
            LASSIMILATE(JL) = .FALSE.
            NINA=NINA+1
            INDEX_NON_ASSI(NINA)=JL
            !
    ENDIF !LDLAND

  ENDDO !JL

 !
 !
 ! 2.4 Loop over LAI discretized levels to estimate ZASSIMIlation and conductance
 !
 !  The calculate(KLON) array is of type logical to indicate whether 
 !  we have to sum over this LAI fixed level (the LAI of
 !  the point for the PFT is lower or equal to the LAI level value). 
 !  The number of such points is incremented in nic and the 
 !  corresponding point is indexed in the index_calc array.
    
 DO JLEVEL = 1, NLAI
        !
        NIC=0
        !
        IF (NIA .GT. 0) then
          DO INIA=1,NIA
            LCALCULATE(INDEX_ASSI(INIA)) = (ZLAITAB(JLEVEL) .LE. PLAI(INDEX_ASSI(INIA)) )
            IF ( LCALCULATE(INDEX_ASSI(INIA))) THEN
               NIC=NIC+1
               INDEX_CALC(NIC)=INDEX_ASSI(INIA)
            ENDIF
          ENDDO
        ENDIF
        !
        ! 2.4.2 ZASSIMIlation for C4 plants (Collatz et al., 1992)
        !
        IF (NIC .GT. 0) THEN
          DO INIC=1,NIC

            ! Analytical resolution of the ZASSIMIlation based Yin et al. (2009)
            ICINIC=INDEX_CALC(INIC)
                   
            IF (ZCTYPE(ICINIC) .EQ. 2_JPIM) THEN
              !C4                

              ! Eq. 28 of Yin et al. (2009)
              ZFCYC = 1._JPRD - ( 4._JPRD *(1._JPRD -REAL(RFPSIR,KIND=JPRD))* &
                  & (1._JPRD + REAL(RFQ,KIND=JPRD))+3._JPRD*REAL(RH_PROTONS,KIND=JPRD)* &
                  & REAL(RFPSEUDO,KIND=JPRD))/(3._JPRD*REAL(RH_PROTONS) - &
                  & 4._JPRD*(1._JPRD-REAL(RFPSIR,KIND=JPRD)))
                                    
              ! See paragraph after eq. (20b) of Yin et al.
              ZRM=REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD)/2._JPRD
                                
              ! We assume that cs_star equals ZCI_STAR (see Comment in legend of Fig. 6 of Yin et al. (2009)
              ! Equation 26 of Yin et al. (2009)
              ZCS_STAR = (REAL(RGBS,KIND=JPRD)*REAL(ZLOW_GAMMA_STAR(ICINIC),KIND=JPRD)* & 
                       &  REAL(ROI,KIND=JPRD)-(1._JPRD+REAL(ZLOW_GAMMA_STAR(ICINIC),KIND=JPRD)* &
                       &  REAL(RALPHA,KIND=JPRD)/0.047_JPRD)*REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD)+ &
                       &  REAL(ZRM,KIND=JPRD))  / ( REAL(RGBS,KIND=JPRD) + REAL(RKP,KIND=JPRD)) 

              ! eq. 11 of Yin et al (2009)
              ZJ2 = REAL(ZJJ(ICINIC,JLEVEL),KIND=JPRD)/ &
                  & (1._JPRD-REAL(RFPSEUDO,KIND=JPRD)/(1._JPRD-ZFCYC))

              ! Equation right after eq. (20d) of Yin et al. (2009)
              ZZ=(2._JPRD+REAL(RFQ,KIND=JPRD)-ZFCYC)/(REAL(RH_PROTONS,KIND=JPRD)*(1._JPRD-ZFCYC))

              ZVPJ2 = REAL(RFPSIR,KIND=JPRD) * ZJ2 * ZZ/ 2._JPRD

              ZA_3=9999._JPRD



              ! See eq. right after eq. 18 of Yin et al. (2009)
              DO LIMIT_PHOTO=1,2
                ! Is Vc limiting the ZASSIMIlation
                IF ( limit_photo .EQ. 1 ) THEN
                  ZA = 1._JPRD + REAL(RKP,KIND=JPRD) / REAL(RGBS,KIND=JPRD)
                  ZB = 0._JPRD
                  ZX1 = REAL(ZVC2(ICINIC,JLEVEL),KIND=JPRD)
                  ZX2 = REAL(ZKMC(ICINIC)/ZKMO(icinic),KIND=JPRD)
                  ZX3 = REAL(ZKMC(ICINIC),KIND=JPRD)
               ! Is J limiting the ZASSIMIlation
                ELSE
                  ZA = 1._JPRD
                  ZB = ZVPJ2
                  ZX1 = (1._JPRD- REAL(RFPSIR,KIND=JPRD)) * ZJ2 * ZZ/ 3._JPRD
                  ZX2 = 7._JPRD * REAL(ZLOW_GAMMA_STAR(ICINIC),KIND=JPRD) / 3._JPRD
                  ZX3 = 0._JPRD
                ENDIF

                ZM = REAL(ZFVPD(ICINIC),KIND=JPRD)-REAL(ZG0VAR(ICINIC),KIND=JPRD)/ &
                   & REAL(ZGBCO2(ICINIC),KIND=JPRD)
                ZD = REAL(ZG0VAR(ICINIC),KIND=JPRD)*(REAL(ZCO2(ICINIC),KIND=JPRD)-ZCS_STAR) + &
                   & REAL(ZFVPD(ICINIC),KIND=JPRD)*REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD)
                ZF = (ZB-ZRM-REAL(ZLOW_GAMMA_STAR(ICINIC),KIND=JPRD)*REAL(ROI,KIND=JPRD)* &
                   & REAL(RGBS,KIND=JPRD))*ZX1*ZD + ZA*REAL(RGBS,KIND=JPRD)* &
                   & ZX1*REAL(ZCO2(ICINIC),KIND=JPRD)*ZD
                ZJ = (ZB-ZRM+REAL(RGBS,KIND=JPRD)*ZX3 + ZX2*REAL(RGBS,KIND=JPRD)* &
                   & REAL(ROI,KIND=JPRD))*ZM + (REAL(RALPHA,KIND=JPRD)* &
                   & ZX2/0.047_JPRD-1._JPRD)*ZD + ZA*REAL(RGBS,KIND=JPRD)* &
                   & (REAL(ZCO2(ICINIC),KIND=JPRD)*ZM - ZD/REAL(ZGBCO2(ICINIC),KIND=JPRD)- &
                   & (REAL(ZCO2(ICINIC),KIND=JPRD)- ZCS_STAR ))
                ZG = (ZB-ZRM-REAL(ZLOW_GAMMA_STAR(ICINIC),KIND=JPRD)*REAL(ROI,KIND=JPRD)* &
                   & REAL(RGBS,KIND=JPRD))*ZX1*ZM-(REAL(RALPHA,KIND=JPRD)* &
                   & REAL(ZLOW_GAMMA_STAR(ICINIC),KIND=JPRD)/0.047_JPRD+1._JPRD)*ZX1*ZD+ &
                   & ZA*REAL(RGBS,KIND=JPRD)*ZX1*(REAL(ZCO2(ICINIC),KIND=JPRD)*ZM- &
                   & ZD/REAL(ZGBCO2(ICINIC),KIND=JPRD)-(REAL(ZCO2(ICINIC),KIND=JPRD)-ZCS_STAR))
                ZH = -((REAL(RALPHA,KIND=JPRD)*REAL(ZLOW_GAMMA_STAR(ICINIC),KIND=JPRD)/ &
                   & 0.047_JPRD+1._JPRD)*ZX1*ZM + (ZA*REAL(RGBS,KIND=JPRD)* &
                   & ZX1*(ZM-1._JPRD))/REAL(ZGBCO2(ICINIC),KIND=JPRD))

                ZI = (ZB-ZRM + REAL(RGBS,KIND=JPRD)*ZX3 + ZX2*REAL(RGBS,KIND=JPRD)* &
                   & REAL(ROI,KIND=JPRD))*ZD + ZA*RGBS*REAL(ZCO2(ICINIC),KIND=JPRD)*ZD
                ZL = (REAL(RALPHA,KIND=JPRD)*ZX2/0.047_JPRD-1._JPRD)*ZM- &
                   & (ZA*REAL(RGBS,KIND=JPRD)*(ZM-1._JPRD))/REAL(ZGBCO2(ICINIC),KIND=JPRD)
   
                ZP = (ZJ-(ZH-ZL*REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD))) / ZL
                ZQ = (ZI+ZJ*REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD)-ZG) / ZL
                ZR = -(ZF-ZI*REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD)) / ZL 

                ! See Yin et al. (2009) and  Baldocchi (1994)
                ZQQ = ( (ZP**2._JPRD) - 3._JPRD * ZQ) / 9._JPRD
                ZUU = ( 2._JPRD * (ZP**3._JPRD) - 9._JPRD *ZP*ZQ + 27._JPRD *ZR) /54._JPRD




                IF ( (ZQQ .GE. 0._JPRD) .AND. (ABS(ZUU/(ABS(ZQQ)**1.5_JPRD) ) .LE. 1._JPRD) ) THEN
                  ZPSI = ACOS(MAX(-1._JPRD,MIN(ZUU/(ZQQ**1.5_JPRD),1._JPRD)))
                  ZA_3_tmp = -2._JPRD*SQRT(ZQQ)*COS((ZPSI+4._JPRD*REAL(RPI,KIND=JPRD))/3._JPRD) &
                          &  - ZP / 3._JPRD

                  IF (( ZA_3_tmp .LT. ZA_3 )) THEN
                    ZA_3 = ZA_3_tmp
                    INFO_LIMITPHOTO(ICINIC,JLEVEL)=2
                  ELSE
                    ! In case, J is not limiting the ZASSIMIlation
                    ! we have to re-initialise a, b, x1, x2 and x3 values
                    ! in agreement with a Vc-limited ZASSIMIlation 
                    ZA = 1._JPRD + REAL(RKP,KIND=JPRD) / REAL(RGBS,KIND=JPRD)
                    ZB = 0._JPRD
                    ZX1 = REAL(ZVC2(ICINIC,JLEVEL),KIND=JPRD)
                    ZX2 = REAL(ZKMC(ICINIC),KIND=JPRD)/REAL(ZKMO(ICINIC),KIND=JPRD)
                    ZX3 = REAL(ZKMC(ICINIC),KIND=JPRD)
                    INFO_LIMITPHOTO(ICINIC,JLEVEL)=1_JPIM
                  ENDIF
                ENDIF
  
                IF ( ( ZA_3 .EQ. 9999._JPRD ) .OR. &
                   & ( ZA_3 .LT. (-REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD)) ) ) THEN
                  ZA_3 = -REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD)
                ENDIF
                ZASSIMI(ICINIC,JLEVEL) = REAL(ZA_3,KIND=JPRB)

                IF ( ABS( ZASSIMI(ICINIC,JLEVEL) + ZRD(ICINIC,JLEVEL) ) .LT. EPSILON(1._JPRM) ) THEN
                   ZGS(ICINIC,JLEVEL) = ZG0VAR(ICINIC)
                   !leaf_ci keeps its initial value (Ca).
                ELSE
                   ! Eq. 24 of Yin et al. (2009) 
                   ZOBS = ( REAL(RALPHA,KIND=JPRD) * REAL(ZASSIMI(ICINIC,JLEVEL),KIND=JPRD) ) / &
                        & ( 0.047_JPRD * REAL(RGBS,KIND=JPRD) ) + REAL(ROI,KIND=JPRD)
                   ! Eq. 23 of Yin et al. (2009)
!ORIG                   ZCC(ICINIC,JLEVEL) = ( ( ZASSIMI(ICINIC,JLEVEL) + ZRD(ICINIC,JLEVEL) ) * &
!ORIG                                     &   ( ZX2 * ZOBS + ZX3 ) + ZLOW_GAMMA_STAR(ICINIC) &
!ORIG                                     & * ZOBS * ZX1 ) &
!ORIG                                     & / MAX(EPSILON(1._JPRM), ZX1 - ( ZASSIMI(ICINIC,JLEVEL) + &
!ORIG                                     &  ZRD(ICINIC,JLEVEL) ))
!ORIG                   ! Eq. 22 of Yin et al. (2009)
!ORIG                   ZLEAF_CI(ICINIC,JLEVEL) = (ZCC(ICINIC,JLEVEL)-(ZB - ZASSIMI(ICINIC,JLEVEL)-ZRM) &
!ORIG                                                 &    / RGBS ) / ZA

                   ZCC = ( ( REAL(ZASSIMI(ICINIC,JLEVEL),KIND=JPRD) + REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD) ) * &
                       &   ( ZX2 * ZOBS + ZX3 ) + REAL(ZLOW_GAMMA_STAR(ICINIC),KIND=JPRD) &
                       & * ZOBS * ZX1 ) &
                       & / MAX(EPSILON(1._JPRM), ZX1 - ( REAL(ZASSIMI(ICINIC,JLEVEL),KIND=JPRD) + &
                       &  REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD) ))
                   ! Eq. 22 of Yin et al. (2009)
                   ZLEAF_CI = (ZCC-(ZB - REAL(ZASSIMI(ICINIC,JLEVEL),KIND=JPRD)-REAL(ZRM,KIND=JPRD)) &
                                                 &    / REAL(RGBS,KIND=JPRD) ) / ZA

                   ! Eq. 25 of Yin et al. (2009)
                   ! It should be Cs instead of Ca but it seems that 
                   ! other equations in Appendix C make use of Ca (ZCO2)
                   ZGS(ICINIC,JLEVEL) = REAL(REAL(ZG0VAR(ICINIC),KIND=JPRD) + &
                                      & ( REAL(ZASSIMI(ICINIC,JLEVEL),KIND=JPRD) + &
                                      & REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD) ) / &
                                      & ( REAL(ZCO2(ICINIC),KIND=JPRD) - ZCS_STAR ) * &
                                      & REAL(ZFVPD(ICINIC),KIND=JPRD),KIND=JPRB)             
                ENDIF
              ENDDO !end limit_photo                  

            ELSE IF (ZCTYPE(ICINIC) .EQ. 1_JPIM) THEN


            !
            ! 2.4.3 ZASSIMIlation for C3 plants (Farqhuar et al., 1980)
            !
              ZA_1=9999._JPRD

              ! See eq. right after eq. 18 of Yin et al. (2009)
              DO LIMIT_PHOTO=1,2
                ! Is Vc limiting the ZASSIMIlation
                IF ( LIMIT_PHOTO .EQ. 1 ) THEN
                  ZX1 = REAL(ZVC2(ICINIC,JLEVEL),KIND=JPRD)
                  ! It should be O not OI (comment from Vuichard)
                  ZX2 = REAL(ZKMC(ICINIC),KIND=JPRD)*(1._JPRD+2._JPRD*REAL(ZGAMMA_STAR(ICINIC),KIND=JPRD)* &
                      & REAL(ZSCO(ICINIC),KIND=JPRD) / REAL(ZKMO(ICINIC),KIND=JPRD) )

                ! Is J limiting the ZASSIMIlation
                ELSE
                  ZX1 = REAL(ZJJ(ICINIC,JLEVEL),KIND=JPRD)/4._JPRD
                  ZX2 = 2._JPRD * REAL(ZGAMMA_STAR(ICINIC),KIND=JPRD)
                ENDIF

                ! See Appendix B of Yin et al. (2009)
                ZA = REAL(ZG0VAR(ICINIC),KIND=JPRD) * (ZX2 + REAL(ZGAMMA_STAR(ICINIC),KIND=JPRD)) + &
                   & (REAL(ZG0VAR(ICINIC),KIND=JPRD) / REAL(ZGM(ICINIC),KIND=JPRD) + &
                   & REAL(ZFVPD(ICINIC),KIND=JPRD)) * (ZX1 - REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD))
                ZB = REAL(ZCO2(ICINIC),KIND=JPRD) * (ZX1 - REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD) ) - &
                   & REAL(ZGAMMA_STAR(ICINIC),KIND=JPRD) * ZX1 - REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD) * ZX2
                ZC = REAL(ZCO2(ICINIC),KIND=JPRD) + ZX2 + (1._JPRD/REAL(ZGM(ICINIC),KIND=JPRD) + &
                   & 1._JPRD/REAL(ZGBCO2(ICINIC),KIND=JPRD) ) *(ZX1 - REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD)) 
                ZD = ZX2 + REAL(ZGAMMA_STAR(ICINIC),KIND=JPRD)+ &
                   & (ZX1 - REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD))/ &
                   &  REAL(ZGM(ICINIC),KIND=JPRD)
                ZM = 1._JPRD/REAL(ZGM(ICINIC),KIND=JPRD) + &
                   & ( REAL(ZG0VAR(ICINIC),KIND=JPRD)/REAL(ZGM(ICINIC),KIND=JPRD) + &
                   & REAL(ZFVPD(ICINIC),KIND=JPRD) ) * &
                   &  ( 1._JPRD/REAL(ZGM(ICINIC),KIND=JPRD) + 1._JPRD/REAL(ZGBCO2(ICINIC),KIND=JPRD) )  
   
                ZP = -( ZD + (ZX1 - REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD) ) / &
                   & REAL(ZGM(ICINIC),KIND=JPRD) + ZA * &
                   & (1._JPRD/REAL(ZGM(ICINIC),KIND=JPRD) + 1._JPRD/REAL(ZGBCO2(ICINIC),KIND=JPRD) ) + &
                   & ( REAL(ZG0VAR(ICINIC),KIND=JPRD)/REAL(ZGM(ICINIC),KIND=JPRD) + &
                   & REAL(ZFVPD(ICINIC),KIND=JPRD) ) * ZC ) / ZM
   
                ZQ = (ZD * ( ZX1 - REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD) ) + ZA*ZC + &
                   & ( REAL(ZG0VAR(ICINIC),KIND=JPRD)/REAL(ZGM(ICINIC),KIND=JPRD) + &
                  &  REAL(ZFVPD(ICINIC),KIND=JPRD) ) * ZB ) / ZM
                ZR = - ZA * ZB / ZM

                ! See Yin et al. (2009) 
                ZQQ = ( (ZP**2._JPRD) - 3._JPRD * ZQ) / 9._JPRD
                ZUU = ( 2._JPRD* (ZP**3._JPRD) - 9._JPRD *ZP*ZQ + 27._JPRD *ZR) /54._JPRD

                IF ( (ZQQ .GE. 0._JPRD) .AND. (ABS(ZUU/(ZQQ**1.5_JPRD) ) .LE. 1._JPRD) ) THEN
                   ZPSI = ACOS(ZUU/(ZQQ**1.5_JPRD))
                   ZA_1_tmp = -2._JPRD * SQRT(ZQQ) * COS( ZPSI / 3._JPRD ) - ZP / 3._JPRD

                  
                   IF (( ZA_1_tmp .LT. ZA_1 )) THEN
                      ZA_1 = ZA_1_tmp
                      INFO_LIMITPHOTO(ICINIC,JLEVEL)=2
                   ELSE
                     ! In case, J is not limiting the ZASSIMIlation
                     ! we have to re-initialise x1 and x2 values
                     ! in agreement with a Vc-limited ZASSIMIlation 
                     ZX1 = REAL(ZVC2(ICINIC,JLEVEL),KIND=JPRD)
                     ! It should be O not OI (comment from Vuichard)
                     ZX2 = REAL(ZKMC(ICINIC),KIND=JPRD) * ( 1._JPRD + 2._JPRD* &
                        &  REAL(ZGAMMA_STAR(ICINIC),KIND=JPRD)*REAL(ZSCO(ICINIC),KIND=JPRD) / REAL(ZKMO(ICINIC),KIND=JPRD) )               
                     INFO_LIMITPHOTO(ICINIC,JLEVEL)=1
                   ENDIF
                ENDIF
              ENDDO !LIMIT_PHOTO 

              IF ( (ZA_1 .EQ. 9999._JPRD) .OR. ( ZA_1 .LT. (-REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD)) ) ) THEN
                ZA_1 = -REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD)
              ENDIF
              ZASSIMI(ICINIC,JLEVEL) = REAL(ZA_1,KIND=JPRB)

              IF (ABS(ZASSIMI(ICINIC,JLEVEL) + ZRD(ICINIC,JLEVEL)) .LT. EPSILON(1._JPRM)) THEN
                ZGS(ICINIC,JLEVEL) = ZG0VAR(ICINIC)
              ELSE
                ! Eq. 18 of Yin et al. (2009)
!ORIG                ZCC(ICINIC,JLEVEL) = ( ZGAMMA_STAR(ICINIC) * ZX1 + ( ZASSIMI(ICINIC,JLEVEL) + &
!ORIG                                 &  ZRD(ICINIC,JLEVEL) ) * ZX2 )  &
!ORIG                                 &  / MAX( EPSILON(1._JPRM), ZX1 - ( ZASSIMI(ICINIC,JLEVEL) +&
!ORIG                                 &  ZRD(ICINIC,JLEVEL) ) )
!ORIG                ZLEAF_CI(ICINIC,JLEVEL) = ZCC(ICINIC,JLEVEL) + ZASSIMI(ICINIC,JLEVEL) / ZGM(ICINIC)

                ZCC = ( REAL(ZGAMMA_STAR(ICINIC),KIND=JPRD) * ZX1 + ( REAL(ZASSIMI(ICINIC,JLEVEL),KIND=JPRD) + &
                    &  REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD) ) * ZX2 )  &
                    &  / MAX( EPSILON(1._JPRM), ZX1 - ( REAL(ZASSIMI(ICINIC,JLEVEL),KIND=JPRD) +&
                    &  REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD) ) )
                ! Eq. 17 of Yin et al. (2009)
                ZLEAF_CI = ZCC + REAL(ZASSIMI(ICINIC,JLEVEL),KIND=JPRD) / REAL(ZGM(ICINIC),KIND=JPRD)
                ! See eq. right after eq. 15 of Yin et al. (2009)
                ZCI_STAR = REAL(ZGAMMA_STAR(ICINIC),KIND=JPRD) - &
                         & REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD) / REAL(ZGM(ICINIC),KIND=JPRD)
                ! 
                ! Eq. 15 of Yin et al. (2009)
!ORIG                ZGS(ICINIC,JLEVEL) = REAL(REAL(ZG0VAR(ICINIC),KIND=JPRD)+ &
!ORIG                                   & REAL((ZASSIMI(ICINIC,JLEVEL),KIND=JPRD) + &
!ORIG                                   &  REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD)) / &
!ORIG                                   & ( REAL(ZLEAF_CI(ICINIC,JLEVEL),KIND=JPRD) &
!ORIG                                   &  - ZCI_STAR ) * REAL(ZFVPD(ICINIC),KIND=JPRD),KIND=JPRB)

                ! Eq. 15 of Yin et al. (2009)
                ZGS(ICINIC,JLEVEL) = REAL(REAL(ZG0VAR(ICINIC),KIND=JPRD)+ &
                                   & (REAL(ZASSIMI(ICINIC,JLEVEL),KIND=JPRD) + &
                                   &  REAL(ZRD(ICINIC,JLEVEL),KIND=JPRD)) / &
                                   & (ZLEAF_CI -ZCI_STAR) * REAL(ZFVPD(ICINIC),KIND=JPRD),KIND=JPRB)

              ENDIF
                  
           ELSE
             !STOP "Wrong ZCTYPE in FARQUHAR"
             write (*,*) "Wrong CTYPE (C3/C4) in FARQUHAR",ZCTYPE(ICINIC)
             ZASSIMI(ICINIC,JLEVEL) = 0._JPRB
             ZGS(ICINIC,JLEVEL) = ZG0VAR(ICINIC)
           ENDIF ! C3/C4
                 
           ENDDO !ICINIC
         ENDIF !nic
         !
         IF (NIC .GT. 0) THEN
         !
           DO INIC=1,NIC
           !
           ! 2.4.4 Estimatation of the stomatal conductance (Ball et al., 1987).
           !
           ICINIC=INDEX_CALC(INIC)
           !
           ! keep stomatal conductance of topmost level
           !
           IF ( JLEVEL .EQ. 1 ) THEN
             ZLEAF_GS_TOP(ICINIC) = ZGS(ICINIC,JLEVEL)
           ENDIF
           !
           !! 2.4.5 Integration at the canopy level
           ! total ZASSIMIlation and conductance
!orig           ZASSIMTOT(ICINIC) = ZASSIMTOT(ICINIC) + &
!orig           &                  ZASSIMI(ICINIC,JLEVEL) * (ZLAITAB(JLEVEL+1)-ZLAITAB(JLEVEL))
!orig           ZRDTOT(ICINIC) = ZRDTOT(ICINIC) + &
!orig                ZRD(ICINIC,JLEVEL) * (ZLAITAB(JLEVEL+1)-ZLAITAB(JLEVEL))
           ZASSIMTOT(ICINIC) = ZASSIMTOT(ICINIC) + &
           &                   (ZASSIMI(ICINIC,JLEVEL)+ZRD(ICINIC,JLEVEL))*&
           &                   (ZLAITAB(JLEVEL+1)-ZLAITAB(JLEVEL))

           PGSTOT(ICINIC) = PGSTOT(ICINIC) + &
           & ZGS(ICINIC,JLEVEL) * (ZLAITAB(JLEVEL+1)-ZLAITAB(JLEVEL))
           !
           ILAI(ICINIC) = JLEVEL
           !
          ENDDO
            !
        ENDIF !NIC
        
        ! Integrate the total dark respiration (assumed to be equal to day respiration)
        ! Dark respiration is required as output to add to heterotrophic respiration
        DO JL=KIDIA,KFDIA
            ZRDTOT(JL) = ZRDTOT(JL) + &
                ZRD(JL,JLEVEL) * (ZLAITAB(JLEVEL+1)-ZLAITAB(JLEVEL))
        ENDDO

      ENDDO  ! loop over LAI steps
 
      !! Calculated intercellular CO2 over NLAI needed for the chemistry module
      !cim(:,JVT)=0._JPRB
      !laisum(:)=0._JPRB
      !DO JLEVEL=1,NLAI
      !   WHERE (ZLAITAB(JLEVEL) .LE. lai(:,JVT) )
      !      cim(:,JVT)= cim(:,JVT)+leaf_ci(:,JVT,JLEVEL)*(ZLAITAB(JLEVEL+1)-ZLAITAB(JLEVEL))
      !      laisum(:)=laisum(:)+ (ZLAITAB(JLEVEL+1)-ZLAITAB(JLEVEL))
      !   ENDWHERE
      !ENDDO
      !WHERE (laisum(:)>0._JPRB)
      !   cim(:,JVT)= cim(:,JVT)/laisum(:)
      !ENDWHERE


      !
      !! 2.5 Calculate resistances
      !
      IF (NIA .GT. 0) THEN
         !
         DO INIA=1,NIA
            !
            IAINIA=INDEX_ASSI(INIA)

            !
            ! cimean is the "mean ci" calculated in such a way that ASSIMILATION
            ! calculated in enerbil is equivalent to assimtot
            !
            !IF ( ABS(gsmean(IAINIA,JVT)-ZG0VAR(IAINIA)*laisum(IAINIA)) .GT. min_sechiba) THEN
            !   cimean(IAINIA,JVT) = (fvpd(IAINIA)*(ZASSIMTOT(IAINIA)+ZRDTOT(IAINIA))) /&
            !     (gsmean(IAINIA,JVT)-ZG0VAR(IAINIA)*laisum(IAINIA)) + ZGAMMA_STAR(IAINIA) 
            !ELSE
            !   cimean(IAINIA,JVT) = ZGAMMA_STAR(IAINIA) 
            !ENDIF
                 
            ! conversion from umol m-2 (PFT) s-1 to kgCO2 m-2 s-1(mesh area)
            PGPP(IAINIA) = ZASSIMTOT(IAINIA)*1E-6_JPRB*RMCO2


!old            ! Dark respiration is required as output to add to heterotrophic respiration
!old             PRD(IAINIA) = PAN(IAINIA)*RRDCF !note this is not equivalent to Ag formula Rd=Am/9 (Am is max assim)
            !
            ! conversion from mol/m^2/s to m/s
            !
            ! As in Pearcy, Schulze and Zimmermann
            ! Measurement of transpiration and leaf conductance
            ! Chapter 8 of Plant Physiological Ecology
            ! Field methods and instrumentation, 1991
            ! Editors:
            !
            !    Robert W. Pearcy,
            !    James R. Ehleringer,
            !    Harold A. Mooney,
            !    Philip W. Rundel
            !
            ! ISBN: 978-0-412-40730-7 (Print) 978-94-010-9013-1 (Online)

            PGSTOT(IAINIA) =  RMOL_TO_M_1 *(PTM(IAINIA)/RTP_00)*&
                 (RPB_STD/ZPB(IAINIA))*PGSTOT(IAINIA)*RRATIO_H2O_TO_CO2

            ZGSTOP(IAINIA) =  RMOL_TO_M_1 * (PTM(IAINIA)/RTP_00)*&
                 (RPB_STD/ZPB(IAINIA))*ZLEAF_GS_TOP(IAINIA)*RRATIO_H2O_TO_CO2*&
                 ZLAITAB(ILAI(IAINIA)+1)

            !! Mean stomatal conductance for water vapour
!            PGSMEAN(IAINIA) = PGSTOT(IAINIA)/PLAI(IAINIA)
            !
            IF (ZGSTOP(IAINIA) .LT. EPSILON(1._JPRM)) THEN
		ZRVEGET(IAINIA) = RUNDEF_SECHIBA
            ELSE
	        ZRVEGET(IAINIA) = 1._JPRB/ZGSTOP(IAINIA)
  	    ENDIF
            !
            !
!check            ! rstruct is the difference between rtot (=1./gstot) and rveget
            !
            ! Correction Nathalie - le 27 Mars 2006 - Interdire a rstruct d'etre negatif
!            IF (PGSTOT(IAINIA) .LT. EPSILON(1._JPRM)) THEN
!                ZRSTRUCT(IAINIA) = RUNDEF_SECHIBA
!            ELSE
! !check               ZRSTRUCT(IAINIA) = 1._JPRB/PGSTOT(IAINIA) - ZRVEGET(IAINIA)
!                ZRSTRUCT(IAINIA) = MAX( 1._JPRB/PGSTOT(IAINIA) - &
!                  &                 ZRVEGET(IAINIA), EPSILON(1._JPRM))
!            ENDIF
            !
            !
            !! wind is a global variable of the diffuco module.
            !speed = MAX(min_wind, wind(IAINIA))
            !
            ! beta for transpiration
            !
            ! Corrections Nathalie - 28 March 2006 - on advices of Fred Hourdin
            !! Introduction of a potentiometer RVEG_PFT to settle the rveg+rstruct sum problem in the coupled mode.
            !! RVEG_PFT=1 in the offline mode. RVEG_PFT is a global variable declared in the diffuco module.
            !vbeta3(IAINIA,JVT) = veget_max(IAINIA,JVT) * &
            !  (1._JPRB - zqsvegrap(IAINIA)) * &
            !  (1._JPRB / (1._JPRB + speed * q_cdrag(IAINIA) * (rveget(IAINIA,JVT) + &
            !   rstruct(IAINIA,JVT))))
            !! Global resistance of the canopy to evaporation
            !ZCRESIST=(1._JPRB / (1._JPRB + speed * q_cdrag(IAINIA) * &
            !     veget(IAINIA,JVT)/veget_max(IAINIA,JVT) * &
            !     (RVEG_PFT*(rveget(IAINIA,JVT) + rstruct(IAINIA,JVT)))))

            !IF ( PF2(IAINIA,JVT) >= min_sechiba ) THEN
            !   vbeta3(IAINIA,JVT) = veget(IAINIA,JVT) * &
            !        (1._JPRB - zqsvegrap(IAINIA)) * ZCRESIST + &
            !        MIN( vbeta23(IAINIA,JVT), veget(IAINIA,JVT) * &
            !        zqsvegrap(IAINIA) * ZCRESIST )
            !ELSE
               ! Because of a minimum conductance G0, vbeta3 cannot be 0._JPRB even if PF2=0
               ! in the above equation.
               ! Here, we force transpiration to be 0._JPRB when the soil cannot deliver it
            !   vbeta3(IAINIA,JVT) = 0._JPRB
            !END IF

            ! vbeta3pot for computation of potential transpiration (needed for irrigation)
            !vbeta3pot(IAINIA,JVT) = MAX(0._JPRB, veget(IAINIA,JVT) * ZCRESIST)
            !
            !
         ENDDO
         !
      ENDIF

     ! Net CO2 assimilation 
      DO JL=KIDIA,KFDIA
            PAN(JL) = (ZASSIMTOT(JL)-ZRDTOT(JL))*1E-6_JPRB*RMCO2
            PRD(JL) = ZRDTOT(JL)*1E-6_JPRB*RMCO2
      ENDDO

      
      !
   !END DO         ! loop over vegetation types
   !
   
   IF (ALLOCATED(ZVC2)   )     DEALLOCATE (ZVC2) 
   IF (ALLOCATED(ZVJ2)   )     DEALLOCATE (ZVJ2) 
   IF (ALLOCATED(ZASSIMI)   )  DEALLOCATE (ZASSIMI) 
   IF (ALLOCATED(ZGS)   )      DEALLOCATE (ZGS) 
   IF (ALLOCATED(ZLAITAB)   )  DEALLOCATE (ZLAITAB) 
   IF (ALLOCATED(ZLIGHT)   )   DEALLOCATE (ZLIGHT) 
   IF (ALLOCATED(ZRD)   )      DEALLOCATE (ZRD) 
   IF (ALLOCATED(ZJJ)   )      DEALLOCATE (ZJJ) 

!ORIG   IF (ALLOCATED(ZLEAF_CI)   ) DEALLOCATE (ZLEAF_CI) 
!ORIG   IF (ALLOCATED(ZCC)   )      DEALLOCATE (ZCC) 
      
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('FARQUHAR_MOD:FARQUHAR',1,ZHOOK_HANDLE)
END SUBROUTINE FARQUHAR



    
END MODULE FARQUHAR_MOD
