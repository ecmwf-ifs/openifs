MODULE YOS_AGS
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
!*    ** *AGS* - 
!     -----------------------------------------------------------------

TYPE :: TAGS
REAL(KIND=JPRB) :: RCO2 
REAL(KIND=JPRB) :: RMAIR
REAL(KIND=JPRB) :: RMH2O
REAL(KIND=JPRB) :: RMCO2
REAL(KIND=JPRB) :: RMC
REAL(KIND=JPRB) :: RDMAX_AGS
REAL(KIND=JPRB) :: RPARCF
REAL(KIND=JPRB) :: RRACCF
REAL(KIND=JPRB) :: RPCCO2
REAL(KIND=JPRB) :: RIAOPT
REAL(KIND=JPRB) :: RDSPOPT
REAL(KIND=JPRB) :: RXGT
REAL(KIND=JPRB) :: RDIFRACF
REAL(KIND=JPRB) :: RXBOMEGA
REAL(KIND=JPRB) :: RRDCF
REAL(KIND=JPRB) :: RAMMIN
REAL(KIND=JPRB) :: RCONDCTMIN
REAL(KIND=JPRB) :: RCONDSTMIN
REAL(KIND=JPRB) :: RANFMINIT
REAL(KIND=JPRB) :: RAIRTOH2O
REAL(KIND=JPRB) :: RCO2TOH2O
REAL(KIND=JPRB) :: RAW
REAL(KIND=JPRB) :: RASW
REAL(KIND=JPRB) :: RBW
REAL(KIND=JPRB) :: RDMAXN
REAL(KIND=JPRB) :: RDMAXX
REAL(KIND=JPRB) :: RRESPFACTOR_NIT
REAL(KIND=JPRB) :: RCNS_NIT
REAL(KIND=JPRB) :: RCA_1x_CO2_NIT
REAL(KIND=JPRB) :: RCA_2x_CO2_NIT
REAL(KIND=JPRB) :: RCC_NIT

REAL(KIND=JPRB) :: RABC(3)
REAL(KIND=JPRB) :: RPOI(3) 

REAL(KIND=JPRB) :: RQ10 !Q10 value in respiration parameterization
REAL(KIND=JPRB) :: RTAULIM !percentage of limitation of efolding time

! Parameters depending on latitudinal band (North/Tropics/South)
REAL(KIND=JPRB) :: RVCH4S(3)

! Parameters depending on photosynthesis mechanism (C3 or c4)

REAL(KIND=JPRB),ALLOCATABLE :: RVTOPT(:)     
REAL(KIND=JPRB),ALLOCATABLE :: RVFZERO(:) 
REAL(KIND=JPRB),ALLOCATABLE :: RVFZEROST(:)     
REAL(KIND=JPRB),ALLOCATABLE :: RVEPSO(:)    
REAL(KIND=JPRB),ALLOCATABLE :: RVGAMM(:)       
REAL(KIND=JPRB),ALLOCATABLE :: RVQDGAMM(:)      
REAL(KIND=JPRB),ALLOCATABLE :: RVQDGMES(:)      
REAL(KIND=JPRB),ALLOCATABLE :: RVT1GMES(:)       
REAL(KIND=JPRB),ALLOCATABLE :: RVT2GMES(:) 
REAL(KIND=JPRB),ALLOCATABLE :: RVAMMAX(:)         
REAL(KIND=JPRB),ALLOCATABLE :: RVQDAMMAX(:)      
REAL(KIND=JPRB),ALLOCATABLE :: RVT1AMMAX(:)
REAL(KIND=JPRB),ALLOCATABLE :: RVT2AMMAX(:)      
REAL(KIND=JPRB),ALLOCATABLE :: RVAH(:)         
REAL(KIND=JPRB),ALLOCATABLE :: RVBH(:)        

! Parameters depending on vegetation type

LOGICAL,ALLOCATABLE :: LVSTRESS(:)      
 
REAL(KIND=JPRB),ALLOCATABLE :: RVBSLAI(:)        
REAL(KIND=JPRB),ALLOCATABLE :: RVLAIMIN(:)      
REAL(KIND=JPRB),ALLOCATABLE :: RVSEFOLD(:)       
REAL(KIND=JPRB),ALLOCATABLE :: RVGMES(:)       
REAL(KIND=JPRB),ALLOCATABLE :: RVGC(:)           
REAL(KIND=JPRB),ALLOCATABLE :: RVDMAX(:)         
REAL(KIND=JPRB),ALLOCATABLE :: RVF2I(:)
REAL(KIND=JPRB),ALLOCATABLE :: RVCE(:)         
REAL(KIND=JPRB),ALLOCATABLE :: RVCF(:)         
REAL(KIND=JPRB),ALLOCATABLE :: RVCNA(:)
REAL(KIND=JPRB),ALLOCATABLE :: RVBSLAI_NITRO(:)
REAL(KIND=JPRB),ALLOCATABLE :: RXBOMEGAM(:)

REAL(KIND=JPRB),ALLOCATABLE :: RVR0VT(:,:) !Reference respiration tabulated by vegetation type
REAL(KIND=JPRB),ALLOCATABLE :: RVCH4QVT(:) !CH4 Q10 temperature dependency by vegetation type

!calculated:
REAL(KIND=JPRB),ALLOCATABLE :: RVANMAX(:)      
END TYPE TAGS

! NAME               TYPE     DESCRIPTION
! ----               ----     -----------
! *RCO2*             REAL     atmospheric CO2 concentration (kgCO2 kgAir-1)
! *RMAIR*            REAL     molecular mass of air (kg mol-1)
! *RMH2O*            REAL     molecular mass of water (kg mol-1)
! *RMCO2*            REAL     molecular mass of CO2 (kg mol-1)
! *RMC*              REAL     molecular mass of C (kg mol-1)
! *RDMAX_AGS*        REAL     maximum specific humidity deficit tolerated by
!                             vegetation (kg kg-1)
!                             for AGS and LAI
! *RPARCF*           REAL     coefficient: PAR fraction of incoming solar 
!                             radiation (0-1)
! *RRACCF*           REAL     factor for aerodynamic resistance for CO2
! *RPCCO2*           REAL     proportion of Carbon in dry plant biomass (0-1)
! *RIAOPT*           REAL     optimum value for absorbed global radiation 
!                             (W m-2)
! *RDSPOPT*          REAL     optimum value for specific humidity deficit 
!                             (kg kg-1)
! *RXGT*             REAL     distribution of leaves
! *RDIFRACF*         REAL     coefficient used in computation of fraction 
!                             of diffusion
! *RXBOMEGA*         REAL     foliage scattering coefficient
! *RXBOMEGAM         REAL     foliage scattering coefficient depending on vegetation type
! *RRDCF*            REAL     dark respiration factor/coefficient
! *RAMMIN*           REAL     minimum criteria for maximum net assimilation 
!                             (kgCO2 kgAir-1 m s-1)
! *RCONDCTMIN*       REAL     minimum canopy conductance (m s-1)
! *RCONDSTMIN*       REAL     minimum stomatal conductance for CO2 (m s-1)
! *RANFMINIT*        REAL     initial maximum leaf assimilation (kgCO2 m2 s-1)
! *RAIRTOH2O*        REAL     ratio of molecular masses of air and H2O
! *RCO2TOH2O*        REAL     ratio of the binary diffusivities of CO2 and H2O 
!                             in air
! *RAW*              REAL     coefficient for stress response of woody species 
! *RASW*             REAL     coefficient for stress response of woody species 
! *RBW*              REAL     coefficient for stress response of woody species 
! *RDMAXN*           REAL     minimum air deficit stress parameters for
!                             herbaceous species (kg kg-1)
! *RDMAXX*           REAL     maximum air deficit stress parameters for
!                             herbaceous species (kg kg-1)
! *RRESPFACTOR_NIT*  REAL     maintenance respiration rate (% per day)
!                             of structural biomass (Fauri�, 1994) (s-1)
! *RCNS_NIT*         REAL     Nitrogen concentration of structural biomass (0-1)
! *RCA_1x_CO2_NIT*   REAL     rate of nitrogen dilution of above-ground biomass
!                             at ambiant (1x) [CO2] (Calvet and Soussana 2001)
! *RCA_2x_CO2_NIT*   REAL     rate of nitrogen dilution of above-ground biomass
!                             at doubled (2x) [CO2] (Calvet and Soussana 2001)
! *RCC_NIT*          REAL     proportion of active biomass for 1t ha-1 of total
!                             above-ground biomass (0-1)

! *RABC       REAL     abscissa needed for integration of net assimilation and 
!                      stomatal conductance over canopy depth (-)
! *RPOI*      REAL     Gaussian weights for integration of net assimilation and 
!                      stomatal conductance over canopy depth (-)

! *RVTOPT*    REAL     optimum temperature for evaluating compensation point (C)
! *RVFZERO*   REAL     ideal value of f (-) (no photorespiration or specific 
!                      humidity deficit)
!                      for AGS and LAI
! *RVFZERO*   REAL     ideal value of f (-) (no photorespiration or specific
!                      humidity deficit)
!                      for AST, LST and NIT 
! *RVEPSO*    REAL     maximum initial quantum use efficiency 
!                      (kgCO2 J-1 PAR m3 kgAir-1)
! *RVGAMM*    REAL     CO2 compensation concentration (kgCO2 kgAir-1)
! *RVQDGAMM*  REAL     Q10 function for CO2 compensation concentration (-)
!  RVQDGMES   REAL     Q10 function for mesophyll conductance (-)
!  RVT1GMES   REAL     minimum reference temperature for computing temperature 
!                      response function for mesophyll conductance (C)
!  RVT2GMES   REAL     maximum reference temperature for computing temperature 
!                      response function for mesophyll conductance (C)
!  RVAMMAX    REAL     leaf photosynthetic capacity (kgCO2 kgAir-1 m s-1)
!  RVQDAMMAX  REAL     Q10 function for leaf photosynthetic capacity (-)
!  RVT1AMMAX  REAL     minimum reference temperature for computing temperature
!                      response function for leaf photosynthetic capacity (C)
!  RVT2AMMAX  REAL     maximum reference temperature for computing temperature
!                      response function for leaf photosynthetic capacity (C)
!  RVAH       REAL     coefficient for herbaceous water stress response 
!  RVBH       REAL     coefficient for herbaceous water stress response 

!  LVSTRESS   LOGICAL  vegetation response type to water stress 
!                      (true:defensive false:offensive) 
!  RVBSLAI    REAL     ratio d(biomass)/d(lai) (kg m-2)
!  RVLAIMIN   REAL     minimum LAI (Leaf Area Index) (m2 m-2)
!  RVSEFOLD   REAL     e-folding time for senescence (s)
!  RVGMES     REAL     mesophyll conductance (m s-1)
!  RVGC       REAL     cuticular conductance (m s-1)
!  RVDMAX     REAL     maximum specific humidity deficit tolerated by 
!                      vegetation (kg kg-1)
!                      for AST, LST and NIT
!  RVF2I      REAL     critical normilized soil water content for stress 
!                      parameterisation (m3 m-3)
!  RVCE       REAL     specific leaf area (SLA) sensitivity to nitrogen 
!                      concentration (m2 kg-1 %-1)
!  RVCF       REAL     lethal minimum value of SLA (m2 kg-1)
!  RVCNA      REAL     nitrogen concentration of active biomass (=leaf biomass)
!                      (%: 0-1)
!  RVBSLAI_NITRO REAL  ratio d(biomass)/d(lai) from nitrogen decline theory 
!                      (kg m-2)

!  RVANMAX    REAL     maximum photosynthesis rate (kgCO2 kgAir-1 m s-1) 
!  RVR0VT     REAL     Reference respiration tabulated by vegetation type and C4 type
!     -----------------------------------------------------------------
END MODULE YOS_AGS
