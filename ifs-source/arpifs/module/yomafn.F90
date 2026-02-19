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

MODULE YOMAFN

USE PARKIND1, ONLY : JPIM
USE PARFPOS, ONLY : JPOSDYN, JPOSSCVA, JPOSPHY, JPOSVX2, JPOSFSU, &
 & JPOSAERO, JPOSGHG, JPOSCHEM, JPOSERA40, JPOSNOGW, JPOSCHEMFLX, & 
 & JPOSAEROUT, JPOSAEROCLIM, JPOSUVP, JPOSSFX, JPOSEZDIAG, JPOSAERODIAG, &
 & JPOSAERAOT, JPOSAERLISI, JPOSAERO_WVL_DIAG, JPOSAERO_WVL_DIAG_TYPES, &
 & JPOSEDRP, JPOSEMIS2D, JPOSEMIS2DAUX, JPOSEMIS3D
USE FULLPOS_MIX, ONLY : FULLPOS_TYPE
USE TYPE_FPDSPHYS, ONLY : FPDSPHY

!====================================================================================
! Identifiers, used for example in FULL-POS.                                        |
! Fields are divided into 5 classes: DYN3D, DYN2D, PHYSOL, CFU, XFU                 |
! Some redundancies may occur.                                                      |
! Identifiers are TFP_[rootfield] for DYN3D and DYN2D                               |
! Identifiers are GFP_[rootfield] for PHYSOL, CFU, XFU                              |
!====================================================================================

IMPLICIT NONE

SAVE

!= Individual 3D dynamic (DYN3D) ====================================================:

TYPE ALL_FULLPOS_TYPES

TYPE(FULLPOS_TYPE) :: Z     ! Geopotential height (g*z).
TYPE(FULLPOS_TYPE) :: T     ! Temperature
TYPE(FULLPOS_TYPE) :: U     ! U-momentum
TYPE(FULLPOS_TYPE) :: V     ! V-momentum
TYPE(FULLPOS_TYPE) :: Q     ! Specific Humidity
TYPE(FULLPOS_TYPE) :: HU    ! Relative Humidity
TYPE(FULLPOS_TYPE) :: VV    ! Pressure coordinate vertical velocity (omega=DPi/Dt)
TYPE(FULLPOS_TYPE) :: VOR   ! Vorticity
TYPE(FULLPOS_TYPE) :: DIV   ! Divergence
TYPE(FULLPOS_TYPE) :: PSI   ! Velocity potential
TYPE(FULLPOS_TYPE) :: KHI   ! Stream function
TYPE(FULLPOS_TYPE) :: L     ! Atmospheric liquid water
TYPE(FULLPOS_TYPE) :: I     ! Atmospheric solid water
TYPE(FULLPOS_TYPE) :: LRAD  ! Total atmospheric liquid water for radiation
TYPE(FULLPOS_TYPE) :: IRAD  ! Total atmospheric solid water for radiation
TYPE(FULLPOS_TYPE) :: SN    ! Atmospheric snow
TYPE(FULLPOS_TYPE) :: RR    ! Atmospheric rain (rain + snow in the case LPROCLD)
TYPE(FULLPOS_TYPE) :: GR    ! Atmospheric graupel
TYPE(FULLPOS_TYPE) :: HL    ! Atmospheric hail
TYPE(FULLPOS_TYPE) :: TKE   ! Turbulent kinetic energy
TYPE(FULLPOS_TYPE) :: EFB1  ! First variable for EFB scheme
TYPE(FULLPOS_TYPE) :: EFB2  ! Second variable for EFB scheme
TYPE(FULLPOS_TYPE) :: EFB3  ! Third variable for EFB scheme
TYPE(FULLPOS_TYPE) :: TH    ! Potential temperature
TYPE(FULLPOS_TYPE) :: THPW  ! Moist (irreversible) pseudo-adiabatic potential temperature (Theta'w).
TYPE(FULLPOS_TYPE) :: TPW   ! Wet bulb (Moist pseudo-adiabatic) temperature (T'w)
TYPE(FULLPOS_TYPE) :: CLF   ! Cloud fraction
TYPE(FULLPOS_TYPE) :: O3MX  ! Ozone mixing ratio
TYPE(FULLPOS_TYPE) :: CPF   ! Convective precipitation flux
TYPE(FULLPOS_TYPE) :: SPF   ! Stratiform precipitation flux
TYPE(FULLPOS_TYPE) :: WND   ! Wind velocity
TYPE(FULLPOS_TYPE) :: ETH   ! Equivalent potential temperature (Theta_e).
TYPE(FULLPOS_TYPE) :: ABS   ! Absolute Vorticity
TYPE(FULLPOS_TYPE) :: STD   ! Stretching Deformation
TYPE(FULLPOS_TYPE) :: PV    ! Potential Vorticity
TYPE(FULLPOS_TYPE) :: SHD   ! Shearing Deformation
TYPE(FULLPOS_TYPE) :: P     ! Hydrostatic pressure
TYPE(FULLPOS_TYPE) :: MG    ! Montgomery geopotential
TYPE(FULLPOS_TYPE) :: PD    ! Pressure departure (NH)
TYPE(FULLPOS_TYPE) :: VD    ! Vertical Divergence (in particuliar in NH model)
TYPE(FULLPOS_TYPE) :: VW    ! True vertical velocity w=Dz/Dt
TYPE(FULLPOS_TYPE) :: SRE   ! Simulated reflectivity in mm/h (microphysics)
TYPE(FULLPOS_TYPE) :: SREDB ! Simulated reflectivity in dBZ (microphysics)
TYPE(FULLPOS_TYPE) :: THV   ! Virtual potential temperature Theta_v (microphysics)
TYPE(FULLPOS_TYPE) :: ETAD  ! Etadot=Deta/Dt (trajectory calculations)
TYPE(FULLPOS_TYPE) :: RKTH  ! Rasch-Kristjansson enthalpy tendency
TYPE(FULLPOS_TYPE) :: RKTQV ! Rasch-Kristjansson water vapour tendency
TYPE(FULLPOS_TYPE) :: RKTQC ! Rasch-Kristjansson condensates tendency
TYPE(FULLPOS_TYPE) :: PHYCTY ! Mass change rate from physics included in CTY

! Prognostic convection & turbulence:
TYPE(FULLPOS_TYPE) :: DAL   ! Downdraught mesh fraction
TYPE(FULLPOS_TYPE) :: DOM   ! Downdraught vertical velocity
TYPE(FULLPOS_TYPE) :: UAL   ! Updraught mesh fraction
TYPE(FULLPOS_TYPE) :: UOM   ! Updraught vertical velocity
TYPE(FULLPOS_TYPE) :: UEN   ! Updraught entrainment
TYPE(FULLPOS_TYPE) :: UNEBH ! Pseudo historic convective cloudiness
TYPE(FULLPOS_TYPE) :: LCONV ! Convective liquid water content
TYPE(FULLPOS_TYPE) :: ICONV ! Convective solid water content
TYPE(FULLPOS_TYPE) :: RCONV ! Convective rain
TYPE(FULLPOS_TYPE) :: SCONV ! Convective snow
!---------------------
TYPE(FULLPOS_TYPE) :: TTE   ! Total turbulent energy
TYPE(FULLPOS_TYPE) :: MXL   ! Prognostic mixing length
TYPE(FULLPOS_TYPE) :: SHTUR ! Shear source term for turbulence           
TYPE(FULLPOS_TYPE) :: FQTUR ! Flux form source term for turbulence -moisture 
TYPE(FULLPOS_TYPE) :: FSTUR ! Flux form source term for turbulence -enthalpy 
TYPE(FULLPOS_TYPE) :: EDR   ! Eddy dissipation rate
!---------------------

TYPE(FULLPOS_TYPE) :: FUA(JPOSVX2)       ! Free upper air fields
TYPE(FULLPOS_TYPE) :: EXT(JPOSSCVA)      ! Extra-GFL variables (former passive scalars)
TYPE(FULLPOS_TYPE) :: EZDIAG(JPOSEZDIAG) ! EZDIAG
TYPE(FULLPOS_TYPE) :: GHG(JPOSGHG)       ! Greenhouse Gases
TYPE(FULLPOS_TYPE) :: CHEM(JPOSCHEM)     ! Chemistry fields 
TYPE(FULLPOS_TYPE) :: AERO(JPOSAERO)     ! Aerosols
TYPE(FULLPOS_TYPE) :: LRCH4              ! CH4 (methane) loss rate
TYPE(FULLPOS_TYPE) :: EMIS3D(JPOSEMIS3D) ! 3D emissions for atmospheric composition

!---------------------

TYPE(FULLPOS_TYPE) :: PTB   ! Pressure of iso-T Kelvin
TYPE(FULLPOS_TYPE) :: HTB   ! Altitude of iso-T Kelvin
TYPE(FULLPOS_TYPE) :: RHO   ! Density of humid air
!----------------------

!= Individual 2D dynamic (DYN2D) ====================================================:

TYPE(FULLPOS_TYPE) :: SP    ! Surface hydrostatic pressure
TYPE(FULLPOS_TYPE) :: SPNH  ! Surface total pressure
TYPE(FULLPOS_TYPE) :: MSL   ! Mean sea level hydrostatic pressure
TYPE(FULLPOS_TYPE) :: MSLNH ! Mean sea level total pressure
TYPE(FULLPOS_TYPE) :: CUF1  ! Filtred ln(Ps), NCUFNR=1
TYPE(FULLPOS_TYPE) :: CUF2  ! Filtred ln(Ps), NCUFNR=2
TYPE(FULLPOS_TYPE) :: CUF3  ! Filtred ln(Ps), NCUFNR=3
TYPE(FULLPOS_TYPE) :: CUF4  ! Filtred ln(Ps), NCUFNR=4
TYPE(FULLPOS_TYPE) :: CUF5  ! Filtred ln(Ps), NCUFNR=5
TYPE(FULLPOS_TYPE) :: QNH   ! QNH
TYPE(FULLPOS_TYPE) :: FIS   ! Interpolated (spectral) model orography
TYPE(FULLPOS_TYPE) :: GM    ! Mapping factor
TYPE(FULLPOS_TYPE) :: FOL   ! Tropopause Folding Indicator
TYPE(FULLPOS_TYPE) :: WWS   ! Surface vertical velocity w_s (NH)
TYPE(FULLPOS_TYPE) :: LNSP  ! Log of surface hydrostatic pressure
TYPE(FULLPOS_TYPE) :: UCLS  ! U cls = U-component of wind at 10 meters (pbl)
TYPE(FULLPOS_TYPE) :: VCLS  ! V cls = V-component of wind at 10 meters (pbl)
TYPE(FULLPOS_TYPE) :: TCLS  ! T cls = Temperature at 2 meters (pbl)
TYPE(FULLPOS_TYPE) :: QCLS  ! Q cls = Specific humidity at 2 meters (pbl)
TYPE(FULLPOS_TYPE) :: RCLS  ! HU cls = Relative humidity at 2 meters (pbl)
TYPE(FULLPOS_TYPE) :: FCLS  ! Module of wind velocity at 2 meters (pbl)
TYPE(FULLPOS_TYPE) :: TX    ! Maximum temperature at 2 meters
TYPE(FULLPOS_TYPE) :: TN    ! Minimum temperature at 2 meters
TYPE(FULLPOS_TYPE) :: CAPE  ! CAPE (Convective available potential energy)
TYPE(FULLPOS_TYPE) :: CIEN  ! CIEN (Convective inhibition energy)
TYPE(FULLPOS_TYPE) :: STRMMU! Storm motion - u component
TYPE(FULLPOS_TYPE) :: STRMMV! Storm motion - v component
TYPE(FULLPOS_TYPE) :: SRH   ! Storm relative helicity
TYPE(FULLPOS_TYPE) :: MOCO  ! MOCON (Moisture convergence)
TYPE(FULLPOS_TYPE) :: TWV   ! Total water vapour content in a vertical column
TYPE(FULLPOS_TYPE) :: UGST  ! U gusts
TYPE(FULLPOS_TYPE) :: VGST  ! V gusts
TYPE(FULLPOS_TYPE) :: FGST  ! Module of gusts
TYPE(FULLPOS_TYPE) :: HCLP  ! Height of the PBL
TYPE(FULLPOS_TYPE) :: VEIN  ! ventilation index in PBL
TYPE(FULLPOS_TYPE) :: UJET  ! ICAO jet U
TYPE(FULLPOS_TYPE) :: VJET  ! ICAO jet V
TYPE(FULLPOS_TYPE) :: PJET  ! ICAO jet pressure
TYPE(FULLPOS_TYPE) :: TCAO  ! ICAO tropopause temperature
TYPE(FULLPOS_TYPE) :: PCAO  ! ICAO tropopause pressure
TYPE(FULLPOS_TYPE) :: HTPW  ! Altitude of iso-Theta'w=0 Celsius, from bottom
TYPE(FULLPOS_TYPE) :: HTPW1 ! Altitude of iso-Theta'w=1 Celsius, from bottom
TYPE(FULLPOS_TYPE) :: HTPW2 ! altitude of iso-T'w=1.5 Celsius, from bottom
TYPE(FULLPOS_TYPE) :: HUX   ! Recomputed maximum relative moisture at 2 meters
TYPE(FULLPOS_TYPE) :: HUN   ! Recomputed minimum relative moisture at 2 meters
TYPE(FULLPOS_TYPE) :: SREX  ! Recomputed maximum simulated reflectivities in mm/h
TYPE(FULLPOS_TYPE) :: SREDBX! Recomputed maximum simulated reflectivities in dBZ
TYPE(FULLPOS_TYPE) :: TOPR  ! Pressure of top of reflectivities
TYPE(FULLPOS_TYPE) :: IET   ! Isobaric equivalent temperature (for ALARO)
TYPE(FULLPOS_TYPE) :: SMC   ! Forecast surface moisture convergence
TYPE(FULLPOS_TYPE) :: ASMC  ! Analysed surface moisture convergence
TYPE(FULLPOS_TYPE) :: VSMC  ! varpack-purpose surface moisture convergence

TYPE(FULLPOS_TYPE) :: FSU(JPOSFSU) ! Free surface field

TYPE(FULLPOS_TYPE) :: MSAT7C1
TYPE(FULLPOS_TYPE) :: MSAT7C2                            ! Meteosat 7 MVIRI channel 1 and 2
TYPE(FULLPOS_TYPE) :: MSAT8C1
TYPE(FULLPOS_TYPE) :: MSAT8C2
TYPE(FULLPOS_TYPE) :: MSAT8C3
TYPE(FULLPOS_TYPE) :: MSAT8C4     ! Meteosat 8 SEVIRI channels 1 to 4
TYPE(FULLPOS_TYPE) :: MSAT8C5
TYPE(FULLPOS_TYPE) :: MSAT8C6
TYPE(FULLPOS_TYPE) :: MSAT8C7
TYPE(FULLPOS_TYPE) :: MSAT8C8     ! Meteosat 8 SEVIRI channels 5 to 8
TYPE(FULLPOS_TYPE) :: MSAT9C1
TYPE(FULLPOS_TYPE) :: MSAT9C2
TYPE(FULLPOS_TYPE) :: MSAT9C3
TYPE(FULLPOS_TYPE) :: MSAT9C4     ! Meteosat 9 SEVIRI channels 1 to 4
TYPE(FULLPOS_TYPE) :: MSAT9C5
TYPE(FULLPOS_TYPE) :: MSAT9C6
TYPE(FULLPOS_TYPE) :: MSAT9C7
TYPE(FULLPOS_TYPE) :: MSAT9C8     ! Meteosat 9 SEVIRI channels 5 to 8
TYPE(FULLPOS_TYPE) :: GOES11C1
TYPE(FULLPOS_TYPE) :: GOES11C2
TYPE(FULLPOS_TYPE) :: GOES11C3
TYPE(FULLPOS_TYPE) :: GOES11C4 ! GOES 11 Imager channels 1 to 4
TYPE(FULLPOS_TYPE) :: GOES12C1
TYPE(FULLPOS_TYPE) :: GOES12C2
TYPE(FULLPOS_TYPE) :: GOES12C3
TYPE(FULLPOS_TYPE) :: GOES12C4 ! GOES 12 Imager channels 1 to 4
TYPE(FULLPOS_TYPE) :: MTSAT1C1
TYPE(FULLPOS_TYPE) :: MTSAT1C2
TYPE(FULLPOS_TYPE) :: MTSAT1C3
TYPE(FULLPOS_TYPE) :: MTSAT1C4 ! MTSAT-1R Imager channels 1 to 4
                                                                   
TYPE(FULLPOS_TYPE) :: LCL  ! Lifting Condensation level
TYPE(FULLPOS_TYPE) :: FCL  ! Free convection level
TYPE(FULLPOS_TYPE) :: EL   ! Equilibrium level
TYPE(FULLPOS_TYPE) :: TCVS ! Temperature of convection

TYPE(FULLPOS_TYPE) :: NOGW(JPOSNOGW)     ! Diagnostic fields for NORO GWD scheme
TYPE(FULLPOS_TYPE) :: UVP(JPOSUVP)       ! UV-processor output fields
TYPE(FULLPOS_TYPE) :: MSK                ! Mask extra domain
TYPE(FULLPOS_TYPE) :: ERA40(JPOSERA40)   ! ERA40 diagnostic fields
TYPE(FULLPOS_TYPE) :: EDRP(JPOSEDRP)     ! Diagnostic fields of EDR Parameters for CAT and MWT

TYPE(FULLPOS_TYPE) :: AEROUT(JPOSAEROUT)     ! Aerosol output fields
TYPE(FULLPOS_TYPE) :: AEROCLIM(JPOSAEROCLIM) ! Aerosol climatology
TYPE(FULLPOS_TYPE) :: AERAOT(JPOSAERAOT)     ! Aerosol optical thicknesses
TYPE(FULLPOS_TYPE) :: AERLISI(JPOSAERLISI)   ! Aerosol lidar simulator

END TYPE ALL_FULLPOS_TYPES

! Individual Surface physical fields  =================================================:

TYPE ALL_FPDSPHY_TYPES

TYPE(FPDSPHY) :: ALUVP                ! MODIS-derived parallel albedo UV-vis  (ECMWF)
TYPE(FPDSPHY) :: ALUVD                ! MODIS-derived diffuse albedo UV-vis   (ECMWF)
TYPE(FPDSPHY) :: ALNIP                ! MODIS-derived parallel albedo Near-IR (ECMWF)
TYPE(FPDSPHY) :: ALNID                ! MODIS-derived diffuse albedo Near-IR  (ECMWF)
TYPE(FPDSPHY) :: ALUVI                ! MODIS-derived isotropic UV-vis albedo component (ECMWF)
TYPE(FPDSPHY) :: ALUVV                ! MODIS-derived volumetric UV-vis albedo component (ECMWF)
TYPE(FPDSPHY) :: ALUVG                ! MODIS-derived geometric UV-vis albedo component (ECMWF)
TYPE(FPDSPHY) :: ALNII                ! MODIS-derived isotropic near-IR albedo component (ECMWF)
TYPE(FPDSPHY) :: ALNIV                ! MODIS-derived volumetric near-IR albedo component (ECMWF)
TYPE(FPDSPHY) :: ALNIG                ! MODIS-derived geometric near-IR albedo component (ECMWF)
TYPE(FPDSPHY) :: SDFOR                ! Standard deviation of filtered subgrid orography (units m)
TYPE(FPDSPHY) :: CGPP                 ! GPP flux adjustment coefficient
TYPE(FPDSPHY) :: CREC                 ! REC flux adjustment coefficient
TYPE(FPDSPHY) :: LSM                  ! Land/sea mask
TYPE(FPDSPHY) :: GFIS                 ! OUTPUT Grid-point orography (times g)
TYPE(FPDSPHY) :: SFIS                 ! INPUT Grid-point orography (times g)
TYPE(FPDSPHY) :: ST                   ! Surface temperature
TYPE(FPDSPHY) :: DST                  ! Deep soil temperature
TYPE(FPDSPHY) :: RDST                 ! INTERPOLATED surface temperature
TYPE(FPDSPHY) :: SSW                  ! Surface soil wetness
TYPE(FPDSPHY) :: DSW                  ! Deep soil wetness
TYPE(FPDSPHY) :: FSSW                 ! Frozen superficial soil wetness
TYPE(FPDSPHY) :: FDSW                 ! Frozen deep soil wetness
TYPE(FPDSPHY) :: RDSW                 ! Climatology relaxation relative soil wetness
TYPE(FPDSPHY) :: CSSW                 ! Climatology relative surface soil wetness
TYPE(FPDSPHY) :: CDSW                 ! Climatology relative deep soil wetness
TYPE(FPDSPHY) :: SD                   ! Snow depth
TYPE(FPDSPHY) :: SDSL                 ! Snow depth total
TYPE(FPDSPHY) :: SR                   ! Surface roughness (times g)
TYPE(FPDSPHY) :: BSR                  ! Roughness length of the bare surface (times g)
TYPE(FPDSPHY) :: AL                   ! Albedo
TYPE(FPDSPHY) :: EMIS                 ! Emissivity
TYPE(FPDSPHY) :: SDOG                 ! Standard deviation of orography (times g)
TYPE(FPDSPHY) :: VEG                  ! Percentage of vegetation
TYPE(FPDSPHY) :: SOTY                 ! Soil type 
TYPE(FPDSPHY) :: LAN                  ! Percentage of land
TYPE(FPDSPHY) :: CLK                  ! Lake cover
TYPE(FPDSPHY) :: DL                   ! Lake depth
TYPE(FPDSPHY) :: LMLT                 ! Lake mix layer temperature
TYPE(FPDSPHY) :: LMLD                 ! Lake mix layer depth
TYPE(FPDSPHY) :: LBLT                 ! Lake bottom layer temperature
TYPE(FPDSPHY) :: LTLT                 ! Lake total layer temperature
TYPE(FPDSPHY) :: LSHF                 ! Lake shape factor
TYPE(FPDSPHY) :: LICT                 ! Lake ice temperature
TYPE(FPDSPHY) :: LICD                 ! Lake ice thickness
TYPE(FPDSPHY) :: ACOT                 ! Anisotropy coefficient of topography
TYPE(FPDSPHY) :: DPAT                 ! Direction of the principal axis of the topography
TYPE(FPDSPHY) :: IVEG                 ! Index of vegetation
TYPE(FPDSPHY) :: RSMIN                ! Stomatal minimum resistance
TYPE(FPDSPHY) :: ARG                  ! Percentage of clay within soil
TYPE(FPDSPHY) :: SAB                  ! Percentage of sand within soil
TYPE(FPDSPHY) :: D2                   ! Soil Depth
TYPE(FPDSPHY) :: LAI                  ! Leaf area index
TYPE(FPDSPHY) :: HV                   ! Resistance to evapotranspiration
TYPE(FPDSPHY) :: Z0H                  ! Roughness length for heat (times g)
TYPE(FPDSPHY) :: IC                   ! Interception content
TYPE(FPDSPHY) :: ALSN                 ! Surface snow albedo
TYPE(FPDSPHY) :: SNDE                 ! Surface snow density
TYPE(FPDSPHY) :: BAAL                 ! Albedo of the bare surface (old)
TYPE(FPDSPHY) :: ALBHIS               ! Albedo total
TYPE(FPDSPHY) :: ALS                  ! Albedo of bare ground
TYPE(FPDSPHY) :: ALV                  ! Albedo of vegetation
TYPE(FPDSPHY) :: PSRHU                ! Surface relative moisture
TYPE(FPDSPHY) :: PADOU                ! U-component of vector anisotropy
TYPE(FPDSPHY) :: PADOV                ! V-component of vector anisotropy
TYPE(FPDSPHY) :: PCAAG                ! Analysed RMS of geopotential (out of CANARI)
TYPE(FPDSPHY) :: PCAPG                ! Forecasted RMS of geopotential (out of CANARI)
TYPE(FPDSPHY) :: IDZ0                 ! Interpolated dynamic surface roughness length (times g) 
TYPE(FPDSPHY) :: ITZ0                 ! Interpolated thermal surface roughness length (times g)
TYPE(FPDSPHY) :: PVGMX                ! Maximum proportion of vegetation
TYPE(FPDSPHY) :: Z0V                  ! Roughness length for vegetation (times g)
TYPE(FPDSPHY) :: PURB                 ! Proportion of urbanisation
TYPE(FPDSPHY) :: D2MX                 ! Maximum soil depth
TYPE(FPDSPHY) :: STL1                 ! Soil first level temperature (ECMWF)
TYPE(FPDSPHY) :: STL2                 ! Soil second level temperature (ECMWF)
TYPE(FPDSPHY) :: STL3                 ! Soil third level temperature (ECMWF)
TYPE(FPDSPHY) :: STL4                 ! Soil fourth level temperature (ECMWF)
TYPE(FPDSPHY) :: SWL1                 ! Soil first level wetness (ECMWF)
TYPE(FPDSPHY) :: SWL2                 ! Soil second level wetness (ECMWF)
TYPE(FPDSPHY) :: SWL3                 ! Soil third level wetness (ECMWF)
TYPE(FPDSPHY) :: SWL4                 ! Soil fourth level wetness (ECMWF)
TYPE(FPDSPHY) :: TSN                  ! Temperature of snow layer (ECMWF)
TYPE(FPDSPHY) :: ISOR                 ! Anisotropy of surface orography (ECMWF)
TYPE(FPDSPHY) :: ANOR                 ! Angle of surface orography (ECMWF)
TYPE(FPDSPHY) :: SLOR                 ! Slope of surface orography (ECMWF)
TYPE(FPDSPHY) :: LSRH                 ! Logarithm of surface roughness (ECMWF)
TYPE(FPDSPHY) :: SRC                  ! Skin wetness (ECMWF)
TYPE(FPDSPHY) :: SKT                  ! Skin temperature (ECMWF)
TYPE(FPDSPHY) :: LSP                  ! Large scale precipitation (ECMWF)
TYPE(FPDSPHY) :: CP                   ! Convective precipitation (ECMWF)
TYPE(FPDSPHY) :: TP                   ! Total precipitation (ECMWF)
TYPE(FPDSPHY) :: SF                   ! Snowfall (ECMWF)
TYPE(FPDSPHY) :: FZRA                 ! Freezing rain (ECMWF)
TYPE(FPDSPHY) :: BLD                  ! Boundary layer dissipation (ECMWF)
TYPE(FPDSPHY) :: SSHF                 ! Surface sensible heat flux (ECMWF)
TYPE(FPDSPHY) :: SLHF                 ! Surface latent heat flux (ECMWF)
TYPE(FPDSPHY) :: NEE                  ! Surface net ecosystem exchange of CO2 (ECMWF)
TYPE(FPDSPHY) :: GPP                  ! Surface gross primary production of CO2 (ECMWF)
TYPE(FPDSPHY) :: REC                  ! Surface ecosystem respiration of CO2 (ECMWF)
TYPE(FPDSPHY) :: MSLD                 ! Mean sea level pressure (ECMWF)
TYPE(FPDSPHY) :: SP                   ! Surface pressure (ECMWF)
TYPE(FPDSPHY) :: TCC                  ! Total cloud cover (ECMWF)
TYPE(FPDSPHY) :: EC10U                ! U-wind at 10 m (ECMWF)
TYPE(FPDSPHY) :: EC10V                ! V-wind at 10 m (ECMWF)
TYPE(FPDSPHY) :: EC10SI               ! Wind speed at 10 m (ECMWF)
TYPE(FPDSPHY) :: EC2T                 ! Temperature at 2 m (ECMWF)
TYPE(FPDSPHY) :: EC2D                 ! Dewpoint at 2 m (ECMWF)
TYPE(FPDSPHY) :: EC2SH                ! Specific humidity at 2 m (ECMWF)
TYPE(FPDSPHY) :: SSR                  ! Surface solar radiation (ECMWF)
TYPE(FPDSPHY) :: STR                  ! Surface thermal radiation (ECMWF)
TYPE(FPDSPHY) :: TSR                  ! Top solar radiation (ECMWF)
TYPE(FPDSPHY) :: TTR                  ! Top thermal radiation (ECMWF)
TYPE(FPDSPHY) :: EWSS                 ! U-wind stress (ECMWF)
TYPE(FPDSPHY) :: NSSS                 ! V-wind stress (ECMWF)
TYPE(FPDSPHY) :: E                    ! Water evaporation (ECMWF)
TYPE(FPDSPHY) :: PEV                  ! Potential water evaporation (ECMWF)
TYPE(FPDSPHY) :: CCC                  ! Convective cloud cover (ECMWF)
TYPE(FPDSPHY) :: LCC                  ! Low cloud cover (ECMWF)
TYPE(FPDSPHY) :: MCC                  ! Medium cloud cover (ECMWF)
TYPE(FPDSPHY) :: HCC                  ! High cloud cover (ECMWF)
TYPE(FPDSPHY) :: LGWS                 ! Zonal gravity wave stress (ECMWF)
TYPE(FPDSPHY) :: MGWS                 ! Meridian gravity wave stress (ECMWF)
TYPE(FPDSPHY) :: GWD                  ! Gravity wave dissipation (ECMWF)
TYPE(FPDSPHY) :: MX2T                 ! Maximum temperature at 2 m since last p.p. (ECMWF)
TYPE(FPDSPHY) :: MN2T                 ! Minimum temperature at 2 m since last p.p. (ECMWF)
TYPE(FPDSPHY) :: MX2T3                ! Maximum temperature at 2 m since last 3 hours (ECMWF)
TYPE(FPDSPHY) :: MN2T3                ! Minimum temperature at 2 m since last 3 hours (ECMWF)
TYPE(FPDSPHY) :: MX2T6                ! Maximum temperature at 2 m since last 6 hours (ECMWF)
TYPE(FPDSPHY) :: MN2T6                ! Minimum temperature at 2 m since last 6 hours (ECMWF)
TYPE(FPDSPHY) :: MXTPR                ! Maximum total precip since last p.p. (ECMWF)    
TYPE(FPDSPHY) :: MNTPR                ! Minimum total precip since last p.p. (ECMWF)    
TYPE(FPDSPHY) :: MXTPR3               ! Maximum total precip since last 3 hours (ECMWF) 
TYPE(FPDSPHY) :: MNTPR3               ! Minimum total precip since last 3 hours (ECMWF) 
TYPE(FPDSPHY) :: MXTPR6               ! Maximum total precip since last 6 hours (ECMWF) 
TYPE(FPDSPHY) :: MNTPR6               ! Minimum total precip since last 6 hours (ECMWF) 
TYPE(FPDSPHY) :: LSRR                 ! Large-scale rainfall rate (ECMWF)
TYPE(FPDSPHY) :: CRR                  ! Convective rainfall rate (ECMWF)
TYPE(FPDSPHY) :: LSSFR                ! Large-scale snowfall rate (ECMWF)
TYPE(FPDSPHY) :: CSFR                 ! Convective snowfall rate (ECMWF)
TYPE(FPDSPHY) :: PTYPE                ! Precipitation type (ECMWF)
TYPE(FPDSPHY) :: ILSPF                ! Instantaneous large-scale precipitation fraction (ECMWF)
TYPE(FPDSPHY) :: RO                   ! Runoff (total) (ECMWF)
TYPE(FPDSPHY) :: SRO                  ! Surface Runoff (ECMWF)
TYPE(FPDSPHY) :: SSRO                 ! Sub-Surface Runoff (ECMWF)
TYPE(FPDSPHY) :: ALB                  ! Albedo (ECMWF)
TYPE(FPDSPHY) :: IEWSS                ! Instantaneous surface zonal component of stress
TYPE(FPDSPHY) :: INSSS                ! Instantaneous surface meridian component of stress
TYPE(FPDSPHY) :: ISSHF                ! Instantaneous surface heat flux
TYPE(FPDSPHY) :: IE                   ! Instantaneous surface moisture flux (ECMWF)
TYPE(FPDSPHY) :: INEE                 ! Instantaneous net ecosystem exchange of CO2 (ECMWF)
TYPE(FPDSPHY) :: IGPP                 ! Instantaneous gross primary production of CO2 (ECMWF)
TYPE(FPDSPHY) :: IREC                 ! Instantaneous ecosystem respiration of CO2 (ECMWF)
TYPE(FPDSPHY) :: ICH4                 ! Instantaneous wetland CH4 (ECMWF)
TYPE(FPDSPHY) :: CSF                  ! Convective snow fall (ECMWF)
TYPE(FPDSPHY) :: LSF                  ! Large scale snowfall (ECMWF)
TYPE(FPDSPHY) :: Z0F                  ! Surface roughness (ECMWF)
TYPE(FPDSPHY) :: LZ0H                 ! Logarithm of z0 times heat flux (ECMWF)
TYPE(FPDSPHY) :: TCW                  ! Total water content in a vertical column (ECMWF)
TYPE(FPDSPHY) :: TCWV                 ! Total water vapor content in a vertical column (ECMWF)
TYPE(FPDSPHY) :: SSRD                 ! Downward surface solar radiation (ECMWF)
TYPE(FPDSPHY) :: STRD                 ! Downward surface thermic radiation (ECMWF)
TYPE(FPDSPHY) :: SSRDC                ! Clear-sky downward surface solar radiation (ECMWF)
TYPE(FPDSPHY) :: STRDC                ! Clear-sky downward surface thermic radiation (ECMWF)
TYPE(FPDSPHY) :: TCO3                 ! Total ozone content in a vertical column (ECMWF)
TYPE(FPDSPHY) :: TCGHG(JPOSGHG)       ! Total column greenhouse Gases
TYPE(FPDSPHY) :: TCCHEM(JPOSCHEM)     ! Total column chemical fields
TYPE(FPDSPHY) :: CHEMFLXO(JPOSCHEMFLX)! Surface total flux chemical fields (emissions + deposition) 
TYPE(FPDSPHY) :: CHEMWDFLX(JPOSCHEMFLX)! Dry deposition flux chemical fields  
TYPE(FPDSPHY) :: CHEMDDFLX(JPOSCHEMFLX)! Wet deposition flux chemical fields  

TYPE(FPDSPHY) :: SUND                 ! Sunshine duration (ECMWF)
TYPE(FPDSPHY) :: CHAR                 ! Charnock parameter (ECMWF)
TYPE(FPDSPHY) :: BLH                  ! Height of boundary layer (ECMWF)
TYPE(FPDSPHY) :: BV                   ! Budget values
TYPE(FPDSPHY) :: TCLW                 ! Total column liquid water
TYPE(FPDSPHY) :: TCIW                 ! Total column ice water
TYPE(FPDSPHY) :: TCRW                 ! Total column rain water
TYPE(FPDSPHY) :: TCSW                 ! Total column snow water
TYPE(FPDSPHY) :: TCSLW                ! Total column supercooled liquid water
TYPE(FPDSPHY) :: VX2(JPOSVX2)         ! Extra surface fields
TYPE(FPDSPHY) :: FSU(JPOSFSU)         ! Free surface fields
TYPE(FPDSPHY) :: CVL                  ! Low vegetation cover
TYPE(FPDSPHY) :: CO2TYP               ! CO2 photosynthesis type (c3/c4) for low vegetation cover
TYPE(FPDSPHY) :: CVH                  ! High vegetation cover
TYPE(FPDSPHY) :: FWET                 ! Wetland fraction
TYPE(FPDSPHY) :: CUR                  ! Urban cover
TYPE(FPDSPHY) :: TVL                  ! Low vegetation type
TYPE(FPDSPHY) :: TVH                  ! High vegetation type
TYPE(FPDSPHY) :: LAIL                 ! Low vegetation LAI
TYPE(FPDSPHY) :: LAIH                 ! High vegetation LAI
TYPE(FPDSPHY) :: CI                   ! Sea ice cover
TYPE(FPDSPHY) :: ASN                  ! Snow albedo
TYPE(FPDSPHY) :: RSN                  ! Snow density
TYPE(FPDSPHY) :: WSN                  ! Snow Liquid water content
TYPE(FPDSPHY) :: SST                  ! Sea surface temperature
TYPE(FPDSPHY) :: ISTL1                ! Ice surface temperature - layer 1
TYPE(FPDSPHY) :: ISTL2                ! Ice surface temperature - layer 2
TYPE(FPDSPHY) :: ISTL3                ! Ice surface temperature - layer 3
TYPE(FPDSPHY) :: ISTL4                ! Ice surface temperature - layer 4
TYPE(FPDSPHY) :: TSRC                 ! Top solar radiation clear sky
TYPE(FPDSPHY) :: TTRC                 ! Top thermal radiation clear sky
TYPE(FPDSPHY) :: SSRC                 ! Surface solar radiation clear sky
TYPE(FPDSPHY) :: STRC                 ! Surface thermal radiation clear sky
TYPE(FPDSPHY) :: ES                   ! Evaporation of snow
TYPE(FPDSPHY) :: SMLT                 ! Snow melt
TYPE(FPDSPHY) :: LSPF                 ! Large scale precipitation fraction
TYPE(FPDSPHY) :: EC10FG               ! Wind gust at 10 m (max since previous p.p.)
TYPE(FPDSPHY) :: EC10FG3              ! Wind gust at 10 m (max since last 3 hours)
TYPE(FPDSPHY) :: EC10FG6              ! Wind gust at 10 m (max since last 6 hours)
TYPE(FPDSPHY) :: I10FG                ! Wind gust at 10 m ("instantaneous")
TYPE(FPDSPHY) :: ASEA                 ! Marine aerosols
TYPE(FPDSPHY) :: ALAN                 ! Continental aerosols
TYPE(FPDSPHY) :: ASOO                 ! Carbone aerosols
TYPE(FPDSPHY) :: ADES                 ! Desert aerosols
TYPE(FPDSPHY) :: ASUL                 ! Sulfate aerosols
TYPE(FPDSPHY) :: AVOL                 ! Volcano aerosols
TYPE(FPDSPHY) :: O3A                  ! First ozone profile (A one)
TYPE(FPDSPHY) :: O3B                  ! Second ozone profile (B one)
TYPE(FPDSPHY) :: O3C                  ! Third ozone profile (C one)
TYPE(FPDSPHY) :: SPAR                 ! Surface photo active radiation
TYPE(FPDSPHY) :: SUVB                 ! Surface UV-B radiation
TYPE(FPDSPHY) :: CAPE                 ! Convective available potential energy (CAPE)
TYPE(FPDSPHY) :: MUCAPE               ! Maximum unstable CAPE
TYPE(FPDSPHY) :: MLCAPE50             ! CAPE from 50 hPa mixed layer
TYPE(FPDSPHY) :: MLCAPE100            ! CAPE from 100 hPa mixed layer
TYPE(FPDSPHY) :: PDEPL                ! Pressure at parcel departure level for MUCAPE
TYPE(FPDSPHY) :: CAPES                ! CAPE-shear
TYPE(FPDSPHY) :: MXCAP6               ! Maximum Conv.avail.potential energy in last 6h
TYPE(FPDSPHY) :: MXCAPS6              ! Maximum CAPE-Shear in last 6h
TYPE(FPDSPHY) :: TROPOTP              ! Pressure of thermal tropopause
TYPE(FPDSPHY) :: SPARC                ! Surface clear-sky parallel radiation
TYPE(FPDSPHY) :: STINC                ! TOA incident solar radiation
TYPE(FPDSPHY) :: SFDIR                ! Surface total sky direct shortwave radiation
TYPE(FPDSPHY) :: SCDIR                ! Surface clear sky direct shortwave radiation
TYPE(FPDSPHY) :: SDSRP                ! Surface total sky direct beam shortwave radiation
TYPE(FPDSPHY) :: VIMD                 ! Vertically integrated mass divergence
TYPE(FPDSPHY) :: CBASE                ! Cloud base level
TYPE(FPDSPHY) :: EC0DEGL              ! Zero degree Celsius level
TYPE(FPDSPHY) :: ECM10DEGL            ! Minus 10 degree Celsius level
TYPE(FPDSPHY) :: VISIH                ! Horizontal visibility
TYPE(FPDSPHY) :: CIN                  ! Convective Inhibition
TYPE(FPDSPHY) :: MLCIN50              ! CIN from 50 hPa mixed layer
TYPE(FPDSPHY) :: MLCIN100             ! CIN from 100 hPa mixed layer
TYPE(FPDSPHY) :: KINDEX               ! Convective K-Index 
TYPE(FPDSPHY) :: TTINDEX              ! Convective TT-Index
TYPE(FPDSPHY) :: CBASEA               ! Cloud base level aviation
TYPE(FPDSPHY) :: CTOPC                ! Cloud top  level convection
TYPE(FPDSPHY) :: ZTWETB0              ! Wet bulb 0-degree Celsius level
TYPE(FPDSPHY) :: ZTWETB1              ! Wet bulb 1-degree Celsius level
TYPE(FPDSPHY) :: TPR                  ! Total precipitation rate (ECMWF)
TYPE(FPDSPHY) :: AERODIAG(JPOSAERO,JPOSAERODIAG)                          ! Per-aerosol-type diagnostic fields
TYPE(FPDSPHY) :: AERO_WVL_DIAG(JPOSAERO_WVL_DIAG,JPOSAERO_WVL_DIAG_TYPES) ! Per-aerosol-type diagnostic fields
TYPE(FPDSPHY) :: EMIS2D(JPOSEMIS2D)   ! 2D emission fluxes for composition
TYPE(FPDSPHY) :: EMIS2DAUX(JPOSEMIS2DAUX) ! 2D emission auxiliary fields for composition

TYPE(FPDSPHY) :: SO2DD                ! SO2 dry deposition velocity       (ECMWF)
TYPE(FPDSPHY) :: VIWVE                ! Vertical integral of eastward water vapour flux (ECMWF)
TYPE(FPDSPHY) :: VIWVN                ! Vertical integral of northward water vapour flux (ECMWF)
TYPE(FPDSPHY) :: URBF                 ! Urban fraction                    (ECMWF)
TYPE(FPDSPHY) :: ODSS                 ! Optical depth sea salt aerosols   (ECMWF)
TYPE(FPDSPHY) :: ODDU                 ! Optical depth dust aerosols       (ECMWF)
TYPE(FPDSPHY) :: ODOM                 ! Optical depth organic matter aerosols (ECMWF)
TYPE(FPDSPHY) :: ODBC                 ! Optical depth black C aerosols    (ECMWF)
TYPE(FPDSPHY) :: ODSU                 ! Optical depth sulphate aerosols   (ECMWF)
TYPE(FPDSPHY) :: ODNI                 ! Optical depth nitrate aerosols    (ECMWF)
TYPE(FPDSPHY) :: ODAM                 ! Optical depth ammonium aerosols   (ECMWF)
TYPE(FPDSPHY) :: ODSOA                ! Optical depth secondary org. aer. (ECMWF)
TYPE(FPDSPHY) :: ODVFA                ! Optical depth volc.fly ash aeros. (ECMWF)
TYPE(FPDSPHY) :: ODVSU                ! Optical depth volc.sulphate aeros.(ECMWF)
TYPE(FPDSPHY) :: ODTOACC              ! Optical depth total aerosol acc.  (ECMWF)
TYPE(FPDSPHY) :: AEPM1                ! Partic.Matter le 1 um             (ECMWF)
TYPE(FPDSPHY) :: AEPM25               ! Partic.Matter le 2.5 um           (ECMWF)
TYPE(FPDSPHY) :: AEPM10               ! Partic.Matter le 10 um            (ECMWF)
TYPE(FPDSPHY) :: UVBED                ! UV Biologically Effective Dose    (ECMWF)
TYPE(FPDSPHY) :: UVBEDCS              ! UV Biologically Effective Dose Clear Sky (ECMWF)
TYPE(FPDSPHY) :: DMSO                 ! Oceanic DMS (dimethylsulfide)     (ECMWF)
TYPE(FPDSPHY) :: AERDEP               ! Dust emission potential           (ECMWF)
TYPE(FPDSPHY) :: AERLTS               ! Lifting threshold speed           (ECMWF)
TYPE(FPDSPHY) :: AERSCC               ! Soil clay content                 (ECMWF)
TYPE(FPDSPHY) :: DSF                  ! Dust Source Function              (ECMWF)
TYPE(FPDSPHY) :: DSZ                  ! Dust Size distribution modulation (ECMWF)
TYPE(FPDSPHY) :: UCUR                 ! Ocean current U component         (ECMWF)
TYPE(FPDSPHY) :: VCUR                 ! Ocean current V component         (ECMWF)
TYPE(FPDSPHY) :: EC100U               ! U-wind at 100 m                   (ECMWF)
TYPE(FPDSPHY) :: EC100V               ! V-wind at 100 m                   (ECMWF)
TYPE(FPDSPHY) :: EC100SI              ! Wind speed at 100 m               (ECMWF)
TYPE(FPDSPHY) :: EC200U               ! U-wind at 200 m                   (ECMWF)
TYPE(FPDSPHY) :: EC200V               ! V-wind at 200 m                   (ECMWF)
TYPE(FPDSPHY) :: EC200SI              ! Wind speed at 200 m               (ECMWF)
TYPE(FPDSPHY) :: ZUST                 ! Friction velocity                 (ECMWF)
TYPE(FPDSPHY) :: EC10NU               ! Neutral U-wind at 10 m            (ECMWF)
TYPE(FPDSPHY) :: EC10NV               ! Neutral V-wind at 10 m            (ECMWF)
TYPE(FPDSPHY) :: DNDZN                ! Minimum vertical refractivity gradient inside trapping layer (ECMWF)
TYPE(FPDSPHY) :: DNDZA                ! Mean vertical refractivity gradient inside trapping layer    (ECMWF)
TYPE(FPDSPHY) :: DCTB                 ! Duct base height                   (ECMWF)
TYPE(FPDSPHY) :: TPLB                 ! Trapping layer base height         (ECMWF)
TYPE(FPDSPHY) :: TPLT                 ! Trapping layer top height          (ECMWF)
TYPE(FPDSPHY) :: LITOTI               ! Total lightning flash density (instantaneous) (ECMWF)
TYPE(FPDSPHY) :: LITOTA1              ! Total lightning flash density (averaged over past 1h) (ECMWF)
TYPE(FPDSPHY) :: LITOTA3              ! Total lightning flash density (averaged over past 3h) (ECMWF)
TYPE(FPDSPHY) :: LITOTA6              ! Total lightning flash density (averaged over past 6h) (ECMWF)
TYPE(FPDSPHY) :: LICGI                ! Cloud-to-ground flash density (instantaneous) (ECMWF)
TYPE(FPDSPHY) :: LICGA1               ! Cloud-to-ground lightning flash density (averaged over past 1h) (ECMWF)
TYPE(FPDSPHY) :: LICGA3               ! Cloud-to-ground lightning flash density (averaged over past 3h) (ECMWF)
TYPE(FPDSPHY) :: LICGA6               ! Cloud-to-ground lightning flash density (averaged over past 6h) (ECMWF)
TYPE(FPDSPHY) :: PTYPEMODE1           ! Most frequent precip type in last 1h (mode)
TYPE(FPDSPHY) :: PTYPEMODE3           ! Most frequent precip type in last 3h (mode)
TYPE(FPDSPHY) :: PTYPEMODE6           ! Most frequent precip type in last 6h (mode)
TYPE(FPDSPHY) :: PTYPESEVR1           ! Most severe precip type in last 1h
TYPE(FPDSPHY) :: PTYPESEVR3           ! Most severe precip type in last 3h
TYPE(FPDSPHY) :: PTYPESEVR6           ! Most severe precip type in last 6h
TYPE(FPDSPHY) :: CLBT                 ! Cloudy brightness temperature     (ECMWF)
TYPE(FPDSPHY) :: CSBT                 ! Clear-sky brightness temperature  (ECMWF)
TYPE(FPDSPHY) :: FCA1                 ! Fraction of calcite over dust 1st bin   (ECMWF)
TYPE(FPDSPHY) :: FCA2                 ! Fraction of calcite over dust 2nd bin   (ECMWF)
TYPE(FPDSPHY) :: ICTH                 ! Ice thickess          (ECMWF)
TYPE(FPDSPHY) :: EC20D                ! Depth of 20d isotherm (ECMWF)
TYPE(FPDSPHY) :: SSH                  ! Sea surface height    (ECMWF)
TYPE(FPDSPHY) :: MLD                  ! Mixed layer depth     (ECMWF)
TYPE(FPDSPHY) :: SSS                  ! Sea surface salinity  (ECMWF)
TYPE(FPDSPHY) :: TEM3                 ! Pottem avg 300m       (ECMWF)
TYPE(FPDSPHY) :: SAL3                 ! Sal avg 300m          (ECMWF)

! Individual surface cumulated fluxes (CFU) ========================================:

TYPE(FPDSPHY) :: CLSP    ! Large Scale liquid precipitation
TYPE(FPDSPHY) :: CCP     ! Convective liquid precipitation
TYPE(FPDSPHY) :: CLSS    ! Large Scale Snow fall
TYPE(FPDSPHY) :: CCSF    ! Convective Snow Fall
TYPE(FPDSPHY) :: CLSG    ! Large Scale Graupel fall
TYPE(FPDSPHY) :: CCSG    ! Convective Graupel Fall
TYPE(FPDSPHY) :: CLSH    ! Large Scale Hail fall
TYPE(FPDSPHY) :: CCSH    ! Convective Hail Fall
TYPE(FPDSPHY) :: CUSS    ! U-stress
TYPE(FPDSPHY) :: CVSS    ! V-stress
TYPE(FPDSPHY) :: CSSH    ! Surface Sensible Heat Flux
TYPE(FPDSPHY) :: CSLH    ! Surface Latent Heat Flux
TYPE(FPDSPHY) :: CTSP    ! Tendency of Surface pressure
TYPE(FPDSPHY) :: CTCC    ! Total Cloud cover
TYPE(FPDSPHY) :: CBLD    ! Boundary Layer Dissipation
TYPE(FPDSPHY) :: CSMR    ! Surface downward moon radiation
TYPE(FPDSPHY) :: CSSR    ! Surface solar radiation
TYPE(FPDSPHY) :: CSTR    ! Surface Thermal radiation
TYPE(FPDSPHY) :: CTSR    ! Top Solar radiation
TYPE(FPDSPHY) :: CTTR    ! Top Thermal radiation
TYPE(FPDSPHY) :: CCCC    ! Convective Cloud Cover
TYPE(FPDSPHY) :: CHCC    ! High Cloud Cover
TYPE(FPDSPHY) :: CMCC    ! Medium Cloud ;over
TYPE(FPDSPHY) :: CLCC    ! Low Cloud Cover
TYPE(FPDSPHY) :: CUGW    ! U-Gravity-Wave Stress
TYPE(FPDSPHY) :: CVGW    ! V-Gravity-Wave Stress
TYPE(FPDSPHY) :: CUTO    ! U-Total Stress = U-stress + U-Gravity-Wave Stress
TYPE(FPDSPHY) :: CVTO    ! V-Total Stress = V-stress + V-Gravity-Wave Stress
TYPE(FPDSPHY) :: CE      ! Water Evaporation
TYPE(FPDSPHY) :: CS      ! Snow Sublimation
TYPE(FPDSPHY) :: CT      ! Water & Snow Sublimation
TYPE(FPDSPHY) :: CLHE    ! Latent Heat Evaporation
TYPE(FPDSPHY) :: CLHS    ! Latent Heat Sublimation
TYPE(FPDSPHY) :: CLHT    ! Total Latent Heat 
TYPE(FPDSPHY) :: CWS     ! Soil Moisture
TYPE(FPDSPHY) :: CSNS    ! Snow mass
TYPE(FPDSPHY) :: CQTO    ! Total precipitable water
TYPE(FPDSPHY) :: CTO3    ! Total Ozone
TYPE(FPDSPHY) :: CTME    ! Top mesospheric enthalpy
TYPE(FPDSPHY) :: CICE    ! Solid specific moisture
TYPE(FPDSPHY) :: CLI     ! Liquid specific moisture
TYPE(FPDSPHY) :: CCVU    ! Contribution of Convection to U
TYPE(FPDSPHY) :: CCVV    ! Contribution of Convection to V
TYPE(FPDSPHY) :: CCVQ    ! Contribution of Convection to Q
TYPE(FPDSPHY) :: CCVS    ! Contribution of Convection to Cp.T
TYPE(FPDSPHY) :: CTUQ    ! Contribution of Turbulence to Q
TYPE(FPDSPHY) :: CTUS    ! Contribution of Turbulence to Cp.T
TYPE(FPDSPHY) :: CSOC    ! Surface clear sky shortwave radiative flux
TYPE(FPDSPHY) :: CTSOC   ! Top clear sky shortwave radiative flux
TYPE(FPDSPHY) :: CTHC    ! Surface clear sky longwave radiative flux
TYPE(FPDSPHY) :: CTTHC   ! Top clear sky longwave radiative flux
TYPE(FPDSPHY) :: CSOP    ! Surface parallel solar flux
TYPE(FPDSPHY) :: CTOP    ! Top parallel solar flux
TYPE(FPDSPHY) :: CSOD    ! Surface down solar flux
TYPE(FPDSPHY) :: CTHD    ! Surface down thermic flux
TYPE(FPDSPHY) :: CDNI    ! Surface direct normal irradiance
TYPE(FPDSPHY) :: CGNI    ! Surface global normal irradiance
TYPE(FPDSPHY) :: CFON    ! Melt snow
TYPE(FPDSPHY) :: CCHS    ! Heat flux in soil
TYPE(FPDSPHY) :: CEAS    ! Water flux in soil
TYPE(FPDSPHY) :: CSRU    ! Surface soil runoff
TYPE(FPDSPHY) :: CDRU    ! Deep soil runoff
TYPE(FPDSPHY) :: CIRU    ! Interception soil layer runoff
TYPE(FPDSPHY) :: CETP    ! Evapotranspiration flux
TYPE(FPDSPHY) :: CTP     ! Transpiration flux
TYPE(FPDSPHY) :: CDUTP   ! Duration of total precipitations
TYPE(FPDSPHY) :: CFLASH  ! Lightning density (cumulative)

! Individual surface instantaneous fluxes/diagnostics (XFU) ========================:

TYPE(FPDSPHY) :: XTCC    ! Total Cloud cover
TYPE(FPDSPHY) :: X10U    ! 10 meters U-component of wind
TYPE(FPDSPHY) :: X10V    ! 10 meters V-component of wind
TYPE(FPDSPHY) :: X10FF   ! 10 meters wind
TYPE(FPDSPHY) :: X2T     ! 2 meters Temerature
TYPE(FPDSPHY) :: X2TPW   !2 meters Wet Bulb Temperature
TYPE(FPDSPHY) :: XMRT    ! Mean Radiant Temperature
TYPE(FPDSPHY) :: X2SH    ! 2 meters Specific Humidity
TYPE(FPDSPHY) :: X2RH    ! 2 meters Relative Humidity
TYPE(FPDSPHY) :: XCCC    ! Convective Cloud Cover
TYPE(FPDSPHY) :: XHCC    ! High Cloud Cover
TYPE(FPDSPHY) :: XMCC    ! Medium Cloud Cover
TYPE(FPDSPHY) :: XLCC    ! Low Cloud Cover
TYPE(FPDSPHY) :: XX2T    ! Maximum Temperature at 2 meters since previous p.p
TYPE(FPDSPHY) :: XN2T    ! Minimum Temperature at 2 meters since previous p.p.
TYPE(FPDSPHY) :: XCVU    ! Contribution of Convection to U
TYPE(FPDSPHY) :: XCVV    ! Contribution of Convection to V
TYPE(FPDSPHY) :: XCVQ    ! Contribution of Convection to Q
TYPE(FPDSPHY) :: XCVS    ! Contribution of Convection to Cp.T
TYPE(FPDSPHY) :: XTUU    ! Contribution of Turbulence to U
TYPE(FPDSPHY) :: XTUV    ! Contribution of Turbulence to V
TYPE(FPDSPHY) :: XTUQ    ! Contribution of Turbulence to Q
TYPE(FPDSPHY) :: XTUS    ! Contribution of Turbulence to Cp.T
TYPE(FPDSPHY) :: XGDU    ! Contribution of Gravity Wave drag to U
TYPE(FPDSPHY) :: XGDV    ! Contribution of Gravity Wave drag to V
TYPE(FPDSPHY) :: XLSP    ! Large Scale Precipitation
TYPE(FPDSPHY) :: XCP     ! Convective precipitation
TYPE(FPDSPHY) :: XLSS    ! Large Scale Snow fall
TYPE(FPDSPHY) :: XCSF    ! Convective Snow Fall
TYPE(FPDSPHY) :: XLSG    ! Large Scale Graupel fall
TYPE(FPDSPHY) :: XCSG    ! Convective Graupel Fall
TYPE(FPDSPHY) :: XLSH    ! Large Scale Hail fall
TYPE(FPDSPHY) :: XCSH    ! Convective Hail Fall
TYPE(FPDSPHY) :: XSSR    ! Surface solar radiation
TYPE(FPDSPHY) :: XSTR    ! Surface Thermal radiation
TYPE(FPDSPHY) :: XTSR    ! Top Solar radiation
TYPE(FPDSPHY) :: XTTR    ! Top Thermal radiation
TYPE(FPDSPHY) :: XCAPE   ! CAPE out of the model
TYPE(FPDSPHY) :: XCTOP   ! Pressure of top of deep convection
TYPE(FPDSPHY) :: XMOCO   ! MOCON out of the model
TYPE(FPDSPHY) :: XCLPH   ! Height (in meters) of the PBL out of the model
TYPE(FPDSPHY) :: XVEIN   ! Ventilation index in PBL
TYPE(FPDSPHY) :: XGUST   ! Gusts out of the model
TYPE(FPDSPHY) :: XUGST   ! U-momentum of gusts out of the model
TYPE(FPDSPHY) :: XVGST   ! V-momentum of gusts out of the model
TYPE(FPDSPHY) :: XX2HU   ! Maximum relative moisture at 2 meters
TYPE(FPDSPHY) :: XN2HU   ! Minimum relative moisture at 2 meters
TYPE(FPDSPHY) :: INCTX   ! Increment to maxi temperature (ECMWF)
TYPE(FPDSPHY) :: INCTN   ! Increment to mini temperature (ECMWF)
TYPE(FPDSPHY) :: INCHX   ! Increment to maxi relative moisture (ECMWF)
TYPE(FPDSPHY) :: INCHN   ! Increment to mini relative moisture (ECMWF)
TYPE(FPDSPHY) :: XTHW    ! Theta'w' surface flux (AROME)
TYPE(FPDSPHY) :: XXDIAGH ! Hail diagnostic (AROME)
TYPE(FPDSPHY) :: ACCGREL ! Dummy hail field (AROME)
TYPE(FPDSPHY) :: X10NU   ! 10 meters U-neutral wind
TYPE(FPDSPHY) :: X10NV   ! 10 meters V-neutral wind
TYPE(FPDSPHY) :: XUGST2  ! U-momentum of gusts2 out of the model
TYPE(FPDSPHY) :: XVGST2  ! V-momentum of gusts2 out of the model
TYPE(FPDSPHY) :: VISICLD ! Visibility due to fog
TYPE(FPDSPHY) :: VISIHYD ! Visibility due to precipitations
TYPE(FPDSPHY) :: MXCLWC  ! Maximum of CLWC
TYPE(FPDSPHY) :: VISICLD2! Visibility due to fog
TYPE(FPDSPHY) :: VISIHYD2! Visibility due to precipitations
TYPE(FPDSPHY) :: MXCLWC2 ! Maximum of CLWC
TYPE(FPDSPHY) :: XPTYPE  ! Frequent Precipitation Type
TYPE(FPDSPHY) :: XPTYPESEV ! Severe Precipitation Type
TYPE(FPDSPHY) :: XPTYPE2 ! Frequent Precipitation Type
TYPE(FPDSPHY) :: XPTYPESEV2! Severe Precipitation Type

! PREP field descriptors
TYPE(FPDSPHY) :: SFXPRE (JPOSSFX) 

END TYPE ALL_FPDSPHY_TYPES


TYPE TAFN

! All fields for dynamics ==========================================================

TYPE(FULLPOS_TYPE) :: TFP_DYNDS(JPOSDYN)

TYPE(ALL_FULLPOS_TYPES) :: TFP

! All fields for physics ===========================================================

TYPE(FPDSPHY) :: GFP_PHYDS(JPOSPHY)

TYPE(ALL_FPDSPHY_TYPES) :: GFP

! Miscellaneous controls of physical fields  =======================================

INTEGER(KIND=JPIM) :: MSTARTGFP_PHY ! Start index for physical fields in GFP_PHYDS
INTEGER(KIND=JPIM) :: MSTARTGFP_CFU ! Start index for cumulated fluxes in GFP_PHYDS
INTEGER(KIND=JPIM) :: MSTARTGFP_XFU ! Start index for instantaneous fluxes in GFP_PHYDS
INTEGER(KIND=JPIM) :: MENDGFP_PHY   ! End index for physical fields in GFP_PHYDS
INTEGER(KIND=JPIM) :: MENDGFP_CFU   ! End index for cumulated fluxes in GFP_PHYDS
INTEGER(KIND=JPIM) :: MENDGFP_XFU   ! End index for instantaneous fluxes in GFP_PHYDS

! PREP maski numbers, start from 101 to 104
INTEGER (KIND=JPIM) :: NSFXMSK_CNT
INTEGER (KIND=JPIM) :: NSFXMSK_BAS
! PREP masks
CHARACTER (LEN=8), ALLOCATABLE :: CLSFXMASK (:)
! NSFXPRE_CNT : number of Surfex PREP fields to remember
INTEGER(KIND=JPIM) :: NSFXPRE_CNT

!     ------------------------------------------------------------------

END TYPE TAFN

END MODULE YOMAFN


