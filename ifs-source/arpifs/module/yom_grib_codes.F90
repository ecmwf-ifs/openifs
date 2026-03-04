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

MODULE YOM_GRIB_CODES

USE PARKIND1  ,ONLY : JPIM
USE PAR_GFL , ONLY : JPGHG

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    GRIB CODING DESCRIPTORS


! NGRB...  - Fixed GRIB codes according to ECMWF local Code Table 2
! The characters after NGRB should be the same as the MARS short name !!!!!!
! NGRBd.. where d.. are digits means not assigned (yet) in MARS
! Some of the descriptions below might not agree with MARS manual (they should)

! NGRBSTRF  -  1 Stream function
! NGRBVP    -  2 Velocity potential
! NGRBPT    -  3 Potential Temperature
! NGRBSRO   -  8 Surface Runoff
! NGRBSSRO  -  9 Sub-Surface Runoff 
! NGRBALUVP - 15 MODIS albedo UV-vis parallel radiation
! NGRBALUVD - 16 MODIS albedo UV-vis diffuse radiation
! NGRBALNIP - 17 MODIS albedo Near-IR parallel radiation
! NGRBALNID - 18 MODIS albedo Near-IR diffuse radiation

! NGRBALUVI - 210186 MODIS albedo UV-vis parallel radiation (isom.)
! NGRBALNII - 210189 MODIS albedo Near-IR parallel radiation (isom.)
! NGRBALUVV - 210187 MODIS albedo UV-vis parallel radiation (volu.)
! NGRBALNIV - 210190 MODIS albedo Near-IR parallel radiation (volu.)
! NGRBALUVG - 210188 MODIS albedo UV-vis parallel radiation (geom.)
! NGRBALNIG - 210191 MODIS albedo Near-IR parallel radiation (geom.)

! NGRBPARCS - 20 surface clear-sky PARadiation
! NGRBCL   -  26 Lake cover
! NGRBCVL  -  27 Low vegetation cover
! NGRBCVH  -  28 High vegetation cover
! NGRBFWET -  200026 Wetland Fraction (experimental)
! NGRBCUR  -  200199 Urban Cover (experimental)
! NGRBTVL  -  29 Low vegetation type
! NGRBTVH  -  30 High vegetation type
! NGRBCI   -  31 Sea ice cover
! NGRBASN  -  32 Snow albedo
! NGRBRSN  -  33 Snow density
! NGRBSSTK  - 34 Sea surface temperature
! NGRBISTL1 - 35 Ice surface temperature layer 1
! NGRBISTL2 - 36 Ice surface temperature layer 2
! NGRBISTL3 - 37 Ice surface temperature layer 3
! NGRBISTL4 - 38 Ice surface temperature layer 4
! NGRBSWVL1 - 39 Volumetric soil water content layer 1
! NGRBSWVL2 - 40 Volumetric soil water content layer 2
! NGRBSWVL3 - 41 Volumetric soil water content layer 3      
! NGRBSWVL4 - 42 Volumetric soil water content layer 4 
! NGRBSLT   - 43 Soil type
! NGRBES   -  44 Evaporation of snow
! NGRBSMLT -  45 Snow melt
! NGRBDSRP -  47 Direct solar radiation 
!                Incident on a plane perpendicular to the Sun's direction
! NGRB10FG -  49 gust at 10 m level
! NGRBLSPF -  50 large scale precipitation fraction 

! NGRBMONT  - 53 Montgomery Geopotential
! NGRBPRES  - 54 Pressure on Theta and PV surfaces

! NGRBUVB  -  57 surface UV-B radiation
! NGRBPAR  -  58 surface PARadiation
! NGRBCAPE -  59 convect.avail.potential energy
! NGRBMLCAPE50  - 228231 Convective Inhibition for near surface 50 hPa mixed layer parcel
! NGRBMLCAPE100 - 228233 Convective Inhibition for near surface 100 hPa mixed layer parcel
! NGRBMUCAPE - 228235 Most Unstable CAPE (using Tv)
! NGRBMUDEPL  - 228237 Departure level (Pa) of Most Unstable CAPE
! NGRBCAPES-  228044 wind shear*sqrt(cape)
! NGRBMXCAP6- 228035 maximum CAPE in the last 6 hours
! NGRBMXCAPS6-228036 maximum CAPES in the last 6 hours
 
! NGRBPV   -  60 Potential Vorticity
! NGRBLAIL -  66 Leaf Area Index Low vegitation
! NGRBLAIH -  67 Leaf Area Index High vegitation

! NGRBSDFOR - 74 standard deviation of a filtered orography

! NGRBCRWC -  75 Precipitating rain water content
! NGRBCSWC -  76 Precipitating snow water content

! NGRBTCLW -  78 Total column liquid water
! NGRBTCIW -  79 Total column ice water
! NGRBTCRW -  228089 Total column rain water
! NGRBTCSW -  228090 Total column snow water
! NGRBTCSLW - 228088 Total column supercooled liquid water
! NGRBSPD  -  80 !! 80 and 81 extra grib code introduced to 
! NGRBSVD  -  81 !! introduce extra fields for NH Not MARS codes!!!!!!
! NGRBALPHA-  80 Alpha control variable

! NGRB082 to NGRB117 reserved for extra fields. Do not use for permanent post-processed fields

! Codes for the lake model
! NGRBDL    - 228007 Lake depth
! NGRBLMLT  - 228008 Lake mix-layer temperature
! NGRBLMLD  - 228009 Lake mix-layer depth
! NGRBLBLT  - 228010 Lake bottom layer temperature
! NGRBLTLT  - 228011 Lake total layer temperature
! NGRBLSHF  - 228012 Lake shape factor
! NGRBLICT  - 228013 Lake ice temperature
! NGRBLICD  - 228014 Lake ice depth

! Codes for 2D and 3D extra fields
! NGRBMINXTRA  to NGRBMAXXTRA

! NGRBMX2T3 - 228026 Maximum temperature at 2 m since last 3 hours 
! NGRBMN2T3 - 228027 Minimum temperature at 2 m since last 3 hours 
! NGRB10FG3 - 228028 Wind gust at 10 metres since last 3 hours
! NGRBI10FG - 228029 Wind gust at 10 metres ("instantaneous")
! NGRBMX2T6 - 121 Maximum temperature at 2 m since last 6 hours 
! NGRBMN2T6 - 122 Minimum temperature at 2 m since last 6 hours 
! NGRB10FG6 - 123 Wind gust at 10 metres since last 6 hours 
! NGRBEMIS  - 124 Surface Longwave emissivity # has replaced NGRB212

! NGRBETADOT - 077 Etadotdpdeta

! NGRBAT   - 127 Atmospheric tide
! NGRBBV   - 128 Budget values
! NGRBZ    - 129 Geopotential (at the surface orography)
! NGRBT    - 130 Temperature
! NGRBU    - 131 U-velocity
! NGRBV    - 132 V-velocity
! NGRBUCUR - 131 U-velocity, using (ocean) table 151
! NGRBVCUR - 132 V-velocity, using (ocean) table 151
! NGRBSSS  - 151130 sea surface salinity, (ocean table 151)

! NGRBQ    - 133 Specific humidity
! NGRBSP   - 134 Surface pressure
! NGRBW    - 135 Vertical velocity
! NGRBTCW  - 136 Total column water
! NGRBTCWV - 137 Total column water vapour
! NGRBVO   - 138 Vorticity (relative)
! NGRBSTL1 - 139 Surface temperature level 1
! NGRBSDSL - 141 Snow depth (total)
! NGRBSD   - 228141 Snow depth (multi-layer)
! NGRBWSN  - 228038 Snow liquid water (multi-layer)
! NGRBLSP  - 142 Large scale precipitation
! NGRBCP   - 143 Convective precipitation
! NGRBSF   - 144 Snow fall
! NGRBFZRA - 228216 Freezing rain accumulation
! NGRBBLD  - 145 Boundary layer dissipation
! NGRBSSHF - 146 Surface sensible heat flux
! NGRBSLHF - 147 Surface latent heat flux
! NGRBCHNK - 148 Charnock parameter 
! NGRBSNR  - 149 Surface net radiation
! NGRBTNR  - 150 Top net radiation
! NGRBMSL  - 151 Mean sea level pressure
! NGRBLNSP - 152 Log surface pressure
! NGRBSWHR - 153 Not used?
! NGRBLWHR - 154 Not used?
! NGRBD    - 155 Divergence
! NGRBGH   - 156 Height (geopotential)
! NGRBR    - 157 Relative humidity
! NGRBTSP  - 158 Tendency of surface pressure
! NGRBBLH  - 159 Boundary layer height
! NGRBSDOR - 160 Standard deviation of orography
! NGRBISOR - 161 Anisotropy of subgrid scale orography
! NGRBANOR - 162 Angle of subgrid scale orography
! NGRBSLOR - 163 Slope of subgrid scale orography
! NGRBTCC  - 164 Total cloud cover
! NGRB10U  - 165 10 metre u wind
! NGRB10V  - 166 10 metre v wind
! NGRBZUST - 228003 Friction velocity
! NGRBFDIR - 228021 Surface total sky direct SW
! NGRBCDIR - 228022 Surface clear sky direct SW
! NGRBCBASE- 228023 Cloud base level
! NGRB0DEGL- 228024 Zero deg. level
! NGRBM10DEGL- 228020 -10 deg. level
! NGRBVISIH- 3020 Visibility ! Changed from 228025 to be WMO compliant
! NGRBCIN  - 228001 Convective Inhibition
! NGRBMLCIN50  - 228232 Convective Inhibition for near surface 50 hPa mixed layer parcel
! NGRBMLCIN100 - 228234 Convective Inhibition for near surface 100 hPa mixed layer parcel
! NGRBKINDEX - 260121 Convective K-Index
! NGRBTTINDEX- 260123 Convective TT-Index
! NGRBCBASEA-260109 Ceiling = Cloud base level aviation
! NGRBCTOPC- 228046 Cloud top level convective
! NGRBZTWETB0-228047 Zero Deg Wet Bulb temperature
! NGRBZTWETB1-228048 One  Deg Wet Bulb temperature
! NGRBTROPOTP-228045 Pressure at thermal tropopause
! NGRB10NU - 228131 10 metre u neutral wind
! NGRB10NV - 228132 10 metre v neutral wind
! NGRB2T   - 167 2 metre temperature
! NGRB2D   - 168 2 metre dewpoint temperature
! NGRB2SH  - 174096 2 metre specific humidity
! NGRBSSRD - 169 Surface solar radiation downwards
! NGRBSTL2 - 170 Soil temperature level 2
! NGRBLSM  - 172 Land/sea mask
! NGRBSR   - 173 Surface roughness
! NGRBAL   - 174 Albedo
! NGRBSTRD - 175 Surface thermal radiation downwards
! NGRBSSR  - 176 Surface solar radiation
! NGRBSTR  - 177 Surface thermal radiation
! NGRBTSR  - 178 Top solar radiation
! NGRBTTR  - 179 Top thermal radiation
! NGRBEWSS - 180 U-stress
! NGRBNSSS - 181 V-stress
! NGRBE    - 182 Evaporation
! NGRBPEV  - 228251 Potential evaporation
! NGRBSTL3 - 183 Soil temperature level 3
! NGRBCCC  - 185 Convective cloud cocer
! NGRBLCC  - 186 Low cloud cover
! NGRBMCC  - 187 Medium cloud cover
! NGRBHCC  - 188 High cloud cover
! NGRBSUND - 189 Sunshine duration
! NGRBEWOV - 190 EW component of sub-grid scale orographic variance
! NGRBNSOV - 191 NS component of sub-grid scale orographic variance
! NGRBNWOV - 192 NWSE component of sub-grid scale orographic variance
! NGRBNEOV - 193 NESW component of sub-grid scale orographic variance
! NGRBBTMP - 194 Brightness temperature (K)
! NGRBCLBT - 260510 Cloudy brightness temperature
! NGRBCSBT - 260511 Clear-sky brightness temperature
! NGRBLGWS - 195 Latitudinal component of gravity wave stress
! NGRBMGWS - 196 Meridional component of gravity wave stress
! NGRBGWD  - 197 Gravity wave dissipation
! NGRBSRC  - 198 Skin reservoir content
! NGRBVEG  - 199 Percentage of vegetation
! NGRBVSO  - 200 variance of sub-grid scale orogrophy
! NGRBMX2T - 201 Maximum temperature at 2m since last post-processing
! NGRBMN2T - 202 Minimum temperature at 2m since last post-processing
! NGRBO3   - 203 Ozone mixing ratio (EC prognostic ozone)
! NGRBPAW  - 204 Precipitation analysis weights
! NGRBRO   - 205 Runoff
! NGRBTCO3 - 206 Total column ozone 
! NGRB10SI - 207 10m wind speed
! NGRBTSRC - 208 Top solar radiation clear sky 
! NGRBTTRC - 209 Top thermal radiation clear sky
! NGRBSSRC - 210 Surface solar radiation clear sky
! NGRBSTRC - 211 Surface thermal radiation clear sky

! NGRBTISR - 212 TOA incident solar radiation


!-- bunch of codes, confusing in their use ...
! NGRBDHR  - 214 Diabatic heating by radiation ?
! NGRBDHVD - 215 Diabatic heating by vertical diffusion?
! NGRBDHCC - 216 Diabatic heating by cumulus convection?
! NGRBDHLC - 217 Diabatic heating by large-scale condensation?
! NGRBVDZW - 218 INTSURFTEMPERATU (vertical diffusion of zonal wind?)
! NGRBVDMW - 219 PROFTEMPERATURE (vertical diffusion of meridional wind?) 

! NGRB221  - 221 Not used
! NGRBCTZW - 222 Convective tendency of zonal wind
! NGRBCTMW - 223 Convective tendency of meridional wind
! NGRBVDH  - 224 Not used (Vertical diffusion of humidity)
! NGRBHTCC - 225 Not used (Humidity tendency by cumulus convection)
! NGRBHTLC - 226 Not used (Humidity tendency by large-scale condensation)
! NGRBCRNH - 227 Not used (Change from removal of negative humidity)
! NGRBTP   - 228 Total precipitation
! NGRBIEWS - 229 Istantaneous X-surface stress
! NGRBINSS - 230 Istantaneous Y-surface stress
! NGRBISHF - 231 Istantaneous surface heat flux
! NGRBIE   - 232 Istantaneous moisture flux (evaporation)
! NGRBLSRH - 234 Logarithm of surface roughness length for heat
! NGRBSKT  - 235 Skin temperature
! NGRBSTL4 - 236 Soil temperature level 4
! NGRBTSN  - 238 Temperature of snow layer
! NGRBCSF  - 239 Convective snow-fall
! NGRBLSF  - 240 Large scale snow-fall
! NGRBTPR    - 260048 Total precipitation rate
! NGRBILSPF  - 228217 Large-scale precipitation fraction
! NGRBCRR    - 228218 Convective rain rate
! NGRBLSRR   - 228219 Large scale rain rate
! NGRBCSFR   - 228220 Convective snowfall rate water equivalent
! NGRBLSSFR  - 228221 Large scale snowfall rate water equivalent
! NGRBMXTPR3 - 228222 Max precip rate in last 3 hours
! NGRBMNTPR3 - 228223 Min precip rate in last 3 hours
! NGRBMXTPR6 - 228224 Max precip rate in last 6 hours
! NGRBMNTPR6 - 228225 Min precip rate in last 6 hours
! NGRBMXTPR  - 228226 Max precip rate since last post-processing
! NGRBMNTPR  - 228227 Min precip rate since last post-processing
! NGRBPTYPE  - 260015 Precipitation type

! NGRBACF  - 241 Not used
! NGRBALW  - 242 Not used
! NGRBFAL  - 243 Forecast albedo
! NGRBFSR  - 244 Forecast surface roughness
! NGRBFLSR - 245 Forecast logarithm of surface roughness for heat
! NGRBCLWC - 246 Cloud liquid water content
! NGRBCIWC - 247 Cloud ice water content
! NGRBCC   - 248 Cloud cover
! NGRB100U - 228246 100m u wind
! NGRB100V - 228247 100m v wind
! NGRB100SI - 228249 100m wind speed
! NGRB200U - 228239 200m u wind
! NGRB200V - 228240 200m v wind
! NGRB200SI - 228241 200m wind speed
! NGRBDNDZN - 228015 Minimum refractivity gradient inside trapping layer
! NGRBDNDZA - 228016 Mean refractivity gradient inside trapping layer
! NGRBDCTB  - 228017 Duct base height
! NGRBTPLB  - 228018 Trapping layer base height
! NGRBTPLT  - 228019 Trapping layer top height
! NGRBSSRDC - 228129 Surface clear-sky downward solar radiation 
! NGRBSTRDC - 228130 Surface clear-sky downward thermal radiation

!-- land carbon dioxide fields in Table 228 ---------------------
! NGRBFASGPPCOEF - 228078 Gross Primary Production flux adjustment factor
! NGRBFASRECCOEF - 228079 Ecosystem respiration flux adjustment factor
! NGRBNEE  - 228080 Net ecosysten exchange for CO2
! NGRBGPP  - 228081 Gross primary production for CO2
! NGRBREC  - 228082 Ecosystem respiration for CO2
! NGRBINEE - 228083 Istantaneous net ecosysten exchange for CO2
! NGRBIGPP - 228084 Istantaneous gross primary production for CO2
! NGRBIREC - 228085 Istantaneous ecosysten respiration for CO2

! NGRBCO2TYP- 129172 CO2 photosynthesis type (C3/C4) 

!-- ocean fields in Table 151 ---------------------
! NGRBOCT   - 128 *In situ* ocean temperature
! NGRBOCS   - 130 salinity
! NGRBOCU   - 131 zonal velocity U
! NGRBOCV   - 132 meridional velocity V
! NGRBOCVVS - 135 viscosity
! NGRBOCVDF - 136 diffusibity
! NGRBOCDEP - 137 bathymetry (bottom layer depth)
! NGRBOCLDP - 176 layer thickness (scalar & vector)
! NGRBOCLZ  - 213 layer depth
! NGRBADVT  - 214 correction term for temperature 
! NGRBADVS  - 215 correction term for salinity

!-- aerosols in Table 210 -------------------------
! NGRBAERMR01 - 001 aerosol mixing ratio 1
! NGRBAERMR02 - 002 aerosol mixing ratio 2
! NGRBAERMR03 - 003 aerosol mixing ratio 3
! NGRBAERMR04 - 004 aerosol mixing ratio 4
! NGRBAERMR05 - 005 aerosol mixing ratio 5
! NGRBAERMR06 - 006 aerosol mixing ratio 6
! NGRBAERMR07 - 007 aerosol mixing ratio 7
! NGRBAERMR08 - 008 aerosol mixing ratio 8
! NGRBAERMR09 - 009 aerosol mixing ratio 9
! NGRBAERMR10 - 010 aerosol mixing ratio 10
! NGRBAERMR11 - 011 aerosol mixing ratio 11
! NGRBAERMR12 - 012 aerosol mixing ratio 12
! NGRBAERMR13 - 013 aerosol mixing ratio 13
! NGRBAERMR14 - 014 aerosol mixing ratio 14
! NGRBAERMR15 - 015 aerosol mixing ratio 15
! NGRBAERGN01 - 016 aerosol gain acc. 1 2D
! NGRBAERGN02 - 017 aerosol gain acc. 2 2D
! NGRBAERGN03 - 018 aerosol gain acc. 3 2D
! NGRBAERGN04 - 019 aerosol gain acc. 4 2D
! NGRBAERGN05 - 020 aerosol gain acc. 5 2D
! NGRBAERGN06 - 021 aerosol gain acc. 6 2D
! NGRBAERGN07 - 022 aerosol gain acc. 7 2D
! NGRBAERGN08 - 023 aerosol gain acc. 8 2D
! NGRBAERGN09 - 024 aerosol gain acc. 9 2D
! NGRBAERGN10 - 025 aerosol gain acc. 10 2D
! NGRBAERGN11 - 026 aerosol gain acc. 11 2D
! NGRBAERGN12 - 027 aerosol gain acc. 12 2D
! NGRBAERGN13 - 028 aerosol gain acc. 13 2D
! NGRBAERGN14 - 029 aerosol gain acc. 14 2D
! NGRBAERGN15 - 030 aerosol gain acc. 15 2D

! NGRBAERLS01 - 031 black carbon biog.      clim2D 
! NGRBAERLS02 - 032 black carbon fossil     clim2D
! NGRBAERLS03 - 033 organic matter biog.    clim2D
! NGRBAERLS04 - 034 organic amtter fossil   clim2D
! NGRBAERLS05 - 035 sulphate low-level      clim2D
! NGRBAERLS06 - 036 sulphate high-level     clim2D
! NGRBAERLS07 - 037 volcanic continuous     clim2D
! NGRBAERLS08 - 038 volcanic explosive      clim2D
! NGRBAERLS09 - 039 secondary organic       clim2D
! NGRBAERLS10 - 040 black carbon GFED           2D
! NGRBAERLS11 - 041 organic matter GFED         2D
! NGRBAERLS12 - 042 sulphate GFED               2D
! NGRBAERLS13 - 043 oceanic DMS             clim2D

! NGRBAERLS14 - 044 aerosol loss acc. 14        2D
! NGRBAERLS15 - 045 aerosol loss acc. 15        2D
! NGRBAERPR   - 046 aerosol precursor mixing ratio 
! NGRBAERSM   - 047 small aerosols mixing ratio  
! NGRBAERLG   - 048 large aerosols mixing ratio ! This is used for the total mixing 
!                                                 ratio if JB_STRUCT%JB_DATA%NAEROCV=1
! NGRBAODPR   - 049 aerosol precursor opt.depth 2D
! NGRBAODSM   - 050 small aerosols opt. depth   2D
! NGRBAODLG   - 051 large aerosols opt. depth   2D
! NGRBAERDEP  - 052 dust emission potential clim2D 
! NGRBAERLTS  - 053 lifting threshold speed clim2D
! NGRBAERSCC  - 054 soli clay content       clim2D

! NGRBAEPM1   - 072 PM1 particulate matter <= 1um   2D
! NGRBAEPM25  - 073 PM2.5 particulate matter <= 2.5um 2D
! NGRBAEPM10  - 074 PM10 particulate matter <= 10um  2D

! NGRBAEODTO  - 207 aerosol total optical depth 2D
! NGRBAEODSS  - 208 optical depth sea salt      2D
! NGRBAEODDU  - 209 optical depth dust          2D
! NGRBAEODOM  - 210 optical depth organic matter2D
! NGRBAEODBC  - 211 optical depth black carbon  2D
! NGRBAEODSU  - 212 optical depth sulphate      2D
! NGRBAEODNI  - 250 optical depth Nitrate      2D
! NGRBAEODAM  - 251 optical depth Ammonium      2D

! NGRBAEODTO469  - 213 aerosol total optical depth 469 nm 2D
! NGRBAEODTO670  - 214 aerosol total optical depth 670 nm 2D
! NGRBAEODTO865  - 215 aerosol total optical depth 865 nm 2D
! NGRBAEODTO1240  - 216 aerosol total optical depth 1240 nm 2D
 
! NGRBAEODVSU - 243 optical depth volc.sulphate 2D
! NGRBAEODVFA - 244 optical depth volc.fly ash  2D

! NGRBAERMR16 - 247 aerosol mixing ratio 16
! NGRBAERMR17 - 248 aerosol mixing ratio 17
! NGRBAERMR18 - 249 aerosol mixing ratio 18
! NGRBAERMR19 - 252 aerosol mixing ratio 19
! NGRBAERMR20 - 253 aerosol mixing ratio 20

!-- aerosols in Table 215 -------------------------
! NGRBAEODSOA  - 215226 optical depth secondary organics     2D

!-- aerosols in Table 216 -------------------------

! NGRBAERSO2DD - 216006 dry deposition velocities from SUMO for SO2 2D
! NGRBAERFCA1  - 216043 calcite fraction of dust bin 1
! NGRBAERFCA2  - 216044 calcite fraction of dust bin 2
! NGRBAERDSF   - 216046 Dust source function       clim2D
! NGRBAERDSZ   - 216048 Dust size modulation at emission clim2D
! NGRBAERURBF  - 216122 Urban fraction

! NGRBUVBED     - 214002 UV Biologically Effective Dose 2D
! NGRBUVBEDCS   - 214003 UV Biologically Effective Dose - Clear Sky 2D

! ERA40 DIAGNOSTIC FIELDS
! NGRBMINERA    - 162100  Minimum identifier for ERA40 diagnostic fields
! NGRBMAXERA    - 162113  Maximum identifier for ERA40 diagnostic fields
! NGRBVIWVE     - 162071  Vertical integral of eastward water vapour flux
! NGRBVIWVN     - 162072  Vertical integral of northward water vapour flux

! NGRBCO2OF     - 210067 CO2 - ocean flux
! NGRBCO2NBF    - 210068 CO2 - biosphere flux
! NGRBCO2APF    - 210069 CO2 - anthropogenic emissions
! NGRBCO2FIRE   - 210080 CO2 - biomass burning
! NGRBCH4F      - 210070 CH4 surface fluxes - aggregated field
! NGRBCH4FIRE   - 210082 CH4 - fire emissions
! NGRBCH4       - 219004 CH4 emissions (wetlands)
!
! Lightning fields
! NGRBLITOTI  - 228050 Instantaneous total lightning flash density
! NGRBLITOTA1 - 228051 1h averaged total lightning flash density
! NGRBLITOTA3 - 228057 3h averaged total lightning flash density
! NGRBLITOTA6 - 228058 6h averaged total lightning flash density
! NGRBLICGI   - 228052 Instantaneous cloud-to-ground lightning flash density
! NGRBLICGA1  - 228053 1h averaged cloud-to-ground lightning flash density
! NGRBLICGA3  - 228059 3h averaged cloud-to-ground lightning flash density
! NGRBLICGA6  - 228060 6h averaged cloud-to-ground lightning flash density

! Precipitation type most frequent (mode) and most severe
! NGRBPTYPESEVR1 - 260318 Precipitation type (most severe) in the last 1 hour 
! NGRBPTYPESEVR3 - 260319 Precipitation type (most severe) in the last 3 hours
! NGRBPTYPESEVR6 - 260338 Precipitation type (most severe) in the last 6 hours
! NGRBPTYPEMODE1 - 260320 Precipitation type (most frequent) in the last 1 hour 
! NGRBPTYPEMODE3 - 260321 Precipitation type (most frequent) in the last 3 hours
! NGRBPTYPEMODE6 - 260339 Precipitation type (most frequent) in the last 6 hours

! Ocean fields
! NGRBICETK  - 174098 Ice thickness
! NGRBMLD    - 151148 Mixed layer depth
! NGRBSL     - 151145 Sea lavel
! NGRB20D    - 151163 20 degree isotherm depth
! NGRBSSS0   - 151130 Sea surface salinity
! NGRBTEM300 - 151164 Mean temparature over 300m
! NGRBSAL300 - 151175 Mean salinity content over 300m

!---------------------------------------------------
! NGRBGHG(JPGHG)   - 210061 GHG1: Carbon dioxide
!                  - 210062 GHG2: Methane
!                  - 210063 GHG3: Nitrous oxide
! NGRBTCGHG(JPGHG) - 210064 Total column GHG1: Carbon Dioxide
!                  - 210065 Total column GHG2: Methane
!                  - 210066 Total column GHG3: Nitrous Oxide

!---------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NGRBSTRF  =  1
INTEGER(KIND=JPIM), PARAMETER :: NGRBVP    =  2
INTEGER(KIND=JPIM), PARAMETER :: NGRBPT    =  3
INTEGER(KIND=JPIM), PARAMETER :: NGRBSRO   =  8
INTEGER(KIND=JPIM), PARAMETER :: NGRBSSRO  =  9
INTEGER(KIND=JPIM), PARAMETER :: NGRBALUVP = 15
INTEGER(KIND=JPIM), PARAMETER :: NGRBALUVD = 16
INTEGER(KIND=JPIM), PARAMETER :: NGRBALNIP = 17
INTEGER(KIND=JPIM), PARAMETER :: NGRBALNID = 18

INTEGER(KIND=JPIM), PARAMETER :: NGRBALUVI = 210186
INTEGER(KIND=JPIM), PARAMETER :: NGRBALNII = 210189
INTEGER(KIND=JPIM), PARAMETER :: NGRBALUVV = 210187
INTEGER(KIND=JPIM), PARAMETER :: NGRBALNIV = 210190
INTEGER(KIND=JPIM), PARAMETER :: NGRBALUVG = 210188
INTEGER(KIND=JPIM), PARAMETER :: NGRBALNIG = 210191

INTEGER(KIND=JPIM), PARAMETER :: NGRBPARCS = 20
INTEGER(KIND=JPIM), PARAMETER :: NGRBUCTP  = 21
INTEGER(KIND=JPIM), PARAMETER :: NGRBUCLN  = 22
INTEGER(KIND=JPIM), PARAMETER :: NGRBUCDV  = 23
INTEGER(KIND=JPIM), PARAMETER :: NGRBCL    = 26
INTEGER(KIND=JPIM), PARAMETER :: NGRBCVL   = 27
INTEGER(KIND=JPIM), PARAMETER :: NGRBCO2TYP= 129172
INTEGER(KIND=JPIM), PARAMETER :: NGRBCVH   = 28
INTEGER(KIND=JPIM), PARAMETER :: NGRBFWET  = 200026
INTEGER(KIND=JPIM), PARAMETER :: NGRBCUR   = 200199
INTEGER(KIND=JPIM), PARAMETER :: NGRBTVL   = 29
INTEGER(KIND=JPIM), PARAMETER :: NGRBTVH   = 30
INTEGER(KIND=JPIM), PARAMETER :: NGRBCI    = 31
INTEGER(KIND=JPIM), PARAMETER :: NGRBASN   = 32
INTEGER(KIND=JPIM), PARAMETER :: NGRBRSN   = 33
INTEGER(KIND=JPIM), PARAMETER :: NGRBSSTK  = 34
INTEGER(KIND=JPIM), PARAMETER :: NGRBISTL1 = 35
INTEGER(KIND=JPIM), PARAMETER :: NGRBISTL2 = 36
INTEGER(KIND=JPIM), PARAMETER :: NGRBISTL3 = 37
INTEGER(KIND=JPIM), PARAMETER :: NGRBISTL4 = 38
INTEGER(KIND=JPIM), PARAMETER :: NGRBSWVL1 = 39
INTEGER(KIND=JPIM), PARAMETER :: NGRBSWVL2 = 40
INTEGER(KIND=JPIM), PARAMETER :: NGRBSWVL3 = 41
INTEGER(KIND=JPIM), PARAMETER :: NGRBSWVL4 = 42
INTEGER(KIND=JPIM), PARAMETER :: NGRBSLT   = 43
INTEGER(KIND=JPIM), PARAMETER :: NGRBES    = 44
INTEGER(KIND=JPIM), PARAMETER :: NGRBSMLT  = 45
INTEGER(KIND=JPIM), PARAMETER :: NGRB10FG  = 49
INTEGER(KIND=JPIM), PARAMETER :: NGRBLSPF  = 50

INTEGER(KIND=JPIM), PARAMETER :: NGRBMONT  = 53
INTEGER(KIND=JPIM), PARAMETER :: NGRBPRES  = 54

INTEGER(KIND=JPIM), PARAMETER :: NGRBUVB   = 57
INTEGER(KIND=JPIM), PARAMETER :: NGRBPAR   = 58
INTEGER(KIND=JPIM), PARAMETER :: NGRBCAPE  = 59
INTEGER(KIND=JPIM), PARAMETER :: NGRBMUCAPE= 228235
INTEGER(KIND=JPIM), PARAMETER :: NGRBMLCAPE50 = 228231
INTEGER(KIND=JPIM), PARAMETER :: NGRBMLCAPE100= 228233
INTEGER(KIND=JPIM), PARAMETER :: NGRBMUDEPL = 228237
INTEGER(KIND=JPIM), PARAMETER :: NGRBCAPES = 228044
INTEGER(KIND=JPIM), PARAMETER :: NGRBMXCAP6= 228035
INTEGER(KIND=JPIM), PARAMETER :: NGRBMXCAPS6=228036
INTEGER(KIND=JPIM), PARAMETER :: NGRBPV    = 60
INTEGER(KIND=JPIM), PARAMETER :: NGRBLAIL  = 66
INTEGER(KIND=JPIM), PARAMETER :: NGRBLAIH  = 67

INTEGER(KIND=JPIM), PARAMETER :: NGRBSDFOR = 74

INTEGER(KIND=JPIM), PARAMETER :: NGRBCRWC  = 75
INTEGER(KIND=JPIM), PARAMETER :: NGRBCSWC  = 76

INTEGER(KIND=JPIM), PARAMETER :: NGRBETADOT= 77
INTEGER(KIND=JPIM), PARAMETER :: NGRBTCLW  = 78
INTEGER(KIND=JPIM), PARAMETER :: NGRBTCIW  = 79
INTEGER(KIND=JPIM), PARAMETER :: NGRBTCSLW = 228088
INTEGER(KIND=JPIM), PARAMETER :: NGRBTCRW  = 228089
INTEGER(KIND=JPIM), PARAMETER :: NGRBTCSW  = 228090

INTEGER(KIND=JPIM), PARAMETER :: NGRBTROPOTP= 228045

! LARPEGE
INTEGER(KIND=JPIM), PARAMETER :: NGRBSPD   = 80
INTEGER(KIND=JPIM), PARAMETER :: NGRBSVD   = 81
! LECMWF
INTEGER(KIND=JPIM), PARAMETER :: NGRB080   = 80
INTEGER(KIND=JPIM), PARAMETER :: NGRB081   = 81
! LECV
INTEGER(KIND=JPIM), PARAMETER :: NGRBALPHA= 90

INTEGER(KIND=JPIM), PARAMETER :: NGRBMINXTRA  = 082
INTEGER(KIND=JPIM), PARAMETER :: NGRBMAXXTRA  = 117

INTEGER(KIND=JPIM), PARAMETER :: NGRB082  = 082
INTEGER(KIND=JPIM), PARAMETER :: NGRB083  = 083
INTEGER(KIND=JPIM), PARAMETER :: NGRB084  = 084
INTEGER(KIND=JPIM), PARAMETER :: NGRB085  = 085
INTEGER(KIND=JPIM), PARAMETER :: NGRB086  = 086
INTEGER(KIND=JPIM), PARAMETER :: NGRB087  = 087
INTEGER(KIND=JPIM), PARAMETER :: NGRB088  = 088
INTEGER(KIND=JPIM), PARAMETER :: NGRB089  = 089
INTEGER(KIND=JPIM), PARAMETER :: NGRB090  = 090
INTEGER(KIND=JPIM), PARAMETER :: NGRB091  = 091
INTEGER(KIND=JPIM), PARAMETER :: NGRB092  = 092
INTEGER(KIND=JPIM), PARAMETER :: NGRB093  = 093
INTEGER(KIND=JPIM), PARAMETER :: NGRB094  = 094
INTEGER(KIND=JPIM), PARAMETER :: NGRB095  = 095
INTEGER(KIND=JPIM), PARAMETER :: NGRB096  = 096
INTEGER(KIND=JPIM), PARAMETER :: NGRB097  = 097
INTEGER(KIND=JPIM), PARAMETER :: NGRB098  = 098
INTEGER(KIND=JPIM), PARAMETER :: NGRB099  = 099
INTEGER(KIND=JPIM), PARAMETER :: NGRB100  = 100
INTEGER(KIND=JPIM), PARAMETER :: NGRB101  = 101
INTEGER(KIND=JPIM), PARAMETER :: NGRB102  = 102
INTEGER(KIND=JPIM), PARAMETER :: NGRB103  = 103
INTEGER(KIND=JPIM), PARAMETER :: NGRB104  = 104
INTEGER(KIND=JPIM), PARAMETER :: NGRB105  = 105
INTEGER(KIND=JPIM), PARAMETER :: NGRB106  = 106
INTEGER(KIND=JPIM), PARAMETER :: NGRB107  = 107
INTEGER(KIND=JPIM), PARAMETER :: NGRB108  = 108
INTEGER(KIND=JPIM), PARAMETER :: NGRB109  = 109
INTEGER(KIND=JPIM), PARAMETER :: NGRB110  = 110
INTEGER(KIND=JPIM), PARAMETER :: NGRB111  = 111
INTEGER(KIND=JPIM), PARAMETER :: NGRB112  = 112
INTEGER(KIND=JPIM), PARAMETER :: NGRB113  = 113
INTEGER(KIND=JPIM), PARAMETER :: NGRB114  = 114
INTEGER(KIND=JPIM), PARAMETER :: NGRB115  = 115
INTEGER(KIND=JPIM), PARAMETER :: NGRB116  = 116
INTEGER(KIND=JPIM), PARAMETER :: NGRB117  = 117

INTEGER(KIND=JPIM), PARAMETER :: NGRB118  = 118
INTEGER(KIND=JPIM), PARAMETER :: NGRB119  = 119
INTEGER(KIND=JPIM), PARAMETER :: NGRB120  = 120

INTEGER(KIND=JPIM), PARAMETER :: NGRBMX2T3 = 228026
INTEGER(KIND=JPIM), PARAMETER :: NGRBMN2T3 = 228027
INTEGER(KIND=JPIM), PARAMETER :: NGRB10FG3 = 228028
INTEGER(KIND=JPIM), PARAMETER :: NGRBI10FG = 228029
INTEGER(KIND=JPIM), PARAMETER :: NGRBMX2T6 = 121
INTEGER(KIND=JPIM), PARAMETER :: NGRBMN2T6 = 122
INTEGER(KIND=JPIM), PARAMETER :: NGRB10FG6 = 123
INTEGER(KIND=JPIM), PARAMETER :: NGRBEMIS  = 124

INTEGER(KIND=JPIM), PARAMETER :: NGRBAT   = 127
INTEGER(KIND=JPIM), PARAMETER :: NGRBBV   = 128
INTEGER(KIND=JPIM), PARAMETER :: NGRBZ    = 129
INTEGER(KIND=JPIM), PARAMETER :: NGRBT    = 130
INTEGER(KIND=JPIM), PARAMETER :: NGRBU    = 131
INTEGER(KIND=JPIM), PARAMETER :: NGRBV    = 132
INTEGER(KIND=JPIM), PARAMETER :: NGRBUCUR = 151131
INTEGER(KIND=JPIM), PARAMETER :: NGRBVCUR = 151132
INTEGER(KIND=JPIM), PARAMETER :: NGRBSSS  = 151130
INTEGER(KIND=JPIM), PARAMETER :: NGRBQ    = 133
INTEGER(KIND=JPIM), PARAMETER :: NGRBSP   = 134
INTEGER(KIND=JPIM), PARAMETER :: NGRBW    = 135
INTEGER(KIND=JPIM), PARAMETER :: NGRBTCW  = 136
INTEGER(KIND=JPIM), PARAMETER :: NGRBTCWV = 137
INTEGER(KIND=JPIM), PARAMETER :: NGRBVO   = 138
INTEGER(KIND=JPIM), PARAMETER :: NGRBSTL1 = 139

INTEGER(KIND=JPIM), PARAMETER :: NGRBSDSL = 141  ! back-comp single-layer
!!INTEGER(KIND=JPIM), PARAMETER :: NGRBSD   = 3066 ! grib2 multi-layer
INTEGER(KIND=JPIM), PARAMETER :: NGRBSD   = 228141 ! grib2 multi-layer
INTEGER(KIND=JPIM), PARAMETER :: NGRBWSN  = 228038

INTEGER(KIND=JPIM), PARAMETER :: NGRBLSP  = 142
INTEGER(KIND=JPIM), PARAMETER :: NGRBCP   = 143
INTEGER(KIND=JPIM), PARAMETER :: NGRBSF   = 144
INTEGER(KIND=JPIM), PARAMETER :: NGRBFZRA = 228216
INTEGER(KIND=JPIM), PARAMETER :: NGRBBLD  = 145
INTEGER(KIND=JPIM), PARAMETER :: NGRBSSHF = 146
INTEGER(KIND=JPIM), PARAMETER :: NGRBSLHF = 147
INTEGER(KIND=JPIM), PARAMETER :: NGRBCHNK = 148
INTEGER(KIND=JPIM), PARAMETER :: NGRBSNR  = 149
INTEGER(KIND=JPIM), PARAMETER :: NGRBTNR  = 150
INTEGER(KIND=JPIM), PARAMETER :: NGRBMSL  = 151
INTEGER(KIND=JPIM), PARAMETER :: NGRBLNSP = 152
INTEGER(KIND=JPIM), PARAMETER :: NGRBSWHR = 153
INTEGER(KIND=JPIM), PARAMETER :: NGRBLWHR = 154
INTEGER(KIND=JPIM), PARAMETER :: NGRBD    = 155
INTEGER(KIND=JPIM), PARAMETER :: NGRBGH   = 156
INTEGER(KIND=JPIM), PARAMETER :: NGRBR    = 157
INTEGER(KIND=JPIM), PARAMETER :: NGRBTSP  = 158
INTEGER(KIND=JPIM), PARAMETER :: NGRBBLH  = 159
INTEGER(KIND=JPIM), PARAMETER :: NGRBSDOR = 160
INTEGER(KIND=JPIM), PARAMETER :: NGRBISOR = 161
INTEGER(KIND=JPIM), PARAMETER :: NGRBANOR = 162
INTEGER(KIND=JPIM), PARAMETER :: NGRBSLOR = 163
INTEGER(KIND=JPIM), PARAMETER :: NGRBTCC  = 164
INTEGER(KIND=JPIM), PARAMETER :: NGRB10U  = 165
INTEGER(KIND=JPIM), PARAMETER :: NGRB10V  = 166
INTEGER(KIND=JPIM), PARAMETER :: NGRBZUST = 228003
INTEGER(KIND=JPIM), PARAMETER :: NGRBFDIR = 228021
INTEGER(KIND=JPIM), PARAMETER :: NGRBCDIR = 228022
INTEGER(KIND=JPIM), PARAMETER :: NGRBDSRP = 47
INTEGER(KIND=JPIM), PARAMETER :: NGRBCBASE= 228023
INTEGER(KIND=JPIM), PARAMETER :: NGRB0DEGL= 228024
INTEGER(KIND=JPIM), PARAMETER :: NGRBM10DEGL= 228020
INTEGER(KIND=JPIM), PARAMETER :: NGRBCIN  = 228001
INTEGER(KIND=JPIM), PARAMETER :: NGRBMLCIN50 = 228232 
INTEGER(KIND=JPIM), PARAMETER :: NGRBMLCIN100= 228234
INTEGER(KIND=JPIM), PARAMETER :: NGRBKINDEX = 260121
INTEGER(KIND=JPIM), PARAMETER :: NGRBTTINDEX= 260123
INTEGER(KIND=JPIM), PARAMETER :: NGRBCBASEA= 260109
INTEGER(KIND=JPIM), PARAMETER :: NGRBCTOPC = 228046
INTEGER(KIND=JPIM), PARAMETER :: NGRBZTWETB0 = 228047
INTEGER(KIND=JPIM), PARAMETER :: NGRBZTWETB1 = 228048
INTEGER(KIND=JPIM), PARAMETER :: NGRBVISIH= 3020
INTEGER(KIND=JPIM), PARAMETER :: NGRB10NU = 228131
INTEGER(KIND=JPIM), PARAMETER :: NGRB10NV = 228132
INTEGER(KIND=JPIM), PARAMETER :: NGRB100U = 228246
INTEGER(KIND=JPIM), PARAMETER :: NGRB100V = 228247
INTEGER(KIND=JPIM), PARAMETER :: NGRB100SI= 228249
INTEGER(KIND=JPIM), PARAMETER :: NGRB200U = 228239
INTEGER(KIND=JPIM), PARAMETER :: NGRB200V = 228240
INTEGER(KIND=JPIM), PARAMETER :: NGRB200SI= 228241
INTEGER(KIND=JPIM), PARAMETER :: NGRB2T   = 167
INTEGER(KIND=JPIM), PARAMETER :: NGRB2D   = 168
INTEGER(KIND=JPIM), PARAMETER :: NGRB2SH  = 174096
INTEGER(KIND=JPIM), PARAMETER :: NGRBSSRD = 169
INTEGER(KIND=JPIM), PARAMETER :: NGRBSTL2 = 170
INTEGER(KIND=JPIM), PARAMETER :: NGRBLSM  = 172
INTEGER(KIND=JPIM), PARAMETER :: NGRBSR   = 173
INTEGER(KIND=JPIM), PARAMETER :: NGRBAL   = 174
INTEGER(KIND=JPIM), PARAMETER :: NGRBSTRD = 175
INTEGER(KIND=JPIM), PARAMETER :: NGRBSSR  = 176
INTEGER(KIND=JPIM), PARAMETER :: NGRBSTR  = 177
INTEGER(KIND=JPIM), PARAMETER :: NGRBTSR  = 178
INTEGER(KIND=JPIM), PARAMETER :: NGRBTTR  = 179
INTEGER(KIND=JPIM), PARAMETER :: NGRBEWSS = 180
INTEGER(KIND=JPIM), PARAMETER :: NGRBNSSS = 181
INTEGER(KIND=JPIM), PARAMETER :: NGRBE    = 182
INTEGER(KIND=JPIM), PARAMETER :: NGRBPEV  = 228251
INTEGER(KIND=JPIM), PARAMETER :: NGRBSTL3 = 183
INTEGER(KIND=JPIM), PARAMETER :: NGRBCCC  = 185
INTEGER(KIND=JPIM), PARAMETER :: NGRBLCC  = 186
INTEGER(KIND=JPIM), PARAMETER :: NGRBMCC  = 187
INTEGER(KIND=JPIM), PARAMETER :: NGRBHCC  = 188
INTEGER(KIND=JPIM), PARAMETER :: NGRBSUND = 189
INTEGER(KIND=JPIM), PARAMETER :: NGRBEWOV = 190
INTEGER(KIND=JPIM), PARAMETER :: NGRBNSOV = 191
INTEGER(KIND=JPIM), PARAMETER :: NGRBNWOV = 192
INTEGER(KIND=JPIM), PARAMETER :: NGRBNEOV = 193
INTEGER(KIND=JPIM), PARAMETER :: NGRBBTMP = 194
INTEGER(KIND=JPIM), PARAMETER :: NGRBCLBT = 260510
INTEGER(KIND=JPIM), PARAMETER :: NGRBCSBT = 260511
INTEGER(KIND=JPIM), PARAMETER :: NGRBLGWS = 195
INTEGER(KIND=JPIM), PARAMETER :: NGRBMGWS = 196
INTEGER(KIND=JPIM), PARAMETER :: NGRBGWD  = 197
INTEGER(KIND=JPIM), PARAMETER :: NGRBSRC  = 198
INTEGER(KIND=JPIM), PARAMETER :: NGRBVEG  = 199
INTEGER(KIND=JPIM), PARAMETER :: NGRBVSO  = 200
INTEGER(KIND=JPIM), PARAMETER :: NGRBMX2T = 201
INTEGER(KIND=JPIM), PARAMETER :: NGRBMN2T = 202
INTEGER(KIND=JPIM), PARAMETER :: NGRBO3   = 203
INTEGER(KIND=JPIM), PARAMETER :: NGRBPAW  = 204
INTEGER(KIND=JPIM), PARAMETER :: NGRBRO   = 205
INTEGER(KIND=JPIM), PARAMETER :: NGRBTCO3 = 206
INTEGER(KIND=JPIM), PARAMETER :: NGRB10SI = 207
INTEGER(KIND=JPIM), PARAMETER :: NGRBTSRC = 208
INTEGER(KIND=JPIM), PARAMETER :: NGRBTTRC = 209
INTEGER(KIND=JPIM), PARAMETER :: NGRBSSRC = 210
INTEGER(KIND=JPIM), PARAMETER :: NGRBSTRC = 211

INTEGER(KIND=JPIM), PARAMETER :: NGRBTISR = 212
INTEGER(KIND=JPIM), PARAMETER :: NGRBVIMD = 213

INTEGER(KIND=JPIM), PARAMETER :: NGRBDHR  = 214
INTEGER(KIND=JPIM), PARAMETER :: NGRBDHVD = 215
INTEGER(KIND=JPIM), PARAMETER :: NGRBDHCC = 216
INTEGER(KIND=JPIM), PARAMETER :: NGRBDHLC = 217
INTEGER(KIND=JPIM), PARAMETER :: NGRBVDZW = 218
INTEGER(KIND=JPIM), PARAMETER :: NGRBVDMW = 219
INTEGER(KIND=JPIM), PARAMETER :: NGRBCTZW = 222
INTEGER(KIND=JPIM), PARAMETER :: NGRBCTMW = 223
INTEGER(KIND=JPIM), PARAMETER :: NGRBVDH  = 224
INTEGER(KIND=JPIM), PARAMETER :: NGRBHTCC = 225
INTEGER(KIND=JPIM), PARAMETER :: NGRBHTLC = 226
INTEGER(KIND=JPIM), PARAMETER :: NGRBCRNH = 227
INTEGER(KIND=JPIM), PARAMETER :: NGRBTP   = 228
INTEGER(KIND=JPIM), PARAMETER :: NGRBIEWS = 229
INTEGER(KIND=JPIM), PARAMETER :: NGRBINSS = 230
INTEGER(KIND=JPIM), PARAMETER :: NGRBISHF = 231
INTEGER(KIND=JPIM), PARAMETER :: NGRBIE   = 232
INTEGER(KIND=JPIM), PARAMETER :: NGRBLSRH = 234
INTEGER(KIND=JPIM), PARAMETER :: NGRBSKT  = 235
INTEGER(KIND=JPIM), PARAMETER :: NGRBSTL4 = 236
INTEGER(KIND=JPIM), PARAMETER :: NGRBTSN  = 238
INTEGER(KIND=JPIM), PARAMETER :: NGRBCSF  = 239
INTEGER(KIND=JPIM), PARAMETER :: NGRBLSF  = 240
INTEGER(KIND=JPIM), PARAMETER :: NGRBTPR    = 260048
INTEGER(KIND=JPIM), PARAMETER :: NGRBILSPF  = 228217
INTEGER(KIND=JPIM), PARAMETER :: NGRBCRR    = 228218
INTEGER(KIND=JPIM), PARAMETER :: NGRBLSRR   = 228219
INTEGER(KIND=JPIM), PARAMETER :: NGRBCSFR   = 228220
INTEGER(KIND=JPIM), PARAMETER :: NGRBLSSFR  = 228221
INTEGER(KIND=JPIM), PARAMETER :: NGRBMXTPR3 = 228222
INTEGER(KIND=JPIM), PARAMETER :: NGRBMNTPR3 = 228223
INTEGER(KIND=JPIM), PARAMETER :: NGRBMXTPR6 = 228224
INTEGER(KIND=JPIM), PARAMETER :: NGRBMNTPR6 = 228225
INTEGER(KIND=JPIM), PARAMETER :: NGRBMXTPR  = 228226
INTEGER(KIND=JPIM), PARAMETER :: NGRBMNTPR  = 228227
INTEGER(KIND=JPIM), PARAMETER :: NGRBPTYPE  = 260015

INTEGER(KIND=JPIM), PARAMETER :: NGRBACF  = 241
INTEGER(KIND=JPIM), PARAMETER :: NGRBALW  = 242
INTEGER(KIND=JPIM), PARAMETER :: NGRBFAL  = 243
INTEGER(KIND=JPIM), PARAMETER :: NGRBFSR  = 244
INTEGER(KIND=JPIM), PARAMETER :: NGRBFLSR = 245
INTEGER(KIND=JPIM), PARAMETER :: NGRBCLWC = 246
INTEGER(KIND=JPIM), PARAMETER :: NGRBCIWC = 247
INTEGER(KIND=JPIM), PARAMETER :: NGRBCC   = 248
INTEGER(KIND=JPIM), PARAMETER :: NGRBAIW  = 249
INTEGER(KIND=JPIM), PARAMETER :: NGRBICE  = 250
INTEGER(KIND=JPIM), PARAMETER :: NGRBATTE = 251
INTEGER(KIND=JPIM), PARAMETER :: NGRBATHE = 252
INTEGER(KIND=JPIM), PARAMETER :: NGRBATZE = 253
INTEGER(KIND=JPIM), PARAMETER :: NGRBATMW = 254
INTEGER(KIND=JPIM), PARAMETER :: NGRB255  = 255
!-- lake depth ancillary + lake prognostics -- Table 228 ---
INTEGER(KIND=JPIM) :: NGRBDL=228007
INTEGER(KIND=JPIM) :: NGRBLMLT=228008
INTEGER(KIND=JPIM) :: NGRBLMLD=228009
INTEGER(KIND=JPIM) :: NGRBLBLT=228010
INTEGER(KIND=JPIM) :: NGRBLTLT=228011
INTEGER(KIND=JPIM) :: NGRBLSHF=228012
INTEGER(KIND=JPIM) :: NGRBLICT=228013
INTEGER(KIND=JPIM) :: NGRBLICD=228014
!-- land carbon dioxide diagnostics -- Table 228 ---
INTEGER(KIND=JPIM) :: NGRBNEE  = 228080
INTEGER(KIND=JPIM) :: NGRBGPP  = 228081
INTEGER(KIND=JPIM) :: NGRBREC  = 228082
INTEGER(KIND=JPIM) :: NGRBINEE = 228083
INTEGER(KIND=JPIM) :: NGRBIGPP = 228084
INTEGER(KIND=JPIM) :: NGRBIREC = 228085
!-- ducting diagnostics -- Table 228 ---------------
INTEGER(KIND=JPIM), PARAMETER :: NGRBDNDZN = 228015
INTEGER(KIND=JPIM), PARAMETER :: NGRBDNDZA = 228016
INTEGER(KIND=JPIM), PARAMETER :: NGRBDCTB  = 228017
INTEGER(KIND=JPIM), PARAMETER :: NGRBTPLB  = 228018
INTEGER(KIND=JPIM), PARAMETER :: NGRBTPLT  = 228019
INTEGER(KIND=JPIM), PARAMETER :: NGRBASCAT_SM_CDFA  = 228253
INTEGER(KIND=JPIM), PARAMETER :: NGRBASCAT_SM_CDFB  = 228254
!-- clearsky downward surface fluxes -- Table 228 ---------------
INTEGER(KIND=JPIM), PARAMETER :: NGRBSSRDC= 228129 
INTEGER(KIND=JPIM), PARAMETER :: NGRBSTRDC= 228130 
!-- ocean mixed layer -- Table 151 ----------------
INTEGER(KIND=JPIM), PARAMETER :: NGRBOCT   = 151129 
INTEGER(KIND=JPIM), PARAMETER :: NGRBOCS   = 151130
INTEGER(KIND=JPIM), PARAMETER :: NGRBOCU   = 151131
INTEGER(KIND=JPIM), PARAMETER :: NGRBOCV   = 151132
INTEGER(KIND=JPIM), PARAMETER :: NGRBOCVVS = 105 
INTEGER(KIND=JPIM), PARAMETER :: NGRBOCVDF = 106 
INTEGER(KIND=JPIM), PARAMETER :: NGRBOCDEP = 107 
INTEGER(KIND=JPIM), PARAMETER :: NGRBOCLDP = 108 
INTEGER(KIND=JPIM), PARAMETER :: NGRBOCLZ  = 109 
INTEGER(KIND=JPIM), PARAMETER :: NGRBADVT  = 110 
INTEGER(KIND=JPIM), PARAMETER :: NGRBADVS  = 111
INTEGER(KIND=JPIM), PARAMETER :: NGRBOCUC  = 112
INTEGER(KIND=JPIM), PARAMETER :: NGRBOCVC  = 113
INTEGER(KIND=JPIM), PARAMETER :: NGRBUSTRC = 114
INTEGER(KIND=JPIM), PARAMETER :: NGRBVSTRC = 115

!-- aerosols -- Table 210 --------------------------
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR01=210001
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR02=210002
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR03=210003
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR04=210004
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR05=210005
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR06=210006
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR07=210007
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR08=210008
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR09=210009
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR10=210010
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR11=210011
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR12=210012
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR13=210013
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR14=210014
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR15=210015

INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN01=210016
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN02=210017
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN03=210018
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN04=210019
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN05=210020
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN06=210021
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN07=210022
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN08=210023
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN09=210024
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN10=210025
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN11=210026
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN12=210027
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN13=210028
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN14=210029
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERGN15=210030

INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS01  =210031
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS02  =210032
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS03  =210033
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS04  =210034
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS05  =210035
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS06  =210036
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS07  =210037
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS08  =210038
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS09  =210039
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS10  =210040
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS11  =210041
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS12  =210042
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS13  =210043

INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS14=210044
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLS15=210045

INTEGER(KIND=JPIM), PARAMETER :: NGRBAERPR  =210046
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERSM  =210047
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLG  =210048
INTEGER(KIND=JPIM), PARAMETER :: NGRBAODPR  =210049
INTEGER(KIND=JPIM), PARAMETER :: NGRBAODSM  =210050
INTEGER(KIND=JPIM), PARAMETER :: NGRBAODLG  =210051
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERDEP =210052
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERLTS =210053
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERSCC =210054

INTEGER(KIND=JPIM), PARAMETER :: NGRBAEPM1  =210072
INTEGER(KIND=JPIM), PARAMETER :: NGRBAEPM25 =210073
INTEGER(KIND=JPIM), PARAMETER :: NGRBAEPM10 =210074

INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODTO =210207
INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODSS =210208
INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODDU =210209
INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODOM =210210
INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODBC =210211
INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODSU =210212

INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODTO469=210213
INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODTO670=210214
INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODTO865=210215
INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODTO1240=210216

INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODVSU=210243
INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODVFA=210244

INTEGER(KIND=JPIM), PARAMETER :: NGRBTAEDEC550=210245
INTEGER(KIND=JPIM), PARAMETER :: NGRBTAEDAB550=210246

INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR16=210247
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR17=210248
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR18=210249

INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODNI =210250
INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODAM =210251

INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR19=210252
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERMR20=210253

!--Injection height for biomass burning emissions---
INTEGER(KIND=JPIM), PARAMETER :: NGRBINJFIRE = 210119

!-- aerosols -- Table 216 --------------------------
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERSO2DD =216006
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERFCA1  =216043
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERFCA2  =216044
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERDSF =216046
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERDSZ =216048
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERURBF  =216122

!---------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NGRBUVBED   = 214002
INTEGER(KIND=JPIM), PARAMETER :: NGRBUVBEDCS = 214003

INTEGER(KIND=JPIM), PARAMETER :: NGRBMINERA  = 162100 
INTEGER(KIND=JPIM), PARAMETER :: NGRBMAXERA  = 162113
INTEGER(KIND=JPIM), PARAMETER :: NGRBVIWVE   = 162071
INTEGER(KIND=JPIM), PARAMETER :: NGRBVIWVN   = 162072

!---------------------------------------------------
INTEGER(KIND=JPIM), PARAMETER :: NGRBCO2OF   = 210067
INTEGER(KIND=JPIM), PARAMETER :: NGRBCO2NBF  = 210068
INTEGER(KIND=JPIM), PARAMETER :: NGRBCO2APF  = 210069
INTEGER(KIND=JPIM), PARAMETER :: NGRBCO2FIRE = 210080
INTEGER(KIND=JPIM), PARAMETER :: NGRBCH4F    = 210070
INTEGER(KIND=JPIM), PARAMETER :: NGRBCH4FIRE = 210082
INTEGER(KIND=JPIM), PARAMETER :: NGRBCH4WET  = 228104
!---------------------------------------------------
INTEGER(KIND=JPIM), PARAMETER :: NGRBFASGPPCOEF = 228078 
INTEGER(KIND=JPIM), PARAMETER :: NGRBFASRECCOEF = 228079
!---------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NGRBNOXLOG = 210121    ! use NO2 for C-IFS
!INTEGER(KIND=JPIM), PARAMETER :: NGRBNOXLOG = 210129    ! use NOX for coupled system

INTEGER(KIND=JPIM), PARAMETER :: NGRBLITOTI  = 228050
INTEGER(KIND=JPIM), PARAMETER :: NGRBLITOTA1 = 228051
INTEGER(KIND=JPIM), PARAMETER :: NGRBLITOTA3 = 228057
INTEGER(KIND=JPIM), PARAMETER :: NGRBLITOTA6 = 228058
INTEGER(KIND=JPIM), PARAMETER :: NGRBLICGI   = 228052
INTEGER(KIND=JPIM), PARAMETER :: NGRBLICGA1  = 228053
INTEGER(KIND=JPIM), PARAMETER :: NGRBLICGA3  = 228059
INTEGER(KIND=JPIM), PARAMETER :: NGRBLICGA6  = 228060

INTEGER(KIND=JPIM), PARAMETER :: NGRBPTYPESEVR1 = 260318
INTEGER(KIND=JPIM), PARAMETER :: NGRBPTYPESEVR3 = 260319
INTEGER(KIND=JPIM), PARAMETER :: NGRBPTYPESEVR6 = 260338
INTEGER(KIND=JPIM), PARAMETER :: NGRBPTYPEMODE1 = 260320
INTEGER(KIND=JPIM), PARAMETER :: NGRBPTYPEMODE3 = 260321
INTEGER(KIND=JPIM), PARAMETER :: NGRBPTYPEMODE6 = 260339

!--Further aerosol diagnostics -- Table 215 --------
INTEGER(KIND=JPIM), PARAMETER :: NGRBACCAOD550 = 215089

INTEGER(KIND=JPIM), PARAMETER :: NGRBAOT532  = 215093
INTEGER(KIND=JPIM), PARAMETER :: NGRBNAOT532 = 215094
INTEGER(KIND=JPIM), PARAMETER :: NGRBAAOT532 = 215095

INTEGER(KIND=JPIM), PARAMETER :: NGRBAEREXT355  = 215180
INTEGER(KIND=JPIM), PARAMETER :: NGRBAEREXT532  = 215181
INTEGER(KIND=JPIM), PARAMETER :: NGRBAEREXT1064 = 215182

INTEGER(KIND=JPIM), PARAMETER :: NGRBAERBACKSCATTOA355  = 215183
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERBACKSCATTOA532  = 215184
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERBACKSCATTOA1064 = 215185

INTEGER(KIND=JPIM), PARAMETER :: NGRBAERBACKSCATGND355  = 215186
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERBACKSCATGND532  = 215187
INTEGER(KIND=JPIM), PARAMETER :: NGRBAERBACKSCATGND1064 = 215188

INTEGER(KIND=JPIM), PARAMETER :: NGRBAEODSOA =215226

!-- Ocean model output on IFS grid ---

INTEGER(KIND=JPIM), PARAMETER :: NGRBICETK = 174098
INTEGER(KIND=JPIM), PARAMETER :: NGRBMLD   = 151148
INTEGER(KIND=JPIM), PARAMETER :: NGRBSL    = 151145
INTEGER(KIND=JPIM), PARAMETER :: NGRB20D   = 151163
INTEGER(KIND=JPIM), PARAMETER :: NGRBSSSO  = 151130
INTEGER(KIND=JPIM), PARAMETER :: NGRBTEM300= 151164
INTEGER(KIND=JPIM), PARAMETER :: NGRBSAL300= 151175

INTEGER(KIND=JPIM), PARAMETER, DIMENSION(JPGHG) :: NGRBGHG = (/&
 & 210061, 210062, 210063/)
INTEGER(KIND=JPIM),  PARAMETER,DIMENSION(JPGHG) :: NGRBTCGHG = (/&
 & 210064, 210065, 210066/)
!     ------------------------------------------------------------------

END MODULE YOM_GRIB_CODES
