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

module varno_module
use parkind1,   only : jpim
implicit none

type odb_varno
   integer(kind=jpim) :: u                                    = 3   !    upper air u component     
   integer(kind=jpim) :: v                                    = 4   !    upper air v component     
   integer(kind=jpim) :: z                                    = 1   !    geopotential              
   integer(kind=jpim) :: dz                                   = 57  !    thickness                 
   integer(kind=jpim) :: rh                                   = 29  !    upper air rel. humidity   
   integer(kind=jpim) :: pwc                                  = 9   !    precipitable water content      
   integer(kind=jpim) :: rh2m                                 = 58  !    2m rel. humidity          
   integer(kind=jpim) :: q2m                                  = 281 !    2m specific humidity (q, kg/kg)         
   integer(kind=jpim) :: t                                    = 2   !    upper air temperature (K)
   integer(kind=jpim) :: td                                   = 59  !    upper air dew point (K) 
   integer(kind=jpim) :: t2m                                  = 39  !    2m temperature (K)           
   integer(kind=jpim) :: td2m                                 = 40  !    2m dew point (K)             
   integer(kind=jpim) :: ts                                   = 11  !    surface temperature (K)      
   integer(kind=jpim) :: ptend                                = 30  !    pressure tendency         
   integer(kind=jpim) :: w                                    = 60  !    past weather (w)          
   integer(kind=jpim) :: ww                                   = 61  !    present weather (ww)      
   integer(kind=jpim) :: vv                                   = 62  !    visibility                
   integer(kind=jpim) :: ch                                   = 63  !    type of high clouds (ch)  
   integer(kind=jpim) :: cm                                   = 64  !    type of middle clouds (cm)   
   integer(kind=jpim) :: cl                                   = 65  !    type of low clouds (cl)   
   integer(kind=jpim) :: nh                                   = 66  !    cloud base height (nh)    
   integer(kind=jpim) :: nn                                   = 67  !    low cloud amount (n)      
   integer(kind=jpim) :: hshs                                 = 68  !    additional cloud group height (hh)   
   integer(kind=jpim) :: c                                    = 69  !    additional cloud group type (c)   
   integer(kind=jpim) :: ns                                   = 70  !    additional cloud group amount (ns)   
   integer(kind=jpim) :: sdepth                               = 71  !    snow depth                
   integer(kind=jpim) :: e                                    = 72  !    state of ground (e)       
   integer(kind=jpim) :: tgtg                                 = 73  !    ground temperature (tgtg)   
   integer(kind=jpim) :: spsp1                                = 74  !    special phenomena (spsp)#1   
   integer(kind=jpim) :: spsp2                                = 75  !    special phenomena (spsp)#2   
   integer(kind=jpim) :: rs                                   = 76  !    ice code type (rs)        
   integer(kind=jpim) :: eses                                 = 77  !    ice thickness (eses)      
   integer(kind=jpim) :: is                                   = 78  !    ice (is)                  
   integer(kind=jpim) :: trtr                                 = 79  !    original time period of rain obs. (trtr)   
   integer(kind=jpim) :: rr                                   = 80  !    6hr rain (liquid part)    
   integer(kind=jpim) :: jj                                   = 81  !    max. temperature (jj)     
   integer(kind=jpim) :: vs                                   = 82  !    ship speed (vs)           
   integer(kind=jpim) :: ds                                   = 83  !    ship direction (ds)       
   integer(kind=jpim) :: hwhw                                 = 84  !    wave height               
   integer(kind=jpim) :: pwpw                                 = 85  !    wave period               
   integer(kind=jpim) :: dwdw                                 = 86  !    wave direction            
   integer(kind=jpim) :: gclg                                 = 87  !    general cloud group       
   integer(kind=jpim) :: rhlc                                 = 88  !    rel. humidity from low clouds      
   integer(kind=jpim) :: rhmc                                 = 89  !    rel. humidity from middle clouds   
   integer(kind=jpim) :: rhhc                                 = 90  !    rel. humidity from high clouds     
   integer(kind=jpim) :: n                                    = 91  !    total amount of clouds    
   integer(kind=jpim) :: sfall                                = 92  !    6hr snowfall (solid part of rain)   
   integer(kind=jpim) :: ps                                   = 110 !    surface pressure          
   integer(kind=jpim) :: dd                                   = 111 !    wind direction            
   integer(kind=jpim) :: ff                                   = 112 !    wind force                
   integer(kind=jpim) :: rawbt                                = 119 !    brightness temperature (K)
   integer(kind=jpim) :: rawra                                = 120 !    raw radiance              
   integer(kind=jpim) :: satcl                                = 121 !    cloud amount from satellite   
   integer(kind=jpim) :: scatss                               = 122 !    sigma 0   
   integer(kind=jpim) :: du                                   = 5   !    wind shear (du)   
   integer(kind=jpim) :: dv                                   = 6   !    wind shear (dv)   
   integer(kind=jpim) :: u10m                                 = 41  !    10m u component (m/s)
   integer(kind=jpim) :: v10m                                 = 42  !    10m v component (m/s)  
   integer(kind=jpim) :: rhlay                                = 19  !    layer rel. humidity   
   integer(kind=jpim) :: cllqw                                = 123 !    cloud liquid water   
   integer(kind=jpim) :: scatv                                = 124 !    ambiguous v component   
   integer(kind=jpim) :: scatu                                = 125 !    ambiguous u component   
   integer(kind=jpim) :: q                                    = 7   !    specific humidity (q)   
   integer(kind=jpim) :: scatwd                               = 126 !    ambiguous wind direction   
   integer(kind=jpim) :: scatws                               = 127 !    ambiguous wind speed       
   integer(kind=jpim) :: vsp                                  = 8   !    vertical speed   
   integer(kind=jpim) :: vt                                   = 56  !    virtual temperature   
   integer(kind=jpim) :: o3lay                                = 206 !    layer ozone   
   integer(kind=jpim) :: height                               = 156 !    height   
   integer(kind=jpim) :: onedvar                              = 215 !    1d-var model level (pseudo)-variable   
   integer(kind=jpim) :: w2                                   = 160 !    past weather 2 (used in synoptic maps)   
   integer(kind=jpim) :: cpt                                  = 130 !    characteristic of pressure tendency (used in synoptic maps)   
   integer(kind=jpim) :: tsts                                 = 12  !    sea water temperature (used in synoptic maps)   
   integer(kind=jpim) :: refl                                 = 192 !    radar reflectivity  
   integer(kind=jpim) :: apdss                                = 128 !    atmospheric path delay in satellite signal   
   integer(kind=jpim) :: bend_angle                           = 162 !    radio occultation bending angle   
   integer(kind=jpim) :: los                                  = 187 !    horizontal line-of-sight wind component     
   integer(kind=jpim) :: aerod                                = 174 !    aerosol optical depth at 0.55 microns  
   integer(kind=jpim) :: limb_radiance                        = 163 !    Limb Radiances  
   integer(kind=jpim) :: chem1                                = 181 !    chem1: no2/nox   
   integer(kind=jpim) :: chem2                                = 182 !    chem2: so2   
   integer(kind=jpim) :: chem3                                = 183 !    chem3: co   
   integer(kind=jpim) :: chem4                                = 184 !    chem4: hcho   
   integer(kind=jpim) :: chem5                                = 185 !    chem5: go3   
   integer(kind=jpim) :: chem6                                = 284 !    chem6: so2volc
   integer(kind=jpim) :: cod                                  = 175 !    cloud optical depth   
   integer(kind=jpim) :: rao                                  = 176 !    Ratio of fine mode to total aerosol optical depth at 0.55 microns   
   integer(kind=jpim) :: od                                   = 177 !    optical depth   
   integer(kind=jpim) :: rfltnc                               = 178 !    Aerosol reflectance multi-channel   
   integer(kind=jpim) :: nsoilm                               = 179 !    normalized soil moisture  (0-100%) 
   integer(kind=jpim) :: soilm                                = 180 !    soil moisture   
   integer(kind=jpim) :: flgt_phase                           = 201 !    phase of aircraft flight  
   integer(kind=jpim) :: height_assignment_method             = 211 !   Height assignment method   
   integer(kind=jpim) :: dopp                                 = 195 !    radar doppler wind   
   integer(kind=jpim) :: ghg1                                 = 186 !    ghg1: carbon dioxide   
   integer(kind=jpim) :: ghg2                                 = 188 !    ghg2: methane   
   integer(kind=jpim) :: ghg3                                 = 189 !    ghg3: nitrous oxide   
   integer(kind=jpim) :: bt_real                              = 190 !    brightness temperature real part  
   integer(kind=jpim) :: bt_imaginary                         = 191 !    brightness temperature imaginary part  
   integer(kind=jpim) :: prc                                  = 202 !    radar rain rate   
   integer(kind=jpim) :: lnprc                                = 203 !    log(radar rain rate mm/h + epsilon)   
   integer(kind=jpim) :: libksc                               = 222 !    lidar backscattering  
   integer(kind=jpim) :: ralt_swh                             = 220 !    significant wave height (m)           
   integer(kind=jpim) :: ralt_sws                             = 221 !    surface wind speed (m/s)              
   integer(kind=jpim) :: rawbt_clear                          = 193 !    brightness temperature for clear  (K) 
   integer(kind=jpim) :: rawbt_cloudy                         = 194 !    brightness temperature for cloudy (K)
   integer(kind=jpim) :: binary_snow_cover                    = 223 !    binary snow cover (0: no snow / 1: presence of snow)
   integer(kind=jpim) :: salinity                             = 224 !    ocean salinity (PSU)
   integer(kind=jpim) :: potential_temp                       = 225 !    potential temperature (Kelvin)
   integer(kind=jpim) :: humidity_mixing_ratio                = 226 !  humidity mixing ratio (kg/kg)
   integer(kind=jpim) :: airframe_icing                       = 227 !  airframe icing
   integer(kind=jpim) :: turbulence_index                     = 228 !  turbulence index
   integer(kind=jpim) :: lidar_aerosol_extinction             = 236 !  lidar aerosol extinction
   integer(kind=jpim) :: lidar_cloud_backscatter              = 237 !  lidar cloud backscatter
   integer(kind=jpim) :: lidar_cloud_extinction               = 238 !  lidar cloud extinction
   integer(kind=jpim) :: cloud_radar_reflectivity             = 239 !  cloud radar reflectivity
   integer(kind=jpim) :: pmsl                                 = 108 !  Mean sea-level pressure (Pa)
   integer(kind=jpim) :: liatbk                               = 280 !    lidar attenuated backscatter
end type
type(odb_varno)    :: varno

INTEGER(KIND=JPIM), PARAMETER :: NVNUMAX = 284

end module varno_module
