SUBROUTINE SUCDH1S

USE PARKIND1  ,ONLY : JPIM     ,JPRB,   JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMDPHY   , ONLY : NCSS
USE YOMCDH1S   ,ONLY : NLEVI   , &
     &  NDHVTLS,NDHFTLS,NDHVTSS,NDHFTSS, &
     &  NDHVTTS,NDHFTTS,NDHVTIS,NDHFTIS, &
     &  NDHVSSS,NDHFSSS,NDHVIIS,NDHFIIS, &
     &  NDHVWLS,NDHFWLS,NDHVRESS,NDHFRESS, &
     &  NDHVCO2S,NDHFCO2S, &
     &  NDHVBIOS,NDHFBIOS,NDHVVEGS,NDHFVEGS

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUCDH1S',0,ZHOOK_HANDLE)

NLEVI   = NCSS  ! nr of ice layers
! ! tiles
! NDHVTLS = 3  ! fraction, tskin, albedo
! NDHFTLS = 8  ! swdn, swup, lwdn, lwup, sens, lat, ground, evap


! tiles
NDHVTLS = 3  ! fraction, tskin, albedo
NDHFTLS = 9  ! swdn, swup, lwdn, lwup, sens, lat, ground, evap
! 1- fraction 
! 2 -  skin temp
! 3 - albedo
! 4 - swdown
! 5 - swup 
! 6 - lwdown
! 7 - lwup 
! 8 - sens
! 9 - latent
! 10 - ground
! 11 - evap
! 12 - skin conductivity # new 

! ! snow temp
! NDHVTSS = 6  ! heatcap, temp, energy/surface, thermal depth, density, alb
! NDHFTSS = 9  ! swdn, swup, lwdn, lwup, sens, lat, ground, phase change,
!              ! snow heat content change 
! ! snow mass
! NDHVSSS = 1  ! depth
! NDHFSSS = 4  ! ls sfall, cv sfall, evap, melt

NDHVTSS = 6 
NDHFTSS = 9 
!  1 - heat capacity
!  2 - temperature
!  3 - energy/surface
!  4 - thermal depth  (old scheme)   / number of snow active layers in multi-layer 
!  5 - density
!  6 - Albedo
!  7 - SWdown
!  8 - SWup
!  9 - lwdown
!  10 - lwup
!  11 - sensible 
!  12 - latent 
!  13 - ground flux 
!  14 - phase changes
!  15 - snow heat content change 

NDHVSSS = 2    
NDHFSSS = 6 
! 1 - snow mass
! 2 - large-scale snowfall
! 3 - convective snowfall 
! 4 - snow evaporatin 
! 5 - snow phase changes 
! 6- liquid water # new
! 7- intercepted rainfall # new
! 8 - snow runoff # new

! soil temp
NDHVTTS = 4  ! heat cap, temp, energy/surface, depth
NDHFTTS = 11 ! swdn, swup, lwdn, lwup, sens, lat, (ground)"soil-layers", 
              ! phase change, snow-soil flux, ground, soil heat content change
! sea ice
NDHVTIS = 4  ! heat cap, temp, energy/surface, depth
NDHFTIS = 9  ! swdn, swup, lwdn, lwup, sens, lat, ground, phase change, ice-soil
! interception
NDHVIIS = 1  ! depth
NDHFIIS = 3  ! ls net interc, cv net interc, evap
! soil water
NDHVWLS = 2  ! total soil water, frozen soil water
NDHFWLS = 8  ! ls tfallm, cv tfall, melt tfall, runoff, root extr, gw flux,
             ! bare ground evap incl mismatches, tile 7 condensation
NDHVRESS = 2 ! Canopy resistance, Aerodynamic resisitance
NDHFRESS = 0

! CO2 (per vegetation type, not normalized to grid square)
NDHVCO2S = 0
NDHFCO2S = 7 ! gross assimilation, dark respiration, net assimilation,
             ! soil+structural biomass respiration, ecosystem respiration, net CO2 flux, soil+str respiration without Q10 dependency
! Biomass (per vegetation type, not normalized to grid square)
NDHVBIOS = 0
NDHFBIOS = 5 ! leaf biomass, leaf biomass loss, leaf biomass gain, above ground
             ! structural biomass, below ground structural biomass
! Variables per vegetation type, not normalized to grid square	   
NDHVVEGS = 7 ! fraction, LAI (m2/m2), 
             ! Canopy conductance, Aerodynamic conductance, stress function F2,
             ! Ds (spec. humidity deficit at canopy level),Dmax 	
NDHFVEGS = 1 ! Latent heat flux (W/m2) 

IF (LHOOK) CALL DR_HOOK('SUCDH1S',1,ZHOOK_HANDLE)
 
RETURN
END SUBROUTINE SUCDH1S
