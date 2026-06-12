MODULE YOMGPD1S
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE
SAVE

!*    Grid point array for one timelevel physics fields

! GPD    : generic for one timelevel physics fields

! VFMASK  - catchment mask
! VFZ0F   - roughness length for momentum (constant) (NB: units M)
! VFALBF  - surface shortwave albedo
! VFALUVP - MODIS ALBEDO UV-VIS PARALLEL (DIRECT) RADIATION
! VFALUVD - MODIS ALBEDO UV-VIS DIFFUSE RADIATION
! VFALNIP - MODIS ALBEDO NEAR IR PARALLEL (DIRECT) RADIATION
! VFALNID - MODIS ALBEDO NEAR IR DIFFUSE RADIATION
! VFAL[UV|NI][I|V|G] - 6-component MODIS albedo coefficients
! VFITM   - land-sea mask
! VFGEO   - surface geopotential
! VFZ0H   - roughness length for heat (constant)  (NB: units M)
! VFCVL   - low vegetation cover
! VFCVH   - high vegetation cover
! VFCUR   - urban cover (PASSIVE)
! VFTVL   - low vegetation type
! VFTVH   - high vegetation type
! VFLAIL  - low vegetation lai
! VFLAIH  - high vegetation lai
! VFCO2TYP- photosynthesis pathway (C3/C4) (only applicable to low vegetation types)
! VFRSML  - low  vegetation minimum stomatal resistance
! VFRSMH  - high vegetation minimum stomatal resistance
! VFCOVVT - vegetation coverage (Cveg) of each vegetation type
! VFR0VT  - reference respiration at 25 

! VFVAROR - orographic variance

! VFSOTY  - soil type
! VFSDOR  - orographic standard deviation
! VFSST   - Sea Surface temperature
! VFCI    - sea ice fraction
! VFLDEPTH- LAKE DEPTH                                           !FLAKE
! VFCLAKE - LAKE FRACTION                                        !FLAKE
! VFZO    - vertical layer depth of ocean mixed layer model      !KPP
! VFHO    - vertical layer thickness of ocean mixed layer model  !KPP
! VFDO    - interface layer depth of ocean mixed layer model     !KPP
! VFOCDEPTH - ocean depth for ocean mixed layer model            !KPP
! VFADVT  - temperature advection correction term of KPP model   !KPP
! VFADVS  - salinity advection correction term of KPP model      !KPP
! VFTRI0  - trigonal matrix for diff. eq. of KPP model           !KPP
! VFTRI1  - trigonal matrix for diff. eq. of KPP model           !KPP
! VFHO_INV - 1/VFHO                                              !KPP 
! VFSWDK_SAVE - radiative coefficient for KPP model              !KPP
! VDLSP   - large scale precipitation             
! VDCP    - convective precipitation              
! VDSF    - snow fall                             
! VDSSHF  - surface sensible heat flux            
! VDSLHF  - surface latent heat flux              
! VDMSL   - mean sea level pressure               
! VDSSR   - surface solar radiation               
! VDSTR   - surface thermal radiation             
! VDEWSS  - U stress                 
! VDNSSS  - V stress
! VDE     - evaporation                            
! VDRO    - runoff                                 
! VDROS   - surface runoff                                 
! VDMLT   - snow melt                                 
! VDALB   - surface shortwave albedo as seen by model 
! VDIEWSS - instantaneous U stress                 
! VDINSSS - instantaneous V stress
! VDISSHF - instantaneous sensible heat flux
! VDIE    - instantaneous evaporation
! VDCSF   - convective snow fall
! VDLSSF  - large scale snow fall
! VDZ0F   - roughness length for momentum (varying over sea)  (NB: units M)
! VDZ0H   - roughness length for heat (varying over sea)  (NB: units M)
! VDSSRD  - surface solar radiation downwards               
! VDSTRD  - surface thermal radiation downwards            
! VDIEWSSTL- Instantaneous U stress for each tile                 
! VDINSSSTL- Intantaneous V stress for each tile
! VDISSHFTL- Instantaneous sensible heat flux for each tile
! VDIETL  - instantaneous evaporation for each tile
! VDTSKTL - skin temperature for each tile

! VDANDAYVT- daily net CO2 assimilation (sum) for each vegetation type
! VDANFMVT- daily maximum leaf assimilation (highest value) for each vegetation 
!           type
! VDRESPBSTR        - respiration of above ground structural biomass per VT
! VDRESPBSTR2       - respiration of below ground structural biomass per VT
! VDBIOMASS_LAST    - (active) leaf biomass of previous day per VT
! NB: value only after nitro_decline, not after laigain!!!!
! VDBLOSSVT- active biomass loss per VT
! VDBGAINVT- active biomass gain per VT

! VFPGLOB - global grid indices gaussian reduced

REAL(KIND=JPRB),ALLOCATABLE,TARGET:: GPD(:,:)

REAL(KIND=JPRB),POINTER:: VFMASK(:)
REAL(KIND=JPRB),POINTER:: VFZ0F(:)
REAL(KIND=JPRB),POINTER:: VFALBF(:)
REAL(KIND=JPRB),POINTER:: VFALUVP(:)
REAL(KIND=JPRB),POINTER:: VFALUVD(:)
REAL(KIND=JPRB),POINTER:: VFALNIP(:)
REAL(KIND=JPRB),POINTER:: VFALNID(:)
REAL(KIND=JPRB),POINTER:: VFALUVI(:)
REAL(KIND=JPRB),POINTER:: VFALUVV(:)
REAL(KIND=JPRB),POINTER:: VFALUVG(:)
REAL(KIND=JPRB),POINTER:: VFALNII(:)
REAL(KIND=JPRB),POINTER:: VFALNIV(:)
REAL(KIND=JPRB),POINTER:: VFALNIG(:)
REAL(KIND=JPRB),POINTER:: VFITM(:)
REAL(KIND=JPRB),POINTER:: VFGEO(:)
REAL(KIND=JPRB),POINTER:: VFZ0H(:)
REAL(KIND=JPRB),POINTER:: VFFWET(:)
REAL(KIND=JPRB),POINTER:: VFCVL(:)
REAL(KIND=JPRB),POINTER:: VFCVH(:)
REAL(KIND=JPRB),POINTER:: VFCUR(:)
REAL(KIND=JPRB),POINTER:: VFTVL(:)
REAL(KIND=JPRB),POINTER:: VFTVH(:)
REAL(KIND=JPRB),POINTER:: VFLAIL(:)
REAL(KIND=JPRB),POINTER:: VFLAIH(:)
REAL(KIND=JPRB),POINTER:: VFCO2TYP(:)
REAL(KIND=JPRB),POINTER:: VFRSML(:)
REAL(KIND=JPRB),POINTER:: VFRSMH(:)
REAL(KIND=JPRB),POINTER:: VFCVT(:,:)
REAL(KIND=JPRB),POINTER:: VFLAIVT(:,:)
REAL(KIND=JPRB),POINTER:: VFCOVVT(:,:)
REAL(KIND=JPRB),POINTER:: VFR0VT(:,:)


REAL(KIND=JPRB),POINTER:: VFSOTY(:)
REAL(KIND=JPRB),POINTER:: VFSDOR(:)
REAL(KIND=JPRB),POINTER:: VFSST(:)
REAL(KIND=JPRB),POINTER:: VFCI(:)
REAL(KIND=JPRB),POINTER:: VDLSP(:)
REAL(KIND=JPRB),POINTER:: VDCP(:)
REAL(KIND=JPRB),POINTER:: VDSF(:)
REAL(KIND=JPRB),POINTER:: VDSSHF(:)
REAL(KIND=JPRB),POINTER:: VDSLHF(:)
REAL(KIND=JPRB),POINTER:: VDSSR(:)
REAL(KIND=JPRB),POINTER:: VDSTR(:)
REAL(KIND=JPRB),POINTER:: VDEWSS(:)
REAL(KIND=JPRB),POINTER:: VDNSSS(:)
REAL(KIND=JPRB),POINTER:: VDE(:)
REAL(KIND=JPRB),POINTER:: VDRO(:)
REAL(KIND=JPRB),POINTER:: VDROS(:)
REAL(KIND=JPRB),POINTER:: VDMLT(:)
REAL(KIND=JPRB),POINTER:: VDALB(:)
REAL(KIND=JPRB),POINTER:: VDIEWSS(:)
REAL(KIND=JPRB),POINTER:: VDINSSS(:)
REAL(KIND=JPRB),POINTER:: VDISSHF(:)
REAL(KIND=JPRB),POINTER:: VDIE(:)
REAL(KIND=JPRB),POINTER:: VDCSF(:)
REAL(KIND=JPRB),POINTER:: VDLSSF(:)
REAL(KIND=JPRB),POINTER:: VDZ0F(:)
REAL(KIND=JPRB),POINTER:: VDZ0H(:)
REAL(KIND=JPRB),POINTER:: VDSSRD(:)
REAL(KIND=JPRB),POINTER:: VDSTRD(:)
REAL(KIND=JPRB),POINTER:: VDIEWSSTL(:,:)
REAL(KIND=JPRB),POINTER:: VDINSSSTL(:,:)
REAL(KIND=JPRB),POINTER:: VDISSHFTL(:,:)
REAL(KIND=JPRB),POINTER:: VDIETL(:,:)
REAL(KIND=JPRB),POINTER:: VDTSKTL(:,:)

REAL(KIND=JPRB),POINTER:: VFLDEPTH(:)       !FLAKE
REAL(KIND=JPRB),POINTER:: VFCLAKE(:)        !FLAKE 
REAL(KIND=JPRB),POINTER:: VFCLAKEF(:)        !FLAKE + Flood

REAL(KIND=JPRB),POINTER:: VFZO(:,:)         !KPP
REAL(KIND=JPRB),POINTER:: VFHO(:,:)         !KPP
REAL(KIND=JPRB),POINTER:: VFDO(:,:)         !KPP
REAL(KIND=JPRB),POINTER:: VFOCDEPTH(:)      !KPP
REAL(KIND=JPRB),POINTER:: VFADVT(:,:)       !KPP
REAL(KIND=JPRB),POINTER:: VFADVS(:,:)       !KPP
REAL(KIND=JPRB),POINTER:: VFTRI0(:,:)       !KPP
REAL(KIND=JPRB),POINTER:: VFTRI1(:,:)       !KPP
REAL(KIND=JPRB),POINTER:: VFHO_INV(:,:)     !KPP
REAL(KIND=JPRB),POINTER:: VFSWDK_SAVE(:,:)  !KPP

REAL(KIND=JPRB),POINTER:: VDANDAYVT(:,:)
REAL(KIND=JPRB),POINTER:: VDANFMVT(:,:)
REAL(KIND=JPRB),POINTER:: VDRESPBSTR(:,:)
REAL(KIND=JPRB),POINTER:: VDRESPBSTR2(:,:)
REAL(KIND=JPRB),POINTER:: VDBIOMASS_LAST(:,:)
REAL(KIND=JPRB),POINTER:: VDBLOSSVT(:,:)
REAL(KIND=JPRB),POINTER:: VDBGAINVT(:,:)

REAL(KIND=JPRB),POINTER:: VFPGLOB(:) 


END MODULE YOMGPD1S
