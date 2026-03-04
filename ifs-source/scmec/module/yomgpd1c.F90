! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

MODULE YOMGPD1C

!*    Grid point array for one time-level physics fields

! VEXTRA  - generic for extra 3-d fields
! VEXTR2  - generic for extra 2-d fields

! VUSTRTI - E-W (instantaneous) surface stress for each tile
! VVSTRTI - N-S (instantaneous) surface stress for each tile
! VAHFSTI - (instantaneous) surface sensible heat flux for each tile
! VEVAPTI - (instantaneous) evaporation for each tile
! VTSKTI  - skin temperature for each tile

! VFZ0F   - g * surface roughness length for momentum (callpar in - varying over sea)
! VFLZ0H  - logarithm of roughness length for heat    (callpar in)
! VDZ0F   - g * surface roughness length for momentum (callpar out - current)
! VDLZ0H  - logarithm of roughness length for heat    (callpar out - current)
! VFALBF  - surface shortwave albedo                  (callpar in)
! VDALB   - surface shortwave albedo as seen by model (callpar out)
! VFEMISF - surface longwave emmisivity               (callpar in)
! VFEMIS  - surface longwave emmisivity               (callpar out)
! VFITM   - land-sea mask
! VFVAROG - directional orographic variances                     
! VCVL    - low vegetation cover
! VCVH    - high vegetation cover
! VTVL    - low vegetation dominant type
! VTVH    - high vegetation dominant type
! VCI     - sea-ice fraction
! VSST    - (open) sea surface temperature

! VNTOP   - index of convective cloud top         
! VNBAS   - index of convective cloud base         
! VNACPR  - averaged convective precipitation rate 
! VNTYPE  - convection type (0: none, 1: deep, 2: shallow, 3: mid-level)

! VREMTU  - upward longwave emmisivity             

! VDLSP   - large scale precipitation (sfc)           
! VDCP    - convective precipitation  (sfc)       
! VDCSF   - convective snow fall      (sfc)        
! VDLSSF  - large scale snow fall     (sfc)        
! VDBLD   - boundary layer dissipitaion           
! VDSSHF  - surface sensible heat flux            
! VDSLHF  - surface latent heat flux              
! VDMSL   - mean sea level pressure               
! VDTCC   - total cloud cover                     
! VD10U   - 10 metre u                            
! VD10V   - 10 metre v                            
! VD2T    - 2 metre temperature                   
! VD2D    - 2 metre dewpoint temperature          
! VDSSR   - surface solar radiation               
! VDSTR   - surface thermal radiation             
! VDTSR   - top solar radiation                    
! VDTTR   - top thermal radiation                  
! VDEWSS  - U-sterss                               
! VDNSSS  - V-stress  
! VDE     - evaporation                            
! VBLH    - boundary layer height
! VISUND  - sunshine duration (fraction of time step)
! VDCCC   - convective cloud cover                 
! VDLCC   - low cloud cover                        
! VDMCC   - medium cloud cover                     
! VDHCC   - high cloud cover                       
! VDLGWS  - lat. comp. of gravity wave stress      
! VDMGWS  - mer. comp. of gravity wave stress      
! VDGWD   - gravity wave dissipation               
! VDMX2T  - max temp. at 2 m since prev. p.p.      
! VDMN2T  - min temp. at 2 m since prev. p.p.      
! VDRO    - runoff                                 
! VDIEWSS - instantaneous U stress                 
! VDINSSS - instantaneous V stress                 
! VDISSHF - instantaneous surface sensible heat flux       
! VDIE    - instantaneous surface evaporation              
! VLAIL   - leaf area index low vegitation
! VLAIH   - leaf area index high vegitation
! VSOTY   - soil type               
! VALUVP  - UV parallel albedo
! VALUVD  - UV diffuse albedo      
! VALNIP  - nearIR parallel albedo        
! VALNID  - nearIR diffuse albedo  



! VI10FG  - wind gust at 10 m for current time level (m/s)
! VILSPF  - large scale precipitation fraction (0-1)

! VRCLFR  - cloud fractional cover for radiation
! VRCLLW  - cloud liquid water     for radiation
! VRCLIW  - cloud ice water        for radiation
! VQSAT   - specific humidity at saturation

! NOGWFIELD - Increments of horizontal wind from NOOGW scheme
! SPFIELDS  - Prognostic quantities for SP-CRM scheme

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE PARDIM1C , ONLY : NTILESMX ,JPCEXTR  ,JPFLEV

IMPLICIT NONE

SAVE

INTEGER(KIND=JPIM), PARAMETER :: KLON = 1

REAL(KIND=JPRB), pointer :: VEXTRA(:,:,:)
REAL(KIND=JPRB), pointer :: VEXTR2(:,:)

REAL(KIND=JPRB), target :: VUSTRTI(1,NTILESMX)
REAL(KIND=JPRB), target :: VVSTRTI(1,NTILESMX)
REAL(KIND=JPRB), target :: VAHFSTI(1,NTILESMX)
REAL(KIND=JPRB), target :: VEVAPTI(1,NTILESMX)
REAL(KIND=JPRB), target :: VTSKTI (1,NTILESMX)


INTEGER         :: INTYPE(KLON)
INTEGER         :: IPBLTYPE(KLON)

REAL(KIND=JPRB), target  :: VFEMIS(KLON) 
REAL(KIND=JPRB), target  :: VISUND(KLON) !-> PISUND -> PRAD%PISUND

REAL(KIND=JPRB) :: VDLSP(KLON)
REAL(KIND=JPRB) :: VDCP(KLON)
REAL(KIND=JPRB) :: VDBLD(KLON)
REAL(KIND=JPRB) :: VDSSHF(KLON)
REAL(KIND=JPRB) :: VDSLHF(KLON)
REAL(KIND=JPRB) :: VDMSL(KLON)

REAL(KIND=JPRB) :: VDSSR(KLON)
REAL(KIND=JPRB) :: VDSTR(KLON)
REAL(KIND=JPRB) :: VDTSR(KLON)
REAL(KIND=JPRB) :: VDTTR(KLON)
REAL(KIND=JPRB) :: VDEWSS(KLON)
REAL(KIND=JPRB) :: VDNSSS(KLON)
REAL(KIND=JPRB) :: VDE(KLON)


REAL(KIND=JPRB) :: VDLGWS(KLON)
REAL(KIND=JPRB) :: VDMGWS(KLON)
REAL(KIND=JPRB) :: VDGWD(KLON)
REAL(KIND=JPRB) :: VDMX2T(KLON)
REAL(KIND=JPRB) :: VDMN2T(KLON)
REAL(KIND=JPRB) :: VDRO(KLON)
REAL(KIND=JPRB) :: VDIEWSS(KLON)
REAL(KIND=JPRB) :: VDINSSS(KLON)
REAL(KIND=JPRB) :: VDISSHF(KLON)
REAL(KIND=JPRB) :: VDIE(KLON)
REAL(KIND=JPRB) :: VDCSF(KLON)
REAL(KIND=JPRB) :: VDLSSF(KLON)

REAL(KIND=JPRB) :: VI10FG(KLON) ! -> PI10FG -> PDIAG%PI10FG 
REAL(KIND=JPRB) :: VILSPF(KLON,JPFLEV) ! -> PILSPF -> PDIAG%PCOVPTOT(JL,JK)

REAL(KIND=JPRB) :: VQSAT(JPFLEV)  ! Beware not to used it as target!!!

REAL(KIND=JPRB) :: NOGWFIELD(JPFLEV,2)

! Quantities for SETTLS extrapolation scheme in SL advection
!  NOTE: In this case also GFL fields should be treated by this
!        technique. The stored LS horizontal advection should
!        be considered as RHS.
REAL(KIND=JPRB) :: RHST9(JPFLEV), RHSQ9(JPFLEV)
REAL(KIND=JPRB) :: RHSU9(JPFLEV), RHSV9(JPFLEV)
!  Quantity used by SETTLST scheme for trajectory research
REAL(KIND=JPRB) :: WRL9(JPFLEV)
!  Quantity used for SL physics on GFL (LSLPHY=true)
REAL(KIND=JPRB), allocatable :: GFLSLP(:,:,:)

!     ------------------------------------------------------------------

END MODULE YOMGPD1C
