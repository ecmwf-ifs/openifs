INTERFACE
SUBROUTINE SURFEXCDRIVER    (YDSURF, CDCONF &
 & , KIDIA, KFDIA, KLON, KLEVS, KTILES, KVTYPES, KDIAG, KSTEP &
 & , KLEVSN, KLEVI, KDHVTLS, KDHFTLS, KDHVTSS, KDHFTSS &
 & , KDHVTTS, KDHFTTS, KDHVTIS, KDHFTIS, K_VMASS &
 & , KDHVCO2S,KDHFCO2S,KDHVVEGS,KDHFVEGS &
 & , PTSTEP,PTSTEPF &
! input data, non-tiled
 & , KTVL, KCO2TYP, KTVH, PCVL, PCVH, PCUR &
 & , PLAIL, PLAIH, PFWET, PLAT &
 & , PSNM , PRSN &
 & , PMU0 , PCARDI &
 & , PUMLEV, PVMLEV, PTMLEV, PQMLEV, PCMLEV,PAPHMS, PGEOMLEV, PCPTGZLEV &
 & , PSST, PTSKM1M, PCHAR, PCHARHQ, PSSRFL, PSLRFL, PEMIS, PTICE, PTSN &
 & , PHLICE,PTLICE,PTLWML & 
 & , PTHKICE,PSNTICE &
 & , PWLMX, PUCURR, PVCURR, PI10FGCV &
! input data, soil
 & , PTSAM1M, PWSAM1M, KSOTY &
! input data, tiled
 & , PFRTI, PALBTI &
! updated data, tiled
 & , PUSTRTI, PVSTRTI, PAHFSTI, PEVAPTI, PTSKTI &
 & , PANDAYVT,PANFMVT &
! updated data, non-tiled
 & , PZ0M, PZ0H &
! output data, tiled
 & , PSSRFLTI, PQSTI, PDQSTI, PCPTSTI, PCFHTI, PCFQTI, PCSATTI, PCAIRTI &
 & , PCPTSTIU, PCSATTIU, PCAIRTIU, PRAQTI, PTSRF, PLAMSK &
! output data, non-tiled
 & , PKHLEV, PKCLEV, PCFMLEV, PKMFL, PKHFL, PKQFL, PEVAPSNW &
 & , PZ0MW, PZ0HW, PZ0QW, PBLENDPP, PCPTSPP, PQSAPP, PBUOMPP, PZDLPP &
! output data, non-tiled CO2
 & , PAN,PAG,PRD,PRSOIL_STR,PRECO,PCO2FLUX,PCH4FLUX&
 & , PWETB, PWETL, PWETLU, PWETH, PWETHS &
! o.utput data, diagnostics
 & , PDHTLS, PDHTSS, PDHTTS, PDHTIS &
 & , PDHVEGS, PEXDIAG, PDHCO2S &
 & , PRPLRG &
 & , LSICOUP &
 & )

USE PARKIND1, ONLY : JPIM, JPRB
USE ISO_C_BINDING


! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!------------------------------------------------------------------------

!  PURPOSE:
!    Routine SURFEXCDRIVER controls the ensemble of routines that prepare
!    the surface exchange coefficients and associated surface quantities
!    needed for the solution of the vertical diffusion equations. 

!  SURFEXCDRIVER is called by VDFMAIN

!  METHOD:
!    This routine is only a shell needed by the surface library
!    externalisation.

!  AUTHOR:
!    P. Viterbo       ECMWF May 2005   

!  REVISION HISTORY:
!    E. Dutra/G. Balsamo May 2008 Add lake tile
!    S. Boussetta/G.Balsamo May 2009 Add lai
!    S. Boussetta/G.Balsamo May 2010 Add CTESSEL
!    N. Semane  04-10-2012  Include small planet PRPLRG
!    Linus Magnusson     10-09-28 Sea-ice

!  INTERFACE: 

!    Characters (In):
!      CDCONF   :    IFS Configuration

!    Integers (In):
!      KIDIA    :    Begin point in arrays
!      KFDIA    :    End point in arrays
!      KLON     :    Length of arrays
!      KLEVS    :    Number of soil layers
!      KTILES   :    Number of tiles
!      KVTYPES  :    Number of biomes for carbon
!      KDIAG    :    Number of diagnostic parameters
!      KSTEP    :    Time step index
!      KLEVSN   :    Number of snow layers (diagnostics) 
!      KLEVI    :    Number of sea ice layers (diagnostics)
!      KDHVTLS  :    Number of variables for individual tiles
!      KDHFTLS  :    Number of fluxes for individual tiles
!      KDHVTSS  :    Number of variables for snow energy budget
!      KDHFTSS  :    Number of fluxes for snow energy budget
!      KDHVTTS  :    Number of variables for soil energy budget
!      KDHFTTS  :    Number of fluxes for soil energy budget
!      KDHVTIS  :    Number of variables for sea ice energy budget
!      KDHFTIS  :    Number of fluxes for sea ice energy budget
!      K_VMASS  :    Controls the use of vector functions in the IBM scientific
!                     library. Set K_VMASS=0 to use standard functions
!      KTVL     :    Dominant low vegetation type
!      KCO2TYP  :    Photosynthetic pathway for low vegetation type (C3/C4)
!      KTVH     :    Dominant high vegetation type
!      KSOTY    :    SOIL TYPE                                        (1-7)

!    *KDHVCO2S*     Number of variables for CO2
!    *KDHFCO2S*     Number of fluxes for CO2
!    *KDHVVEGS*     Number of variables for vegetation
!    *KDHFVEGS*     Number of fluxes for vegetation


!    Reals (In):
!      PTSTEP   :    Timestep
!      PCVL     :    Low vegetation fraction
!      PCVH     :    High vegetation fraction
!      PCUR     :    Urban cover                                     (0-1)
!      PLAIL    :    Low vegetation LAI
!      PLAIH    :    High vegetation LAI

!     PSNM      :       SNOW MASS (per unit area)                      kg/m**2
!     PRSN      :      SNOW DENSITY                                   kg/m**3

!      PMU0          : COS SOLAR angle
!      PCARDI        : CONCENTRATION ATMOSPHERIC CO2
!      PRPLRG    :   ! GRAVITY SMALL PLANET FACTOR 

!    Reals with tile index (In): 
!      PFRTI    :    TILE FRACTIONS                                   (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!      PALBTI   :    Tile albedo                                      (0-1)

!    Reals independent of tiles (In):
!      PUMLEV   :    X-VELOCITY COMPONENT, lowest atmospheric level   m/s
!      PVMLEV   :    Y-VELOCITY COMPONENT, lowest atmospheric level   m/s
!      PTMLEV   :    TEMPERATURE,   lowest atmospheric level          K
!      PQMLEV   :    SPECIFIC HUMIDITY                                kg/kg
!      PCMLEV   :    ATMOSPHERIC CO2 (TRACER)                         kg/kg
!      PAPHMS   :    Surface pressure                                 Pa
!      PGEOMLEV :    Geopotential, lowest atmospehric level           m2/s2
!      PCPTGZLEV:    Geopotential, lowest atmospehric level           J/kg
!      PSST     :    (OPEN) SEA SURFACE TEMPERATURE                   K
!      PTSKM1M  :    SKIN TEMPERATURE                                 K
!      PCHAR    :    CHARNOCK PARAMETER                               -
!      PCHARHQ  :    CHARNOCK PARAMETER FOR HEAT AND MOISTURE         -
!      PSSRFL   :    NET SHORTWAVE RADIATION FLUX AT SURFACE          W/m2
!      PSLRFL   :    NET LONGWAVE RADIATION FLUX AT SURFACE           W/m2
!      PEMIS    :    MODEL SURFACE LONGWAVE EMISSIVITY
!      PTSAM1M  :    SURFACE TEMPERATURE                              K
!      PWSAM1M  :    SOIL MOISTURE ALL LAYERS                         m**3/m**3
!      PTICE    :    Ice temperature, top slab                        K
!      PTSN     :    Snow temperature                                 K
!      PHLICE   :    Lake ice thickness                               m
!      PTLICE   :    Lake ice temperature                             K
!      PTLWML   :    Lake mean water temperature                      K
!      PTHKICE  :    Sea-ice thickness                                m
!      PSNTICE  :    Snow thickness on sea-ice                        m
!      PWLMX    :    Maximum interception layer capacity              kg/m**2
!      PUCURR   :    u component of ocean surface current             m/s
!      PVCURR   :    v component of ocean surface current             m/s
!      PI10FGCV :    gust velocity from deep convection               m/s

!    Reals with tile index (In/Out):
!      PUSTRTI  :    SURFACE U-STRESS                                 N/m2 
!      PVSTRTI  :    SURFACE V-STRESS                                 N/m2
!      PAHFSTI  :    SURFACE SENSIBLE HEAT FLUX                       W/m2
!      PEVAPTI  :    SURFACE MOISTURE FLUX                            KG/m2/s
!      PTSKTI   :    SKIN TEMPERATURE                                 K

!    	UPDATED PARAMETERS FOR VEGETATION TYPES (REAL):
!      PANDAYVT :    DAILY NET CO2 ASSIMILATION OVER CANOPY           kg_CO2/m**2
!      PANFMVT  :    MAXIMUM LEAF ASSIMILATION                        kg_CO2/kg_Air m/s

!    Reals independent of tiles (In/Out):
!      PZ0M     :    AERODYNAMIC ROUGHNESS LENGTH                     m
!      PZ0H     :    ROUGHNESS LENGTH FOR HEAT                        m

!    Reals with tile index (Out):
!      PSSRFLTI :    Tiled NET SHORTWAVE RADIATION FLUX AT SURFACE    W/m2
!      PQSTI    :    Tiled SATURATION Q AT SURFACE                    kg/kg
!      PDQSTI   :    Tiled DERIVATIVE OF SATURATION Q-CURVE           kg/kg/K
!      PCPTSTI  :    Tiled DRY STATIC ENERGY AT SURFACE               J/kg
!      PCFHTI   :    Tiled EXCHANGE COEFFICIENT AT THE SURFACE        ????
!      PCFQTI   :    Tiled EXCHANGE COEFFICIENT AT THE SURFACE        ????
!      PCSATTI  :    MULTIPLICATION FACTOR FOR QS AT SURFACE          -
!                      FOR SURFACE FLUX COMPUTATION
!      PCAIRTI  :    MULTIPLICATION FACTOR FOR Q AT  LOWEST MODEL     - 
!                      LEVEL FOR SURFACE FLUX COMPUTATION
!      PCPTSTIU :    AS PCPTSTI FOR UNSTRESSED LOW VEGETATION         J/kg
!      PCSATTIU :    AS PCSATTI FOR UNSTRESSED LOW VEGETATION         -
!      PCAIRTIU :    AS PCAIRTI FOR UNSTRESSED LOW VEGETATION         - 
!      PRAQTI   :    Aerodynamic resistance                           s/m 
!      PTSRF    :    Tiled surface temperature for each tile 
!                    Boundary condition in surfseb                    K
!      PLAMSK   :    Tiled skin layer conductivity                    W m-2 K-1

!    Reals independent of tiles (Out):
!      PKHLEV   :    SURFACE LAYER: CH*U                              m/s
!      PKCLEV   :    SURFACE LAYER: Cc*U, tracer transfer             m/s
!      PCFMLEV  :    PROP. TO EXCH. COEFF. FOR MOMENTUM               ????
!                     (C-STAR IN DOC.) (SURFACE LAYER ONLY)
!      PKMFL    :    Kinematic momentum flux                          ????
!      PKHFL    :    Kinematic heat flux                              ????
!      PKQFL    :    Kinematic moisture flux                          ????
!      PEVAPSNW :    Evaporation from snow under forest               kgm-2s-1
!      PZ0MW    :    Roughness length for momentum, WMO station       m
!      PZ0HW    :    Roughness length for heat, WMO station           m
!      PZ0QW    :    Roughness length for moisture, WMO station       m
!      PBLENDPP :    Blending weight for 10 m wind postprocessing     m
!      PCPTSPP  :    Cp*Ts for post-processing of weather parameters  J/kg
!      PQSAPP   :    Apparent surface humidity for post-processing    kg/kg
!                     of weather parameters
!      PBUOMPP  :    Buoyancy flux, for post-processing of gustiness  ???? 
!      PZDLPP   :    z/L for post-processing of weather parameters    -
!      PDHTLS   :    Diagnostic array for tiles (see module yomcdh)
!                      (Wm-2 for energy fluxes, kg/(m2s) for water fluxes)
!      PDHTSS   :    Diagnostic array for snow T (see module yomcdh)
!                      (Wm-2 for fluxes)
!      PDHTTS   :    Diagnostic array for soil T (see module yomcdh)
!                      (Wm-2 for fluxes)
!      PDHTIS   :    Diagnostic array for ice T (see module yomcdh)
!                      (Wm-2 for fluxes)

!     *PDHVEGS*      Diagnostic array for vegetation (see module yomcdh) 
!     *PEXDIAG*      Diagnostic array for optional pp of canopy resistances
!     *PDHCO2S*      Diagnostic array for CO2 (see module yomcdh)

!     EXTERNALS.
!     ----------

!     ** SURFEXCDRIVER_CTL CALLS SUCCESSIVELY:
!         *VUPDZ0*
!         *VSURF*
!         *CO2* 
!         *VEXCS*
!         *VEVAP*
!         *VSFLX*

!  DOCUMENTATION:
!    See Physics Volume of IFS documentation

!------------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

TYPE(C_PTR),                  INTENT(IN)  :: YDSURF

CHARACTER(LEN=1)  ,INTENT(IN)    :: CDCONF 

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM),INTENT(IN)    :: KVTYPES
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIAG
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVSN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVI 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTLS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTLS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTTS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTTS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_VMASS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTEPF
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCO2TYP(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTY(:)
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVCO2S
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFCO2S
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVVEGS
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFVEGS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCUR(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAIL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAIH(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFWET(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNM(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSN(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCARDI
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCMLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOMLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICE(:)   
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLICE(:)  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLWML(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTHKICE(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNTICE(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHAR(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHARHQ(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTICE(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSN(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWLMX(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI10FGCV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSAM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSAM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSTRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSTRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PANDAYVT(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PANFMVT(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0M(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0H(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSSRFLTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQSTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDQSTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCPTSTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCFHTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCFQTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCSATTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAIRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCPTSTIU(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCSATTIU(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAIRTIU(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAQTI(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTSRF(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLAMSK(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKHLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKCLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCFMLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKMFL(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKHFL(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKQFL(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEVAPSNW(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ0MW(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ0HW(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ0QW(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBLENDPP(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCPTSPP(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQSAPP(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBUOMPP(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZDLPP(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAN(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAG(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRD(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRSOIL_STR(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRECO(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCO2FLUX(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCH4FLUX(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWETB(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWETL(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWETLU(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWETH(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWETHS(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHVEGS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXDIAG(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHCO2S(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTLS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTSS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTTS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTIS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPLRG

LOGICAL           ,INTENT(IN)    :: LSICOUP

!------------------------------------------------------------------------

END SUBROUTINE SURFEXCDRIVER
END INTERFACE
