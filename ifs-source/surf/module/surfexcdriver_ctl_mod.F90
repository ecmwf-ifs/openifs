MODULE SURFEXCDRIVER_CTL_MOD
CONTAINS
SUBROUTINE SURFEXCDRIVER_CTL(CDCONF &
 & , KIDIA, KFDIA, KLON, KLEVS, KTILES, KVTYPES, KDIAG, KSTEP &
 & , KLEVSN, KLEVI, KDHVTLS, KDHFTLS, KDHVTSS, KDHFTSS &
 & , KDHVTTS, KDHFTTS, KDHVTIS, KDHFTIS, K_VMASS &
 & , PTSTEP,PTSTEPF &
! input data, non-tiled
 & , KTVL, KCO2TYP, KTVH, PCVL, PCVH, PCUR &
 & , PLAIL, PLAIH, PFWET, PLAT &
 & , PSNM , PRSN &
 & , PMU0 , PCARDI &
 & , PUMLEV, PVMLEV, PTMLEV, PQMLEV, PCMLEV, PAPHMS, PGEOMLEV, PCPTGZLEV &
 & , PSST, PTSKM1M, PCHAR, PCHARHQ, PSSRFL, PSLRFL, PEMIS, PTICE, PTSN &
 & , PHLICE,PTLICE,PTLWML &   
 & , PTHKICE, PSNTICE &
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
 & , PCPTSTIU,PCSATTIU, PCAIRTIU,PRAQTI,PTSRF,PLAMSK &
! output data, non-tiled
 & , PKHLEV, PKCLEV, PCFMLEV, PKMFL, PKHFL, PKQFL, PEVAPSNW &
 & , PZ0MW, PZ0HW, PZ0QW, PBLENDPP, PCPTSPP, PQSAPP, PBUOMPP, PZDLPP &
! output data, non-tiled CO2
 & , PAN,PAG,PRD,PRSOIL_STR,PRECO,PCO2FLUX,PCH4FLUX&
! output canopy resistance  
 & , PWETB, PWETL, PWETLU, PWETH, PWETHS &
! output data, diagnostics
 & , PDHTLS, PDHTSS, PDHTTS, PDHTIS &
 & , PDHVEGS, PEXDIAG, PDHCO2S &
 & , PRPLRG &
! LIM switch
 & , LSICOUP &
 & , YDCST, YDEXC, YDVEG, YDAGS, YDAGF, YDSOIL, YDFLAKE, YDURB & 
 & )

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_EXC   , ONLY : TEXC
USE YOS_CST   , ONLY : TCST
USE YOS_VEG   , ONLY : TVEG
USE YOS_AGS   , ONLY : TAGS
USE YOS_AGF   , ONLY : TAGF
USE YOS_SOIL  , ONLY : TSOIL
USE YOS_FLAKE , ONLY : TFLAKE
USE YOS_URB   , ONLY : TURB
USE YOS_THF   , ONLY : R4LES, R5LES, R2ES, R4IES, R3LES, R3IES, R5IES
USE VUPDZ0_MOD
USE VSURF_MOD
USE VEXCS_MOD
USE VEVAP_MOD
USE SURFSEB_CTL_MOD
USE SRFCOTWO_MOD
USE VSFLX_MOD
USE VLAMSK_MOD 

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
!    A. Beljaars      10/12/2005  TOFD
!    A. Beljaars      17/05/2007  clean-up of roughness length initialization
!    G. Balsamo       15/11/2007  Use aggregated Z0M for drag and dominant low
!                                 for post-processing of 2m T/TD.
!    E. Dutra/G. Balsamo  01/05/2008  lake tile
!    A. Beljaars/M.Koehler 14/01/2009 Surfcae flux bugfix for stability
!    S. Boussetta/G.Balsamo May 2009 Add lai
!    G. Balsamo       15/09/2009  Fix lake tile temperature initialization
!    S. Boussetta/G.Balsamo May 2010 Add CTESSEL
!    I. Sandu        May 2010 Reduce blending height for post-processing
!    N.Semane+P.Bechtold 04-10-2012 Add PRPLRG factor for small planet
!    Linus Magnusson      10-09-28 Sea-ice
!    I. Sandu    24-02-2014  Lambda skin values by vegetation type instead of tile
!    A. Beljaars      26/02/2014  Compute unstressed evaporation
!    A. Agusti-Panareda 26/10/2017 Call CTESSEL at step 0
!    G. Balsamo+I. Sandu 22/03/2018 - correction 2T for wet skin tile and lakes
!    A. Agusti-Panareda Nov 2020  couple atm CO2 tracer (new input) with photosynthesis  
!    V.Bastrikov,F.Maignan,P.Peylin,A.Agusti-Panareda/S.Boussetta Feb 2021 Add Farquhar photosynthesis model

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
!      KDIAG    :    Number of diagnostic fields
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
!      KCO2TYP :    Type of photosynthetic pathway for low vegetation type (c3/c4)
!      KTVH     :    Dominant high vegetation type
!      KSOTY    :    SOIL TYPE                                        (1-6)

!    Reals (In):
!      PTSTEP   :    Timestep 
!      PTSTEPF  :    Full actual model timestep (PTSTEP can be a sub-step)

!    Reals with tile index (In): 
!      PFRTI    :    TILE FRACTIONS                                   (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!            9 : LAKE                  10 : URBAN
!      PALBTI   :    Tile albedo                                      (0-1)

!    Reals independent of tiles (In):
!      PCVL     :    LOW VEGETATION COVER (CLIMATOLOGICAL)
!      PCVH     :    HIGH VEGETATION COVER (CLIMATOLOGICAL)
!      PCUR     :    URBAN COVER                                     (0-1)
!      PLAIL    :    LOW VEGETATION LAI                              m2/m2
!      PLAIH    :    HIGH VEGETATION LAI                             m2/m2

!     PSNM      :       SNOW MASS (per unit area)                      kg/m**2
!     PRSN      :      SNOW DENSITY                                   kg/m**3

!      PMU0     : COS SOLAR
!      PCARDI   : CONCENTRATION ATMOSPHERIC CO2


!      PUMLEV   :    X-VELOCITY COMPONENT, lowest atmospheric level   m/s
!      PVMLEV   :    Y-VELOCITY COMPONENT, lowest atmospheric level   m/s
!      PTMLEV   :    TEMPERATURE,   lowest atmospheric level          K
!      PQMLEV   :    SPECIFIC HUMIDITY                                kg/kg
!      PCMLEV   :    ATMOSPHERIC CO2                                  kg/kg
!      PAPHMS   :    Surface pressure                                 Pa
!      PGEOMLEV :    Geopotential, lowest atmospehric level           m2/s2
!      PCPTGZLEV:    Geopotential, lowest atmospehric level           J/kg
!      PSST     :    (OPEN) SEA SURFACE TEMPERATURE                   K
!      PUSTRTI  :    X-STRESS                                         N/m2
!      PVSTRTI  :    Y-STRESS                                         N/m2
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
!      PTHKICE  :    Sea-ice thickness
!      PSNTICE  :    Snow thickness on sea-ice
!      PWLMX    :    Maximum interception layer capacity              kg/m**2
!      PUCURR   :    u component of ocean surface current             m/s
!      PVCURR   :    v component of ocean surface current             m/s
!      PI10FGCV :    gust velocity from deep convcetion               m/s

!    Reals with tile index (In/Out):
!      PUSTRTI  :    SURFACE U-STRESS                                 N/m2 
!      PVSTRTI  :    SURFACE V-STRESS                                 N/m2
!      PAHFSTI  :    SURFACE SENSIBLE HEAT FLUX                       W/m2
!      PEVAPTI  :    SURFACE MOISTURE FLUX                            KG/m2/s
!      PTSKTI   :    SKIN TEMPERATURE                                 K

!    UPDATED PARAMETERS FOR VEGETATION TYPES (REAL):
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
!      PCFHTI   :    Tiled transfer coefficient heat     Rho*Ch*U     kgm-2s-1
!      PCFQTI   :    Tiled transfer coefficient moisture Rho*Cq*U     kgm-2s-1
!      PCSATTI  :    MULTIPLICATION FACTOR FOR QS AT SURFACE          -
!                      FOR SURFACE FLUX COMPUTATION
!      PCAIRTI  :    MULTIPLICATION FACTOR FOR Q AT  LOWEST MODEL     - 
!                      LEVEL FOR SURFACE FLUX COMPUTATION
!      PCPTSTIU :    AS PCPTSTI FOR UNSTRESSED LOW VEGETATION         J/kg
!      PCSATTIU :    AS PCSATTI FOR UNSTRESSED LOW VEGETATION         -
!      PCAIRTIU :    AS PCAIRTI FOR UNSTRESSED LOW VEGETATION         -
!      PRAQTI   :    AERODYNAMIC RESISTANCE, tiled                    s/m
!      PTSRF    :    Tiled surface temperature for each tile 
!                    Boundary condition in surfseb                    K
!      PLAMSK   :    Tiled skin layer conductivity                    W m-2 K-1

!    Reals independent of tiles (Out):
!      PKHLEV   :    Scaled transfer Coeff. heat         Ch*U         m/s
!      PKCLEV   :    Scaled transfer Coeff. for tracers  Ch*U         m/s
!      PCFMLEV  :    Transfer coeff. for momentum    Rho*Cm*U         kgm-2s-1
!      PKMFL    :    Kinematic momentum flux                          m2s-2
!      PKHFL    :    Kinematic heat flux                              Kms-1
!      PKQFL    :    Kinematic moisture flux                          kgkg-1ms-1
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
!     *PEXDIAG*      Extra diagnostics for pp of canopy stresses
!     *PDHCO2S*      Diagnostic array for CO2 (see module yomcdh)
!       PWETB        Surface resistance bare soild 
!       PWETL        Canopy resistance of low vegetation  
!       PWETLU       Canopy resistance of low vegetation , unstressed 
!       PWETH        Canopy resistance of high vegetation  
!       PWETHS       Canopy resistance of high vegetation snow cover  



!     EXTERNALS.
!     ----------

!     ** SURFEXCDRIVER_CTL CALLS SUCCESSIVELY:
!         *VUPDZ0*
!         *VSURF*
!         *CO2* 
!         *VEXCS*
!         *VEVAP*
!         *SURFSEB_CTL*

!         *VSFLX* (may be not needed as fluxes are computed within this routine)

!  DOCUMENTATION:
!    See Physics Volume of IFS documentation

!------------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

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

!CO2
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAN(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAG(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRD(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRSOIL_STR(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRECO(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCO2FLUX(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCH4FLUX(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHVEGS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHCO2S(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXDIAG(:,:) 

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTLS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTSS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTTS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTIS(:,:,:)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPLRG
LOGICAL           ,INTENT(IN)    :: LSICOUP  

!canopy / bare soild resistances 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWETB(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWETL(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWETLU(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWETH(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWETHS(:)

TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TEXC)        ,INTENT(IN)    :: YDEXC
TYPE(TVEG)        ,INTENT(IN)    :: YDVEG
TYPE(TAGS)        ,INTENT(IN)    :: YDAGS
TYPE(TAGF)        ,INTENT(IN)    :: YDAGF
TYPE(TSOIL )      ,INTENT(IN)    :: YDSOIL
TYPE(TFLAKE)      ,INTENT(IN)    :: YDFLAKE
TYPE(TURB)        ,INTENT(IN)    :: YDURB

! Local variables
!ZLIQ is passed to compute soil moisture scaling factor in Reco (CO2 routine) CTESSEL
REAL(KIND=JPRB) :: PWLIQ(KLON,KLEVS,KTILES)

INTEGER(KIND=JPIM) :: IFRMAX(KLON),IFRLMAX(KLON)

REAL(KIND=JPRB) :: ZZ0MTI(KLON,KTILES) , ZZ0HTI(KLON,KTILES) ,&
                 & ZZ0QTI(KLON,KTILES) , ZBUOMTI(KLON,KTILES),&
                 & ZZDLTI(KLON,KTILES) , &
                 & ZQSATI(KLON,KTILES) , ZCFMTI(KLON,KTILES) ,&
                 & ZKMFLTI(KLON,KTILES), ZKHFLTI(KLON,KTILES),&
                 & ZKQFLTI(KLON,KTILES), ZZQSATI(KLON,KTILES),&
                 & ZJS(KLON,KTILES)    , ZJQ(KLON,KTILES)    ,&
                 & ZSSK(KLON,KTILES)   , ZTSK(KLON,KTILES)   ,&
                 & ZSSH(KLON,KTILES)   , ZSLH(KLON,KTILES)   ,&
                 & ZSTR(KLON,KTILES)   , ZG0(KLON,KTILES)    

REAL(KIND=JPRB) :: ZANTI(KLON,KTILES)  , ZAGTI(KLON,KTILES),&
                 & ZRDTI(KLON,KTILES)  , ZLAMSK(KLON,KTILES)

REAL(KIND=JPRB) :: ZLAIVT(KLON,2+1),& ! number vegetation tiles
                   ZCVT(KLON,2+1)
INTEGER(KIND=JPIM) :: KVEG(KLON,KTILES)
INTEGER(KIND=JPIM)  :: KVTTL(KTILES) !link tile number/veg type



REAL(KIND=JPRB) :: ZFRMAX(KLON)   , ZFRLMAX(KLON)  , ZALB(KLON)     , &
                 & ZSRFD(KLON)    , ZWETL(KLON)    , ZWETH(KLON)    , &
                 & ZWETHS(KLON)   , ZWETB(KLON)    , ZKHLEV(KLON)   , &
                 & ZTSA(KLON)     , ZCSNW(KLON)    , ZSSRFL1(KLON)  , &
                 & ZCBLENDM(KLON) , ZCBLENDH(KLON) , ZSL(KLON)      , &
                 & ZQL(KLON)      , ZASL(KLON)     , ZBSL(KLON)     , &
                 & ZAQL(KLON)     , ZBQL(KLON)     , ZRHO(KLON)     , &
                 & ZWETLU(KLON)   , ZDSN(KLON)     , ZKCLEV(KLON) 



INTEGER(KIND=JPIM) :: JL, JTILE, IITT,JK
LOGICAL :: LLINIT

REAL(KIND=JPRB) :: ZDUA, ZZCDN, ZQSSN, ZCOR, ZRG, &
                 & ZZ0MWMO, ZBLENDWMO, ZBLENDZ0, ZCOEF1
REAL(KIND=JPRB) :: ZRVTRSR

LOGICAL         :: LLAND, LLSICE, LLHISSR(KLON)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fcsttre.h"

!*         1.     Set up of general quantities
!                 ----------------------------

IF (LHOOK) CALL DR_HOOK('SURFEXCDRIVER_CTL_MOD:SURFEXCDRIVER_CTL',0,ZHOOK_HANDLE)
ASSOCIATE(LEOCWA=>YDEXC%LEOCWA, LEOCCO=>YDEXC%LEOCCO, REPDU2=>YDEXC%REPDU2, &
 & RKAP=>YDEXC%RKAP, RZ0ICE=>YDEXC%RZ0ICE, &
 & RALFMAXSN=>YDSOIL%RALFMAXSN, &
 & LEFLAKE=>YDFLAKE%LEFLAKE, RH_ICE_MIN_FLK=>YDFLAKE%RH_ICE_MIN_FLK, &
 & RTT=>YDCST%RTT, RCPD=>YDCST%RCPD, RD=>YDCST%RD, RSIGMA=>YDCST%RSIGMA, &
 & RG=>YDCST%RG, RETV=>YDCST%RETV, &
 & RVZ0M=>YDVEG%RVZ0M, LECTESSEL=>YDVEG%LECTESSEL, RVTRSR=>YDVEG%RVTRSR)

ZRG         = 1.0_JPRB/RG        !     -"-
DO JL=KIDIA,KFDIA
  ZRHO(JL)=PAPHMS(JL)/( RD*PTMLEV(JL)*(1.0_JPRB+RETV*PQMLEV(JL)) )
ENDDO


! initialisation of carbon related array : 

ZANTI(KIDIA:KFDIA,:)=0._JPRB
ZAGTI(KIDIA:KFDIA,:)=0._JPRB
ZRDTI(KIDIA:KFDIA,:)=0._JPRB

!mapping betwen tiles and vegetation TYPE
KVEG(KIDIA:KFDIA,:)=0_JPIM
KVEG(KIDIA:KFDIA,4)=KTVL(KIDIA:KFDIA) !type low veg
KVEG(KIDIA:KFDIA,6)=KTVH(KIDIA:KFDIA) !type high veg KVEG
KVEG(KIDIA:KFDIA,7)=KTVH(KIDIA:KFDIA) !shaded snow same type high veg !why ?

!initialisation LAI
ZLAIVT(KIDIA:KFDIA,:)=0.0_JPRB
ZLAIVT(KIDIA:KFDIA,1)=PLAIL(KIDIA:KFDIA)
ZLAIVT(KIDIA:KFDIA,2)=PLAIH(KIDIA:KFDIA)

!intialisation vegetation type fraction
ZCVT(KIDIA:KFDIA,:)=0.0_JPRB
ZCVT(KIDIA:KFDIA,1)=PCVL(KIDIA:KFDIA)
ZCVT(KIDIA:KFDIA,2)=PCVH(KIDIA:KFDIA)

! Total snow depth (m) 
ZDSN(KIDIA:KFDIA) = SUM( PSNM(KIDIA:KFDIA,:) / PRSN(KIDIA:KFDIA,:),DIM=2) 

!*         1.1  ESTIMATE SURF.FL. FOR STEP 0
!*              (ASSUME NEUTRAL STRATIFICATION)
IF ( KSTEP == 0 ) THEN
  DO JTILE=1,KTILES
    DO JL=KIDIA,KFDIA
      PTSKTI(JL,JTILE)=PTSKM1M(JL)
    ENDDO
  ENDDO
  IF ((.NOT. LEOCWA) .AND. (.NOT. LEOCCO)) THEN
    DO JL=KIDIA,KFDIA
      PTSKTI(JL,1)=PSST(JL)
    ENDDO
  ENDIF
ENDIF

!*         1.2  UPDATE Z0
CALL VUPDZ0(KIDIA,KFDIA,KLON,KTILES,KSTEP,CDCONF,&
   & KTVL,KTVH,PCVL,PCVH,PCUR,PUMLEV,PVMLEV,&
   & PTMLEV,PQMLEV,PAPHMS,PGEOMLEV,ZDSN,&
   & PUSTRTI,PVSTRTI,PAHFSTI,PEVAPTI,&
   & PHLICE,& 
   & PTSKTI,PCHAR,PCHARHQ,PFRTI,PUCURR,PVCURR,&
   & YDCST,YDEXC,YDVEG,YDFLAKE,YDURB,&
   & ZZ0MTI,ZZ0HTI,ZZ0QTI,ZBUOMTI,ZZDLTI,PRAQTI)

!*         1.3  FIND DOMINANT SURFACE TYPE and DOMINANT LOW
!*              parameters for postprocessing
ZFRMAX(KIDIA:KFDIA)=PFRTI(KIDIA:KFDIA,1)
ZFRLMAX(KIDIA:KFDIA)=PFRTI(KIDIA:KFDIA,1)
IFRMAX(KIDIA:KFDIA)=1
IFRLMAX(KIDIA:KFDIA)=1
DO JTILE=2,KTILES
  DO JL=KIDIA,KFDIA
    IF (PFRTI(JL,JTILE)  >  ZFRMAX(JL)) THEN
      ZFRMAX(JL)=PFRTI(JL,JTILE)
      IFRMAX(JL)=JTILE
    ENDIF
    IF (PFRTI(JL,JTILE)  >  ZFRLMAX(JL) .AND. &
      JTILE.NE.6 .AND. JTILE.NE.7) THEN
      ZFRLMAX(JL)=PFRTI(JL,JTILE)
      IFRLMAX(JL)=JTILE
      IF (JTILE.EQ.3.OR.JTILE.EQ.9) THEN
!* for tiles wet-skin or lakes attribute if present
!* low-vegetation (4) if present or bare soil (8) 
         IF (PFRTI(JL,8).GT.0.0_JPRB) IFRLMAX(JL)=8
         IF (PFRTI(JL,4).GT.0.0_JPRB) IFRLMAX(JL)=4
!        IF (PFRTI(JL,8).GT.0.0_JPRB) IFRLMAX(JL)=8
!        IF (PFRTI(JL,4).GT.PFRTI(JL,8)) IFRLMAX(JL)=4
!        IF (PFRTI(JL,5).GT.PFRTI(JL,4) .AND. PFRTI(JL,5).GT.PFRTI(JL,8) ) IFRLMAX(JL)=5
      ENDIF
    ENDIF
  ENDDO
ENDDO

!*         Use tile average (log) Z0 for M and H
ZBLENDZ0=10._JPRB/PRPLRG
ZCBLENDM(KIDIA:KFDIA)=PFRTI(KIDIA:KFDIA,1)&
           &/(LOG(ZBLENDZ0/ZZ0MTI(KIDIA:KFDIA,1)))**2
ZCBLENDH(KIDIA:KFDIA)=PFRTI(KIDIA:KFDIA,1)&
           &/(LOG(ZBLENDZ0/ZZ0HTI(KIDIA:KFDIA,1)))**2
DO JTILE=2,KTILES
  DO JL=KIDIA,KFDIA
    ZCBLENDM(JL)=ZCBLENDM(JL)&
           &+PFRTI(JL,JTILE)/(LOG(ZBLENDZ0/ZZ0MTI(JL,JTILE)))**2
    ZCBLENDH(JL)=ZCBLENDH(JL)&
           &+PFRTI(JL,JTILE)/(LOG(ZBLENDZ0/ZZ0HTI(JL,JTILE)))**2
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  PZ0M(JL)=ZBLENDZ0*EXP(-1._JPRB/SQRT(ZCBLENDM(JL)))
  PZ0H(JL)=ZBLENDZ0*EXP(-1._JPRB/SQRT(ZCBLENDH(JL)))
ENDDO

!*         Put dominant tile Z0M on all tiles
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    ZZ0MTI(JL,JTILE)=PZ0M(JL)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         2.     SURFACE BOUNDARY CONDITIONS FOR T AND Q
!                 ---------------------------------------

!    2.1 Albedo
ZALB(KIDIA:KFDIA)=&
 &  PFRTI(KIDIA:KFDIA,1)*PALBTI(KIDIA:KFDIA,1)&
 & +PFRTI(KIDIA:KFDIA,2)*PALBTI(KIDIA:KFDIA,2)&
 & +PFRTI(KIDIA:KFDIA,3)*PALBTI(KIDIA:KFDIA,3)&
 & +PFRTI(KIDIA:KFDIA,4)*PALBTI(KIDIA:KFDIA,4)&
 & +PFRTI(KIDIA:KFDIA,5)*PALBTI(KIDIA:KFDIA,5)&
 & +PFRTI(KIDIA:KFDIA,6)*PALBTI(KIDIA:KFDIA,6)&
 & +PFRTI(KIDIA:KFDIA,7)*PALBTI(KIDIA:KFDIA,7)&
 & +PFRTI(KIDIA:KFDIA,8)*PALBTI(KIDIA:KFDIA,8) 

IF (LEFLAKE) THEN
  ZALB(KIDIA:KFDIA)=ZALB(KIDIA:KFDIA)&
 & +PFRTI(KIDIA:KFDIA,9)*PALBTI(KIDIA:KFDIA,9)   
ENDIF

IF (KTILES .GT. 9) THEN
  ZALB(KIDIA:KFDIA)=ZALB(KIDIA:KFDIA)+PFRTI(KIDIA:KFDIA,10)*PALBTI(KIDIA:KFDIA,10)
ENDIF     

ZSSRFL1(KIDIA:KFDIA)=0._JPRB

LLHISSR(KIDIA:KFDIA)=.FALSE.
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
! Disaggregate solar flux but limit to 700 W/m2 (due to inconsistency
!  with albedo)
    PSSRFLTI(JL,JTILE)=((1.0_JPRB-PALBTI(JL,JTILE))/&
   & (1.0_JPRB-ZALB(JL)))*PSSRFL(JL)
    IF (PSSRFLTI(JL,JTILE) > 700._JPRB) THEN
      LLHISSR(JL)=.TRUE.
      PSSRFLTI(JL,JTILE)=700._JPRB
    ENDIF

! Compute averaged net solar flux after limiting to 700 W/m2
    ZSSRFL1(JL)=ZSSRFL1(JL)+PFRTI(JL,JTILE)*PSSRFLTI(JL,JTILE) 
  ENDDO
ENDDO

DO JTILE=1,KTILES
 DO JL=KIDIA,KFDIA  
   DO JK=1,KLEVS
PWLIQ(JL,JK,JTILE)=0.0_JPRB
   ENDDO
 ENDDO
ENDDO

DO JTILE=1,KTILES


IF     (JTILE  ==  4) THEN
KVTTL(JTILE)=1 !type low veg
ELSEIF (JTILE  ==  6) THEN
KVTTL(JTILE)=2 !type high veg 
ELSEIF (JTILE  ==  7) THEN
KVTTL(JTILE)=2 !type high veg 
ELSE
KVTTL(JTILE)=3
ENDIF

ENDDO

DO JTILE=1,KTILES

  IF (JTILE  ==  1.OR. JTILE  ==  2.OR. JTILE  ==  3.OR. JTILE  ==  5 .OR. JTILE == 9) THEN 
    DO JL=KIDIA,KFDIA  
      DO JK=1,KLEVS
        PWLIQ(JL,JK,JTILE)=0.0_JPRB
      ENDDO
    ENDDO
  ENDIF 
  DO JL=KIDIA,KFDIA
    IF (LLHISSR(JL)) THEN
      PSSRFLTI(JL,JTILE)=PSSRFLTI(JL,JTILE)*PSSRFL(JL)/ZSSRFL1(JL)
    ENDIF
    ZSRFD(JL)=PSSRFLTI(JL,JTILE)/(1.0_JPRB-PALBTI(JL,JTILE))  
  ENDDO
  CALL VSURF(KIDIA,KFDIA,KLON,KTILES,KLEVS,JTILE,&
   & KTVL,KCO2TYP,KTVH,&
   & KVTTL(JTILE),KVEG,&
   & ZLAIVT(:,KVTTL(JTILE)),& 
   & PMU0,PCO2FLUX,&
   & PFRTI, PLAIL, PLAIH,&
   & PTMLEV  ,PQMLEV  , PCMLEV, PAPHMS,&
   & PTSKTI(:,JTILE),PWSAM1M,PTSAM1M,KSOTY,&
   & ZSRFD,PRAQTI(:,JTILE),ZQSATI(:,JTILE),&
   & PQSTI(:,JTILE)  ,PDQSTI(:,JTILE)  ,&
   & ZWETB ,PCPTSTI(:,JTILE) ,ZWETL, ZWETLU, ZWETH, ZWETHS ,&
   & PEVAPTI(:,JTILE) ,&
   & ZANTI(:,JTILE),ZAGTI(:,JTILE),ZRDTI(:,JTILE) ,&
   & PWLIQ(:,:,JTILE), PDHVEGS , PEXDIAG,&
   & YDCST, YDVEG, YDEXC, YDAGS, YDAGF, YDSOIL, YDFLAKE, YDURB)
ENDDO

! cp arraysi, later repalce local copy  
PWETB(KIDIA:KFDIA)=ZWETB(KIDIA:KFDIA) 
PWETL(KIDIA:KFDIA)=ZWETL(KIDIA:KFDIA) 
PWETLU(KIDIA:KFDIA)=ZWETLU(KIDIA:KFDIA) 
PWETH(KIDIA:KFDIA)=ZWETH(KIDIA:KFDIA) 
PWETHS(KIDIA:KFDIA)=ZWETHS(KIDIA:KFDIA) 


    IF (LECTESSEL) THEN  ! usage of CTESSEL

      IF ( KSTEP == 0 ) THEN
!orig        PCO2FLUX(KIDIA:KFDIA)=0._JPRB
!orig        PCH4FLUX(KIDIA:KFDIA)=0._JPRB
!orig        PAG(KIDIA:KFDIA)=0._JPRB
!orig        PRD(KIDIA:KFDIA)=0._JPRB
!orig        PAN(KIDIA:KFDIA)=0._JPRB
!orig        PRSOIL_STR(KIDIA:KFDIA)=0._JPRB
!orig        PRECO(KIDIA:KFDIA)=0._JPRB

        PANDAYVT(KIDIA:KFDIA,:)=0._JPRB
        PANFMVT(KIDIA:KFDIA,:)=0._JPRB
        IF (SIZE(PDHVEGS) > 0) PDHVEGS(KIDIA:KFDIA,:,:)=0._JPRB
        IF (SIZE(PDHCO2S) > 0) PDHCO2S(KIDIA:KFDIA,:,:)=0._JPRB
      ENDIF

!orig      ELSE
        CALL SRFCOTWO(KIDIA,KFDIA,KLON,KLEVS,KTILES,&
   & KVEG,KVTTL,KCO2TYP,KSOTY,& 
   & PTSTEP ,&
   & PTMLEV,PQMLEV,PAPHMS,&
   & ZCVT,PFRTI,ZLAIVT,PWLIQ,&
!for respiration the soil temperature of the second layer is used 
   & PWSAM1M(:,:),PTSAM1M(:,:),&
     & ZDSN,PFWET,PLAT,&
   & ZANTI,ZAGTI,ZRDTI,&
   & YDCST,YDVEG,YDSOIL,YDAGS,&
   & PANDAYVT,PANFMVT,&
   & PAG,PRD,PAN,PRSOIL_STR,&
   & PRECO,PCO2FLUX,PCH4FLUX,&
   & PDHCO2S)

!orig      ENDIF
    ENDIF


!*          3.3x Surface temperature and Skin conductivity 
!               -------------------------
IF (LEOCWA .OR. LEOCCO) THEN
  PTSRF(KIDIA:KFDIA,1)=PTSKTI(KIDIA:KFDIA,1)
ELSE
  PTSRF(KIDIA:KFDIA,1)=PSST(KIDIA:KFDIA)
ENDIF
PTSRF(KIDIA:KFDIA,2)=PTICE(KIDIA:KFDIA)
PTSRF(KIDIA:KFDIA,3)=PTSAM1M(KIDIA:KFDIA,1)
PTSRF(KIDIA:KFDIA,4)=PTSAM1M(KIDIA:KFDIA,1)
PTSRF(KIDIA:KFDIA,5)=PTSN(KIDIA:KFDIA,1)
PTSRF(KIDIA:KFDIA,6)=PTSAM1M(KIDIA:KFDIA,1)
PTSRF(KIDIA:KFDIA,7)=PTSN(KIDIA:KFDIA,1)
PTSRF(KIDIA:KFDIA,8)=PTSAM1M(KIDIA:KFDIA,1)
IF (KTILES>=9) THEN
  WHERE(PHLICE(KIDIA:KFDIA) > RH_ICE_MIN_FLK )
    PTSRF(KIDIA:KFDIA,9) = PTLICE(KIDIA:KFDIA)
  ELSEWHERE
    PTSRF(KIDIA:KFDIA,9) = PTLWML(KIDIA:KFDIA)
  ENDWHERE
ENDIF
IF (KTILES>=10) THEN
PTSRF(KIDIA:KFDIA,10)=PTSAM1M(KIDIA:KFDIA,1)
ENDIF
!! PLAMSK with Full time-step (this is wrong but should is kept for now 
!! to allow bit identical results 
CALL VLAMSK(KIDIA,KFDIA,KLON,KTILES,KTVL,KTVH,&
            PTSTEPF,PTSKTI,PTSRF,&
            PSNM,PRSN,PSNTICE,&
            PWSAM1M,KSOTY,&
            YDCST,YDVEG,YDSOIL,YDURB,LSICOUP,&
            PLAMSK)

!! PLAMSK with actual VDF* time-step 
CALL VLAMSK(KIDIA,KFDIA,KLON,KTILES,KTVL,KTVH,&
            PTSTEP,PTSKTI,PTSRF,&
            PSNM,PRSN,PSNTICE,&
            PWSAM1M,KSOTY,&
            YDCST,YDVEG,YDSOIL,YDURB,LSICOUP,&
            ZLAMSK)

! DDH diagnostics

IF (SIZE(PDHTLS) > 0 .AND. SIZE(PDHTSS) > 0 .AND. SIZE(PDHTTS) > 0 .AND. &
  & SIZE(PDHTIS) > 0 .AND. SIZE(PDHVEGS) > 0 .AND. SIZE(PDHCO2S) > 0) &
  & CALL COMPUTE_DDH

!*         3.     EXCHANGE COEFFICIENTS
!                 ---------------------

!*         3.1  SURFACE EXCHANGE COEFFICIENTS


LLINIT=KSTEP == 0
IF (KSTEP <= 3) THEN
  IITT=3
ELSE
  IITT=1
ENDIF

DO JTILE=1,KTILES

  CALL VEXCS(KIDIA,KFDIA,KLON,IITT,K_VMASS,LLINIT,&
   & PUMLEV,PVMLEV,PTMLEV,PQMLEV,PAPHMS,PGEOMLEV,PCPTGZLEV,&
   & PCPTSTI(:,JTILE),ZQSATI(:,JTILE),&
   & ZZ0MTI(:,JTILE),ZZ0HTI(:,JTILE),&
   & ZZ0QTI(:,JTILE),ZZDLTI(:,JTILE),ZBUOMTI(:,JTILE),&
   & PUCURR,PVCURR,PI10FGCV,&
   & YDCST,YDEXC,&
   & ZCFMTI(:,JTILE),PCFHTI(:,JTILE),&
   & PCFQTI(:,JTILE),ZKHLEV,ZKCLEV)

  DO JL=KIDIA,KFDIA
    IF (JTILE == IFRMAX(JL)) THEN 
      PKHLEV(JL)=ZKHLEV(JL)
      PKCLEV(JL)=ZKCLEV(JL)
    ENDIF
  ENDDO
ENDDO

!*         3.2  EQUIVALENT EVAPOTRANSPIRATION EFFICIENCY COEFFICIENT

DO JTILE=1,KTILES
  IF     (JTILE == 1) THEN
    ZTSA(KIDIA:KFDIA)=PSST(KIDIA:KFDIA)
  ELSEIF (JTILE == 2) THEN
    ZTSA(KIDIA:KFDIA)=PTICE(KIDIA:KFDIA)
  ELSEIF (JTILE == 5 .OR. JTILE == 7) THEN
    ZTSA(KIDIA:KFDIA)=PTSN(KIDIA:KFDIA,1)

  ELSEIF (JTILE == 9) THEN
    DO JL=KIDIA,KFDIA
      IF(PHLICE(JL) > RH_ICE_MIN_FLK ) THEN
        ZTSA(JL)=PTLICE(JL)
      ELSE
        ZTSA(JL)=PTLWML(JL)
      ENDIF
    ENDDO
  ELSE
    ZTSA(KIDIA:KFDIA)=PTSAM1M(KIDIA:KFDIA,1)
  ENDIF
  CALL VEVAP(KIDIA,KFDIA,KLON,PTSTEP,JTILE,&
   & PWLMX ,PTMLEV  ,PQMLEV  ,PAPHMS,PTSKTI(:,JTILE),ZTSA,&
   & PQSTI(:,JTILE),PCFQTI(:,JTILE),ZWETB,ZWETL,ZWETLU,ZWETH,ZWETHS,&
   & YDCST,YDVEG,YDURB,&
   & PCPTSTI(:,JTILE),PCSATTI(:,JTILE),PCAIRTI(:,JTILE),&
   & PCPTSTIU(:,JTILE),PCSATTIU(:,JTILE),PCAIRTIU(:,JTILE),ZCSNW,PDHVEGS)
ENDDO

!          COMPUTE SNOW EVAPORATION FROM BELOW TREES i.e. TILE 7

! Note the use of qsat(Tsnow), rather than tile 7 skin. Skin T7 is a
! canopy temperature, definitely not what is desirable. Skin T5 can go
! up (and down ..) freely, not really what we want. The use of
! qsat (Tsnow) is tantamount to neglecting the skin effect there.

DO JL=KIDIA,KFDIA
  IF (PFRTI(JL,7) > 0.0_JPRB) THEN
    !Only compute snow evap. when there is a snow cover under trees 
    !(just to skip funct computation)
    ZQSSN=FOEEW(PTSN(JL,1))/PAPHMS(JL)
    ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSSN)
    ZQSSN=ZQSSN*ZCOR
    PEVAPSNW(JL)=PCFQTI(JL,7)*ZCSNW(JL)*(PQMLEV(JL)-ZQSSN)
  ELSE
    PEVAPSNW(JL)=0.0_JPRB
  END IF

ENDDO

!*         3.3  COMPUTE SURFACE FLUXES FOR TILES

DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
!   use previous times tep fluxes for heat and moisture
    ZKHFLTI(JL,JTILE)=PAHFSTI(JL,JTILE)/(ZRHO(JL)*RCPD)
    ZKQFLTI(JL,JTILE)=PEVAPTI(JL,JTILE)/ZRHO(JL)

    ZKMFLTI(JL,JTILE)=ZCFMTI(JL,JTILE)*SQRT((PUMLEV(JL)-PUCURR(JL))**2&
   & +(PVMLEV(JL)-PVCURR(JL))**2)/ZRHO(JL)
  ENDDO
ENDDO


!*         3.3a   PREPARE ARRAY'S FOR CALL TO SURFACE ENERGY
!                 BALANCE ROUTINE
IF (KSTEP == 0) THEN 

  

  ZASL(KIDIA:KFDIA)=0.0_JPRB
  ZBSL(KIDIA:KFDIA)=PCPTGZLEV(KIDIA:KFDIA)
  ZAQL(KIDIA:KFDIA)=0.0_JPRB
  ZBQL(KIDIA:KFDIA)=PQMLEV(KIDIA:KFDIA)


!*         3.3b   CALL TO SURFACE ENERGY BALANCE ROUTINE

  CALL SURFSEB_CTL(KIDIA,KFDIA,KLON,KTILES,KTVL,KTVH,&
   & PTSTEP,PCPTSTI,PTSKTI,PQSTI,&
   & PDQSTI,PCFHTI,PCFQTI,&
   & PCAIRTI,PCSATTI,&
   & PSSRFLTI,PFRTI,PTSRF,ZLAMSK,&
   & PSNM(:,1),PRSN(:,1),PHLICE,&
   & PSLRFL,PTSKM1M,PEMIS,&
   & ZASL,ZBSL,ZAQL,ZBQL,&
   & PTHKICE, PSNTICE, &
   & YDCST,YDEXC,YDVEG,YDFLAKE,YDSOIL,YDURB,&
   !out
   & ZJS,ZJQ,ZSSK,ZTSK,&
   & ZSSH,ZSLH,ZSTR,ZG0,&
   & ZSL,ZQL,&
   & LSICOUP)  

  DO JTILE=1,KTILES
    DO JL=KIDIA,KFDIA
      ZKHFLTI(JL,JTILE)=ZSSH(JL,JTILE)/(ZRHO(JL)*RCPD)
      ZKQFLTI(JL,JTILE)=ZJQ(JL,JTILE)/ZRHO(JL)
    ENDDO
  ENDDO

ENDIF


!          ADD SNOW EVAPORATION FROM BELOW TREES i.e. TILE 7

ZKQFLTI(KIDIA:KFDIA,7)=ZKQFLTI(KIDIA:KFDIA,7)+&
 & ZCSNW(KIDIA:KFDIA)*ZKQFLTI(KIDIA:KFDIA,5)  

!*         3.4  COMPUTE SURFACE FLUXES, WEIGHTED AVERAGE OVER TILES

PKMFL(KIDIA:KFDIA)=0.0_JPRB
PKHFL(KIDIA:KFDIA)=0.0_JPRB
PKQFL(KIDIA:KFDIA)=0.0_JPRB
PCFMLEV(KIDIA:KFDIA)=0.0_JPRB
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    PKMFL(JL)=PKMFL(JL)+PFRTI(JL,JTILE)*ZKMFLTI(JL,JTILE)
    PKHFL(JL)=PKHFL(JL)+PFRTI(JL,JTILE)*ZKHFLTI(JL,JTILE)
    PKQFL(JL)=PKQFL(JL)+PFRTI(JL,JTILE)*ZKQFLTI(JL,JTILE)
    PCFMLEV(JL)=PCFMLEV(JL)+PFRTI(JL,JTILE)*ZCFMTI(JL,JTILE)
  ENDDO
ENDDO

!*         4.  Preparation for "POST-PROCESSING" of surface weather parameters

!          POST-PROCESSING WITH MINIMUM OF LOCAL AND EFFECTIVE
!          SURFACE ROUGHNESS LENGTH. THE LOCAL ONES ARE FOR
!          WMO-TYPE WIND STATIONS I.E. OPEN TERRAIN WITH GRASS

ZBLENDWMO=40._JPRB/PRPLRG
ZZ0MWMO=0.03_JPRB/PRPLRG
DO JL=KIDIA,KFDIA
  IF (PZ0M(JL)  >  ZZ0MWMO) THEN
    PZ0MW(JL)=ZZ0MWMO
    PBLENDPP(JL)=ZBLENDWMO
  ELSE
    PZ0MW(JL)=PZ0M(JL)
    PBLENDPP(JL)=PGEOMLEV(JL)*ZRG
  ENDIF
ENDDO

DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    ZZQSATI(JL,JTILE)=PQMLEV(JL)*(1.0_JPRB-PCAIRTI(JL,JTILE))&
     & +PCSATTI(JL,JTILE)*PQSTI(JL,JTILE)  
    ZZQSATI(JL,JTILE)=MAX(1.0E-12_JPRB,ZZQSATI(JL,JTILE))
  ENDDO
ENDDO

!          ROUGHNESS LENGTH FOR HEAT and MOISTURE ARE TAKEN
!          FROM THE DOMINANT LOW-VEG. TYPE
DO JL=KIDIA,KFDIA
  JTILE=IFRLMAX(JL)
  PZ0HW(JL)=ZZ0HTI(JL,JTILE)
  PZ0QW(JL)=ZZ0QTI(JL,JTILE)
  PCPTSPP(JL)=PCPTSTI(JL,JTILE)
  PQSAPP(JL)=ZZQSATI(JL,JTILE)
  PBUOMPP(JL)=ZBUOMTI(JL,JTILE)
  PZDLPP(JL)=ZZDLTI(JL,JTILE)
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFEXCDRIVER_CTL_MOD:SURFEXCDRIVER_CTL',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE COMPUTE_DDH

! DDH diagnostics computation, mostly radiation related
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFEXCDRIVER_CTL:COMPUTE_DDH',0,ZHOOK_HANDLE)

! radiation related tiled quantities

DO JTILE=1,KTILES
  ZRVTRSR=0._JPRB
  IF (JTILE == 3 .OR. JTILE == 4) THEN
    ZRVTRSR=YDVEG%RVTRSR(1)
  ENDIF
  IF (JTILE == 6 .OR. JTILE == 7) THEN
    ZRVTRSR=YDVEG%RVTRSR(3)
  ENDIF
!   ZRVTRSR=0._JPRB  !! Added E.Dutra 1/4/2016 to allow snow EB check in offline using tiled output 
  PDHTLS(KIDIA:KFDIA,JTILE,4)=PFRTI(KIDIA:KFDIA,JTILE)*&
   & (1.-ZRVTRSR)*PSSRFL(KIDIA:KFDIA)/(1.0_JPRB-ZALB(KIDIA:KFDIA))
  PDHTLS(KIDIA:KFDIA,JTILE,5)=PFRTI(KIDIA:KFDIA,JTILE)*&
   & (1.-ZRVTRSR)*&
   & PSSRFLTI(KIDIA:KFDIA,JTILE)-PDHTLS(KIDIA:KFDIA,JTILE,4)  
  PDHTLS(KIDIA:KFDIA,JTILE,6)=PFRTI(KIDIA:KFDIA,JTILE)*&
   & (PSLRFL(KIDIA:KFDIA)/PEMIS(KIDIA:KFDIA)+&
   & YDCST%RSIGMA*PTSKTI(KIDIA:KFDIA,JTILE)**4)  
  PDHTLS(KIDIA:KFDIA,JTILE,7)=PFRTI(KIDIA:KFDIA,JTILE)*&
   & PSLRFL(KIDIA:KFDIA)/PEMIS(KIDIA:KFDIA)-PDHTLS(KIDIA:KFDIA,JTILE,6)
ENDDO

! tiles: fraction (1), temperature (2), albedo (3)
PDHTLS(KIDIA:KFDIA,:,1)=PFRTI(KIDIA:KFDIA,:)
PDHTLS(KIDIA:KFDIA,:,2)=PTSKTI(KIDIA:KFDIA,:)
PDHTLS(KIDIA:KFDIA,:,3)=PALBTI(KIDIA:KFDIA,:)
PDHTLS(KIDIA:KFDIA,:,12)=PLAMSK(KIDIA:KFDIA,:)
! snow radiative fluxes (7-10) and albedo (6) 
PDHTSS(KIDIA:KFDIA,1,7)=(PFRTI(KIDIA:KFDIA,5)+PFRTI(KIDIA:KFDIA,7))*&
 & PSSRFL(KIDIA:KFDIA)/(1.0_JPRB-ZALB(KIDIA:KFDIA))  
PDHTSS(KIDIA:KFDIA,1,8)=PFRTI(KIDIA:KFDIA,5)*PSSRFLTI(KIDIA:KFDIA,5)+&
 & PFRTI(KIDIA:KFDIA,7)*PSSRFLTI(KIDIA:KFDIA,7)-&
 & PDHTSS(KIDIA:KFDIA,1,7)  
PDHTSS(KIDIA:KFDIA,1,9)=(PFRTI(KIDIA:KFDIA,5)+PFRTI(KIDIA:KFDIA,7))*&
 & PSLRFL(KIDIA:KFDIA)/PEMIS(KIDIA:KFDIA)+&
 & YDCST%RSIGMA*(PFRTI(KIDIA:KFDIA,5)*PTSKTI(KIDIA:KFDIA,5)**4+&
 & PFRTI(KIDIA:KFDIA,7)*PTSKTI(KIDIA:KFDIA,7)**4)  
PDHTSS(KIDIA:KFDIA,1,10)=(PFRTI(KIDIA:KFDIA,5)+PFRTI(KIDIA:KFDIA,7))*&
 & PSLRFL(KIDIA:KFDIA)/PEMIS(KIDIA:KFDIA)-&
 & PDHTSS(KIDIA:KFDIA,1,9)  
WHERE (PFRTI(KIDIA:KFDIA,5)+PFRTI(KIDIA:KFDIA,7) > 0.001)
  PDHTSS(KIDIA:KFDIA,1,6)=(PFRTI(KIDIA:KFDIA,5)*PALBTI(KIDIA:KFDIA,5)+&
   & PFRTI(KIDIA:KFDIA,7)*PALBTI(KIDIA:KFDIA,7))/&
   & (PFRTI(KIDIA:KFDIA,5)+PFRTI(KIDIA:KFDIA,7))  
ELSEWHERE
  PDHTSS(KIDIA:KFDIA,1,6)=YDSOIL%RALFMAXSN
ENDWHERE
IF (KLEVSN>1) THEN
  PDHTSS(KIDIA:KFDIA,2:KLEVSN,6:10)=0.0_JPRB
ENDIF

! soil radiative fluxes (5-8)
PDHTTS(KIDIA:KFDIA,1,5)=(PFRTI(KIDIA:KFDIA,3)+PFRTI(KIDIA:KFDIA,4)+&
 & PFRTI(KIDIA:KFDIA,6)+PFRTI(KIDIA:KFDIA,8))*&
 & PSSRFL(KIDIA:KFDIA)/(1.0_JPRB-ZALB(KIDIA:KFDIA))  
PDHTTS(KIDIA:KFDIA,1,6)=PFRTI(KIDIA:KFDIA,3)*PSSRFLTI(KIDIA:KFDIA,3)+&
 & PFRTI(KIDIA:KFDIA,4)*PSSRFLTI(KIDIA:KFDIA,4)+&
 & PFRTI(KIDIA:KFDIA,6)*PSSRFLTI(KIDIA:KFDIA,6)+&
 & PFRTI(KIDIA:KFDIA,8)*PSSRFLTI(KIDIA:KFDIA,8)-&
 & PDHTTS(KIDIA:KFDIA,1,5)  
PDHTTS(KIDIA:KFDIA,1,7)=(PFRTI(KIDIA:KFDIA,3)+PFRTI(KIDIA:KFDIA,4)+&
 & PFRTI(KIDIA:KFDIA,6)+PFRTI(KIDIA:KFDIA,8))*&
 & PSLRFL(KIDIA:KFDIA)/PEMIS(KIDIA:KFDIA)+&
 & YDCST%RSIGMA*(PFRTI(KIDIA:KFDIA,3)*PTSKTI(KIDIA:KFDIA,3)**4+&
 & PFRTI(KIDIA:KFDIA,4)*PTSKTI(KIDIA:KFDIA,4)**4+&
 & PFRTI(KIDIA:KFDIA,6)*PTSKTI(KIDIA:KFDIA,6)**4+&
 & PFRTI(KIDIA:KFDIA,8)*PTSKTI(KIDIA:KFDIA,8)**4)  
PDHTTS(KIDIA:KFDIA,1,8)=(PFRTI(KIDIA:KFDIA,3)+PFRTI(KIDIA:KFDIA,4)+&
 & PFRTI(KIDIA:KFDIA,6)+PFRTI(KIDIA:KFDIA,8))*&
 & PSLRFL(KIDIA:KFDIA)/PEMIS(KIDIA:KFDIA)-&
 & PDHTTS(KIDIA:KFDIA,1,7)
! if urban
IF (KTILES>=10) THEN
 PDHTTS(KIDIA:KFDIA,1,5)=(PFRTI(KIDIA:KFDIA,3)+PFRTI(KIDIA:KFDIA,4)+&
  & PFRTI(KIDIA:KFDIA,6)+PFRTI(KIDIA:KFDIA,8)+PFRTI(KIDIA:KFDIA,10))*&
  & PSSRFL(KIDIA:KFDIA)/(1.0_JPRB-ZALB(KIDIA:KFDIA))
 PDHTTS(KIDIA:KFDIA,1,6)=PFRTI(KIDIA:KFDIA,3)*PSSRFLTI(KIDIA:KFDIA,3)+&
  & PFRTI(KIDIA:KFDIA,4)*PSSRFLTI(KIDIA:KFDIA,4)+&
  & PFRTI(KIDIA:KFDIA,6)*PSSRFLTI(KIDIA:KFDIA,6)+&
  & PFRTI(KIDIA:KFDIA,8)*PSSRFLTI(KIDIA:KFDIA,8)+&
  & PFRTI(KIDIA:KFDIA,10)*PSSRFLTI(KIDIA:KFDIA,10)-&
  & PDHTTS(KIDIA:KFDIA,1,5)
 PDHTTS(KIDIA:KFDIA,1,7)=(PFRTI(KIDIA:KFDIA,3)+PFRTI(KIDIA:KFDIA,4)+&
  & PFRTI(KIDIA:KFDIA,6)+PFRTI(KIDIA:KFDIA,8)+PFRTI(KIDIA:KFDIA,10))*&
  & PSLRFL(KIDIA:KFDIA)/PEMIS(KIDIA:KFDIA)+&
  & YDCST%RSIGMA*(PFRTI(KIDIA:KFDIA,3)*PTSKTI(KIDIA:KFDIA,3)**4+&
  & PFRTI(KIDIA:KFDIA,4)*PTSKTI(KIDIA:KFDIA,4)**4+&
  & PFRTI(KIDIA:KFDIA,6)*PTSKTI(KIDIA:KFDIA,6)**4+&
  & PFRTI(KIDIA:KFDIA,8)*PTSKTI(KIDIA:KFDIA,8)**4+&
  & PFRTI(KIDIA:KFDIA,10)*PTSKTI(KIDIA:KFDIA,10)**4)
 PDHTTS(KIDIA:KFDIA,1,8)=(PFRTI(KIDIA:KFDIA,3)+PFRTI(KIDIA:KFDIA,4)+&
  & PFRTI(KIDIA:KFDIA,6)+PFRTI(KIDIA:KFDIA,8)+PFRTI(KIDIA:KFDIA,10))*&
  & PSLRFL(KIDIA:KFDIA)/PEMIS(KIDIA:KFDIA)-&
  & PDHTTS(KIDIA:KFDIA,1,7)
ENDIF
 
PDHTTS(KIDIA:KFDIA,2:KLEVS,5:8)=0.0_JPRB
! ice radiative fluxes (5-8)
PDHTIS(KIDIA:KFDIA,1,5)=PFRTI(KIDIA:KFDIA,2)*PSSRFL(KIDIA:KFDIA)/&
 & (1.0_JPRB-ZALB(KIDIA:KFDIA))  
PDHTIS(KIDIA:KFDIA,1,6)=PFRTI(KIDIA:KFDIA,2)*PSSRFLTI(KIDIA:KFDIA,2)-&
 & PDHTIS(KIDIA:KFDIA,1,5)  
PDHTIS(KIDIA:KFDIA,1,7)=PFRTI(KIDIA:KFDIA,2)*&
 & PSLRFL(KIDIA:KFDIA)/PEMIS(KIDIA:KFDIA)+&
 & YDCST%RSIGMA*PFRTI(KIDIA:KFDIA,2)*PTSKTI(KIDIA:KFDIA,2)**4  
PDHTIS(KIDIA:KFDIA,1,8)=PFRTI(KIDIA:KFDIA,2)*&
 & PSLRFL(KIDIA:KFDIA)/PEMIS(KIDIA:KFDIA)-&
 & PDHTIS(KIDIA:KFDIA,1,7)  
PDHTIS(KIDIA:KFDIA,2:KLEVS,5:8)=0.0_JPRB


PDHVEGS(KIDIA:KFDIA,1,1)=PCVL(KIDIA:KFDIA)
PDHVEGS(KIDIA:KFDIA,2,1)=PCVH(KIDIA:KFDIA)


! latent heat for vegetation tiles without snow
DO JL=KIDIA,KFDIA
   DO JTILE=4,4
      IF (PFRTI(JL,JTILE)/= 0._JPRB) THEN
         PDHVEGS(JL,KVTTL(JTILE),8)=PEVAPTI(JL,JTILE)
      ELSE
         PDHVEGS(JL,KVTTL(JTILE),8)=0._JPRB
      ENDIF
   ENDDO
   DO JTILE=6,6
      IF (PFRTI(JL,JTILE)/= 0._JPRB) THEN    
         PDHVEGS(JL,KVTTL(JTILE),8)=PEVAPTI(JL,JTILE)
      ELSE
         PDHVEGS(JL,KVTTL(JTILE),8)=0._JPRB
      ENDIF
   ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('SURFEXCDRIVER_CTL:COMPUTE_DDH',1,ZHOOK_HANDLE)

END SUBROUTINE COMPUTE_DDH

END SUBROUTINE SURFEXCDRIVER_CTL
END MODULE SURFEXCDRIVER_CTL_MOD
