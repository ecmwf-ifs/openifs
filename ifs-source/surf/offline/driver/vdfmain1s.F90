!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFMAIN1S  ( CDCONF, &
 & KIDIA  , KFDIA  , KLON   , KLEV   , KLEVS  , KSTEP  , KTILES , KVTYPES, KDIAG, &
 & KTRAC  , KLEVSN , KLEVI  , KDHVTLS, KDHFTLS, KDHVTSS, KDHFTSS, &
 & KDHVTTS, KDHFTTS, KDHVTIS, KDHFTIS, &
 & KDHVCO2S,KDHFCO2S, KDHVVEGS,KDHFVEGS,  &
 & PMU0   , &
 & PTSPHY , KTVL   , KCO2TYP, KTVH   , PCVL   , PCVH   , PCUR, PFWET  , PLAT  , PLAIL  , PLAIH  , &
 & PSNM1M , PRSNM1M , &
 & PUM1   , PVM1   , PTM1   , PQM1   , PCM1   , &
 & PAPHM1 , PGEOM1 , PTSKM1M, PTSAM1M, PWSAM1M, &
 & PSSRFL , PSLRFL , PEMIS  , &
 & PTSN , PTICE  , &
 & PHLICE , PTLICE , PTLWML , &
 & PSST   , KSOTY  , PFRTI  , PALBTI , PWLMX  , &
 & PCHAR  , PCHARHQ, PUCURR , PVCURR , PTSKRAD, PCFLX  , &
 ! OUTPUT
 & PZ0M   , PZ0H   , &
 & PVDIS  , PAHFLEV, PAHFLSB, PFWSB  , &
 & PU10M  , PV10M  , PT2M   , PD2M   , PQ2M   , PZINV  , &
 & PSSRFLTI,PEVAPSNW,PEXDIAG,PGUST  , PZIDLWV, &
 ! OUTPUT TENDENCIES
 & PTE    , PQE    , PVOM   , PVOL   , &
 & PTENC  , PTSKE1 , &
 ! UPDATED FIELDS FOR TILES
 & PUSTRTI, PVSTRTI, PAHFSTI, PEVAPTI, PTSKTI , PSLRFLTI,&
  !-UPDATED FIELDS FOR VEGETATION TYPES
 & PANDAYVT,PANFMVT,&
 ! OUTPUT FLUXES
 & PEVAPTIU,PDIFTS , PDIFTQ , PSTRTU , PSTRTV , PKH, &
  & PAN,PAG,PRD,PRSOIL_STR,PRECO,PCO2FLUX,PCH4FLUX,& 
 ! DDH OUTPUTS
 & PDHTLS , PDHTSS , PDHTTS , PDHTIS, PDHCO2S,PDHVEGS )

! (C) Copyright 1982- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!***

!**   *VDFMAIN1S* - DOES THE VERTICAL EXCHANGE OF U,V,SLG,QT BY TURBULENCE.

!     J.F.GELEYN       20/04/82   Original  
!     C.A.BLONDIN      18/12/86
!     A.C.M. BELJAARS  20/10/89   IFS-VERSION (TECHNICAL REVISION OF CY34)
!     A.C.M. BELJAARS  26/03/90   OBUKHOV-L UPDATE 
!     A.C.M. BELJAARS  30/09/98   SURFACE TILES 
!     P. Viterbo       17/05/2000 Surface DDH for TILES
!     D. Salmond       15/10/2001 FULLIMP mods
!     S. Abdalla       27/11/2001 Passing Zi/L to waves
!     A. Beljaars       2/05/2003 New tile coupling     
!     P.Viterbo        24/05/2004 Change surface units
!     M. Ko"hler        3/12/2004 Moist Advection-Diffusion
!     P. Viterbo       17/06/2005 surf external library
!     G. Balsamo       03/07/2006 add soil type
!     G. Balsamo       25/08/2009 add lai clim
!     S. Boussetta/G.Balsamo May 2010 Add CTESSEL
!     E. Dutra         10/10/2014  net longwave tiled 
!     A. Agusti-Panareda 09/04/2021 atmospheric CO2 passed to land surface

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE FOUR
!     PROGNOSTIC VARIABLES U,V,T AND Q DUE TO THE VERTICAL EXCHANGE BY
!     TURBULENT (= NON-MOIST CONVECTIVE) PROCESSES. THESE TENDENCIES ARE
!     OBTAINED AS THE DIFFERENCE BETWEEN THE RESULT OF AN IMPLICIT
!     TIME-STEP STARTING FROM VALUES AT T-1 AND THESE T-1 VALUES. ALL
!     THE DIAGNOSTIC COMPUTATIONS (EXCHANGE COEFFICIENTS, ...) ARE DONE
!      FROM THE T-1 VALUES. AS A BY-PRODUCT THE ROUGHNESS LENGTH OVER SEA
!     IS UPDATED ACCORDINGLY TO THE *CHARNOCK FORMULA. HEAT AND MOISTURE
!     SURFACE FLUXES AND THEIR DERIVATIVES AGAINST TS, WS AND WL
!     (THE LATTER WILL BE LATER WEIGHTED WITH THE SNOW FACTOR IN
!     *VDIFF*), LATER TO BE USED FOR SOIL PROCESSES TREATMENT, ARE ALSO
!     COMPUTED AS WELL AS A STABILITY VALUE TO BE USED AS A DIAGNOSTIC
!     OF THE DEPTH OF THE WELL MIXED LAYER IN CONVECTIVE COMPUTATIONS.

!     INTERFACE.
!     ----------
!          *VDIFF* TAKES THE MODEL VARIABLES AT T-1 AND RETURNS THE VALUES
!     FOR THE PROGNOSTIC TIME T+1 DUE TO VERTICAL DIFFUSION.
!     THE MODEL VARIABLES, THE MODEL DIMENSIONS AND THE DIAGNOSTICS DATA
!     ARE PASSED AS SUBROUTINE ARGUMENTS. CONSTANTS THAT DO NOT CHANGE
!     DURING A MODEL RUN (E.G. PHYSICAL CONSTANTS, SWITCHES ETC.) ARE
!     STORED IN A SINGLE COMMON BLOCK *YOMVDF*, WHICH IS INITIALIZED
!     BY SET-UP ROUTINE *SUVDF*.

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLEV*         NUMBER OF LEVELS
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEVS*        NUMBER OF SOIL LAYERS
!    *KSTEP*        CURRENT TIME STEP INDEX
!    *KTILES*       NUMBER OF TILES (I.E. SUBGRID AREAS WITH DIFFERENT 
!                   OF SURFACE BOUNDARY CONDITION)
!    *KVTYPES*      NUMBER OF biomes for land carbon
!    *KDIAG*        NUMBER of diagnostic parameters
!    *KTRAC*        Number of tracers
!    *KLEVSN*       Number of snow layers (diagnostics) 
!    *KLEVI*        Number of sea ice layers (diagnostics)
!    *KDHVTLS*      Number of variables for individual tiles
!    *KDHFTLS*      Number of fluxes for individual tiles
!    *KDHVTSS*      Number of variables for snow energy budget
!    *KDHFTSS*      Number of fluxes for snow energy budget
!    *KDHVTTS*      Number of variables for soil energy budget
!    *KDHFTTS*      Number of fluxes for soil energy budget
!    *KDHVTIS*      Number of variables for sea ice energy budget
!    *KDHFTIS*      Number of fluxes for sea ice energy budget

!    *KTVL*         VEGETATION TYPE FOR LOW VEGETATION FRACTION
!    *KCO2TYP*       TYPE OF PHOTOSYNTHESIS PATHWAY FOR LOW VEGETATION (C3/C4)
!    *KTVH*         VEGETATION TYPE FOR HIGH VEGETATION FRACTION
!    *KSOTY*        SOIL TYPE                                     (1-6)

!    *KCNT*         Index of vdf sub steps.

!     INPUT PARAMETERS (LOGICAL)

!     INPUT PARAMETERS (REAL)

!    *PTSPHY*       TIME STEP FOR THE PHYSICS

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):

!    *PCVL*         LOW VEGETATION COVER                          -  
!    *PCVH*         HIGH VEGETATION COVER                         -  
!    *PCUR*         URBAN COVER (PASSIVE)                        (0-1)
!    *PLAIL*        LOW VEGETATION LAI                           m2/m2
!    *PLAIH*        HIGH VEGETATION LAI                          m2/m2

!     *PSNM1M*       SNOW MASS (per unit area)                      kg/m**2
!     *PRSNM1M*      SNOW DENSITY                                   kg/m**3


!    *PUM1*         X-VELOCITY COMPONENT                          M/S
!    *PVM1*         Y-VELOCITY COMPONENT                          M/S
!    *PTM1*         TEMPERATURE                                   K
!    *PQM1*         SPECIFIC HUMIDITY                             KG/KG
!    *PLM1*         SPECIFIC CLOUD LIQUID WATER                   KG/KG
!    *PIM1*         SPECIFIC CLOUD ICE                            KG/KG
!    *PAM1*         CLOUD FRACTION                                1
!    *PCM1*         TRACER CONCENTRATION                          KG/KG
!    *PAPM1*        PRESSURE ON FULL LEVELS                       PA
!    *PAPHM1*       PRESSURE ON HALF LEVELS                       PA
!    *PGEOM1*       GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL AT HALF LEVELS                   M2/S2
!    *PTSKM1M*      SKIN TEMPERATURE                              K
!    *PTSAM1M*      SURFACE TEMPERATURE                           K
!    *PWSAM1M*      SOIL MOISTURE ALL LAYERS                      M**3/M**3
!    *PSSRFL*       NET SHORTWAVE RADIATION FLUX AT SURFACE       W/M2
!    *PSLRFL*       NET LONGWAVE RADIATION FLUX AT SURFACE        W/M2
!    *PEMIS*        MODEL SURFACE LONGWAVE EMISSIVITY
!    *PHRLW*        LONGWAVE HEATING RATE                         K/s
!    *PHRSW*        SHORTWAVE HEATING RATE                        K/s
!    *PTSN*         SNOW TEMPERATURE                              K
!    *PTICE*        ICE TEMPERATURE (TOP SLAB)                    K
!    *PHLICE*       LAKE ICE THICKNESS                            m
!    *PTLICE*       LAKE ICE TEMPERATURE                          K
!    *PTLWML*       LAKE MEAN WATER TEMPERATURE                   K
!    *PSST*         (OPEN) SEA SURFACE TEMPERATURE                K
!    *PFRTI*        TILE FRACTIONS                                (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!    *PALBTI*       BROADBAND ALBEDO FOR TILE FRACTIONS
!    *PWLMX*        MAXIMUM SKIN RESERVOIR CAPACITY               kg/m**2
!    *PCHAR*        CHARNOCK PARAMETER                            -
!    *PCHARHQ       CHARNOCK PARAMETER FOR HEAT AND MOISTURE      -
!    *PUCURR*       OCEAN CURRENT X_COMPONENT
!    *PVCURR*       OCEAN CURRENT Y_COMPONENT
!    *PUSTOKES*     SURFACE STOKES VELOCITY X_COMPONENT
!    *PVSTOKES*     SURFACE STOKES VELOCITY Y_COMPONENT
!    *PTSKRAD*      SKIN TEMPERATURE OF LATEST FULL RADIATION
!                      TIMESTEP                                   K
!    *PCFLX*        TRACER SURFACE FLUX                           kg/(m2 s)

!    *PMU0*         LOCAL COSINE OF INSTANTANEOUS MEAN SOLAR ZENITH ANGLE


!     INPUT PARAMETERS (LOGICAL):

!     CONTRIBUTIONS TO BUDGETS (OUTPUT,REAL):

!    *PVDIS*        DISSIPATION                                   W/M2
!    *PAHFLEV*      LATENT HEAT FLUX  (SNOW/ICE FREE PART)        W/M2
!    *PAHFLSB*      LATENT HEAT FLUX  (SNOW/ICE COVERED PART)     W/M2

!     UPDATED PARAMETERS (REAL):

!    *PTE*          TEMPERATURE TENDENCY                          K/S
!    *PQE*          MOISTURE TENDENCY                             KG/(KG S)
!    *PLE*          LIQUID WATER TENDENCY                         KG/(KG S)
!    *PIE*          ICE WATER TENDENCY                            KG/(KG S)
!    *PAE*          CLOUD FRACTION TENDENCY                       1/S)
!    *PVOM*         MERIODINAL VELOCITY TENDENCY (DU/DT)          M/S2
!    *PVOL*         LATITUDE TENDENCY            (DV/DT)          M/S2
!    *PTENC*        TRACER TENDENCY                               KG/(KG S)
!    *PTSKE1*       SKIN TEMPERATURE TENDENCY                     K/S
!    *PZ0M*         AERODYNAMIC ROUGHNESS LENGTH                  M
!    *PZ0H*         ROUGHNESS LENGTH FOR HEAT                     M

!     UPDATED PARAMETERS FOR TILES (REAL): 

!    *PUSTRTI*      SURFACE U-STRESS                              N/M2 
!    *PVSTRTI*      SURFACE V-STRESS                              N/M2
!    *PAHFSTI*      SURFACE SENSIBLE HEAT FLUX                    W/M2
!    *PEVAPTI*      SURFACE MOISTURE FLUX                         KG/M2/S
!    *PTSKTI*       SKIN TEMPERATURE                              K
!    *PSLRFLTI*     Net longwave radiation                        W/m2 (only output here)

!    UPDATED PARAMETERS FOR VEGETATION TYPES (REAL): 

!    *PANDAYVT*     DAILY NET CO2 ASSIMILATION OVER CANOPY    KG_CO2/M2
!    *PANFMVT*      MAXIMUM LEAF ASSIMILATION                KG_CO2/KG_AIR M/S  

!     OUTPUT PARAMETERS (REAL):

!    *PEVAPTIU*      SURFACE MOISTURE FLUX (unstressed low veg)    KG/M2/S
!    *PFWSB*        EVAPORATION OF SNOW                           KG/(M**2*S)
!    *PU10M*        U-COMPONENT WIND AT 10 M                      M/S
!    *PV10M*        V-COMPONENT WIND AT 10 M                      M/S
!    *P10NU*        U-COMPONENT NEUTRAL WIND AT 10 M              M/S
!    *P10NV*        V-COMPONENT NEUTRAL WIND AT 10 M              M/S
!    *PUST*         FRICTION VELOCITY                             M/S
!    *PT2M*         TEMPERATURE AT 2M                             K
!    *PD2M*         DEW POINT TEMPERATURE AT 2M                   K
!    *PQ2M*         SPECIFIC HUMIDITY AT 2M                       KG/KG
!    *PEXDIAG*      Extra diagnostic fields
!    *PGUST*        GUST AT 10 M                                  M/S
!    *PZIDLWV*      Zi/L used for gustiness in wave model         M/M
!                   (NOTE: Positive values of Zi/L are set to ZERO)
!    *PBLH*         PBL HEIGHT (dry diagnostic based on Ri#)      M
!    *PZINV*        PBL HEIGHT (moist parcel, not for stable PBL) M
!    *PSSRFLTI*     NET SHORTWAVE RADIATION FLUX AT SURFACE, FOR
!                      EACH TILE                                  W/M2
!    *PEVAPSNW*     EVAPORATION FROM SNOW UNDER FOREST            KG/(M2*S)
!    *PSTRTU*       TURBULENT FLUX OF U-MOMEMTUM            KG*(M/S)/(M2*S)
!    *PSTRTV*       TURBULENT FLUX OF V-MOMEMTUM            KG*(M/S)/(M2*S)
!    *PDIFTS*       TURBULENT FLUX OF HEAT                         J/(M2*S)
!    *PDIFTQ*       TURBULENT FLUX OF SPECIFIC HUMIDITY           KG/(M2*S)
!    *PDIFTL*       TURBULENT FLUX OF LIQUID WATER                KG/(M2*S)
!    *PDIFTI*       TURBULENT FLUX OF ICE WATER                   KG/(M2*S)

!    *PKH*          TURB. DIFF. COEFF. FOR HEAT ABOVE SURF. LAY.  (M2/S)
!                   IN SURFACE LAYER: CH*U                        (M/S)
!    *PDHTLS*       Diagnostic array for tiles (see module yomcdh)
!                      (Wm-2 for energy fluxes, kg/(m2s) for water fluxes)
!    *PDHTSS*       Diagnostic array for snow T (see module yomcdh)
!                      (Wm-2 for fluxes)
!    *PDHTTS*       Diagnostic array for soil T (see module yomcdh)
!                      (Wm-2 for fluxes)
!    *PDHTIS*       Diagnostic array for ice T (see module yomcdh)
!                      (Wm-2 for fluxes)


!    *PAN*          NET CO2 ASSIMILATION OVER CANOPY          KG_CO2/M2/S
!    *PAG*          GROSS CO2 ASSIMILATION OVER CANOPY        KG_CO2/M2/S
!    *PRD*          DARK RESPIRATION                          KG_CO2/M2/S
!    *PRSOIL_STR*   RESPIRATION FROM SOIL AND STRUCTURAL BIOMASS KG_CO2/M2/S
!    *PRECO*        ECOSYSTEM RESPIRATION                     KG_CO2/M2/S
!    *PCO2FLUX*     CO2 FLUX                                  KG_CO2/M2/S
!    *PCH4FLUX*     CH4 FLUX                                  KG_CO2/M2/S

!     *PSNM1M*       SNOW MASS (per unit area)                      kg/m**2
!     *PRSNM1M*      SNOW DENSITY                                   kg/m**3
!     Additional parameters for flux boundary condtion (in 1D model):

!    *LLSFCFLX*     If .TRUE. flux boundary condtion is used 
!    *ZFSH1D*       Specified sensible heat flux (W/m2)
!    *ZFLH1D*       Specified latent heat flux (W/m2)

!     METHOD.
!     -------

!          FIRST AN AUXIALIARY VARIABLE CP(Q)T+GZ IS CREATED ON WHICH
!     THE VERTICAL DIFFUSION PROCESS WILL WORK LIKE ON U,V AND Q. THEN
!     ALONG THE VERTICAL AND AT THE SURFACE, EXCHANGE COEFFICIENTS (WITH
!     THE DIMENSION OF A PRESSURE THICKNESS) ARE COMPUTED FOR MOMENTUM
!     AND FOR HEAT (SENSIBLE PLUS LATENT). THE LETTERS M AND H ARE USED
!     TO DISTINGUISH THEM AND THE COMPUTATION IS THE RESULT OF A
!     CONDITIONAL MERGE BETWEEN THE STABLE AND THE UNSTABLE CASE
!     (DEPENDING ON THE SIGN OF THE *RICHARDSON BULK NUMBER).
!          IN THE SECOND PART OF THE ROUTINE THE IMPLICIT LINEAR
!     SYSTEMS FOR U,V FIRST AND T,Q SECOND ARE SOLVED BY A *GAUSSIAN
!     ELIMINATION BACK-SUBSTITUTION METHOD. FOR T AND Q THE LOWER
!     BOUNDARY CONDITION DEPENDS ON THE SURFACE STATE.
!     OVER LAND, TWO DIFFERENT REGIMES OF EVAPORATION PREVAIL:
!     A STOMATAL RESISTANCE DEPENDENT ONE OVER THE VEGETATED PART
!     AND A SOIL RELATIVE HUMIDITY DEPENDENT ONE OVER THE
!     BARE SOIL PART OF THE GRID MESH.
!     POTENTIAL EVAPORATION TAKES PLACE OVER THE SEA, THE SNOW
!     COVERED PART AND THE LIQUID WATER COVERED PART OF THE
!     GRID MESH AS WELL AS IN CASE OF DEW DEPOSITION.
!          FINALLY ONE RETURNS TO THE VARIABLE TEMPERATURE TO COMPUTE
!     ITS TENDENCY AND THE LATER IS MODIFIED BY THE DISSIPATION'S EFFECT
!     (ONE ASSUMES NO STORAGE IN THE TURBULENT KINETIC ENERGY RANGE) AND
!     THE EFFECT OF MOISTURE DIFFUSION ON CP. Z0 IS UPDATED AND THE
!     SURFACE FLUXES OF T AND Q AND THEIR DERIVATIVES ARE PREPARED AND
!     STORED LIKE THE DIFFERENCE BETWEEN THE IMPLICITELY OBTAINED
!     CP(Q)T+GZ AND CP(Q)T AT THE SURFACE.

!     EXTERNALS.
!     ----------

!     *VDFMAIN1S* CALLS SUCESSIVELY:
!         *SURFEXCDRIVER*
!         *VDFEXCU*
!         *VDFDIFM*
!         *VDFDIFH*
!         *VDFDIFC*
!         *VDFINCR*
!         *VDFSDRV*
!         *VDFPPCFL*
!         *VDFUPDZ0*

!     REFERENCE.
!     ----------

!          SEE VERTICAL DIFFUSION'S PART OF THE MODEL'S DOCUMENTATION
!     FOR DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB       ,JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

USE YOEVDF   , ONLY : RVDIFTS
USE YOMCST   , ONLY : RG       ,RD       ,&
                    & RCPD     ,RETV     ,RLVTT    ,RLSTT    ,RTT 
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
                    & R4IES    ,R5LES    ,R5IES    ,RVTMP2   ,R5ALVCP  ,&
                    & R5ALSCP  ,RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,&
                    & RTICECU  ,RTWAT_RTICE_R      ,RTWAT_RTICECU_R  
USE YOMJFH   , ONLY : N_VMASS
USE YOEPHY   , ONLY : LVDFTRAC 
USE YOMGF1S  , ONLY : RALT

USE YOERDI   , ONLY : RCARDI ! C-TESSEL
USE YOMDPHY  , ONLY : YDSURF
USE YOMCST   , ONLY : RSIGMA


IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVTYPES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIAG
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
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
CHARACTER(LEN=1)  ,INTENT(IN)    :: CDCONF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCO2TYP(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(KLON) 

INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVCO2S 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFCO2S
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVVEGS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFVEGS


INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTY(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCUR(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAIL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAIH(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFWET(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT(KLON)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(KLON)

REAL(KIND=JPRB),    INTENT(IN)  :: PSNM1M(KLON,KLEVSN)
REAL(KIND=JPRB),    INTENT(IN)  :: PRSNM1M(KLON,KLEVSN)
 
!1s INTEGER(KIND=JPIM),INTENT(IN)    :: KCNT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCM1(KLON,KLEV,KTRAC) 
!1s REAL(KIND=JPRB)   ,INTENT(IN)    :: PLM1(KLON,KLEV)
!1s REAL(KIND=JPRB)   ,INTENT(IN)    :: PIM1(KLON,KLEV)
!1s REAL(KIND=JPRB)   ,INTENT(IN)    :: PAM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV)
!1s REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV)
!1s REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSAM1M(KLON,KLEVS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSAM1M(KLON,KLEVS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) 
!1s REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRLW(KLON,KLEV) 
!1s REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRSW(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSN(KLON,KLEVSN) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTICE(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICE(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLICE(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLWML(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWLMX(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHAR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHARHQ(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKRAD(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFLX(KLON,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0M(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0H(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAHFLEV(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAHFLSB(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFWSB(KLON) 
!1s REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBIR(KLON)
!1s REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVAR(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU10M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PV10M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT2M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD2M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQ2M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZINV(KLON)
!1s REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBLH(KLON) 
!1s INTEGER(KIND=JPIM),INTENT(OUT)   :: KHPBLN(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSSRFLTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEVAPSNW(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXDIAG(KLON,KDIAG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGUST(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZIDLWV(KLON) 
!1s REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWUAVG(KLON)
!1s LOGICAL           ,INTENT(IN)    :: LDNODECP(KLON)
!1s INTEGER(KIND=JPIM),INTENT(OUT)   :: KPBLTYPE(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQE(KLON,KLEV) 
!1s REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLE(KLON,KLEV) 
!1s REAL(KIND=JPRB)   ,INTENT(INOUT) :: PIE(KLON,KLEV) 
!1s REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKE1(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSTRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSTRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEVAPTIU(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKTI(KLON,KTILES)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSLRFLTI(KLON,KTILES)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PANDAYVT(KLON,KVTYPES)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PANFMVT(KLON,KVTYPES)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTS(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTQ(KLON,0:KLEV) 
!1s REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTL(KLON,0:KLEV) 
!1s REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTI(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRTU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRTV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTLS(KLON,KTILES,KDHVTLS+KDHFTLS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTSS(KLON,KLEVSN,KDHVTSS+KDHFTSS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTTS(KLON,KLEVS,KDHVTTS+KDHFTTS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTIS(KLON,KLEVI,KDHVTIS+KDHFTIS) 


!# O PAN is defined as a vector it is a surface quantity
!#not a quantity of the lowest atmospherci level.
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAN(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAG(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRD(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRSOIL_STR(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRECO(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCO2FLUX(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCH4FLUX(KLON)

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHVEGS(KLON,KVTYPES,KDHVVEGS+KDHFVEGS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHCO2S(KLON,KVTYPES,KDHVCO2S+KDHFCO2S)

!*         0.2    LOCAL VARIABLES

REAL(KIND=JPRB) ::    ZRAQTI(KLON,KTILES)
REAL(KIND=JPRB) ::    ZKCLEV(KLON)

REAL(KIND=JPRB) ::    ZWETB(KLON)
REAL(KIND=JPRB) ::    ZWETL(KLON)
REAL(KIND=JPRB) ::    ZWETLU(KLON)
REAL(KIND=JPRB) ::    ZWETH(KLON)
REAL(KIND=JPRB) ::    ZWETHS(KLON)

REAL(KIND=JPRB) ::    ZDIFTQT(KLON,0:KLEV), ZDIFTSLG(KLON,0:KLEV) 

REAL(KIND=JPRB) ::    ZCPTGZ(KLON,KLEV) , ZCFM(KLON,KLEV)   , ZCFH(KLON,KLEV)   ,&
                    & ZUDIF(KLON,KLEV)  , ZVDIF(KLON,KLEV)  ,&
                    & ZQTDIF(KLON,KLEV) , ZSLGDIF(KLON,KLEV),&
                    & ZSLGM1(KLON,KLEV) , ZQTM1(KLON,KLEV)  , ZQTE(KLON,KLEV)   ,&
                    & ZSLGE(KLON,KLEV)    
REAL(KIND=JPRB) ::    ZKHFL(KLON)       , ZKQFL(KLON)       , ZKMFL(KLON)  
REAL(KIND=JPRB) ::    ZQEA(KLON,KLEV)   , ZLEA(KLON,KLEV)   , ZIEA(KLON,KLEV)   ,&
                    & ZQTEA(KLON,KLEV)  , ZSLGEA(KLON,KLEV) , ZAEA(KLON,KLEV)   ,&
                    & ZTEA(KLON,KLEV)   , ZUEA(KLON,KLEV)   , ZVEA(KLON,KLEV)   ,&
                    & ZSLGEWODIS(KLON,KLEV)  
REAL(KIND=JPRB) ::    ZZ0MW(KLON)       , ZZ0HW(KLON)       , ZZ0QW(KLON)       ,&
                    & ZBLEND(KLON)      , ZFBLEND(KLON)
REAL(KIND=JPRB) ::    ZZCPTS(KLON)      , ZZQSA(KLON)       , ZZBUOM(KLON)      ,&
                    & ZZZDL(KLON)
REAL(KIND=JPRB) ::    ZTUPD(KLON,KLEV)  , ZQUPD(KLON,KLEV)  , ZLUPD(KLON,KLEV)  ,&
                    & ZIUPD(KLON,KLEV)  , ZQTUPD(KLON,KLEV) , ZLIUPD(KLON,KLEV) ,&
                    & ZSLGUPD(KLON,KLEV), ZAUPD(KLON,KLEV)  
REAL(KIND=JPRB) ::    ZTINI(KLON,KLEV)  , ZQINI(KLON,KLEV)  , ZLINI(KLON,KLEV)  ,&
                    & ZIINI(KLON,KLEV)  , ZAINI(KLON,KLEV)  
REAL(KIND=JPRB) ::    ZVARGEN           , ZTAU(KLON)        ,&
                    & ZQSVAR(KLON,KLEV) , ZANEW(KLON,KLEV)  , ZLNEW(KLON,KLEV)  
REAL(KIND=JPRB) ::    ZSVFLUXCLD(KLON,0:KLEV)               , ZSVFLUXSUB(KLON,0:KLEV),&
                    & ZSVFLUX(KLON,0:KLEV),ZBUOYPOS(KLON)   , ZBUOYNEG(KLON)    ,&
                    & ZDZH(KLON,0:KLEV)
REAL(KIND=JPRB) ::    ZALFA1            , ZALFA2            , ZDELQ             ,&
                    & ZCORQS(KLON,KLEV) , ZDQSDTEMP(KLON,KLEV)
REAL(KIND=JPRB) ::    ZMFLX(KLON,0:KLEV),&
                    & ZQTUH(KLON,0:KLEV), ZSLGUH(KLON,0:KLEV), ZCLDBASE(KLON)
REAL(KIND=JPRB) ::    ZFSH1D(KLON)      , ZFLH1D(KLON)

REAL(KIND=JPRB) ::    ZCPTSTI(KLON,KTILES), ZQSTI(KLON,KTILES)  ,&
                    & ZDQSTI(KLON,KTILES) , ZCSATTI(KLON,KTILES),&
                    & ZCAIRTI(KLON,KTILES), ZCFHTI(KLON,KTILES) ,&
                    & ZCFQTI(KLON,KTILES) , ZAHFLTI(KLON,KTILES),&
                    & ZTSKTIP1(KLON,KTILES),ZCPTSTIU(KLON,KTILES),&
                    & ZCSATTIU(KLON,KTILES),ZCAIRTIU(KLON,KTILES),&
                    & ZTSRF(KLON,KTILES),ZLAMSK(KLON,KTILES)
REAL(KIND=JPRB) ::    ZG0(KLON,KTILES),&
                    & Z10NU(KLON), Z10NV(KLON), ZUST(KLON),&
                    & ZUSTOKES(KLON), ZVSTOKES(KLON) , PI10FGCV(KLON)

LOGICAL ::            LLRUNDRY(KLON)

LOGICAL ::            LDEBUGPRINT

INTEGER(KIND=JPIM) :: ITOP, JK, JL, JT

REAL(KIND=JPRB) ::    ZGDPH, ZRHO, ZTMST, ZRG, ZRTMST,ZWORK
LOGICAL ::            LLSFCFLX
REAL(KIND=JPRB) ::    ZEXTSHF, ZEXTLHF
REAL(KIND=JPRB) ::    ZALFAW(KLON,KLEV), ZFACW, ZFACI, ZFAC, ZESDP, ZCOR
REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE

REAL(KIND=JPRB) ::    ZTHKICE(KLON),ZSNTICE(KLON),PRPLRG
REAL(KIND=JPRB) ::    ZATMCO2(KLON)

!LOGICAL         ::    LNEMOLIMTHK
LOGICAL         ::    LSICOUP

#include "surfexcdriver.h"
#include "surfpp.h"

!1s #include "cloudvar.intfb.h"
!1s #include "cover.intfb.h"
#include "vdfdifh1s.intfb.h"
#include "vdfdifm1s.intfb.h"
#include "vdfdifc.intfb.h"
!gc #include "vdfdpbl.intfb.h"
!gc #include "vdfexcu.intfb.h"
!1s #include "vdfhghtn.intfb.h"
#include "vdfincr.intfb.h"
!1s #include "vdffblend.intfb.h"

#include "fcttre.h"

LDEBUGPRINT=.FALSE.
!     ------------------------------------------------------------------

!*         1.     INITIALIZE CONSTANTS
!                 --------------------

IF (LHOOK) CALL DR_HOOK('VDFMAIN1S',0,ZHOOK_HANDLE)

ZTMST       = PTSPHY
ZRTMST      = 1.0_JPRB/PTSPHY    ! optimization
ZRG         = 1.0_JPRB/RG        !     -"-
LLRUNDRY(:) = .FALSE.  ! option to run dry updrafts with no condensation
!DO JL=KIDIA,KFDIA
! IF ( .NOT. LDNODECP(JL) )  LLRUNDRY(JL) = .TRUE. ! run dry for real vdfmain's
! IF ( .NOT. LDNODECP(JL) )  PBIR(JL) = 1.0        ! always decouple for real vdfmain's
!ENDDO


! CTESSEL initialization 
PANDAYVT(KIDIA:KFDIA,:)=0.0_JPRB
PANFMVT(KIDIA:KFDIA,:)=0.0_JPRB
PDHVEGS(KIDIA:KFDIA,:,:)=0.0_JPRB
PDHCO2S(KIDIA:KFDIA,:,:)=0.0_JPRB

PAN(KIDIA:KFDIA)=0.0_JPRB
PAG(KIDIA:KFDIA)=0.0_JPRB
PRD(KIDIA:KFDIA)=0.0_JPRB
PRSOIL_STR(KIDIA:KFDIA)=0.0_JPRB
PRECO(KIDIA:KFDIA)=0.0_JPRB
PCO2FLUX(KIDIA:KFDIA)=0.0_JPRB
PCH4FLUX(KIDIA:KFDIA)=0.0_JPRB

!*         1.0  Fixed fluxes for flux boundary condition ([W/m^2] downward)

LLSFCFLX = .FALSE.
ZEXTSHF  = 0.0   ! SPECIFY SURFACE FLUX SENSIBLE
ZEXTLHF  = 0.0   ! SPECIFY SURFACE FLUX LATENT
Z10NU(KIDIA:KFDIA)=0. 
Z10NV(KIDIA:KFDIA)=0.
ZUST(KIDIA:KFDIA)=0.
ZUSTOKES(KIDIA:KFDIA)=0.
ZVSTOKES(KIDIA:KFDIA)=0.

IF (LLSFCFLX) THEN
  ZFSH1D(KIDIA:KFDIA) = ZEXTSHF
  ZFLH1D(KIDIA:KFDIA) = ZEXTLHF
ELSE
  ZFSH1D(KIDIA:KFDIA) = 0.0_JPRB
  ZFLH1D(KIDIA:KFDIA) = 0.0_JPRB
ENDIF


!*         1.1  Store initial tendencies for flux calculation
!*              and initialize variable.

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZQEA(JL,JK)=PQE(JL,JK)
!1s    ZLEA(JL,JK)=PLE(JL,JK)
!1s    ZIEA(JL,JK)=PIE(JL,JK)
!1s    ZAEA(JL,JK)=PAE(JL,JK)
    ZTEA(JL,JK)=PTE(JL,JK)
    ZUEA(JL,JK)=PVOM(JL,JK)
    ZVEA(JL,JK)=PVOL(JL,JK)
  ENDDO
ENDDO

DO JK=0,KLEV
  DO JL=KIDIA,KFDIA
    ZMFLX(JL,JK)  = 0.0_JPRB
    ZSLGUH(JL,JK) = 0.0_JPRB
    ZQTUH(JL,JK)  = 0.0_JPRB
  ENDDO
ENDDO

!*         1.2  dry static energy cp(q)*T + gz

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZCPTGZ(JL,JK)  =PGEOM1(JL,JK)+PTM1(JL,JK)*RCPD*(1.0_JPRB+RVTMP2*PQM1(JL,JK))
    IF (LDEBUGPRINT) THEN
      PRINT *,"DEBUGGING GEOPOTENTIAL"
      PRINT *,"PGEOM1(JL,JK),PTM1(JL,JK),RCPD,RVTMP2,PQM1(JL,JK)",PGEOM1(JL,JK),PTM1(JL,JK),RCPD,RVTMP2,PQM1(JL,JK)
      PRINT *,"ZCPTGZ(JL,JK)",ZCPTGZ(JL,JK)
    ENDIF
  ENDDO
ENDDO

!*         2.  Compute all surface related quantities
!          ------------------------------------------

!! OFFLINE only
LSICOUP=.FALSE.
!LNEMOLIMTHK=.FALSE.
ZSNTICE(:) = 0.0_JPRB
ZTHKICE(:) = 273._JPRB
PRPLRG=1._JPRB
PI10FGCV(:) = 0.0_JPRB 

!Use first tracer as atmospheric CO2 to be used in CTESSEL
ZATMCO2=PCM1(:,KLEV,1)

CALL SURFEXCDRIVER(YDSURF,CDCONF=CDCONF, &
 & KIDIA=KIDIA, KFDIA=KFDIA, KLON=KLON, KLEVS=KLEVS, KTILES=KTILES, KVTYPES=KVTYPES, KDIAG=KDIAG, &
 & KSTEP=KSTEP, KLEVSN=KLEVSN, KLEVI=KLEVI, KDHVTLS=KDHVTLS, KDHFTLS=KDHFTLS, &
 & KDHVTSS=KDHVTSS, KDHFTSS=KDHFTSS, KDHVTTS=KDHVTTS, KDHFTTS=KDHFTTS, &
 & KDHVTIS=KDHVTIS, KDHFTIS=KDHFTIS, K_VMASS=N_VMASS, &
 & KDHVCO2S=KDHVCO2S, KDHFCO2S=KDHFCO2S,KDHVVEGS=KDHVVEGS,KDHFVEGS=KDHFVEGS, &
 & PTSTEP=PTSPHY,PTSTEPF=PTSPHY, &
! input data, non-tiled
 & KTVL=KTVL, KCO2TYP=KCO2TYP, KTVH=KTVH, PCVL=PCVL, PCVH=PCVH, PCUR=PCUR, &
 & PLAIL=PLAIL, PLAIH=PLAIH, PFWET=PFWET, PLAT=PLAT,&
 & PSNM=PSNM1M, PRSN=PRSNM1M, &
 & PMU0=PMU0 ,PCARDI=RCARDI, &
 & PUMLEV=PUM1(:,KLEV), PVMLEV=PVM1(:,KLEV), PTMLEV=PTM1(:,KLEV), &
 & PQMLEV=PQM1(:,KLEV), PCMLEV=ZATMCO2, PAPHMS=PAPHM1(:,KLEV), PGEOMLEV=PGEOM1(:,KLEV), &
 & PCPTGZLEV=ZCPTGZ(:,KLEV), PSST=PSST, PTSKM1M=PTSKM1M, PCHAR=PCHAR, PCHARHQ=PCHARHQ, &
 & PSSRFL=PSSRFL, PSLRFL=PSLRFL, PEMIS=PEMIS, PTICE=PTICE, PTSN=PTSN, &
 & PHLICE=PHLICE,PTLICE=PTLICE,PTLWML=PTLWML, &
 & PTHKICE=ZTHKICE, PSNTICE=ZSNTICE, &
 & PWLMX=PWLMX, PUCURR=PUCURR, PVCURR=PVCURR,PI10FGCV=PI10FGCV, &
! input data, soil
 & PTSAM1M=PTSAM1M, PWSAM1M=PWSAM1M, KSOTY=KSOTY,&
! input data, tiled
 & PFRTI=PFRTI, PALBTI=PALBTI, &
! updated data, tiled
 & PUSTRTI=PUSTRTI, PVSTRTI=PVSTRTI, PAHFSTI=PAHFSTI, PEVAPTI=PEVAPTI, &
 & PTSKTI=PTSKTI, &
 & PANDAYVT=PANDAYVT,PANFMVT=PANFMVT,&
! updated data, non-tiled
 & PZ0M=PZ0M, PZ0H=PZ0H, &
! output data, tiled
 & PSSRFLTI=PSSRFLTI, PQSTI=ZQSTI, PDQSTI=ZDQSTI, PCPTSTI=ZCPTSTI, &
 & PCFHTI=ZCFHTI, PCFQTI=ZCFQTI, PCSATTI=ZCSATTI, PCAIRTI=ZCAIRTI, &
 & PCPTSTIU=ZCPTSTIU, PCSATTIU=ZCSATTIU, PCAIRTIU=ZCAIRTIU, & 
 & PRAQTI=ZRAQTI, PTSRF=ZTSRF,PLAMSK=ZLAMSK,&
! output data, non-tiled
 & PKHLEV=PKH(:,KLEV), PKCLEV=ZKCLEV, PCFMLEV=ZCFM(:,KLEV), PKMFL=ZKMFL, PKHFL=ZKHFL, &
 & PKQFL=ZKQFL, PEVAPSNW=PEVAPSNW, PZ0MW=ZZ0MW, PZ0HW=ZZ0HW, PZ0QW=ZZ0QW, &
 & PBLENDPP=ZBLEND, PCPTSPP=ZZCPTS, PQSAPP=ZZQSA, PBUOMPP=ZZBUOM, &
 & PZDLPP=ZZZDL, &
! output data, non-tiled CO2
 & PAN=PAN,PAG=PAG,PRD=PRD,PRSOIL_STR=PRSOIL_STR,PRECO=PRECO,PCO2FLUX=PCO2FLUX,PCH4FLUX=PCH4FLUX,&
 & PWETB=ZWETB, PWETL=ZWETL, PWETLU=ZWETLU, PWETH=ZWETH, PWETHS=ZWETHS, & 
! output data, diagnostics
 & PDHTLS=PDHTLS, PDHTSS=PDHTSS, PDHTTS=PDHTTS, PDHTIS=PDHTIS, &
 & PDHVEGS=PDHVEGS,PEXDIAG=PEXDIAG,PDHCO2S=PDHCO2S, &
 & PRPLRG=PRPLRG, &
! LIM switch
 & LSICOUP=LSICOUP &
! & LNEMOLIMTHK=LNEMOLIMTHK &
 & )

! Re-scale dussion coefficients for code clean-up in SURF
! Further clean-up still needed in OSM  (A. Beljaars)

DO JT=1,KTILES
  DO JL=KIDIA,KFDIA
    ZCFHTI(JL,JT)=ZCFHTI(JL,JT)*RG*PTSPHY*RVDIFTS
    ZCFQTI(JL,JT)=ZCFQTI(JL,JT)*RG*PTSPHY*RVDIFTS
  ENDDO
ENDDO
DO JL=KIDIA,KFDIA
  ZCFM(JL,KLEV)=ZCFM(JL,KLEV)*RG*PTSPHY*RVDIFTS
ENDDO

!     ------------------------------------------------------------------

!*         3.     NEW VARIABLES S, SLG, QT
!*                (at initial time level)
!                 ------------------------

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA

!*             total water and generalized liquid water static energy 
!*              slg = cp*T + gz - Lcond*ql - Ldep*qi

!1s    ZSLGM1(JL,JK) = ZCPTGZ(JL,JK) - RLVTT * PLM1(JL,JK) - RLSTT * PIM1(JL,JK)
!1s    ZSLGE(JL,JK)  = RCPD * ( ( 1.0_JPRB + RVTMP2 * PQM1(JL,JK) ) * PTE(JL,JK) &!dcpT/dt
!1s                & + RVTMP2 * PTM1(JL,JK)   * PQE(JL,JK) ) &                    !  -"-
!1s                & - RLVTT * PLE(JL,JK) - RLSTT * PIE(JL,JK)                    !dLqli/dt  
!1s    ZQTM1(JL,JK ) = PQM1(JL,JK) + PLM1(JL,JK) + PIM1(JL,JK)
!1s    ZQTE(JL,JK)   = PQE(JL,JK)  + PLE(JL,JK)  + PIE(JL,JK)             !dyn. qt tendency

!1s             slg = cp*T + gz
    ZSLGM1(JL,JK) = ZCPTGZ(JL,JK)
    ZSLGE(JL,JK)  = RCPD * ( ( 1.0_JPRB + RVTMP2 * PQM1(JL,JK) ) * PTE(JL,JK) &!dcpT/dt
                & + RVTMP2 * PTM1(JL,JK)   * PQE(JL,JK) )                      !  -"-
    ZQTM1(JL,JK ) = PQM1(JL,JK)             !dyn. qt tendency (no PLM1 PIM1)
    ZQTE(JL,JK)   = PQE(JL,JK)              !dyn. qt tendency (no PLE PIE in offline)
    ZSLGEA(JL,JK) = ZSLGE(JL,JK)
    ZQTEA(JL,JK)  = ZQTE(JL,JK)
  ENDDO
ENDDO


!     ------------------------------------------------------------------
 
!*         4.     CALCULATE QT VARIANCE AFTER DYN+RAD (for var. equ.)
!*                AND CONSISTENT T, QV, QL, QI AND CLOUD FRACTION
!*                (for tendency calculation at end) 
!                 ---------------------------------------------------

!*         4.0  state after dynamics and radiation

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
!1s    ZLUPD(JL,JK)  = PLM1(JL,JK) + PLE(JL,JK) * ZTMST
!1s    ZIUPD(JL,JK)  = PIM1(JL,JK) + PIE(JL,JK) * ZTMST
    ZQUPD(JL,JK)  = PQM1(JL,JK) + PQE(JL,JK) * ZTMST
    ZTUPD(JL,JK)  = PTM1(JL,JK) + PTE(JL,JK) * ZTMST
!          total condensate 
!     ZLIUPD(JL,JK) = ZLUPD(JL,JK) + ZIUPD(JL,JK)
  ENDDO
ENDDO


!*         4.1  qsat, dqsat/dT and alfa (= liquid fraction)
!*              (after dyn+rad, kept constant through vdfmain)

!          qsat

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
!1s    ZQSVAR(JL,JK) = FOEEWM(ZTUPD(JL,JK))/PAPM1(JL,JK)
!     ZQSVAR(JL,JK) = MIN(0.5_JPRB,ZQSVAR(JL,JK))
!     ZQSVAR(JL,JK) = ZQSVAR(JL,JK)/(1.0_JPRB-RETV*ZQSVAR(JL,JK))
!     ZQSVAR(JL,JK) = MAX(ZQSVAR(JL,JK),ZQUPD(JL,JK)) !don't allow input supersat.

!          dqsat/dT correction factor (1+L/cp*dqsat/dT) & alfa

    ZALFAW(JL,JK)=FOEALFA(ZTUPD(JL,JK))
    ZFACW=R5LES/((ZTUPD(JL,JK)-R4LES)**2)
    ZFACI=R5IES/((ZTUPD(JL,JK)-R4IES)**2)
    ZFAC=ZALFAW(JL,JK)*ZFACW+(1.0_JPRB-ZALFAW(JL,JK))*ZFACI
    ZESDP=0. ! Emanuel Dutra ZESDP must be defined to calculate next ZCOR 
    ZCOR=1.0_JPRB/(1.0_JPRB-RETV*ZESDP)
!     ZDQSDTEMP(JL,JK)=ZFAC*ZCOR*ZQSVAR(JL,JK)  !dqsat/dT
!     ZCORQS(JL,JK)=MAX(1.0_JPRB,1.0_JPRB+FOELDCPM(ZTUPD(JL,JK))*ZDQSDTEMP(JL,JK))
  ENDDO
ENDDO


!1s !*         4.2  diagnose total water variance
!1s !*              (qv, ql+qi, qsat -> variance;  cloud fraction ignored)
!1s 
!1s IF ( KCNT == 1 ) THEN  ! prognostic qt variance within vdf interations
!1s   CALL CLOUDVAR &
!1s !---input
!1s  & ( KIDIA, KFDIA, KLON  , KLEV  , 1    , KLEV, &
!1s  &   ZTUPD, ZQUPD, ZQSVAR, ZLIUPD, PAPM1, &
!1s !---output
!1s  &   PVAR , ZANEW, ZLNEW ) !last two are dummy args  
!1s ENDIF


!1s !*         4.3  diagnose cloud fraction and T, qv, ql, qi
!1s !*              (T, qv, ql qi should be approximately unchanged;
!1s !*              zdelq ~ 0)  
!1s 
!1s !          From total water and its variance calculate 
!1s !          cloud cover (zanew) and conversion of qv to qc (zdelq).
!1s 
!1s CALL COVER &
!1s !---input
!1s  & ( KIDIA, KFDIA , KLON  , KLEV , 1   , KLEV, &
!1s  &   ZQUPD, ZQSVAR, ZLIUPD, PAPM1, PVAR, &
!1s !---output
!1s  &   ZLNEW, ZANEW )  

!          Add modifications to estimate of initial state.

! DO JK=1,KLEV
!   DO JL=KIDIA,KFDIA
!     ZDELQ        = ( ZLNEW(JL,JK) - ZLIUPD(JL,JK) ) / ZCORQS(JL,JK)
!     ZQINI(JL,JK) = ZQUPD(JL,JK) -                            ZDELQ
!     ZLINI(JL,JK) = ZLUPD(JL,JK) +       ZALFAW(JL,JK)      * ZDELQ
!     ZIINI(JL,JK) = ZIUPD(JL,JK) + ( 1.0_JPRB - ZALFAW(JL,JK)) * ZDELQ
!     ZTINI(JL,JK) = ZTUPD(JL,JK) + FOELDCPM(ZTUPD(JL,JK))   * ZDELQ
!     ZAINI(JL,JK) = ZANEW(JL,JK)
!   ENDDO
! ENDDO

!     ------------------------------------------------------------------

!*         5.     EXCHANGE COEFFICIENTS
!                 ---------------------

!1s !*         5.4  COMPUTATION OF THE PBL EXTENSION
!1s 
!1s !          SET PBL HEIGHT-INDEX TO 1
!1s 
!1s DO JL=KIDIA,KFDIA
!1s   KHPBLN(JL)=1
!1s ENDDO
ITOP=1

!          FLUX BOUNDARY CONDITION
IF (LLSFCFLX) THEN
  DO JL=KIDIA,KFDIA
    ZRHO = PAPHM1(JL,KLEV)/( RD*PTM1(JL,KLEV)*(1.0_JPRB+RETV*PQM1(JL,KLEV)) )
    ZKHFL(JL) = ZFSH1D(JL) / ( RCPD*(1.0_JPRB+RVTMP2*PQM1(JL,KLEV)) ) / ZRHO
    ZKQFL(JL) = ZFLH1D(JL) / RLVTT / ZRHO
  ENDDO
ENDIF


!1s !*         5.5  BOUNDARY LAYER HEIGHT FOR DIANOSTICS ONLY
!1s 
!1s CALL VDFDPBL(KIDIA,KFDIA,KLON,KLEV,&
!1s  & PUM1,PVM1,PTM1,PQM1,PGEOM1,&
!1s  & ZKMFL,ZKHFL,ZKQFL,PBLH)  
!1s 
!1s 
!1s !*         5.6  PARCEL UPDRAFT
!1s 
!1s !          iterate inversion height (PZINV)
!1s  
!1s CALL VDFHGHTN (KIDIA   , KFDIA   , KLON    , KLEV    , KHPBLN   , ZTMST,&
!1s              & PTM1    , PQM1    , PLM1    , PIM1    , PAM1,&
!1s              & PAPHM1  , PAPM1   , PGEOM1  , PGEOH,&
!1s              & ZKMFL   , ZKHFL   , ZKQFL   , ZMFLX,&
!1s              & ZSLGUH  , ZQTUH   , PZINV   , PWUAVG  , ZCLDBASE,&
!1s              & PBIR    , LDNODECP, LLRUNDRY, KPBLTYPE)  
!1s Assign PZINV to dummy value (to comply with intent OUT)
PZINV=0.0_JPRB
!1s 
!1s 
!1s !*         5.7  EXCHANGE COEFFICIENTS ABOVE THE SURFACE LAYER
!1s 
!1s CALL VDFEXCU(KIDIA  , KFDIA  , KLON   , KLEV   , ZTMST  , PZ0M   , &
!1s            & PHRLW  , PHRSW  , PUM1   , PVM1   , PTM1   , PQM1   , &
!1s            & PAPHM1 , PAPM1  , PGEOM1 , PGEOH  , ZCPTGZ , &
!1s            & ZKMFL  , ZKHFL  , ZKQFL  , ZCFM   , ZCFH   , &
!1s            & PZINV  , KHPBLN , PKH    , ZCLDBASE        , KPBLTYPE)  
!1s 
!1s 
!1s !*         5.8  MASS FLUX MODIFICATIONS
!1s 
!1s !          Selective mflux=0 for single mflux layers
!1s DO JL=KIDIA,KFDIA
!1s   IF ( ZMFLX(JL,KLEV-2) < 1.E-40 ) THEN
!1s     ZMFLX(JL,KLEV-1) = 0.0
!1s   ENDIF
!1s ENDDO
!1s !          Turn off mass-flux below convection
!1s DO JL=KIDIA,KFDIA
!1s   IF ( KPBLTYPE(JL) == 3 ) THEN
!1s     ZMFLX(JL,:)  = 0.0_JPRB
!1s     ZSLGUH(JL,:) = 0.0_JPRB
!1s     ZQTUH (JL,:) = 0.0_JPRB
!1s   ENDIF
!1s ENDDO


!     ------------------------------------------------------------------

!*         6.     SOLVE ADVECTION-DIFFUSION EQUATION
!                 ----------------------------------

!*         6.1  MOMENTUM
! 

CALL VDFDIFM1S (KIDIA, KFDIA, KLON , KLEV  , ITOP, &
            & ZTMST, PUM1 , PVM1 , PAPHM1, ZCFM, &
            & PVOM , PVOL , ZUDIF, ZVDIF)  


!*         6.2  GENERALIZED LIQUID WATER STATIC ENERGY AND TOTAL WATER
IF (LDEBUGPRINT) THEN
PRINT *,'DEBUGGING VDFDIFH1S --->'
PRINT *,'KIDIA',KIDIA  
PRINT *,'KFDIA', KFDIA  
PRINT *,'KLON', KLON   
PRINT *,'KLEV', KLEV   
PRINT *,'ITOP', ITOP   
PRINT *,'KTILES', KTILES
PRINT *,'ZTMST', ZTMST  
PRINT *,'ZFSH1D', ZFSH1D 
PRINT *,'ZFLH1D', ZFLH1D 
PRINT *,'LLSFCFLX', LLSFCFLX
PRINT *,'KSOTY', KSOTY  
PRINT *,'PFRTI', PFRTI  
PRINT *,'PSSRFLTI', PSSRFLTI
PRINT *,'PSLRFL',PSLRFL 
PRINT *,'PEMIS', PEMIS  
PRINT *,'PEVAPSNW', PEVAPSNW
PRINT *,'ZSLGM1', ZSLGM1 
PRINT *,'PTM1', PTM1   
PRINT *,'PQM1', PQM1   
PRINT *,'ZQTM1', ZQTM1  
PRINT *,'PAPHM1', PAPHM1 
! PRINT *,'ZCFH', ZCFH   
PRINT *,'ZCFHTI', ZCFHTI 
PRINT *,'ZCFQTI', ZCFQTI 
PRINT *,'ZMFLX', ZMFLX  
PRINT *,'ZSLGUH', ZSLGUH 
PRINT *,'ZQTUH', ZQTUH  
! PRINT *,'ZSLGDIF',ZSLGDIF
! PRINT *,'ZQTDIF', ZQTDIF 
PRINT *,'ZCPTSTI', ZCPTSTI
PRINT *,'ZQSTI', ZQSTI  
PRINT *,'ZCAIRTI', ZCAIRTI
PRINT *,'ZCSATTI', ZCSATTI
PRINT *,'ZDQSTI', ZDQSTI 
PRINT *,'PTSKTI', PTSKTI 
PRINT *,'PTSKRAD', PTSKRAD
PRINT *,'PTSAM1M(1,1)', PTSAM1M(1,1)    
PRINT *,'PTSN', PTSN 
PRINT *,'PTICE', PTICE  
PRINT *,'PSST', PSST
PRINT *,'ZTSKTIP1', ZTSKTIP1
PRINT *,'ZSLGE',ZSLGE  
PRINT *,'PTE', PTE    
PRINT *,'ZQTE', ZQTE
PRINT *,'PEVAPTI', PEVAPTI
PRINT *,'PAHFSTI', PAHFSTI
PRINT *,'ZAHFLTI', ZAHFLTI
! PRINT *,'ZSTR', PSLRFLTI  
PRINT *,'ZG0', ZG0
ENDIF
CALL VDFDIFH1S (KIDIA  , KFDIA  , KLON   , KLEV   , ITOP   , KTILES, KTVL, KTVH, KLEVSN, &
            & ZTMST  , ZFSH1D , ZFLH1D , LLSFCFLX, &
            & PFRTI  , PSSRFLTI,PSLRFL , PEMIS  , PEVAPSNW, &
            & PHLICE , PTLICE , PTLWML , &
            & ZSLGM1 , PTM1   , PQM1   , ZQTM1  , PAPHM1 , &
            & ZCFH   , ZCFHTI , ZCFQTI , ZMFLX  , ZSLGUH , ZQTUH  , &
            & ZSLGDIF, ZQTDIF , ZCPTSTI, ZQSTI  , ZCAIRTI, ZCSATTI, &
            & ZCPTSTIU,ZCAIRTIU,ZCSATTIU,ZTSRF,ZLAMSK,PSNM1M(1,1),PRSNM1M(1,1), &
            & ZDQSTI , PTSKTI , PTSKRAD, PTSAM1M(1,1)    , &
            & PTSN , PTICE  , PSST, &
            & ZTSKTIP1,ZSLGE  , PTE    , ZQTE, &
            & PEVAPTI, PAHFSTI, ZAHFLTI, PSLRFLTI   , ZG0,PEVAPTIU)

! ! Update Tiled LW components: what is done in surfexcdriver does not acconnts for the LWtiling ....
! DO JT=1,KTILES
!   ! upwad 
!   PDHTLS(KIDIA:KFDIA,JT,7)=PFRTI(KIDIA:KFDIA,JT)*&
!    & PEMIS(KIDIA:KFDIA)*RSIGMA*PTSKTI(KIDIA:KFDIA,JT)**4
!   ! lwdown
!   PDHTLS(KIDIA:KFDIA,JT,6)=PFRTI(KIDIA:KFDIA,JT)*&
!    & PSLRFLTI(KIDIA:KFDIA,JT)-PDHTLS(KIDIA:KFDIA,JT,7)
! ENDDO

IF (LDEBUGPRINT) THEN
PRINT *,'DEBUGGING VDFDIFH1S <---'
PRINT *,'KIDIA',KIDIA  
PRINT *,'KFDIA', KFDIA  
PRINT *,'KLON', KLON   
PRINT *,'KLEV', KLEV   
PRINT *,'ITOP', ITOP   
PRINT *,'KTILES', KTILES
PRINT *,'ZTMST', ZTMST  
PRINT *,'ZFSH1D', ZFSH1D 
PRINT *,'ZFLH1D', ZFLH1D 
PRINT *,'LLSFCFLX', LLSFCFLX
PRINT *,'KSOTY', KSOTY  
PRINT *,'PFRTI', PFRTI  
PRINT *,'PSSRFLTI', PSSRFLTI
PRINT *,'PSLRFL',PSLRFL 
PRINT *,'PEMIS', PEMIS  
PRINT *,'PEVAPSNW', PEVAPSNW
PRINT *,'ZSLGM1', ZSLGM1 
PRINT *,'PTM1', PTM1   
PRINT *,'PQM1', PQM1   
PRINT *,'ZQTM1', ZQTM1  
PRINT *,'PAPHM1', PAPHM1 
! PRINT *,'ZCFH', ZCFH   
PRINT *,'ZCFHTI', ZCFHTI 
PRINT *,'ZCFQTI', ZCFQTI 
PRINT *,'ZMFLX', ZMFLX  
PRINT *,'ZSLGUH', ZSLGUH 
PRINT *,'ZQTUH', ZQTUH  
PRINT *,'ZSLGDIF',ZSLGDIF
PRINT *,'ZQTDIF', ZQTDIF 
PRINT *,'ZCPTSTI', ZCPTSTI
PRINT *,'ZQSTI', ZQSTI  
PRINT *,'ZCAIRTI', ZCAIRTI
PRINT *,'ZCSATTI', ZCSATTI
PRINT *,'ZDQSTI', ZDQSTI 
PRINT *,'PTSKTI', PTSKTI 
PRINT *,'PTSKRAD', PTSKRAD
PRINT *,'PTSAM1M(1,1)', PTSAM1M(1,1)    
PRINT *,'PTSN', PTSN 
PRINT *,'PTICE', PTICE  
PRINT *,'PSST', PSST
PRINT *,'ZTSKTIP1', ZTSKTIP1
PRINT *,'ZSLGE',ZSLGE  
PRINT *,'PTE', PTE    
PRINT *,'ZQTE', ZQTE
PRINT *,'PEVAPTI', PEVAPTI
PRINT *,'PAHFSTI', PAHFSTI
PRINT *,'ZAHFLTI', ZAHFLTI
PRINT *,'ZSTR', PSLRFLTI
PRINT *,'ZG0', ZG0
ENDIF
 

!*         6.3  INCREMENTATION OF U AND V TENDENCIES, STORAGE OF
!*              THE DISSIPATION, COMPUTATION OF MULTILEVEL FLUXES.

CALL VDFINCR (KIDIA  , KFDIA  , KLON   , KLEV   , ITOP   , ZTMST  , &
            & PUM1   , PVM1   , ZSLGM1 , PTM1   , ZQTM1  , PAPHM1 , PGEOM1 , &
            & ZCFM   , ZUDIF  , ZVDIF  , ZSLGDIF, ZQTDIF , &
            & PVOM   , PVOL   , ZSLGE  , ZQTE   , ZSLGEWODIS, &
            & PVDIS  , PSTRTU , PSTRTV)  


!          6.4  Solve for tracers

IF (LVDFTRAC .AND. KTRAC > 0) THEN 
  CALL VDFDIFC(KIDIA,KFDIA,KLON,KLEV,ITOP,KTRAC,&
             & ZTMST,PCM1,PTENC,PAPHM1,ZCFH,PCFLX)
ENDIF


!     ------------------------------------------------------------------

!*         7.     SURFACE FLUXES - TILES
!                 ----------------------
!*         AND    COMPUTE 2M TEMPERATURE AND HUMIDITY, 10M WIND,
!*                  and gustiness

!1s !  Compute wind speed at blending height
!1s 
!1s CALL VDFFBLEND(KIDIA,KFDIA,KLON,KLEV, &
!1s  & PUM1, PVM1, PGEOM1, ZBLEND, &
!1s  & ZFBLEND)
!1s The wind speed is simply 
DO JL=KIDIA,KFDIA
  ZFBLEND(JL)=SQRT(PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2)
  ZBLEND(JL)=RALT
ENDDO

! Wrap-up computations for the surface and 2T/2D/10U/10V/gustiness computation

CALL SURFPP( YDSURF,KIDIA=KIDIA,KFDIA=KFDIA,KLON=KLON,KTILES=KTILES, &
 & KDHVTLS=KDHVTLS,KDHFTLS=KDHFTLS, &
 & PTSTEP=PTSPHY, &
! input
 & PFRTI=PFRTI, PAHFLTI=ZAHFLTI, PG0TI=ZG0, &
 & PSTRTULEV=PSTRTU(:,KLEV), PSTRTVLEV=PSTRTV(:,KLEV), PTSKM1M=PTSKM1M, &
 & PUMLEV=PUM1(:,KLEV), PVMLEV=PVM1(:,KLEV), PQMLEV=PQM1(:,KLEV), &
 & PGEOMLEV=PGEOM1(:,KLEV), PCPTSPP=ZZCPTS, PCPTGZLEV=ZCPTGZ(:,KLEV), &
 & PAPHMS=PAPHM1(:,KLEV), PZ0MW=ZZ0MW, PZ0HW=ZZ0HW, PZ0QW=ZZ0QW, &
 & PZDL=ZZZDL, PQSAPP=ZZQSA, PBLEND=ZBLEND, PFBLEND=ZFBLEND, PBUOM=ZZBUOM, &
 & PZ0M=PZ0M, PEVAPSNW=PEVAPSNW,PSSRFLTI=PSSRFLTI, PSLRFL=PSLRFL, PSST=PSST, &
 & PUCURR=PUCURR, PVCURR=PVCURR, PUSTOKES=ZUSTOKES, PVSTOKES=ZVSTOKES, &
! updated
 & PAHFSTI=PAHFSTI, PEVAPTI=PEVAPTI, PTSKE1=PTSKE1,PTSKTIP1=ZTSKTIP1, &
! output
 & PDIFTSLEV=PDIFTS(:,KLEV), PDIFTQLEV=PDIFTQ(:,KLEV), PUSTRTI=PUSTRTI, &
 & PVSTRTI=PVSTRTI,  PTSKTI=PTSKTI, PAHFLEV=PAHFLEV, PAHFLSB=PAHFLSB, &
 & PFWSB=PFWSB, PU10M=PU10M, PV10M=PV10M, PT2M=PT2M, PD2M=PD2M, PQ2M=PQ2M, &
 & PGUST=PGUST, P10NU=Z10NU, P10NV=Z10NV, PUST=ZUST, &
! output DDH
 & PDHTLS=PDHTLS, &
 & PRPLRG=PRPLRG &
 & )

!1s PDIFTL  (KIDIA:KFDIA,KLEV) = 0.0_JPRB
!1s PDIFTI  (KIDIA:KFDIA,KLEV) = 0.0_JPRB
ZDIFTQT (KIDIA:KFDIA,KLEV) = 0.0_JPRB
ZDIFTSLG(KIDIA:KFDIA,KLEV) = 0.0_JPRB
ZDIFTQT (KIDIA:KFDIA,KLEV) = PDIFTQ(KIDIA:KFDIA,KLEV)
ZDIFTSLG(KIDIA:KFDIA,KLEV) = PDIFTS(KIDIA:KFDIA,KLEV)

!     ------------------------------------------------------------------

!*         9.     SLG, QT, U, V FLUX COMPUTATIONS AND T,SKIN TENDENCY
!                 ---------------------------------------------------

DO JL=KIDIA,KFDIA
  ZDIFTQT (JL,0) = 0.0_JPRB
  PDIFTQ  (JL,0) = 0.0_JPRB
!1s   PDIFTL  (JL,0) = 0.0_JPRB
!1s   PDIFTI  (JL,0) = 0.0_JPRB
  PDIFTS  (JL,0) = 0.0_JPRB
  ZDIFTSLG(JL,0) = 0.0_JPRB
  PSTRTU  (JL,0) = 0.0_JPRB
  PSTRTV  (JL,0) = 0.0_JPRB
ENDDO

DO JK=KLEV-1,1,-1
  DO JL=KIDIA,KFDIA
    ZGDPH = - (PAPHM1(JL,JK)-PAPHM1(JL,JK+1)) * ZRG
!...change in slg,qt,u,v tendencies are converted to fluxes
    ZDIFTSLG(JL,JK) = ( ZSLGEWODIS(JL,JK+1) - ZSLGEA(JL,JK+1) ) * ZGDPH &
                    & + ZDIFTSLG(JL,JK+1)  
    ZDIFTQT(JL,JK)  = (ZQTE (JL,JK+1)-ZQTEA(JL,JK+1))*ZGDPH + ZDIFTQT(JL,JK+1)
    PSTRTU(JL,JK)   = (PVOM(JL,JK+1) -ZUEA(JL,JK+1)) *ZGDPH + PSTRTU(JL,JK+1)
    PSTRTV(JL,JK)   = (PVOL(JL,JK+1) -ZVEA(JL,JK+1)) *ZGDPH + PSTRTV(JL,JK+1)
  ENDDO
ENDDO


!     ------------------------------------------------------------------

!*         10.    VARIANCE EVOLUTION EQUATION
!                 ---------------------------

!          use vdf updated variables for cloud variance

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZQTUPD(JL,JK)  = ZQTM1(JL,JK)  + ZQTE(JL,JK)  * ZTMST
    ZSLGUPD(JL,JK) = ZSLGM1(JL,JK) + ZSLGE(JL,JK) * ZTMST
  ENDDO
ENDDO

!1s !          turbulence source and sink terms
!1s !          solve d(sigma^2)/dt = - 2 w'q' dq/dz - sigma^2/tau - d(w'sigma^2')/dz
!1s !          (analytically)
!1s 
!1s DO JL=KIDIA,KFDIA
!1s   ZTAU(JL) = PZINV(JL) / MAX(PWUAVG(JL),0.01_JPRB)
!1s   ZTAU(JL) = MAX(ZTAU(JL), 100.0_JPRB)  ! prevent tau=0 for stable boundary layer (zi=0)
!1s ENDDO
!1s 
!1s DO JK=1,KLEV
!1s   DO JL=KIDIA,KFDIA
!1s     IF ( JK >= KHPBLN(JL) ) THEN        ! only within PBL - for now
!1s       IF ( JK < KLEV ) THEN             ! var. gener. not for KLEV - upstream impossible
!1s !...centered differencing:
!1s !       ZVARGEN= ZDIFTQT(JL,JK)   / PAPHM1(JL,JK)   * RD * (PTM1(JL,JK)+PTM1(JL,JK+1)) / 2.0 &
!1s !            & * (ZQTUPD(JL,JK)-ZQTUPD(JL,JK+1)) / (PGEOM1(JL,JK)-PGEOM1(JL,JK+1)) * RG      &
!1s !            & + ZDIFTQT(JL,JK-1) / PAPHM1(JL,JK-1) * RD * (PTM1(JL,JK-1)+PTM1(JL,JK)) / 2.0 &
!1s !            & * (ZQTUPD(JL,JK-1)-ZQTUPD(JL,JK)) / (PGEOM1(JL,JK-1)-PGEOM1(JL,JK)) * RG
!1s !...upstream differencing:
!1s         ZVARGEN     = ZDIFTQT(JL,JK)                                   &! - rho*w'qt'
!1s                   & / PAPHM1(JL,JK) * RD * (PTM1(JL,JK)+PTM1(JL,JK+1)) &! * 2 / rho
!1s                   & * (ZQTUPD(JL,JK)-ZQTUPD(JL,JK+1))                  &! * dqt/dz
!1s                   & / (PGEOM1(JL,JK)-PGEOM1(JL,JK+1)) * RG  
!1s         ZVARGEN     = MAX(ZVARGEN,0.0_JPRB)            ! exclude countergradient flow
!1s         PVAR(JL,JK) = ZVARGEN * ZTAU(JL) + (PVAR(JL,JK)- ZVARGEN * ZTAU(JL) ) &
!1s                   & * EXP( - ZTMST / ZTAU(JL) )  
!1s       ELSE
!1s         PVAR(JL,JK) = PVAR(JL,JK) * EXP( -ZTMST / ZTAU(JL) )    ! decay only for KLEV
!1s       ENDIF
!1s       PVAR(JL,JK) = MAX(PVAR(JL,JK),0.0_JPRB)
!1s     ENDIF
!1s   ENDDO
!1s ENDDO


!     ------------------------------------------------------------------

!*         11.    CONVERT VARIANCE -> QV, QL, QI AND CLOUD FRACTION
!*                (on final state; only for PBL, above PBL tendencies
!*                will be overwritten)
!                 ---------------------------------------------------

!          Guess an (arbitrary) state consistent with final qt and slg.
!          (choice: Use ql and qi state after dyn+rad.  Total water 
!          conservation then dictates qv and energy conservation 
!          dictates T.  Z*UPD values are then estimates, before cover,
!          of final profiles.)

! DO JK=1,KLEV
!   DO JL=KIDIA,KFDIA
! !1s     ZLUPD(JL,JK)  = PLM1(JL,JK) + PLE(JL,JK) * ZTMST   ! arbitrary
! !1s     ZIUPD(JL,JK)  = PIM1(JL,JK) + PIE(JL,JK) * ZTMST   ! choice
! !     ZQUPD(JL,JK)  = ZQTUPD(JL,JK) - ZLUPD(JL,JK) - ZIUPD(JL,JK)
! !     ZTUPD(JL,JK)  = ( ZSLGUPD(JL,JK) - PGEOM1(JL,JK) &
! !                 & + RLVTT * ZLUPD(JL,JK) + RLSTT * ZIUPD(JL,JK) &
! !                 & ) / ( RCPD * ( 1.0_JPRB + RVTMP2 * ZQUPD(JL,JK) ) )  
! !     ZLIUPD(JL,JK) = ZLUPD(JL,JK) + ZIUPD(JL,JK) ! total condensate
!   ENDDO
! ENDDO

!1s !          From total water and its variance calculate 
!1s !          cloud cover (ZANEW) and conversion of qv to qc (ZDELQ).
!1s !          (This produces modifications to qv, ql, qi and T. Cloud
!1s !          fraction is diagnosed.)
!1s 
!1s CALL COVER &
!1s !---input
!1s  & ( KIDIA, KFDIA , KLON  , KLEV , 1   , KLEV, &
!1s  &   ZQUPD, ZQSVAR, ZLIUPD, PAPM1, PVAR, &
!1s !---output
!1s  &   ZLNEW, ZANEW )

! DO JK=1,KLEV
!   DO JL=KIDIA,KFDIA
!     
! !          add modifications to estimate of final state
!   
! !     ZDELQ        = ( ZLNEW(JL,JK) - ZLIUPD(JL,JK) ) / ZCORQS(JL,JK)
! !     ZQUPD(JL,JK) = ZQUPD(JL,JK) -                            ZDELQ
! !     ZLUPD(JL,JK) = ZLUPD(JL,JK) +       ZALFAW(JL,JK)      * ZDELQ
! !     ZIUPD(JL,JK) = ZIUPD(JL,JK) + ( 1.0_JPRB - ZALFAW(JL,JK)) * ZDELQ
! !     ZTUPD(JL,JK) = ZTUPD(JL,JK) + FOELDCPM(ZTUPD(JL,JK))   * ZDELQ
! !     ZAUPD(JL,JK) = ZANEW(JL,JK)
! 
! !          vdf tendencies from final and initial state, where the
! !          initial state is taken consistent with the beta distribution
! 
! !     PQE(JL,JK) = ( ZQUPD(JL,JK) - ZQINI(JL,JK) ) * ZRTMST + PQE(JL,JK)
! ! !1s     PLE(JL,JK) = ( ZLUPD(JL,JK) - ZLINI(JL,JK) ) * ZRTMST + PLE(JL,JK)
! ! !1s     PIE(JL,JK) = ( ZIUPD(JL,JK) - ZIINI(JL,JK) ) * ZRTMST + PIE(JL,JK)
! !     PTE(JL,JK) = ( ZTUPD(JL,JK) - ZTINI(JL,JK) ) * ZRTMST + PTE(JL,JK)
! !1s     PAE(JL,JK) = ( ZAUPD(JL,JK) - ZAINI(JL,JK) ) * ZRTMST + PAE(JL,JK)
! 
!   ENDDO
! ENDDO

!     ------------------------------------------------------------------

!*         12.    Q, QL, QI AND S FLUX COMPUTATIONS
!                 ---------------------------------
 
DO JK=KLEV-1,1,-1
  DO JL=KIDIA,KFDIA
    ZGDPH = - (PAPHM1(JL,JK)-PAPHM1(JL,JK+1)) * ZRG
!...changes in q,l,i tendencies are converted to fluxes
    PDIFTQ(JL,JK) = (PQE (JL,JK+1) - ZQEA(JL,JK+1)) * ZGDPH + PDIFTQ(JL,JK+1)
!1s   slg=s (for offline)
    PDIFTS(JL,JK) = ZDIFTSLG(JL,JK) 
!1s   PDIFTL(JL,JK) = (PLE (JL,JK+1) - ZLEA(JL,JK+1)) * ZGDPH + PDIFTL(JL,JK+1)
!1s   PDIFTI(JL,JK) = (PIE (JL,JK+1) - ZIEA(JL,JK+1)) * ZGDPH + PDIFTI(JL,JK+1)
!1s !...slg=s-Lc*ql-Ld*qi (same for fluxes)
!1s    PDIFTS(JL,JK) = ZDIFTSLG(JL,JK) 
!1s                & + RLVTT * PDIFTL(JL,JK) + RLSTT * PDIFTI(JL,JK)  
  ENDDO
ENDDO


!     ------------------------------------------------------------------

!1s !*         13.    TKE buoyancy production/consumption decoupling criteria
!1s !                 -------------------------------------------------------
!1s 
!1s !          buoyancy flux: rho*w'sv' = alfa1 * rho*w'slg'  + alfa2 * L * rho*w'qt'
!1s !          (Martin Koehler, used here)
!1s !          (half level!; convert downward to upward flux; not ice generalized)
!1s !          attention: Stull 1988, p.551 misreads Moeng & Randall 1984 who
!1s !                     originally use slgv: 
!1s !          buoyancy flux: rho*w'sv' = alfa1 * rho*w'slgv' + alfa2 * rho*w'qt' 
!1s 
!1s IF ( .FALSE. ) THEN
!1s 
!1s   ZSVFLUXCLD(:,:) = 0.0_JPRB      !initialize (attention PBL top)
!1s   ZSVFLUXSUB(:,:) = 0.0_JPRB
!1s   ZSVFLUX   (:,:) = 0.0_JPRB
!1s   DO JK=1,KLEV
!1s     DO JL=KIDIA,KFDIA
!1s 
!1s !          cloud
!1s 
!1s       ZALFA1 = ( 1.0_JPRB + (1.0_JPRB+RETV) * ZQSVAR(JL,JK) - ZQTM1(JL,JK)   &
!1s        & + (1.0_JPRB+RETV) * PTM1(JL,JK) * ZDQSDTEMP(JL,JK)      ) &
!1s        & / ( 1.0_JPRB + RLVTT/RCPD * ZDQSDTEMP(JL,JK) )  
!1s       ZALFA2 = RLVTT * ZALFA1 - RCPD * PTM1(JL,JK)
!1s       ZSVFLUXCLD(JL,JK) = - ZALFA1 * ZDIFTSLG(JL,JK) - ZALFA2 * ZDIFTQT(JL,JK)
!1s 
!1s !          sub-cloud
!1s 
!1s       ZALFA1 = 1.0_JPRB + RETV *  ZQTM1(JL,JK)
!1s       ZALFA2 = RETV * RCPD * PTM1(JL,JK)
!1s       ZSVFLUXSUB(JL,JK) = - ZALFA1 * ZDIFTSLG(JL,JK) - ZALFA2 * ZDIFTQT(JL,JK)
!1s 
!1s       IF ( PGEOH(JL,JK)*ZRG < PZINV(JL) ) THEN     ! only define sv-flux within PBL
!1s         IF ( PGEOH(JL,JK)*ZRG > ZCLDBASE(JL) .AND. &
!1s            & KPBLTYPE(JL) == 2 ) THEN              ! to prevent cldbase=-100  
!1s           ZSVFLUX(JL,JK) = ZSVFLUXCLD(JL,JK)
!1s         ELSE
!1s           ZSVFLUX(JL,JK) = ZSVFLUXSUB(JL,JK)
!1s         ENDIF
!1s       ENDIF
!1s 
!1s     ENDDO
!1s   ENDDO
!1s 
!1s !          buoyancy flux integral (pos/neg): rho * w'sv' * dz integral
!1s !          ... simple constant flux layers between two full levels
!1s !          (future: extend to linear assumption between half levels, 
!1s !                   continuous zi & zcldbase)
!1s 
!1s   DO JK=1,KLEV-1
!1s     DO JL=KIDIA,KFDIA
!1s       ZDZH(JL,JK) = ( PGEOM1(JL,JK) - PGEOM1(JL,JK+1) ) * ZRG ! layer thickness
!1s     ENDDO
!1s   ENDDO
!1s   DO JL=KIDIA,KFDIA
!1s     ZDZH(JL,KLEV) = PGEOM1(JL,KLEV) * ZRG
!1s     ZDZH(JL,0)    = 0.0_JPRB                                  ! approximation
!1s 
!1s !          pos & neg flux within PBL and it's ratio (BIR)
!1s 
!1s     ZBUOYPOS(JL) = SUM( ZSVFLUX(JL,KHPBLN(JL):KLEV)*ZDZH(JL,KHPBLN(JL):KLEV), &
!1s      & MASK = ZSVFLUX(JL,KHPBLN(JL):KLEV) > 0.0_JPRB )  
!1s     ZBUOYNEG(JL) = SUM( ZSVFLUX(JL,KHPBLN(JL):KLEV)*ZDZH(JL,KHPBLN(JL):KLEV), &
!1s      & MASK = ZSVFLUX(JL,KHPBLN(JL):KLEV) < 0.0_JPRB )  
!1s     IF ( LDNODECP(JL) ) THEN          ! test vdfmain to get decoupling criteria
!1s       PBIR(JL)   = - ZBUOYNEG(JL) / MAX(ZBUOYPOS(JL), 1e-10_JPRB)
!1s       PBIR(JL)   = MIN(PBIR(JL), 10.0_JPRB)
!1s     ENDIF
!1s   ENDDO
!1s 
!1s ENDIF

IF (LHOOK) CALL DR_HOOK('VDFMAIN1S',1,ZHOOK_HANDLE)
END SUBROUTINE VDFMAIN1S
