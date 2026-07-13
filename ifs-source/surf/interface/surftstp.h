INTERFACE
SUBROUTINE SURFTSTP    (YDSURF, KIDIA , KFDIA , KLON  , KLEVS , KTILES,&
 & KLEVSN, KLEVO , KLEVI ,  KSTART , KSTEP,&
 & KDHVTSS , KDHFTSS,&
 & KDHVTTS , KDHFTTS,&
 & KDHVTIS , KDHFTIS,&
 & KDHVSSS , KDHFSSS,&
 & KDHVIIS , KDHFIIS,&
 & KDHVWLS , KDHFWLS,&
 & KDHVBIOS, KDHFBIOS, KDHVVEGS, KDHFVEGS,&                  !CTESSEL
 & KTVL , KTVH , KVEG , KSOTY,&
 & PTSPHY , PSDOR , PFRTI,&
 & PSST, PUSTRTI , PVSTRTI,&
 & PAHFSTI, PEVAPTI, PSSRFLTI,PSLRFLTI,&
 & PLAT , PANFMVT , PANDAYVT , PMU0,&                        !CTESSEL
 & PCVT, PLAIVT, PLAIL, PLAIH,&                              !CTESSEL
 & PLAILC,  PLAIHC,&
 & LDLAND, LDSICE, LDSI, LDNH, LDOCN_KPP,&
 & PSNM1M  ,PTSNM1M,PASNM1M,PRSNM1M,PWSNM1M,&
 & PAPRS,  PTSKM1M ,PTSAM1M,PTIAM1M,&
 & PWLM1M  ,PWSAM1M,&
 & PTLICEM1M,PTLMNWM1M,PTLWMLM1M,PTLBOTM1M,PTLSFM1M,& 
 & PHLICEM1M,PHLMLM1M,PGEMU,PLDEPTH,LDLAKE,PCLAKE, &
 & PRSFC   ,PRSFL,& 
 & PSLRFL  ,PSSFC   ,PSSFL,& 
 & PCVL    ,PCVH    ,PCUR    ,PWLMX   ,PEVAPSNW,&
 & PUSRF   ,PVSRF   ,PTSRF,& 
 & PZO     ,PHO     ,PHO_INV ,PDO     ,POCDEPTH ,&
 & PUO0    ,PVO0    ,PUOC    ,PVOC    ,PTO0     ,&
 & PSO0    ,PADVT   ,PADVS   ,PTRI0   ,PTRI1    ,&
 & PSWDK_SAVE, PUSTRC ,PVSTRC,&
 & PUSTOKES,PVSTOKES,PTAUOCX ,PTAUOCY ,PPHIOC   ,&
 & PWSEMEAN,PWSFMEAN  ,&
!-DIAGNOSTICS OUTPUT
 & PTSDFL  , PROFD , PROFS,&
 & PWFSD   , PMELT , PFWEV, PENES,&
 & PDIFM   , PDIFT , PDIFS, POTKE,&
 & PRESPBSTR,PRESPBSTR2,PBIOMASS_LAST,&                      !CTESSEL
 & PBIOMASSTR_LAST,PBIOMASSTR2_LAST,&                        !CTESSEL
 & PBLOSSVT, PBGAINVT, &                                     !CTESSEL
 & PLAI    , PBIOM , PBLOSS, PBGAIN, PBIOMSTR, PBIOMSTR2, &  !CTESSEL
!-TENDENCIES OUTPUT
 & PSNE1   , PTSNE1 , PASNE1,&
 & PRSNE1  , PWSNE1 , PTSAE1 , PTIAE1,&
 & PWLE1   , PWSAE1,&
 & PTLICEE1,PTLMNWE1,PTLWMLE1,&        
 & PTLBOTE1,PTLSFE1,PHLICEE1,PHLMLE1,&  
 & PUOE1   , PVOE1  , PTOE1  ,PSOE1,& 
 & PLAIE1  , PBSTRE1, PBSTR2E1,&                             !CTESSEL 
!-DDH OUTPUTS
 & PDHTSS  , PDHTTS , PDHTIS,&
 & PDHSSS  , PDHIIS , PDHWLS,&
 & PDHBIOS , PDHVEGS)                                        !CTESSEL

USE PARKIND1, ONLY : JPIM, JPRB
USE ISO_C_BINDING

! (C) Copyright 1993- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!     ------------------------------------------------------------------
!**** *SURFTSTP* - UPDATES LAND VALUES OF TEMPERATURE, MOISTURE AND SNOW.

!     PURPOSE.
!     --------
!          This routine updates the sea ice values of temperature
!     and the land values of soil temperature, skin soil water,
!     moisture in soil layers (in M scaled to the first reservoir depth),
!     snow mass (in M water equivalent), snow temperature, snow density
!     and snow albedo. For temperature, this is done via a forward time
!     step damped with some implicit linear considerations: as if all
!     fluxes that explicitely depend on the variable had only a linear
!     variation around the t-1 value. For moisture: (a) the evaporation
!     of the interception layer is treated implicitely, via a
!     linearization, while the other sources are treated explicitely;
!     (b) The soil transfer is "fully semi-implicit" similar to vertical
!     diffusion. Lower boundary conditons are no heat flux and free
!     drainage.

!**   INTERFACE.
!     ----------
!          *SURFTSTP* IS CALLED FROM *CALLPAR*.
!          THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE:
!     TSA,WL,WSA,SN AT T-1,TSK,SURFACE FLUXES COMPUTED IN OTHER PARTS OF
!     THE PHYSICS, AND W AND LAND-SEA MASK. IT RETURNS AS AN OUTPUT
!     TENDENCIES TO THE SAME VARIABLES (TSA,WL,WSA,SN).

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLEV*       NUMBER OF LEVELS
!    *KLON*       NUMBER OF GRID POINTS PER PACKET
!    *KLEVS*      NUMBER OF SOIL LAYERS
!    *KTILES*     NUMBER OF TILES (I.E. SUBGRID AREAS WITH DIFFERENT 
!                 SURFACE BOUNDARY CONDITION)
!    *KTVL*       VEGETATION TYPE FOR LOW VEGETATION FRACTION
!    *KTVH*       VEGETATION TYPE FOR HIGH VEGETATION FRACTION
!    *KSOTY*      SOIL TYPE                                   (1-7)
!    *KLEVSN*     Number of snow layers  
!    *KLEVI*      Number of sea ice layers (diagnostics)
!    *KLEVO*      NUMBER OF LAYERS OF OCEAN MIXED LAYER MODEL       
!    *KSTART*     FIRST TIMESTEP                                    
!    *KSTEP*      CURRENT TIMESTEP                                  

!     INPUT PARAMETERS (REAL):
!    *PSICE*      .LT.0.5 FOR LATITUDE ROWS WITH NO SEA ICE POINTS
!    *PEPS*       COEFFICIENT FOR TIME FILTER
!    *PTSPHY*     TIME STEP FOR THE PHYSICS                      S
!    *PSDOR*      OROGRAPHIC PARAMETER                         (m)
!    *PFRTI*      TILE FRACTIONS                              (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!    *PAHFSTI*    SURFACE SENSIBLE HEAT FLUX, FOR EACH TILE    W/M2
!    *PEVAPTI*    SURFACE MOISTURE FLUX, FOR EACH TILE         KG/M2/S
!    *PSSRFLTI*   NET SHORTWAVE RADIATION FLUX AT SURFACE, FOR
!                 EACH TILE                                    W/M2
!    *PSLRFLTI*   NET LONGWAVE  RADIATION AT THE SURFACE, FOR
!                 EACH TILE                                W/M**2
!    *PUO0*       HORIZONTAL VELOCITY OF OCEAN MIXED LAYER MODEL  
!    *PVO0*       HORIZONTAL VELOCITY OF OCEAN MIXED LAYER MODEL  
!    *PTO0*       TEMPERATURE OF OCEAN MIXED LAYER MODEL          
!    *PSO0*       SALINITY OF OCEAN MIXED LAYER MODEL             
!    *PADVT*      TEMPERATURE ADVECTION TERM FOR OCEAN MIXED LAYER MODEL 
!    *PADVS*      SALINITY ADVECTION TERM FOR OCEAN MIXED LAYER MODEL    
!    *PUSTOKES*   X-COMPONENT STOKES DRIFT
!    *PVSTOKES*   Y-COMPONENT STOKES DRIFT
!    *PTAUOCX*    X-COMPONENT SURFACE MOMENTUM FLUX TO OCEANS
!    *PTAUOCY*    Y-COMPONENT SURFACE MOMENTUM FLUX TO OCEANS
!    *PPHIOC*     ENERGY FLUX INTO OCEAN (DIMENSIONAL)
!    *PWSEMEAN*   WINDSEA VARIANCE
!    *PWSFMEAN*   WINDSEA MEAN FREQUENCY

!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE)
!    *LDSICE*     SEA ICE MASK (.T. OVER SEA ICE)
!    *LDLAKE*     LAKE MASK (.T. OVER LAKE)
!    *LDSI*       TRUE IF THERMALLY RESPONSIVE SEA-ICE
!    *LDNH*       TRUE FOR NORTHERN HEMISPHERE LATITUDE ROW
!    *LDOCN_KPP*  TRUE IF OCEAN MIXED LAYER GRID                    
!    *PCLAKE*     REAL LAKE COVER 

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!    *PSNM1M*     SNOW MASS (per unit area)                      kg/m**2
!    *PWSNM1M*    SNOW LIQUID WATER CONTENT (per unit area)      kg/m**2
!    *PTSAM1M*    SOIL TEMPERATURE                               K
!    *PTIAM1M*    SEA-ICE TEMPERATURE                            K
!          (NB: REPRESENTS THE FRACTION OF ICE ONLY, NOT THE GRID-BOX)
!    *PWLM1M*     SKIN RESERVOIR WATER CONTENT                   kg/m**2
!    *PWSAM1M*    SOIL MOISTURE                                m**3/m**3
!    *PTLICEM1M*  LAKE ICE TEMPERATURE                           K
!    *PTLMNWM1M*  LAKE MEAN WATER TEMPERATURE                    K
!    *PTLWMLM1M*  LAKE MIX-LAYER TEMPERATURE                     K
!    *PTLBOTM1M*  LAKE BOTTOM TEMPERATURE                        K
!    *PTLSFM1M*   LAKE SHAPE FACTOR (THERMOCLINE)                -
!    *PHLICEM1M*  LAKE ICE THICKNESS                             m
!    *PHLMLM1M*   LAKE MIX-LAYER THICKNESS                       m
!    *PGEMU*      SIN OF LATITUDE                                -
!    *PLDEPTH*    LAKE DEPTH                                     m
!    *PRSFC*      CONVECTIVE RAIN FLUX AT THE SURFACE          KG/M**2/S
!    *PRSFL*      LARGE SCALE RAIN FLUX AT THE SURFACE         KG/M**2/S
!    *PSLRFL*     NET LONGWAVE  RADIATION AT THE SURFACE         W/M**2
!    *PSSFC*      CONVECTIVE  SNOW FLUX AT THE SURFACE         KG/M**2/S
!    *PSSFL*      LARGE SCALE SNOW FLUX AT THE SURFACE         KG/M**2/S
!    *PCVL*       LOW VEGETATION COVER  (CORRECTED)              (0-1)
!    *PCVH*       HIGH VEGETATION COVER (CORRECTED)              (0-1)
!    *PCUR*       URBAN COVER (PASSIVE)                          (0-1)
!    *PWLMX*      MAXIMUM SKIN RESERVOIR CAPACITY                kg/m**2
!    *PEVAPSNW*   EVAPORATION FROM SNOW UNDER FOREST           KG/M**2/S
!    *POCDEPTH*     OCEAN DEPTH FOR OCEAN MIXED LAYER MODEL    (m)  
!    *PZO*          VERTICAL LAYER FOR OCEAN MIXED LAYER MODEL (m)  
!    *PDO*          VERTICAL LAYER FOR OCEAN MIXED LAYER MODEL (m)  
!    *PUSRF*      LOWEST LEVEL WIND U COMPONENT                   M/S
!    *PVSRF*      LOWEST LEVEL WIND V COMPONENT                   M/S
!    *PTSRF*      LOWEST LEVEL AIR TEMPERATURE                    K
!    *PAPRS*      LOWEST LEVEL AIR PRESSURE                       Pa
!    *POCDEPTH*     OCEAN DEPTH FOR OCEAN MIXED LAYER MODEL    (m)  
!    *PZO*          VERTICAL LAYER FOR OCEAN MIXED LAYER MODEL (m)  
!    *PDO*          VERTICAL LAYER FOR OCEAN MIXED LAYER MODEL (m)  
!     OUTPUT PARAMETERS (TENDENCIES):
!    *PSNE1*      SNOW MASS (per unit area) TENDENCY           kg/m**2/S
!    *PWSNE1*     SNOW LIQUD WATER CONTENT (per unit area) TENDENCY      kg/m**2/S
!    *PTSNE1*     SNOW TEMPERATURE                              K/S
!    *PASNE1*     SNOW ALBEDO                                   -/S
!    *PRSNE1*     SNOW DENSITY                                KG/M**3/S
!    *PTSAE1*     SOIL TEMPERATURE TENDENCY                     K/S
!    *PTIAE1*     SEA-ICE TEMPERATURE TENDENCY                  K/S
!          (NB: REPRESENTS THE FRACTION OF ICE ONLY, NOT THE GRID-BOX)
!    *PWLE1*      SKIN RESERVOIR WATER CONTENT TENDENCY       kg/m**2/s
!    *PWSAE1*     SOIL MOISTURE TENDENCY                   (m**3/m**3)/S
!    *PTLICEE1*   LAKE ICE TEMPERATURE TENDENCY                 K/S
!    *PTLMNWE1*   LAKE MEAN WATER TEMPERATURE TENDENCY          K/S
!    *PTLWMLE1*   LAKE MIX-LAYER TEMPERATURE TENDENCY           K/S
!    *PTLBOTE1*   LAKE BOTTOM TEMPERATURE TENDENCY              K/S
!    *PTLSFE1*    LAKE SHAPE FACTOR (THERMOCLINE) TENDENCY      -/S
!    *PHLICEE1*   LAKE ICE THICKNESS TENDENCY                   m/S
!    *PHLMLE1*    LAKE MIX-LAYER THICKNESS TENDENCY             m/S
!    *PWSAE1M*    SOIL MOISTURE TENDENCY                   (m**3/m**3)/S
!    *PUOE1*      VELOCITY TENDENCY OF OCEAN MIXED LAYER MODEL M/S**2
!    *PVOE1*      VELOCITY TENDENCY OF OCEAN MIXED LAYER MODEL M/S**2
!    *PTOE1*      TEMPERATURE TENDENCY OF OCEAN MIXED LAYER MODEL  K/S
!    *PSOE1*      SALINITY TENDENCY OF OCEAN MIXED LAYER MODEL  psu/S

!     OUTPUT PARAMETERS (DIAGNOSTIC):
!    *PTSDFL*     UPWARD FLUX BETWEEN SURFACE AND DEEP LAYER   W/M**2
!    *PROFD*      DEEP LAYER RUN-OFF                          kg/m**2/s
!    *PROFS*      SURFACE RUN-OFF                             kg/m**2/s
!    *PWFSD*      WATER FLUX BETWEEN LAYER 1 AND 2            kg/m**2/s
!    *PMELT*      WATER FLUX CORRESPONDING TO SNOW MELT       kg/m**2/s
!    *PFWEV*      EVAPORATION OVER LAND SURFACE               kg/m**2/s
!    *PENES*      SOIL ENERGY per unit area                   J/M**2
!    *PDHTSS*     Diagnostic array for snow T (see module yomcdh)
!    *PDHTTS*     Diagnostic array for soil T (see module yomcdh)
!    *PDHTIS*     Diagnostic array for ice T (see module yomcdh)
!    *PDHSSS*     Diagnostic array for snow mass (see module yomcdh)
!    *PDHIIS*     Diagnostic array for interception layer (see module yomcdh)
!    *PDHWLS*     Diagnostic array for soil water (see module yomcdh)
!    *PDIFM*      VISCOSITY OF OCEAN MIXED LAYER MODEL                
!    *PDIFT*      TEMPERATURE DIFFUSIVITY OF OCEAN MIXED LAYER MODEL  
!    *PDIFS*      SCALAR DIFFUSIVITY OF OCEAN MIXED LAYER MODEL       
!    *POTKE*      TURBULENT KINETIC ENERGY IN THE OCEAN MIXED LAYER   !OCEAN TKE


!     METHOD.
!     -------
!          STRAIGHTFORWARD ONCE THE DEFINITION OF THE CONSTANTS IS
!     UNDERSTOOD. FOR THIS REFER TO DOCUMENTATION. FOR THE TIME FILTER
!     SEE CORRESPONDING PART OF THE DOCUMENTATION OF THE ADIABATIC CODE.

!     EXTERNALS.
!     ----------
!          *SRFT*   COMPUTES THE TEMPERATURE CHANGES BEFORE THE SNOW
!                   MELTING.
!          *SRFWL*  COMPUTES THE SKIN RESERVOIR CHANGES.
!          *SRFWEXC*COMPUTES THE RUN-OFF AND SETS THE COEFFICIENTS OF
!                   THE TRIDIAGONAL MOISTURE SYSTEM OF EQUATIONS.
!          *SRFWDIF*SOLUTION OF THE TRIDIAGONAL MOISTURE SYSTEM OF EQUATIO
!          *SRFWINC*INCRMENTS OF SOIL MOSITURE, COMPUTATION OF DIAGNOSTICS
!          *SRFSN*  COMPUTES SNOW DEPTH CHANGES BEFORE SNOW MELTING.
!          *SRFSML* COMPUTES TEMPERATURE, SKIN RESERVOIR, SOIL WATER
!                   AND SNOW DEPTH CHANGES DUE TO SNOW MELTING.
!          *SRFWNG* CORRECTIONS TO AVOID NEGATIVE SOIL MOISTURE.

!     REFERENCE.
!     ----------
!          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
!     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.
!     ------------------------------------------------------------------

!     P.VITERBO       E.C.M.W.F.     16/03/93.
!     MODIFIED BY
!     SURFACE TILES   P.VITERBO/ACMB          8/03/99.
!     Surface DDH for TILES P. Viterbo        17/05/2000.
!     J.F. Estrade *ECMWF* 03-10-01 move in surf vob
!     P. Viterbo   *ECMWF* 04-05-24 New argument PENES; change surface units
!     G. Balsamo   *ECMWF* 06-07-03 Add soil type and orographic var. (runoff)
!     V. Stepanenko        08-01-22 Stress components, sea surface temperature
!                                   are added as arguments
!     E. Dutra/G. Balsamo  08-05-01 Add lake tile
!     E. Dutra/G. Balsamo  08-06-08 Add revised snow
!     Y. Takaya    *ECMWF* 08-10-07 Implement an ocean mixed layer model
!     E. Dutra             09-11-16 snow 2009 cleaning 
!     E. Dutra             10/10/2014      net longwave tiled 

IMPLICIT NONE

! Declaration of arguments

TYPE(C_PTR)       ,INTENT(IN)    :: YDSURF
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART       
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP        
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVSN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVI 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTTS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTTS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVSSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFSSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVIIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFIIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVWLS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFWLS 
!CTESSEL specific
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVBIOS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFBIOS
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVVEGS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFVEGS
INTEGER(KIND=JPIM),INTENT(INOUT) :: KVEG(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCVT(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAIVT(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAIL(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAIH(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PANFMVT(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PANDAYVT(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRESPBSTR(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRESPBSTR2(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBIOMASS_LAST(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBIOMASSTR_LAST(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBIOMASSTR2_LAST(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBLOSSVT(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBGAINVT(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAIE1(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBSTRE1(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBSTR2E1(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLAI(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBIOM(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBLOSS(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBGAIN(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBIOMSTR(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBIOMSTR2(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHBIOS(:,:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHVEGS(:,:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAILC(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAIHC(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(KLON)
!End CTESSEL specific
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVO         
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTY(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSDOR(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSTRTI(:,:)  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSTRTI(:,:)  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAHFSTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVAPTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFLTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFLTI(:,:)
LOGICAL           ,INTENT(IN)    :: LDLAND(:) 
LOGICAL           ,INTENT(IN)    :: LDSICE(:) 
LOGICAL           ,INTENT(IN)    :: LDSI 
LOGICAL           ,INTENT(IN)    :: LDNH(:) 
LOGICAL           ,INTENT(IN)    :: LDLAKE(:)  
LOGICAL           ,INTENT(IN)    :: LDOCN_KPP(:)  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLAKE(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSNM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PASNM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSNM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSNM1M(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSAM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTIAM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWLM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSAM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSFC(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSFL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSFC(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSFL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCUR(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWLMX(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVAPSNW(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSRF(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSRF(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSRF(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTSDFL(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PROFD(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PROFS(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWFSD(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMELT(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFWEV(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PENES(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSNE1(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSNE1(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PASNE1(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRSNE1(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWSNE1(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSAE1(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTIAE1(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWLE1(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWSAE1(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTSS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTTS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTIS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHSSS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHIIS(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDHWLS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLICEM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLMNWM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLWMLM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLBOTM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLSFM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICEM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLMLM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLDEPTH(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTLICEE1(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTLMNWE1(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTLWMLE1(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTLBOTE1(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTLSFE1(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHLICEE1(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHLMLE1(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZO(:,:)     
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PHO(:,:)     
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PHO_INV(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDO(:,:)     
REAL(KIND=JPRB)   ,INTENT(IN)    :: POCDEPTH(:)  
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PADVT(:,:)   
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PADVS(:,:)   
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUO0(:,:)    
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVO0(:,:)    
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUOC(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOC(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSTRC(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSTRC(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTO0(:,:)    
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSO0(:,:)    
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUOE1(:,:)   
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOE1(:,:)   
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTOE1(:,:)   
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSOE1(:,:)   
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTRI0(:,:)   
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTRI1(:,:)   
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSWDK_SAVE(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFM(:,:)   
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFT(:,:)   
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFS(:,:)   
REAL(KIND=JPRB)   ,INTENT(INOUT) :: POTKE(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSTOKES(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSTOKES(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUOCX(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUOCY(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPHIOC(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSEMEAN(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSFMEAN(:)

!     ------------------------------------------------------------------

END SUBROUTINE SURFTSTP
!
! 
! this space is left for visual comparison to SURFTSTP* and CALLPAR call
! The variables KLEVSN and KDH* are used in the external routine only
! The number of arguments per line is rearranged for readibility
!
!
!
!
!
!
!
!
!
!
!
!
!
END INTERFACE
