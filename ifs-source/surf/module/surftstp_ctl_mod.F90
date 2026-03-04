! (C) Copyright 1993- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
MODULE SURFTSTP_CTL_MOD
CONTAINS
SUBROUTINE SURFTSTP_CTL(KIDIA , KFDIA , KLON  , KLEVS , KTILES,&
 & KLEVSN , KLEVO , KLEVI , KSTART , KSTEP,&
 !
 !
 !
 ! this space is left for visual comparison to SURFTSTP.* and CALLPAR call
 ! The variables KLEVSN and KDH* are used in the external routine only
 ! The number of arguments per line is rearranged for readibility
 !
 !
 !
 & KTVL  , KTVH , KVEG , KSOTY,&
 & PTSPHY , PSDOR , PFRTI,&
 & PSST, PUSTRTI , PVSTRTI,&
 & PAHFSTI, PEVAPTI, PSSRFLTI,PSLRFLTI,&
 & PLAT, PANFMVT, PANDAYVT,PMU0,&                                  !CTESSEL
 & PCVT, PLAIVT, PLAIL, PLAIH,&                                    !CTESSEL
 & PLAILC,  PLAIHC, &
 & LDLAND,  LDSICE, LDSI, LDNH, LDOCN_KPP,&
 & PSNM1M  ,PTSNM1M, PASNM1M, PRSNM1M,PWSNM1M,&
 & PAPRS,   PTSKM1M ,PTSAM1M, PTIAM1M,&
 & PWLM1M  ,PWSAM1M,&
 & PTLICEM1M,PTLMNWM1M,PTLWMLM1M,PTLBOTM1M,PTLSFM1M,&
 & PHLICEM1M,PHLMLM1M,PGEMU,PLDEPTH,LDLAKE,PCLAKE, &
 & PRSFC   ,PRSFL,&
 & PSLRFL  ,PSSFC   ,PSSFL,&
 & PCVL    ,PCVH    ,PCUR    ,PWLMX   ,PEVAPSNW,&
 & PUSRF   ,PVSRF   ,PTSRF, &
 & PZO     ,PHO     ,PHO_INV ,PDO     ,POCDEPTH ,&
 & PUO0    ,PVO0    ,PUOC    ,PVOC    ,PTO0     ,&
 & PSO0    ,PADVT   ,PADVS   ,PTRI0   ,PTRI1    ,&
 & PSWDK_SAVE, PUSTRC ,PVSTRC, &
 & PUSTOKES,PVSTOKES,PTAUOCX , PTAUOCY, PPHIOC  ,&
 & PWSEMEAN,PWSFMEAN  ,&
 & YDCST   ,YDEXC   ,YDSOIL  ,YDVEG   ,YDFLAKE  ,&
 & YDURB   ,YDAGS   ,YDMLM   ,YDOCEAN_ML,&
!-DIAGNOSTICS OUTPUT
 & PTSDFL  , PROFD , PROFS,&
 & PWFSD   , PMELT , PFWEV, PENES,&
 & PDIFM   , PDIFT , PDIFS, POTKE,&
 & PRESPBSTR,PRESPBSTR2,PBIOMASS_LAST,&                            !CTESSEL
 & PBIOMASSTR_LAST,PBIOMASSTR2_LAST,&                              !CTESSEL
 & PBLOSSVT, PBGAINVT,&                                            !CTESSEL
 & PLAI    , PBIOM , PBLOSS, PBGAIN, PBIOMSTR, PBIOMSTR2,&         !CTESSEL
!-TENDENCIES OUTPUT
 & PSNE1   , PTSNE1 , PASNE1,&
 & PRSNE1  , PWSNE1 ,PTSAE1 , PTIAE1,&
 & PWLE1   , PWSAE1,&
 & PTLICEE1,PTLMNWE1,PTLWMLE1,&
 & PTLBOTE1,PTLSFE1,PHLICEE1,PHLMLE1,&
 & PUOE1   ,PVOE1   ,PTOE1  ,PSOE1,&
 & PLAIE1  ,PBSTRE1 ,PBSTR2E1, &                                   !CTESSEL
!-DDH OUTPUTS
 & PDHTSS  , PDHTTS , PDHTIS,&
 & PDHSSS  , PDHIIS , PDHWLS,&
 & PDHBIOS , PDHVEGS)                                              !CTESSEL

!-END OF CALL SURFTSTP

USE PARKIND1  , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_THF   , ONLY : RHOH2O
USE YOS_CST   , ONLY : TCST
USE YOS_EXC   , ONLY : TEXC
USE YOS_SOIL  , ONLY : TSOIL
USE YOS_VEG   , ONLY : TVEG
USE YOS_FLAKE , ONLY : TFLAKE
USE YOS_URB   , ONLY : TURB
USE YOS_AGS   , ONLY : TAGS
USE YOS_MLM   , ONLY : TMLM
USE YOS_OCEAN_ML , ONLY : TOCEAN_ML
USE SRFENE_MOD
USE SRFSN_MOD
USE SRFSN_LWIMP_MOD
USE SRFSN_DRIVER_MOD
USE SRFRCG_MOD
USE SRFT_MOD
USE SRFI_MOD
USE SRFWL_MOD
USE SRFWEXC_MOD
USE SRFWEXC_VG_MOD
USE SRFWDIF_MOD
USE SRFWINC_MOD
USE SRFWNG_MOD
USE SRFVEGEVOL_MOD
USE FLAKE_DRIVER_MOD



!!OIFS_rm USE OCEAN_ML_DRIVER_MOD   !KPP

USE OCEAN_ML_DRIVER_V2_MOD   !MIXED_LAYER_PJ

USE ABORT_SURF_MOD

!**** *SURFTSTP* - UPDATES LAND VALUES OF TEMPERATURE, MOISTURE AND SNOW.

!     PURPOSE.
!     --------
!          This routine updates the sea ice values of temperature
!     and the land values of soil temperature, interception layer
!     water (kg/m**2), multilayer soil water (in m**3/m**3),
!     snow mass (in kg/m**2), snow temperature, snow density
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
!    *KLEVI*      Number of sea ice layers (diagnostics)
!    *KLEVO*      NUMBER OF LAYERS OF OCEAN MIXED LAYER MODEL       !KPP
!    *KSTART*     FIRST TIMESTEP                                    !KPP
!    *KSTEP*      CURRENT TIMESTEP                                  !KPP

!     INPUT PARAMETERS (REAL):
!    *PSICE*      .LT.0.5 FOR LATITUDE ROWS WITH NO SEA ICE POINTS
!    *PEPS*       COEFFICIENT FOR TIME FILTER
!    *PTSPHY*     TIME STEP FOR THE PHYSICS                      S
!    *PSDOR*      OROGRAPHIC PARAMETER                        (m)
!    *PFRTI*      TILE FRACTIONS                              (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!            9 : LAKE                  10 : URBAN
!    *PAHFSTI*      SURFACE SENSIBLE HEAT FLUX, FOR EACH TILE  W/M2
!    *PSST*         SEA SURFACE TEMPERATURE                    K
!    *PUSTRTI*      X-STRESS                                   N/m2
!    *PVSTRTI*      Y-STRESS                                   N/m2
!    *PEVAPTI*      SURFACE MOISTURE FLUX, FOR EACH TILE      KG/M2/S
!    *PSSRFLTI*     NET SHORTWAVE RADIATION FLUX AT SURFACE, FOR
!                       EACH TILE                                 W/M2
!    *PSLRFLTI*       NET LONGWAVE  RADIATION AT THE SURFACE, FOR
!                       EACH TILE                                W/M**2
!    *PUO0*         HORIZONTAL VELOCITY OF OCEAN MIXED LAYER MODEL  !KPP
!    *PVO0*         HORIZONTAL VELOCITY OF OCEAN MIXED LAYER MODEL  !KPP
!    *PTO0*         TEMPERATURE OF OCEAN MIXED LAYER MODEL          !KPP
!    *PSO0*         SALINITY OF OCEAN MIXED LAYER MODEL             !KPP
!    *PADVT*        TEMPERATURE ADVECTION TERM FOR OCEAN MIXED LAYER MODEL !KPP
!    *PADVS*        SALINITY ADVECTION TERM FOR OCEAN MIXED LAYER MODEL    !KPP
!    *PUSTOKES*     X-COMPONENT STOKES DRIFT
!    *PVSTOKES*     Y-COMPONENT STOKES DRIFT
!    *PTAUOCX*      SURFACE MOMENTUM FLUX TO OCEANS X-DIRECTION
!    *PTAUOCY*      SURFACE MOMENTUM FLUX TO OCEANS Y-DIRECTION
!    *PPHIOC*       ENERGY FLUX INTO OCEANS (DIMENSIONAL)
!    *PWSEMEAN*     WINDSEA VARIANCE
!    *PWSFMEAN*     WINDSEA MEAN FREQUENCY

!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE)
!    *LDSICE*     SEA ICE MASK (.T. OVER SEA ICE)
!    *LDLAKE*     LAKE MASK (.T. OVER LAKE)
!    *LDSI*       TRUE IF THERMALLY RESPONSIVE SEA-ICE
!    *LDNH*       TRUE FOR NORTHERN HEMISPHERE LATITUDE ROW
!    *LDOCN_KPP   TRUE IF OCEAN MIXED LAYER GRID                    !KPP
!    *PCLAKE      REAL - lake cover

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!    *PSNM1M*     SNOW MASS (per unit area)                      KG/M**2
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
!    *PCUR*       URBAN COVER                                    (0-1)
!    *PWLMX*      MAXIMUM SKIN RESERVOIR CAPACITY                kg/m**2
!    *PEVAPSNW*   EVAPORATION FROM SNOW UNDER FOREST           KG/M**2/S
!    *POCDEPTH*   OCEAN DEPTH FOR OCEAN MIXED LAYER MODEL    (m)  !KPP
!    *PZO*        VERTICAL LAYER FOR OCEAN MIXED LAYER MODEL (m)  !KPP
!    *PDO*        VERTICAL LAYER FOR OCEAN MIXED LAYER MODEL (m)  !KPP
!    *PUSRF*      LOWEST LEVEL WIND U COMPONENT                   M/S
!    *PVSRF*      LOWEST LEVEL WIND V COMPONENT                   M/S
!    *PTSRF*      LOWEST LEVEL AIR TEMPERATURE                    K
!    *PAPRS*      LOWEST LEVEL AIR PRESSURE                       Pa
!    *POCDEPTH*   OCEAN DEPTH FOR OCEAN MIXED LAYER MODEL    (m)  !KPP
!    *PZO*        VERTICAL LAYER FOR OCEAN MIXED LAYER MODEL (m)  !KPP
!    *PDO*        VERTICAL LAYER FOR OCEAN MIXED LAYER MODEL (m)  !KPP
!     OUTPUT PARAMETERS (TENDENCIES):
!    *PSNE1*      SNOW MASS (per unit area) TENDENCY           KG/M**2/S
!    *PWSNE1*     SNOW LIQUD WATER CONTENT (per unit area) TENDENCY      kg/m**2/S
!    *PTSNE1*     SNOW TEMPERATURE                              K/S
!    *PASNE1*     SNOW ALBEDO                                   -/S
!    *PRSNE1*     SNOW DENSITY                                KG/M**3/S
!    *PTSAE1*     SOIL TEMPERATURE TENDENCY                     K/S
!    *PTIAE1*     SEA-ICE TEMPERATURE TENDENCY                  K/S
!          (NB: REPRESENTS THE FRACTION OF ICE ONLY, NOT THE GRID-BOX)
!    *PWLE1*      SKIN RESERVOIR WATER CONTENT TENDENCY       kg/m**2/S
!    *PWSAE1*     SOIL MOISTURE TENDENCY                   (m**3/m**3)/S
!    *PTLICEE1*   LAKE ICE TEMPERATURE TENDENCY                 K/S
!    *PTLMNWE1*   LAKE MEAN WATER TEMPERATURE TENDENCY          K/S
!    *PTLWMLE1*   LAKE MIX-LAYER TEMPERATURE TENDENCY           K/S
!    *PTLBOTE1*   LAKE BOTTOM TEMPERATURE TENDENCY              K/S
!    *PTLSFE1*    LAKE SHAPE FACTOR (THERMOCLINE) TENDENCY      -/S
!    *PHLICEE1*   LAKE ICE THICKNESS TENDENCY                   m/S
!    *PHLMLE1*    LAKE MIX-LAYER THICKNESS TENDENCY             m/S
!    *PTOE1*      TEMPERATURE TENDENCY OF OCEAN MIXED LAYER MODEL  K/S
!    *PTOE1*      SALINITY TENDENCY OF OCEAN MIXED LAYER MODEL  psu/S
!    *PUOE1*      VELOCITY TENDENCY OF OCEAN MIXED LAYER MODEL M/S**2
!    *PVOE1*      VELOCITY TENDENCY OF OCEAN MIXED LAYER MODEL M/S**2

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
!    *PDIFM*      VISCOSITY OF OCEAN MIXED LAYER MODEL                !KPP
!    *PDIFT*      TEMPERATURE DIFFUSIVITY OF OCEAN MIXED LAYER MODEL  !KPP
!    *PDIFS*      SCALAR DIFFUSIVITY OF OCEAN MIXED LAYER MODEL       !KPP
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
!          *SRFWDIF*SOLUTION OF THE TRIDIAGONAL MOISTURE SYSTEM OF EQUATIONS
!          *SRFWINC*INCREMENTS OF SOIL MOSITURE, COMPUTATION OF DIAGNOSTICS
!          *SRFSN*  COMPUTES SNOW DEPTH CHANGES BEFORE SNOW MELTING.
!          *SRFSN_LWIMP*  REVISED SNOW SCHEME W. DIAG. LIQ. WATER.
!          *SRFSML* COMPUTES TEMPERATURE, SKIN RESERVOIR, SOIL WATER
!                   AND SNOW DEPTH CHANGES DUE TO SNOW MELTING.
!          *SRFWNG* CORRECTIONS TO AVOID NEGATIVE SOIL MOISTURE.

!     REFERENCE.
!     ----------
!          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
!     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.

!     P.VITERBO       E.C.M.W.F.     16/03/93.
!     MODIFIED BY
!     SURFACE TILES   P.VITERBO/ACMB          8/03/99.
!     Surface DDH for TILES P. Viterbo        17/05/2000.
!     J.F. Estrade *ECMWF* 03-10-01 move in surf vob
!     P. Viterbo    24-05-2004      New argument PENES; Change surface units
!     G. Balsamo    03-07-2004      Add soil type and orographic var (runoff)
!     V. Stepanenko/G. Balsamo 01-05-2008   Add lake parameterization (FLAKE)
!     Y. Takaya     07-10-2008      Add ocean mixed layer model (KPP scheme)
!     E. Dutra      12-11.2008      Include new snow parameterization
!     E. Dutra      16-11-2009      snow 2009 cleaning
!     P. JANSSEN    02-08-2011      ADD OCEAN MIXED LAYER MODEL THAT INCLUDES
!                                   SEA STATE EFFECTS
!     E. Dutra      10/10/2014      net longwave tiled
!     G. Balsamo    05-02-2015      soil and lake tendencies updates for fractional ice
!     ------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVSN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVI

!CTESSEL specific
INTEGER(KIND=JPIM),INTENT(INOUT) :: KVEG(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCVT(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:)
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(:)
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
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTLICEE1(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTLMNWE1(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTLWMLE1(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTLBOTE1(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTLSFE1(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PHLICEE1(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PHLMLE1(:)
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
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TEXC)        ,INTENT(IN)    :: YDEXC
TYPE(TSOIL)       ,INTENT(IN)    :: YDSOIL
TYPE(TVEG)        ,INTENT(IN)    :: YDVEG
TYPE(TFLAKE)      ,INTENT(IN)    :: YDFLAKE
TYPE(TURB)        ,INTENT(IN)    :: YDURB
TYPE(TAGS)        ,INTENT(IN)    :: YDAGS
TYPE(TMLM)        ,INTENT(IN)    :: YDMLM
TYPE(TOCEAN_ML)   ,INTENT(INOUT) :: YDOCEAN_ML

!*         0.2    DECLARATION OF LOCAL VARIABLES.
!                 ----------- -- ----- ----------

REAL(KIND=JPRB) :: ZSN(KLON,KLEVSN),ZTSN(KLON,KLEVSN),ZASN(KLON),&
 & ZRSN(KLON,KLEVSN),ZTSA(KLON,KLEVS)      ,ZTIA(KLON,KLEVS)   ,&
 & ZWL(KLON)      ,ZWSA(KLON,KLEVS),ZWSN(KLON,KLEVSN)
REAL(KIND=JPRB) :: ZCFW(KLON,KLEVS),ZRHSW(KLON,KLEVS), ZWSADIF(KLON,KLEVS),&
 & ZSAWGFL(KLON,KLEVS),ZCTSA(KLON,KLEVS)
REAL(KIND=JPRB) :: ZTSFC(KLON),ZTSFL(KLON)
REAL(KIND=JPRB) :: ZRSFC(KLON)    ,&
                 & ZRSFL(KLON)    ,&
                 & ZSSFC(KLON)    ,ZSSFL(KLON)
REAL(KIND=JPRB) :: ZLW(KLON), ZPMSNINT(KLON)
! Test ad
REAL(KIND=JPRB) :: ZTSFC_TAD(KLON),ZTSFL_TAD(KLON), ZSSFC_TAD(KLON), ZSSFL_TAD(KLON)

REAL(KIND=JPRB) :: ZFWEL1(KLON)   ,ZFWE234(KLON)
REAL(KIND=JPRB) :: ZGSN(KLON) ,ZMSN(KLON),ZEMSSN(KLON),ZEINTTI(KLON,KTILES)
REAL(KIND=JPRB) :: ZCDAWZ(KLON,KLEVS)
REAL(KIND=JPRB) :: ZSLRFLTI(KLON,KTILES)
REAL(KIND=JPRD) :: ZTSPHY
REAL(KIND=JPRB) :: ZHOH2O,ZSNPERT
REAL(KIND=JPRB) :: ZRSNM1M(KLON,KLEVSN)

INTEGER(KIND=JPIM) :: JK, JL, JT,KLMAX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

REAL(KIND=JPRD) :: ZTLICE(KLON),ZTLMNW(KLON),ZTLWML(KLON),&
                   &ZTLBOT(KLON),ZTLSF(KLON),ZHLICE(KLON),ZHLML(KLON)
REAL(KIND=JPRB) :: ZTSFCIN(KLON),ZTSFLIN(KLON),ZROFS(KLON)

!     ------------------------------------------------------------------

!*         1.1   SET UP SOME CONSTANTS, INITIALISE ARRAYS.
!                --- -- ---- -----------------------------
IF (LHOOK) CALL DR_HOOK('SURFTSTP_CTL_MOD:SURFTSTP_CTL',0,ZHOOK_HANDLE)
ASSOCIATE(RTT=>YDCST%RTT, &
 & LELWTL=>YDEXC%LELWTL, &
 & LEFLAKE=>YDFLAKE%LEFLAKE, RH_ICE_MIN_FLK=>YDFLAKE%RH_ICE_MIN_FLK, &
 & LOCMLTKE=>YDMLM%LOCMLTKE, &
 & LEOCML=>YDOCEAN_ML%LEOCML, &
 & LESN09=>YDSOIL%LESN09,LESNML=>YDSOIL%LESNML, LEVGEN=>YDSOIL%LEVGEN, RALFMINPSN=>YDSOIL%RALFMINPSN, &
 & RDAW=>YDSOIL%RDAW, RHOMAXSN=>YDSOIL%RHOMAXSN, &
 & LECTESSEL=>YDVEG%LECTESSEL)

ZTSPHY=(1.0_JPRD)/PTSPHY
ZHOH2O=1.0_JPRB/RHOH2O

DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA
    ZCDAWZ(JL,JK)=RDAW(JK)
  ENDDO
ENDDO
DO JL=KIDIA,KFDIA

  ! decide if longwave is tiled or not, if not set the same value for all tiles
  IF (LELWTL) THEN
    DO JT=1,KTILES
      ZSLRFLTI(JL,JT) = PSLRFLTI(JL,JT)
    ENDDO
  ELSE
    DO JT=1,KTILES
      ZSLRFLTI(JL,JT) = PSLRFL(JL)
    ENDDO
  ENDIF


  ZFWEL1(JL)=0.0_JPRB
  ZFWE234(JL)=0.0_JPRB
  PROFS(JL)=0.0_JPRB

! protect against negative precip

  ZRSFC(JL)=MAX(0.0_JPRB,PRSFC(JL))
  ZRSFL(JL)=MAX(0.0_JPRB,PRSFL(JL))
  ZSSFC(JL)=MAX(0.0_JPRB,PSSFC(JL))
  ZSSFL(JL)=MAX(0.0_JPRB,PSSFL(JL))
ENDDO

CALL SRFENE(KIDIA, KFDIA  , KLON  , KLEVS,&
 & LDLAND , LDSICE ,&
 & PTSAM1M, KSOTY, PCVL , PCVH ,&
 & YDCST  , YDSOIL ,&
 & PENES)


!     ------------------------------------------------------------------
!*         0.1     INTERCEPTION RESERVOIR AND ASSOCIATED RUNOFF CHANGES.

CALL SRFWL(KIDIA,KFDIA,KTILES,&
 & PTSPHY,PWLM1M,PCVL,PCVH,PWLMX,&
 & PFRTI  , PEVAPTI ,ZRSFC,ZRSFL,PEVAPSNW,&
 & ZWL,ZFWEL1,ZTSFC,ZTSFL,ZEINTTI,&
 & LDLAND,&
 & YDSOIL,YDVEG,&
 & PDHIIS)

!     ------------------------------------------------------------------

!*         2.     SNOW CHANGES.
!                 ---- --------
! IF ( LESNML ) THEN
IF ( LESNML ) THEN
!
  CALL SRFSN_DRIVER(KIDIA   ,KFDIA   ,KLON   ,KLEVSN, PTSPHY, LDLAND,&
  & PSDOR, &
  ! input at T-1 prognostics
  & PSNM1M, PTSNM1M, PRSNM1M  ,PWSNM1M, PASNM1M, &
  ! input at T-1 fluxes or constants
  & PFRTI, PTSAM1M , PUSRF, PVSRF, PTSRF ,&
  & ZSSFC, ZSSFL   , ZTSFC, ZTSFL, &
  & PSLRFLTI, PSSRFLTI , PAHFSTI, PEVAPTI, PEVAPSNW, &
  & PWSAM1M , KSOTY, PAPRS, &
  ! input derived types (constants)
  & YDSOIL , YDCST, &
  ! output prognostics at T
  & ZSN,  ZTSN, ZRSN, ZWSN, ZASN, &
  ! output fluxes
  & ZGSN, ZMSN, ZEMSSN, &
  ! Diagnostics
  & PDHTSS, PDHSSS)

ELSE


  IF ( LESN09 ) THEN
  !* LIQUID WATER DIAGNOSTIC
    CALL SRFSN_LWIMP(KIDIA  ,KFDIA  ,KLON   ,KTILES   ,PTSPHY, LDLAND, &
  & PSNM1M(:,1) ,PTSNM1M(:,1) ,PASNM1M,PRSNM1M(:,1),PTSAM1M,PHLICEM1M,  &
  & ZSLRFLTI, PSSRFLTI,PFRTI  ,PAHFSTI,PEVAPTI,            &
  & ZSSFC  ,ZSSFL   ,PEVAPSNW,                           &
  & ZTSFC  ,ZTSFL   ,PUSRF  ,PVSRF   ,PTSRF,             &
  & YDCST  ,YDVEG   ,YDSOIL ,YDFLAKE, YDEXC,             &
  & ZSN(:,1)    ,ZTSN(:,1)    ,ZASN   ,ZRSN(:,1)   ,ZGSN   ,ZMSN,       &
  & ZEMSSN ,ZTSFCIN,ZTSFLIN,                             &
  & PDHTSS , PDHSSS )

    ! CHANGE THROUGHFALL AT THE SURFACE
    DO JL=KIDIA,KFDIA
      ZTSFC(JL)=ZTSFC(JL)-ZTSFCIN(JL)
      ZTSFL(JL)=ZTSFL(JL)-ZTSFLIN(JL)
    ENDDO

  ELSE
  !* BASE ORIGINAL SNOW SCHEME
    CALL SRFSN(      KIDIA  ,KFDIA  ,KLON   ,KTILES, PTSPHY,          &
  & PSNM1M(:,1) ,PTSNM1M(:,1) ,PASNM1M,PRSNM1M(:,1),PTSAM1M,PHLICEM1M,       &
  & PSLRFL ,PSSRFLTI,PFRTI  ,PAHFSTI,PEVAPTI,                 &
  & ZSSFC  ,ZSSFL   ,PEVAPSNW ,                               &
  & YDCST  ,YDVEG   ,YDSOIL ,YDFLAKE,                         &
  & ZSN(:,1)    ,ZTSN(:,1)    ,ZASN   ,ZRSN(:,1)   ,ZGSN   ,ZMSN   ,         &
  & ZEMSSN ,                                                  &
  & PDHTSS , PDHSSS )


  ENDIF
  !! If multi-layer is bypassed :
  ZWSN(KIDIA:KFDIA,:) = 0._JPRB
  IF (KLEVSN > 1 ) THEN
    DO JK=2,KLEVSN
      ZSN(KIDIA:KFDIA,JK) = 0._JPRB
      ZTSN(KIDIA:KFDIA,JK) = ZTSN(KIDIA:KFDIA,1)
      ZRSN(KIDIA:KFDIA,JK) = ZRSN(KIDIA:KFDIA,1)
    ENDDO
  ENDIF


ENDIF
ZROFS(KIDIA:KFDIA)=0._JPRB
! Take care of permanent snow areas!
!Permanent snow reset was potentially dangerous at 1m threshold and 9m is used instead.
ZSNPERT=9000.0_JPRB ! permanent snow threshold
KLMAX=KLEVSN ! permanent snow max index
IF ( .NOT. LESNML ) KLMAX=1
DO JL=KIDIA,KFDIA
  IF ( SUM(PSNM1M(JL,:)) >=ZSNPERT ) THEN
    ZASN(JL)=RALFMINPSN
    IF (.NOT. YDSOIL%LESNWBCON ) THEN
      !=======***NOTE***==========
       ! SETTING ZSN IN PERMANENT SNOW AREAS DOES NOT CONSERVE MASS !!!!
       ! IF MASS IS TO BE CONSERVED THE FOLLOWING LINE SHOULD BE COMMENTED !
       IF ( .NOT. LESNML ) THEN
         ZRSN(JL,KLMAX)=RHOMAXSN
         ZSN(JL,KLMAX) = 10000.0_JPRB !reset to glaciers value of 10000 kg/m2 (SWE=10m)
       ELSE ! RSN is not interactive over glacier, fixed to 300
         ZRSN(JL,1:KLEVSN) =RHOMAXSN
         ZWSN(JL,1:KLEVSN) = 0._JPRB
         IF (KLMAX /= KLEVSN) THEN
           ZSN(JL,KLMAX)  = 10000.0_JPRB - SUM(ZSN(JL,1:KLMAX-1)) - SUM(ZSN(JL,KLMAX+1:KLEVSN))  !reset to glaciers value of 10000 kg/m2 (SWE=10m)
         ELSE
           ZSN(JL,KLMAX)  = 10000.0_JPRB - SUM(ZSN(JL,1:KLMAX-1))                        !reset to glaciers value of 10000 kg/m2 (SWE=10m)
         ENDIF
       ENDIF
    ELSE
      ! try to keep water balance by removing snow mass > 10000 as calving - runoff
      ! this only applies to the single layer
      IF ( .NOT. LESNML ) THEN
        ZRSN(JL,KLMAX)=RHOMAXSN
        ZROFS(JL) = MAX(0._JPRB,ZSN(JL,KLMAX)-SUM(PSNM1M(JL,:)))*ZTSPHY
        ZSN(JL,KLMAX)=10000.0_JPRB
      ENDIF
    ENDIF
  ENDIF
ENDDO

!*         2b.    URBAN CHANGES.
!                 ---- --------
  ! INTERCEPTION BY RUN-OFF ASSOCIATED WITH URBAN FRACTION
  ! CHANGE THROUGHFALL AT THE SURFACE (SRFWEXC_VG ACCOUNTS FOR ROF FROM SUB-URB LAYER)
  !   ROUGH ESTIMATE RUNOFF OF 0.3 - Drainage estimate will need to be revised or updated based on mapping
  ! Likely to be an underestimate - Paul & Meyer 2001 suggest 55% runoff on fully urban surfaces
IF (KTILES .GT. 9) THEN
 DO JL=KIDIA,KFDIA
  ZTSFC(JL)=ZTSFC(JL)*(1.0_JPRB-(0.3_JPRB*PFRTI(JL,10)))
  ZTSFL(JL)=ZTSFL(JL)*(1.0_JPRB-(0.3_JPRB*PFRTI(JL,10)))
 ENDDO
ENDIF

!     ------------------------------------------------------------------
!*         3.     SOIL HEAT BUDGET.
!                 ---- ---- -------
!*         3.1     COMPUTE SOIL HEAT CAPACITY.

CALL SRFRCG(KIDIA  , KFDIA  , KLON   ,KTILES, KLEVS ,&
 & LDLAND , LDSICE ,&
 & PTSAM1M, KSOTY  , PCVL   , PCVH   ,PCUR  ,&
 & YDCST, YDSOIL   , YDURB  ,&
 & ZCTSA)

!*         3.2     TEMPERATURE CHANGES.
! soil
CALL SRFT(KIDIA  , KFDIA  , KLON   , KTILES, KLEVS  , &
 & PTSPHY , PTSAM1M, PWSAM1M, KSOTY, &
 & PFRTI  , PAHFSTI,PEVAPTI , &
 & ZSLRFLTI , PSSRFLTI, ZGSN  , &
 & ZCTSA  , LDLAND , &
 & YDCST  , YDSOIL , YDFLAKE, YDURB,&
 & ZTSA   ,PTSDFL  , PDHTTS)
! sea-ice
IF (LDSI) THEN
  CALL SRFI(&
   & KIDIA  , KFDIA  , KLON   , KLEVS  , KLEVI  , &
   & PTSPHY , PTIAM1M, PFRTI  , PAHFSTI, PEVAPTI, &
   & ZSLRFLTI(:,2) ,PSSRFLTI, LDSICE , LDNH   , &
   & YDCST  ,YDSOIL  ,&
   & ZTIA   ,PDHTIS)
ELSE
  DO JK=1,KLEVS
    DO JL=KIDIA,KFDIA
      ZTIA(JL,JK)=PTIAM1M(JL,JK)
      IF (SIZE(PDHTIS) > 0) THEN
        PDHTIS(JL,JK,1:4)=0.0_JPRB
        PDHTIS(JL,JK,9:12)=0.0_JPRB
      ENDIF
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------
!*         4.     SOIL WATER BUDGET.
!                 ---- ----- -------
!!*         4.1     INTERCEPTION RESERVOIR AND ASSOCIATED RUNOFF CHANGES.
!
!CALL SRFWL(KIDIA,KFDIA,KTILES,&
! & PTSPHY,PWLM1M,PCVL,PCVH,PWLMX,&
! & PFRTI  , PEVAPTI ,ZRSFC,ZRSFL,PEVAPSNW,&
! & ZWL,ZFWEL1,ZTSFC,ZTSFL,ZEINTTI,&
! & LDLAND,PDHIIS)

!*         4.2  COMPUTATION OF SOIL MOISTURE DIFFUSIVITY AND RIGHT HAND
!                 SIDES
IF (LEVGEN) THEN
  CALL SRFWEXC_VG(KIDIA,KFDIA,KLON,KLEVS,KTILES,&
 & PTSPHY,KTVL,KTVH,KSOTY,PSDOR,PFRTI,PEVAPTI,&
 & PWSAM1M,PTSAM1M,PCUR,&
 & ZTSFC,ZTSFL,ZMSN,ZEMSSN,ZEINTTI,PEVAPSNW,&
 & YDSOIL,YDVEG,YDURB,&
 & PROFS,ZCFW,ZRHSW,&
 & ZSAWGFL,ZFWEL1,ZFWE234,&
 & LDLAND,PDHWLS)
ELSE
  CALL SRFWEXC(KIDIA,KFDIA,KLON,KLEVS,KTILES,&
 & PTSPHY,KTVL,KTVH,PFRTI,PEVAPTI,&
 & PWSAM1M,PTSAM1M,PCUR,&
 & ZTSFC,ZTSFL,ZMSN,ZEMSSN,ZEINTTI,PEVAPSNW,&
 & YDSOIL,YDVEG,YDURB,&
 & PROFS,ZCFW,ZRHSW,&
 & ZSAWGFL,ZFWEL1,ZFWE234,&
 & LDLAND,PDHWLS)
ENDIF


!*         4.3  SOLUTION OF TRIDIAGONAL SYSTEM OF EQUATIONS

CALL SRFWDIF(KIDIA,KFDIA,KLON,KLEVS,PWSAM1M,ZCFW,ZRHSW,ZCDAWZ,YDSOIL,ZWSADIF,LDLAND)

!*         4.4  UPDATE SOIL MOISTURE VALUES

CALL SRFWINC(KIDIA,KFDIA,KLEVS,&
 & PTSPHY,PWSAM1M,ZWSADIF,ZSAWGFL,ZCFW,&
 & YDSOIL,&
 & ZWSA,PWFSD,LDLAND,PDHWLS)

!*         4.5     WATER CORRECTIONS.

CALL SRFWNG(KIDIA,KFDIA,KLEVS,PTSPHY,KSOTY,&
 & ZWL,PWLMX,ZWSA,&
 & YDSOIL,&
 & PROFS,PROFD,PWFSD,&
 & LDLAND,PDHWLS)

!*         4.6     CONTRIBUTION OF BOTTOM DRAINAGE TO TOTAL RUNOFF.

WHERE (LDLAND(KIDIA:KFDIA)) PROFD(KIDIA:KFDIA)=&
 & PROFD(KIDIA:KFDIA)+ZSAWGFL(KIDIA:KFDIA,KLEVS)*PTSPHY


!*        4.7     Lake model FLAKE

IF (LEFLAKE) THEN
  CALL FLAKE_DRIVER( KIDIA, KFDIA, KLON,LDLAKE,&
                 & PSSRFLTI(:,9),ZSLRFLTI(:,9),PAHFSTI(:,9),PEVAPTI(:,9), &
                 & PUSTRTI(:,9),PVSTRTI(:,9),PFRTI,&
                 & PLDEPTH,PGEMU,PTSPHY,&
                 & YDCST,YDFLAKE,&
                 & PTLICEM1M,PTLMNWM1M,PTLWMLM1M,PTLBOTM1M,PTLSFM1M, &
                 & PHLICEM1M,PHLMLM1M,&
                 & ZTLICE,ZTLMNW,ZTLWML,&
                 & ZTLBOT,ZTLSF,ZHLICE,ZHLML )

ELSE ! The prognostic variables are set equal to previous timestep

  ZTLICE  = PTLICEM1M
  ZTLMNW  = PTLMNWM1M
  ZTLWML  = PTLWMLM1M
  ZTLBOT  = PTLBOTM1M
  ZTLSF   = PTLSFM1M
  ZHLICE  = PHLICEM1M
  ZHLML   = PHLMLM1M

ENDIF

! Set meaningful value for FLAKE prognostic variables on ocean points (using proxy sea temperatures/ice conditions)
DO JL =KIDIA,KFDIA
! Over land LEFLAKE is set to true even if lakes have zero fraction to get meaningful
! lake prognostic values (this is useful in FULLPOS horizontal interpolation)
    IF (.NOT.LDLAND(JL)  .AND. .NOT. LDLAKE(JL) ) THEN
!sea points that have no lake fraction are set to SST
      ZTLWML(JL)=MAX(PSST(JL),RTT)
      ZTLMNW(JL)=MAX(PSST(JL),RTT)
      ZTLBOT(JL)=MAX(PSST(JL),RTT)
      ZTLSF(JL) = 0.65_JPRB ! The lake shape factor is kept constant to 0.65
      ZTLICE(JL)= ZTIA(JL,1)
      IF (LDSICE(JL)) THEN  ! The sea-ice presence is used to set a ice-depth over ocean (capped to 0.1m)
         ZHLICE(JL)= MIN(PFRTI(JL,2)/10.0_JPRB,0.1_JPRB)
      ELSE
         ZHLICE(JL)= 0.0_JPRB
      ENDIF
      ZHLML(JL) = MIN(PLDEPTH(JL),50.0_JPRB)
    ENDIF
    IF (.NOT.LDLAND(JL)  .AND. LDLAKE(JL) ) THEN
      ZTSA(JL,:)=ZTLWML(JL) ! The lake temperature is put into the soil temperature over resolved lakes.
    ENDIF
ENDDO

!*        4.8     Ocean Mixed layer
! ! OIFS_rm IF(LEOCML) THEN
! ! OIFS_rm   CALL OCEAN_ML_DRIVER &
! ! OIFS_rm  & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KSTART   ,&
! ! OIFS_rm  &   KSTEP    ,LDOCN_KPP,LDSICE   ,PSSRFLTI(:,1),ZSLRFLTI(:,1),&
! ! OIFS_rm  &   PAHFSTI(:,1), PEVAPTI(:,1),PUSTRTI(:,1), PVSTRTI(:,1),&
! ! OIFS_rm  &   PUSTRC   ,PVSTRC   ,ZRSFC    ,ZRSFL    ,ZSSFC    ,&
! ! OIFS_rm  &   ZSSFL    ,PGEMU    ,PZO      ,PHO      ,PHO_INV  ,&
! ! OIFS_rm  &   PDO      ,POCDEPTH ,PUO0     ,PVO0     ,PUOC     ,&
! ! OIFS_rm  &   PVOC     ,PTO0     ,PSO0     ,PUOE1    ,PVOE1    ,&
! ! OIFS_rm  &   PTOE1    ,PSOE1    ,PADVT    ,PADVS    ,PTSPHY   ,&
! ! OIFS_rm  &   PDIFM    ,PDIFT    ,PDIFS    ,PTRI0    ,PTRI1    ,&
! ! OIFS_rm  &   PSWDK_SAVE, YDCST  ,YDOCEAN_ML )
! ! OIFS_rm ENDIF


IF(LOCMLTKE) THEN
  CALL OCEAN_ML_DRIVER_V2 &
 & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,PTSPHY   ,&
 &   PGEMU    ,PSSRFLTI(:,1),ZSLRFLTI(:,1),PAHFSTI(:,1),&
 &   PEVAPTI(:,1),  PUSTOKES  ,PVSTOKES ,&
 &   PTAUOCX  ,PTAUOCY   ,PPHIOC  ,PWSEMEAN ,PWSFMEAN ,&
 &   PUO0     ,PVO0     ,PTO0     ,PSO0     ,&
 &   YDCST    ,YDEXC    ,YDMLM    ,&
 &   POTKE    ,PUOE1    ,PVOE1    ,PTOE1    ,PSOE1    )
ENDIF
!     ------------------------------------------------------------------
!*         5.     COMPUTE TENDENCIES (TEND. FOR TSK COMP. IN VDIFF)
!                 ------- ---------- ------------------------------
!*         5.1   TENDENCIES

DO JK=1,KLEVSN
  DO JL=KIDIA,KFDIA
    PSNE1(JL,JK)   =PSNE1(JL,JK)+&
   & (ZSN(JL,JK)-PSNM1M(JL,JK))*ZTSPHY
    PTSNE1(JL,JK)  =PTSNE1(JL,JK)+&
   & (ZTSN(JL,JK)-PTSNM1M(JL,JK))*ZTSPHY
    PRSNE1(JL,JK)   =PRSNE1(JL,JK)+&
   & (ZRSN(JL,JK)-PRSNM1M(JL,JK))*ZTSPHY
    PWSNE1(JL,JK)   =PWSNE1(JL,JK)+&
   & (ZWSN(JL,JK)-PWSNM1M(JL,JK))*ZTSPHY
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  PASNE1(JL)  =PASNE1(JL)+&
   & (ZASN(JL)-PASNM1M(JL))*ZTSPHY
  PTLICEE1(JL)   =PTLICEE1(JL)+&
   & (ZTLICE(JL)-PTLICEM1M(JL))*ZTSPHY
  PTLMNWE1(JL)   =PTLMNWE1(JL)+&
   & (ZTLMNW(JL)-PTLMNWM1M(JL))*ZTSPHY
  PTLWMLE1(JL)   =PTLWMLE1(JL)+&
   & (ZTLWML(JL)-PTLWMLM1M(JL))*ZTSPHY
  PTLBOTE1(JL)   =PTLBOTE1(JL)+&
   & (ZTLBOT(JL)-PTLBOTM1M(JL))*ZTSPHY
  PTLSFE1(JL)   =PTLSFE1(JL)+&
   & (ZTLSF(JL)-PTLSFM1M(JL))*ZTSPHY
  PHLICEE1(JL)   =PHLICEE1(JL)+&
   & (ZHLICE(JL)-PHLICEM1M(JL))*ZTSPHY
  PHLMLE1(JL)   =PHLMLE1(JL)+&
   & (ZHLML(JL)-PHLMLM1M(JL))*ZTSPHY
ENDDO

DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA
    PTIAE1(JL,JK)=PTIAE1(JL,JK)+&
     & (ZTIA(JL,JK)-PTIAM1M(JL,JK))*ZTSPHY
  ENDDO
ENDDO
! PTSA represents the grid-box average T; the sea-ice fraction
!    has a T representative of the fraction
!    this is calculated only if the fraction of lakes is 0
!    (otherwise lake temperature can contaminate soil temperature tendency)
DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA
    IF (PFRTI(JL,9).GT.0.0_JPRB) THEN
       PTSAE1(JL,JK)=PTSAE1(JL,JK)+&
       & ((ZTSA(JL,JK)-PTSAM1M(JL,JK)))*ZTSPHY
    ELSE
       PTSAE1(JL,JK)=PTSAE1(JL,JK)+&
       & (PFRTI(JL,2)*&
       & (ZTIA(JL,JK)-PTIAM1M(JL,JK))+&
       & (1.0_JPRB-PFRTI(JL,1)-PFRTI(JL,2))*&
       & (ZTSA(JL,JK)-PTSAM1M(JL,JK)))*ZTSPHY
    ENDIF
  ENDDO
ENDDO
DO JL=KIDIA,KFDIA
  PWLE1(JL)   =PWLE1(JL)+&
   & (ZWL(JL)-PWLM1M(JL))*ZTSPHY
ENDDO
DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA
    PWSAE1(JL,JK)=PWSAE1(JL,JK)+&
     & (ZWSA(JL,JK)-PWSAM1M(JL,JK))*ZTSPHY
  ENDDO
ENDDO

!*         5.2  DIAGNOSTICS

DO JL=KIDIA,KFDIA
  PFWEV(JL) =ZFWEL1(JL)+ZFWE234(JL)
  PROFS(JL) =PROFS(JL)*ZTSPHY
  PROFD(JL) =PROFD(JL)*ZTSPHY
  PTSDFL(JL)=-PTSDFL(JL)
  PMELT(JL) =ZMSN(JL)
ENDDO

!*         5.3 precipitaiton over resolved lakes as surface runoff
IF (YDSOIL%LEROLAKE) THEN
  DO JL=KIDIA,KFDIA
    IF ( .NOT. LDLAND(JL)  .AND. LDLAKE(JL) ) THEN
      PROFS(JL)=PROFS(JL)+ZRSFC(JL)+ZRSFL(JL)+ZSSFC(JL)+ZSSFL(JL)
    ENDIF
    !* Add ice calving to surface runoff
    PROFS(JL)=PROFS(JL)+ZROFS(JL)
  ENDDO
ENDIF

IF (LECTESSEL) THEN


!*         6.     INTERACTIVE VEGETATION.

CALL SRFVEGEVOL(KIDIA,KFDIA,KLON,KLEVS,KSTEP, &
 & PTSPHY, PLAT,PMU0,&
 & PCVT,PCVL,PCVH, KVEG,KTVL,KTVH, &
 & PTSKM1M,PTSAM1M, &
 & PANFMVT, PANDAYVT, &
 & PLAIVT,PLAIL,PLAIH, &
 & PLAILC,  PLAIHC, &
 & PRESPBSTR,PRESPBSTR2,PBIOMASS_LAST,PBIOMASSTR_LAST,PBIOMASSTR2_LAST, &
 & PBLOSSVT, PBGAINVT, &
 & PLAIE1, PBSTRE1, PBSTR2E1, &
 & PLAI, PBIOM, PBLOSS, PBGAIN, PBIOMSTR, PBIOMSTR2, &
 & PDHBIOS,PDHVEGS, &
 & YDCST,YDVEG,YDAGS)

ELSE
  ! Compute PLAI if this routine is skiped
  DO JL=KIDIA,KFDIA
    PLAI(JL)=PCVL(JL)*PLAIL(JL)+PCVH(JL)*PLAIH(JL)
  ENDDO
  PBIOM(KIDIA:KFDIA)=0._JPRB
  PBLOSS(KIDIA:KFDIA)=0._JPRB
  PBGAIN(KIDIA:KFDIA)=0._JPRB
  PBIOMSTR(KIDIA:KFDIA)=0._JPRB
  PBIOMSTR2(KIDIA:KFDIA)=0._JPRB
  IF (SIZE(PDHBIOS) > 0) THEN
    PDHBIOS(KIDIA:KFDIA,1:YDVEG%NVTILES,1:5)=0._JPRB
  ENDIF
ENDIF

!    ------------------------------------------------------------------
!*   7. Budget computations : TESTING ONLY !!
!*   Should be comment out when running coupled to IFS !
!    ------------------------------------------------------------------

IF (YDSOIL%LEWBCHECK) THEN
  CALL BUDGET_MASS_DDH
ENDIF



END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFTSTP_CTL_MOD:SURFTSTP_CTL',1,ZHOOK_HANDLE)

CONTAINS
SUBROUTINE BUDGET_MASS_DDH
IMPLICIT NONE
INTEGER(KIND=JPRB) :: JL
REAL(KIND=JPRB)    :: ZFLUX,ZSTORAGE,ZRES,ZSNE,ZSOILE,ZEVAP,ZLAKEE
REAL(KIND=JPRB)    :: ZWTRH,ZEPSILON,ZECanopMISS
REAL(KIND=JPRB)    :: ZSOILEL(KLEVS)

ZEPSILON=100._JPRB*EPSILON(ZEPSILON)
ZWTRH=MAX(1e-10_JPRB,ZEPSILON) !10  ! snow mip 1e-7 !
DO JL=KIDIA,KFDIA
!   IF (LDLAND(JL) .AND. (PFRTI(JL,9) == 0._JPRB) .AND. PFRTI(JL,1) == 0._JPRB ) THEN
  IF (LDLAND(JL) ) THEN
    ZSTORAGE=PWLE1(JL)
    ZSNE=0._JPRB
    DO JK=1,KLEVSN
      ZSNE=ZSNE+PSNE1(JL,JK)
    ENDDO
    ZSOILE=0._JPRB
    ZSOILEL(:)=0._JPRB
    DO JK=1,KLEVS
      ZSOILEL(JK)=PWSAE1(JL,JK)*RHOH2O*YDSOIL%RDAW(JK)
      ZSOILE=ZSOILE+ZSOILEL(JK)
    ENDDO
    ZLAKEE=PFRTI(JL,9)*PEVAPTI(JL,9)+PFRTI(JL,1)*PEVAPTI(JL,1)
    ZSTORAGE=PWLE1(JL)+ZSNE+ZSOILE+ZLAKEE
    ZECanopMISS=(ZEINTTI(JL,3)+ZEINTTI(JL,4)+&
      & ZEINTTI(JL,6)+ZEINTTI(JL,7)+ZEINTTI(JL,8) )
    ZEVAP=0._JPRB
    DO JT=1,KTILES
      ZEVAP=ZEVAP+PFRTI(JL,JT)*PEVAPTI(JL,JT)
    ENDDO
    ZEVAP=ZEVAP !+ ZEINTTI(JL,3)-MAX(0._JPRB,(PFRTI(JL,3)*PEVAPTI(JL,3)))
    ZFLUX=PRSFC(JL)+PRSFL(JL)+PSSFC(JL)+PSSFL(JL)-PROFS(JL)-PROFD(JL)+ZEVAP
    ZRES=ZSTORAGE-ZFLUX
    ZWTRH=MAX(ZEPSILON,MAX(ABS(ZSTORAGE),ABS(ZFLUX))*ZEPSILON)

    if (ZEPSILON .ne. ZWTRH ) then
      print*,ZEPSILON, ZWTRH
      read(*,*)
    endif
    IF ( ABS(ZRES) > ZWTRH ) THEN
      write(*,'("DDH TWBAL:storage,flux,residual,threshold,EPSILON (kg/m2/s):",1X,I10,5(1X,E14.6E3))') &
            JL,ZSTORAGE,ZFLUX,ZRES,ZWTRH,ZEPSILON
      write(*,'("DDH WBAL1:Int,t,t+1 tend",1X,2(1X,E14.6E3))') PWLM1M(JL),ZWL(JL)
      write(*,'("DDH WBAL1:Int. tend",1X,1(1X,E14.6E3))') PWLE1(JL)
      write(*,'("DDH WBAL1:Snow t,t+1",1X,2(1X,E14.6E3))')SUM(PSNM1M(JL,:)),SUM(ZSN(JL,:))
      write(*,'("DDH WBAL1:Snow tend",1X,1(1X,E14.6E3))') ZSNE
      write(*,'("DDH WBAL1:soil t,t+1",1X,8(1X,E14.6E3))') PWSAM1M(JL,:)*RHOH2O*YDSOIL%RDAW(:),ZWSA(JL,:)*RHOH2O*YDSOIL%RDAW(:)
      write(*,'("DDH WBAL1:soil tend",1X,1(1X,E14.6E3))') ZSOILE
      write(*,'("DDH WBAL1:soil tend Levels",1X,4(1X,E14.6E3))') ZSOILEL(:)
      write(*,'("DDH WBAL1:soil flux ZSAWGFL",1X,4(1X,E14.6E3))')ZSAWGFL(JL,:)
      write(*,'("DDH WBAL1:Lake tend",1X,1(1X,E14.6E3))') ZLAKEE
      write(*,'("DDH WBAL1:total rainf",1X,1(1X,E14.6E3))') PRSFC(JL)+PRSFL(JL)
      write(*,'("DDH WBAL1:total THROUGHFALL",1X,1(1X,E14.6E3))') ZTSFC(JL)+ZTSFL(JL)
      write(*,'("DDH WBAL1:PEINTTI",1X,9(1X,E14.6E3))')ZEINTTI(JL,:)
      write(*,'("DDH WBAL1:total snowfall",1X,1(1X,E14.6E3))') PSSFC(JL)+PSSFL(JL)
      write(*,'("DDH WBAL1:surf ro",1X,1(1X,E14.6E3))') -PROFS(JL)
      write(*,'("DDH WBAL1:sub surf ro",1X,1(1X,E14.6E3))') -PROFD(JL)
      write(*,'("DDH WBAL1:tot evap",1X,1(1X,E14.6E3))') ZEVAP
      write(*,'("DDH WBAL1:evap snow",1X,1(1X,E14.6E3))') PEVAPSNW(JL)*PFRTI(JL,7)
      write(*,'("DDH WBAL1:ZECanopMISS,",1X,1(1X,E14.6E3))')ZECanopMISS
      write(*,'("DDH WBAL1:ZFWEL1,ZFWE234,",1X,2(1X,E14.6E3))') ZFWEL1(JL),ZFWE234(JL)
      write(*,'("DDH WBAL1:ZFWEL1-E3&8,",1X,2(1X,E14.6E3))') ZFWEL1(JL)-(PFRTI(JL,3)*PEVAPTI(JL,3)+PFRTI(JL,8)*PEVAPTI(JL,8))
      write(*,'("DDH WBAL1:ZFWEL234-E467,",1X,2(1X,E14.6E3))') ZFWE234(JL)-(PFRTI(JL,4)*PEVAPTI(JL,4)+PFRTI(JL,6) &
          & *PEVAPTI(JL,6)+PFRTI(JL,7)*PEVAPTI(JL,7))
      write(*,'("DDH WBAL1:tile evap",1X,9(1X,E14.6E3))') PFRTI(JL,:)*PEVAPTI(JL,:)
      write(*,'("DDH WBAL1:tile fractions",1X,9(1X,E14.6E3))') PFRTI(JL,:)
      write(*,*)'----END---'
      write(*,*)''
      IF ( YDSOIL%LEWBCHECKAbort ) THEN
        CALL ABORT_SURF('ALL WATER BALANCE CLOSURE ERROR')
      ENDIF
    ENDIF



  ENDIF
ENDDO

END SUBROUTINE BUDGET_MASS_DDH

END SUBROUTINE SURFTSTP_CTL
END MODULE SURFTSTP_CTL_MOD
