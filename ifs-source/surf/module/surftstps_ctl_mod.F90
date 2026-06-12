MODULE SURFTSTPS_CTL_MOD
CONTAINS
SUBROUTINE SURFTSTPS_CTL(KIDIA , KFDIA , KLON  , KLEVS , KLEVSN, &
 & KTILES, KSOTY,&
 & PTSPHY , PSDOR, PFRTI,&
 & PAHFSTI, PEVAPTI, PSSRFLTI,&
 & LDLAND,  LDSICE, LDSI, LDNH,&
 & PSNM1M  ,PTSNM1M,PRSNM1M,&
 & PTSAM1M, PTIAM1M,&
 & PWLM1M  ,PWSAM1M,&
 & PHLICEM1M, PAPRS,&
 & PRSFC   ,PRSFL,&
 & PSLRFL  ,PSSFC  ,PSSFL,&
 & PCVL    ,PCVH   ,PWLMX   ,PEVAPSNW,&
 & YDCST   ,YDVEG  ,YDSOIL  ,YDFLAKE,&
!-TENDENCIES OUTPUT
 & PTSNE1 , PTSAE1 , PTIAE1)

!-END OF CALL SURFTSTPS

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST   ,ONLY : TCST
USE YOS_VEG   ,ONLY : TVEG
USE YOS_SOIL  ,ONLY : TSOIL
USE YOS_FLAKE ,ONLY : TFLAKE

USE SRFSN_LWIMPS_MOD
USE SRFSN_LWIMPMLS_MOD
USE SRFSN_VGRID_MOD 
USE SRFSN_REGRID_MOD
USE SRFSN_SSRABSS_MOD
USE SRFSN_WEBALS_MOD
USE SRFTS_MOD
USE SRFIS_MOD
USE SRFWLS_MOD
USE SRFRCGS_MOD

#ifdef DOC

! (C) Copyright 2011- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SURFTSTPS* - UPDATES LAND VALUES OF TEMPERATURE AND SNOW.

!     PURPOSE.
!     --------
!          This routine updates the sea ice values of temperature
!     and the land values of soil temperature and snow temperature. 
!     For temperature, this is done via a forward time
!     step damped with some implicit linear considerations: as if all
!     fluxes that explicitely depend on the variable had only a linear
!     variation around the t-1 value. 

!**   INTERFACE.
!     ----------
!          *SURFTSTPS* IS CALLED FROM *CALLPAR*.
!          THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE:
!     TSA,WL,SN AT T-1,TSK,SURFACE FLUXES COMPUTED IN OTHER PARTS OF
!     THE PHYSICS, AND W AND LAND-SEA MASK. IT RETURNS AS AN OUTPUT
!     TENDENCIES TO THE SAME VARIABLES (TSA,WL,SN).


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
!    *KSOTY*      SOIL TYPE                                   (1-7)

!     INPUT PARAMETERS (REAL):
!    *PTSPHY*     TIME STEP FOR THE PHYSICS                      S
!    *PFRTI*      TILE FRACTIONS                              (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!    *PAHFSTI*      SURFACE SENSIBLE HEAT FLUX, FOR EACH TILE  W/M2
!    *PEVAPTI*      SURFACE MOISTURE FLUX, FOR EACH TILE      KG/M2/S
!    *PSSRFLTI*     NET SHORTWAVE RADIATION FLUX AT SURFACE, FOR
!                       EACH TILE                                 W/M2

!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE)
!    *LDSICE*     SEA ICE MASK (.T. OVER SEA ICE)
!    *LDSI*       TRUE IF THERMALLY RESPONSIVE SEA-ICE
!    *LDNH*       TRUE FOR NORTHERN HEMISPHERE LATITUDE ROW

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!    *PSNM1M*     SNOW MASS (per unit area)                      KG/M**2
!    *PRSNM1M*    SNOW DENSITY                                 KG/M3
!    *PTSNM1M*    SNOW TEMPERATURE                               K
!    *PTSAM1M*    SOIL TEMPERATURE                               K
!    *PTIAM1M*    SEA-ICE TEMPERATURE                            K
!          (NB: REPRESENTS THE FRACTION OF ICE ONLY, NOT THE GRID-BOX)
!    *PWLM1M*     SKIN RESERVOIR WATER CONTENT                   kg/m**2
!    *PWSAM1M*    SOIL MOISTURE                                m**3/m**3
!    *PHLICEM1M*  LAKE ICE THICKNESS                             m
!    *PRSFC*      CONVECTIVE RAIN FLUX AT THE SURFACE          KG/M**2/S
!    *PRSFL*      LARGE SCALE RAIN FLUX AT THE SURFACE         KG/M**2/S
!    *PSLRFL*     NET LONGWAVE  RADIATION AT THE SURFACE         W/M**2
!    *PSSFC*      CONVECTIVE  SNOW FLUX AT THE SURFACE         KG/M**2/S
!    *PSSFL*      LARGE SCALE SNOW FLUX AT THE SURFACE         KG/M**2/S
!    *PCVL*       LOW VEGETATION COVER  (CORRECTED)              (0-1)
!    *PCVH*       HIGH VEGETATION COVER (CORRECTED)              (0-1)
!    *PWLMX*      MAXIMUM SKIN RESERVOIR CAPACITY                kg/m**2
!    *PEVAPSNW*   EVAPORATION FROM SNOW UNDER FOREST           KG/M**2/S


!     OUTPUT PARAMETERS (TENDENCIES):
!    *PTSNE1*     SNOW TEMPERATURE                              K/S
!    *PTSAE1*     SOIL TEMPERATURE TENDENCY                     K/S
!    *PTIAE1*     SEA-ICE TEMPERATURE TENDENCY                  K/S

!     METHOD.
!     -------
!          STRAIGHTFORWARD ONCE THE DEFINITION OF THE CONSTANTS IS
!     UNDERSTOOD. FOR THIS REFER TO DOCUMENTATION. FOR THE TIME FILTER
!     SEE CORRESPONDING PART OF THE DOCUMENTATION OF THE ADIABATIC CODE.

!     EXTERNALS.
!     ----------
!          *SRFWLS*        COMPUTES THE SKIN RESERVOIR CHANGES.
!          *SRFSN_LWIMPS*  REVISED SNOW SCHEME W. DIAG. LIQ. WATER.
!          *SRFRCGS*       COMPUTE SOIL HEAT CAPACITY.
!          *SRFTS*         COMPUTES THE TEMPERATURE CHANGES BEFORE THE SNOW
!                            MELTING.
!          *SRFIS*         COMPUTES TEMPERATURE EVOLUTION OF SEA ICE.

!     REFERENCE.
!     ----------
!          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
!     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.

!     Original   
!     --------
!          Simplified version based on SURFTSTP
!     M. Janiskova              E.C.M.W.F.     27-07-2011  

!     Modifications
!     -------------

!     ------------------------------------------------------------------
#endif

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVSN
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAHFSTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVAPTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFLTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTY(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSDOR(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:)
LOGICAL           ,INTENT(IN)    :: LDLAND(:)
LOGICAL           ,INTENT(IN)    :: LDSICE(:)
LOGICAL           ,INTENT(IN)    :: LDSI
LOGICAL           ,INTENT(IN)    :: LDNH(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNM1M(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSNM1M(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSNM1M(:,:)
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWLMX(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVAPSNW(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICEM1M(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(:)
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TVEG)        ,INTENT(IN)    :: YDVEG
TYPE(TSOIL)       ,INTENT(IN)    :: YDSOIL
TYPE(TFLAKE)      ,INTENT(IN)    :: YDFLAKE

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSNE1(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSAE1(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTIAE1(:,:)

!*         0.2    DECLARATION OF LOCAL VARIABLES.
!

REAL(KIND=JPRB) :: ZTSN(KLON,YDSOIL%NLEVSN)
REAL(KIND=JPRB) :: ZTSA(KLON,KLEVS),ZTIA(KLON,KLEVS)
REAL(KIND=JPRB) :: ZCTSA(KLON,KLEVS)
REAL(KIND=JPRB) :: ZTSFC(KLON),ZTSFL(KLON)
REAL(KIND=JPRB) :: ZRSFC(KLON),ZRSFL(KLON),ZSSFC(KLON),ZSSFL(KLON)
REAL(KIND=JPRB) :: ZGSN(KLON)

REAL(KIND=JPRB) :: ZSNOTRS(KLON, YDSOIL%NLEVSN+1) 
REAL(KIND=JPRB) :: ZWSNM1M(KLON, YDSOIL%NLEVSN) 

LOGICAL         :: LLNOSNOW(KLON)
REAL(KIND=JPRB) :: ZDSNOUT(KLON,YDSOIL%NLEVSN)
REAL(KIND=JPRB) :: ZSSNM1M_REG(KLON,YDSOIL%NLEVSN),&
                 & ZRSNM1M_REG(KLON,YDSOIL%NLEVSN),&
                 & ZWSNM1M_REG(KLON,YDSOIL%NLEVSN),&
                 & ZTSNM1M_REG(KLON,YDSOIL%NLEVSN)
INTEGER(KIND=JPIM) :: KLEVSNA(KLON)
REAL(KIND=JPRB) :: ZFRSN

REAL(KIND=JPRB) :: ZTSPHY

INTEGER(KIND=JPIM) :: JK, JL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*         1.1   SET UP SOME CONSTANTS, INITIALISE ARRAYS.
!                --- -- ---- -----------------------------
IF (LHOOK) CALL DR_HOOK('SURFTSTPS_CTL_MOD:SURFTSTPS_CTL',0,ZHOOK_HANDLE)

ZTSPHY=1.0_JPRB/PTSPHY

DO JL=KIDIA,KFDIA
! protect against negative precip

  ZRSFC(JL)=MAX(0.0_JPRB,PRSFC(JL))
  ZRSFL(JL)=MAX(0.0_JPRB,PRSFL(JL))
  ZSSFC(JL)=MAX(0.0_JPRB,PSSFC(JL))
  ZSSFL(JL)=MAX(0.0_JPRB,PSSFL(JL))
ENDDO

!     ------------------------------------------------------------------
!*         0.1     INTERCEPTION RESERVOIR AND ASSOCIATED RUNOFF CHANGES.

CALL SRFWLS(KIDIA,KFDIA,KLON,&
 & PTSPHY,PWLM1M,PCVL,PCVH,PWLMX,&
 & PFRTI  , PEVAPTI ,ZRSFC,ZRSFL,&
 & LDLAND ,&
 & YDSOIL ,YDVEG,&
 & ZTSFC  ,ZTSFL)

!     ------------------------------------------------------------------

!*         2.     SNOW CHANGES.

!* LIQUID WATER DIAGNOSTIC  (LESN09)
IF (YDSOIL%LESNML) THEN
!!CALL SRFSN_LWIMPMLS(KIDIA  ,KFDIA  ,KLON   ,PTSPHY, &
!! & PSNM1M ,PTSNM1M ,PRSNM1M,PTSAM1M,PHLICEM1M,    &
!! & PSLRFL ,PSSRFLTI,PFRTI  ,PAHFSTI,PEVAPTI,      &
!! & ZSSFC  ,ZSSFL   ,PEVAPSNW,                     &
!! & ZTSFC  ,ZTSFL   ,                              &
!! & PAPRS   ,PWSAM1M , KSOTY   ,                   &
!! & YDCST  ,YDVEG   ,YDSOIL ,YDFLAKE ,             &
!! & ZTSN   ,ZGSN )

! Initialise some snowML related working arrays
  DO JL=KIDIA,KFDIA
    ZFRSN=MAX(PFRTI(JL,5)+PFRTI(JL,7),YDSOIL%RFRTINY)
    IF (ZFRSN < YDSOIL%RFRSMALL) THEN
      LLNOSNOW(JL)=.TRUE.
    ELSE
      LLNOSNOW(JL)=.FALSE.
    ENDIF
  ENDDO
  ZWSNM1M(KIDIA:KFDIA, 1:YDSOIL%NLEVSN)     = 0._JPRB
  ZSSNM1M_REG(KIDIA:KFDIA, 1:YDSOIL%NLEVSN) = 0._JPRB 
  ZRSNM1M_REG(KIDIA:KFDIA, 1:YDSOIL%NLEVSN) = 0._JPRB
  ZWSNM1M_REG(KIDIA:KFDIA, 1:YDSOIL%NLEVSN) = 0._JPRB
  ZTSNM1M_REG(KIDIA:KFDIA, 1:YDSOIL%NLEVSN) = 0._JPRB
  ZSNOTRS(KIDIA:KFDIA, 1:YDSOIL%NLEVSN+1) = 0._JPRB
!----
! Note: Steps 2.1 and 2.2 are required only at time step 0 to re-balance 
! the input profiles, as mass budget is not computed in simplified phys.
! In TL/AD these steps are performed during the reading of trajectory fields
! in read_surfgrid_traj so in principle their use in tl/ad is not required 
!----
! 2.1 compute vertical grid from previous time step profiles
  CALL SRFSN_VGRID(KIDIA,KFDIA,KLON,YDSOIL%NLEVSN,LLNOSNOW,&
              &  PSDOR, &
              &  PSNM1M,PRSNM1M,&
              &  YDSOIL%RLEVSNMIN,YDSOIL%RLEVSNMAX,&
              &  YDSOIL%RLEVSNMIN_GL,YDSOIL%RLEVSNMAX_GL,&
              &  ZDSNOUT,KLEVSNA, .TRUE.)

! 2.2 Regrid snowML fields based on vertical grid
  CALL SRFSN_REGRID(KIDIA,KFDIA,KLON,YDSOIL%NLEVSN,LLNOSNOW,&
                  & PSNM1M,PRSNM1M,ZWSNM1M,PTSNM1M,ZDSNOUT,&
                  & ZSSNM1M_REG,ZRSNM1M_REG,ZWSNM1M_REG,ZTSNM1M_REG)

! 2.3 Compute absorption of solar radiation into the snowpack
  ZSNOTRS(KIDIA:KFDIA, 1:YDSOIL%NLEVSN+1)   = 0._JPRB
  CALL SRFSN_SSRABSS(KIDIA,KFDIA,KLON,YDSOIL%NLEVSN,&
                  & LLNOSNOW,PFRTI,PSSRFLTI,&
                  & ZSSNM1M_REG,ZRSNM1M_REG,&
                  & YDSOIL,YDCST,&
                  & ZSNOTRS)

! 2.4 Compute energy balance of the snowpack
  CALL SRFSN_WEBALS(KIDIA,KFDIA,KLON,YDSOIL%NLEVSN,&
   & PTSPHY,PFRTI,LDLAND, LLNOSNOW,  &
   & ZSSNM1M_REG,ZWSNM1M_REG,ZRSNM1M_REG,ZTSNM1M_REG,&
   & PSLRFL,PSSRFLTI,PAHFSTI ,PEVAPTI,&
   & ZSSFC, ZSSFL, PEVAPSNW, ZTSFC  ,ZTSFL   ,&
   & PTSAM1M(KIDIA:KFDIA,1), ZSNOTRS,         &
   & PAPRS, PWSAM1M(KIDIA:KFDIA,1), KSOTY,&
   & YDSOIL,YDCST,&
   & ZTSN, ZGSN)

 
ELSE

  CALL SRFSN_LWIMPS(KIDIA  ,KFDIA  ,KLON   ,PTSPHY, &
   & PSNM1M(KIDIA:KFDIA,1) ,PTSNM1M(KIDIA:KFDIA,1) ,PRSNM1M(KIDIA:KFDIA,1),&
   & PTSAM1M,PHLICEM1M,    &
   & PSLRFL ,PSSRFLTI,PFRTI  ,PAHFSTI,PEVAPTI,      &
   & ZSSFC  ,ZSSFL   ,PEVAPSNW,                     &
   & ZTSFC  ,ZTSFL   ,                              &
   & YDCST  ,YDVEG   ,YDSOIL ,YDFLAKE ,             &
   & ZTSN(KIDIA:KFDIA,1)   ,ZGSN )
  
ENDIF
!     ------------------------------------------------------------------
!*         3.     SOIL HEAT BUDGET.
!                 ---- ---- -------

!*         3.1     COMPUTE SOIL HEAT CAPACITY.

CALL SRFRCGS(KIDIA  , KFDIA  , KLON, KLEVS ,&
 & LDLAND , LDSICE ,&
 & PTSAM1M, KSOTY  , PCVL   , PCVH   ,&
 & YDCST  , YDSOIL ,&
 & ZCTSA)

!*         3.2     TEMPERATURE CHANGES.

! soil
CALL SRFTS(KIDIA  , KFDIA  , KLON   , KLEVS  , &
 & PTSPHY , PTSAM1M, PWSAM1M, KSOTY, &
 & PFRTI  , PAHFSTI,PEVAPTI ,&
 & PSLRFL ,PSSRFLTI, ZGSN   ,&
 & ZCTSA  , ZTSA   , LDLAND ,&
 & YDCST  , YDSOIL , YDFLAKE)  

! sea-ice
IF (LDSI) THEN
  CALL SRFIS(&
   & KIDIA  , KFDIA  , KLON   , KLEVS  , &
   & PTSPHY , PTIAM1M, PAHFSTI, PEVAPTI, &
   & PSLRFL ,PSSRFLTI, ZTIA   , LDSICE , LDNH, &
   & YDCST  , YDSOIL)
ELSE
  DO JK=1,KLEVS
    DO JL=KIDIA,KFDIA
      ZTIA(JL,JK)=PTIAM1M(JL,JK)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------
!*         5.     COMPUTE TENDENCIES (TEND. FOR TSK COMP. IN VDIFF)
!                 ------- ---------- ------------------------------
!*         5.1   TENDENCIES
DO JK=1,YDSOIL%NLEVSN
  DO JL=KIDIA,KFDIA
      PTSNE1(JL,JK) =PTSNE1(JL,JK)+&
     & (ZTSN(JL,JK)-PTSNM1M(JL,JK))*ZTSPHY
  ENDDO
ENDDO

DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA
    PTIAE1(JL,JK)=PTIAE1(JL,JK)+&
     & (ZTIA(JL,JK)-PTIAM1M(JL,JK))*ZTSPHY
  ENDDO
ENDDO
! PTSA represents the grid-box average T; the sea-ice fraction
!    has a T representative of the fraction
DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA
    PTSAE1(JL,JK)=PTSAE1(JL,JK)+&
     & (PFRTI(JL,2)*&
     & (ZTIA(JL,JK)-PTIAM1M(JL,JK))+&
     & (1.0_JPRB-PFRTI(JL,1)-PFRTI(JL,2))*&
     & (ZTSA(JL,JK)-PTSAM1M(JL,JK)))*ZTSPHY
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('SURFTSTPS_CTL_MOD:SURFTSTPS_CTL',1,ZHOOK_HANDLE)

END SUBROUTINE SURFTSTPS_CTL
END MODULE SURFTSTPS_CTL_MOD
