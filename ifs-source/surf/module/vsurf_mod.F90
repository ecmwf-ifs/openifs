MODULE VSURF_MOD
CONTAINS
SUBROUTINE VSURF(KIDIA,KFDIA,KLON,KTILES,KLEVS,KTILE,&
 & KTVL, KCO2TYP, KTVH,&
 & JVTTL,KVEG, &
 & PLAI,&
 & PMU0,PCO2FLUX, &
 & PFRTI, PLAIL, PLAIH,&
 & PTMLEV, PQMLEV  , PCMLEV, PAPHMS,&
 & PTSKM1M,PWSAM1M,PTSAM1M,KSOTY,&
 & PSRFD ,PRAQ  ,PQSAM ,&
 & PQS   ,PDQS  ,&
 & PWETB ,PCPTS ,PWETL, PWETLU, PWETH, PWETHS , &  
 & PEVAP ,&
 & PAN,PAG,PRD ,PPWLIQ ,&
 & PDHVEGS, PEXDIAG, &
 & YDCST, YDVEG, YDEXC, YDAGS, YDAGF, YDSOIL, YDFLAKE, YDURB)
  
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_THF  , ONLY : R4LES, R5LES, R2ES, R4IES, R3LES, R3IES, R5IES, RVTMP2
USE YOS_CST  , ONLY : TCST
USE YOS_VEG  , ONLY : TVEG
USE YOS_EXC  , ONLY : TEXC
USE YOS_AGS  , ONLY : TAGS
USE YOS_AGF  , ONLY : TAGF
USE YOS_SOIL , ONLY : TSOIL
USE YOS_FLAKE, ONLY : TFLAKE
USE YOS_URB  , ONLY : TURB
USE COTWORESTRESS_MOD

! (C) Copyright 1990- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!     ------------------------------------------------------------------

!**   *VSURF* - PREPARES SURFACE BOUNDARY CONDITION FOR T AND Q

!     DERIVED FROM VDIFF (CY34) BY
!     A.C.M. BELJAARS       E.C.M.W.F.    18-1-90
!     Modified P.VITERBO AND A.C.M. BELJAARS  E.C.M.W.F.    16-3-93
!     Modified ACM Beljaars  26-03-99  Tiling of the surface
!     P. Viterbo     24-05-2004     Change surface units
!     P. Viterbo ECMWF 12/05/2005 Externalize SURF
!                     (based on VDFSURF)
!     G. Balsamo ECMWF 22/05/2006   Evaporative fraction f(soil)
!     G. Balsamo ECMWF 03/07/2006   Add soil type
!     E. Dutra/G. Balsamo 01/05/2008   Add lake tile
!     G. Balsamo ECMWF 9/3/2010 Bare ground evaporation
!     S. Boussetta/G.Balsamo May 2009 Add lai
!     S. Boussetta/G.Balsamo May 2010 Add CTESSEL based on:
!     M.H. Voogt (KNMI) "C-Tessel"  09/2005 
!     S. Lafont "C-TESSEL" 18/05/2006
!     S. Boussetta/G.Balsamo June 2010 Add soil moisture scaling factor for Reco
!     S. Boussetta/G.Balsamo June 2011 modularisation of Ags call
!     A. Beljaars      26/02/2014   compute unstressed evaporation
!     A. Agusti-Panareda Nov 2020  couple atm CO2 tracer with photosynthesis 
!     A. Agusti-Panareda May 2021  Pass soil temperature to photosynthesis 
!     A. Agusti-Panareda June 2021 Pass photosynthetic pathway for low vegetation (c3/c4)

!     PURPOSE
!     -------

!     PREPARE SURFACE BOUNDARY CONDITION FOR Q AND T, E.G. FRACTIONAL
!     SURFACE COVER (SNOW AND VEGETATION), SATURATION SPECIFIC HUMIDITY
!     AT THE SURFACE, RELATIVE HUMIDITY OVER BARE LAND AND THE STOMATAL
!     RESISTANCE.

!     INTERFACE
!     ---------

!     *VSURF* IS CALLED BY *SURFEXCDRIVER*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KTILES*       NUMBER OF TILES
!     *KLEVS*        NUMBER OF SOIL LAYERS
!     *KTILE*        TILE INDEX
!     *KTVL*         VEGETATION TYPE FOR LOW VEGETATION FRACTION
!     *KCO2TYP*     TYPE OF PHOTOSYNTHETIC PATHWAY FOR LOW VEGETATION (C3/C4)
!     *KTVH*         VEGETATION TYPE FOR HIGH VEGETATION FRACTION

!     *JVTTL*        CROSS-REFERENCE BETWEEN TILES AND VEGETATION TYPES 


!     *KSOTY*        SOIL TYPE (1-7) 

!     INPUT PARAMETERS (REAL):

!     *PFRTI*      TILE FRACTIONS                                   (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!            9 : LAKE                  10 : URBAN
!     *PLAIL*        LAI OF LOW VEGETATION
!     *PLAIH*        LAI OF HIGH VEGETATION

!     *PLAI*         LEAF AREA INDEX for both low and high                 (-)
!     *PMU0*        LOCAL COSINE OF INSTANTANEOUS MEAN SOLAR ZENITH ANGLE

!     *PTMLEV*      TEMPERATURE AT T-1, lowest model level
!     *PQMLEV*      SPECIFIC HUMIDITY AT T-1, lowest model level
!     *PCMLEV*      ATMOSPHERIC CO2 AT T-1, lowest model level
!     *PAPHMS*      PRESSURE AT T-1, surface
!     *PTSKM1M*      SURFACE TEMPERATURE
!     *PWSAM1M*      SOIL MOISTURE ALL LAYERS                   M**3/M**3
!     *PTSAM1M*      SOIL TEMPERATURE ALL LAYERS  
!     *PSRFD*        DOWNWARD SHORT WAVE RADIATION FLUX AT SURFACE
!     *PRAQ*         PRELIMINARY AERODYNAMIC RESISTANCE

!     OUTPUT PARAMETERS (REAL):

!     *PQSAM*        SPECIFIC HUMIDITY AT THE SURFACE
!     *PQS*          SATURATION Q AT SURFACE
!     *PDQS*         DERIVATIVE OF SATURATION Q-CURVE AT SURFACE T
!     *PWETB*        BARE SOIL RESISTANCE
!     *PCPTS*        DRY STATIC ENRGY AT SURFACE
!     *PWETL*        CANOPY RESISTANCE LOW VEGETATION
!     *PWETLU*       CANOPY RESISTANCE LOW VEGETATION IN UNSTRESSED CONDITIONS
!     *PWETH*        CANOPY RESISTANCE HIGH VEGETATION, SNOW FREE
!     *PWETHS*       CANOPY RESISTANCE HIGH VEGETATION WITH SNOW
!     *PEXDIAG*      Extra fields for potential pp of canopy


!     *PDHVEGS*      Diagnostic array for vegetation (see module yomcdh)
!     *PAN*          NET CO2 ASSIMILATION OVER CANOPY          KG_CO2/M2/S
!                    positive downwards, to be changed for diagnostic output
!     *PAG*          GROSS CO2 ASSIMILATION OVER CANOPY        KG_CO2/M2/S
!                    positive downwards, to be changed for diagnostic output
!     *PRD*          DARK RESPIRATION                          KG_CO2/M2/S
!                    positive upwards

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILE 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCO2TYP(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(:)

INTEGER(KIND=JPIM),INTENT(IN)    :: JVTTL
INTEGER(KIND=JPIM),INTENT(IN)    :: KVEG(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAI(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVAP(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAN(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAG(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRD(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHVEGS(:,:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXDIAG(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCO2FLUX(:)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAIL(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAIH(:)
INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTY(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSAM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSAM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSRFD(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRAQ(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQSAM(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQS(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDQS(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWETB(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCPTS(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWETL(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWETLU(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWETH(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWETHS(:) 

!ZLIQ is passed to compute soil moisture scaling factor in Reco (CO2 routine) CTESSEL
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPWLIQ(KLON,KLEVS)

TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TVEG)        ,INTENT(IN)    :: YDVEG
TYPE(TEXC)        ,INTENT(IN)    :: YDEXC
TYPE(TAGS)        ,INTENT(IN)    :: YDAGS
TYPE(TAGF)        ,INTENT(IN)    :: YDAGF
TYPE(TSOIL)       ,INTENT(IN)    :: YDSOIL
TYPE(TFLAKE)      ,INTENT(IN)    :: YDFLAKE
TYPE(TURB)        ,INTENT(IN)    :: YDURB

!*    LOCAL STORAGE
!     ----- -------

INTEGER(KIND=JPIM) :: JK, JL, JS
REAL(KIND=JPRB) ::  ZLIQ(KLON,KLEVS), ZLIQR(KLON,KLEVS),ZF2(KLON),ZWROOT(KLON)
REAL(KIND=JPRB) ::  ZDSP(KLON),ZDMAXT(KLON)
REAL(KIND=JPRB) ::  ZWET(KLON)
REAL(KIND=JPRB) ::  ZTSK(KLON)
REAL(KIND=JPRB) ::  ZEPSILON

REAL(KIND=JPRB) ::  ZTSOIL(KLON)
REAL(KIND=JPRB) ::  ZCOR, ZEPSF3, ZF, ZF1H, ZF1L, ZF2H, ZF2L, ZF2B, &
 & ZF3H, ZF3L, ZHSTRH, ZHSTRL, ZLAIH, ZLAIL, ZQSAIR, &
 & ZRSMINH, ZRSMINL, ZRSMINB, ZRVA, ZRVB, ZSRFL, ZWROOTH, ZWROOTL, &
 & ZQWEVAP, ZWPWP, ZQWEVAPBARE, ZBARE, ZWPBARE,&
 & ZSALIN, ZRCLU
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
#include "fcsttre.h"
!     ------------------------------------------------------------------

!*       1.     INITIALIZE CONSTANTS
!               ---------- ----------

IF (LHOOK) CALL DR_HOOK('VSURF_MOD:VSURF',0,ZHOOK_HANDLE)
ASSOCIATE(RCPD=>YDCST%RCPD, RETV=>YDCST%RETV, RLSTT=>YDCST%RLSTT, &
 & RLVTT=>YDCST%RLVTT, RTT=>YDCST%RTT, &
 & LEOCSA=>YDEXC%LEOCSA, &
 & LEVGEN=>YDSOIL%LEVGEN, RQWEVAP=>YDSOIL%RQWEVAP, RQWEVAPM=>YDSOIL%RQWEVAPM, &
 & RTF1=>YDSOIL%RTF1, RTF2=>YDSOIL%RTF2, RTF3=>YDSOIL%RTF3, RTF4=>YDSOIL%RTF4, &
 & RWCAP=>YDSOIL%RWCAP, RWCAPM=>YDSOIL%RWCAPM, RWPWP=>YDSOIL%RWPWP, &
 & RWPWPM=>YDSOIL%RWPWPM, RWRESTM=>YDSOIL%RWRESTM, &
 & LEAGS=>YDVEG%LEAGS, LECTESSEL=>YDVEG%LECTESSEL, RCEPSW=>YDVEG%RCEPSW, &
 & LEFARQUHAR=>YDVEG%LEFARQUHAR, LEAIRCO2COUP=>YDVEG%LEAIRCO2COUP, &
 & RVHSTR=>YDVEG%RVHSTR, RVLAI=>YDVEG%RVLAI, RVROOTSA=>YDVEG%RVROOTSA, &
 & RVRSMIN=>YDVEG%RVRSMIN, RURBRES=>YDURB%RURBRES)


! This is needed as unitialized values end up being passed around otherwise
ZDSP(:)=0.0_JPRB
ZDMAXT(:)=0.0_JPRB

ZRVA=5000._JPRB
ZRVB=10._JPRB
ZEPSF3=0.00001_JPRB ! security value for exponential sat-deficit dependence
ZRSMINB=50._JPRB  ! bare soil minimum resistance
ZRCLU=50._JPRB    ! unstressed canopy resis. for pot. evap. over grassland (50 s/m)

ZEPSILON=EPSILON(ZEPSILON)
!     ------------------------------------------------------------------


!          2.    PREPARE SURFACE BOUNDARY CONDITION
!                ------- ------- -------- ---------

!*         2.1   RELATIVE HUMIDITY OVER THE BARE LAND PART

!                BARE SOIL RESISTANCE IS COMPUTED FOR KTILE=4

!*         2.2   SATURATION PARAMETERS,

IF (LEOCSA .AND. KTILE  ==  1) THEN
  ZSALIN=0.98_JPRB
ELSE
  ZSALIN=1.0_JPRB
ENDIF

! Fix for single precision and spurios cases in which Tsk tile 7 overshoots to very large values.
! This can occur for non-physical tile fractions (<epsilon), but it should be 
! re-evaluated in the future, looking at the surface energy balance in sp.
! This also avoids unphysical Tsk entering cotworestress later in vsurf when pfrti_7<epsilon
ZTSK(KIDIA:KFDIA)=PTSKM1M(KIDIA:KFDIA)
IF (KTILE==7)THEN
  WHERE (PFRTI(KIDIA:KFDIA,KTILE) <= ZEPSILON .AND.ZTSK(KIDIA:KFDIA)>RTT+100._JPRB*ZEPSILON )
    ZTSK(KIDIA:KFDIA)=RTT
  ENDWHERE
ENDIF
DO JL=KIDIA,KFDIA
  !*PQS(JL)=FOEEW(PTSKM1M(JL))/PAPHMS(JL)
  PQS(JL)=FOEEW(ZTSK(JL))/PAPHMS(JL)
  ZCOR=ZSALIN/(1.0_JPRB-RETV  *PQS(JL))
  PQS(JL)=PQS(JL)*ZCOR
!*  PDQS(JL)=PQS(JL)*ZCOR*FOEDESU(PTSKM1M(JL))
  PDQS(JL)=PQS(JL)*ZCOR*FOEDESU(ZTSK(JL))
ENDDO




!*         2.3   DEFINITION OF THE STOMATAL RESISTANCE AND BARE SOIL RES
!*               DOES WORK FOR TYPE 4, 6 AND 8 WHEN ROUTINE IS CALLED FOR 
!*               TYPE 4

 IF (KTILE == 4 ) THEN  !for tiles 6,7, and 8 variables will be already computed
!IF (KTILE==4 .OR. KTILE==6 .OR. KTILE==7 .OR. KTILE==8 ) THEN
!                Compute first liquid fraction of soil water to 
!                be used later in stress functions
!          CONTRIBUTION TO APPARENT ENERGY, TAKING INTO ACCOUNT
!          FREEZING/MELTING OF SOIL WATER.

  DO JK=1,KLEVS
    DO JL=KIDIA,KFDIA
      IF(PTSAM1M(JL,JK) < RTF1.AND.PTSAM1M(JL,JK) > RTF2) THEN
        ZF=0.5_JPRB*(1.0_JPRB-SIN(RTF4*(PTSAM1M(JL,JK)-RTF3)))
      ELSEIF (PTSAM1M(JL,JK) <= RTF2) THEN
        ZF=1.0_JPRB
      ELSE
        ZF=0.0_JPRB
      ENDIF
      IF (LEVGEN) THEN
        JS=KSOTY(JL)
        ZLIQ(JL,JK)=MAX(RWPWPM(JS),MIN(RWCAPM(JS),PWSAM1M(JL,JK)*(1.0_JPRB-ZF)))
        ZLIQR(JL,JK)=MAX(RWRESTM(JS),MIN(RWCAPM(JS),PWSAM1M(JL,JK)*(1.0_JPRB-ZF)))
      ELSE
        ZLIQ(JL,JK)=MAX(RWPWP,MIN(RWCAP,PWSAM1M(JL,JK)*(1.0_JPRB-ZF)))
        ZLIQR(JL,JK)=MAX(0.05_JPRB,MIN(RWCAP,PWSAM1M(JL,JK)*(1.0_JPRB-ZF)))
      ENDIF
    ENDDO
  ENDDO



  DO JL=KIDIA,KFDIA

    ZRSMINL=RVRSMIN(KTVL(JL))
    ZRSMINH=RVRSMIN(KTVH(JL))

!           leaf area index  : ZLAI
    ZLAIL=PLAIL(JL)
    ZLAIH=PLAIH(JL)

!           bare ground fraction
    ZBARE=PFRTI(JL,8)
    IF (KTILES .GT. 9) THEN
     ZBARE=PFRTI(JL,8)+PFRTI(JL,9)
    ENDIF
!           soil moisture stress function : F2
    ZWROOTL=0.0_JPRB
    ZWROOTH=0.0_JPRB
    DO JK=1,KLEVS
      ZWROOTL=ZWROOTL+ZLIQ(JL,JK)*RVROOTSA(JK,KTVL(JL))
      ZWROOTH=ZWROOTH+ZLIQ(JL,JK)*RVROOTSA(JK,KTVH(JL))
    ENDDO
    IF (LEVGEN) THEN
       JS=KSOTY(JL)
       ZWPWP=RWPWPM(JS)
       ZQWEVAP=RQWEVAPM(JS)
!      bare ground evaporation stress is calculated with the weighted average of
!      residual and wilting point soil moisture (since it is common soil)
       ZWPBARE=(RWPWPM(JS)*(1.0_JPRB-ZBARE)+RWRESTM(JS)*ZBARE)
       IF (JS >=1 ) THEN
          ZQWEVAPBARE=1._JPRB/(RWCAPM(JS)-ZWPBARE)
       ELSE
          ZQWEVAPBARE=0._JPRB
       ENDIF
       ZF2B=MAX(RCEPSW,MIN(1.0_JPRB,(ZLIQR(JL,1)-ZWPBARE)*ZQWEVAPBARE))
    ELSE
       ZWPWP=RWPWP
       ZQWEVAP=RQWEVAP
       ZWPBARE=(RWPWP*(1.0_JPRB-ZBARE)+0.05_JPRB*ZBARE)
       ZQWEVAPBARE=1._JPRB/(RWCAP-ZWPBARE)
       ZF2B=MAX(RCEPSW,MIN(1.0_JPRB,(ZLIQ(JL,1)-ZWPWP)*ZQWEVAP))
    ENDIF

    IF (PTSAM1M(JL,1) <= RTT ) THEN
!   if first soil layer temperature is freezing then shutdown transpiration
       ZF2L=RCEPSW
       ZF2H=RCEPSW
    ELSE
       ZF2L=MAX(RCEPSW,MIN(1.0_JPRB,(ZWROOTL-ZWPWP)*ZQWEVAP))
       ZF2H=MAX(RCEPSW,MIN(1.0_JPRB,(ZWROOTH-ZWPWP)*ZQWEVAP))
    ENDIF


!    ZF2L=MAX(RCEPSW,MIN(1.0_JPRB,(ZWROOTL-ZWPWP)*ZQWEVAP))
!    ZF2H=MAX(RCEPSW,MIN(1.0_JPRB,(ZWROOTH-ZWPWP)*ZQWEVAP))

!           radiation stress function (proposed by Alan Betts): ZF1 
    ZSRFL=PSRFD(JL)/250._JPRB
    ZF1L=1.0_JPRB/MAX(1.0_JPRB,0.81_JPRB*(1.+ZSRFL)/(ZSRFL+0.05_JPRB))
    ZF1H=ZF1L

!           atmospheric moisture deficit stress function : F3
    ZHSTRL=RVHSTR(KTVL(JL))
    ZHSTRH=RVHSTR(KTVH(JL))
    ZQSAIR=FOEEW(PTMLEV(JL))/PAPHMS(JL)
    ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAIR)
    ZQSAIR=ZQSAIR*ZCOR
    ZF3L=EXP(-ZHSTRL*(ZQSAIR-PQMLEV(JL)))
    ZF3H=EXP(-ZHSTRH*(ZQSAIR-PQMLEV(JL)))
    ZF3L=MAX(ZEPSF3,MIN(1.0_JPRB,ZF3L))
    ZF3H=MAX(ZEPSF3,MIN(1.0_JPRB,ZF3H))

    IF(ZLAIL /= 0.0_JPRB) THEN
      PWETL(JL)=ZRSMINL/(ZLAIL*ZF1L*ZF2L*ZF3L)
    ELSE
      PWETL(JL) =1.0E+6_JPRB
    ENDIF


    PWETLU(JL)=ZRCLU/(ZF1L*ZF3L) !Pot. Evap. Canopy resist. calc.

    IF(ZLAIH /= 0.0_JPRB) THEN
      PWETH(JL)=ZRSMINH/(ZLAIH*ZF1H*ZF2H*ZF3H)
    ELSE
      PWETH(JL)=1.0E+6_JPRB
    ENDIF   

    PWETHS(JL)=PWETH(JL)
    PWETB(JL)=ZRSMINB/ZF2B

    PEXDIAG(JL,1)=PWETL(JL)
    PEXDIAG(JL,2)=ZF1L
    PEXDIAG(JL,3)=ZF2L
    PEXDIAG(JL,4)=ZF3L
    PEXDIAG(JL,5)=ZLAIL
    PEXDIAG(JL,6)=PWETH(JL)
    PEXDIAG(JL,7)=ZF1H
    PEXDIAG(JL,8)=ZF2H
    PEXDIAG(JL,9)=ZF3H
    PEXDIAG(JL,10)=ZLAIH
  ENDDO

ENDIF !tiles 4 

! Provide also for dew conditions because they are needed in depvel_wes.F90
! introduce zero ing for PQSAM 
!IF (KTILE == 4) THEN
!  DO JL=KIDIA,KFDIA
!    IF (PQMLEV(JL) > PQS(JL)) THEN
!      PWETL(JL)=0.0_JPRB
!    ENDIF
!  ENDDO
!ELSEIF (KTILE == 6) THEN
!  DO JL=KIDIA,KFDIA
!    IF (PQMLEV(JL) > PQS(JL)) THEN
!      PWETH(JL)=0.0_JPRB
!    ENDIF
!  ENDDO
!ELSEIF (KTILE == 7) THEN
!  DO JL=KIDIA,KFDIA
!    IF (PQMLEV(JL) > PQS(JL)) THEN
!      PWETHS(JL)=0.0_JPRB
!    ENDIF
!  ENDDO
!ELSEIF (KTILE == 8) THEN
!  DO JL=KIDIA,KFDIA
!    IF (PQMLEV(JL) > PQS(JL)) THEN
!      PWETB(JL)=0.0_JPRB
!    ENDIF
!  ENDDO
!ENDIF


!*         2.4   APPARENT SURFACE HUMIDITY
!!
IF (KTILE  ==  1.OR. KTILE  ==  2.OR. KTILE  ==  3.OR. KTILE  ==  5 .OR. KTILE == 9 ) THEN   
  DO JL=KIDIA,KFDIA
    PQSAM(JL)=PQS(JL)
  ENDDO
ELSEIF (KTILE  ==  8) THEN
  DO JL=KIDIA,KFDIA
    PQSAM(JL)=PQS(JL)+MIN(0.0_JPRB,(PQMLEV(JL)-PQS(JL)))*PWETB(JL)/(PWETB(JL)+PRAQ(JL))
  ENDDO
ELSEIF (KTILE  ==  4) THEN
  DO JL=KIDIA,KFDIA
    PQSAM(JL)=PQS(JL)+MIN(0.0_JPRB,(PQMLEV(JL)-PQS(JL)))*PWETL(JL)/(PWETL(JL)+PRAQ(JL))
  ENDDO
ELSEIF (KTILE == 6) THEN ! I.E. HIGH VEGETATION, SNOW FREE
  DO JL=KIDIA,KFDIA
    PQSAM(JL)=PQS(JL)+MIN(0.0_JPRB,(PQMLEV(JL)-PQS(JL)))*PWETH(JL)/(PWETH(JL)+PRAQ(JL))
  ENDDO
ELSEIF (KTILE == 10) THEN
  DO JL=KIDIA,KFDIA
    PQSAM(JL)=PQS(JL)+MIN(0.0_JPRB,(PQMLEV(JL)-PQS(JL)))*RURBRES/(RURBRES+PRAQ(JL))
  ENDDO
ELSE ! I.E. HIGH VEGETATION WITH SNOW (7)
  DO JL=KIDIA,KFDIA
    PQSAM(JL)=PQS(JL)+MIN(0.0_JPRB,(PQMLEV(JL)-PQS(JL)))*PWETHS(JL)/(PWETHS(JL)+PRAQ(JL))
  ENDDO
ENDIF


IF (LECTESSEL) THEN  !test on usage of CTESSEL
!=========================================
IF (KTILE==4 .OR. KTILE==6 .OR. KTILE==7 .OR. KTILE==8 .OR. KTILE==10) THEN
!                Compute first liquid fraction of soil water to 
!                be used later in stress functions
!          CONTRIBUTION TO APPARENT ENERGY, TAKING INTO ACCOUNT
!          FREEZING/MELTING OF SOIL WATER.

  DO JK=1,KLEVS
    DO JL=KIDIA,KFDIA
      IF(PTSAM1M(JL,JK) < RTF1.AND.PTSAM1M(JL,JK) > RTF2) THEN
        ZF=0.5_JPRB*(1.0_JPRB-SIN(RTF4*(PTSAM1M(JL,JK)-RTF3)))
      ELSEIF (PTSAM1M(JL,JK) <= RTF2) THEN
        ZF=1.0_JPRB
      ELSE
        ZF=0.0_JPRB
      ENDIF
      IF (LEVGEN) THEN
        JS=KSOTY(JL)
        ZLIQ(JL,JK)=MAX(RWPWPM(JS),MIN(RWCAPM(JS),PWSAM1M(JL,JK)*(1.0_JPRB-ZF)))
        ZLIQR(JL,JK)=MAX(RWRESTM(JS),MIN(RWCAPM(JS),PWSAM1M(JL,JK)*(1.0_JPRB-ZF)))
      ELSE
        ZLIQ(JL,JK)=MAX(RWPWP,MIN(RWCAP,PWSAM1M(JL,JK)*(1.0_JPRB-ZF)))
        ZLIQR(JL,JK)=MAX(0.05_JPRB,MIN(RWCAP,PWSAM1M(JL,JK)*(1.0_JPRB-ZF)))
      ENDIF
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
!           bare ground fraction
    ZBARE=PFRTI(JL,8)
    IF (KTILES .GT. 9) THEN
     ZBARE=PFRTI(JL,8)+PFRTI(JL,9)
    ENDIF
!           soil moisture stress function : F2
    ZWROOTL=0.0_JPRB
    ZWROOTH=0.0_JPRB
    DO JK=1,KLEVS
      ZWROOTL=ZWROOTL+ZLIQ(JL,JK)*RVROOTSA(JK,KTVL(JL))
      ZWROOTH=ZWROOTH+ZLIQ(JL,JK)*RVROOTSA(JK,KTVH(JL))
    ENDDO
    IF (LEVGEN) THEN
       JS=KSOTY(JL)
       ZWPWP=RWPWPM(JS)
       ZQWEVAP=RQWEVAPM(JS)
!      bare ground evaporation stress is calculated with the weighted average of
!      residual and wilting point soil moisture (since it is common soil)
       ZWPBARE=(RWPWPM(JS)*(1.0_JPRB-ZBARE)+RWRESTM(JS)*ZBARE)
       IF (JS >=1 ) THEN
          ZQWEVAPBARE=1._JPRB/(RWCAPM(JS)-ZWPBARE)
       ELSE
          ZQWEVAPBARE=0._JPRB
       ENDIF
       ZF2B=MAX(RCEPSW,MIN(1.0_JPRB,(ZLIQR(JL,1)-ZWPBARE)*ZQWEVAPBARE))
    ELSE
       ZWPWP=RWPWP
       ZQWEVAP=RQWEVAP
       ZWPBARE=(RWPWP*(1.0_JPRB-ZBARE)+0.05_JPRB*ZBARE)
       ZQWEVAPBARE=1._JPRB/(RWCAP-ZWPBARE)
       ZF2B=MAX(RCEPSW,MIN(1.0_JPRB,(ZLIQ(JL,1)-ZWPWP)*ZQWEVAP))
    ENDIF

    IF (PTSAM1M(JL,1) <= RTT ) THEN
!   if first soil layer temperature is freezing then shutdown transpiration
       ZF2L=RCEPSW
       ZF2H=RCEPSW
    ELSE
       ZF2L=MAX(RCEPSW,MIN(1.0_JPRB,(ZWROOTL-ZWPWP)*ZQWEVAP))
       ZF2H=MAX(RCEPSW,MIN(1.0_JPRB,(ZWROOTH-ZWPWP)*ZQWEVAP))
    ENDIF

!    ZF2L=MAX(RCEPSW,MIN(1.0_JPRB,(ZWROOTL-ZWPWP)*ZQWEVAP))
!    ZF2H=MAX(RCEPSW,MIN(1.0_JPRB,(ZWROOTH-ZWPWP)*ZQWEVAP))

!           radiation stress function (proposed by Alan Betts): ZF1 
    ZSRFL=PSRFD(JL)/250._JPRB
    ZF1L=1.0_JPRB/MAX(1.0_JPRB,0.81_JPRB*(1.+ZSRFL)/(ZSRFL+0.05_JPRB))
    ZF1H=ZF1L
!           atmospheric moisture deficit stress function : F3
    ZHSTRL=RVHSTR(KTVL(JL))
    ZHSTRH=RVHSTR(KTVH(JL))
    ZQSAIR=FOEEW(PTMLEV(JL))/PAPHMS(JL)
    ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSAIR)
    ZQSAIR=ZQSAIR*ZCOR
    ZF3L=EXP(-ZHSTRL*(ZQSAIR-PQMLEV(JL)))
    ZF3H=EXP(-ZHSTRH*(ZQSAIR-PQMLEV(JL)))
    ZF3L=MAX(ZEPSF3,MIN(1.0_JPRB,ZF3L))
    ZF3H=MAX(ZEPSF3,MIN(1.0_JPRB,ZF3H))
  ENDDO
ENDIF !tiles 4 

!==========================================

IF (KTILE  ==  1.OR. KTILE  ==  2.OR. KTILE  ==  3.OR. KTILE  ==  5 .OR. KTILE == 9 ) THEN 
  DO JL=KIDIA,KFDIA  
    DO JK=1,KLEVS
      PPWLIQ(JL,JK)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF 

IF (KTILE==4 .OR. KTILE==6 .OR. KTILE==7 .OR. KTILE==8 .OR. KTILE==10) THEN

! default for THE STOMATAL RESISTANCE 
  DO JL=KIDIA,KFDIA
    ZWET(JL)=0.0_JPRB
    ZF2(JL)=1._JPRB
  ENDDO
! The CO2/canopy resistance routine is only called for vegetation tiles
! and the shaded snow tile.

  IF (KTILE==4 .OR. KTILE==6 .OR. KTILE==7 ) THEN
!!                Compute weighted average of unfrozen soil water and F2
    DO JL=KIDIA,KFDIA  
      ZWROOT(JL)=0._JPRB
      DO JK=1,KLEVS
        ZWROOT(JL)=ZWROOT(JL)+ZLIQ(JL,JK)*RVROOTSA(JK,KVEG(JL,KTILE))
      ENDDO
      IF (LEVGEN) THEN
         JS=KSOTY(JL)
         ZWPWP=RWPWPM(JS)
         ZQWEVAP=RQWEVAPM(JS)
      ELSE
         ZWPWP=RWPWP
         ZQWEVAP=RQWEVAP
      ENDIF

      IF (PTSAM1M(JL,1) <= RTT ) THEN
!   if first soil layer temperature is freezing then shutdown transpiration
         ZF2(JL)=RCEPSW
      ELSE
         ZF2(JL)=MAX(RCEPSW,MIN(1.0_JPRB,(ZWROOT(JL)-ZWPWP)*ZQWEVAP))
      ENDIF

!      ZF2(JL)=MAX(RCEPSW,MIN(1.0_JPRB,(ZWROOT(JL)-ZWPWP)*ZQWEVAP))

      IF (KTILE==4 .OR. KTILE==6 ) THEN
         PDHVEGS(JL,JVTTL,5)=ZF2(JL)
      ENDIF
    ENDDO

    DO JL=KIDIA,KFDIA  
      DO JK=1,KLEVS
        PPWLIQ(JL,JK)=ZLIQ(JL,JK)
      ENDDO
    ENDDO

!   Set the soil temperature to be used for acclimation of photosynthetic traits
!   Soil layer 3  (28 - 100cm) as it matches best with past 30-day mean of air T
    DO JL=KIDIA,KFDIA  
       ZTSOIL(JL) = PTSAM1M(JL,3)
    ENDDO

!*         2.3b   STOMATAL RESISTANCE

! COTWORESTRESS (stress parameterisation from JC Calvet) is called. 
! ZWET is not computed for shaded snow, it will be taken from the high veg value
    CALL COTWORESTRESS(KIDIA,KFDIA,KLON,KVEG(:,KTILE),KCO2TYP,PFRTI,&
         & PTMLEV,PQMLEV,PCMLEV,PAPHMS,&
         & ZTSK, ZTSOIL,&
         & PEVAP,PLAI,&
         & PSRFD,PRAQ,PMU0,&
         & ZF2,PQS,&
         & YDCST,YDAGS,YDAGF,YDVEG,YDFLAKE, &
         & PAN,PAG,PRD,&
         & ZWET,ZDSP,ZDMAXT)

! check on dew-fall conditions
    DO JL=KIDIA,KFDIA  
      IF (PQMLEV(JL) > PQS(JL)) ZWET(JL)=0.0_JPRB
      PDHVEGS(JL,JVTTL,6)=ZDSP(JL)
      PDHVEGS(JL,JVTTL,7)=ZDMAXT(JL)
    ENDDO
  ENDIF !vegetation tiles, shaded snow tile

! bare soil (& urban)
  IF( KTILE == 8 .OR. KTILE == 10 ) THEN
    DO JL=KIDIA,KFDIA 
!           bare ground fraction
      ZBARE=PFRTI(JL,8)
      IF ( KTILES .GT. 9 ) THEN
       ZBARE=PFRTI(JL,8)+PFRTI(JL,9)
      ENDIF
      IF (LEVGEN) THEN
         JS=KSOTY(JL)
         ZWPWP=RWPWPM(JS)
         ZQWEVAP=RQWEVAPM(JS)
!      bare ground evaporation stress is calculated with the weighted average of
!      residual and wilting point soil moisture (since it is common soil)
         ZWPBARE=(RWPWPM(JS)*(1.0_JPRB-ZBARE)+RWRESTM(JS)*ZBARE)
         IF (JS >=1 ) THEN
            ZQWEVAPBARE=1._JPRB/(RWCAPM(JS)-ZWPBARE)
         ELSE
            ZQWEVAPBARE=0._JPRB
         ENDIF
         ZF2B=MAX(RCEPSW,MIN(1.0_JPRB,(ZLIQR(JL,1)-ZWPBARE)*ZQWEVAPBARE))
      ELSE
         ZWPWP=RWPWP
         ZQWEVAP=RQWEVAP
         ZWPBARE=(RWPWP*(1.0_JPRB-ZBARE)+0.05_JPRB*ZBARE)
         ZQWEVAPBARE=1._JPRB/(RWCAP-ZWPBARE)
         ZF2B=MAX(RCEPSW,MIN(1.0_JPRB,(ZLIQ(JL,1)-ZWPWP)*ZQWEVAP))
      ENDIF
      ZWET(JL)=ZRSMINB/ZF2B
! check on dew-fall conditions
! we want ZET output 
!      IF (PQMLEV(JL) > PQS(JL)) ZWET(JL)=0.0_JPRB
    ENDDO

  ENDIF
  IF (LEAGS) THEN
! put Rsmin value from ctessel to the actual system variables
    IF (KTILE == 4) THEN
      DO JL=KIDIA,KFDIA
        PWETL(JL)=ZWET(JL)
      ENDDO
    ELSEIF (KTILE == 6) THEN
      DO JL=KIDIA,KFDIA
        PWETH(JL)=ZWET(JL)
      ENDDO
    ELSEIF (KTILE == 7) THEN
      DO JL=KIDIA,KFDIA
        PWETHS(JL)=ZWET(JL)
      ENDDO
    ELSEIF (KTILE == 8) THEN
      DO JL=KIDIA,KFDIA
        PWETB(JL)=ZWET(JL)
      ENDDO
    ELSEIF (KTILE == 10) THEN
      DO JL=KIDIA,KFDIA
        PWETB(JL)=ZWET(JL)
      ENDDO
    ENDIF
  ENDIF

ENDIF !tiles 4, 6, 7 , 8, 10


IF (LEAGS) THEN
!*         2.4   APPARENT SURFACE HUMIDITY
! for bare soil and vegetation tiles
  DO JL=KIDIA,KFDIA
! check only here
    IF (PQMLEV(JL) > PQS(JL)) ZWET(JL)=0.0_JPRB
    IF ( KTILE == 8 .OR. KTILE == 4 .OR. KTILE==6 .OR. KTILE==7 .OR. KTILE==10 ) THEN 
      PQSAM(JL)=PQS(JL)+(PQMLEV(JL)-PQS(JL))*ZWET(JL)/(ZWET(JL)+PRAQ(JL))
    ELSE
      PQSAM(JL)=PQS(JL)
    ENDIF
  ENDDO
ENDIF
ENDIF ! LECTESSEL test

!*         2.5   DRY STATIC ENERGY AT THE SURFACE

DO JL=KIDIA,KFDIA
  PCPTS(JL)=PTSKM1M(JL)*RCPD*(1.0_JPRB+RVTMP2*PQSAM(JL))
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('VSURF_MOD:VSURF',1,ZHOOK_HANDLE)
END SUBROUTINE VSURF
END MODULE VSURF_MOD
