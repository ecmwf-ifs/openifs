MODULE SURFEXCDRIVERS_CTL_MOD
CONTAINS
SUBROUTINE SURFEXCDRIVERS_CTL( &
 &   KIDIA, KFDIA, KLON, KLEVS, KTILES, KSTEP &
 & , PTSTEP, PRVDIFTS &
 & , LDSURF2 &
! input data, non-tiled
 & , KTVL, KTVH, PCVL, PCVH &
 & , PLAIL, PLAIH &
 & , PUMLEV, PVMLEV, PTMLEV, PQMLEV, PAPHMS, PGEOMLEV, PCPTGZLEV &
 & , PSST, PTSKM1M, PCHAR, PSSRFL, PTICE, PTSNOW &
 & , PWLMX &
! input data, soil
 & , PTSAM1M, PWSAM1M, KSOTY &
! input data, tiled
 & , PFRTI, PALBTI &
!
 & , YDCST, YDEXC, YDVEG, YDSOIL, YDFLAKE & 
! updated data, tiled
 & , PUSTRTI, PVSTRTI, PAHFSTI, PEVAPTI, PTSKTI &
! updated data, non-tiled
 & , PZ0M, PZ0H &
! output data, tiled
 & , PSSRFLTI, PQSTI, PDQSTI, PCPTSTI, PCFHTI, PCFQTI, PCSATTI, PCAIRTI &
! output data, non-tiled
 & , PCFMLEV, PKMFL, PKHFL, PKQFL, PEVAPSNW &
 & , PZ0MW, PZ0HW, PZ0QW, PCPTSPP, PQSAPP, PBUOMPP &
 & )

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_THF  , ONLY : R4LES, R5LES, R2ES, R4IES, R3LES, R3IES, R5IES
USE YOS_CST  , ONLY : TCST
USE YOS_EXC  , ONLY : TEXC
USE YOS_VEG  , ONLY : TVEG
USE YOS_SOIL , ONLY : TSOIL
USE YOS_FLAKE, ONLY : TFLAKE
USE VUPDZ0S_MOD
USE VSURFS_MOD
USE VEXCSS_MOD
USE VEVAPS_MOD

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!------------------------------------------------------------------------

!  PURPOSE:
!    Routine SURFEXCDRIVERS controls the ensemble of routines that prepare
!    the surface exchange coefficients and associated surface quantities
!    needed for the solution of the vertical diffusion equations. 

!  SURFEXCDRIVERS is called by VDFMAINS

!  METHOD:
!    This routine is only a shell needed by the surface library
!    externalisation.

!  AUTHOR:
!    P. Viterbo       ECMWF May 2005   

!  REVISION HISTORY:
!    M. Janiskova     27/06/2005 removed option for MASS vector functions
!                                not use in correcsponding TL/AD
!    A.Beljaars       10/12/2005 TOFD
!    M. Janiskova     10/03/2006 call for simplified routines (suffix s)
!                                instead of full NL routines
!    G. Balsamo       03/07/2006 Add soil type
!    M. Janiskova     21/05/2007 clean-up of roughness length initialization
!    S. Boussetta/G.Balsamo May 2009 Add lai
!    M. Janiskova     July 2011->2013  modified computation of snow evaporation
!    M. Janiskova     Jan 2015   use previous time step fluxes for heat&momentum

!  INTERFACE: 

!    Integers (In):
!      KIDIA    :    Begin point in arrays
!      KFDIA    :    End point in arrays
!      KLON     :    Length of arrays
!      KLEVS    :    Number of soil layers
!      KTILES   :    Number of tiles
!      KSTEP    :    Time step index
!      KTVL     :    Dominant low vegetation type
!      KTVH     :    Dominant high vegetation type
!      KSOTY    :    SOIL TYPE                                        (1-7)

!    Reals (In):
!      PTSTEP   :    Timestep
!      PRVDIFTS :    Semi-implicit factor for vertical diffusion discretization

!    Reals with tile index (In): 
!      PFRTI    :    TILE FRACTIONS                                   (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!      PALBTI   :    Tile albedo                                      (0-1)

!    Reals independent of tiles (In):
!      PCVL     :    LOW VEGETATION COVER                             -  
!      PCVH     :    HIGH VEGETATION COVER                            -  
!      PLAIL    :    LOW VEGETATION LAI
!      PLAIH    :    HIGH VEGETATION LAI

!  Logical:
!      LDSURF2  :    TRUE when simplified surface scheme called

!      PUMLEV   :    X-VELOCITY COMPONENT, lowest atmospheric level   m/s
!      PVMLEV   :    Y-VELOCITY COMPONENT, lowest atmospheric level   m/s
!      PTMLEV   :    TEMPERATURE,   lowest atmospheric level          K
!      PQMLEV   :    SPECIFIC HUMIDITY                                kg/kg
!      PAPHMS   :    Surface pressure                                 Pa
!      PGEOMLEV :    Geopotential, lowest atmospehric level           m2/s2
!      PCPTGZLEV:    Geopotential, lowest atmospehric level           J/kg
!      PSST     :    (OPEN) SEA SURFACE TEMPERATURE                   K
!      PTSKM1M  :    SKIN TEMPERATURE                                 K
!      PCHAR    :    "EQUIVALENT" CHARNOCK PARAMETER                  -
!      PSSRFL   :    NET SHORTWAVE RADIATION FLUX AT SURFACE          W/m2
!      PTSAM1M  :    SURFACE TEMPERATURE                              K
!      PWSAM1M  :    SOIL MOISTURE ALL LAYERS                         m**3/m**3
!      PTICE    :    Ice temperature, top slab                        K
!      PTSNOW   :    Snow temperature                                 K
!      PWLMX    :    Maximum interception layer capacity              kg/m**2

!    Reals with tile index (In/Out):
!      PUSTRTI  :    SURFACE U-STRESS                                 N/m2 
!      PVSTRTI  :    SURFACE V-STRESS                                 N/m2
!      PAHFSTI  :    SURFACE SENSIBLE HEAT FLUX                       W/m2
!      PEVAPTI  :    SURFACE MOISTURE FLUX                            KG/m2/s
!      PTSKTI   :    SKIN TEMPERATURE                                 K

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

!    Reals independent of tiles (Out):
!      PCFMLEV  :    PROP. TO EXCH. COEFF. FOR MOMENTUM               ????
!                     (C-STAR IN DOC.) (SURFACE LAYER ONLY)
!      PKMFL    :    Kinematic momentum flux                          ????
!      PKHFL    :    Kinematic heat flux                              ????
!      PKQFL    :    Kinematic moisture flux                          ????
!      PEVAPSNW :    Evaporation from snow under forest               kgm-2s-1
!      PZ0MW    :    Roughness length for momentum, WMO station       m
!      PZ0HW    :    Roughness length for heat, WMO station           m
!      PZ0QW    :    Roughness length for moisture, WMO station       m
!      PCPTSPP  :    Cp*Ts for post-processing of weather parameters  J/kg
!      PQSAPP   :    Apparent surface humidity for post-processing    kg/kg
!                     of weather parameters
!      PBUOMPP  :    Buoyancy flux, for post-processing of gustiness  ???? 


!     EXTERNALS.
!     ----------

!     ** SURFEXCDRIVERS_CTL CALLS SUCCESSIVELY:
!         *VUPDZ0S*
!         *VSURFS*
!         *VEXCSS*
!         *VEVAPS*

!  DOCUMENTATION:
!    See Physics Volume of IFS documentation

!------------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRVDIFTS
LOGICAL           ,INTENT(IN)    :: LDSURF2

INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTY(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAIL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAIH(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOMLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHAR(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTICE(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSNOW(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWLMX(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSAM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSAM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBTI(:,:) 
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TEXC)        ,INTENT(IN)    :: YDEXC
TYPE(TVEG)        ,INTENT(IN)    :: YDVEG
TYPE(TSOIL)       ,INTENT(IN)    :: YDSOIL
TYPE(TFLAKE)      ,INTENT(IN)    :: YDFLAKE
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSTRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSTRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKTI(:,:) 
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
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCFMLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKMFL(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKHFL(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKQFL(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEVAPSNW(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ0MW(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ0HW(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ0QW(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCPTSPP(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQSAPP(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBUOMPP(:)

! Local variables

INTEGER(KIND=JPIM) :: IFRMAX(KLON)

REAL(KIND=JPRB) :: ZZ0MTI(KLON,KTILES) , ZZ0HTI(KLON,KTILES) ,&
                 & ZZ0QTI(KLON,KTILES) , ZBUOMTI(KLON,KTILES),&
                 & ZZDLTI(KLON,KTILES) , ZRAQTI(KLON,KTILES) ,&
                 & ZQSATI(KLON,KTILES) , ZCFMTI(KLON,KTILES) ,&
                 & ZKMFLTI(KLON,KTILES), ZKHFLTI(KLON,KTILES),&
                 & ZKQFLTI(KLON,KTILES), ZZQSATI(KLON,KTILES)

REAL(KIND=JPRB) :: ZFRMAX(KLON)   , ZALB(KLON)     , ZSSRFL1(KLON)  , &
                 & ZSRFD(KLON)    , ZWETL(KLON)    , ZWETH(KLON)    , &
                 & ZWETHS(KLON)   , ZWETB(KLON)    , &
                 & ZTSA(KLON)     , ZCSNW(KLON)

INTEGER(KIND=JPIM) :: JL, JTILE, KTILE

REAL(KIND=JPRB) :: ZQSSN, ZCOR, ZCONS1, ZZ0MWMO, ZZ0HWMO, ZCDRO 
REAL(KIND=JPRB) :: ZDIV1, ZDIV2, Z3S, Z4S
REAL(KIND=JPRB) :: ZCONS2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

LOGICAL         :: LLAND, LLSICE, LLHISSR(KLON)

#include "fcsttre.h"

!*         1.     Set up of general quantities
!                 ----------------------------

IF (LHOOK) CALL DR_HOOK('SURFEXCDRIVERS_CTL_MOD:SURFEXCDRIVERS_CTL',0,ZHOOK_HANDLE)
ASSOCIATE(RCPD=>YDCST%RCPD, RD=>YDCST%RD, RETV=>YDCST%RETV, RG=>YDCST%RG, &
 & RSIGMA=>YDCST%RSIGMA, RTT=>YDCST%RTT, &
 & REPDU2=>YDEXC%REPDU2, RKAP=>YDEXC%RKAP, RZ0ICE=>YDEXC%RZ0ICE, &
 & RALFMAXSN=>YDSOIL%RALFMAXSN, &
 & RVTRSR=>YDVEG%RVTRSR, RVZ0M=>YDVEG%RVZ0M)

ZCONS1=1./(RG*PTSTEP)

!*         1.1  ESTIMATE SURF.FL. FOR STEP 0
!*              (ASSUME NEUTRAL STRATIFICATION)

IF ( KSTEP == 0) THEN
  DO JTILE=2,KTILES
    DO JL=KIDIA,KFDIA
      PTSKTI(JL,JTILE)=PTSKM1M(JL)
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    PTSKTI(JL,1)=PSST(JL)
  ENDDO
ENDIF

!*         1.2  UPDATE Z0

CALL VUPDZ0S(KIDIA,KFDIA,KLON,KTILES,KSTEP,&
     & KTVL,KTVH,PCVL,PCVH,PUMLEV, PVMLEV,&
     & PTMLEV,PQMLEV,PAPHMS,PGEOMLEV,&
     & PUSTRTI,PVSTRTI,PAHFSTI,PEVAPTI,&
     & PTSKTI,PCHAR,PFRTI, &
     & YDCST,YDEXC,YDVEG,YDFLAKE, &
     & ZZ0MTI,ZZ0HTI,ZZ0QTI,ZBUOMTI,ZZDLTI,ZRAQTI)

!*         1.3  FIND DOMINANT SURFACE TYPE parameters for postprocessing

ZFRMAX(KIDIA:KFDIA)=PFRTI(KIDIA:KFDIA,1)
IFRMAX(KIDIA:KFDIA)=1
DO JTILE=2,KTILES
  DO JL=KIDIA,KFDIA
    IF (PFRTI(JL,JTILE)  >  ZFRMAX(JL)) THEN
      ZFRMAX(JL)=PFRTI(JL,JTILE)
      IFRMAX(JL)=JTILE
    ENDIF
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  JTILE=IFRMAX(JL)
  PZ0M(JL)=ZZ0MTI(JL,JTILE)
  PZ0H(JL)=ZZ0HTI(JL,JTILE)
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
    IF (LLHISSR(JL)) THEN
      PSSRFLTI(JL,JTILE)=PSSRFLTI(JL,JTILE)*PSSRFL(JL)/ZSSRFL1(JL)
    ENDIF
    ZSRFD(JL)=PSSRFLTI(JL,JTILE)/(1.0_JPRB-PALBTI(JL,JTILE))  
  ENDDO

  CALL VSURFS(KIDIA,KFDIA,KLON,KLEVS,JTILE,&
   & KTVL, KTVH,&
   & PLAIL, PLAIH,&
   & PTMLEV,PQMLEV,PAPHMS,&
   & PTSKTI(:,JTILE),PWSAM1M,PTSAM1M,KSOTY,&
   & ZSRFD,ZRAQTI(:,JTILE),&
   & YDCST,YDVEG,YDSOIL,&
   & ZQSATI(:,JTILE),PQSTI(:,JTILE),PDQSTI(:,JTILE),&
   & ZWETB,PCPTSTI(:,JTILE),ZWETL,ZWETH,ZWETHS)
ENDDO

!*         3.     EXCHANGE COEFFICIENTS
!                 ---------------------

!*         3.1  SURFACE EXCHANGE COEFFICIENTS
DO JTILE=1,KTILES

  CALL VEXCSS(KIDIA,KFDIA,KLON,PTSTEP,PRVDIFTS,&
   & PUMLEV,PVMLEV,PTMLEV,PQMLEV,PAPHMS,PGEOMLEV,PCPTGZLEV,&
   & PCPTSTI(:,JTILE),ZQSATI(:,JTILE),&
   & ZZ0MTI(:,JTILE),ZZ0HTI(:,JTILE),&
   & ZZ0QTI(:,JTILE),ZBUOMTI(:,JTILE),&
   & YDCST,YDEXC,&
   & ZCFMTI(:,JTILE),PCFHTI(:,JTILE),&
   & PCFQTI(:,JTILE))
ENDDO

!*         3.2  EQUIVALENT EVAPOTRANSPIRATION EFFICIENCY COEFFICIENT

DO JTILE=1,KTILES
  IF     (JTILE == 1) THEN
    ZTSA(KIDIA:KFDIA)=PSST(KIDIA:KFDIA)
  ELSEIF (JTILE == 2) THEN
    ZTSA(KIDIA:KFDIA)=PTICE(KIDIA:KFDIA)
  ELSEIF (JTILE == 5 .OR. JTILE == 7) THEN
    ZTSA(KIDIA:KFDIA)=PTSNOW(KIDIA:KFDIA)
  ELSE
    ZTSA(KIDIA:KFDIA)=PTSAM1M(KIDIA:KFDIA,1)
  ENDIF
  CALL VEVAPS(KIDIA,KFDIA,KLON,PTSTEP,PRVDIFTS,JTILE,&
   & PWLMX ,PTMLEV  ,PQMLEV  ,PAPHMS,PTSKTI(:,JTILE),ZTSA,&
   & PQSTI(:,JTILE),PCFQTI(:,JTILE),ZWETB,ZWETL,ZWETH,ZWETHS,&
   & YDCST,YDVEG,&
   & PCPTSTI(:,JTILE),PCSATTI(:,JTILE),PCAIRTI(:,JTILE),&
   & ZCSNW)  
ENDDO

!          COMPUTE SNOW EVAPORATION FROM BELOW TREES i.e. TILE 7

IF (LDSURF2) THEN

! Note the use of qsat(Tsnow), rather than tile 7 skin. Skin T7 is a
! canopy temperature, definitely not what is desirable. Skin T5 can go
! up (and down ..) freely, not really what we want. The use of
! qsat (Tsnow) is tantamount to neglecting the skin effect there.

  DO JL=KIDIA,KFDIA
    IF (PFRTI(JL,7) > 0.0_JPRB) THEN
!Only compute snow evap. when there is a snow cover under trees (just to skip funct computation)
!      ZQSSN=FOEEW(PTSNOW(JL))/PAPHMS(JL)
      IF (PTSNOW(JL) > RTT) THEN
        ZDIV1 = 1.0_JPRB/(PTSNOW(JL)-R4LES)
        Z3S = R3LES*(PTSNOW(JL)-RTT)*ZDIV1
      ELSE
        ZDIV1 = 1.0_JPRB/(PTSNOW(JL)-R4IES)
        Z3S = R3IES*(PTSNOW(JL)-RTT)*ZDIV1
      ENDIF
      Z4S = EXP(Z3S)
      ZDIV2 = 1.0_JPRB/PAPHMS(JL)
      ZQSSN = (R2ES*Z4S)*ZDIV2

      ZCOR = 1.0_JPRB/(1.0_JPRB-RETV  *ZQSSN)
      ZQSSN= ZQSSN*ZCOR
      PEVAPSNW(JL)=ZCONS1*PCFQTI(JL,7)*ZCSNW(JL)*&
       & (PQMLEV(JL)-ZQSSN)
    ELSE
      PEVAPSNW(JL)=0.0_JPRB
    ENDIF
  ENDDO
ELSE
  DO JL=KIDIA,KFDIA
    ZQSSN=FOEEW(PTSNOW(JL))/PAPHMS(JL)
    ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSSN)
    ZQSSN=ZQSSN*ZCOR
    PEVAPSNW(JL)=ZCONS1*PCFQTI(JL,7)*ZCSNW(JL)*&
     & (PQMLEV(JL)-ZQSSN)
  ENDDO
ENDIF

!*         3.3  COMPUTE SURFACE FLUXES FOR TILES
!               (replaces vsflx routine without currents)

ZCONS2 = 1.0_JPRB/(RG*PTSTEP*PRVDIFTS)
DO JL=KIDIA,KFDIA
  DO JTILE=1,KTILES
    ZCDRO = ( RD*PTMLEV(JL)*(1.0_JPRB+RETV*PQMLEV(JL)) )/PAPHMS(JL)
    ZKMFLTI(JL,JTILE) = ZCDRO * ZCONS2 * ZCFMTI(JL,JTILE) &
        & * SQRT(PUMLEV(JL)**2+PVMLEV(JL)**2)

!   use previous times step fluxes for heat and moisture
    ZKHFLTI(JL,JTILE) = ZCDRO * PAHFSTI(JL,JTILE) / RCPD
    ZKQFLTI(JL,JTILE) = ZCDRO * PEVAPTI(JL,JTILE) 
  ENDDO
ENDDO

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

!          POST-PROCESSING WITH LOCAL INSTEAD OF EFFECTIVE
!          SURFACE ROUGHNESS LENGTH. THE LOCAL ONES ARE FOR
!          WMO-TYPE WIND STATIONS I.E. OPEN TERRAIN WITH GRASS

ZZ0MWMO=0.03_JPRB
ZZ0HWMO=0.003_JPRB
DO JL=KIDIA,KFDIA
  JTILE=IFRMAX(JL)
  IF (JTILE  >  2.AND. ZZ0MTI(JL,JTILE)  >  ZZ0MWMO) THEN
    PZ0MW(JL)=ZZ0MWMO
    PZ0HW(JL)=ZZ0HWMO
    PZ0QW(JL)=ZZ0HWMO
  ELSE
    PZ0MW(JL)=ZZ0MTI(JL,JTILE)
    PZ0HW(JL)=ZZ0HTI(JL,JTILE)
    PZ0QW(JL)=ZZ0QTI(JL,JTILE)
  ENDIF
ENDDO

DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    ZZQSATI(JL,JTILE)=PQMLEV(JL)*(1.0_JPRB-PCAIRTI(JL,JTILE))&
     & +PCSATTI(JL,JTILE)*PQSTI(JL,JTILE)  
    ZZQSATI(JL,JTILE)=MAX(1.0E-12_JPRB,ZZQSATI(JL,JTILE))
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  JTILE=IFRMAX(JL)
  PCPTSPP(JL)=PCPTSTI(JL,JTILE)
  PQSAPP(JL)=ZZQSATI(JL,JTILE)
  PBUOMPP(JL)=ZBUOMTI(JL,JTILE)
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFEXCDRIVERS_CTL_MOD:SURFEXCDRIVERS_CTL',1,ZHOOK_HANDLE)

END SUBROUTINE SURFEXCDRIVERS_CTL
END MODULE SURFEXCDRIVERS_CTL_MOD
