MODULE SURFPP_CTL_MOD
CONTAINS
SUBROUTINE SURFPP_CTL( KIDIA,KFDIA,KLON,KTILES, KDHVTLS, KDHFTLS &
 & , PTSTEP &
! input
 & , PFRTI, PAHFLTI, PG0TI, PSTRTULEV, PSTRTVLEV, PTSKM1M &
 & , PUMLEV, PVMLEV, PQMLEV, PGEOMLEV, PCPTSPP ,PCPTGZLEV &
 & , PAPHMS, PZ0MW, PZ0HW, PZ0QW, PZDL, PQSAPP, PBLEND, PFBLEND, PBUOM &
 & , PZ0M, PEVAPSNW, PSSRFLTI, PSLRFL, PSST &
 & , PUCURR, PVCURR, PUSTOKES, PVSTOKES &
 & , YDCST, YDEXC, YDFLAKE &
! updated
 & , PAHFSTI, PEVAPTI, PTSKE1, PTSKTIP1 &
! output
 & , PDIFTSLEV, PDIFTQLEV, PUSTRTI, PVSTRTI, PTSKTI, PAHFLEV, PAHFLSB, PFWSB  &
 & , PU10M, PV10M, PT2M, PD2M, PQ2M &
 & , PGUST, P10NU, P10NV, PUST &
! output DDH
 & , PDHTLS &
 & , PRPLRG)

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST  , ONLY : TCST
USE YOS_EXC  , ONLY : TEXC
USE YOS_FLAKE, ONLY : TFLAKE

USE SPPCFL_MOD
USE SPPGUST_MOD
USE VOSKIN_MOD

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!------------------------------------------------------------------------

!  PURPOSE:
!    Routine SURFPP controls the computation of quantities at the end of
!     vertical diffusion, including routines to post-process weather elements
!     and gustiness.

!  SURFPP is called by VDFMAIN

!  METHOD:
!    This routine is a shell needed by the surface library  externalisation.

!  AUTHOR:
!    P. Viterbo       ECMWF May 2005

!  REVISION HISTORY:
!    05-01-2006    T. Stockdale   ocean surface currents
!    A. Beljaars      ECMWF Feb 2006  Revised gust to accomodate stochastic physics
!    E. Dutra/G. Balsamo    May 2008  Add lake tile
!    N.Semane+P.Bechtold 04-10-2012 Add PRPLRG factor for small planet

!  INTERFACE: 

!    Integers (In):
!      KIDIA    :    Begin point in arrays
!      KFDIA    :    End point in arrays
!      KLON     :    Length of arrays
!      KTILES   :    Number of files
!      KDHVTLS  :    Number of variables for individual tiles
!      KDHFTLS  :    Number of fluxes for individual tiles


!    Reals (In):
!      PTSTEP    :  Timestep                                          s
!      PFRTI     :  TILE FRACTIONS                                   (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!            9 : LAKE                  10 : URBAN
!      PAHFLTI   :  Surface latent heat flux                         Wm-2
!      PG0TI     :  Surface ground heat flux                         W/m2
!      PSTRTULEV :  TURBULENT FLUX OF U-MOMEMTUM                     kg/(m*s2)
!      PSTRTVLEV :  TURBULENT FLUX OF V-MOMEMTUM                     kg/(m*s2)
!      PTSKM1M   :  Skin temperature, t                              K
!      PUMLEV    :  X-VELOCITY COMPONENT, lowest atmospheric level   m/s
!      PVMLEV    :  Y-VELOCITY COMPONENT, lowest atmospheric level   m/s
!      PQMLEV    :  SPECIFIC HUMIDITY                                kg/kg
!      PGEOMLEV  :  Geopotential, lowest atmospehric level           m2/s2
!      PCPTSPP   :  Cp*Ts for post-processing of weather parameters  J/kg
!      PCPTGZLEV :  Geopotential, lowest atmospehric level           J/kg
!      PAPHMS    :  Surface pressure                                 Pa
!      PZ0MW     :  Roughness length for momentum, WMO station       m
!      PZ0HW     :  Roughness length for heat, WMO station           m
!      PZ0QW     :  Roughness length for moisture, WMO station       m
!      PZDL      :  z/L                                              -
!      PQSAPP    :  Apparent surface humidity                        kg/kg
!      PBLEND    :  Blending weight for 10 m wind postprocessing     m
!      PFBLEND   :  Wind speed at blending weight for 10 m wind PP   m/s
!      PBUOM     :  Buoyancy flux, for post-processing of gustiness  ????
!      PZ0M     :    AERODYNAMIC ROUGHNESS LENGTH                    m
!      PEVAPSNW :    Evaporation from snow under forest              kgm-2s-1
!      PSSRFLTI  :  NET SOLAR RADIATION AT THE SURFACE, TILED        Wm-2
!      PSLRFL    :  NET THERMAL RADIATION AT THE SURFACE             Wm-2
!      PSST      :  Sea surface temperatute                          K
!      PUCURR    :  U component of ocean surface current             m/s
!      PVCURR    :  V component of ocean surface current             m/s
!      PUSTOKES  :  U component of surface Stokes velocity           m/s
!      PVSTOKES  :  V component of surface Stokes velocity           m/s

!    Reals (Updated):
!      PAHFSTI   :  SURFACE SENSIBLE HEAT FLUX                       W/m2
!      PEVAPTI   :  SURFACE MOISTURE FLUX                            kg/m2/s
!      PTSKE1    :  SKIN TEMPERATURE TENDENCY                        K/s
!      PTSKTIP1  :  Tile skin temperature, t+1                       K

!    Reals (Out):
!      PDIFTSLEV :  TURBULENT FLUX OF HEAT                           J/(m2*s)
!      PDIFTQLEV :  TURBULENT FLUX OF SPECIFIC HUMIDITY              kg/(m2*s)
!      PUSTRTI   :  SURFACE U-STRESS                                 N/m2 
!      PVSTRTI   :  SURFACE V-STRESS                                 N/m2 
!      PTSKTI    :  SKIN TEMPERATURE                                 K
!      PAHFLEV   :  LATENT HEAT FLUX  (SNOW/ICE FREE PART)           W/m2
!      PAHFLSB   :  LATENT HEAT FLUX  (SNOW/ICE COVERED PART)        W/m2
!      PFWSB     :  EVAPORATION OF SNOW                              kg/(m**2*s)
!      PU10M     :  U-COMPONENT WIND AT 10 M                         m/s
!      PV10M     :  V-COMPONENT WIND AT 10 M                         m/s
!      P10NU     :  U-COMPONENT NEUTRAL WIND AT 10 M                 m/s
!      P10NV     :  V-COMPONENT NEUTRAL WIND AT 10 M                 m/s
!      PUST      :  FRICTION VELOCITY                                m/s
!      PT2M      :  TEMPERATURE AT 2M                                K
!      PD2M      :  DEW POINT TEMPERATURE AT 2M                      K
!      PQ2M      :  SPECIFIC HUMIDITY AT 2M                          kg/kg
!      PGUST     :  GUST AT 10 M                                     m/s
!      PDHTLS    :  Diagnostic array for tiles (see module yomcdh)
!                      (Wm-2 for energy fluxes, kg/(m2s) for water fluxes)

!     EXTERNALS.
!     ----------

!     ** SURFPP_CTL CALLS SUCCESSIVELY:
!         *SPPCFL*
!         *SPPGUST*
!         *VOSKIN*

!  DOCUMENTATION:
!    See Physics Volume of IFS documentation

!------------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTLS
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTLS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAHFLTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PG0TI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTRTULEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTRTVLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOMLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTSPP(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0MW(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0HW(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0QW(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZDL(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSAPP(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBLEND(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFBLEND(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBUOM(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0M(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVAPSNW(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFLTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSTOKES(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSTOKES(:)
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TEXC)        ,INTENT(IN)    :: YDEXC
TYPE(TFLAKE)      ,INTENT(IN)    :: YDFLAKE
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKE1(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKTIP1(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTSLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTQLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUSTRTI(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVSTRTI(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PTSKTI(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAHFLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAHFLSB(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFWSB(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU10M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PV10M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P10NU(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P10NV(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUST(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT2M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD2M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQ2M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGUST(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTLS(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPLRG

! Local variables

INTEGER(KIND=JPIM) :: JTILE, JL

REAL(KIND=JPRB) :: ZTSK(KLON),ZAHFSM(KLON),ZEVAPM(KLON),ZUSTAR(KLON)
REAL(KIND=JPRB) :: ZRTMST,ZRHO
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFPP_CTL_MOD:SURFPP_CTL',0,ZHOOK_HANDLE)
ASSOCIATE(LEFLAKE=>YDFLAKE%LEFLAKE, &
 & RD=>YDCST%RD, RETV=>YDCST%RETV, RLSTT=>YDCST%RLSTT, &
 & LEOCWA=>YDEXC%LEOCWA, LEOCCO=>YDEXC%LEOCCO, REPUST=>YDEXC%REPUST)

ZRTMST      = 1.0_JPRB/PTSTEP    ! optimization

!*         1.     SURFACE FLUXES - TILES
!                 ----------------------

!*         1.1  SURFACE FLUXES OF HEAT AND MOISTURE FOR THE 
!*              DIFFERENT TILES AND THE MEAN OVER TILES

ZAHFSM(KIDIA:KFDIA) = 0.0_JPRB
ZEVAPM(KIDIA:KFDIA) = 0.0_JPRB
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    ZAHFSM(JL)=ZAHFSM(JL)+PFRTI(JL,JTILE)*PAHFSTI(JL,JTILE)
    ZEVAPM(JL)=ZEVAPM(JL)+PFRTI(JL,JTILE)*PEVAPTI(JL,JTILE)
  ENDDO
ENDDO

PDIFTSLEV  (KIDIA:KFDIA) = 0.0_JPRB
PDIFTQLEV  (KIDIA:KFDIA) = 0.0_JPRB
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    PDIFTSLEV(JL)=PDIFTSLEV(JL)+PFRTI(JL,JTILE)*PAHFSTI(JL,JTILE)
    PDIFTQLEV(JL)=PDIFTQLEV(JL)+PFRTI(JL,JTILE)*PEVAPTI(JL,JTILE)

    PUSTRTI(JL,JTILE)=PSTRTULEV(JL)
    PVSTRTI(JL,JTILE)=PSTRTVLEV(JL)
    IF (PFRTI(JL,JTILE) == 0._JPRB .AND. JTILE <= 8 ) THEN
      PAHFSTI(JL,JTILE)=ZAHFSM(JL)
      PEVAPTI(JL,JTILE)=ZEVAPM(JL)
    ENDIF
  ENDDO
ENDDO

IF (SIZE(PDHTLS) > 0) CALL COMPUTE_DDH

!*         1.2  PARAMETERS AND DERIVATIVES (SET TO 0) FOR LAND 
!*              SURFACE SCHEME

DO JL=KIDIA,KFDIA
  PAHFLEV(JL)=PFRTI(JL,1)*PAHFLTI(JL,1)&
   & +PFRTI(JL,3)*PAHFLTI(JL,3)&
   & +PFRTI(JL,4)*PAHFLTI(JL,4)&
   & +PFRTI(JL,6)*PAHFLTI(JL,6)&
   & +PFRTI(JL,7)*(PAHFLTI(JL,7)-RLSTT*PEVAPSNW(JL))&
   & +PFRTI(JL,8)*PAHFLTI(JL,8) 
  IF (LEFLAKE) THEN
    PAHFLEV(JL)=PAHFLEV(JL)+PFRTI(JL,9)*PAHFLTI(JL,9)     
  ENDIF
  PAHFLSB(JL)=PFRTI(JL,2)*PAHFLTI(JL,2)+PFRTI(JL,5)*PAHFLTI(JL,5)&
   & +PFRTI(JL,7)*RLSTT*PEVAPSNW(JL)  
  PFWSB(JL) = PFRTI(JL,5)*PEVAPTI(JL,5)&
   & +PFRTI(JL,7)*PEVAPSNW(JL)  
ENDDO

DO JL=KIDIA,KFDIA
  ZRHO = PAPHMS(JL)/( RD*PTSKM1M(JL)*(1.0_JPRB+RETV*PQMLEV(JL)) )
  ZUSTAR(JL)=MAX(REPUST,SQRT(SQRT(PSTRTULEV(JL)**2+PSTRTVLEV(JL)**2)/ZRHO))
ENDDO

!      2. Post-processing of weather parameters
!         -------------------------------------

CALL SPPCFL(KIDIA,KFDIA,KLON &
 & , PUMLEV, PVMLEV, PQMLEV, PGEOMLEV, PCPTSPP, PCPTGZLEV &
 & , PAPHMS, PZ0MW, PZ0HW, PZ0QW, PZDL, PQSAPP &
 & , PBLEND, PFBLEND, PUCURR, PVCURR &
 & , YDCST, YDEXC &
 & , PU10M, PV10M, P10NU, P10NV, PUST, PT2M, PD2M, PQ2M, PRPLRG)

!      3. Post-processing of wind gusts
!         -----------------------------

CALL SPPGUST(KIDIA,KFDIA,KLON &
 & , PZ0M, PBUOM, ZUSTAR, PU10M, PV10M &
 & , YDCST, YDEXC &
 & , PGUST )

!         4. SKIN LAYER 
!            ---- -----

!         4.1 OCEAN SKIN EFFECTS I.E. TILE 1 ONLY
 
IF (LEOCWA .OR. LEOCCO) THEN                 
  CALL VOSKIN(KIDIA,KFDIA,KLON,&
   & PTSTEP,&
   & PSSRFLTI(:,1) ,PSLRFL ,PAHFSTI(:,1),PAHFLTI(:,1),&
   & PUSTRTI(:,1),PVSTRTI(:,1),&
   & PUMLEV,PVMLEV,PTSKTI(:,1),PSST,PUSTOKES,PVSTOKES,&
   & YDCST,YDEXC,&
   & PTSKTIP1(:,1),PRPLRG)  
ENDIF

!         4.2 SKIN TEMPERATURE AVERAGING OVER TILES 

ZTSK(KIDIA:KFDIA)=0.0_JPRB
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    ZTSK(JL)=ZTSK(JL)+PFRTI(JL,JTILE)*PTSKTIP1(JL,JTILE)
  ENDDO
ENDDO

DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    IF (PFRTI(JL,JTILE) == 0._JPRB .AND. JTILE <= 8 ) THEN
      PTSKTI(JL,JTILE)=ZTSK(JL)
    ELSE
      PTSKTI(JL,JTILE)=PTSKTIP1(JL,JTILE)
    ENDIF
  ENDDO
ENDDO


!         4.3 SKIN LAYER TENDENCY 

DO JL=KIDIA,KFDIA
  PTSKE1(JL) = PTSKE1(JL) + ( ZTSK(JL) - PTSKM1M(JL) ) * ZRTMST
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFPP_CTL_MOD:SURFPP_CTL',1,ZHOOK_HANDLE)
CONTAINS

SUBROUTINE COMPUTE_DDH

! DDH diagnostics computation, skin temperature
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFPP_CTL:COMPUTE_DDH',0,ZHOOK_HANDLE)

DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA

    PDHTLS(JL,JTILE,8)=PFRTI(JL,JTILE)*PAHFSTI(JL,JTILE)
    PDHTLS(JL,JTILE,9)=PFRTI(JL,JTILE)*PAHFLTI(JL,JTILE)
    IF (JTILE == 1) THEN
      PDHTLS(JL,JTILE,10)=0.0_JPRB
    ELSE
      PDHTLS(JL,JTILE,10)=PFRTI(JL,JTILE)*PG0TI(JL,JTILE)
    ENDIF

    PDHTLS(JL,JTILE,11)=PFRTI(JL,JTILE)*PEVAPTI(JL,JTILE)
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('SURFPP_CTL:COMPUTE_DDH',1,ZHOOK_HANDLE)

END SUBROUTINE COMPUTE_DDH

END SUBROUTINE SURFPP_CTL
END MODULE SURFPP_CTL_MOD
