
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!------------------------------------------------------------------------

!  PURPOSE:
!    Routine SURFPPS controls the computation of quantities at the end of
!     vertical diffusion, includting routines to post-process weather elements
!     and gustiness.

!  SURFPPS is called by VDFMAINS

!  METHOD:
!    This routine is a shell needed by the surface library  externalisation.

!  AUTHOR:
!    P. Viterbo       ECMWF May 2005

!  REVISION HISTORY:


!  INTERFACE: 

!    Integers (In):
!      KIDIA    :    Begin point in arrays
!      KFDIA    :    End point in arrays
!      KLON     :    Length of arrays
!      KTILES   :    Number of files


!    Reals (In):
!      PTSTEP    :  Timestep                                          s
!      PFRTI     :  TILE FRACTIONS                                   (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!      PTSKTIP1  :  Tile skin temperature, t+1                       K
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
!      PQSAPP    :  Apparent surface humidity                        kg/kg
!      PBUOM     :  Buoyancy flux, for post-processing of gustiness  ????

!    Reals (Updated):
!      PAHFSTI   :  SURFACE SENSIBLE HEAT FLUX                       W/m2
!      PEVAPTI   :  SURFACE MOISTURE FLUX                            kg/m2/s
!      PTSKE1    :  SKIN TEMPERATURE TENDENCY                        K/s

!    Reals (Out):
!      PDIFTSLEV :  TURBULENT FLUX OF HEAT                           J/(m2*s)
!      PDIFTQLEV :  TURBULENT FLUX OF SPECIFIC HUMIDITY              kg/(m2*s)
!      PUSTRTI   :  SURFACE U-STRESS                                 N/m2 
!      PVSTRTI   :  SURFACE V-STRESS                                 N/m2 
!      PTSKTI    :  SKIN TEMPERATURE                                 K
!      PU10M     :  U-COMPONENT WIND AT 10 M                         m/s
!      PV10M     :  V-COMPONENT WIND AT 10 M                         m/s
!      P10NU     :  U-COMPONENT NEUTRAL WIND AT 10 M                 m/s
!      P10NV     :  V-COMPONENT NEUTRAL WIND AT 10 M                 m/s
!      PT2M      :  TEMPERATURE AT 2M                                K
!      PD2M      :  DEW POINT TEMPERATURE AT 2M                      K
!      PQ2M      :  SPECIFIC HUMIDITY AT 2M                          kg/kg

!     EXTERNALS.
!     ----------

!     ** SURFPPS_CTL CALLS SUCCESSIVELY:
!         *SPPCFL*

!  DOCUMENTATION:
!    See Physics Volume of IFS documentation

!------------------------------------------------------------------------

MODULE SURFPPS_CTL_MOD
CONTAINS
SUBROUTINE SURFPPS_CTL( KIDIA,KFDIA,KLON,KTILES,LDPPCFLS &
 & , PTSTEP &
! input
 & , PFRTI, PTSKTIP1, PAHFLTI, PG0TI, PSTRTULEV, PSTRTVLEV, PTSKM1M &
 & , PUMLEV, PVMLEV, PQMLEV, PGEOMLEV, PCPTSPP ,PCPTGZLEV &
 & , PAPHMS, PZ0MW, PZ0HW, PZ0QW, PQSAPP, PBUOM &
 & , YDCST, YDEXC &
! input, tile dependent
 & , PZ0MTIW, PZ0HTIW, PZ0QTIW, PQSAPPTI, PCPTSPPTI, PBUOMTI &
! updated
 & , PAHFSTI, PEVAPTI, PTSKE1 &
! output
 & , PDIFTSLEV, PDIFTQLEV, PUSTRTI, PVSTRTI, PTSKTI &
 & , PU10M, PV10M, PT2M, PD2M, PQ2M , P10NU, P10NV &
 & )

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST  , ONLY : TCST
USE YOS_EXC  , ONLY : TEXC

USE SPPCFLS_MOD

!------------------------------------------------------------------------

!  PURPOSE:
!    Routine SURFPPS controls the computation of quantities at the end of
!     vertical diffusion, includting routines to post-process weather elements
!     and gustiness.

!  SURFPPS is called by VDFMAINS

!  METHOD:
!    This routine is a shell needed by the surface library  externalisation.

!  AUTHOR:
!    P. Viterbo       ECMWF May 2005

!  REVISION HISTORY:


!  INTERFACE: 

!    Integers (In):
!      KIDIA    :    Begin point in arrays
!      KFDIA    :    End point in arrays
!      KLON     :    Length of arrays
!      KTILES   :    Number of files


!    Reals (In):
!      PTSTEP    :  Timestep                                          s
!      PFRTI     :  TILE FRACTIONS                                   (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!      PTSKTIP1  :  Tile skin temperature, t+1                       K
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
!      PQSAPP    :  Apparent surface humidity                        kg/kg
!      PBUOM     :  Buoyancy flux, for post-processing of gustiness  ????

!    Reals (Updated):
!      PAHFSTI   :  SURFACE SENSIBLE HEAT FLUX                       W/m2
!      PEVAPTI   :  SURFACE MOISTURE FLUX                            kg/m2/s
!      PTSKE1    :  SKIN TEMPERATURE TENDENCY                        K/s

!    Reals (Out):
!      PDIFTSLEV :  TURBULENT FLUX OF HEAT                           J/(m2*s)
!      PDIFTQLEV :  TURBULENT FLUX OF SPECIFIC HUMIDITY              kg/(m2*s)
!      PUSTRTI   :  SURFACE U-STRESS                                 N/m2 
!      PVSTRTI   :  SURFACE V-STRESS                                 N/m2 
!      PTSKTI    :  SKIN TEMPERATURE                                 K
!      PU10M     :  U-COMPONENT WIND AT 10 M                         m/s
!      PV10M     :  V-COMPONENT WIND AT 10 M                         m/s
!      P10NU     :  U-COMPONENT NEUTRAL WIND AT 10 M                 m/s
!      P10NV     :  V-COMPONENT NEUTRAL WIND AT 10 M                 m/s
!      PT2M      :  TEMPERATURE AT 2M                                K
!      PD2M      :  DEW POINT TEMPERATURE AT 2M                      K
!      PQ2M      :  SPECIFIC HUMIDITY AT 2M                          kg/kg

!     EXTERNALS.
!     ----------

!     ** SURFPPS_CTL CALLS SUCCESSIVELY:
!         *SPPCFL*

!  DOCUMENTATION:
!    See Physics Volume of IFS documentation

!------------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES
LOGICAL           ,INTENT(IN)    :: LDPPCFLS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKTIP1(:,:)
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSAPP(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBUOM(:)
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TEXC)        ,INTENT(IN)    :: YDEXC
! Tile dependent pp
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0MTIW(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0HTIW(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0QTIW(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSAPPTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTSPPTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBUOMTI(:,:)

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKE1(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTSLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTQLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUSTRTI(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVSTRTI(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTSKTI(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU10M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PV10M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P10NU(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P10NV(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT2M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD2M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQ2M(:) 

! Local variables
LOGICAL :: LT2MTILE
LOGICAL :: LWIND

INTEGER(KIND=JPIM) :: JTILE, JL

REAL(KIND=JPRB) :: ZTSK(KLON), ZAHFSM(KLON),ZEVAPM(KLON)
REAL(KIND=JPRB) :: ZT2M(KLON,KTILES), ZD2M(KLON,KTILES), ZQ2M(KLON,KTILES)

REAL(KIND=JPRB) :: ZZ0MW(KLON)
REAL(KIND=JPRB) :: ZT2M_DL(KLON), ZQ2M_DL(KLON), ZD2M_DL(KLON)
REAL(KIND=JPRB) :: ZU10M_DUMMY(KLON), ZV10M_DUMMY(KLON), Z10NU_DUMMY(KLON),&
                 & Z10NV_DUMMY(KLON)

REAL(KIND=JPRB) :: ZRTMST
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFPPS_CTL_MOD:SURFPPS_CTL',0,ZHOOK_HANDLE)
ASSOCIATE(RLSTT=>YDCST%RLSTT)

ZRTMST      = 1.0_JPRB/PTSTEP    ! optimization

!*         1.     SURFACE FLUXES - TILES
!                 ----------------------

!*         1.1  SURFACE FLUXES OF HEAT AND MOISTURE FOR THE 
!*              DIFFERENT TILES AND THE MEAN OVER TILES

ZTSK  (KIDIA:KFDIA) = 0.0_JPRB
ZAHFSM(KIDIA:KFDIA) = 0.0_JPRB
ZEVAPM(KIDIA:KFDIA) = 0.0_JPRB
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    ZTSK(JL)=ZTSK(JL)+PFRTI(JL,JTILE)*PTSKTIP1(JL,JTILE)
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
    IF (PFRTI(JL,JTILE) == 0._JPRB) THEN
      PAHFSTI(JL,JTILE)=ZAHFSM(JL)
      PEVAPTI(JL,JTILE)=ZEVAPM(JL)
      PTSKTI(JL,JTILE)=ZTSK(JL)
    ELSE
      PTSKTI(JL,JTILE)=PTSKTIP1(JL,JTILE)
    ENDIF
  ENDDO
ENDDO

!*         1.2  Skin temperature tendencies

DO JL=KIDIA,KFDIA
  PTSKE1(JL) = PTSKE1(JL) + ( ZTSK(JL) - PTSKM1M(JL) ) * ZRTMST
ENDDO

!      2. Post-processing of weather parameters
!         -------------------------------------

IF (LDPPCFLS) THEN
! 2.1 2M DIAGNOSTICS (t2m,d2m,q2m) AVERAGING OVER TILES 
! For ocean tiles (1 and 2), we still use dominant .
  LT2MTILE=.TRUE. 
  IF (LT2MTILE)THEN
    ! Initialise just for safety
    ZT2M(KIDIA:KFDIA,1:KTILES)=0._JPRB
    ZD2M(KIDIA:KFDIA,1:KTILES)=0._JPRB
    ZQ2M(KIDIA:KFDIA,1:KTILES)=0._JPRB
    ! Dummy fields for second call to sppcfl for wind calculation
    ZT2M_DL(KIDIA:KFDIA)=0._JPRB
    ZD2M_DL(KIDIA:KFDIA)=0._JPRB
    ZQ2M_DL(KIDIA:KFDIA)=0._JPRB
    ! Dummy fields for first call to sppcfl for winds
    ZU10M_DUMMY(KIDIA:KFDIA)=0._JPRB
    ZV10M_DUMMY(KIDIA:KFDIA)=0._JPRB
    Z10NU_DUMMY(KIDIA:KFDIA)=0._JPRB
    Z10NV_DUMMY(KIDIA:KFDIA)=0._JPRB
 
    LWIND=.FALSE.
    DO JTILE=1,KTILES
      CALL SPPCFLS(KIDIA,KFDIA,KLON &
       & , PUMLEV, PVMLEV, PQMLEV, PAPHMS, PGEOMLEV, PCPTGZLEV &
       & , PCPTSPPTI(KIDIA:KFDIA,JTILE), PQSAPPTI(KIDIA:KFDIA,JTILE) &
       & , PZ0MTIW(KIDIA:KFDIA,JTILE), PZ0HTIW(KIDIA:KFDIA,JTILE) &
       & , PZ0QTIW(KIDIA:KFDIA,JTILE), PBUOMTI(KIDIA:KFDIA,JTILE) &
       & , YDCST, YDEXC &
       & , ZU10M_DUMMY, ZV10M_DUMMY, Z10NU_DUMMY, Z10NV_DUMMY &
       & , ZT2M(KIDIA:KFDIA,JTILE), ZD2M(KIDIA:KFDIA,JTILE), ZQ2M(KIDIA:KFDIA,JTILE) &
       & , LWIND)
    ENDDO
    ! Compute postprocessed wind over dominant low
    LWIND=.TRUE.
    ZZ0MW(KIDIA:KFDIA)=MAX(PZ0MW(KIDIA:KFDIA),PZ0HW(KIDIA:KFDIA))
    CALL SPPCFLS(KIDIA,KFDIA,KLON &
       & , PUMLEV, PVMLEV, PQMLEV, PAPHMS, PGEOMLEV, PCPTGZLEV &
       & , PCPTSPP, PQSAPP &
       & , ZZ0MW, PZ0HW        &
       & , PZ0QW, PBUOM               &
       & , YDCST, YDEXC               &
       & , PU10M, PV10M, P10NU, P10NV &
       & , ZT2M_DL, ZD2M_DL, ZQ2M_DL  &
       & , LWIND )

    PT2M(KIDIA:KFDIA)=0._JPRB
    PD2M(KIDIA:KFDIA)=0._JPRB
    PQ2M(KIDIA:KFDIA)=0._JPRB
    DO JTILE=1,KTILES
      DO JL=KIDIA,KFDIA
        PT2M(JL)=PT2M(JL)+&
                         &PFRTI(JL,JTILE)*ZT2M(JL,JTILE)
        PD2M(JL)=PD2M(JL)+&
                         &PFRTI(JL,JTILE)*ZD2M(JL,JTILE)
        PQ2M(JL)=PQ2M(JL)+&
                         &PFRTI(JL,JTILE)*ZQ2M(JL,JTILE)
      ENDDO
    ENDDO
    ! Overwrite where dominant water or sea-ice. 
    ! Using PLSM would be cleaner, but this should be eq.
    DO JL=KIDIA,KFDIA
    IF ((PFRTI(JL,3)+PFRTI(JL,4)+PFRTI(JL,5)+PFRTI(JL,6)+PFRTI(JL,7)+PFRTI(JL,8)+PFRTI(JL,9))<1.0_JPRB)THEN
        PT2M(JL)=ZT2M_DL(JL)
        PD2M(JL)=ZD2M_DL(JL)
        PQ2M(JL)=ZQ2M_DL(JL)
      ENDIF
    ENDDO

 
  ELSE ! Business as usual
    ZZ0MW(KIDIA:KFDIA)=MAX(PZ0MW(KIDIA:KFDIA),PZ0HW(KIDIA:KFDIA))
    CALL SPPCFLS(KIDIA,KFDIA,KLON &
     & , PUMLEV, PVMLEV, PQMLEV, PAPHMS, PGEOMLEV, PCPTGZLEV &
     & , PCPTSPP, PQSAPP, ZZ0MW, PZ0HW, PZ0QW, PBUOM &
     !*& , PCPTSPP, PQSAPP, PZ0MW, PZ0MW, PZ0HW, PZ0QW, PBUOM &
     & , YDCST, YDEXC &
     & , PU10M, PV10M, P10NU, P10NV, PT2M, PD2M, PQ2M &
     & , LWIND)
  ENDIF
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFPPS_CTL_MOD:SURFPPS_CTL',1,ZHOOK_HANDLE)

END SUBROUTINE SURFPPS_CTL
END MODULE SURFPPS_CTL_MOD
