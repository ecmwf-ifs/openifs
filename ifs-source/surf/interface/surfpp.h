INTERFACE
SUBROUTINE SURFPP( YDSURF,KIDIA,KFDIA,KLON,KTILES, KDHVTLS, KDHFTLS &
 & , PTSTEP &
! input
 & , PFRTI, PAHFLTI, PG0TI, PSTRTULEV, PSTRTVLEV, PTSKM1M &
 & , PUMLEV, PVMLEV, PQMLEV, PGEOMLEV, PCPTSPP ,PCPTGZLEV &
 & , PAPHMS, PZ0MW, PZ0HW, PZ0QW, PZDL, PQSAPP, PBLEND, PFBLEND, PBUOM &
 & , PZ0M, PEVAPSNW, PSSRFLTI, PSLRFL, PSST &
 & , PUCURR, PVCURR, PUSTOKES, PVSTOKES &
! updated
 & , PAHFSTI, PEVAPTI, PTSKE1, PTSKTIP1 &
! output
 & , PDIFTSLEV, PDIFTQLEV, PUSTRTI, PVSTRTI, PTSKTI, PAHFLEV, PAHFLSB, PFWSB  &
 & , PU10M, PV10M, PT2M, PD2M, PQ2M &
 & , PGUST, P10NU, P10NV, PUST &
! output DDH
 & , PDHTLS &
 & , PRPLRG)

USE PARKIND1, ONLY : JPIM, JPRB
USE, INTRINSIC :: ISO_C_BINDING

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
!      PFRTI     :    TILE FRACTIONS                                   (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
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
!      PUCURR    :  U-comp of ocean surface current                  m/s
!      PVCURR    :  V-comp of ocean surface current                  m/s
!      PUSTOKES  :  U-comp of surface Stokes velocity                m/s
!      PVSTOKES  :  V-comp of surface Stokes velocity                m/s

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

!  DOCUMENTATION:
!    See Physics Volume of IFS documentation

!------------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

TYPE(C_PTR)       ,INTENT(IN)    :: YDSURF
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
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKE1(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKTIP1(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTSLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTQLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUSTRTI(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVSTRTI(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTSKTI(:,:)
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

END SUBROUTINE SURFPP
END INTERFACE
