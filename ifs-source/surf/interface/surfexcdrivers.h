INTERFACE
SUBROUTINE SURFEXCDRIVERS    ( YDSURF, &
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
!                                not use in corresponding TL/AD
!    G. Balsamo       03/07/2005 Add soil type
!    S. Boussetta/G.Balsamo May 2009 Add lai
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
!      PCVL     :    Low vegetation fraction
!      PCVH     :    High vegetation fraction
!      PLAIL    :    Low vegetation LAI
!      PLAIH    :    High vegetation LAI

!  Logical:
!      LDSURF2  :    TRUE when simplified surface scheme called

!    Reals with tile index (In): 
!      PFRTI    :    TILE FRACTIONS                                   (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!      PALBTI   :    Tile albedo                                      (0-1)

!    Reals independent of tiles (In):
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
!         *VUPDZ0*
!         *VSURF*
!         *VEXCS*
!         *VEVAP*

!  DOCUMENTATION:
!    See Physics Volume of IFS documentation

!------------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

TYPE(C_PTR)       ,INTENT(IN)    :: YDSURF
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAIL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAIH(:)
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


!------------------------------------------------------------------------

END SUBROUTINE SURFEXCDRIVERS

END INTERFACE
