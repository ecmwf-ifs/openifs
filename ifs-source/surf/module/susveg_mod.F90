MODULE SUSVEG_MOD
CONTAINS
SUBROUTINE SUSVEG(LD_LELAIV,LD_LECTESSEL,LD_LEAGS,LD_LEFARQUHAR,&
               & LD_LEAIRCO2COUP,PRLAIINT, &
               & YDDIM,YDCST,YDSOIL,YDVEG)

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_DIM   , ONLY : TDIM
USE YOS_CST   , ONLY : TCST
USE YOS_SOIL  , ONLY : TSOIL
USE YOS_VEG   , ONLY : TVEG

USE SRFROOTFR_MOD

! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**   *SUSVEG* IS THE SET-UP ROUTINE FOR COMMON BLOCK *YOS_VEG*

!     PURPOSE
!     -------
!          THIS ROUTINE INITIALIZES THE CONSTANTS IN COMMON BLOCK
!     *YOS_VEG*

!     INTERFACE.
!     ----------
!     CALLLED FRPM *SUSURF*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     A.C.M. BELJAARS         E.C.M.W.F.       2/11/89

!     MODIFICATIONS
!     -------------
!     J.-J. MORCRETTE  ECMWF      91/07/14
!     A. BELJAARS      ECMWF      99/01/05    SURFACE TILES
!     M. HAMRUD        ECMWF      10/11/00    FIX FOR POT. OVERWRITE
!     M.Hamrud                    01/10/2003  CY28 Cleaning
!     P.Viterbo                   24/10/2004  surf library
!     P. Viterbo       ECMWF      09/06/2005  move in surf vob
!     A. BELJAARS      ECMWF      03/12/05    roughness length tables for TOFD
!     E. Dutra/G.Balsamo          01/05/08    add lake tile
!     S. Boussetta/G.Balsamo May 2009 Add lai
!     S. Boussetta/G.Balsamo May 2010 Add CTESSEL for CO2
!     I. Sandu                    May 2010   modify momentum and heat roughness length tables 
!     S. Boussetta/G.Balsamo June 2011 Add LEAGS for modularity of CO2&Evap
!     I. Sandu                    March 2013   strength of coupling in unstable cond for high veg 
!     I. Sandu    24-02-2014  Lambda skin values by vegetation type instead of tile
!     G. Balsamo/I. Sandu    February 2016    reduced coupling over forests
!     A. Agusti-Panareda Nov 2020 Add LEAIRCO2COUP for use of variable atm CO2 in photosynthesis
!     A. Agusti-Panareda Jul 2021 Add LEFARQUHAR for use of Farquhar photosynthesis model
!     ------------------------------------------------------------------


IMPLICIT NONE

LOGICAL,         INTENT(IN)    :: LD_LELAIV
LOGICAL,         INTENT(IN)    :: LD_LECTESSEL
LOGICAL,         INTENT(IN)    :: LD_LEAGS
LOGICAL,         INTENT(IN)    :: LD_LEFARQUHAR
LOGICAL,         INTENT(IN)    :: LD_LEAIRCO2COUP
REAL(KIND=JPRB), INTENT(IN)    :: PRLAIINT
TYPE(TDIM),      INTENT(IN)    :: YDDIM
TYPE(TCST),      INTENT(IN)    :: YDCST
TYPE(TSOIL),     INTENT(IN)    :: YDSOIL
TYPE(TVEG),      INTENT(INOUT) :: YDVEG

INTEGER(KIND=JPIM) ::  ITILES, IVTYPES, JS
REAL(KIND=JPRB)    :: ZLARGE, ZSNOW, ZCONV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSVEG_MOD:SUSVEG',0,ZHOOK_HANDLE)
ASSOCIATE(RETV=>YDCST%RETV, &
 & NCSS=>YDDIM%NCSS, NTILES=>YDDIM%NTILES, &
 & RDAW=>YDSOIL%RDAW, &
 & LEAGS=>YDVEG%LEAGS,LEAIRCO2COUP=>YDVEG%LEAIRCO2COUP,LECTESSEL=>YDVEG%LECTESSEL, LEFARQUHAR=>YDVEG%LEFARQUHAR, LELAIV=>YDVEG%LELAIV, &
 & NVTILES=>YDVEG%NVTILES, NVTYPES=>YDVEG%NVTYPES, RCEPSW=>YDVEG%RCEPSW, &
 & RCVC=>YDVEG%RCVC, REPEVAP=>YDVEG%REPEVAP, REPSR=>YDVEG%REPSR, &
 & RLAIINT=>YDVEG%RLAIINT, RLHAERO=>YDVEG%RLHAERO, RLHAEROS=>YDVEG%RLHAEROS, &
 & RVINTER=>YDVEG%RVINTER, RVLT=>YDVEG%RVLT, RVRAD=>YDVEG%RVRAD)

! Number of vegetation types
NVTYPES=20

! Number of vegetation TILES h/l
NVTILES=2

IVTYPES = MAX(NVTYPES,20)

! Other constants
REPEVAP=1.E-10_JPRB
RVINTER=0.5_JPRB

! Vegetation cover
 
! 0 value acts as default for vegetation type, simplifying code ...
IF(.NOT.ALLOCATED(YDVEG%RVCOV)) ALLOCATE(YDVEG%RVCOV(0:IVTYPES))
YDVEG%RVCOV(1)=0.9_JPRB     ! Crops, Mixed Farming
YDVEG%RVCOV(2)=0.85_JPRB    ! Short Grass
YDVEG%RVCOV(3)=0.9_JPRB     ! Evergreen Needleleaf Trees
YDVEG%RVCOV(4)=0.9_JPRB     ! Deciduous Needleleaf Trees
YDVEG%RVCOV(5)=0.9_JPRB     ! Deciduous Broadleaf Trees
YDVEG%RVCOV(6)=0.99_JPRB    ! Evergreen Broadleaf Trees
YDVEG%RVCOV(7)=0.7_JPRB     ! Tall Grass
YDVEG%RVCOV(8)=0.0_JPRB       ! Desert
YDVEG%RVCOV(9)=0.5_JPRB       ! Tundra
YDVEG%RVCOV(10)=0.9_JPRB    ! Irrigated Crops
YDVEG%RVCOV(11)=0.1_JPRB    ! Semidesert
YDVEG%RVCOV(12)=0.0_JPRB      ! Ice Caps and Glaciers
YDVEG%RVCOV(13)=0.6_JPRB    ! Bogs and Marshes
YDVEG%RVCOV(14)=0.0_JPRB      ! Inland Water
YDVEG%RVCOV(15)=0.0_JPRB      ! Ocean
YDVEG%RVCOV(16)=0.5_JPRB      ! Evergreen Shrubs
YDVEG%RVCOV(17)=0.5_JPRB      ! Deciduous Shrubs
YDVEG%RVCOV(18)=0.9_JPRB    ! Mixed Forest/woodland
YDVEG%RVCOV(19)=0.9_JPRB    ! Interrupted Forest
YDVEG%RVCOV(20)=0.6_JPRB    ! Water and Land Mixtures
YDVEG%RVCOV(0)=YDVEG%RVCOV(8)

! Leaf area index

IF(.NOT.ALLOCATED(YDVEG%RVLAI)) ALLOCATE(YDVEG%RVLAI(0:IVTYPES))
YDVEG%RVLAI(1)=3._JPRB     ! Crops, Mixed Farming
YDVEG%RVLAI(2)=2._JPRB     ! Short Grass
YDVEG%RVLAI(3)=5._JPRB     ! Evergreen Needleleaf Trees
YDVEG%RVLAI(4)=5._JPRB     ! Deciduous Needleleaf Trees
YDVEG%RVLAI(5)=5._JPRB     ! Deciduous Broadleaf Trees
YDVEG%RVLAI(6)=6._JPRB     ! Evergreen Broadleaf Trees
YDVEG%RVLAI(7)=2._JPRB     ! Tall Grass
YDVEG%RVLAI(8)=0.5_JPRB      ! Desert
YDVEG%RVLAI(9)=1.0_JPRB       ! Tundra
YDVEG%RVLAI(10)=3._JPRB    ! Irrigated Crops
YDVEG%RVLAI(11)=0.5_JPRB     ! Semidesert
YDVEG%RVLAI(12)=0.0_JPRB     ! Ice Caps and Glaciers
YDVEG%RVLAI(13)=4._JPRB    ! Bogs and Marshes
YDVEG%RVLAI(14)=0.0_JPRB     ! Inland Water
YDVEG%RVLAI(15)=0.0_JPRB     ! Ocean
YDVEG%RVLAI(16)=3._JPRB    ! Evergreen Shrubs
YDVEG%RVLAI(17)=1.5_JPRB   ! Deciduous Shrubs
YDVEG%RVLAI(18)=5._JPRB    ! Mixed Forest/woodland
YDVEG%RVLAI(19)=2.5_JPRB   ! Interrupted Forest
YDVEG%RVLAI(20)=4._JPRB    ! Water and Land Mixtures
YDVEG%RVLAI(0)=YDVEG%RVLAI(8)

LELAIV = LD_LELAIV
LECTESSEL = LD_LECTESSEL
LEAGS = LD_LEAGS
LEFARQUHAR = LD_LEFARQUHAR
LEAIRCO2COUP = LD_LEAIRCO2COUP
RLAIINT=PRLAIINT

! Skin layer conductivity

ZLARGE=1.E10_JPRB
!attention these values are still used in vevap for tile 7
RLHAERO=15._JPRB
RLHAEROS=10._JPRB

! Unstable skin layer conductivity
IF(.NOT.ALLOCATED(YDVEG%RVLAMSK)) ALLOCATE(YDVEG%RVLAMSK(0:IVTYPES))
YDVEG%RVLAMSK(1)=10._JPRB    ! Crops, Mixed Farming
YDVEG%RVLAMSK(2)=10._JPRB    ! Short Grass
YDVEG%RVLAMSK(3)=10._JPRB    ! Evergreen Needleleaf Trees
YDVEG%RVLAMSK(4)=10._JPRB    ! Deciduous Needleleaf Trees
YDVEG%RVLAMSK(5)=10._JPRB    ! Deciduous Broadleaf Trees
YDVEG%RVLAMSK(6)=10._JPRB    ! Evergreen Broadleaf Trees
YDVEG%RVLAMSK(7)=10._JPRB    ! Tall Grass
YDVEG%RVLAMSK(8)=15._JPRB    ! Desert
YDVEG%RVLAMSK(9)=10._JPRB     ! Tundra
YDVEG%RVLAMSK(10)=10._JPRB   ! Irrigated Crops
YDVEG%RVLAMSK(11)=10._JPRB   ! Semidesert
YDVEG%RVLAMSK(12)=58._JPRB      ! Ice Caps and Glaciers - the value we had for the snow tile
YDVEG%RVLAMSK(13)=10._JPRB   ! Bogs and Marshes
YDVEG%RVLAMSK(14)=ZLARGE      ! Inland Water
YDVEG%RVLAMSK(15)=ZLARGE      ! Ocean
YDVEG%RVLAMSK(16)=10._JPRB   ! Evergreen Shrubs
YDVEG%RVLAMSK(17)=10._JPRB   ! Deciduous Shrubs
YDVEG%RVLAMSK(18)=10._JPRB   ! Mixed Forest/woodland
YDVEG%RVLAMSK(19)=10._JPRB   ! Interrupted Forest
YDVEG%RVLAMSK(20)=ZLARGE   ! Water and Land Mixtures
YDVEG%RVLAMSK(0)=YDVEG%RVLAMSK(8)

! Stable skin layer conductivity
IF(.NOT.ALLOCATED(YDVEG%RVLAMSKS)) ALLOCATE(YDVEG%RVLAMSKS(0:IVTYPES))
YDVEG%RVLAMSKS(1)=10._JPRB    ! Crops, Mixed Farming
YDVEG%RVLAMSKS(2)=10._JPRB    ! Short Grass
YDVEG%RVLAMSKS(3)=10._JPRB    ! Evergreen Needleleaf Trees
YDVEG%RVLAMSKS(4)=10._JPRB    ! Deciduous Needleleaf Trees
YDVEG%RVLAMSKS(5)=10._JPRB    ! Deciduous Broadleaf Trees
YDVEG%RVLAMSKS(6)=10._JPRB    ! Evergreen Broadleaf Trees
YDVEG%RVLAMSKS(7)=10._JPRB    ! Tall Grass
YDVEG%RVLAMSKS(8)=15._JPRB    ! Desert
YDVEG%RVLAMSKS(9)=10._JPRB     ! Tundra
YDVEG%RVLAMSKS(10)=10._JPRB   ! Irrigated Crops
YDVEG%RVLAMSKS(11)=10._JPRB   ! Semidesert
YDVEG%RVLAMSKS(12)=58._JPRB      ! Ice Caps and Glaciers - the value we had for the ice tile
YDVEG%RVLAMSKS(13)=10._JPRB   ! Bogs and Marshes
YDVEG%RVLAMSKS(14)=ZLARGE      ! Inland Water
YDVEG%RVLAMSKS(15)=ZLARGE      ! Ocean
YDVEG%RVLAMSKS(16)=10._JPRB   ! Evergreen Shrubs
YDVEG%RVLAMSKS(17)=10._JPRB   ! Deciduous Shrubs
YDVEG%RVLAMSKS(18)=10._JPRB   ! Mixed Forest/woodland
YDVEG%RVLAMSKS(19)=10._JPRB   ! Interrupted Forest
YDVEG%RVLAMSKS(20)=ZLARGE   ! Water and Land Mixtures
YDVEG%RVLAMSKS(0)=YDVEG%RVLAMSKS(8)

! Transmission of net solar rad. through vegetation
IF(.NOT.ALLOCATED(YDVEG%RVTRSR)) ALLOCATE(YDVEG%RVTRSR(0:IVTYPES))
YDVEG%RVTRSR(1)=0.05_JPRB    ! Crops, Mixed Farming
YDVEG%RVTRSR(2)=0.05_JPRB    ! Short Grass
YDVEG%RVTRSR(3)=0.03_JPRB    ! Evergreen Needleleaf Trees
YDVEG%RVTRSR(4)=0.03_JPRB    ! Deciduous Needleleaf Trees
YDVEG%RVTRSR(5)=0.03_JPRB    ! Deciduous Broadleaf Trees
YDVEG%RVTRSR(6)=0.03_JPRB    ! Evergreen Broadleaf Trees
YDVEG%RVTRSR(7)=0.05_JPRB    ! Tall Grass
YDVEG%RVTRSR(8)=0.00_JPRB    ! Desert
YDVEG%RVTRSR(9)=0.05_JPRB     ! Tundra
YDVEG%RVTRSR(10)=0.05_JPRB   ! Irrigated Crops
YDVEG%RVTRSR(11)=0.05_JPRB   ! Semidesert
YDVEG%RVTRSR(12)=0.00_JPRB      ! Ice Caps and Glaciers - the value we had for the snow tile
YDVEG%RVTRSR(13)=0.05_JPRB   ! Bogs and Marshes
YDVEG%RVTRSR(14)=0.00_JPRB     ! Inland Water
YDVEG%RVTRSR(15)=0.00_JPRB      ! Ocean
YDVEG%RVTRSR(16)=0.05_JPRB   ! Evergreen Shrubs
YDVEG%RVTRSR(17)=0.05_JPRB   ! Deciduous Shrubs
YDVEG%RVTRSR(18)=0.03_JPRB   ! Mixed Forest/woodland
YDVEG%RVTRSR(19)=0.03_JPRB   ! Interrupted Forest
YDVEG%RVTRSR(20)=0.00_JPRB   ! Water and Land Mixtures
YDVEG%RVTRSR(0)=YDVEG%RVTRSR(8)



! Root fraction

IF(.NOT.ALLOCATED(YDVEG%RVROOTSA)) ALLOCATE(YDVEG%RVROOTSA(NCSS,0:NVTYPES))
IF (NCSS >= 1) CALL SRFROOTFR(NCSS,NVTYPES,RDAW,YDVEG%RVROOTSA(:,1:NVTYPES))
DO JS=1,NCSS
  YDVEG%RVROOTSA(JS,0)=YDVEG%RVROOTSA(JS,8)
ENDDO

! Set other constants

RCEPSW  =1.E-3_JPRB

! Minimum stomatal resitance for each vegetation type (s/m)

!                 The additinal correction coefficient (1.30 etc.) accounts
!                 for the decreased radiation stress function at full 
!                 light saturation.

IF(.NOT.ALLOCATED(YDVEG%RVRSMIN)) ALLOCATE (YDVEG%RVRSMIN(0:IVTYPES))
!RVRSMIN(1)=180._JPRB    ! Crops, Mixed Farming
!RVRSMIN(2)=110._JPRB    ! Short Grass
YDVEG%RVRSMIN(1)=100._JPRB    ! Crops, Mixed Farming
YDVEG%RVRSMIN(2)=100._JPRB    ! Short Grass
!RVRSMIN(3)=500._JPRB    ! Evergreen Needleleaf Trees
!RVRSMIN(4)=500._JPRB    ! Deciduous Needleleaf Trees
YDVEG%RVRSMIN(3)=250._JPRB    ! Evergreen Needleleaf Trees
YDVEG%RVRSMIN(4)=250._JPRB    ! Deciduous Needleleaf Trees
YDVEG%RVRSMIN(5)=175._JPRB    ! Deciduous Broadleaf Trees
YDVEG%RVRSMIN(6)=240._JPRB    ! Evergreen Broadleaf Trees
YDVEG%RVRSMIN(7)=100._JPRB    ! Tall Grass
YDVEG%RVRSMIN(8)=250._JPRB    ! Desert
YDVEG%RVRSMIN(9)=80._JPRB     ! Tundra
!RVRSMIN(10)=180._JPRB   ! Irrigated Crops
YDVEG%RVRSMIN(10)=100._JPRB   ! Irrigated Crops
YDVEG%RVRSMIN(11)=150._JPRB   ! Semidesert
YDVEG%RVRSMIN(12)=0.0_JPRB      ! Ice Caps and Glaciers
YDVEG%RVRSMIN(13)=240._JPRB   ! Bogs and Marshes
YDVEG%RVRSMIN(14)=0.0_JPRB      ! Inland Water
YDVEG%RVRSMIN(15)=0.0_JPRB      ! Ocean
YDVEG%RVRSMIN(16)=225._JPRB   ! Evergreen Shrubs
YDVEG%RVRSMIN(17)=225._JPRB   ! Deciduous Shrubs
YDVEG%RVRSMIN(18)=250._JPRB   ! Mixed Forest/woodland
YDVEG%RVRSMIN(19)=175._JPRB   ! Interrupted Forest
YDVEG%RVRSMIN(20)=150._JPRB   ! Water and Land Mixtures
YDVEG%RVRSMIN(0)=YDVEG%RVRSMIN(8)

! Parameter in humidity stress function (m/s mbar, converted to M/S kgkg-1)
IF(.NOT.ALLOCATED(YDVEG%RVHSTR)) ALLOCATE (YDVEG%RVHSTR(0:IVTYPES))
ZCONV=1013.25_JPRB*(RETV+1)
YDVEG%RVHSTR(1)=0.0_JPRB     ! Crops, Mixed Farming
YDVEG%RVHSTR(2)=0.0_JPRB     ! Short Grass
YDVEG%RVHSTR(3)=0.03_JPRB  ! Evergreen Needleleaf Trees
YDVEG%RVHSTR(4)=0.03_JPRB  ! Deciduous Needleleaf Trees
YDVEG%RVHSTR(5)=0.03_JPRB  ! Deciduous Broadleaf Trees
YDVEG%RVHSTR(6)=0.03_JPRB  ! Evergreen Broadleaf Trees
YDVEG%RVHSTR(7)=0.0_JPRB     ! Tall Grass
YDVEG%RVHSTR(8)=0.0_JPRB     ! Desert
YDVEG%RVHSTR(9)=0.0_JPRB     ! Tundra
YDVEG%RVHSTR(10)=0.0_JPRB    ! Irrigated Crops
YDVEG%RVHSTR(11)=0.0_JPRB    ! Semidesert
YDVEG%RVHSTR(12)=0.0_JPRB    ! Ice Caps and Glaciers
YDVEG%RVHSTR(13)=0.0_JPRB    ! Bogs and Marshes
YDVEG%RVHSTR(14)=0.0_JPRB    ! Inland Water
YDVEG%RVHSTR(15)=0.0_JPRB    ! Ocean
YDVEG%RVHSTR(16)=0.0_JPRB    ! Evergreen Shrubs
YDVEG%RVHSTR(17)=0.0_JPRB    ! Deciduous Shrubs
YDVEG%RVHSTR(18)=0.03_JPRB ! Mixed Forest/woodland
YDVEG%RVHSTR(19)=0.03_JPRB ! Interrupted Forest
YDVEG%RVHSTR(20)=0.0_JPRB    ! Water and Land Mixtures
YDVEG%RVHSTR(0)=YDVEG%RVHSTR(8)
YDVEG%RVHSTR(:)=ZCONV*YDVEG%RVHSTR(:)

! Roughness length for momentum (Mahfouf et al. 1995)

IF(.NOT.ALLOCATED(YDVEG%RVZ0M)) ALLOCATE (YDVEG%RVZ0M(0:IVTYPES))
YDVEG%RVZ0M(1)=0.25_JPRB     ! Crops, Mixed Farming
YDVEG%RVZ0M(2)=0.1_JPRB     ! Short Grass
YDVEG%RVZ0M(3)=2.00_JPRB     ! Evergreen Needleleaf Trees
YDVEG%RVZ0M(4)=2.00_JPRB     ! Deciduous Needleleaf Trees
YDVEG%RVZ0M(5)=2.00_JPRB     ! Deciduous Broadleaf Trees
YDVEG%RVZ0M(6)=2.00_JPRB     ! Evergreen Broadleaf Trees
YDVEG%RVZ0M(7)=0.47_JPRB     ! Tall Grass
YDVEG%RVZ0M(8)=0.013_JPRB    ! Desert                    # Masson et al.
YDVEG%RVZ0M(9)=0.034_JPRB     ! Tundra
YDVEG%RVZ0M(10)=0.5_JPRB    ! Irrigated Crops           # Crops type 1
YDVEG%RVZ0M(11)=0.17_JPRB    ! Semidesert 
YDVEG%RVZ0M(12)=0.0013_JPRB  ! Ice Caps and Glaciers     # Mason et al. 
YDVEG%RVZ0M(13)=0.5_JPRB    ! Bogs and Marshes
YDVEG%RVZ0M(14)=0.0001_JPRB  ! Inland Water              # Not used but needs value here
YDVEG%RVZ0M(15)=0.0001_JPRB  ! Ocean                     # Not used but needs value here
YDVEG%RVZ0M(16)=0.1_JPRB    ! Evergreen Shrubs
YDVEG%RVZ0M(17)=0.25_JPRB    ! Deciduous Shrubs
YDVEG%RVZ0M(18)=2.00_JPRB    ! Mixed Forest/woodland
YDVEG%RVZ0M(19)=1.1_JPRB    ! Interrupted Forest        # New value invented here
YDVEG%RVZ0M(20)=0.02_JPRB    ! Water and Land Mixtures   # Not used but needs value here
YDVEG%RVZ0M(0)=YDVEG%RVZ0M(8)      !                           # Bare soil value

!
! Roughness length for heat

IF(.NOT.ALLOCATED(YDVEG%RVZ0H)) ALLOCATE (YDVEG%RVZ0H(0:IVTYPES))
YDVEG%RVZ0H(1)=YDVEG%RVZ0M( 1)/100._JPRB     ! Crops, Mixed Farming
YDVEG%RVZ0H(2)=YDVEG%RVZ0M( 2)/100._JPRB     ! Short Grass
YDVEG%RVZ0H(3)=YDVEG%RVZ0M( 3)              ! Evergreen Needleleaf Trees
YDVEG%RVZ0H(4)=YDVEG%RVZ0M( 4)              ! Deciduous Needleleaf Trees
YDVEG%RVZ0H(5)=YDVEG%RVZ0M( 5)              ! Deciduous Broadleaf Trees
YDVEG%RVZ0H(6)=YDVEG%RVZ0M( 6)              ! Evergreen Broadleaf Trees
YDVEG%RVZ0H(7)=YDVEG%RVZ0M( 7)/100._JPRB     ! Tall Grass
YDVEG%RVZ0H(8)=YDVEG%RVZ0M( 8)/100._JPRB     ! Desert 
YDVEG%RVZ0H(9)=YDVEG%RVZ0M( 9)/100._JPRB     ! Tundra
YDVEG%RVZ0H(10)=YDVEG%RVZ0M(10)/100._JPRB    ! Irrigated Crops 
YDVEG%RVZ0H(11)=YDVEG%RVZ0M(11)/100._JPRB    ! Semidesert 
YDVEG%RVZ0H(12)=YDVEG%RVZ0M(12)/10._JPRB    ! Ice Caps and Glaciers 
YDVEG%RVZ0H(13)=YDVEG%RVZ0M(13)/100._JPRB    ! Bogs and Marshes
YDVEG%RVZ0H(14)=YDVEG%RVZ0M(14)/10._JPRB    ! Inland Water  
YDVEG%RVZ0H(15)=YDVEG%RVZ0M(15)/10._JPRB    ! Ocean 
YDVEG%RVZ0H(16)=YDVEG%RVZ0M(16)/100._JPRB    ! Evergreen Shrubs
YDVEG%RVZ0H(17)=YDVEG%RVZ0M(17)/100._JPRB    ! Deciduous Shrubs
YDVEG%RVZ0H(18)=YDVEG%RVZ0M(18)             ! Mixed Forest/woodland
YDVEG%RVZ0H(19)=YDVEG%RVZ0M(19)/100._JPRB    ! Interrupted Forest 
YDVEG%RVZ0H(20)=YDVEG%RVZ0M(20)/10._JPRB    ! Water and Land Mixtures 
YDVEG%RVZ0H(0)=YDVEG%RVZ0H(8)        

! Simplified physics constants

RCVC=240._JPRB
RVLT=4._JPRB
RVRAD  =0.5_JPRB
REPSR  =1.E-10_JPRB

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSVEG_MOD:SUSVEG',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE SUSVEG
END MODULE SUSVEG_MOD
