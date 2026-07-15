! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DEPVEL_GC (YDMODEL, KIDIA, KFDIA, KLON, KTRAC, KCHEM,  KTILES, PFRTI , PGLAT, PSNS, &
                                  &  PCVL , PCVH ,PLAIL, PLAIH,  KTVL, KTVH, &
                                  &  PCRB, PCRL, PCRLU , PCRH, PCRHS, &
                                  &  PMU0, PRSF,PGEOM1,PGEOH,PHFLUX, &
                                  &  PTS, PT2M, PD2M, PFRSO, PRAQTI, PKCLEV, &  
                                  &  PUSTAR, PUSTARGUST,  PZ0M, PDEPVELCLIM, PDEPVEL )
!!    PURPOSE
!!    -------
!!    The purpose of this routine is to cumpute the deposition velocities of 
!!    the trace gases
!!    Computation is done for a given date, forecast start time, 
!! 
!!   
!!**  METHOD
!!    ------
!!    See [Seinfeld et Pandis, 1998, "Atmospheric Chemistry and
!!    Physics", chap. 19, pp.958-996] based primarily on [Wesely, 1989]. 
!!    Note also the following differences with the Wesely proposed method:
!!    For IFS:   
!!    The aerodynamic resistance is taken from IFS surface scheme  
!!    The stomatal resistance is taken from IFS surface scheme  
!!    Mapping of Wesley classes to IFS surface classes    
!!
!!
!!    REFERENCE
!!    ---------
!!    CAMS_42 Deliverable report D42.3.1.1  (Dec. 2020)
!!
!
!!    AUTHOR
!!    ------
!!   Taken from Geos-Chem code by Doug Finch, November 2020
!!
!!
!!    MODIFICATIONS
!!    -------------
!!    SR 07/2022  Integrate aerosol species
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE DRYDEP_PAR, ONLY : RSMAX 
 
USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD
USE YOMRIP0   ,ONLY : NINDAT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TYPE_MODEL , ONLY : MODEL
! USE YOMLUN, ONLY: NULERR
USE YOMCST, ONLY: RD,RG

IMPLICIT NONE

TYPE(MODEL)        , INTENT(INOUT)       :: YDMODEL
INTEGER(KIND=JPIM) , INTENT(IN)          :: KLON, KTRAC, KCHEM   
INTEGER(KIND=JPIM) , INTENT(IN)          :: KIDIA, KFDIA    
INTEGER(KIND=JPIM) , INTENT(IN)          :: KTILES 


! Tile fraction
REAL(KIND=JPRB) , DIMENSION (KLON,KTILES), INTENT(IN)   ::   PFRTI ! Tile fraction
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL

! surface resistance , tiled  
REAL(KIND=JPRB) , DIMENSION (KLON,KTILES), INTENT(IN)   ::   PRAQTI 
! lowest level transfer coefficient for tracers (reverse of aerodynamic resistance)  
REAL(KIND=JPRB) , DIMENSION (KLON), INTENT(IN)   ::   PKCLEV 

REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PTS   ! surface temperature
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PT2M   ! 2m temperature
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PD2M   ! 2m dew point
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PMU0  ! COS of SZA
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PRSF  ! Pressure at middle of bottom model level
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PGEOM1 ! Geopotential height at middle of bottom model level.
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PGEOH  ! Geopotential height at bottom of bottom model level.
REAL(KIND=JPRB) , DIMENSION (KLON,KTILES),INTENT(IN)     :: PHFLUX  ! Sensible heat flux

REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PCVL ! LOW VEGETATION COVER 
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PCVH !  HIGH VEGETATION COVER 
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PLAIL ! LOW VEGETAION LAI  
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PLAIH ! HIGH VEGETAION LAI  
INTEGER(KIND=JPIM) , DIMENSION (KLON),INTENT(IN)     :: KTVL  ! LOW VEGETAION TYPE   
INTEGER(KIND=JPIM) , DIMENSION (KLON),INTENT(IN)     :: KTVH  ! HIGH VEGETAION TYPE  

REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PCRB ! Bare soil canopy resistance as calculated by surf   
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PCRL ! low vegetation  soil canopy resistance as calculated by surf   
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PCRLU ! low vegetation unstressed  canopy resistance as calculated by surf   
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PCRH !  high vegetation canopy resistance as calculated by surf   
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PCRHS ! high vegetation iunder snow canopy resistance as calculated by surf   


! instantanneous solar radiation flux at the surface
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PFRSO
!  

REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PGLAT  !  Latidtude rad
! REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PGLAM  !  LONGITUDE rad

REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PSNS    ! SNOW
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PUSTAR  ! friction velocity
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PUSTARGUST  ! friction velocity with gusts
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PZ0M    ! Roughness length for momentum M
REAL(KIND=JPRB), DIMENSION (KLON, KTRAC), INTENT(IN)     :: PDEPVELCLIM

! Output 
REAL(KIND=JPRB), DIMENSION (KLON, KTRAC), INTENT(OUT)     :: PDEPVEL

! local 
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZVEG0 ! vegetation cover 
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZRHCL   ! relative humidity
REAL(KIND=JPRB)                         :: ZRSTO1  ! stomatal resistances for water vapor from ifs/surf
REAL(KIND=JPRB)                         :: ZRAERO ! aerodynamic resistances
REAL(KIND=JPRB)                         :: ZWRB   ! laminar resistances
REAL(KIND=JPRB)                         :: ZWRC   ! surface resistances for GEOS-Chem veg types
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZRAIN  ! soil wet or not with rain
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZDEW   ! soil wet or not with dew
! REAL(KIND=JPRB)                         :: ZLON, ZLAT

! not KVTYPES seems to be 2 not 20 ???
REAL(KIND=JPRB), DIMENSION (0:20)        :: ZRSMIN_MOD ! Minimum stomatal resistances according to yom_veg rsmin
REAL(KIND=JPRB), DIMENSION (12)          :: ZRHOP, ZWETD  ! Aerosol density and wet diameter
REAL(KIND=JPRB)                          :: ZRES_TILE ! Temporary resistances
REAL(KIND=JPRB)                          :: ZDEPTILE ! Proportion of veg on tile
REAL(KIND=JPRB)                          :: ZLAI ! Leaf area index as input for surface resistance mod
REAL(KIND=JPRB)                          :: ZXM ! Molecular weight
REAL(KIND=JPRB)                          :: ZAIRDEN ! Air density [kg/m3]
REAL(KIND=JPRB)                          :: ZHEIGHT ! Height  [m]
REAL(KIND=JPRB)                          :: ZSIGMA, ZTMP1, ZTMP2
REAL(KIND=JPRB), DIMENSION (12)          :: ZVFRAC  ! Volume fraction

INTEGER(KIND=JPIM) , DIMENSION (20)      :: IVEG_GC_TYPES! Reference table to map iveg to GEOS-Chem types
INTEGER(KIND=JPIM) :: IVEG_ZH, ISEASON_WE
INTEGER(KIND=JPIM) :: IRH(KLON)
INTEGER(KIND=JPIM)                       :: ILOW_VEG_NUM, IHIGH_VEG_NUM ! Reference number for veg types
INTEGER(KIND=JPIM)                       :: IVEG_GC   ! type of vegetation in reference to GEOS-Chem
INTEGER(KIND=JPIM)                       :: ID_NOWET  ! Reconstructed dry tile 
INTEGER(KIND=JPIM)                       :: IVEGI ! TYPE OF VEGETATION IN IFS  1... 20
INTEGER(KIND=JPIM)  :: ITRCHEM,  ITR1, ITRAERO, JL,  ID , ITILE , IDEBUG, JTAB, IMM, IBIN, ITYP
REAL(KIND=JPHOOK) ::  ZHOOK_HANDLE 

!
!-------------------------------------------------------------------------------

#include "fcttim.func.h"
! #include "ddr_aero_res_gc.intfb.h"
#include "ddr_laminar_res_gc.intfb.h"
#include "ddr_surf_res_gc.intfb.h"
#include "ddr_surf_res_gc_v2.intfb.h"
#include "aer_drydepvel.intfb.h"
#include "aer_drydepvelzh14.intfb.h"
#include "aer_drydepvelem20.intfb.h"
#include "ddr_zh_seasonone.intfb.h"


IF (LHOOK) CALL DR_HOOK('DEPVEL_GC',0,ZHOOK_HANDLE)

ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL, YCHEM=>YDMODEL%YRML_GCONF%YGFL%YCHEM, &
   & KCHEM_DRYDEP=>YDMODEL%YRML_CHEM%YRCHEM%KCHEM_DRYDEP, &
   & YDEAERSNK=>YDMODEL%YRML_PHY_AER%YREAERSNK, &
   & YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM)
ASSOCIATE(YAERO_DESC=>YDEAERATM%YAERO_DESC, &
          & NDRYDEPVEL_DYN => YDEAERSNK%NDRYDEPVEL_DYN, &
          & YDRYDEP => YDMODEL%YRML_CHEM%YRDRYDEP, &
          & RRHTAB=>YDEAERSNK%RRHTAB,                               &
          & LAERDUST_NEWBIN=>YDEAERATM%LAERDUST_NEWBIN, &
          & RSSDENS_RHTAB=>YDEAERSNK%RSSDENS_RHTAB,                 &
           & RRHO_DD=>YDEAERSNK%RRHO_DD, RRHO_SS=>YDEAERSNK%RRHO_SS, &
          & RSSGROWTH_RHTAB=>YDEAERSNK%RSSGROWTH_RHTAB,             &
          & RSOAGROWTH_RHTAB=>YDEAERSNK%RSOAGROWTH_RHTAB,             &
          & RSO4GROWTH_RHTAB=>YDEAERSNK%RSO4GROWTH_RHTAB,             &
          & ROMGROWTH_RHTAB=>YDEAERSNK%ROMGROWTH_RHTAB,             &
          & RNIGROWTH_RHTAB=>YDEAERSNK%RNIGROWTH_RHTAB,             &
          & RAMGROWTH_RHTAB=>YDEAERSNK%RAMGROWTH_RHTAB,             &
          & RRHO_WAT=>YDEAERSNK%RHO_WAT, &
          & RRHO_SO4=>YDEAERSNK%RRHO_SO4, &
          & RRHO_SOA=>YDEAERSNK%RRHO_SOA, &
          & RRHO_NI=>YDEAERSNK%RRHO_NI, &
          & RRHO_AM=>YDEAERSNK%RRHO_AM, &
          & RRHO_ASH=>YDEAERSNK%RRHO_ASH, &
          & RRHO_BC=>YDEAERSNK%RRHO_BC, &
          & RRHO_OM=>YDEAERSNK%RRHO_OM, &
          & RMMD_SO4=>YDEAERSNK%RMMD_SO4, &
          & RMMD_SOA=>YDEAERSNK%RMMD_SOA, &
          & RMMD_NI=>YDEAERSNK%RMMD_NI, &
          & RMMD_AM=>YDEAERSNK%RMMD_AM, &
          & RMMD_ASH=>YDEAERSNK%RMMD_ASH, &
          & RMMD_BC=>YDEAERSNK%RMMD_BC, &
          & RMMD_OM=>YDEAERSNK%RMMD_OM, &
          & RMMD_DD=>YDEAERSNK%RMMD_DD, &
          & RMMD_SS=>YDEAERSNK%RMMD_SS, &
          & RAERDUST_REBOUND=>YDEAERATM%RAERDUST_REBOUND)

!KCHEM_DRYDEP: Switch between dry deposition options:
!KCHEM_DRYDEP=0: no online dry deposition velocities
!KCHEM_DRYDEP=1: Original, Wesely-based SUMO configuration (see routine depvel_wsl)
!KCHEM_DRYDEP=2: New GEOS-Chem-type mapping of land-use clases, but more standard configuration
!KCHEM_DRYDEP=3: Following GEOS-Chem type parameterization in ddr_surf_res_gc_V2, 
!KCHEM_DRYDEP=4: Same as '3', but use the stomatal resitance from CTESSEL, rather than GEOS-Chem type parameterization


!RVRSMIN(2)=110._JPRB    ! Short Grass
ZRSMIN_MOD(1)=100._JPRB    ! Crops, Mixed Farming
ZRSMIN_MOD(2)=100._JPRB    ! Short Grass
!RVRSMIN(3)=500._JPRB    ! Evergreen Needleleaf TrL, KLON, Pees
!RVRSMIN(4)=500._JPRB    ! Deciduous Needleleaf Trees
ZRSMIN_MOD(3)=250._JPRB    ! Evergreen Needleleaf Trees
ZRSMIN_MOD(4)=250._JPRB    ! Deciduous Needleleaf Trees
ZRSMIN_MOD(5)=175._JPRB    ! Deciduous Broadleaf Trees
ZRSMIN_MOD(6)=240._JPRB    ! Evergreen Broadleaf Trees
ZRSMIN_MOD(7)=100._JPRB    ! Tall Grass
ZRSMIN_MOD(8)=250._JPRB    ! Desert
ZRSMIN_MOD(9)=80._JPRB     ! Tundra
!RVRSMIN(10)=180._JPRB   ! Irrigated Crops
ZRSMIN_MOD(10)=100._JPRB   ! Irrigated Crops
ZRSMIN_MOD(11)=150._JPRB   ! Semidesert
ZRSMIN_MOD(12)=0.0_JPRB      ! Ice Caps and Glaciers
ZRSMIN_MOD(13)=240._JPRB   ! Bogs and Marshes
ZRSMIN_MOD(14)=0.0_JPRB      ! Inland Water
ZRSMIN_MOD(15)=0.0_JPRB      ! Ocean
ZRSMIN_MOD(16)=225._JPRB   ! Evergreen Shrubs
ZRSMIN_MOD(17)=225._JPRB   ! Deciduous Shrubs
ZRSMIN_MOD(18)=250._JPRB   ! Mixed Forest/woodland
ZRSMIN_MOD(19)=175._JPRB   ! Interrupted Forest
ZRSMIN_MOD(20)=150._JPRB   ! Water and Land Mixtures
ZRSMIN_MOD(0)=ZRSMIN_MOD(8)

PDEPVEL(KIDIA:KFDIA,1:KTRAC)=PDEPVELCLIM(KIDIA:KFDIA,1:KTRAC)

! IVEG_GC -> THE IVEG VEGETATION TYPES MAPPED TO 11 GEOS-CHEM DRY DEPOSITION TYPES
! REFERENCE TABLE FOR 20 VEG TYPES TO 11 GEOS-Chem DEP TYPES
!1)  LOW - Crops, Mixed Farming         ->    4: Agricultural land
!2)  LOW - Short Grass                  ->    5: Shurb/Grassland
!3)  HIGH -  Evergreen Needleleaf Trees ->    3: Coniferous forest
!4)  HIGH - Deciduous Needleleaf Trees  ->    3: Coniferous forest
!5)  HIGH - Deciduous Broadleaf Trees   ->    2: Deciduous forest
!6)  HIGH - Evergreen Broadleaf Trees   ->    6: Amazon forest
!7)  LOW - Tall Grass                   ->    5: Shrub/Grassland
!8)      - Desert                       ->    8: Desert
!9)  LOW - Tundra                       ->    7: Tundra
!10) LOW -  Irrigated Crops             ->    4: Agricultural
!11) LOW - Semidesert                   ->    5: Shrub/Grassland   - OR?: -
!11) LOW - Semidesert                   ->    8: Desert
!12)     - Ice Caps and Glaciers        ->    1: Snow/Ice
!13) LOW - Bogs and Marshes             ->    9: Wetland
!14)     - Inland Water                 ->   11: Water
!15)     - Ocean                        ->   11: Water
!16) LOW - Evergreen Shrubs             ->    5: Shrub/Grassland
!17) LOW - Deciduous Shrubs             ->    5: Shrub/Grassland
!18) HIGH - Mixed Forest/woodland       ->    2: Deciduous forest - This may be incorrect assumption
!19) HIGH - Interrupted Forest          ->    2: Deciduous forest - This may be incorrect assumption
!20) LOW - Water and Land Mixtures      ->    9: Wetland
IVEG_GC_TYPES = (/4,5,3,3,2,6,5,8,7,4,5,1,9,11,11,5,5,2,2,9/)

! Very simplistic mapping of IFS vegitation  and land types in aperge  classes - one vegeation type per grid box   
! IFS
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!            9 : Lakes

ZRAIN(:)=0.0_JPRB
ZDEW(:)=0.0_JPRB
IRH(:)=1
IMM=NMM(NINDAT)
ZSIGMA=2.0_JPRB

! Trace gas dry deposition
DO ITRCHEM = 1, KCHEM
    IF (YCHEM(ITRCHEM)%IGRIBDV <= 0 ) CYCLE

    ! find postion in KTRAC array
    ITR1=YGFL%NGHG+YGFL%NAERO+ITRCHEM
    DO JL=KIDIA, KFDIA

    ! Debug output
    ! ZLON=PGLAM(JL)/RPI * 180._JPRB
    ! ZLAT=PGLAT(JL)/RPI * 180._JPRB

    IDEBUG=0

    ! IF ( &
    !N-america     &  ((ABS(ZLON - (360. -92.5) ) < 0.4 ) .AND. (ABS(ZLAT - ( 34.) ) < 0.4)) .OR. &
    !Amazon        &  ((ABS(ZLON - (360. -62.5) ) < 0.4 ) .AND. (ABS(ZLAT - (-2.) ) < 0.4)))  &
    !Antarctic     &  ((ABS(ZLON - (   0.0) ) < 0.4 ) .AND. (ABS(ZLAT - (-74.) ) < 0.4)) .OR. &
    !Australian    &  ((ABS(ZLON - ( 130.0) ) < 0.4 ) .AND. (ABS(ZLAT - (-22.) ) < 0.4))) &
    !     &  THEN 
    !      IF (YCHEM(ITRCHEM)%CNAME == "O3") THEN
    !        IDEBUG=1
    !        WRITE(NULERR,"(a30,f8.2, f8.2)")'DEBUG O3 DDEP LON/LAT',ZLON,ZLAT
    !      ENDIF
    !    ENDIF

    ! Initialize deposition velocity
    PDEPVEL(JL, ITR1) = 0.0_JPRB 


    !*    Compute aerodynamic resistance
    ! aerodynamic resistance from reverse exchange coefficient
    ZRAERO = 1.0_JPRB/ PKCLEV(JL)

    !! Resistances limited to a maximum value
    ! IF (ZRAERO > PPRMAX) ZRAERO = PPRMAX
    
    ZAIRDEN =PRSF(JL)/(RD*PTS(JL))  ! Air density, kg/m3
    ZHEIGHT =  (PGEOM1(JL)-PGEOH(JL))/RG ! Height of middle (?) of lowest model level [m]
    ZXM = YCHEM(ITRCHEM)%RMOLMASS * 1E-3 ! molar mass, in units kg/mole

    !*          Compute laminar resistances based on temp,press & molecular weight
    CALL DDR_LAMINAR_RES_GC(PTS(JL), PRSF(JL),ZXM, PUSTAR(JL), YDRYDEP%RDIMO(ITRCHEM), ZWRB)
    ! Resistances limited to a maximum value
    ZVEG0(JL)=PCVL(JL) + PCVH(JL) ! Get vegetation cover

    ! IF (IDEBUG==1) WRITE(NULERR,'(a10,es12.5)')'ZSWRB=',ZWRB

    DO ID = 1, KTILES ! Loop over the number of tiles in each grid box (8)

        ZDEPTILE = PFRTI(JL,ID) ! Get the tile fraction land cover
        ! IF THE FRACTION OF TILE TYPE IS LOW (<0.01) THEN INGORE FOR EFFICIENCY
        IF (ZDEPTILE < 0.01) CYCLE

        ID_NOWET = ID

        !*    Compute aerodynamic resistance, now depending on tile. Currently depreciated.
        ! CALL DDR_AERO_RES_GC(PTS(JL), PZ0M(JL), ZHEIGHT, PUSTAR(JL), ZAIRDEN, PHFLUX(JL,ID), ZRAERO)
        ! IF (IDEBUG==1) WRITE(NULERR,'(a10,es12.5)')'ZRAERO=',ZRAERO

        !COULD PUT A SNOW/ICE FLAG HERE FOR SEPERATE CALCULATION
        ! IF (ID == (2,5,7)) THEN ...

        ILOW_VEG_NUM = KTVL(JL)
        IHIGH_VEG_NUM = KTVH(JL)
        ! INITILIASE ZLAI AS ZERO <- MIGHT NEED TO ADJUST THIS
        ZLAI = 0.0_JPRB

        ! CASE SELECT FOR TILE TYPE - MAPS TO GEOS-CHEM DRY DEP TYPE (IVEG_GC)
        SELECT CASE(ID)
            CASE(1,9) ! WATER / Lakes
                IVEG_GC = 11_JPIM
            CASE(2) ! ICE
                IVEG_GC = 1_JPIM
            CASE(3) ! WET SKIN
                ! IF WET SKIN MORE/LESS THAN 50% THEN ASSIGN TO LOW VEG/BARE SOIL
                ! THIS CURRENTLY DOESN'T MAKE SENSE WITH THE TILE FRACTION ITERATIONS
                IF ( ZVEG0(JL) < 0.5 ) THEN
                  IVEG_GC = 8_JPIM   ! re-assign to bare ground
                  ID_NOWET = 8_JPIM   ! re-assign to bare ground
                ELSE
                  IF(ILOW_VEG_NUM > 0) IVEG_GC = IVEG_GC_TYPES(ILOW_VEG_NUM) ! re-assign to low veg
                  ZLAI = PLAIL(JL)
                  ID_NOWET = 4_JPIM   ! re-assign to low veg
                  IF ( PCVL(JL) <=  PCVH(JL) ) THEN
                    IVEG_GC = IVEG_GC_TYPES(IHIGH_VEG_NUM) ! re-assign to high veg
                    ZLAI = PLAIH(JL)
                    ID_NOWET = 6_JPIM   ! re-assign to high veg
                  ENDIF
                ENDIF
            CASE(4) ! LOW VEG, NO SNOW
                IVEG_GC = IVEG_GC_TYPES(ILOW_VEG_NUM)
                ZLAI = PLAIL(JL)
            CASE(5) ! LOW VEG, SNOW ON
                IVEG_GC = 1_JPIM
            CASE(6) ! HIGH VEG, NO SNOW 
                IVEG_GC = IVEG_GC_TYPES(IHIGH_VEG_NUM)
                ZLAI = PLAIH(JL)
            CASE(7) ! HIGH VEG SNOW COVERED
                IVEG_GC = 1_JPIM
            CASE(8,10) ! BARE SOIL
                IVEG_GC = 8_JPIM
        END SELECT

        ! map IFS stomatal resistance according to dominant tile

        ZRSTO1=0.0_JPRB
        !  using RSMAX gives changes but on places with no vegetation
        !   ZRSTO1(JL)=RSMAX
        ITILE = ID
        IF (ID_NOWET == 4) THEN
          ZRSTO1=PCRL(JL)
        ELSEIF (ID_NOWET == 6) THEN
          ZRSTO1=PCRH(JL)
        ELSEIF (ID_NOWET == 7) THEN
          ZRSTO1=PCRHS(JL)
        ! ELSEIF (ITILE == 8) THEN
        !   ZRSTO1(JL)=PCRB(JL)
        ENDIF
        ZRSTO1 = MIN( ZRSTO1, RSMAX)

        !*       Compute surface and canopy resistances, either for (2) original SUMO type, 
        !*       or (3 & 4) Geos-Chem type parameterization.
        IF (KCHEM_DRYDEP == 2) THEN
          CALL  DDR_SURF_RES_GC (    PTS(JL), PFRSO(JL), &
                           &   ZRSTO1, ID, ID_NOWET, IVEG_GC,  ZLAI,  &
                           &   YCHEM(ITRCHEM)%CNAME, YDRYDEP%RCHEN(ITRCHEM), &
                           &   YDRYDEP%RCHENXP(ITRCHEM), YDRYDEP%RDIMO(ITRCHEM), &
                           &   YDRYDEP%RCF0(ITRCHEM), &
                           &   ZWRC)
        ELSEIF (KCHEM_DRYDEP == 3 .OR. KCHEM_DRYDEP == 4) THEN
          CALL  DDR_SURF_RES_GC_V2 (KCHEM_DRYDEP,  PTS(JL), PFRSO(JL), &
                           &   ZRSTO1, ID, ID_NOWET,  IVEG_GC,  ZLAI,  &
                           &   YCHEM(ITRCHEM)%CNAME, YDRYDEP%RCHEN(ITRCHEM), &
                           &   YDRYDEP%RCHENXP(ITRCHEM), YDRYDEP%RDIMO(ITRCHEM), &
                           &   YDRYDEP%RCF0(ITRCHEM), &
                           &   PMU0(JL),ZXM, PRSF(JL),IDEBUG, &
                           &   ZWRC)
        ENDIF

        !*         Compute Deposition Velocities
        !	        -----------------------------

        ZRES_TILE = (ZWRC+ZWRB+ZRAERO)  ! Sum of tile-specific laminar, aerodynamic and surface resistance 
        ! IF (IDEBUG==1) WRITE(NULERR,'(a20,2es12.5)')'ZRES_TILE  and tile frac',ZRES_TILE, ZDEPTILE

        ! add fraction to sum of deposition velocity.
        PDEPVEL(JL, ITR1) = PDEPVEL(JL, ITR1) + ZDEPTILE/(ZRES_TILE)

    ENDDO ! KTILES
    ! IF (IDEBUG==1) WRITE(NULERR,'(a10,2es12.5)')'-----------'
    ! IF (IDEBUG==1) WRITE(NULERR,'(a10,2es12.5)')'PDEPVEL ',PDEPVEL(JL,ITR1)
    ! IF (IDEBUG==1) WRITE(NULERR,'(a10,2es12.5)')'-----------'


  ENDDO !KIDIA, KFDIA loop
ENDDO ! Chemical species loop (ITRCHEM)

! aerosol dry deposition
IF (NDRYDEPVEL_DYN > 0) THEN
    
  ! First identify relative humidity class
  DO JL=KIDIA, KFDIA
    ZRHCL(JL)= EXP((17.625_JPRB*(PD2M(JL)-273.15_JPRB))/(PD2M(JL)-30.15_JPRB)) / &
  	     & EXP((17.625_JPRB*(PT2M(JL)-273.15_JPRB))/(PT2M(JL)-30.15_JPRB))
    ZTMP1=(17.625_JPRB*(PD2M(JL)-273.15_JPRB))/(PD2M(JL)-30.15_JPRB)
    ZTMP2=(17.625_JPRB*(PT2M(JL)-273.15_JPRB))/(PT2M(JL)-30.15_JPRB)
    DO JTAB=1,12
      IF (ZRHCL(JL)*100._JPRB > RRHTAB(JTAB)) THEN
  	IRH(JL)=JTAB
      ENDIF
    ENDDO
  ENDDO !KIDIA, KFDIA loop

  DO ITRAERO = 1, YGFL%NACTAERO
    IF (YAERO_DESC(ITRAERO)%RDDEPVSEA == 0._JPRB .AND. YAERO_DESC(ITRAERO)%RDDEPVLIC == 0._JPRB) CYCLE

    ! find postion in KTRAC array
    ITR1=YGFL%NGHG+ITRAERO

    ! aerosol size and density for each species
    IBIN=YAERO_DESC(ITRAERO)%NBIN
    ITYP=YAERO_DESC(ITRAERO)%NTYP
    SELECT CASE (ITYP)
       CASE(1)
    	 ZRHOP(:)=RSSDENS_RHTAB(:)
    	 ZWETD(:)=RMMD_SS(IBIN)*1.E-6_JPRB*RSSGROWTH_RHTAB(:)
       CASE(2)
    	 ZRHOP(:)=RRHO_DD(IBIN)
    	 ZWETD(:)=RMMD_DD(IBIN)*1.E-6_JPRB
       CASE(3)
    	 ZWETD(:)=RMMD_OM*1.E-6_JPRB*ROMGROWTH_RHTAB(:)
    	 ZVFRAC(:) = 1.0_JPRB / ROMGROWTH_RHTAB(:)**3
    	 ZRHOP(:) = RRHO_WAT*(1.0_JPRB-ZVFRAC(:)) + ZVFRAC(:)*RRHO_OM
       CASE(4)
    	 ZWETD(:)=RMMD_BC*1.E-6_JPRB
    	 ZRHOP(:) = RRHO_BC
       CASE(5)
    	 ZWETD(:)=RMMD_SO4*1.E-6_JPRB*RSO4GROWTH_RHTAB(:)
    	 ZVFRAC(:) = 1.0_JPRB / RSO4GROWTH_RHTAB(:)**3
    	 ZRHOP(:) = RRHO_WAT*(1.0_JPRB-ZVFRAC(:)) + ZVFRAC(:)*RRHO_SO4
       CASE(6)
    	 ZWETD(:)=RMMD_NI(IBIN)*1.E-6_JPRB*RNIGROWTH_RHTAB(:)
    	 ZVFRAC(:) = 1.0_JPRB / RNIGROWTH_RHTAB(:)**3
    	 ZRHOP(:) = RRHO_WAT*(1.0_JPRB-ZVFRAC(:)) + ZVFRAC(:)*RRHO_NI(IBIN)
       CASE(7)
    	 ZWETD(:)=RMMD_AM*1.E-6_JPRB*RAMGROWTH_RHTAB(:)
    	 ZVFRAC(:) = 1.0_JPRB / RAMGROWTH_RHTAB(:)**3
    	 ZRHOP(:) = RRHO_WAT*(1.0_JPRB-ZVFRAC(:)) + ZVFRAC(:)*RRHO_AM
       CASE(8)
    	 ZWETD(:)=RMMD_SOA*1.E-6_JPRB*RSOAGROWTH_RHTAB(:)
    	 ZVFRAC(:) = 1.0_JPRB / RSOAGROWTH_RHTAB(:)**3
    	 ZRHOP(:) = RRHO_WAT*(1.0_JPRB-ZVFRAC(:)) + ZVFRAC(:)*RRHO_SOA
       CASE(9)
    	 ZWETD(:)= RMMD_ASH*1.E-6_JPRB
    	 ZRHOP(:) = RRHO_ASH
       CASE(10)
    	 ZWETD(:)=RMMD_SO4*1.E-6_JPRB*RSO4GROWTH_RHTAB(:)
    	 ZVFRAC(:) = 1.0_JPRB / RSO4GROWTH_RHTAB(:)**3
    	 ZRHOP(:) = RRHO_WAT*(1.0_JPRB-ZVFRAC(:)) + ZVFRAC(:)*RRHO_SO4
    	CASE DEFAULT
    	  CALL ABOR1('DEPVEL_GC: unsupported aerosol type')
    END SELECT

    DO JL=KIDIA, KFDIA
      ! Initialize deposition velocity
      PDEPVEL(JL, ITR1) = 0.0_JPRB 

      !*    Compute aerodynamic resistance
      ! aerodynamic resistance from reverse exchange coefficient
      ZRAERO = 1.0_JPRB/ PKCLEV(JL)

      !! Resistances limited to a maximum value
      ! IF (ZRAERO > PPRMAX) ZRAERO = PPRMAX

      ZAIRDEN =PRSF(JL)/(RD*PTS(JL))  ! Air density, kg/m3
      ZVEG0(JL)=PCVL(JL) + PCVH(JL) ! Get vegetation cover

      DO ID = 1, KTILES ! Loop over the number of tiles in each grid box (8)

        ZDEPTILE = PFRTI(JL,ID) ! Get the tile fraction land cover
        ! IF THE FRACTION OF TILE TYPE IS LOW (<0.01) THEN INGORE FOR EFFICIENCY
        IF (ZDEPTILE < 0.01) CYCLE

        ID_NOWET=ID
        IVEGI=0
        ZLAI=0.0_JPRB
        SELECT CASE (ID)
           CASE(1,9,0) ! WATER / Lakes
             IVEGI=15
           CASE(2) ! ICE
             IVEGI=12
           CASE(3) ! WET SKIN
             IF ( ZVEG0(JL) < 0.5 ) THEN
               ! re-assign to bare ground
               IVEGI = 8_JPIM	  
               ID_NOWET = 8_JPIM  
             ELSE
               ! re-assign to low veg
               ID_NOWET = 4_JPIM   
               ZLAI = PLAIL(JL)    
               IVEGI=KTVL(JL)
               IF ( PCVL(JL) <=  PCVH(JL) ) THEN
             	 ! re-assign to high veg
             	 ID_NOWET = 6 
             	 ZLAI = PLAIH(JL)
             	 IVEGI=KTVH(JL)
               ENDIF
             ENDIF
           CASE (4) ! LOW VEG, NO SNOW
             ZLAI= PLAIL(JL)
             IVEGI=KTVL(JL)
           CASE (5) ! LOW VEG, SNOW ON
            ! assign ifs snow on low tile to arpege ice land cover
            IVEGI=8
           CASE (6,7) ! HIGH VEG NO SNOW (6) and SNOW COVERED (7)
             ZLAI= PLAIH(JL)
             IVEGI=KTVH(JL)
           CASE (8,10) ! BARE SOIL
             IVEGI=8
           CASE DEFAULT
             CALL ABOR1('DEPVEL_GC: unsupported land class')
        END SELECT

        !  Determines Wesely type 
        CALL DDR_ZH_SEASONONE(IMM, PGLAT(JL), ID_NOWET, IVEGI,ISEASON_WE, IVEG_ZH )
 
        SELECT CASE (NDRYDEPVEL_DYN) 
          CASE(1)
            ! compute deposition velocity following Zhang et al 2001
            CALL AER_DRYDEPVEL(ISEASON_WE, IVEG_ZH, ZRHOP(IRH(JL)),ZWETD(IRH(JL)),ZSIGMA,&
            & PZ0M(JL),PUSTARGUST(JL),PTS(JL),ZAIRDEN,RAERDUST_REBOUND,ZWRC)
          CASE(2)
            CALL AER_DRYDEPVELZH14(ISEASON_WE,IVEGI, IVEG_ZH, ZLAI, ZWETD(IRH(JL)), &
            & ZSIGMA,ZRHOP(IRH(JL)),&
            & PZ0M(JL),PUSTARGUST(JL),PTS(JL),RAERDUST_REBOUND,ZWRC)
          CASE(3)
            CALL AER_DRYDEPVELEM20(ISEASON_WE, IVEG_ZH,ZRHOP(IRH(JL)),ZWETD(IRH(JL)),ZSIGMA, &
            & PZ0M(JL),PUSTARGUST(JL),PTS(JL),ZAIRDEN,ZWRC)
        END SELECT
        ZRES_TILE = (ZWRC+ZRAERO)  ! Sum of tile-specific surface and aerodynamic resistance

        !write (*,*) "DEPVELFINALCHEM",ITRAERO,ITR1,ID,ZRES_TILE,ZRAERO,ZWRC,ZDEPTILE

        ! add fraction to sum of deposition velocity
        PDEPVEL(JL, ITR1) = PDEPVEL(JL, ITR1) + ZDEPTILE/(ZRES_TILE)
        !   SR 03/2018 reduce dry deposition over smoother snow surfaces
        IF (PSNS(JL) > 1.E-3_JPRB ) THEN
          PDEPVEL(JL, ITR1)=PDEPVEL(JL, ITR1)/2.5_JPRB
          PDEPVEL(JL, ITR1)=MIN(PDEPVEL(JL, ITR1),3.E-4_JPRB)
        ENDIF

      ENDDO ! KTILES

    ENDDO !KIDIA, KFDIA loop
  ENDDO ! Aerosol tracers loop (ITRAERO)
ENDIF ! (NDRYDEPVEL_DYN > 0)



END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('DEPVEL_GC',1,ZHOOK_HANDLE)

END SUBROUTINE  DEPVEL_GC 


