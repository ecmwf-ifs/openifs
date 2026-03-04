! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DEPVEL_GC (YDMODEL, KIDIA, KFDIA, KLON, KTRAC, KCHEM,  KTILES, PFRTI , PGLAT, &
                                  &  PCVL , PCVH ,PLAIL, PLAIH,  KTVL, KTVH, &
                                  &  PCRB, PCRL, PCRLU , PCRH, PCRHS, &
                                  &  PMU0, PCFRAC, PRSF,PGEOM1,PGEOH,PHFLUX, &
                                  &  PTS, PFRSO, PRAQTI, PKCLEV, &  
                                  &  PUSTAR,  PZ0M, PDEPVELCLIM, PDEPVEL )
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
!!
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE DRYDEP_PAR, ONLY : RSMAX
 
USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD
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
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PMU0  ! COS of SZA
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PCFRAC ! Cloud fraction
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

REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PUSTAR  ! friction velocity
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     ::  PZ0M   ! iRoughness length for momentum M
REAL(KIND=JPRB), DIMENSION (KLON, KTRAC), INTENT(IN)     :: PDEPVELCLIM

! Output 
REAL(KIND=JPRB), DIMENSION (KLON, KTRAC), INTENT(OUT)     :: PDEPVEL

! local 
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZVEG0 ! vegetation cover 
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZITM   ! simple land sea mask   
REAL(KIND=JPRB)                         :: ZRSTO1  ! stomatal resistances for water vapor from ifs/surf
REAL(KIND=JPRB)                         :: ZRAERO ! aerodynamic resistances
REAL(KIND=JPRB)                         :: ZWRB   ! laminar resistances
REAL(KIND=JPRB)                         :: ZWRC   ! surface resistances for GEOS-Chem veg types
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZRAIN  ! soil wet or not with rain
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZDEW   ! soil wet or not with dew
! REAL(KIND=JPRB)                         :: ZLON, ZLAT

! not KVTYPES seems to be 2 not 20 ???
REAL(KIND=JPRB), DIMENSION (0:20)        :: ZRSMIN_MOD ! Minimum stomatal resistances according to yom_veg rsmin
REAL(KIND=JPRB)                          :: ZRES_TILE ! Temporary resistances
REAL(KIND=JPRB)                          :: ZDEPTILE ! Proportion of veg on tile
REAL(KIND=JPRB)                          :: ZLAI ! Leaf area index as input for surface resistance mod
REAL(KIND=JPRB)                          :: ZXM ! Molecular weight
REAL(KIND=JPRB)                          :: ZAIRDEN ! Air density [kg/m3]
REAL(KIND=JPRB)                          :: ZHEIGHT ! Height  [m]


INTEGER(KIND=JPIM) , DIMENSION (20)      :: IVEG_GC_TYPES! Reference table to map iveg to GEOS-Chem types
INTEGER(KIND=JPIM)                       :: ILOW_VEG_NUM, IHIGH_VEG_NUM ! Reference number for veg types
INTEGER(KIND=JPIM)                       :: IVEG_GC   ! type of vegetation in reference to GEOS-Chem
INTEGER(KIND=JPIM)                       :: ID_NOWET  ! Reconstructed dry tile 
INTEGER(KIND=JPIM)  :: ITR,  ITR1, JL,  ID , ITILE , IDEBUG
REAL(KIND=JPHOOK) ::  ZHOOK_HANDLE 

!
!-------------------------------------------------------------------------------

#include "fcttim.func.h"
! #include "ddr_aero_res_gc.intfb.h"
#include "ddr_laminar_res_gc.intfb.h"
#include "ddr_surf_res_gc.intfb.h"
#include "ddr_surf_res_gc_v2.intfb.h"

IF (LHOOK) CALL DR_HOOK('DEPVEL_GC',0,ZHOOK_HANDLE)

ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL, YCHEM=>YDMODEL%YRML_GCONF%YGFL%YCHEM, &
   & KCHEM_DRYDEP=>YDMODEL%YRML_CHEM%YRCHEM%KCHEM_DRYDEP,YDRYDEP=>YDMODEL%YRML_CHEM%YRDRYDEP)

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
DO ITR = 1, KCHEM
    IF (YCHEM(ITR)%IGRIBDV <= 0 ) CYCLE
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
    !      IF (YCHEM(ITR)%CNAME == "O3") THEN
    !        IDEBUG=1
    !        WRITE(NULERR,"(a30,f8.2, f8.2)")'DEBUG O3 DDEP LON/LAT',ZLON,ZLAT
    !      ENDIF
    !    ENDIF

    ! find postion in KTRAC array (perhaps improve on that by checking CLNAME)
    ITR1=YGFL%NGHG+YGFL%NAERO+ITR
    PDEPVEL(JL, ITR1) = 0.0_JPRB 


    !*    Compute aerodynamic resistance
    ! aerodynamic resistance from reverse exchange coefficient

    ZRAERO = 1.0_JPRB/ PKCLEV(JL)
    !! Resistances limited to a maximum value
    ! IF (ZRAERO > PPRMAX) ZRAERO = PPRMAX
    
    ZAIRDEN =PRSF(JL)/(RD*PTS(JL))  ! Air density, kg/m3
    ZHEIGHT =  (PGEOM1(JL)-PGEOH(JL))/RG ! Height of middle (?) of lowest model level [m]
    ZXM = YCHEM(ITR)%RMOLMASS * 1E-3 ! molar mass, in units kg/mole

    !*          Compute laminar resistances based on temp,press & molecular weight
    CALL DDR_LAMINAR_RES_GC(PTS(JL), PRSF(JL),ZXM, PUSTAR(JL), YDRYDEP%RDIMO(ITR), ZWRB)
    ! Resistances limited to a maximum value
    ZVEG0(JL)=PCVL(JL) + PCVH(JL) ! Get vegetation cover
    ZITM(JL)=1.0 ! Simple land/sea mask

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
            CASE(8) ! BARE SOIL
                IVEG_GC = 8_JPIM
        END SELECT

       !    IF (IDEBUG==1) WRITE(NULERR,'(a10,5i5)')'ID, IDVEG_GC=',ID,IVEG_GC
       !    IF (IDEBUG==1) WRITE(NULERR,'(a10,es12.5)')'LAI=',ZLAI


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
    !    ELSEIF (ITILE == 8) THEN
    !      ZRSTO1(JL)=PCRB(JL)
        ENDIF
        ZRSTO1 = MIN( ZRSTO1, RSMAX)

      !*       Compute surface and canopy resistances, either for (2) original SUMO type, 
      !*       or (3 & 4) Geos-Chem type parameterization.
      IF (KCHEM_DRYDEP == 2) THEN
      CALL  DDR_SURF_RES_GC (    PTS(JL), ZITM(JL), PFRSO(JL), &
                           &   ZRSTO1, ID, ID_NOWET, IVEG_GC,  ZLAI,  &
                           &   YCHEM(ITR)%CNAME, YDRYDEP%RCHEN(ITR), YDRYDEP%RCHENXP(ITR), YDRYDEP%RDIMO(ITR), YDRYDEP%RCF0(ITR), &
                           &   ZWRC)
      ELSEIF (KCHEM_DRYDEP == 3 .OR. KCHEM_DRYDEP == 4) THEN
      CALL  DDR_SURF_RES_GC_V2 (KCHEM_DRYDEP,  PTS(JL), ZITM(JL), PFRSO(JL), &
                           &   ZRSTO1, ID, ID_NOWET,  IVEG_GC,  ZLAI,  &
                           &   YCHEM(ITR)%CNAME, YDRYDEP%RCHEN(ITR), YDRYDEP%RCHENXP(ITR), YDRYDEP%RDIMO(ITR), YDRYDEP%RCF0(ITR), &
                           &   PMU0(JL),PCFRAC(JL),ZXM, PRSF(JL),IDEBUG, &
                           &   ZWRC)
      ENDIF


      ! IF (IDEBUG==1) WRITE(NULERR,'(a10,1es12.5)')'ZRSTO1=',ZRSTO1
      ! IF (IDEBUG==1) WRITE(NULERR,'(a10,2es12.5)')'T,rad=',PTS(JL),PFRSO(JL)
      ! IF (IDEBUG==1) WRITE(NULERR,'(a10,es12.5)')'ZWRC',ZWRC


    !*         Compute Deposition Velocities
    !	        -----------------------------

    ZRES_TILE = (ZWRC+ZWRB+ZRAERO)  ! Sum of (laminar resistance and surface resistance) * tile fraction
    ! IF (IDEBUG==1) WRITE(NULERR,'(a20,2es12.5)')'ZRES_TILE  and tile frac',ZRES_TILE, ZDEPTILE


    ! add fraction to sum of deposition velocity.
    PDEPVEL(JL, ITR1) = PDEPVEL(JL, ITR1) + ZDEPTILE/(ZRES_TILE)

    ENDDO ! KTILES
    ! IF (IDEBUG==1) WRITE(NULERR,'(a10,2es12.5)')'-----------'
    ! IF (IDEBUG==1) WRITE(NULERR,'(a10,2es12.5)')'PDEPVEL ',PDEPVEL(JL,ITR1)
    ! IF (IDEBUG==1) WRITE(NULERR,'(a10,2es12.5)')'-----------'


    ENDDO !KIDIA, KFDIA loop
ENDDO ! Chemical species loop (KCHEM)

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('DEPVEL_GC',1,ZHOOK_HANDLE)

END SUBROUTINE  DEPVEL_GC 


