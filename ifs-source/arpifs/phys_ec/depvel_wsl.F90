! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DEPVEL_WSL (YDMODEL, KIDIA, KFDIA, KLON, KTRAC, KCHEM, KAERO,  KTILES, PFRTI , PGLAT, &
                                  &  PCVL , PCVH ,PLAIL, PLAIH,  KTVL, KTVH, &
                                  &  PCRB, PCRL, PCRLU , PCRH, PCRHS, &
                                  &  PTS, PFRSO, PRAQTI, PKCLEV, &  
                                  &  PUSTAR,  PZ0M, PDEPVELCLIM, PDEPVEL )
!!    PURPOSE
!!    -------
!!    The purpose of this routine is to cumpute the deposition velocities of 
!!    the chemical species of CAMS/MOCAGE 
!!    Computation is done for a given date, forecast start time, 
!! 
!!   
!!**  METHOD
!!    ------
!!    See [Seinfeld et Pandis, 1998, "Atmospheric Chemistry and
!!    Physics", chap. 19, pp.958-996] based primarily on [Wesely, 1989]. 
!!    Note that in [Wesely,1989] only
!!    14 gaseous species are listed, against 17 in [Seinfeld et Pandis, 1998].
!!    Note also the following differences with the Wesely proposed method:
!!    For IFS:   
!!    The aerodynamic resistance is taken from IFS surface scheme  
!!    The stomatal resistance is taken from IFS surface scheme  
!!    Mapping of Wesley classes to IFS surface classes    
!!
!!
!!    REFERENCE
!!    ---------
!!    [Seinfeld et Pandis, 1998, "Atmospheric Chemistry and
!!    Physics", chap. 19, pp.958-996]. 
!!
!
!!    AUTHOR
!!    ------
!!    M. Michou (original for MOCAGE) 
!!
!!
!!    MODIFICATIONS
!!    -------------
!!    Adapated to IFS - J. Flemming                        23 April 2018
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE DRYDEP_PAR, ONLY : PPRMAX, RSMAX
 
USE parkind1 ,ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMRIP0   , ONLY : NINDAT 
USE TYPE_MODEL , ONLY : MODEL


IMPLICIT NONE

TYPE(MODEL)        , INTENT(INOUT)       :: YDMODEL
INTEGER(KIND=JPIM) , INTENT(IN)          :: KLON, KTRAC
INTEGER(KIND=JPIM),INTENT(IN)    :: KCHEM(YDMODEL%YRML_GCONF%YGFL%NCHEM)
INTEGER(KIND=JPIM),INTENT(IN)    :: KAERO(YDMODEL%YRML_GCONF%YGFL%NAERO)

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

REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PUSTAR  ! friction velocity
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     ::  PZ0M   ! iRoughness length for momentum M
REAL(KIND=JPRB), DIMENSION (KLON, KTRAC), INTENT(IN)     :: PDEPVELCLIM

! Output 
REAL(KIND=JPRB), DIMENSION (KLON, KTRAC), INTENT(OUT)     :: PDEPVEL

! local 
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZRSMIN  ! minimum stomatal resistance
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZLAI  ! Lead Area Index
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZVEG0 ! vegetation cover 
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZITM   ! simple land sea mask   
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZRSTO1  ! stomatal resistances for water vapor from ifs/surf
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZRAERO ! aerodynamic resistances
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZWRB   ! laminar resistances 
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZWRC   ! surface resistances 
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZRAIN  ! soil wet or not with rain
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZDEW   ! soil wet or not with dew

! not KVTYPES seems to be 2 not 20 ???
REAL(KIND=JPRB), DIMENSION (0:20)            :: ZRSMIN_MOD ! Minimum stomatal resistances according to yom_veg rsmin  

INTEGER(KIND=JPIM) , DIMENSION (KLON)    :: IVEGA ! TYPE OF VEGETATION OF ARPEGE 1 ... 4 (WATER 1, ICE 2, LOW 3 AND HIGH VEG 4) 
INTEGER(KIND=JPIM) , DIMENSION (KLON)    :: IVEGI ! TYPE OF VEGETATION IN IFS  1 ... 20  
INTEGER(KIND=JPIM) , DIMENSION (KLON)    ::        ISEASON_WE  ! SEASON 
INTEGER(KIND=JPIM) , DIMENSION (KLON)    ::        IVEG_WE     ! VEGEATION TYP, WHICH NUMBER ??
INTEGER(KIND=JPIM) , DIMENSION (KLON)    ::  IDOMTILE, IDOMTILE_ORG   ! DOMINANT TILE NUMBER
INTEGER(KIND=JPIM)  :: ITR,  ITR1, JL,  ID(1) , ITILE, IDDC, IMM  
INTEGER(KIND=JPIM) :: ISHNO3_C, JAER
REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE 

!
!-------------------------------------------------------------------------------

#include "fcttim.func.h"
#include "ddr_we_season.intfb.h"
#include "ddr_laminar_res.intfb.h"
#include "ddr_surf_res1.intfb.h"

IF (LHOOK) CALL DR_HOOK('DEPVEL_WSL',0,ZHOOK_HANDLE)

ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL,YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM)
ASSOCIATE( NCHEM=>YGFL%NCHEM, YCHEM=>YGFL%YCHEM, NAERO=>YGFL%NAERO , &
 & NACTAERO=>YGFL%NACTAERO,YDRYDEP=>YDMODEL%YRML_CHEM%YRDRYDEP,&
 & YAERO_DESC=>YDEAERATM%YAERO_DESC)


IMM=NMM(NINDAT)

!RVRSMIN(2)=110._JPRB    ! Short Grass
ZRSMIN_MOD(1)=100._JPRB    ! Crops, Mixed Farming
ZRSMIN_MOD(2)=100._JPRB    ! Short Grass
!RVRSMIN(3)=500._JPRB    ! Evergreen Needleleaf Trees
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

! Very simplistic mapping of IFS vegitation  and land types in aperge  classes - one vegeation type per grid box   
! IFS
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL

!  arpege  jpvegnb_ar=4   ! Number of vegetation types: 
!                                      ! 1: sea; 2: ice; 3: low vegetation;
!                                      ! 4: forests for the vegetation

ZRAIN(:)=0.0_JPRB
ZDEW(:)=0.0_JPRB

DO JL=KIDIA, KFDIA 
   ZVEG0(JL)=PCVL(JL) + PCVH(JL)
! USE DOMINAT TYPE 
   ID=MAXLOC(PFRTI(JL,:))
   
   IDOMTILE_ORG(JL) = ID(1)
   IDOMTILE(JL) = ID(1)
   ZITM(JL)=1.0
! reassign  tile properties in case of wet skin (3) which could be hi lo veg or bare ground  
   IF ( IDOMTILE_ORG(JL) == 3 ) THEN 
        ZRAIN(JL)=1.0_JPRB
        ZDEW(JL)=1.0_JPRB
        IF ( ZVEG0(JL) < 0.5 ) THEN
             IDOMTILE(JL)=8   ! re-assign to bare ground 
        ELSE
           IDOMTILE(JL)= 4 ! re-assign to low veg   
           IF ( PCVL(JL) <=  PCVH(JL) )  IDOMTILE(JL)= 6 ! re-assign to high veg 
        ENDIF
   ENDIF   
   SELECT CASE (IDOMTILE(JL))  
    CASE(1,9,0)   
     IVEGA(JL)=1
     IVEGI(JL)=15 
     ZLAI(JL)=0.0_JPRB
     ZRSMIN(JL)=0.0_JPRB
     ZITM(JL)=0.0
    CASE(2)
     IVEGA(JL)=2
     IVEGI(JL)=12 
     ZLAI(JL)=0.0_JPRB
     ZRSMIN(JL)=0.0_JPRB
    CASE (6,7) 
     IVEGA(JL)=4
     ZLAI(JL)= PLAIH(JL) 
     IVEGI(JL)=KTVH(JL) 
     ZRSMIN(JL) = ZRSMIN_MOD(KTVH(JL))
    CASE (4)
     IVEGA(JL)=3
     IVEGI(JL)=KTVL(JL) 
     ZRSMIN(JL) = ZRSMIN_MOD(KTVL(JL))
     ZLAI(JL)= PLAIL(JL) 
    CASE (8)
     IVEGA(JL)=3
     IVEGI(JL)=8 
     ZRSMIN(JL)=ZRSMIN_MOD(8)
     ZLAI(JL)= 0.0_JPRB 
! assign ifs snow on low tile to arpege ice land cover
     CASE (5)
     IVEGA(JL)=2
     IVEGI(JL)=8 
     ZRSMIN(JL)=ZRSMIN_MOD(8)
     ZLAI(JL)= 0.0_JPRB
   END SELECT 
ENDDO


! Mapping of IFS classes into 11 regional Weseley classes 
CALL DDR_WE_SEASON(KIDIA, KFDIA, KLON, IMM, PGLAT, IDOMTILE, IVEGI, ISEASON_WE, IVEG_WE ) 

!*    Compute aerodynamic resistance
! aerodynamic resistance from reverse exchange coefficient  
DO JL = KIDIA, KFDIA 
! zraero(jl) = PRAQTI(JL,idomtile_org(jl)) 
! test new output of VEXCS 
  ZRAERO(JL) = 1.0_JPRB/ PKCLEV(JL) 
ENDDO 
!! Resistances limited to a maximum value
 WHERE (ZRAERO(KIDIA:KFDIA) > PPRMAX) ZRAERO(KIDIA:KFDIA) = PPRMAX   

! map IFS stomatal resistance according to dominant tile     
  DO JL = KIDIA, KFDIA  
    ZRSTO1(JL)=0.0_JPRB 
!  using RSMAX gives changes but on places with no vegetation
!   ZRSTO1(JL)=RSMAX 
    ITILE = IDOMTILE(JL) 
    IF (ITILE == 4) THEN
      ZRSTO1(JL)=PCRL(JL)
    ELSEIF (ITILE == 6) THEN
      ZRSTO1(JL)=PCRH(JL)
    ELSEIF (ITILE == 7) THEN
      ZRSTO1(JL)=PCRHS(JL)
!    ELSEIF (ITILE == 8) THEN
!      ZRSTO1(JL)=PCRB(JL)
    ENDIF
! VM    IF ( ZRSTO1(JL) >= 1.0E+6_JPRB )  ZRSTO1(JL) = 0.0_JPRB
    ZRSTO1(JL) = MIN( ZRSTO1(JL), RSMAX) 
  ENDDO 


IDDC=0
ISHNO3_C=0
DO ITR = 1, NCHEM  
  IF (YCHEM(ITR)%IGRIBDV <= 0 ) CYCLE    

    IDDC=IDDC+1
!*          Compute laminar resistances
  CALL DDR_LAMINAR_RES(KIDIA, KFDIA, KLON, PUSTAR, YDRYDEP%RDIMO(ITR), ZWRB)
! Resistances limited to a maximum value
  WHERE (ZWRB(KIDIA:KFDIA) > PPRMAX) ZWRB(KIDIA:KFDIA) = PPRMAX  

!*       Compute surface and canopy resistances
  CALL  DDR_SURF_RES1 ( KIDIA, KFDIA, KLON, &
                       &   PTS, ZITM, PFRSO, &
                       &   ZRSTO1, IDOMTILE_ORG, IDOMTILE, ISEASON_WE, IVEG_WE, &
                       &   YCHEM(ITR)%CNAME, YDRYDEP%RCHEN(ITR), YDRYDEP%RCHENXP(ITR), YDRYDEP%RDIMO(ITR), YDRYDEP%RCF0(ITR), &
                       &   ZWRC)

!*         Compute Deposition Velocities
!	        -----------------------------

! find postion in KTRAC array (perhaps imporve on that by checking CLNAME) 
   ITR1=YGFL%NGHG+YGFL%NAERO+ITR 
   DO JL = KIDIA, KFDIA
         PDEPVEL(JL, ITR1) = 1.0_JPRB / (ZRAERO(JL) + ZWRB(JL) + ZWRC(JL))
   ENDDO
    IF ((TRIM(YCHEM(ITR)%CNAME) == 'HNO3' )) THEN
      ISHNO3_C=ITR1
    ENDIF
ENDDO 
IF (ISHNO3_C /= 0) THEN
  !    Use HNO3 ddep velocities for NO3, if available
  DO JAER=1,NACTAERO
    ITR1=YGFL%NGHG+JAER
    IF (YAERO_DESC(JAER)%CNAME == 'Nitrate_1' ) THEN
       PDEPVEL(KIDIA:KFDIA,ITR1)=PDEPVEL(KIDIA:KFDIA,ISHNO3_C)
    ENDIF
  ENDDO
ENDIF

END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('DEPVEL_WSL',1,ZHOOK_HANDLE)

END SUBROUTINE  DEPVEL_WSL 

