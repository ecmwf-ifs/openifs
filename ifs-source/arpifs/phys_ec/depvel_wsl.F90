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
                                  &  PTS, PT2M, PD2M, PRSF, PSNS, PFRSO, PRAQTI, PKCLEV, &  
                                  &  PUSTAR, PUSTARGUST,  PZ0M, PDEPVELCLIM, PDEPVEL )
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
!!    Add Urban tile    J. McNorton                        24 August 2022
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
USE YOMCST, ONLY: RD


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
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PT2M   ! 2m temperature
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PD2M   ! 2m dew point
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PRSF  ! Pressure at middle of bottom model level
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PSNS    ! SNOW


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
REAL(KIND=JPRB) , DIMENSION (KLON),INTENT(IN)     :: PUSTARGUST  ! friction velocity with gusts
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
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZRHCL   ! relative humidity
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZRAIN  ! soil wet or not with rain
REAL(KIND=JPRB) , DIMENSION (KLON)      :: ZDEW   ! soil wet or not with dew

! not KVTYPES seems to be 2 not 20 ???
REAL(KIND=JPRB), DIMENSION (0:20)            :: ZRSMIN_MOD ! Minimum stomatal resistances according to yom_veg rsmin  

REAL(KIND=JPRB)                          :: ZAIRDEN ! Air density [kg/m3]
REAL(KIND=JPRB)                          :: ZHEIGHT ! Height  [m]
REAL(KIND=JPRB)                          :: ZSIGMA, ZTMP1, ZTMP2
REAL(KIND=JPRB), DIMENSION (12)          :: ZVFRAC  ! Volume fraction
REAL(KIND=JPRB), DIMENSION (12)          :: ZRHOP, ZWETD  ! Aerosol density and wet diameter

INTEGER(KIND=JPIM) :: IRH(KLON)
INTEGER(KIND=JPIM) , DIMENSION (KLON)    :: IVEGA ! TYPE OF VEGETATION OF ARPEGE 1 ... 4 (WATER 1, ICE 2, LOW 3 AND HIGH VEG 4) 
INTEGER(KIND=JPIM) , DIMENSION (KLON)    :: IVEGI ! TYPE OF VEGETATION IN IFS  1 ... 20  
INTEGER(KIND=JPIM) , DIMENSION (KLON)    ::        ISEASON_WE  ! SEASON 
INTEGER(KIND=JPIM) , DIMENSION (KLON)    ::        IVEG_WE     ! VEGEATION TYP, WHICH NUMBER ??
INTEGER(KIND=JPIM) , DIMENSION (KLON)    ::        IVEG_ZH     ! VEGEATION TYP, WHICH NUMBER ??
INTEGER(KIND=JPIM) , DIMENSION (KLON)    ::  IDOMTILE, IDOMTILE_ORG   ! DOMINANT TILE NUMBER
INTEGER(KIND=JPIM)  :: ITR,  ITR1, ITRAERO, JL,  ID(1) , ITILE, IDDC, IMM , IDD, JTAB, IBIN, ITYP


REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE 

!
!-------------------------------------------------------------------------------

#include "fcttim.func.h"
#include "ddr_we_season.intfb.h"
#include "ddr_laminar_res.intfb.h"
#include "ddr_surf_res1.intfb.h"
#include "aer_drydepvel.intfb.h"
#include "aer_drydepvelzh14.intfb.h"
#include "aer_drydepvelem20.intfb.h"
#include "ddr_zh_season.intfb.h"


IF (LHOOK) CALL DR_HOOK('DEPVEL_WSL',0,ZHOOK_HANDLE)

ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL,YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM, &
          & YDEAERSNK=>YDMODEL%YRML_PHY_AER%YREAERSNK)
ASSOCIATE( NCHEM=>YGFL%NCHEM, YCHEM=>YGFL%YCHEM, NAERO=>YGFL%NAERO , &
          & NACTAERO=>YGFL%NACTAERO,&
          & YDRYDEP => YDMODEL%YRML_CHEM%YRDRYDEP, &
          & YAERO_DESC=>YDEAERATM%YAERO_DESC, &
          & NDRYDEPVEL_DYN => YDEAERSNK%NDRYDEPVEL_DYN, &
          & RRHTAB=>YDEAERSNK%RRHTAB,                               &
          & LAERDUST_NEWBIN=>YDEAERATM%LAERDUST_NEWBIN, &
          & RSSDENS_RHTAB=>YDEAERSNK%RSSDENS_RHTAB,                 &
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
          & RRHO_DD=>YDEAERSNK%RRHO_DD, &
          & RRHO_SS=>YDEAERSNK%RRHO_SS, &
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
ZSIGMA=2.0_JPRB

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
    CASE (8,10)
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
                       &   YCHEM(ITR)%CNAME, YDRYDEP%RCHEN(ITR), YDRYDEP%RCHENXP(ITR), &
                       &   YDRYDEP%RDIMO(ITR), YDRYDEP%RCF0(ITR), &
                       &   ZWRC)

!*         Compute Deposition Velocities
!	        -----------------------------

! find postion in KTRAC array (perhaps imporve on that by checking CLNAME) 
   ITR1=YGFL%NGHG+YGFL%NAERO+ITR 
   DO JL = KIDIA, KFDIA
         PDEPVEL(JL, ITR1) = 1.0_JPRB / (ZRAERO(JL) + ZWRB(JL) + ZWRC(JL))
   ENDDO
ENDDO 
! aerosol dry deposition
IF (NDRYDEPVEL_DYN > 0) THEN
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
         CALL ABOR1('DEPVEL_WSL: unsupported aerosol type')
    END SELECT

    CALL DDR_ZH_SEASON(KIDIA, KFDIA, KLON, IMM, IDD,  PGLAT, IDOMTILE, IVEGI, ISEASON_WE, IVEG_ZH)

    DO JL=KIDIA, KFDIA
      ! Initialize deposition velocity
      PDEPVEL(JL, ITR1) = 0.0_JPRB 
      ! Compute RH
      ZRHCL(JL)= EXP((17.625_JPRB*(PD2M(JL)-273.15_JPRB))/(PD2M(JL)-30.15_JPRB))/ &
         & EXP((17.625_JPRB*(PT2M(JL)-273.15_JPRB))/(PT2M(JL)-30.15_JPRB))
      ZTMP1=(17.625_JPRB*(PD2M(JL)-273.15_JPRB))/(PD2M(JL)-30.15_JPRB)
      ZTMP2=(17.625_JPRB*(PT2M(JL)-273.15_JPRB))/(PT2M(JL)-30.15_JPRB)
      DO JTAB=1,12
        IF (ZRHCL(JL)*100._JPRB > RRHTAB(JTAB)) THEN
          IRH(JL)=JTAB
        ENDIF
      ENDDO

      ZAIRDEN =PRSF(JL)/(RD*PTS(JL))  ! Air density, kg/m3
      !  Determines Wesely type 
 
      SELECT CASE (NDRYDEPVEL_DYN) 
        CASE(1)
         ! compute deposition velocity following Zhang et al 2001
          CALL AER_DRYDEPVEL(ISEASON_WE(JL), IVEG_ZH(JL), ZRHOP(IRH(JL)),ZWETD(IRH(JL)),ZSIGMA,&
          & PZ0M(JL),PUSTARGUST(JL),PTS(JL),ZAIRDEN,RAERDUST_REBOUND,ZWRC(JL))
        CASE(2)
          CALL AER_DRYDEPVELZH14(ISEASON_WE(JL),IVEGI(JL), IVEG_ZH(JL), ZLAI(JL), ZWETD(IRH(JL)), &
            & ZSIGMA,ZRHOP(IRH(JL)),&
            & PZ0M(JL),PUSTARGUST(JL),PTS(JL),RAERDUST_REBOUND,ZWRC(JL))
        CASE(3)
          CALL AER_DRYDEPVELEM20(ISEASON_WE(JL), IVEG_ZH(JL),ZRHOP(IRH(JL)),ZWETD(IRH(JL)),ZSIGMA, &
            & PZ0M(JL),PUSTARGUST(JL),PTS(JL),ZAIRDEN,ZWRC(JL))
      END SELECT
      PDEPVEL(JL, ITR1) = 1.0_JPRB / (ZRAERO(JL) + ZWRC(JL))

    !   SR 03/2018 reduce dry deposition over smoother snow surfaces
      IF (PSNS(JL) > 1.E-3_JPRB ) THEN
        PDEPVEL(JL, ITR1)=PDEPVEL(JL, ITR1)/2.5_JPRB
        PDEPVEL(JL, ITR1)=MIN(PDEPVEL(JL, ITR1),3.E-4_JPRB)
      ENDIF
    ENDDO !KIDIA, KFDIA loop
  ENDDO ! Aerosol tracers loop (ITRAERO)
ENDIF ! (NDRYDEPVEL_DYN > 0)



END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('DEPVEL_WSL',1,ZHOOK_HANDLE)

END SUBROUTINE  DEPVEL_WSL 

