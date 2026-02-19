! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DDR_ZH_SEASON(KIDIA, KFDIA, KLON, KMONTH, KDAY,  PLAT,                &
&                        KTILE, KVEG_IFS, KSEASON, KVEG_ZH )
!!                                 snow on the ground. 
!!                                 High and tropical latitudes treated separately.
!!
!!    INTERFACE.
!!    ----------
!!         *DDR_ZH_SEASON* IS CALLED FROM 
!!
!!    PURPOSE
!!    -------
!!   
!!    METHOD
!!    ------
!!
!! INPUTS:
!! -------
!! KMONTH    : Month
!! KDAY      : Day
!! NLON      : Number of grid points in longitude
!! NLAT      : Number of grid points in latitude
!! PLAT      : Latitude of grid point
!! KTILE     : Dominating tile 1:9
!! KVEG_IFS     : Vegetation type 1:20
!!
!! OUTPUTS:
!! -------
!! KSEASON     : Season
!! kveg_ZH     : Vegetation type
!!
!! LOCAL:
!! -------
!!
!!    AUTHOR
!!    ------
!!    M. Michou 
!!
!!
!!    MODIFICATIONS
!!    -------------
!!    Original    July 1999
!!    adapted for IFS land classes , J. Flemming 1.6.2015
!!    adapted for Zhang et al land classes , S.Rémy 17/1/2017
!!
!!-------------------------------------------------------------------------------

USE PARKIND1 ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

  !! Subroutine arguments
  INTEGER (KIND=JPIM), INTENT(IN) :: kmonth 
  INTEGER (KIND=JPIM), INTENT(IN) :: kday
  INTEGER (KIND=JPIM), INTENT(IN) :: klon
  INTEGER (KIND=JPIM), INTENT(IN) :: kidia
  INTEGER (KIND=JPIM), INTENT(IN) :: kfdia
  REAL (KIND=JPRB), DIMENSION(klon), INTENT(IN)     :: plat 
  INTEGER (KIND=JPIM), DIMENSION(klon), INTENT(IN)   :: ktile
  INTEGER (KIND=JPIM), DIMENSION(klon), INTENT(IN)   :: kveg_ifs

  INTEGER (KIND=JPIM), DIMENSION(klon), INTENT(OUT) :: kseason
  INTEGER (KIND=JPIM), DIMENSION(klon), INTENT(OUT) :: kveg_ZH

  !! Local variables
  REAL (KIND=JPRB), PARAMETER  :: zd2r=57.2957795_JPRB
  !REAL (KIND=JPRB), PARAMETER  :: zlat1=65.0_JPRB/zd2r,  zlat2=-60.0_JPRB/zd2r
  REAL (KIND=JPRB), PARAMETER  :: zlat3=30.0_JPRB/zd2r, zlat4=-30.0_JPRB

  REAL (KIND=JPRB)    :: zlat    ! latitude of grid cell considered
  INTEGER (KIND=JPIM) :: jl, iveg_ifs  ! loop controls
  REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DDR_ZH_SEASON',0,ZHOOK_HANDLE)


  DO jl = kidia, kfdia 


        !        1   Determines Wesely's seasonal type according to calendar date,
        !            snow on the ground. High and tropical latitudes treated separately.
        !            -------------------------------------------------------------------
        !              Type #1 : Midsummer with lush vegetation 
        !              Type #2 : Autumn with unharvested cropland
        !              Type #3 : Late automn after frost, no snow
        !              Type #4 : Winter, snow on ground and subfreezing
        !              Type #5 : Transitional spring with partially green short annuals

        zlat = plat(jl)

! assign season 
! tropics always lush summer
           IF ((zlat < zlat3) .AND. (zlat > zlat4)) THEN
              kseason(jl) = 1_JPIM
! season according to date !  - snow / ice cover always winter 4 - see next loop over tiles 
           ELSEIF (zlat >= zlat3 ) THEN
              SELECT CASE(kmonth)         ! NH Dates when change of Wesely season:
                CASE (3,4,5)
                   kseason(jl) = 5_JPIM
                CASE (6,7,8)
                   kseason(jl) = 1_JPIM
                CASE (9,10,11)
                   kseason(jl) = 2_JPIM
                CASE (12,1,2)
!                   kseason(jl) = 4_JPIM
                   kseason(jl) = 3_JPIM
               END SELECT
           ELSEIF ( zlat <= zlat4 ) THEN 
                SELECT CASE(kmonth)         ! SH Dates when change of Wesely season:
                CASE (3,4,5)
                   kseason(jl) = 2_JPIM
                CASE (6,7,8)
                   kseason(jl) = 3_JPIM
                CASE (9,10,11)
                   kseason(jl) = 5_JPIM
                CASE (12,1,2)
!                    kseason(jl) = 4_JPIM
                   kseason(jl) = 1_JPIM
              END SELECT
            ENDIF 

! map to dominat IFS tile 
! IFS
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL,  9: LAKES
! Seinfeld table 19.2
! 1 Urban land, 2: agricultural land , 3 range land 4: deciduous forest, 5 coniferous forest, 6 mixed forest & wet lands 
! 7 water salt and fresh, 8 barren land / desert, 9 nonforested wetland, 10 mixed agricultural and range land
! 11 rocky open areas with low shrubs

! ifs vegetation 
!  
!1)=100._JPRB   L ! Crops, Mixed Farming             2: agricultural land 
!2)=100._JPRB   L ! Short Grass                      3: range land 
!3)=250._JPRB   H ! Evergreen Needleleaf Trees       5: coniferous forest
!4)=250._JPRB   H ! Deciduous Needleleaf Trees       5: coniferous forest
!5)=175._JPRB   H ! Deciduous Broadleaf Trees        4: deciduous forest
!6)=240._JPRB   H ! Evergreen Broadleaf Trees        4: deciduous forest  
!7)=100._JPRB   L ! Tall Grass                       3: range land 
!8)=250._JPRB    ! Desert                            
!9)=80._JPRB    L ! Tundra                          11: rocky open areas with low shrubs
!10)=100._oPRB  L ! Irrigated Crops                  2: agricultural   
!11)=150._JPRB  L ! Semidesert                       8: barren land / desert  
!12)=0.0_JPRB      ! Ice Caps and Glaciers
!13)=240._JPRB  L ! Bogs and Marshes                 9: nonforested wetland 
!14)=0.0_JPRB      ! Inland Water
!15)=0.0_JPRB      ! Ocean
!16)=225._JPRB  L ! Evergreen Shrubs                 3: range land   
!17)=225._JPRB  L ! Deciduous Shrubs                 3: range land
!18)=250._JPRB  H ! Mixed Forest/woodland            6: mixed forest & wet lands 
!19)=175._JPRB  H ! Interrupted Forest               6: mixed forest & wet lands
!20)=150._JPRB  L ! Water and Land Mixtures          9: nonforested wetland


! note wet skin tile has been rassigned to original tile - but rain/dew was set to 1
        iveg_ifs= kveg_ifs(jl)
        SELECT CASE(ktile(jl))
           CASE (1,0)
              kveg_ZH(jl) = 14_JPIM        
           CASE (9)
              kveg_ZH(jl) = 13_JPIM        
           CASE (2)
              kveg_ZH(jl) = 12_JPIM
           CASE (3)
              kveg_ZH(jl) = 11_JPIM
           CASE (4)
              kveg_ZH(jl) = 6_JPIM 
              if (iveg_ifs ==  1_JPIM)  kveg_ZH(jl) =  7_JPIM
              if (iveg_ifs ==  2_JPIM)  kveg_ZH(jl) =  6_JPIM
              if (iveg_ifs ==  7_JPIM)  kveg_ZH(jl) =  6_JPIM
              if (iveg_ifs ==  9_JPIM)  kveg_ZH(jl) =  9_JPIM
              if (iveg_ifs ==  10_JPIM)  kveg_ZH(jl) = 7_JPIM
              if (iveg_ifs ==  11_JPIM)  kveg_ZH(jl) =  8_JPIM
              if (iveg_ifs ==  13_JPIM)  kveg_ZH(jl) =  11_JPIM
              if (iveg_ifs ==  16_JPIM)  kveg_ZH(jl) =  10_JPIM
              if (iveg_ifs ==  17_JPIM)  kveg_ZH(jl) =  10_JPIM
              if (iveg_ifs ==  20_JPIM)  kveg_ZH(jl) =  11_JPIM
           CASE (5)
              kveg_ZH(jl) = 9_JPIM 
              kseason(jl) = 4_JPIM  
           CASE (6)
              kveg_ZH(jl) = 5_JPIM 
              if (iveg_ifs ==  3_JPIM)  kveg_ZH(jl) =  1_JPIM 
              if (iveg_ifs ==  4_JPIM)  kveg_ZH(jl) =  3_JPIM 
              if (iveg_ifs ==  5_JPIM)  kveg_ZH(jl) =  4_JPIM 
              if (iveg_ifs ==  6_JPIM)  kveg_ZH(jl) =  2_JPIM 
              if (iveg_ifs ==  18_JPIM)  kveg_ZH(jl) = 5_JPIM 
              if (iveg_ifs ==  19_JPIM)  kveg_ZH(jl) = 10_JPIM 
           CASE (7)
              kveg_ZH(jl) = 5_JPIM 
              if (iveg_ifs ==  3_JPIM)  kveg_ZH(jl) =  1_JPIM 
              if (iveg_ifs ==  4_JPIM)  kveg_ZH(jl) =  3_JPIM 
              if (iveg_ifs ==  5_JPIM)  kveg_ZH(jl) =  4_JPIM 
              if (iveg_ifs ==  6_JPIM)  kveg_ZH(jl) =  2_JPIM 
              if (iveg_ifs ==  18_JPIM)  kveg_ZH(jl) = 5_JPIM 
              if (iveg_ifs ==  19_JPIM)  kveg_ZH(jl) = 10_JPIM 
              kseason(jl) = 4_JPIM   
           CASE (8)
              kveg_ZH(jl) = 8_JPIM 
        END SELECT


  ENDDO

IF (LHOOK) CALL DR_HOOK('DDR_ZH_SEASON',1,ZHOOK_HANDLE)

END SUBROUTINE DDR_ZH_SEASON 
