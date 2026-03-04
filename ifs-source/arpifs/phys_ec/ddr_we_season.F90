! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DDR_WE_SEASON(kidia, kfdia, klon, kmonth, plat, &
                        &  ktile, kveg_ifs, kseason, kveg_we )
!!
!!    INTERFACE.
!!    ----------
!!         *DDR_WE_SEASON* IS CALLED FROM DEPVEL_WSL 
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
!! PLAT      : Latitude of grid point
!! KTILE     : Dominating tile 1:9
!! KVEG_IFS     : Vegetation type 1:20
!!
!! OUTPUTS:
!! -------
!! KSEASON     : Season
!! kveg_we     : Vegetation type
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
!!
!!-------------------------------------------------------------------------------

USE PARKIND1 ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

  !! Subroutine arguments
  INTEGER (KIND=JPIM), INTENT(IN) :: KMONTH 
  INTEGER (KIND=JPIM), INTENT(IN) :: KLON
  INTEGER (KIND=JPIM), INTENT(IN) :: KIDIA
  INTEGER (KIND=JPIM), INTENT(IN) :: KFDIA
  REAL (KIND=JPRB), DIMENSION(KLON), INTENT(IN)     :: PLAT 
  INTEGER (KIND=JPIM), DIMENSION(KLON), INTENT(IN)   :: KTILE
  INTEGER (KIND=JPIM), DIMENSION(KLON), INTENT(IN)   :: KVEG_IFS

  INTEGER (KIND=JPIM), DIMENSION(KLON), INTENT(OUT) :: KSEASON
  INTEGER (KIND=JPIM), DIMENSION(KLON), INTENT(OUT) :: KVEG_WE

  !! LOCAL VARIABLES
  REAL (KIND=JPRB), PARAMETER  :: ZD2R=57.2957795_JPRB
  REAL (KIND=JPRB), PARAMETER  :: ZLAT3=30.0_JPRB/ZD2R, ZLAT4=-30.0_JPRB/ZD2R


  REAL (KIND=JPRB)    :: ZLAT    ! LATITUDE OF GRID CELL CONSIDERED
  INTEGER (KIND=JPIM) :: JL, IVEG_IFS  ! LOOP CONTROLS
  REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DDR_WE_SEASON',0,ZHOOK_HANDLE)


  DO JL = KIDIA, KFDIA 

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
           IF ((ZLAT < ZLAT3) .AND. (ZLAT > ZLAT4)) THEN
              KSEASON(JL) = 1_JPIM
! season according to date !  - snow / ice cover always winter 4 - see next loop over tiles 
           ELSEIF (ZLAT > ZLAT3 ) THEN
              SELECT CASE(KMONTH)         ! NH DATES WHEN CHANGE OF WESELY SEASON:
                CASE (3,4,5)
                   KSEASON(JL) = 5_JPIM
                CASE (6,7,8)
                   KSEASON(JL) = 1_JPIM
                CASE (9,10,11)
                   KSEASON(JL) = 2_JPIM
                CASE (12,1,2)
!                   KSEASON(JL) = 4_JPIM
                   KSEASON(JL) = 3_JPIM
               END SELECT
           ELSEIF ( ZLAT < ZLAT4 ) THEN 
                SELECT CASE(KMONTH)         ! SH DATES WHEN CHANGE OF WESELY SEASON:
                CASE (3,4,5)
                   KSEASON(JL) = 2_JPIM
                CASE (6,7,8)
                   KSEASON(JL) = 3_JPIM
                CASE (9,10,11)
                   KSEASON(JL) = 5_JPIM
                CASE (12,1,2)
!                    KSEASON(JL) = 4_JPIM
                   KSEASON(JL) = 1_JPIM
              END SELECT
            ENDIF 

! map to dominat IFS tile 
! IFS
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL,  9: LAKES
! Seinfeld table 19.2
! 1 Urban land, 2: agricultural land , 3 range land 4: deciduous forest, 5 coniferous forest, 
! 6 mixed forest & wet lands,  7 water salt and fresh, 8 barren land / desert 
! 9 nonforested wetland 10 mixed agricultural and range land 11 rocky open areas with low shrubs

        iveg_ifs= kveg_ifs(jl)
        SELECT CASE(KTILE(JL))
           CASE (1,9,0)
              KVEG_WE(JL) = 7_JPIM        
           CASE (2)
              KVEG_WE(JL) = 7_JPIM
              KSEASON(JL) = 4_JPIM  
!           CASE (3)
!              KVEG_WE(JL) = 10_JPIM
!              PRAIN(JL)=1_JPIM
           CASE (4)
              KVEG_WE(JL) = 10_JPIM 
              IF (IVEG_IFS ==  1_JPIM)  KVEG_WE(JL) =  2_JPIM
              IF (IVEG_IFS ==  2_JPIM)  KVEG_WE(JL) =  3_JPIM
              IF (IVEG_IFS ==  7_JPIM)  KVEG_WE(JL) =  3_JPIM
              IF (IVEG_IFS ==  9_JPIM)  KVEG_WE(JL) =  11_JPIM
              IF (IVEG_IFS ==  10_JPIM)  KVEG_WE(JL) = 2_JPIM
              IF (IVEG_IFS ==  11_JPIM)  KVEG_WE(JL) =  8_JPIM
              IF (IVEG_IFS ==  13_JPIM)  KVEG_WE(JL) =  9_JPIM
              IF (IVEG_IFS ==  16_JPIM)  KVEG_WE(JL) =  3_JPIM
              IF (IVEG_IFS ==  17_JPIM)  KVEG_WE(JL) =  3_JPIM
              IF (IVEG_IFS ==  20_JPIM)  KVEG_WE(JL) =  9_JPIM
           CASE (5)
              KVEG_WE(JL) = 8_JPIM 
              KSEASON(JL) = 4_JPIM  
           CASE (6)
              KVEG_WE(JL) = 6_JPIM 
              IF (IVEG_IFS ==  3_JPIM)  KVEG_WE(JL) =  5_JPIM 
              IF (IVEG_IFS ==  4_JPIM)  KVEG_WE(JL) =  4_JPIM 
              IF (IVEG_IFS ==  5_JPIM)  KVEG_WE(JL) =  4_JPIM 
              IF (IVEG_IFS ==  6_JPIM)  KVEG_WE(JL) =  5_JPIM 
              IF (IVEG_IFS ==  18_JPIM)  KVEG_WE(JL) = 6_JPIM 
              IF (IVEG_IFS ==  19_JPIM)  KVEG_WE(JL) = 6_JPIM 
           CASE (7)
              KVEG_WE(JL) = 6_JPIM 
              IF (IVEG_IFS ==  3_JPIM)  KVEG_WE(JL) =  5_JPIM 
              IF (IVEG_IFS ==  4_JPIM)  KVEG_WE(JL) =  4_JPIM 
              IF (IVEG_IFS ==  5_JPIM)  KVEG_WE(JL) =  4_JPIM 
              IF (IVEG_IFS ==  6_JPIM)  KVEG_WE(JL) =  5_JPIM 
              IF (IVEG_IFS ==  18_JPIM)  KVEG_WE(JL) =  6_JPIM 
              IF (IVEG_IFS ==  19_JPIM)  KVEG_WE(JL) =  6_JPIM 
              KSEASON(JL) = 4_JPIM   
           CASE (8)
              KVEG_WE(JL) = 8_JPIM 
        END SELECT

  ENDDO

IF (LHOOK) CALL DR_HOOK('DDR_WE_SEASON',1,ZHOOK_HANDLE)

END SUBROUTINE DDR_WE_SEASON 
