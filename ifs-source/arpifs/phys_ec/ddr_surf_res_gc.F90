! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DDR_SURF_RES_GC ( PTS, PITM, PFRSO, &
                        &  PRSTO, KTILE, KTILE_NOWET, KVEG_GC, PLAI, &
                        &  CDNMS, PXCHEN, PXCHENXP, PXDIMO, PXCF0, &
                        &  PWRC)
!!    ---------
!!         *DDR_SURF_RES_GC* IS CALLED FROM *DEPVEL_GC*
!!
!!    PURPOSE
!!    -------
!!    The purpose of this routine is to cumpute the surface resistances for 
!!    the chemical species
!!    Computation is performed for a given date, forecast start time, 
!!    forecast valid time, and cells of a given domain. 
!!   
!!    METHOD
!!    ------
!!    See [Seinfeld et Pandis, 1998, "Atmospheric Chemistry and
!!    Physics", chap. 19, pp.958-996] based primarily on [Wesely, 1989]. 
!!    Note that in [Wesely,1989] only
!!    14 gaseous species are listed, against 17 in [Seinfeld et Pandis, 1998].
!!
!!    REFERENCE
!!    ---------
!!    See [Seinfeld et Pandis, 1998, "Atmospheric Chemistry and
!!    Physics", chap. 19, pp.958-996] based primarily on [Wesely, 1989]. 
!!    and CAMS_42 Deliverable report D42.3.1.1 (Dec. 2020)
!!
!!
!!
!! INPUTS:
!! -------
!! PTS        : Surface temperature
!! PITM       : Land/sea mask
!! PFRSO      : Solar radiation
!! PRSTO      : Stomatal resistances for water vapor(s.m-1)
!! KTILE      : IFS TILE (1... 9) 
!! KTILE_NOWET: IFS TILE (1... 9) no wet skin - reassigned to vegetation or bare ground
!! KVEG_GC     : Vegetation type
!! CDNMS       : Species name
!! PXCHEN      : Henry's specific constant
!! PXCHENXP    : Henry's specific constant exponent
!! PXDIMO      : Molecular diffusivities
!! PXCF0       : Reactivity
!!
!! OUTPUTS:
!! -------
!! ZWRC       : Surface resistances (s.m-1)
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
!!    M. Michou (MF)     : original    July 1999
!!    V. Marecal(MF)/J. Flemming adapted for IFS  in routine ddr_surf_res.F90
!!
!!    V. Huijnen / D. Finch, October 2020: introduction of loop over tiles, 
!!                                         and use of GEOS-Chem specific LUT variables.


USE PARKIND1 ,ONLY : JPIM, JPRB
USE DRYDEP_PAR_GC  ,ONLY : PPRMAX, ZEPS2,ZGSS_W, ZLUSW, & ! GEOS-Chem specific parameters
  &     IRLU, IRAC, IRGSS, IRGSO, IRCLS, IRCLO
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

  !! SUBROUTINE ARGUMENTS
  ! INTEGER (KIND=JPIM), INTENT(IN) :: KL, KLON

  REAL (KIND=JPRB), INTENT(IN)      :: PTS
  REAL (KIND=JPRB), INTENT(IN)      :: PITM
  REAL (KIND=JPRB), INTENT(IN)      :: PFRSO
  REAL (KIND=JPRB), INTENT(IN)      :: PRSTO
  INTEGER (KIND=JPIM), INTENT(IN)   :: KVEG_GC,KTILE,KTILE_NOWET
  REAL (KIND=JPRB), INTENT(IN)      :: PLAI
  CHARACTER (LEN=*), INTENT(IN)     :: CDNMS
  REAL (KIND=JPRB), INTENT(IN)      :: PXCHEN
  REAL (KIND=JPRB), INTENT(IN)      :: PXCHENXP
  REAL (KIND=JPRB), INTENT(IN)      :: PXDIMO
  REAL (KIND=JPRB), INTENT(IN)      :: PXCF0
  REAL (KIND=JPRB), INTENT(OUT)     :: PWRC

  !! LOCAL VARIABLES
  REAL (KIND=JPRB)                  :: ZWRCVEG, ZWRCNVEG
  INTEGER (KIND=JPIM)               :: ITILE  !LOOP CONTROLS
  REAL (KIND=JPRB)                  :: ZLUOW, ZGSSW
  REAL (KIND=JPRB)                  :: ZCLO, ZCLS, ZAC, ZGSO, ZGSS
  REAL (KIND=JPRB)                  :: ZDC, ZLU, ZTCOR, ZTS
  REAL (KIND=JPRB)                  :: ZHENRY, ZLUDRY, ZLUPROD
  REAL (KIND=JPRB)                  :: ZCLX, ZGSX, ZRES
  REAL(KIND=JPRB)                   :: ZTR,ZKA1,ZKW
  REAL(KIND=JPRB),PARAMETER         :: ZHPLUS =3.16227E-6_JPRB ! [H+]. For now assume soil pH=5.5
  LOGICAL                           :: LLWETVEG

  REAL(KIND=JPHOOK) ::  ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DDR_SURF_RES_GC',0,ZHOOK_HANDLE)
  ! In a first step, all cells are considered as water cells to optimise the
  ! execution time for the global domain

ZTCOR = 0.0_JPRB
ZWRCNVEG = 0.0_JPRB
ZWRCVEG = 0.0_JPRB

PWRC=-99999.9_JPRB

! surface independent treatment of HNO3, H2O2 
IF  (  ANY( (/ "HNO3","H2O2" /)   ==  CDNMS ) ) THEN 

  ! no temperature correction over sea / water 
    ZTCOR = 0.0_JPRB
    IF (PTS < 271.0_JPRB .AND. ANY( KTILE == (/ 2 , 4, 5, 6, 7, 8/) ))  THEN
      ZTS=MAX(PTS, 260._JPRB)
      ZTCOR = MIN(PPRMAX, 1000.0_JPRB * EXP(-ZTS + 269.0_JPRB) )
    ENDIF 
    PWRC = MAX(10.0_JPRB, ZTCOR)

ELSE 
! all other species 

! extra resistance for low temperatures - large for below 269 K ! in TM5 only for SO2 etc. 
    IF (PTS < 271.0_JPRB) THEN
!    IF (pts <  280.0_JPRB) THEN
      ZTS=MAX(PTS, 260._JPRB)
      ZTCOR = MIN(PPRMAX, 1000.0_JPRB * EXP(-ZTS + 269.0_JPRB) )
    ELSE
       ZTCOR = 0.0_JPRB
    ENDIF
!      1.2  Adjusts Henry's constant according to species and surface 
!             temperature. Later, could also vary with pH.
    ZHENRY = PXCHEN * EXP(PXCHENXP * ( 1.0_JPRB / PTS  - 1.0_JPRB / 298.0_JPRB   )   )
    IF  ( CDNMS == "NH3" ) THEN
      ZTR=(1.0_JPRB/PTS-1.0_JPRB/298.0_JPRB)
      ZKA1=1.7E-5_JPRB *EXP( -450._JPRB*ZTR )
      ZKW =1.0E-14_JPRB*EXP( -6718._JPRB*ZTR )
      ZHENRY = ZHENRY*(1._JPRB + ZKA1*ZHPLUS/ZKW)
    ENDIF
! IFS SURFACE TILE 
    ITILE=KTILE
    LLWETVEG=.FALSE.

! as wet skin also be vegetation (as given in ktile_nowet) we calculate resistances for wet vegetation,
! which (perhaps) differes from wet bare ground 
    IF (ITILE == 3 .AND. (KTILE_NOWET == 4 .OR. KTILE_NOWET == 6 )) THEN
      ITILE = KTILE_NOWET
      LLWETVEG=.TRUE. 
    ENDIF  

    SELECT CASE (ITILE)  

!        1.  Surface resistance over water
!          Water, Lake 
    CASE (1, 9)  
         !        1.1    Input resistances (see table 19.2)
        ! as they are all the same value through all the seasons (ZGSO = 0, ZGSS = 2000)
        ZGSO = MAX(ZEPS2,IRGSO(11))
        ZGSS = MAX(ZEPS2,IRGSS(11))

        SELECT CASE (TRIM(CDNMS))
            CASE ('O3', 'O3S', 'OX')
                PWRC = ZGSO !GC VERSION
            CASE ('SO2')
                PWRC = ZGSS !GC VERSION
            CASE DEFAULT
                PWRC = 1._JPRB / (1.E-5_JPRB * ZHENRY / ZGSS + PXCF0 / ZGSO)
                PWRC = MIN(PPRMAX, PWRC)
        END SELECT

! wet skin
    CASE (3)
     !   2.3.1  Input resistances (see table 19.2) | is it OK to also  ztcor for the bare ground ? 
     !   special cases 
        ZGSS = ZGSS_W ! UNSURE WHY THIS IS SET AT 50.0 - DF
        ZGSO = MAX(ZEPS2, IRGSO(9)) + ZTCOR ! ALSO USING WETLANDS (AS ABOVE),
        SELECT CASE (TRIM(CDNMS))
            CASE ('O3', 'O3S', 'OX')
                PWRC = ZGSO !GEOS-Chem vesion
            CASE ('SO2')
                PWRC = ZGSS
            CASE DEFAULT
                PWRC = 1._JPRB / (1.E-5_JPRB * ZHENRY / ZGSS + PXCF0 / ZGSO)
                PWRC = MIN(PPRMAX, PWRC)
       END SELECT

! Ice, snow , bare ground, snow on vegetation 
    CASE ( 2, 5, 7, 8 )        
      ! introduce special case for melting - see tm5 code 
  
      !       2.3.1  Input resistances (see table 19.2) | is it OK to also  ztcor for the bare ground ? 
        ZGSO = IRGSO(8) + ZTCOR
        ZGSS = IRGSS(8) + ZTCOR

        SELECT CASE (TRIM(CDNMS))
            CASE ('O3', 'O3S', 'OX')
                PWRC = ZGSO ! GC version
            CASE ('SO2')
                PWRC = ZGSS !GC version
            CASE DEFAULT
                PWRC = 1._JPRB / (1.E-5_JPRB * ZHENRY / ZGSS + PXCF0 / ZGSO)
                PWRC = MIN(PPRMAX, PWRC)
            END SELECT

!        2.4    Computation of the surface resistances for vegetated areas, high and low 
    CASE( 4,6 )
     !     2.4.1    Input resistances (see table 19.2 Seinfeld)
        ZLU  = MAX(ZEPS2, IRLU(KVEG_GC))
        IF (PLAI <= 0._JPRB) THEN
          ZLU = 1E6_JPRB
        ELSE
          ! CORRECT FOR SEASON BY USING LAI - cuticular resistances are per unit area of leaf
          ZLU = ZLU / PLAI
        ENDIF

        ZLUOW = 0.75_JPRB * ZLU + ZTCOR

        ! GEOS-CHEM VERSIONS
        ZCLO = MAX(ZEPS2, IRCLO(KVEG_GC)) + ZTCOR
        ZCLS = MAX(ZEPS2, IRCLS(KVEG_GC)) + ZTCOR
        ZAC = MAX(ZEPS2, IRAC(KVEG_GC))
        ZGSO = MAX(ZEPS2, IRGSO(KVEG_GC)) + ZTCOR
        ZGSS = MAX(ZEPS2, IRGSS(KVEG_GC)) + ZTCOR
        ZGSSW= ZGSS_W + ZTCOR

 !     2.4.2   Bulk canopy stomatal and mesophyll resistances eq. 19.25
        ZWRCVEG = PRSTO * PXDIMO + 1.0_JPRB / (100.0_JPRB * PXCF0 + 3.3E-4_JPRB * ZHENRY)
        ZWRCVEG = 1.0_JPRB / ZWRCVEG

!      2.4.3   Resistance of the outer surfaces in the upper canopy
        ZLUDRY = ZLU / (PXCF0 + 1.E-5_JPRB * ZHENRY) + ZTCOR    ! EQ. 19.26 ( add LAI according to  Gao and Wesely. 1995)
        IF ( LLWETVEG ) THEN
!       the wet skin vegation case has been overly simplified to bare ground wet skin -
!        if tile is wet skin no other info is given
            SELECT CASE (TRIM(CDNMS))
                CASE ('O3', 'O3S', 'OX')
                    ZWRCVEG = ZWRCVEG + 1.0_JPRB / ZLUOW
                CASE ('SO2')
                    ZWRCVEG = ZWRCVEG + 1.0_JPRB / ZLUSW
                CASE DEFAULT
                    ZLUPROD = 1.0_JPRB / (1.0_JPRB / (3.0_JPRB * ZLUDRY) + 1.E-7_JPRB * ZHENRY + PXCF0 / ZLUOW)
                    ZWRCVEG = ZWRCVEG + 1.0_JPRB / ZLUPROD
                END SELECT
        ELSE
            ZWRCVEG = ZWRCVEG + 1.0_JPRB / ZLUDRY
        ENDIF
!        2.4.4  Mixing forced by buoyant convection (ZDC)
!            Later, ZDC could take the slope of the local terrain into account
        ZDC = 100.0_JPRB * (1.0_JPRB + 1000._JPRB / (PFRSO + 10.0_JPRB)) ! FIRST PART OF 19.27 -units of G/pfrso
          !        2.4.5  Lower canopy resistance
        SELECT CASE(TRIM(CDNMS))
            CASE ('O3', 'O3S', 'OX')
                ZWRCVEG = ZWRCVEG + 1.0_JPRB / (ZCLO + ZDC)
            CASE ('SO2' )
                ZWRCVEG = ZWRCVEG + 1.0_JPRB / (ZCLS + ZDC)
            CASE DEFAULT
                ZCLX = 1.0_JPRB / (1.E-5_JPRB * ZHENRY / ZCLS + PXCF0 / ZCLO)  ! 19.28
                ZRES = ZDC + ZCLX
                ZWRCVEG = ZWRCVEG + 1.0_JPRB / ZRES
           END SELECT

!        2.4.6  Resistance due to canopies' height and density (ZAC, see above)
!        2.4.7  "Ground" resistance 
        IF ( LLWETVEG ) THEN
            ZGSS = ZGSSW
        ENDIF
        SELECT CASE(TRIM(CDNMS))
            CASE ('O3', 'O3S', 'OX')
                ZWRCVEG = ZWRCVEG + 1.0_JPRB / (ZGSO + ZAC)
            CASE ('SO2')
                ZWRCVEG = ZWRCVEG + 1.0_JPRB / (ZGSS + ZAC)
            CASE DEFAULT
                ZGSX = 1.0_JPRB / (1.E-5_JPRB * ZHENRY / ZGSS + PXCF0 / ZGSO) ! 19.29
                ZRES = ZAC + ZGSX
                ZWRCVEG = ZWRCVEG + 1.0_JPRB / ZRES
        END SELECT
        PWRC = 1.0_JPRB / ZWRCVEG
     END SELECT

ENDIF

! Resistances limited to a maximum value
IF (PWRC > PPRMAX) PWRC = PPRMAX

IF (LHOOK) CALL DR_HOOK('DDR_SURF_RES_GC',1,ZHOOK_HANDLE)
END SUBROUTINE DDR_SURF_RES_GC

