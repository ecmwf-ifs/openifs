! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE RADVIS &
  &( KIDIA, KFDIA, KLON, KLEV, &
  &  PRSF1, PT, PA, PL, PI, PR, PS, &
  &  PAEREXT, &
  &  PVISALL, PVISAER, PVISCLD )
  
!**** *RADVIS* - ROUTINE COMPUTING THE VISIBILITY

!      J.-J. MORCRETTE , R M Forbes,  ECMWF


!**   INTERFACE.
!     ----------
!          *RADVIS* IS CALLED FROM *RADVIS_LAYER*.

! INPUTS:
! -------

! OUTPUTS:
! --------

!     MODIFICATIONS.
!     -------------
!     Original: JJMorcrette, 20101125  
!     Nov 2022: R Forbes  Updated with improved options as in claervis.F90

!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM , JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RD
USE YOESRTCOP, ONLY : RSASWA, RSASWB, RSFUA0, RSFUA1

!-------------------------------------------------------------------------------

IMPLICIT NONE

! Input 
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV
REAL(KIND=JPRB),INTENT(IN)    :: PRSF1(KLON,KLEV) ! Pressure on full levels (Pa)
REAL(KIND=JPRB),INTENT(IN)    :: PT(KLON,KLEV)    ! Temperature (K)
REAL(KIND=JPRB),INTENT(IN)    :: PA(KLON,KLEV)    ! Cloud fraction (0-1)
REAL(KIND=JPRB),INTENT(IN)    :: PL(KLON,KLEV)    ! Cloud liquid (kg kg-1)
REAL(KIND=JPRB),INTENT(IN)    :: PI(KLON,KLEV)    ! Cloud ice (kg kg-1) 
REAL(KIND=JPRB),INTENT(IN)    :: PR(KLON,KLEV)    ! Rain (kg kg-1)
REAL(KIND=JPRB),INTENT(IN)    :: PS(KLON,KLEV)    ! Snow (kg kg-1)
REAL(KIND=JPRB),INTENT(IN)    :: PAEREXT(KLON)    ! Aerosol extinction coeff (m-1)

! Output
REAL(KIND=JPRB),INTENT(OUT)   :: PVISALL(KLON) ! Rayleigh + cloud + aerosol (m)
REAL(KIND=JPRB),INTENT(OUT), OPTIONAL :: PVISAER(KLON) ! Rayleigh + aerosol (m)
REAL(KIND=JPRB),INTENT(OUT), OPTIONAL :: PVISCLD(KLON) ! Rayleigh + cloud (m)

!-- Local variables --
INTEGER(KIND=JPIM) :: JL, JSW

REAL(KIND=JPRB)    :: ZAIRDENSITY        ! Air density (kg m-3)
REAL(KIND=JPRB)    :: ZLOGNETA           ! Natural log of liminal visual contrast 
REAL(KIND=JPRB)    :: ZRE_LIQ,  ZDE_ICE  ! Cloud liquid and ice hydrometeor sizes
REAL(KIND=JPRB)    :: ZRE_RAIN, ZDE_SNOW ! Rain and snow hydrometeor sizes
REAL(KIND=JPRB)    :: ZND_LIQ            ! Cloud liquid droplet size for Gultepe (2006)
REAL(KIND=JPRB)    :: ZEXTLIQ,  ZEXTICE  ! Extinction coeffs for cloud liquid and ice
REAL(KIND=JPRB)    :: ZEXTRAIN, ZEXTSNOW ! Extinction coeffs for rain and snow
REAL(KIND=JPRB)    :: ZQLWC,    ZQIWC    ! Water contents for cloud liquid and ice
REAL(KIND=JPRB)    :: ZQRWC,    ZQSWC    ! Water contents for rain and snow
REAL(KIND=JPRB)    :: ZRAYEXT            ! Rayleigh extinction coefficient (m-1)
REAL(KIND=JPRB)    :: ZHYDEXT(KLON)      ! Hydrometeor extinction coefficient (m-1)

REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

! Visibility algorithm options: 
! 'IFSRAD'=original radiation eqns
! 'OBSFIT'=empirical fit to observations
CHARACTER(LEN=*), PARAMETER :: CLVISIB_ALGOR='OBSFIT' 

!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RADVIS',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------

! For "IFSRAD" option:
! Particle sizes (in microns) for liquid, ice, rain and snow are fixed
! (in the future, could be linked to a calculated particle effective radius)
ZRE_LIQ  = 10._JPRB
ZDE_ICE  = 60._JPRB
ZRE_RAIN = 1000._JPRB
ZDE_SNOW = 2000._JPRB
! Index for 0.55 um (in fact 0.6250 - 0.4415 um)
JSW=10                

! For "OBSFIT" option:
! For Gultepe cloud liquid - fixed number concfor cloud water droplets (cm-3)
! (in the future, could be linked to a calculated particle effective radius)
ZND_LIQ  = 50._JPRB

! Define -ln(eta) where liminal visual contrast eta=0.02
ZLOGNETA = 3.912023_JPRB


!---------------------------------------------------------------------
!
! Extinction coefficient for clear air (Rayleigh scattering)
!
!---------------------------------------------------------------------
! This could be calculated explicitly with code below:  
!   ASSOCIATE(RNS=>YDERAD%RNS, RSIGAIR=>YDERAD%RSIGAIR)
!   Molecular density in first layer
!   ZNS = RNS*273.15_JPRB/PT(JL,KLEV)
!   Rayleigh scattering cross section (cm2 molec-1)
!   ZSIGAIR = (RSIGAIR/ZNS)/ZNS ! not RSIGAIR/(ZNS*ZNS) because intermediate value busts single precision
!   Rayleigh scattering coefficient (km-1)
!   ZRAYLEIGH = 1.E+05_JPRB * ZSIGAIR * ZNS * PSRF1(JL,KLEV) / 101325._JPRB
! but in practice the extinction coefficient of clean air generally has no 
! practical impact, and has been taken to be equivalent to a visibility of 
! 100km  (=1.E5 m) to ensure that unrealistically high visibilities are 
! never diagnosed following Clark et al. (2008,QJ) and Claxton (2008,QJ).
!---------------------------------------------------------------------
! Calculate extinction coefficient from Rayleigh scattering (m-1) 
ZRAYEXT = ZLOGNETA/1.E5_JPRB


DO JL=KIDIA,KFDIA
 
  !---------------------------------------------------------------------
  !
  ! Extinction coefficient for cloud and precipitation
  !
  !---------------------------------------------------------------------
  ! Calculate air density
  ZAIRDENSITY = PRSF1(JL,KLEV)/(RD*PT(JL,KLEV))

  ! Use grid-box mean values to avoid very low visibilities
  ! Multiply by air_density*1000 to convert from kg kg-1 to g m-3
  ZQLWC   = MAX(0._JPRB, PL(JL,KLEV)*ZAIRDENSITY*1000._JPRB)
  ZQRWC   = MAX(0._JPRB, PR(JL,KLEV)*ZAIRDENSITY*1000._JPRB)
  ZQIWC   = MAX(0._JPRB, PI(JL,KLEV)*ZAIRDENSITY*1000._JPRB)
  ZQSWC   = MAX(0._JPRB, PS(JL,KLEV)*ZAIRDENSITY*1000._JPRB)
  
  ! Calculate extinction coefficient for cloud and precipitation (m-1)
  
  ! Use assumptions from the radiation scheme
  IF (CLVISIB_ALGOR == 'IFSRAD') THEN

    ! IFS radiation scheme extinction coefficient for visible band (m-1)
    ZEXTLIQ  = ZQLWC * (RSASWA(JSW) + RSASWB(JSW) / ZRE_LIQ)
    ZEXTICE  = ZQIWC * (RSFUA0(JSW) + RSFUA1(JSW) / ZDE_ICE)
    ZEXTSNOW = ZQSWC * (RSFUA0(JSW) + RSFUA1(JSW) / ZDE_SNOW)
    ZEXTRAIN = ZQRWC * (RSASWA(JSW) + RSASWB(JSW) / ZRE_RAIN)

  ! Use empirical fit to observations
  ELSEIF (CLVISIB_ALGOR == 'OBSFIT') THEN

    ! Cloud water extinction coefficient (m-1) 
    ! Gultepe (2006) Vis = 1.002*(ZQLWC*Nd)^-0.6473 (*10E-3 to convert km to m)
    ZEXTLIQ = (ZLOGNETA/1002._JPRB)*(ZQLWC*ZND_LIQ)**0.6473_JPRB

    ! Ice extinction coefficient (m-1)
    ! Stoelinga and Warner 1999 (*10E-3 to convert km to m)
    ZEXTICE = 163.9E-3_JPRB*ZQIWC

    ! Snow extinction coefficient (m-1) 
    ! Stallabrass 1985; Stoelinga and Warner 1999: 10.4*SWC**0.78
    ! Modified for IFS to better fit data (*10E-3 to convert km to m)
    ZEXTSNOW = 4.E-3_JPRB*ZQSWC**0.78_JPRB

    ! Rain extinction coefficient (m-1)
    ! Stoelinga and Warner 1999
    ! Modified for IFS to better fit data (*10E-3 to convert km to m)
    ZEXTRAIN = 5.E-3_JPRB*ZQRWC**0.75_JPRB

  ELSE
  
    ZEXTLIQ  = 0.0_JPRB
    ZEXTICE  = 0.0_JPRB
    ZEXTSNOW = 0.0_JPRB
    ZEXTRAIN = 0.0_JPRB

  ENDIF

  ! Calculate total extinction coefficient for cloud+precip hydrometeors (m-1)
  ZHYDEXT(JL) = (ZEXTICE + ZEXTLIQ + ZEXTSNOW + ZEXTRAIN)
  
  !---------------------------------------------------------------------
  !
  ! Total visibility (m) =  clear air (Rayleigh) + aerosols + cloud + precip
  !
  ! Hinkley (1976) Vis = (3.91/Sigma) * (0.55/Lambda)^1.3
  ! For visible wavelength, Lambda = 0.55 um and
  ! Sigma is the total extinction coefficient in m-1
  !---------------------------------------------------------------------
  PVISALL(JL)  = ZLOGNETA / (ZRAYEXT + PAEREXT(JL) + ZHYDEXT(JL))

ENDDO

! Only compute partial visibilities if provided as output arguments to the routine
IF (PRESENT(PVISAER)) THEN
  ! Rayleigh plus aerosols
  DO JL=KIDIA,KFDIA
    PVISAER(JL)  = ZLOGNETA / (ZRAYEXT + PAEREXT(JL))
  ENDDO
ENDIF

IF (PRESENT(PVISCLD)) THEN
  ! Rayleigh plus hydrometeors
  DO JL=KIDIA,KFDIA
    PVISCLD(JL)  = ZLOGNETA / (ZRAYEXT + ZHYDEXT(JL))
  ENDDO
ENDIF

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RADVIS',1,ZHOOK_HANDLE)
END SUBROUTINE RADVIS
