! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CLAERVIS &
  &( KIDIA  , KFDIA , KLON , KLEV, &
  &  PCLAERS, PA, PL, PI, PR, PS , PRSF1  , PT  , &
  &  PVISICL  &
  & )

!**** *CLAERVIS* - ROUTINE COMPUTING VISIBILITY

!      J.-J. MORCRETTE  R.M.FORBES, ECMWF

!**    INTERFACE
!      ---------
!         *CLAERVIS* IS CALLED FROM *CLIMAER_LAYER*

! INPUTS:
! -------

! OUTPUTS:
! --------

!      AUTHOR.
!      -------
!      Original: JJMorcrette, 20120801

!      Modifications:
!      --------------
!      R. Forbes  01-Mar-2014  Corrected input variables
!      K. Yessad (July 2014)   Move some variables.
!      R. Forbes  01-Jun-2015  Changed assumed particle sizes and use grid means
!      R. Forbes  01-Nov-2020  Tidied code
!      R. Forbes  01-Nov-2020  Changed visibility eqns for liq/ice/rain/snow
!-------------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM , JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RD
USE YOESRTCOP, ONLY : RSASWA, RSASWB, RSFUA0, RSFUA1

!-------------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV

! Input
REAL(KIND=JPRB),INTENT(IN)    :: PRSF1(KLON,KLEV) ! pressure on full levels (Pa)
REAL(KIND=JPRB),INTENT(IN)    :: PT(KLON,KLEV)    ! temperature (K)
REAL(KIND=JPRB),INTENT(IN)    :: PA(KLON,KLEV)    ! cloud fraction (0-1)
REAL(KIND=JPRB),INTENT(IN)    :: PL(KLON,KLEV)    ! cloud liquid (kg kg-1)
REAL(KIND=JPRB),INTENT(IN)    :: PI(KLON,KLEV)    ! cloud ice (kg kg-1) 
REAL(KIND=JPRB),INTENT(IN)    :: PR(KLON,KLEV)    ! rain (kg kg-1)
REAL(KIND=JPRB),INTENT(IN)    :: PS(KLON,KLEV)    ! snow (kg kg-1)
REAL(KIND=JPRB),INTENT(IN)    :: PCLAERS(KLON)    ! aerosol extinction coeff (m-1)

! Output                
REAL(KIND=JPRB),INTENT(OUT)   :: PVISICL(KLON)    ! visibility (m)

!-------------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JL, JSW

REAL(KIND=JPRB)    :: ZAIRDENSITY, ZLOGNETA
REAL(KIND=JPRB)    :: ZRE_LIQ ,ZDE_ICE ,ZRE_RAIN ,ZDE_SNOW
REAL(KIND=JPRB)    :: ZND_LIQ
REAL(KIND=JPRB)    :: ZEXTLIQ  ,ZEXTICE  ,ZEXTRAIN  ,ZEXTSNOW
REAL(KIND=JPRB)    :: ZQLWC   ,ZQIWC   ,ZQRWC    ,ZQSWC    
REAL(KIND=JPRB)    :: ZEXTCOEFF_RAYLEIGH
REAL(KIND=JPRB)    :: ZEXTCOEFF_AEROSOL
REAL(KIND=JPRB)    :: ZEXTCOEFF_CONDENSATE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! Visibility algorithm options: 
! 'IFSRAD'=original radiation eqns
! 'OBSFIT'=empirical fit to observations
CHARACTER(LEN=*),PARAMETER :: CLVISIB_ALGOR='OBSFIT' 

!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CLAERVIS',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------

! For IFSRAD option:
! Particle sizes (in microns) for liquid, ice, rain and snow are fixed
ZRE_LIQ  = 10._JPRB !30._JPRB
ZDE_ICE  = 60._JPRB
ZRE_RAIN = 1000._JPRB
ZDE_SNOW = 2000._JPRB

! For Gultepe cloud liquid - fixed number concfor cloud water droplets (cm-3)
ZND_LIQ  = 50._JPRB

! Choose index for 0.55 um (in fact 0.6250 - 0.4415 um)
JSW=10                

! Define -ln(eta) where liminal visual contrast eta=0.02
ZLOGNETA = 3.912023_JPRB


!---------------------------------------------------------------------
!
! Extinction coefficient for clear air (Rayleigh scattering)
!
!---------------------------------------------------------------------
! This could be calculated explicitly with code below:  
!   Molecular density in first layer
!   ZNS = RNS*273.15_JPRB/PT(JL,KLEV)
!   Rayleigh scattering cross section (cm2 molec-1)
!   ZSIGAIR = RSIGAIR/(ZNS*ZNS)
!   Rayleigh scattering coefficient (km-1)
!   ZRAYLEIGH = 1.E+05_JPRB * ZSIGAIR * ZNS * PSRF1(JL,KLEV) / 101325._JPRB
! but in practice the extinction coefficient of clean air generally has no 
! practical impact, and has been taken to be equivalent to a visibility of 
! 100km  (=1.E5 m) to ensure that unrealistically high visibilities are 
! never diagnosed following Clark et al. (2008,QJ) and Claxton (2008,QJ).
!---------------------------------------------------------------------
! Calculate extinction coefficient from Rayleigh scattering (m-1) 
ZEXTCOEFF_RAYLEIGH = ZLOGNETA/1.E5_JPRB


DO JL=KIDIA,KFDIA

  !---------------------------------------------------------------------
  !
  ! Extinction coefficient for aerosol
  !
  !---------------------------------------------------------------------
  ! Convert from km-1 to m-1
  ZEXTCOEFF_AEROSOL = PCLAERS(JL)/1000._JPRB
  
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
  IF (CLVISIB_ALGOR == 'IFSRAD') THEN

    ! IFS radiation scheme extinction coefficient for visible band (m-1)
    ZEXTLIQ  = ZQLWC * (RSASWA(JSW) + RSASWB(JSW) / ZRE_LIQ)
    ZEXTICE  = ZQIWC * (RSFUA0(JSW) + RSFUA1(JSW) / ZDE_ICE)
    ZEXTSNOW = ZQSWC * (RSFUA0(JSW) + RSFUA1(JSW) / ZDE_SNOW)
    ZEXTRAIN = ZQRWC * (RSASWA(JSW) + RSASWB(JSW) / ZRE_RAIN)

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

  ! Calculate total extinction coefficient for cloud+precipitation (m-1)
  ZEXTCOEFF_CONDENSATE = (ZEXTICE + ZEXTLIQ + ZEXTSNOW + ZEXTRAIN)
  
  !---------------------------------------------------------------------
  !
  ! Total visibility (m) =  clear air (Rayleigh) + aerosols + cloud/precip
  !
  ! Hinkley (1976) Vis = (3.91/Sigma) * (0.55/Lambda)^1.3
  ! For visible wavelength, Lambda = 0.55 um and
  ! Sigma is the total extinction coefficient in m-1
  !---------------------------------------------------------------------
  PVISICL(JL)  = ZLOGNETA / &
               & (ZEXTCOEFF_RAYLEIGH + ZEXTCOEFF_AEROSOL + ZEXTCOEFF_CONDENSATE)

ENDDO

!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CLAERVIS',1,ZHOOK_HANDLE)
END SUBROUTINE CLAERVIS
