! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE YOMCLMICST

USE PARKIND1,ONLY : JPRB
USE YOMHOOK ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST  ,ONLY : RPI
USE YOMLUN  ,ONLY : NULOUT

IMPLICIT NONE

SAVE

REAL(KIND=JPRB) :: RXAR,RXBR,RXCCR,RXAI,RXBI,RXAS,RXBS,RXCCS,RXCXS
REAL(KIND=JPRB) :: RXAG,RXBG,RXCCG,RXCXG,RXALPHAR,RXNUR,RXALPHAI,RXNUI
REAL(KIND=JPRB) :: RXALPHAS,RXNUS,RXALPHAG,RXNUG,RXLBEXR,RXLBR,RXLBEXI
REAL(KIND=JPRB) :: RXLBI, RXLBEXS,RXLBS,RXLBEXG,RXLBG
REAL(KIND=JPRB) :: RXRTMIN(6)
REAL(KIND=JPRB) :: RMOMG_RXALPHAR_RXNUR_6
REAL(KIND=JPRB) :: RMOMG_RXALPHAI_RXNUI_ZEXP
REAL(KIND=JPRB) :: RMOMG_RXALPHAS_RXNUS_ZEXP
REAL(KIND=JPRB) :: RMOMG_RXALPHAG_RXNUG_ZEXP

CONTAINS
  REAL FUNCTION MOMG(PALPHA,PNU,PP)
! auxiliary routine used to compute the Pth moment order of the generalized
! gamma law

  IMPLICIT NONE

  REAL(KIND=JPRB),INTENT(IN)     :: PALPHA ! first shape parameter of the dimensionnal distribution
  REAL(KIND=JPRB),INTENT(IN)     :: PNU    ! second shape parameter of the dimensionnal distribution
  REAL(KIND=JPRB),INTENT(IN)     :: PP     ! order of the moment
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  REAL(KIND=JPRB) :: FCGENERALIZED_GAMMA

  IF (LHOOK) CALL DR_HOOK('YOMCLMICST:MOMG',0,ZHOOK_HANDLE)

  MOMG = FCGENERALIZED_GAMMA(PNU+PP/PALPHA)/FCGENERALIZED_GAMMA(PNU)

  IF (LHOOK) CALL DR_HOOK('YOMCLMICST:MOMG',1,ZHOOK_HANDLE)
   
  END FUNCTION MOMG
! ===========================================================================
SUBROUTINE SETUP_CLMICST

!**** *SETUP_CLMICST*  - Initialise constants of microphysics scheme ICE3 of MNH-
!                     this is needed in case of computation of microphysics variables

!     PURPOSE.
!     --------
!        To Initialise constants related to cloud 

!**   INTERFACE.
!     ----------
!       *CALL* *SETUP_CLMICST*

!        EXPLICIT ARGUMENTS
!        --------------------
!            INPUT : None

!            OUTPUT: None

!        IMPLICIT ARGUMENTS
!        --------------------
!           NONE

!     METHOD.
!     -------
!        Duplicated from ini_rain_ice of MNH code

!     EXTERNALS.
!     ----------
   
!     REFERENCE.
!     ----------
!        Meso-Nh scientific documentation 

!     AUTHOR.
!     -------
!        Gwenaelle Hello *METEO-FRANCE*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 06-07-31
!     ------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPRB) :: ZRHOLW,ZEXP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMCLMICST:SETUP_CLMICST',0,ZHOOK_HANDLE)

!*********************************************************
! 0. Initialization  
!*********************************************************


ZRHOLW=1000.0_JPRB ! density of liquid water (kgm-3)

!*********************************************************
! 1. Define Characteristics of the species
!*********************************************************
!*       1.1    Raindrops characteristics

RXAR = (RPI/6.0_JPRB)*ZRHOLW
RXBR = 3.0_JPRB
RXCCR = 8.E6_JPRB    ! N0_r with XCXR = -1, implicitly

!*       1.2    Ice crystal characteristics, Plates case

RXAI = 0.82_JPRB      ! Plates
RXBI = 2.5_JPRB       ! Plates 

!*       1.3    Snowflakes/aggregates characteristics

RXAS = 0.02_JPRB
RXBS = 1.9_JPRB

RXCCS = 5.0_JPRB
RXCXS = 1.0_JPRB

!*       1.4    Heavily rimed crystals characteristics Graupel case

RXAG = 19.6_JPRB  ! Lump graupel case
RXBG = 2.8_JPRB   ! Lump graupel case

RXCCG = 5.E5_JPRB
RXCXG = -0.5_JPRB

!*********************************************************
! 2. Define dimensional distributions of the species
!*********************************************************

!*       2.1    Raindrops distribution

RXALPHAR = 1.0_JPRB  ! Exponential law
RXNUR    = 1.0_JPRB  ! Exponential law

!*       2.2    Ice crystal distribution

RXALPHAI = 3.0_JPRB  ! Gamma law for the ice crystal volume
RXNUI    = 3.0_JPRB  ! Gamma law with little dispersion

RXALPHAS = 1.0_JPRB  ! Exponential law
RXNUS    = 1.0_JPRB  ! Exponential law

RXALPHAG = 1.0_JPRB  ! Exponential law
RXNUG    = 1.0_JPRB  ! Exponential law

!*       2.3    Constants for shape parameter

RXLBEXR = 1.0_JPRB/(-1.0_JPRB-RXBR)
RXLBR   = ( RXAR*RXCCR*MOMG(RXALPHAR,RXNUR,RXBR) )**(-RXLBEXR)

RXLBEXI = 1.0_JPRB/(-RXBI)
RXLBI   = ( RXAI*MOMG(RXALPHAI,RXNUI,RXBI) )**(-RXLBEXI)

RXLBEXS = 1.0_JPRB/(RXCXS-RXBS)
RXLBS   = ( RXAS*RXCCS*MOMG(RXALPHAS,RXNUS,RXBS) )**(-RXLBEXS)

RXLBEXG = 1.0_JPRB/(RXCXG-RXBG)
RXLBG   = ( RXAG*RXCCG*MOMG(RXALPHAG,RXNUG,RXBG))**(-RXLBEXG)

!*       2.4    Minimal values allowed for the mixing ratios

RXRTMIN(1) = 1.0E-20_JPRB
RXRTMIN(2) = 1.0E-20_JPRB
RXRTMIN(3) = 1.0E-10_JPRB
RXRTMIN(4) = 1.0E-10_JPRB
RXRTMIN(5) = 1.0E-15_JPRB
RXRTMIN(6) = 1.0E-15_JPRB

RMOMG_RXALPHAR_RXNUR_6=MOMG(RXALPHAR,RXNUR,6.0_JPRB)
!ZEXP = 7.0_JPRB*RXBI/3.0_JPRB-1.0_JPRB
ZEXP = 2.0_JPRB*RXBI
RMOMG_RXALPHAI_RXNUI_ZEXP=MOMG(RXALPHAI,RXNUI,ZEXP)
!ZEXP = 7.0_JPRB*RXBS/3.0_JPRB-1.0_JPRB
ZEXP = 2.0_JPRB*RXBS
RMOMG_RXALPHAS_RXNUS_ZEXP=MOMG(RXALPHAS,RXNUS,ZEXP)
!ZEXP = 7.0_JPRB*RXBG/3.0_JPRB-1.0_JPRB
ZEXP = 2.0_JPRB*RXBG
RMOMG_RXALPHAG_RXNUG_ZEXP=MOMG(RXALPHAG,RXNUG,ZEXP)

!********************************************************
! 3. Some prints
!********************************************************

WRITE(UNIT=NULOUT,FMT='(" ** SETUP_CLMICST - init of species for use outside aroM physics - ** ")')
WRITE(UNIT=NULOUT,FMT='(" Summary of the ice particule characteristics")')
WRITE(UNIT=NULOUT,FMT='("      PRISTINE ICE")')
WRITE(UNIT=NULOUT,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
                                                  &  RXAI,RXBI
WRITE(UNIT=NULOUT,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
                                                  &  RXALPHAI,RXNUI
WRITE(UNIT=NULOUT,FMT='("              SNOW")')
WRITE(UNIT=NULOUT,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
                                                  &  RXAS,RXBS
WRITE(UNIT=NULOUT,FMT='("           concentration:CC=",E13.6," x=",E13.6)') &
                                                  &  RXCCS,RXCXS
WRITE(UNIT=NULOUT,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
                                                  &  RXALPHAS,RXNUS
WRITE(UNIT=NULOUT,FMT='("            GRAUPEL")')
WRITE(UNIT=NULOUT,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
                                                  &  RXAG,RXBG
WRITE(UNIT=NULOUT,FMT='("           concentration:CC=",E13.6," x=",E13.6)') &
                                                  &  RXCCG,RXCXG
WRITE(UNIT=NULOUT,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
                                                    &  RXALPHAG,RXNUG

IF (LHOOK) CALL DR_HOOK('YOMCLMICST:SETUP_CLMICST',1,ZHOOK_HANDLE)

END SUBROUTINE  SETUP_CLMICST
END MODULE YOMCLMICST
