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

MODULE YOERADGHG

!     -------------------------------------------------------------------
!     *YOERADGHG* - STORAGE OF GAS CLIMATOLOGY INTERPOLATED TO MODEL TIME
!     -------------------------------------------------------------------
! Modifications
!   2019-01-25  R. Hogan  Completely rewrote

USE PARKIND1  ,ONLY : JPRB, JPIM

IMPLICIT NONE

! Type to store greenhouse-gas climatology interpolated to model time
TYPE :: TRADGHG
  ! Sine of latitude where monthly mean climatologies are provided
  REAL(KIND=JPRB), ALLOCATABLE :: SINLAT(:)
  ! Pressure (Pa) at half levels *indexed from 0* up to the number of layers
  REAL(KIND=JPRB), ALLOCATABLE :: PRESSURE_HL(:)
  ! Layer-integrated gas mass concentrations (Pa, i.e. mass mixing
  ! ratio multiplied by pressure difference across the layer),
  ! dimensioned (NLATITUDE,NPRESSURE), having already been
  ! interpolated in time
  REAL(KIND=JPRB), ALLOCATABLE :: MASS_CH4(:,:),  MASS_CO2(:,:),    MASS_N2O(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: MASS_NO2(:,:),  MASS_CFC11(:,:),  MASS_CFC12(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: MASS_CCL4(:,:), MASS_HCFC22(:,:), MASS_O3(:,:)

  ! Number of pressures and latitudes
  INTEGER(KIND=JPIM) :: NLATITUDE = 0, NPRESSURE = 0

CONTAINS

  PROCEDURE :: INIT => INIT_RAD_GHG

END TYPE TRADGHG

TYPE(TRADGHG), POINTER :: YRADGHG => NULL()

CONTAINS

!-----------------------------------------------------------------------
! Initialize pressure and latitude, and allocate space for gases
SUBROUTINE INIT_RAD_GHG(SELF, PLATITUDE, PPRESSURE)

  USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOMCST    ,ONLY : RPI

  CLASS(TRADGHG),  INTENT(INOUT) :: SELF
  REAL(KIND=JPRB), INTENT(IN)    :: PLATITUDE(:), PPRESSURE(:)

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOERADGHG:INIT_RAD_GHG',0,ZHOOK_HANDLE)

  ! We only do anything if the type has not already been initialized
  IF (SELF%NLATITUDE <= 0) THEN
    SELF%NLATITUDE = SIZE(PLATITUDE)
    SELF%NPRESSURE = SIZE(PPRESSURE)
    ALLOCATE(SELF%SINLAT    (SELF%NLATITUDE))
    ALLOCATE(SELF%PRESSURE_HL(0:SELF%NPRESSURE))
    ALLOCATE(SELF%MASS_CH4   (SELF%NLATITUDE,SELF%NPRESSURE))
    ALLOCATE(SELF%MASS_CO2   (SELF%NLATITUDE,SELF%NPRESSURE))
    ALLOCATE(SELF%MASS_N2O   (SELF%NLATITUDE,SELF%NPRESSURE))
    ALLOCATE(SELF%MASS_NO2   (SELF%NLATITUDE,SELF%NPRESSURE))
    ALLOCATE(SELF%MASS_CFC11 (SELF%NLATITUDE,SELF%NPRESSURE))
    ALLOCATE(SELF%MASS_CFC12 (SELF%NLATITUDE,SELF%NPRESSURE))
    ALLOCATE(SELF%MASS_CCL4  (SELF%NLATITUDE,SELF%NPRESSURE))
    ALLOCATE(SELF%MASS_HCFC22(SELF%NLATITUDE,SELF%NPRESSURE))
    ALLOCATE(SELF%MASS_O3    (SELF%NLATITUDE,SELF%NPRESSURE))

    ! Initialize sine of latitude
    SELF%SINLAT = SIN(PLATITUDE * (RPI/180.0_JPRB))

    ! Initialize pressure at half levels
    SELF%PRESSURE_HL(0) = 0.0_JPRB
    SELF%PRESSURE_HL(1:SELF%NPRESSURE-1) = 0.5_JPRB * (PPRESSURE(1:SELF%NPRESSURE-1) + PPRESSURE(2:SELF%NPRESSURE))
    SELF%PRESSURE_HL(SELF%NPRESSURE) = MAX(110000.0_JPRB, 1.5_JPRB*PPRESSURE(SELF%NPRESSURE) &
         &                                               -0.5_JPRB*PPRESSURE(SELF%NPRESSURE-1))
  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOERADGHG:INIT_RAD_GHG',1,ZHOOK_HANDLE)

END SUBROUTINE INIT_RAD_GHG

END MODULE YOERADGHG


