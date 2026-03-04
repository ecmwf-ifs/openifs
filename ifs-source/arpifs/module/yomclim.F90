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

MODULE YOMCLIM
!     ------------------------------------------------------------------
!     *YOMCLIM* - STORAGE OF GAS CLIMATOLOGY FOR USE IN RADIATION SCHEME
!     ------------------------------------------------------------------
! Modifications
!   2019-01-25  R. Hogan  Completely rewrote

USE PARKIND1, ONLY : JPRB, JPIM

IMPLICIT NONE

SAVE

! Type to store greenhouse-gas climatology
TYPE :: TGHGCLIM
  ! Pressure and latitude coordinate variables
  REAL(KIND=JPRB), ALLOCATABLE :: PRESSURE(:) ! Pa
  REAL(KIND=JPRB), ALLOCATABLE :: LATITUDE(:) ! Degrees

  ! Gas volume mixing ratios (unitless) versus latitude, pressure and
  ! month
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_CO2(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_CH4(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_N2O(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_NO2(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_CFC11(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_CFC12(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_O3(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_HCFC22(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_CCL4(:,:,:)

  ! Annual/global surface-mean values, in case we wish to scale to a
  ! prescribed value
  REAL(KIND=JPRB) :: SURF_MEAN_CO2, SURF_MEAN_CH4, SURF_MEAN_N2O, SURF_MEAN_NO2
  REAL(KIND=JPRB) :: SURF_MEAN_CFC11, SURF_MEAN_CFC12, SURF_MEAN_HCFC22, SURF_MEAN_CCL4
  
  ! Number of pressures and latitudes
  INTEGER(KIND=JPIM) :: NPRESSURE = 0, NLATITUDE = 0

CONTAINS

  PROCEDURE :: READ => READ_GHG_CLIMATOLOGY

END TYPE TGHGCLIM

TYPE(TGHGCLIM), POINTER :: YGHGCLIM => NULL()

CONTAINS

!-----------------------------------------------------------------------
! Read greenhouse gas climatology from NetCDF file
SUBROUTINE READ_GHG_CLIMATOLOGY(SELF, CD_FILE_NAME)

  USE EASY_NETCDF_READ_MPI,ONLY : NETCDF_FILE
  USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOMLUN    ,ONLY : NULERR
  USE YOMCST    ,ONLY : RPI

  CLASS(TGHGCLIM), INTENT(INOUT) :: SELF

  ! Input file name: if it starts with "." or "/" then it is treated
  ! as an absolute path, otherwise it is assumed to be in the usual
  ! data directory
  CHARACTER(*), INTENT(IN) :: CD_FILE_NAME

  ! Full name of the gas climatology file
  CHARACTER(LEN=512) :: CL_FILE_NAME
  ! "DATA" directory
  CHARACTER(LEN=512) :: CLDIRECTORY

  ! The NetCDF file containing the input climatology
  TYPE(NETCDF_FILE)  :: FILE

  ! Latitude weight to compute global mean
  REAL(KIND=JPRB), ALLOCATABLE :: ZWEIGHT(:)

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
  !-----------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('YOMCLIM:READ_GHG_CLIMATOLOGY',0,ZHOOK_HANDLE)

  IF (CD_FILE_NAME(1:1) /= "." .AND. CD_FILE_NAME(1:1) /= "/") THEN
    ! Prefix file name by the data directory name from the "DATA"
    ! environment variable, if present
    CALL GET_ENVIRONMENT_VARIABLE("DATA",CLDIRECTORY)
    IF (CLDIRECTORY /= " ") THEN
      CL_FILE_NAME = TRIM(CLDIRECTORY) // "/ifsdata/" // CD_FILE_NAME
    ELSE
      CL_FILE_NAME = TRIM(CD_FILE_NAME)
    ENDIF
  ELSE
    ! File name starts with "." or "/": treat as an absolute path
    CL_FILE_NAME = TRIM(CD_FILE_NAME)
  ENDIF

  CALL FILE%OPEN(TRIM(CL_FILE_NAME), IVERBOSE=4)

  ! Read coordinate variables
  CALL FILE%GET('pressure', SELF%PRESSURE)
  CALL FILE%GET('latitude', SELF%LATITUDE)

  ! Read first gas
  CALL FILE%GET('co2_vmr',  SELF%VMR_CO2)

  ! Check size of first gas array
  SELF%NPRESSURE = SIZE(SELF%PRESSURE)
  SELF%NLATITUDE = SIZE(SELF%LATITUDE)
  IF (SIZE(SELF%VMR_CO2,1) /= SELF%NLATITUDE &
       & .OR. SIZE(SELF%VMR_CO2,2) /= SELF%NPRESSURE &
       & .OR. SIZE(SELF%VMR_CO2,3) /= 12) THEN
    WRITE(NULERR,*) 'Error reading ', TRIM(CL_FILE_NAME), ': expected co2_vmr(', SELF%NLATITUDE, &
         &  ',', SELF%NPRESSURE, ',12) but got co2_vmr(', SIZE(SELF%VMR_CO2,1), &
         &  ',', SIZE(SELF%VMR_CO2,2), ',', SIZE(SELF%VMR_CO2,3), ')'
    CALL ABOR1('YOMCLIM:READ_GHG_CLIMATOLOGY: Inconsistency in gas climatology file')
  ENDIF
  
  ! Read remaining gases
  CALL FILE%GET('ch4_vmr',   SELF%VMR_CH4)
  CALL FILE%GET('n2o_vmr',   SELF%VMR_N2O)
  CALL FILE%GET('cfc11_vmr', SELF%VMR_CFC11)
  CALL FILE%GET('cfc12_vmr', SELF%VMR_CFC12)
  CALL FILE%GET('o3_vmr',    SELF%VMR_O3)
  CALL FILE%GET('hcfc22_vmr',SELF%VMR_HCFC22)
  CALL FILE%GET('ccl4_vmr',  SELF%VMR_CCL4)
  
  ! NO2 is ignored
  IF(ALLOCATED(SELF%VMR_NO2))DEALLOCATE(SELF%VMR_NO2)
  ALLOCATE(SELF%VMR_NO2(SELF%NLATITUDE, SELF%NPRESSURE, 12))
  SELF%VMR_NO2 = 0.0_JPRB
  
  CALL FILE%CLOSE()
  
  ! Compute surface-mean values
  IF(ALLOCATED(ZWEIGHT))DEALLOCATE(ZWEIGHT)
  ALLOCATE(ZWEIGHT(SELF%NLATITUDE))
  ZWEIGHT = COS(SELF%LATITUDE * (RPI / 180.0_JPRB))
  ZWEIGHT = ZWEIGHT / (12.0_JPRB * SUM(ZWEIGHT))

  SELF%SURF_MEAN_CO2   = SUM(SUM(SELF%VMR_CO2(:,SELF%NPRESSURE,:),2)*ZWEIGHT)
  SELF%SURF_MEAN_CH4   = SUM(SUM(SELF%VMR_CH4(:,SELF%NPRESSURE,:),2)*ZWEIGHT)
  SELF%SURF_MEAN_N2O   = SUM(SUM(SELF%VMR_N2O(:,SELF%NPRESSURE,:),2)*ZWEIGHT)
  SELF%SURF_MEAN_NO2   = SUM(SUM(SELF%VMR_NO2(:,SELF%NPRESSURE,:),2)*ZWEIGHT)
  SELF%SURF_MEAN_CFC11 = SUM(SUM(SELF%VMR_CFC11(:,SELF%NPRESSURE,:),2)*ZWEIGHT)
  SELF%SURF_MEAN_CFC12 = SUM(SUM(SELF%VMR_CFC12(:,SELF%NPRESSURE,:),2)*ZWEIGHT)
  !SELF%SURF_MEAN_O3    = SUM(SUM(SELF%VMR_O3(:,SELF%NPRESSURE,:),2)*ZWEIGHT)
  SELF%SURF_MEAN_HCFC22= SUM(SUM(SELF%VMR_HCFC22(:,SELF%NPRESSURE,:),2)*ZWEIGHT)
  SELF%SURF_MEAN_CCL4  = SUM(SUM(SELF%VMR_CCL4(:,SELF%NPRESSURE,:),2)*ZWEIGHT)

  DEALLOCATE(ZWEIGHT)

  IF (LHOOK) CALL DR_HOOK('YOMCLIM:READ_GHG_CLIMATOLOGY',1,ZHOOK_HANDLE)

END SUBROUTINE READ_GHG_CLIMATOLOGY


END MODULE YOMCLIM

