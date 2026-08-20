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

MODULE YOMSOLARIRRADIANCE
! PURPOSE
! -------
!   Derived type for reading in a multi-year time series of total
!   solar irradiance (TSI)
!
! INTERFACE
! ---------
!   YSOLARIRRADIANCE%READ is called from SUECRAD and YSOLARIRRADIANCE%GET
!   is called from UPDRGAS.
!
! AUTHOR
! ------
!   Robin Hogan, ECMWF
!   Original: 2019-03-11
!
! MODIFICATIONS
! -------------
!   2022-12-05  R. Hogan  Read solar cycle number if available

USE PARKIND1, ONLY : JPRB, JPRD, JPIM

IMPLICIT NONE

SAVE

! Type to store total solar irradiance time series
TYPE :: TSOLARIRRADIANCE
  ! Year measured from 00 UTC 1 January 0000 (astronomical numbering).
  ! Many climate datasets contain whole-calendar year averages, so
  ! these are given a date in the middle of the year and linear
  ! interpolation between them is performed.  Thus 2010.5 is treated
  ! as 2 July 2010.
  REAL(KIND=JPRB), ALLOCATABLE :: RYEAR(:)

  ! Total solar irradiance at 1 Astronomical Unit
  REAL(KIND=JPRB), ALLOCATABLE :: RTSI(:)

  ! Solar cycle number according to Rudolf Wolf's convention where 1.0
  ! corresponds to solar minimum in February 1755 and subsequent whole
  ! numbers correspond to the subsequent solar minima. This may be
  ! interpolated to the current date and then used to compute the
  ! solar cycle multiplier, used to scale the amplitude of the
  ! spectral solar variability.
  REAL(KIND=JPRB), ALLOCATABLE :: RSOLARCYCLENUM(:)
  
  ! Number of times
  INTEGER(KIND=JPIM) :: NTIME

CONTAINS

  PROCEDURE :: READ => READ_SOLAR_IRRADIANCE
  PROCEDURE :: GET  => GET_SOLAR_IRRADIANCE

END TYPE TSOLARIRRADIANCE

TYPE(TSOLARIRRADIANCE), POINTER :: YSOLARIRRADIANCE => NULL()

CONTAINS


!-----------------------------------------------------------------------
! Read TSI from NetCDF file
SUBROUTINE READ_SOLAR_IRRADIANCE(SELF, CD_FILE_NAME)

  USE EASY_NETCDF_READ_MPI, ONLY : NETCDF_FILE
  USE YOMHOOK,              ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOMLUN,               ONLY : NULOUT

  CLASS(TSOLARIRRADIANCE), INTENT(INOUT) :: SELF

  ! Input file name: if it starts with "." or "/" then it is treated
  ! as an absolute path, otherwise it is assumed to be in the usual
  ! data directory
  CHARACTER(*), INTENT(IN) :: CD_FILE_NAME

  ! Full name of the TSI file
  CHARACTER(LEN=512) :: CL_FILE_NAME
  ! "DATA" directory
  CHARACTER(LEN=512) :: CLDIRECTORY

  ! The NetCDF file containing TSI
  TYPE(NETCDF_FILE)  :: FILE

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOMSOLARIRRADIANCE:READ_SOLAR_IRRADIANCE',0,ZHOOK_HANDLE)

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

  ! Read year coordinate variable
  CALL FILE%GET('time', SELF%RYEAR)
  SELF%NTIME = SIZE(SELF%RYEAR,1)

  ! Read total solar irradiance
  CALL FILE%GET('tsi', SELF%RTSI)

  IF (FILE%EXISTS('solar_cycle_number')) THEN
    CALL FILE%GET('solar_cycle_number', SELF%RSOLARCYCLENUM)
  ENDIF
  
  CALL FILE%CLOSE()

  IF (LHOOK) CALL DR_HOOK('YOMSOLARIRRADIANCE:READ_SOLAR_IRRADIANCE',1,ZHOOK_HANDLE)

END SUBROUTINE READ_SOLAR_IRRADIANCE


!-----------------------------------------------------------------------
! Get TSI at a particular date, specified in terms of decimal year
SUBROUTINE GET_SOLAR_IRRADIANCE(SELF, PYEAR, PTSI, PSOLARCYCLENUM, PSOLARCYCLEMULT)

  USE YOMHOOK,              ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOMLUN,               ONLY : NULOUT
  USE YOMCST,               ONLY : RPI

  ! For time functions NCCAA, NMM and NDD
#include "fcttim.func.h"

  CLASS(TSOLARIRRADIANCE), INTENT(IN)  :: SELF
  REAL(KIND=JPRB),         INTENT(IN)  :: PYEAR

  ! Output total solar irradiance (W m-2)
  REAL(KIND=JPRB),         INTENT(OUT) :: PTSI

  ! Output solar cycle number and multiplier (1=solar max, -1=min,
  ! 0=mean)
  REAL(KIND=JPRB), INTENT(OUT), OPTIONAL :: PSOLARCYCLENUM, PSOLARCYCLEMULT

  ! Index to first interpolation point
  INTEGER(KIND=JPIM) :: IVAL
  ! Weight of second interpolation point (indexed IVAL+1)
  REAL(KIND=JPRB)    :: ZWEIGHT

  REAL(KIND=JPRB)    :: ZSOLARCYCLENUM, ZSOLARCYCLEMULT
  
  ! Loop index for time
  INTEGER(KIND=JPIM) :: JT

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOMSOLARIRRADIANCE:GET_SOLAR_IRRADIANCE',0,ZHOOK_HANDLE)

  ! Linear interpolation clamped at each end
  IF (PYEAR <= SELF%RYEAR(1)) THEN
    ! Time too early
    IVAL = 1
    ZWEIGHT = 0.0_JPRB
    WRITE(NULOUT,'(A,I8.8)') 'YOMSOLARIRRADIANCE:GET_SOLAR_IRRADIANCE year ', &
         &  PYEAR, ' is before earliest year ', SELF%RYEAR(1)
  ELSEIF (PYEAR >= SELF%RYEAR(SELF%NTIME)) THEN
    ! Time too late
    IVAL = SELF%NTIME-1
    ZWEIGHT = 1.0_JPRB
    WRITE(NULOUT,'(A,I8.8)') 'YOMSOLARIRRADIANCE:GET_SOLAR_IRRADIANCE year ', PYEAR, &
         &  ' is after final year ', SELF%RYEAR(SELF%NTIME)
  ELSE
    ! Time is in range
    IVAL = 1
    ZWEIGHT = 0.0_JPRB
    DO JT = 1,SELF%NTIME-1
      IF (PYEAR >= SELF%RYEAR(JT) .AND. PYEAR < SELF%RYEAR(JT+1)) THEN
        IVAL = JT
        ZWEIGHT = (PYEAR - SELF%RYEAR(JT)) / (SELF%RYEAR(JT+1) - SELF%RYEAR(JT))
        EXIT
      ENDIF
    ENDDO
  ENDIF
  PTSI = SELF%RTSI(JT) + ZWEIGHT * (SELF%RTSI(JT+1) - SELF%RTSI(JT))

  IF (ALLOCATED(SELF%RSOLARCYCLENUM)) THEN
    ZSOLARCYCLENUM = SELF%RSOLARCYCLENUM(JT) &
         &  + ZWEIGHT * (SELF%RSOLARCYCLENUM(JT+1) - SELF%RSOLARCYCLENUM(JT))
    IF ((IVAL == 1 .AND. ZWEIGHT == 0.0_JPRB) &
         &  .OR. (IVAL == SELF%NTIME-1 .AND. ZWEIGHT == 1.0_JPRB)) THEN
      ! Out of range: use solar spectral mean
      ZSOLARCYCLEMULT = 0.0_JPRB
    ELSE
      ! Solar minima occur at whole numbers and lead to a multiplier
      ! of -1.0
      ZSOLARCYCLEMULT = -COS(ZSOLARCYCLENUM*2.0_JPRB*RPI)
    ENDIF
  ELSE
    ! Solar cycle number not available: choose a value corresponding
    ! to a midpoint in the solar cycle
    ZSOLARCYCLENUM  = -0.25_JPRB
    ZSOLARCYCLEMULT = 0.0_JPRB
  ENDIF
  
  IF (PRESENT(PSOLARCYCLENUM)) THEN
    PSOLARCYCLENUM = ZSOLARCYCLENUM
  ENDIF

  IF (PRESENT(PSOLARCYCLEMULT)) THEN
    PSOLARCYCLEMULT = ZSOLARCYCLEMULT
  ENDIF
  
  IF (LHOOK) CALL DR_HOOK('YOMSOLARIRRADIANCE:GET_SOLAR_IRRADIANCE',1,ZHOOK_HANDLE)

END SUBROUTINE GET_SOLAR_IRRADIANCE

END MODULE YOMSOLARIRRADIANCE
