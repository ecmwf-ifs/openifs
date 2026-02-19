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

MODULE YOMGHGTIMESERIES
! PURPOSE
! -------
!   Derived type for reading in a multi-year time series of surface
!   greenhouse gas concentrations, used to scale the monthly
!   climatology for use in the radiation scheme.
!
! INTERFACE
! ---------
!   YGHGTIMESERIES%READ is called from SUECRAD and YGHGTIMESERIES%GET
!   is called from UPDRGAS.
!
! AUTHOR
! ------
!   Robin Hogan, ECMWF
!   Original: 2019-03-05
!
! MODIFICATIONS
! -------------

USE PARKIND1, ONLY : JPRB, JPRD, JPIM

IMPLICIT NONE

SAVE

! Type to store greenhouse-gas time series
TYPE :: TGHGTIMESERIES
  ! Year measured from 00 UTC 1 January 0000 (astronomical numbering).
  ! Many climate datasets contain whole-calendar year averages, so
  ! these are given a date in the middle of the year and linear
  ! interpolation between them is performed.  Thus 2010.5 is treated
  ! as 2 July 2010.
  REAL(KIND=JPRB), ALLOCATABLE :: RYEAR(:)

  ! Volume mixing ratio of the greenhouse gases represented in the
  ! radiation scheme. Note that the input time series may use
  ! "equivalent CFC11" to account for the combined radiative impact of
  ! many minor gases as well, usually including HCFC22 and CCl4, but
  ! not CFC12.  In this case, zero will be used here for the
  ! concentrations of HCFC22 and CCl4.
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_CO2(:)
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_CH4(:)
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_N2O(:)
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_CFC11(:)
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_CFC12(:)
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_HCFC22(:)
  REAL(KIND=JPRB), ALLOCATABLE :: VMR_CCL4(:)

  ! Number of times
  INTEGER(KIND=JPIM) :: NTIME

CONTAINS

  PROCEDURE :: READ => READ_GHG_TIMESERIES
  PROCEDURE :: GET  => GET_GHG_CONCENTRATIONS

END TYPE TGHGTIMESERIES

TYPE(TGHGTIMESERIES), POINTER :: YGHGTIMESERIES => NULL()

CONTAINS


!-----------------------------------------------------------------------
! Read greenhouse gas time series from NetCDF file
SUBROUTINE READ_GHG_TIMESERIES(SELF, CD_FILE_NAME)

  USE EASY_NETCDF_READ_MPI, ONLY : NETCDF_FILE
  USE YOMHOOK,              ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOMLUN,               ONLY : NULOUT

  CLASS(TGHGTIMESERIES), INTENT(INOUT) :: SELF

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

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOMGHGTIMESERIES:READ_GHG_TIMESERIES',0,ZHOOK_HANDLE)

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
  CALL FILE%GET('time',  SELF%RYEAR)
  SELF%NTIME = SIZE(SELF%RYEAR,1)

  ! Read essential gases
  CALL FILE%GET('co2_vmr',   SELF%VMR_CO2)
  CALL FILE%GET('ch4_vmr',   SELF%VMR_CH4)
  CALL FILE%GET('n2o_vmr',   SELF%VMR_N2O)
  CALL FILE%GET('cfc11_vmr', SELF%VMR_CFC11)
  CALL FILE%GET('cfc12_vmr', SELF%VMR_CFC12)

  ! Read optional gases, assumed zero if not present
  IF (FILE%EXISTS('hcfc22_vmr')) THEN
    CALL FILE%GET('hcfc22_vmr',SELF%VMR_HCFC22)
  ELSE
    IF(.NOT.ALLOCATED(SELF%VMR_HCFC22))ALLOCATE(SELF%VMR_HCFC22(SELF%NTIME))
    SELF%VMR_HCFC22(:) = 0.0_JPRB
    IF (FILE%IS_MASTER_TASK) THEN
      WRITE(NULOUT,'(A)') '  hcfc22_vmr not present: setting concentration to zero'
    ENDIF
  ENDIF

  IF (FILE%EXISTS('ccl4_vmr')) THEN
    CALL FILE%GET('ccl4_vmr',SELF%VMR_CCL4)
  ELSE
    IF(.NOT.ALLOCATED(SELF%VMR_CCL4))ALLOCATE(SELF%VMR_CCL4(SELF%NTIME))
    SELF%VMR_CCL4(:) = 0.0_JPRB
    IF (FILE%IS_MASTER_TASK) THEN
      WRITE(NULOUT,'(A)') '  ccl4_vmr not present: setting concentration to zero'
    ENDIF
  ENDIF

  CALL FILE%CLOSE()

  IF (LHOOK) CALL DR_HOOK('YOMGHGTIMESERIES:READ_GHG_TIMESERIES',1,ZHOOK_HANDLE)

END SUBROUTINE READ_GHG_TIMESERIES


!-----------------------------------------------------------------------
! Get greenhouse gas concentrations at a particular date, specified in
! terms of decimal year

SUBROUTINE GET_GHG_CONCENTRATIONS(SELF, PYEAR, PCO2_VMR, PCH4_VMR, PN2O_VMR, &
     &                            PCFC11_VMR, PCFC12_VMR, PHCFC22_VMR, PCCL4_VMR)

  USE YOMHOOK,              ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOMLUN,               ONLY : NULOUT

  CLASS(TGHGTIMESERIES), INTENT(IN)  :: SELF
  REAL(KIND=JPRB),       INTENT(IN)  :: PYEAR
  ! Output gas volume mixing ratios
  REAL(KIND=JPRB),       INTENT(OUT) :: PCO2_VMR, PCH4_VMR, PN2O_VMR
  REAL(KIND=JPRB),       INTENT(OUT) :: PCFC11_VMR, PCFC12_VMR, PHCFC22_VMR, PCCL4_VMR

  ! Index to first interpolation point
  INTEGER(KIND=JPIM) :: IVAL
  ! Weight of second interpolation point (indexed IVAL+1)
  REAL(KIND=JPRB)    :: ZWEIGHT

  ! Loop index for time
  INTEGER(KIND=JPIM) :: JT

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOMGHGTIMESERIES:GET_GHG_CONCENTRATIONS',0,ZHOOK_HANDLE)

  ! Linear interpolation clamped at each end
  IF (PYEAR <= SELF%RYEAR(1)) THEN
    ! Time too early
    IVAL = 1
    ZWEIGHT = 0.0_JPRB
    WRITE(NULOUT,'(A,F0.3)') 'YOMGHGTIMESERIES:GET_GHG_CONCENTRATIONS year ', &
         &  PYEAR, ' is before earliest year ', SELF%RYEAR(1)
  ELSEIF (PYEAR >= SELF%RYEAR(SELF%NTIME)) THEN
    ! Time too late
    IVAL = SELF%NTIME-1
    ZWEIGHT = 1.0_JPRB
    WRITE(NULOUT,'(A,F0.3)') 'YOMGHGTIMESERIES:GET_GHG_CONCENTRATIONS year ', &
         &  PYEAR, ' is after final year ', SELF%RYEAR(SELF%NTIME)
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
  PCO2_VMR    = SELF%VMR_CO2(JT)    + ZWEIGHT * (SELF%VMR_CO2(JT+1)    - SELF%VMR_CO2(JT))
  PCH4_VMR    = SELF%VMR_CH4(JT)    + ZWEIGHT * (SELF%VMR_CH4(JT+1)    - SELF%VMR_CH4(JT))
  PN2O_VMR    = SELF%VMR_N2O(JT)    + ZWEIGHT * (SELF%VMR_N2O(JT+1)    - SELF%VMR_N2O(JT))
  PCFC11_VMR  = SELF%VMR_CFC11(JT)  + ZWEIGHT * (SELF%VMR_CFC11(JT+1)  - SELF%VMR_CFC11(JT))
  PCFC12_VMR  = SELF%VMR_CFC12(JT)  + ZWEIGHT * (SELF%VMR_CFC12(JT+1)  - SELF%VMR_CFC12(JT))
  PHCFC22_VMR = SELF%VMR_HCFC22(JT) + ZWEIGHT * (SELF%VMR_HCFC22(JT+1) - SELF%VMR_HCFC22(JT))
  PCCL4_VMR   = SELF%VMR_CCL4(JT)   + ZWEIGHT * (SELF%VMR_CCL4(JT+1)   - SELF%VMR_CCL4(JT))
  
  IF (LHOOK) CALL DR_HOOK('YOMGHGTIMESERIES:GET_GHG_CONCENTRATIONS',1,ZHOOK_HANDLE)

END SUBROUTINE GET_GHG_CONCENTRATIONS



END MODULE YOMGHGTIMESERIES
