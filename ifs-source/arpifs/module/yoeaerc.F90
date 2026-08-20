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

MODULE YOEAERC

! YOEAERC - Storage of aerosol climatology for radiation scheme
!
! PURPOSE
! -------
!   The TAERC structure stores an aerosol climatology, typically in
!   the form of monthly averages of multiple species with 3x3 degree
!   lon-lat resolution. The original data may be spatially 2D in which
!   case they contain only column integrated mass, and must be
!   distributed in height, or fully 3D. Originally the IFS used an
!   aerosol climatology derived from the MACC project, which was later
!   upgraded to be from free-running CAMS control run (experiment
!   gbst), bias corrected with total aerosol optical depth
!   CAMS-Interim reanalysis 2003-2011 (experiment eac4).
!
! INTERFACE
! ---------
 
!   All information is stored in a TAERC structure. The data are read
!   in using the TAERC%SETUP routine.
!
! AUTHORS
! -------
!   R. Hogan and A. Bozzo, ECMWF 2016-09-28
!
! MODIFICATIONS
! -------------
!   2018-01-22  A. Bozzo  Added 3D climatology
!   2019-02-02  R. Hogan  Open/close file once, rather than for each variable
!   2022-01-26  R. Hogan  Moved Tegen to yoeaerc_tegen.F90
!   2022-01-26  R. Hogan  Generalized the species - no longer hard-coded
!   2022-04-15  R. Hogan  Option to read 4D aerosol climatology with long-term trends.
!   2023-08-28  R. Hogan & R. Checa-Garcia Fixed Bilinear Interp. issue.
!   2024-11-04  R. Hogan  Locate and scale-down unrepresentative fires in climatology
! ------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

! ------------------------------------------------------------------
! Structure for storing aerosol climatology
TYPE :: TEAERC
  ! --- 2D climatology ---
  ! All 12 months are loaded at setup
  ! Vertically integrated mass (kg/m2) dimensioned lon,lat,mon,species
  REAL(KIND=JPRB), ALLOCATABLE :: COLUMN_MASS_MON(:,:,:,:)
  ! Aerosol scale height (m), dimensioned month,species
  REAL(KIND=JPRB), ALLOCATABLE :: SCALE_HEIGHT_MON(:,:)
  ! Vertically integrated mass for the current timestep, dimensioned
  ! lon,lat,species
  REAL(KIND=JPRB), ALLOCATABLE :: COLUMN_MASS(:,:,:)
  ! Cumulative mass fraction at half levels, i.e. the fraction of the
  ! column mass between top-of-atmosphere and a particular half level,
  ! dimensioned lev_model,species
  REAL(KIND=JPRB), ALLOCATABLE :: CUM_FRACTION(:,:)
  ! Keep track of aerosol type (e.g. DU, SU) because scale height is
  ! determined by this
  CHARACTER(2), ALLOCATABLE :: AER_TYPE(:)

  ! --- 4D climatology ---
  ! Only two months are ever held in memory
  ! Aerosol mass mixing ratio (kg/kg) dimensioned lon,lat,lev,species
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:) :: MMR_MON1
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:) :: MMR_MON2
  ! Corresponding pressure, dimensioned lon,lat,lev
  REAL(KIND=JPRB), POINTER, DIMENSION(:,:,:) :: PRESSURE_MON1
  REAL(KIND=JPRB), POINTER, DIMENSION(:,:,:) :: PRESSURE_MON2
  ! Aerosol mass mixing ratio (kg/kg) for the current timestep,
  ! dimensioned lon,lat,lev,species
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:) :: MMR
  ! Corresponding pressure, dimensioned lon,lat,lev
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:)   :: PRESSURE
  ! Some variables have long-term variability, represented by "epochs"
  ! (actually integer years)
  INTEGER(KIND=JPIM), ALLOCATABLE :: EPOCH(:)
  ! Year and month of the currently stored months, noting that the
  ! monthly values correspond ot the 15th of each month
  INTEGER(KIND=JPIM) :: YEAR1, YEAR2, MON1, MON2
  CHARACTER(LEN=512) :: FILE_NAME ! Save for future calls
  CHARACTER(LEN=128), ALLOCATABLE :: VAR_NAME(:)
  LOGICAL, ALLOCATABLE :: IS_SPECIES_4D(:) ! 4D or 3D
  ! 4D climatology also uses MMR and PRESSURE

  ! Climatologies are sometimes calculated from data including intense
  ! forest fires, and even if they are averaged over many years, the
  ! forest fire from one year can be present in the average. This
  ! tends not to be very important in the radiation scheme, but is
  ! apparent in the visibility which uses the concentration at the
  ! surface.  If REMOVE_FIRES==.TRUE. then we identify isolated
  ! high-concentration pixels in the Black Carbon and Organic Matter
  ! fields and work out a scaling that would remove them. Note that
  ! this is only used when the surface properties are requested in
  ! CALC_SURFACE used in the visibility calculation. CALC_PROFILES
  ! used by the radiation scheme does not use this.
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: SURFACE_SCALING
  
  ! --- General ---
  ! Index to output aerosol mass-mixing-ratio array of each
  ! climatological species (in case we need to leave gaps for
  ! prognostic and parametric species)
  INTEGER(KIND=JPIM), ALLOCATABLE :: ISPECIES_INDEX(:)

  ! Number of latitude, longitude, levels, months and aerosol species
  INTEGER(KIND=JPIM) :: NLAT = 0, NLON = 0, NLEV = 0, NMON = 0, NAER = 0

  ! Initial latitude and longitude of look-up table, and increments,
  ! all in radians
  REAL(KIND=JPRB) :: LAT1 = 0.0_JPRB, LON1 = 0.0_JPRB
  REAL(KIND=JPRB) :: DLAT = 0.0_JPRB, DLON = 0.0_JPRB

  LOGICAL :: IS_INITIALISED = .FALSE.
  LOGICAL :: REMOVE_FIRES   = .TRUE.
  INTEGER(KIND=JPIM) :: DIMENSIONALITY = 2

CONTAINS 

  PROCEDURE :: SETUP => SETUP
  PROCEDURE :: UPDATE
  PROCEDURE :: LOAD_MONTH
  PROCEDURE :: CALC_PROFILES
  PROCEDURE :: CALC_SURFACE
  PROCEDURE :: SU_CLIMATOLOGY2D, SU_CLIMATOLOGY4D
  PROCEDURE :: SET_SCALE_HEIGHT, SET_VARIABLE_SCALE_HEIGHT
  PROCEDURE :: LOCATE_FIRES
  PROCEDURE :: GET_INTERP_POINTS

END TYPE TEAERC

! Parameters
REAL(KIND=JPRB), PARAMETER :: RPRESSURE_SCALE_HEIGHT = 8434.0_JPRB


! ------------------------------------------------------------------
#include "abor1.intfb.h"

CONTAINS

! ------------------------------------------------------------------
! Load climatology from a 2D or 3D file (depending on LDAER3D) into
! "SELF", with the species given in the YAER_RAD_DESC structure
SUBROUTINE SETUP(SELF,LDAER3D,YAER_RAD_DESC,CDFILENAME)

  USE YOMHOOK,       ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOMLUN,        ONLY : NULOUT
  USE YOEAERRADDESC, ONLY : TAER_RAD_DESC

  ! Climatology object
  CLASS(TEAERC),       INTENT(INOUT) :: SELF
  ! Description of aerosols seen by radiation
  TYPE(TAER_RAD_DESC), INTENT(IN)    :: YAER_RAD_DESC
  ! True=3D, false=2D climatology
  LOGICAL,             INTENT(IN)    :: LDAER3D
  ! Name of the file to read; if prefixed by "/" treat as an absolute
  ! path, otherwise read from $DATA/ifsdata
  CHARACTER(LEN=*),    INTENT(IN)    :: CDFILENAME

  CHARACTER(LEN=512) :: CLFILENAME, CLDIRECTORY

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOEAERC:SETUP',0,ZHOOK_HANDLE)
  
  IF (CDFILENAME(1:1) /= "." .AND. CDFILENAME(1:1) /= "/") THEN
    ! Prefix file name by the data directory name from the "DATA"
    ! environment variable, if present
    CALL GET_ENVIRONMENT_VARIABLE("DATA",CLDIRECTORY)
    IF (CLDIRECTORY /= " ") THEN
      CLFILENAME = TRIM(CLDIRECTORY) // "/ifsdata/" // CDFILENAME
    ELSE
      CLFILENAME = TRIM(CDFILENAME)
    ENDIF
  ELSE
    ! File name starts with "." or "/": treat as an absolute path
    CLFILENAME = TRIM(CDFILENAME)
  ENDIF

  IF (LDAER3D) THEN
    WRITE(NULOUT,'(A)') 'YOEAERC: Initializing 3D/4D aerosol climatology'
    CALL SU_CLIMATOLOGY4D(SELF, YAER_RAD_DESC, CLFILENAME)
  ELSE
    WRITE(NULOUT,'(A)') 'YOEAERC: Loading 2D aerosol climatology'
    CALL SU_CLIMATOLOGY2D(SELF, YAER_RAD_DESC, CLFILENAME)
  ENDIF
  SELF%IS_INITIALISED = .TRUE.

  IF (LHOOK) CALL DR_HOOK('YOEAERC:SETUP',1,ZHOOK_HANDLE)

END SUBROUTINE SETUP


! ------------------------------------------------------------------
! Load 2D aerosol climatology from file
SUBROUTINE SU_CLIMATOLOGY2D(SELF, YAER_RAD_DESC, CDFILENAME)

  USE EASY_NETCDF_READ_MPI, ONLY : NETCDF_FILE
  USE YOMHOOK             , ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOEAERRADDESC       , ONLY : TAER_RAD_DESC, ICLIMATOLOGICAL
  USE YOMCST              , ONLY : RPI

  ! Arguments
  CLASS(TEAERC),       INTENT(INOUT) :: SELF
  TYPE(TAER_RAD_DESC), INTENT(IN)    :: YAER_RAD_DESC
  ! Name of the file to read; if prefixed by "/" treat as an absolute
  ! path, otherwise read from $DATA/ifsdata
  CHARACTER(LEN=*),    INTENT(IN)    :: CDFILENAME

  ! Local variables

  ! The NetCDF file containing the input climatology
  TYPE(NETCDF_FILE)  :: FILE
  
  ! Temporary storage of an aerosol species
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: ZTEMP_ARR

  ! Longitude and latitude
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:) :: ZLAT, ZLON

  INTEGER(KIND=JPIM) :: JAER  ! Loop index for aerosol type
  INTEGER(KIND=JPIM) :: ICLIM ! Index to current climatological aerosol type

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOEAERC:SU_CLIMATOLOGY2D',0,ZHOOK_HANDLE)
  
  CALL FILE%OPEN(TRIM(CDFILENAME), IVERBOSE=4)

  ICLIM = 0

  ! Read longitude and latitude
  CALL FILE%GET("lat", ZLAT)
  SELF%NLAT = SIZE(ZLAT)
  SELF%LAT1 = ZLAT(1) * RPI / 180.0_JPRB
  SELF%DLAT = (ZLAT(2)-ZLAT(1)) * RPI / 180.0_JPRB
  CALL FILE%GET("lon", ZLON)
  SELF%NLON = SIZE(ZLON)
  SELF%LON1 = ZLON(1) * RPI / 180.0_JPRB
  SELF%DLON = (ZLON(2)-ZLON(1)) * RPI / 180.0_JPRB

  ! Loop over all aerosols in the description structure (prognostic
  ! and diagnostic)
  DO JAER = 1,YAER_RAD_DESC%NAEROSOL
    IF (YAER_RAD_DESC%DATA(JAER)%NSOURCE == ICLIMATOLOGICAL) THEN
      ! Found a climatological aerosol species
      ICLIM = ICLIM + 1
      
      CALL FILE%GET(TRIM(YAER_RAD_DESC%DATA(JAER)%NAME), ZTEMP_ARR)

      IF (ICLIM == 1) THEN
        ! This is the first aerosol type loaded
        SELF%NMON = SIZE(ZTEMP_ARR,3)
        SELF%NAER = YAER_RAD_DESC%NCLIMATOLOGICAL

        ! Allocate main storage array for all the 3D climatological
        ! MACC quantities
        ALLOCATE(SELF%COLUMN_MASS_MON(SELF%NLON,SELF%NLAT,SELF%NMON,SELF%NAER))

        ! Set initial scale height to 2 km
        ALLOCATE(SELF%SCALE_HEIGHT_MON(SELF%NMON,SELF%NAER))
        SELF%SCALE_HEIGHT_MON = 2000.0_JPRB

        ! Space to store aerosol type
        ALLOCATE(SELF%AER_TYPE(SELF%NAER))

        ! Index for where to store each output
        ALLOCATE(SELF%ISPECIES_INDEX(SELF%NAER))
      ENDIF

      ! Store aerosol species
      SELF%COLUMN_MASS_MON(:,:,:,ICLIM) = ZTEMP_ARR

      ! Copy over aerosol type
      SELF%AER_TYPE(ICLIM) = YAER_RAD_DESC%DATA(JAER)%AER_TYPE

      SELF%ISPECIES_INDEX(ICLIM) = JAER
    ENDIF
  ENDDO

  CALL FILE%CLOSE()

  ! Allocate the instantaneous arrays
  ALLOCATE(SELF%COLUMN_MASS(SELF%NLON,SELF%NLAT,ICLIM))

  IF (LHOOK) CALL DR_HOOK('YOEAERC:SU_CLIMATOLOGY2D',1,ZHOOK_HANDLE)

END SUBROUTINE SU_CLIMATOLOGY2D


! ------------------------------------------------------------------
! Initialize 3D/4D aerosol climatology from file
SUBROUTINE SU_CLIMATOLOGY4D(SELF, YAER_RAD_DESC, CDFILENAME)

  USE YOMHOOK             , ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOMLUN              , ONLY : NULERR, NULOUT
  USE YOEAERRADDESC       , ONLY : TAER_RAD_DESC, ICLIMATOLOGICAL
  USE EASY_NETCDF_READ_MPI, ONLY : NETCDF_FILE
  USE YOMCST              , ONLY : RPI

  ! Arguments
  CLASS(TEAERC),       INTENT(INOUT) :: SELF
  TYPE(TAER_RAD_DESC), INTENT(IN)    :: YAER_RAD_DESC
  ! Name of the file to read; if prefixed by "/" treat as an absolute
  ! path, otherwise read from $DATA/ifsdata
  CHARACTER(LEN=*),    INTENT(IN)    :: CDFILENAME

  ! Local variables
  
  ! The NetCDF file containing the input climatology
  TYPE(NETCDF_FILE)  :: FILE

  ! Longitude and latitude
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:) :: ZLAT, ZLON, ZLEV

  INTEGER(KIND=JPIM) :: JAER  ! Loop index for aerosol type
  INTEGER(KIND=JPIM) :: ICLIM ! Index to current climatological aerosol type

  INTEGER(KIND=JPIM) :: NDIMS ! Number of dimensions in NetCDF variable

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


  IF (LHOOK) CALL DR_HOOK('YOEAERC:SU_CLIMATOLOGY4D',0,ZHOOK_HANDLE)

  ! Save for future reads
  SELF%FILE_NAME = CDFILENAME

  CALL FILE%OPEN(TRIM(CDFILENAME), IVERBOSE=4)
  
  ICLIM = 0

  ! Read longitude and latitude
  CALL FILE%GET("lat", ZLAT)
  SELF%NLAT = SIZE(ZLAT)
  SELF%LAT1 = ZLAT(1) * RPI / 180.0_JPRB
  SELF%DLAT = (ZLAT(2)-ZLAT(1)) * RPI / 180.0_JPRB
  CALL FILE%GET("lon", ZLON)
  SELF%NLON = SIZE(ZLON)
  SELF%LON1 = ZLON(1) * RPI / 180.0_JPRB
  SELF%DLON = (ZLON(2)-ZLON(1)) * RPI / 180.0_JPRB
  CALL FILE%GET("lev", ZLEV)
  SELF%NLEV = SIZE(ZLEV)
  SELF%NMON = 12
  
  IF (FILE%EXISTS("epoch")) THEN
    CALL FILE%GET("epoch", SELF%EPOCH)
  ENDIF

  SELF%NAER = YAER_RAD_DESC%NCLIMATOLOGICAL
  IF (SELF%NAER > 0) THEN
    ALLOCATE(SELF%IS_SPECIES_4D(SELF%NAER))
    SELF%IS_SPECIES_4D = .FALSE.
    
    ! Space to store aerosol type
    ALLOCATE(SELF%AER_TYPE(SELF%NAER))

    ! Index for where to store each output
    ALLOCATE(SELF%ISPECIES_INDEX(SELF%NAER))

    ! Storage for variable names
    ALLOCATE(SELF%VAR_NAME(SELF%NAER))

    ! Loop over all aerosols in the description structure (prognostic
    ! and diagnostic)
    DO JAER = 1,YAER_RAD_DESC%NAEROSOL
      IF (YAER_RAD_DESC%DATA(JAER)%NSOURCE == ICLIMATOLOGICAL) THEN
        ! Found a climatological aerosol species
        ICLIM = ICLIM + 1
      
        NDIMS = FILE%GET_RANK(TRIM(YAER_RAD_DESC%DATA(JAER)%NAME))
        SELF%VAR_NAME(ICLIM) = TRIM(YAER_RAD_DESC%DATA(JAER)%NAME)

        IF (NDIMS == 5) THEN
          SELF%IS_SPECIES_4D(ICLIM) = .TRUE.
          WRITE(NULOUT,'(A,I0,A,A,A)') '    Radiatively active aerosol ', JAER, ': ', &
               &  TRIM(SELF%VAR_NAME(ICLIM)), ' (4D climatological)'
        ELSEIF (NDIMS < 4) THEN
          WRITE(NULERR, '(A,A,A)') 'Error: ', TRIM(YAER_RAD_DESC%DATA(JAER)%NAME), &
               &  ' not found or of insufficient dimensions'
          CALL ABOR1('YOEAERC error')
        ELSE
          WRITE(NULOUT,'(A,I0,A,A,A)') '    Radiatively active aerosol ', JAER, ': ', &
               &  TRIM(SELF%VAR_NAME(ICLIM)), ' (3D climatological)'
        ENDIF

        ! Copy over aerosol type
        SELF%AER_TYPE(ICLIM) = YAER_RAD_DESC%DATA(JAER)%AER_TYPE

        SELF%ISPECIES_INDEX(ICLIM) = JAER
      ENDIF
    ENDDO
  
  ENDIF

  CALL FILE%CLOSE()
    
  ! Allocate the interpolation arrays
  ALLOCATE(SELF%MMR_MON1(SELF%NLON,SELF%NLAT,SELF%NLEV,SELF%NAER))
  ALLOCATE(SELF%MMR_MON2(SELF%NLON,SELF%NLAT,SELF%NLEV,SELF%NAER))
  ALLOCATE(SELF%PRESSURE_MON1(SELF%NLON,SELF%NLAT,SELF%NLEV))
  ALLOCATE(SELF%PRESSURE_MON2(SELF%NLON,SELF%NLAT,SELF%NLEV))

  ! Allocate arrays for an instant in time
  ALLOCATE(SELF%MMR(SELF%NLON,SELF%NLAT,SELF%NLEV,SELF%NAER))
  ALLOCATE(SELF%PRESSURE(SELF%NLON,SELF%NLAT,SELF%NLEV))

  IF (SELF%REMOVE_FIRES) THEN
    ALLOCATE(SELF%SURFACE_SCALING(SELF%NLON,SELF%NLAT,SELF%NAER))
    SELF%SURFACE_SCALING = 1.0_JPRB
  ENDIF
  
  SELF%DIMENSIONALITY = 4
  ! Currently no year/month is loaded
  SELF%YEAR1 = -1
  SELF%YEAR2 = -1
  SELF%MON1  = -1
  SELF%MON2  = -1
  
  IF (LHOOK) CALL DR_HOOK('YOEAERC:SU_CLIMATOLOGY4D',1,ZHOOK_HANDLE)

END SUBROUTINE SU_CLIMATOLOGY4D


! ------------------------------------------------------------------
! Set the scale height of all aerosols of type CDAER_TYPE to the
! constant value PSCALE_HEIGHT
SUBROUTINE SET_SCALE_HEIGHT(SELF, CDAER_TYPE, PSCALE_HEIGHT)

  USE YOMHOOK             , ONLY : LHOOK, DR_HOOK, JPHOOK

  ! Arguments
  CLASS(TEAERC),     INTENT(INOUT) :: SELF
  CHARACTER(2),      INTENT(IN)    :: CDAER_TYPE
  REAL(KIND=JPRB),   INTENT(IN)    :: PSCALE_HEIGHT

  ! Local variables
  INTEGER(KIND=JPIM) :: JA

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOEAERC:SET_SCALE_HEIGHT',0,ZHOOK_HANDLE)

  IF (SELF%DIMENSIONALITY == 2) THEN
    DO JA = 1,SELF%NAER
      IF (SELF%AER_TYPE(JA) == CDAER_TYPE) THEN
        SELF%SCALE_HEIGHT_MON(:,JA) = PSCALE_HEIGHT
      ENDIF
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOEAERC:SET_SCALE_HEIGHT',1,ZHOOK_HANDLE)

END SUBROUTINE SET_SCALE_HEIGHT


! ------------------------------------------------------------------
! Set the scale height of all aerosols of type CDAER_TYPE to
! PSCALE_HEIGHT, which varies with month
SUBROUTINE SET_VARIABLE_SCALE_HEIGHT(SELF, CDAER_TYPE, PSCALE_HEIGHT)

  USE YOMHOOK             , ONLY : LHOOK, DR_HOOK, JPHOOK

  ! Arguments
  CLASS(TEAERC),   INTENT(INOUT) :: SELF
  CHARACTER(2),    INTENT(IN)    :: CDAER_TYPE
  REAL(KIND=JPRB), INTENT(IN)    :: PSCALE_HEIGHT(:)

  ! Local variables
  INTEGER(KIND=JPIM) :: JA

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOEAERC:SET_VARIABLE_SCALE_HEIGHT',0,ZHOOK_HANDLE)

  IF (SELF%DIMENSIONALITY == 2) THEN
    DO JA = 1,SELF%NAER
      IF (SELF%AER_TYPE(JA) == CDAER_TYPE) THEN
        SELF%SCALE_HEIGHT_MON(:,JA) = PSCALE_HEIGHT(:)
      ENDIF
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOEAERC:SET_VARIABLE_SCALE_HEIGHT',1,ZHOOK_HANDLE)

END SUBROUTINE SET_VARIABLE_SCALE_HEIGHT


! ------------------------------------------------------------------
! Interpolate the monthly climatology to the current time as stored in
! the real number "PMONTH", where 1.0=January and non-integers
! indicate that linear interpolation in time is required. Values may
! lie between 0 and 13.  Typically the calling function interpolates
! the current time treating the 15th of each month as integer indices.
! The year is provided as an integer in case the input climatology
! also describes long-term trends.
SUBROUTINE UPDATE(SELF, KYEAR, PMONTH, P_ETA_HL)

  USE YOMLUN,  ONLY : NULOUT
  USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

  ! Arguments
  CLASS(TEAERC),      INTENT(INOUT) :: SELF
  INTEGER(KIND=JPIM), INTENT(IN)    :: KYEAR
  REAL(KIND=JPRB),    INTENT(IN)    :: PMONTH
  ! Eta (pressure/surface_pressure) at half levels of a reference
  ! atmosphere
  REAL(KIND=JPRB),   INTENT(IN), OPTIONAL :: P_ETA_HL(:)

  ! Local variables
  ! Indices to first and second month, and corresponding years, in interpolation
  INTEGER(KIND=JPIM) :: IM1, IM2, IY1, IY2
  REAL(KIND=JPRB)    :: ZWEIGHT
  REAL(KIND=JPRB)    :: SCALE_HEIGHT(SELF%NAER)

  INTEGER(KIND=JPIM) :: JA

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOEAERC:UPDATE',0,ZHOOK_HANDLE)

  ! Weight is the decimal part of PMONTH
  ZWEIGHT = PMONTH - INT(PMONTH)

  ! Fold indices into the range 1-12
  IM1 = MOD(INT(PMONTH)-1, 12)+1
  IM2 = MOD(INT(PMONTH),   12)+1

  IY1 = KYEAR
  IY2 = KYEAR
  IF (PMONTH < 1.0_JPRB) THEN
    IY1 = KYEAR-1
  ELSEIF (PMONTH >= 12.0_JPRB) THEN
    IY2 = KYEAR+1
  ENDIF

  IF (SELF%DIMENSIONALITY >= 3) THEN ! 3D/4D climatology

    WRITE(NULOUT,'(A,I4,A,F10.7,A)') 'YOEAERC: Updating 3D/4D aerosol climatology to year=', &
         &  KYEAR, ', month=', PMONTH, ' (Jan15=1.0)'
    
    IF (SELF%MON1 /= IM1 .OR. SELF%YEAR1 /= IY1) THEN
      ! First month in interpolation needs updating - first see if we
      ! can copy over second month
      IF (SELF%MON2 == IM1 .AND. SELF%YEAR2 == IY1) THEN
        WRITE(NULOUT,'(A,I4,A,I2,A)') 'YOEAERC: Copying aerosol data for year=', IY1, ', month=', IM1, &
             & ' from second to first interpolation point'
        SELF%MMR_MON1 = SELF%MMR_MON2
        SELF%PRESSURE_MON1 = SELF%PRESSURE_MON2
        SELF%YEAR1 = IY1
        SELF%MON1  = IM1
      ELSE
        CALL SELF%LOAD_MONTH(IY1, IM1, SELF%MMR_MON1, SELF%PRESSURE_MON1)
        SELF%YEAR1 = IY1
        SELF%MON1  = IM1
      ENDIF
    ENDIF
    IF (SELF%MON2 /= IM2 .OR. SELF%YEAR2 /= IY2) THEN
      ! Second month in interpolation needs updating
      CALL SELF%LOAD_MONTH(IY2, IM2, SELF%MMR_MON2, SELF%PRESSURE_MON2)
      SELF%YEAR2 = IY2
      SELF%MON2  = IM2
    ENDIF

    SELF%PRESSURE = (1.0_JPRB - ZWEIGHT) * SELF%PRESSURE_MON1 &
         &        +             ZWEIGHT  * SELF%PRESSURE_MON2
    SELF%MMR = (1.0_JPRB - ZWEIGHT) * SELF%MMR_MON1 &
         &    +            ZWEIGHT  * SELF%MMR_MON2

    IF (SELF%REMOVE_FIRES) THEN
      CALL SELF%LOCATE_FIRES
    ENDIF
    
  ELSE ! 2D climatology

    WRITE(NULOUT,'(A,I4,A,F10.7,A)') 'YOEAERC: Updating 2D aerosol climatology to year=', &
         &  KYEAR, ', month=', PMONTH, ' (Jan15=1.0)'
    SELF%COLUMN_MASS = (1.0_JPRB - ZWEIGHT) * SELF%COLUMN_MASS_MON(:,:,IM1,:) &
         &           +             ZWEIGHT  * SELF%COLUMN_MASS_MON(:,:,IM2,:)
    IF (.NOT. ALLOCATED(SELF%CUM_FRACTION)) THEN
      ALLOCATE(SELF%CUM_FRACTION(SIZE(P_ETA_HL), SELF%NAER))
    ENDIF
    SCALE_HEIGHT = (1.0_JPRB - ZWEIGHT) * SELF%SCALE_HEIGHT_MON(IM1,:) &
         &       +             ZWEIGHT  * SELF%SCALE_HEIGHT_MON(IM2,:)
    DO JA = 1,SELF%NAER
      SELF%CUM_FRACTION(:,JA) = P_ETA_HL ** (RPRESSURE_SCALE_HEIGHT/SCALE_HEIGHT(JA))
    ENDDO

  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOEAERC:UPDATE',1,ZHOOK_HANDLE)

END SUBROUTINE UPDATE

! ------------------------------------------------------------------
! Load a month snapshot of 3D aerosol climatology from file
SUBROUTINE LOAD_MONTH(SELF, KYEAR, KMON, PMMR_MON, PPRESSURE_MON)

  USE YOMLUN,  ONLY : NULOUT
  USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
  USE EASY_NETCDF_READ_MPI, ONLY : NETCDF_FILE

  ! Arguments
  CLASS(TEAERC),      INTENT(INOUT) :: SELF
  INTEGER(KIND=JPIM), INTENT(IN)    :: KYEAR, KMON

  ! Mass mixing ratio, dimensioned lon,lat,lev,species
  REAL(KIND=JPRB),    INTENT(OUT)   :: PMMR_MON(:,:,:,:)
  ! Pressure (Pa), dimensioned lon,lat,lev
  REAL(KIND=JPRB),    INTENT(OUT)   :: PPRESSURE_MON(:,:,:)

  ! Local variables
  TYPE(NETCDF_FILE)  :: FILE

  REAL(KIND=JPRB)    :: ZWEIGHT

  ! Indices to first and second epoch for interpolation (IY2=0 if no
  ! interpolation)
  INTEGER(KIND=JPIM) :: IY1, IY2

  REAL(KIND=JPRB), ALLOCATABLE :: ZTEMP_ARR1(:,:,:), ZTEMP_ARR2(:,:,:)

  ! Loop indices for epoch and aerosol type
  INTEGER(KIND=JPIM) :: JY, JA

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOEAERC:LOAD_MONTH',0,ZHOOK_HANDLE)

  CALL FILE%OPEN(TRIM(SELF%FILE_NAME), IVERBOSE=4)

  ! For 4D variables we need to find the year indices and weight
  IF (ALLOCATED(SELF%EPOCH)) THEN
    IY1 = 0
    IY2 = 0
    ZWEIGHT = 1.0_JPRB
    DO JY = 1,SIZE(SELF%EPOCH)-1
      IF (KYEAR <= SELF%EPOCH(JY)) THEN
        ! Found an exact year: "<=" rather than "==" because if
        ! requested year is before the first one in the file we clamp
        ! at the first year
        IY1 = JY
        EXIT
      ELSEIF (KYEAR < SELF%EPOCH(JY+1)) THEN
        ! Year falls between two in the file: interpolate
        IY1 = JY
        IY2 = JY+1
        ZWEIGHT = (REAL(SELF%EPOCH(JY+1),JPRB) - REAL(KYEAR,JPRB)) &
             &   /(REAL(SELF%EPOCH(JY+1),JPRB) - REAL(SELF%EPOCH(JY),JPRB))
        EXIT
      ENDIF
    ENDDO
    IF (IY1 == 0) THEN
      ! Year not yet found so must be the last one 
      IY1 = SIZE(SELF%EPOCH)
    ENDIF
  ENDIF

  DO JA = 1,SELF%NAER
    IF (SELF%IS_SPECIES_4D(JA)) THEN
      IF (IY2 == 0) THEN
        ! Exact year required
        CALL FILE%GET(TRIM(SELF%VAR_NAME(JA)), ZTEMP_ARR1, KMON, IY1)
        PMMR_MON(:,:,:,JA) = ZTEMP_ARR1
      ELSE
        ! Interpolate between two years
        CALL FILE%GET(TRIM(SELF%VAR_NAME(JA)), ZTEMP_ARR1, KMON, IY1)
        CALL FILE%GET(TRIM(SELF%VAR_NAME(JA)), ZTEMP_ARR2, KMON, IY2)
        PMMR_MON(:,:,:,JA) = ZWEIGHT  * ZTEMP_ARR1 &
             &   + (1.0_JPRB-ZWEIGHT) * ZTEMP_ARR2
      ENDIF
    ELSE
      ! 3D field
      CALL FILE%GET(TRIM(SELF%VAR_NAME(JA)), ZTEMP_ARR1, KMON)
      PMMR_MON(:,:,:,JA) = ZTEMP_ARR1
    ENDIF
  ENDDO

  CALL FILE%GET("pressure", ZTEMP_ARR1, KMON)
  PPRESSURE_MON = ZTEMP_ARR1

  CALL FILE%CLOSE()

  IF (LHOOK) CALL DR_HOOK('YOEAERC:LOAD_MONTH',1,ZHOOK_HANDLE)

END SUBROUTINE LOAD_MONTH


! ------------------------------------------------------------------
! Calculate profiles of interpolated climatological aerosol mass
! mixing ratio at the requested latitude and longitude
SUBROUTINE CALC_PROFILES(SELF, KSTART, KEND, KCOL, KLEV, PSINLAT, PLON, &
     &                   PPRESSUREH, PAERO_MMR, LDLEAVESPACE)

  USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOMCST,  ONLY : RG
  USE YOMLUN,  ONLY : NULERR
  
  ! Arguments
  CLASS(TEAERC), INTENT(IN) :: SELF
  ! Indices of start and end columns
  INTEGER(KIND=JPIM), INTENT(IN) :: KSTART, KEND
  ! Number of columns
  INTEGER(KIND=JPIM), INTENT(IN) :: KCOL
  ! Number of levels
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  ! Sine of latitude of input columns
  REAL(KIND=JPRB), INTENT(IN) :: PSINLAT(KCOL)
  ! Longitude of input columns in radians
  REAL(KIND=JPRB), INTENT(IN) :: PLON(KCOL)
  ! Half-level pressure
  REAL(KIND=JPRB), INTENT(IN) :: PPRESSUREH(KCOL,KLEV+1)
  ! Output aerosol mass mixing ratio of each species
  REAL(KIND=JPRB), INTENT(INOUT) :: PAERO_MMR(:,:,:)
  ! Do we place climatological aerosol mixing ratios contiguously
  ! (LDLEAVESPACE=.false.) or with space for other types of aerosol
  ! (LDLEAVESPACE=.true.)?
  LOGICAL, OPTIONAL, INTENT(IN) :: LDLEAVESPACE

  ! Local variables
  REAL(KIND=JPRB)    :: ZLATWEIGHT(KSTART:KEND), ZLONWEIGHT(KSTART:KEND)
  INTEGER(KIND=JPIM) :: ILON1(KSTART:KEND), ILON2(KSTART:KEND)
  INTEGER(KIND=JPIM) :: ILAT1(KSTART:KEND), ILAT2(KSTART:KEND)

  ! Interpolation of climatology to model longitude/latitude grid:
  ! mass mixing ratio in kg/kg, plus extra levels at TOA and at very
  ! high pressure to enact clamped extrapolation
  REAL(KIND=JPRB) :: ZAER3D(KSTART:KEND,0:SELF%NLEV+1)
  ! Column-integrated mass in kg/m2
  REAL(KIND=JPRB) :: ZAER2D(KSTART:KEND)
  ! Values in interpolation, in kg/kg
  REAL(KIND=JPRB) :: ZAERA, ZAERB
  ! Pressure on climatology vertical grid (Pa), plus extra levels at
  ! TOA (pressure=0) and at very high pressure to enact clamped
  ! extrapolation
  REAL(KIND=JPRB) :: ZREF_PRS(KSTART:KEND,0:SELF%NLEV+1)
  ! Inverse of mass of air in a layer (m2/kg)
  REAL(KIND=JPRB) :: ZFACTOR(KSTART:KEND,KLEV)

  ! Weight of aerosol in vertical interpolation
  REAL(KIND=JPRB) :: ZWEIGHT

  ! Full-level pressure on model grid (Pa)
  REAL(KIND=JPRB) :: ZPRESSUREF

  ! Index to the output aerosol array
  INTEGER(KIND=JPIM) :: IINDEX(SELF%NAER)

  INTEGER(KIND=JPIM) :: JK, JL, JA, JC, IC

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


  IF (LHOOK) CALL DR_HOOK('YOEAERC:CALC_PROFILES',0,ZHOOK_HANDLE)

  IINDEX = [(JA,JA=1,SELF%NAER)]

  IF (PRESENT(LDLEAVESPACE)) THEN
    IF (LDLEAVESPACE) THEN
      IINDEX = SELF%ISPECIES_INDEX
    ENDIF
  ENDIF

  IF (SIZE(PAERO_MMR,3) < MAXVAL(IINDEX)) THEN
    WRITE(NULERR,'(A,I0,A,I0)') 'Error: PAERO_MMR provided to TEAERC::CALC_PROFILES holds up to ', &
         &  SIZE(PAERO_MMR,3), ' aerosol species, but maximum climatological aerosol index is ', &
         &  MAXVAL(IINDEX)
    CALL ABOR1('YOEAERC error')
  ENDIF
  
  CALL SELF%GET_INTERP_POINTS(KSTART,KEND,PSINLAT(KSTART:KEND), PLON(KSTART:KEND), &
       &                      ZLATWEIGHT, ZLONWEIGHT, ILAT1, ILAT2, ILON1, ILON2)

  IF (SELF%DIMENSIONALITY == 4) THEN ! 3D or 4D climatology
    
    ! Bilinear interpolation of pressure from the climatology to the
    ! lon/lat grid of the model
    ZREF_PRS(KSTART:KEND,0) = 0.0_JPRB
    ZREF_PRS(KSTART:KEND,SELF%NLEV+1) = 1.0E7_JPRB
    DO JK = 1,SELF%NLEV
      DO JL = KSTART,KEND
        ZAERA = SELF%PRESSURE(ILON1(JL),ILAT1(JL),JK) &
             & +ZLATWEIGHT(JL)*( SELF%PRESSURE(ILON1(JL),ILAT2(JL),JK) &
             &                  -SELF%PRESSURE(ILON1(JL),ILAT1(JL),JK))
        ZAERB = SELF%PRESSURE(ILON2(JL),ILAT1(JL),JK) &
             & +ZLATWEIGHT(JL)*( SELF%PRESSURE(ILON2(JL),ILAT2(JL),JK) &
             &                  -SELF%PRESSURE(ILON2(JL),ILAT1(JL),JK))
        ZREF_PRS(JL,JK) = ZAERA + ZLONWEIGHT(JL) * (ZAERB-ZAERA)
      ENDDO
    ENDDO
    ! Do each aerosol species in turn (most efficient in memory access
    ! but possibly not in computation)
    DO JA = 1,SELF%NAER
      ! Bilinear interpolation of aerosol from the climatology to the
      ! lon/lat grid of the model
      ZAER3D(KSTART:KEND,0) = 0.0_JPRB ! Zero at TOA
      DO JK = 1,SELF%NLEV
        DO JL = KSTART,KEND
          ZAERA = SELF%MMR(ILON1(JL),ILAT1(JL),JK,JA) &
               & +ZLATWEIGHT(JL)*( SELF%MMR(ILON1(JL),ILAT2(JL),JK,JA) &
               &                  -SELF%MMR(ILON1(JL),ILAT1(JL),JK,JA))
          ZAERB = SELF%MMR(ILON2(JL),ILAT1(JL),JK,JA) &
               & +ZLATWEIGHT(JL)*( SELF%MMR(ILON2(JL),ILAT2(JL),JK,JA) &
               &                  -SELF%MMR(ILON2(JL),ILAT1(JL),JK,JA))
          ZAER3D(JL,JK) = ZAERA + ZLONWEIGHT(JL) * (ZAERB-ZAERA) 
        ENDDO
      ENDDO
      ! Constant below lowest point in climatology profile
      ZAER3D(KSTART:KEND,SELF%NLEV+1) = ZAER3D(KSTART:KEND,SELF%NLEV) 

      ! Loop over model columns
      DO JL = KSTART,KEND
        ! Linear interpolation in pressure
        IC=0 ! Index to climatology profile
        DO JK = 1,KLEV
          ! Compute full level pressure
          ZPRESSUREF = 0.5_JPRB * (PPRESSUREH(JL,JK)+PPRESSUREH(JL,JK+1))
          ! Seek IC such that ZPRESSUREF is bracketed by climatology
          ! pressures ZREF_PRS(JL,IC) and ZREF_PRS(JL,IC+1)
          DO WHILE (ZREF_PRS(JL,IC+1) < ZPRESSUREF)
            IC = IC+1
          ENDDO
          ZWEIGHT = (ZREF_PRS(JL,IC+1)-ZPRESSUREF)/(ZREF_PRS(JL,IC+1)-ZREF_PRS(JL,IC))
          PAERO_MMR(JL,JK,IINDEX(JA)) = ZWEIGHT * ZAER3D(JL,IC) &
               &           + (1.0_JPRB-ZWEIGHT) * ZAER3D(JL,IC+1)
        ENDDO
      ENDDO
    ENDDO

  ELSE ! 2D climatology

    IF (SIZE(SELF%CUM_FRACTION,1) /= KLEV+1) THEN
      WRITE(NULERR,'(A,I0,A,I0,A)') 'Error: 2D aerosol climatology has cumulative fraction CUM_FRACTION defined on ', &
           &  SIZE(SELF%CUM_FRACTION,1), ' half levels but CALC_PROFILES called with ', KLEV+1, &
           &  ' half levels'
      CALL ABOR1('YOEAERC error')
    ENDIF
    
    ! Factor to convert kg/m2 in a layer to mean mixing ratio in kg/kg
    ZFACTOR = RG / (PPRESSUREH(KSTART:KEND,2:KLEV+1) - PPRESSUREH(KSTART:KEND,1:KLEV))

    DO JA = 1,SELF%NAER
      DO JL = KSTART,KEND
        ! Bilinear interpolation in lon/lat
        ZAERA = SELF%COLUMN_MASS(ILON1(JL),ILAT1(JL),JA) &
             & +ZLATWEIGHT(JL)*( SELF%COLUMN_MASS(ILON1(JL),ILAT2(JL),JA) &
             &                  -SELF%COLUMN_MASS(ILON1(JL),ILAT1(JL),JA))
        ZAERB = SELF%COLUMN_MASS(ILON2(JL),ILAT1(JL),JA) &
             & +ZLATWEIGHT(JL)*( SELF%COLUMN_MASS(ILON2(JL),ILAT2(JL),JA) &
             &                  -SELF%COLUMN_MASS(ILON2(JL),ILAT1(JL),JA))
        ZAER2D(JL) = ZAERA + ZLONWEIGHT(JL) * (ZAERB-ZAERA)
      ENDDO
      DO JK = 1,KLEV
        PAERO_MMR(KSTART:KEND,JK,IINDEX(JA)) = ZAER2D(KSTART:KEND) * ZFACTOR(KSTART:KEND,JK) &
             &  * (SELF%CUM_FRACTION(JK+1,JA) - SELF%CUM_FRACTION(JK,JA))
      ENDDO
    ENDDO

  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOEAERC:CALC_PROFILES',1,ZHOOK_HANDLE)

END SUBROUTINE CALC_PROFILES


! ------------------------------------------------------------------
! Calculate surface values interpolated climatological aerosol mass
! mixing ratio at the requested latitude and longitude
SUBROUTINE CALC_SURFACE(SELF, KSTART, KEND, KCOL, KLEV, PSINLAT, PLON, &
     &                  PPRESSUREH, PAERO_MMR, LDLEAVESPACE)

  USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOMCST,  ONLY : RG
  USE YOMLUN,  ONLY : NULERR

  ! Arguments
  CLASS(TEAERC), INTENT(IN) :: SELF
  ! Indices of start and end columns
  INTEGER(KIND=JPIM), INTENT(IN) :: KSTART, KEND
  ! Number of columns
  INTEGER(KIND=JPIM), INTENT(IN) :: KCOL
  ! Number of levels
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  ! Sine of latitude of input columns
  REAL(KIND=JPRB), INTENT(IN) :: PSINLAT(KCOL)
  ! Longitude of input columns in radians
  REAL(KIND=JPRB), INTENT(IN) :: PLON(KCOL)
  ! Half-level pressure (Pa)
  REAL(KIND=JPRB), INTENT(IN) :: PPRESSUREH(KCOL,KLEV+1)
  ! Output aerosol mass mixing ratio of each species
  REAL(KIND=JPRB), INTENT(OUT) :: PAERO_MMR(KCOL,SELF%NAER)
  ! Do we place climatological aerosol mixing ratios contiguously
  ! (LDLEAVESPACE=.false.) or with space for other types of aerosol
  ! (LDLEAVESPACE=.true.)?
  LOGICAL, OPTIONAL, INTENT(IN) :: LDLEAVESPACE

  ! Local variables
  REAL(KIND=JPRB)    :: ZLATWEIGHT(KSTART:KEND), ZLONWEIGHT(KSTART:KEND)
  INTEGER(KIND=JPIM) :: ILON1(KSTART:KEND), ILON2(KSTART:KEND)
  INTEGER(KIND=JPIM) :: ILAT1(KSTART:KEND), ILAT2(KSTART:KEND)

  ! Interpolation of climatology to model longitude/latitude grid:
  ! mass in layer in kg/m2
  REAL(KIND=JPRB) :: ZAER
  ! Column-integrated mass in kg/m2
  REAL(KIND=JPRB) :: ZAER2D(KSTART:KEND)
  ! Values in interpolation, in kg/m2
  REAL(KIND=JPRB) :: ZAERA, ZAERB
  ! Pressure on climatology vertical grid (Pa), plus TOA and a very
  ! high pressure to ensure the surface is bracketed
  REAL(KIND=JPRB) :: ZREF_PRS(KSTART:KEND,0:SELF%NLEV+1)
  ! Inverse of mass of air in lowest layer (m2/kg)
  REAL(KIND=JPRB) :: ZFACTOR(KSTART:KEND)

  ! Weight of aerosol in vertical interpolation
  REAL(KIND=JPRB) :: ZWEIGHT

  ! Index to the output aerosol array
  INTEGER(KIND=JPIM) :: IINDEX(SELF%NAER)

  INTEGER(KIND=JPIM) :: JK, JL, JA, IC
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


  IF (LHOOK) CALL DR_HOOK('YOEAERC:CALC_SURFACE',0,ZHOOK_HANDLE)

  IINDEX = [(JA,JA=1,SELF%NAER)]

  IF (PRESENT(LDLEAVESPACE)) THEN
    IF (LDLEAVESPACE) THEN
      IINDEX = SELF%ISPECIES_INDEX
    ENDIF
  ENDIF

  CALL SELF%GET_INTERP_POINTS(KSTART,KEND,PSINLAT(KSTART:KEND), PLON(KSTART:KEND), &
       &                      ZLATWEIGHT, ZLONWEIGHT, ILAT1, ILAT2, ILON1, ILON2)

  IF (SELF%DIMENSIONALITY >= 3) THEN ! 3D or 4D climatology

    ! Interpolate reference pressure profile to lon/lat grid of
    ! columns
    ZREF_PRS(KSTART:KEND,0) = 0.0_JPRB
    ZREF_PRS(KSTART:KEND,SELF%NLEV+1) = 1.0E7_JPRB
    DO JK = 1,SELF%NLEV
      DO JL = KSTART,KEND
        ZAERA = SELF%PRESSURE(ILON1(JL),ILAT1(JL),JK) &
             & +ZLATWEIGHT(JL)*( SELF%PRESSURE(ILON1(JL),ILAT2(JL),JK) &
             &                  -SELF%PRESSURE(ILON1(JL),ILAT1(JL),JK))
        ZAERB = SELF%PRESSURE(ILON2(JL),ILAT1(JL),JK) &
             & +ZLATWEIGHT(JL)*( SELF%PRESSURE(ILON2(JL),ILAT2(JL),JK) &
             &                  -SELF%PRESSURE(ILON2(JL),ILAT1(JL),JK))
        ZREF_PRS(JL,JK) = ZAERA + ZLONWEIGHT(JL) * (ZAERB-ZAERA)
      ENDDO
    ENDDO

    ! Interpolate pressure difference across lowest layer of
    ! climatology
    DO JL = KSTART,KEND
      ! Bracket the surface pressure to lie between ZREF_PRS(JL,IC)
      ! and ZREF_PRS(JL,IC+1). Note that MMR's has only NLEV points in
      ! pressure, so if the requested surface pressure is higher (in
      ! Pa, i.e. lower in height) than the lowest (in height)
      ! climatological value, we want to get IC=NLEV-1 so that we
      ! interpolate between NLEV-1 and NLEV, and also a value of
      ! ZWEIGHT=0 so that we clamp to the lowest (in height)
      ! climatological value.
      IC = SELF%NLEV - 1
      DO WHILE (PPRESSUREH(JL,KLEV+1) < ZREF_PRS(JL,IC))
        IC = IC - 1
      ENDDO
      ZWEIGHT = (ZREF_PRS(JL,IC+1)-PPRESSUREH(JL,KLEV+1)) &
           &  / (ZREF_PRS(JL,IC+1)-ZREF_PRS(JL,IC))
      ! Prevent linear extrapolation
      ZWEIGHT = MAX(0.0_JPRB, MIN(ZWEIGHT, 1.0_JPRB))
      
      IF (SELF%REMOVE_FIRES) THEN
        DO JA = 1,SELF%NAER
          ! Interpolate aerosol in height at latitude point 1
          ZAERA = SELF%SURFACE_SCALING(ILON1(JL),ILAT1(JL),JA) &
               &  * (       SELF%MMR(ILON1(JL),ILAT1(JL),IC+1,JA) &
               &  +ZWEIGHT*(SELF%MMR(ILON1(JL),ILAT1(JL),IC,JA)-SELF%MMR(ILON1(JL),ILAT1(JL),IC+1,JA)))
          ZAERB = SELF%SURFACE_SCALING(ILON2(JL),ILAT1(JL),JA) &
               &  * (       SELF%MMR(ILON2(JL),ILAT1(JL),IC+1,JA) &
               &  +ZWEIGHT*(SELF%MMR(ILON2(JL),ILAT1(JL),IC,JA)-SELF%MMR(ILON2(JL),ILAT1(JL),IC+1,JA)))
          ! And interpolate in longitude
          ZAER = ZAERA + ZLONWEIGHT(JL)*(ZAERB-ZAERA)
          ! Interpolate aerosol in height at latitude point 2
          ZAERA = SELF%SURFACE_SCALING(ILON1(JL),ILAT2(JL),JA) &
               &  * (       SELF%MMR(ILON1(JL),ILAT2(JL),IC+1,JA) &
               &  +ZWEIGHT*(SELF%MMR(ILON1(JL),ILAT2(JL),IC,JA)-SELF%MMR(ILON1(JL),ILAT2(JL),IC+1,JA)))
          ZAERB = SELF%SURFACE_SCALING(ILON2(JL),ILAT2(JL),JA) &
               &  * (       SELF%MMR(ILON2(JL),ILAT2(JL),IC+1,JA) &
               &  +ZWEIGHT*(SELF%MMR(ILON2(JL),ILAT2(JL),IC,JA)-SELF%MMR(ILON2(JL),ILAT2(JL),IC+1,JA)))
          ! Final interpolation in longitude and latitude
          PAERO_MMR(JL,IINDEX(JA)) = ZAER &
               &  + ZLATWEIGHT(JL) * (ZAERA + ZLONWEIGHT(JL)*(ZAERB-ZAERA) - ZAER)
        ENDDO
      ELSE
        DO JA = 1,SELF%NAER
          ! Interpolate aerosol in lon/lat at climatology level IC
          ZAERA = SELF%MMR(ILON1(JL),ILAT1(JL),IC,JA) &
               & +ZLATWEIGHT(JL)*( SELF%MMR(ILON1(JL),ILAT2(JL),IC,JA) &
               &                  -SELF%MMR(ILON1(JL),ILAT1(JL),IC,JA))
          ZAERB = SELF%MMR(ILON2(JL),ILAT1(JL),IC,JA) &
               & +ZLATWEIGHT(JL)*( SELF%MMR(ILON2(JL),ILAT2(JL),IC,JA) &
               &                  -SELF%MMR(ILON2(JL),ILAT1(JL),IC,JA))
          ! Weight aerosol concentration from level IC
          ZAER = ZWEIGHT * (ZAERA + ZLONWEIGHT(JL) * (ZAERB-ZAERA))
          ! Interpolate aerosol in lon/lat at climatology level IC+1
          ZAERA = SELF%MMR(ILON1(JL),ILAT1(JL),IC+1,JA) &
               & +ZLATWEIGHT(JL)*( SELF%MMR(ILON1(JL),ILAT2(JL),IC+1,JA) &
               &                  -SELF%MMR(ILON1(JL),ILAT1(JL),IC+1,JA))
          ZAERB = SELF%MMR(ILON2(JL),ILAT1(JL),IC+1,JA) &
               & +ZLATWEIGHT(JL)*( SELF%MMR(ILON2(JL),ILAT2(JL),IC+1,JA) &
               &                  -SELF%MMR(ILON2(JL),ILAT1(JL),IC+1,JA))
          ! Final mass mixing ratio is the weighted average of the
          ! values at IC and IC+1
          PAERO_MMR(JL,IINDEX(JA)) = ZAER &
               &  + (1.0_JPRB-ZWEIGHT) * (ZAERA + ZLONWEIGHT(JL) * (ZAERB-ZAERA))
        ENDDO
      ENDIF
      
    ENDDO

  ELSE ! 2D climatology

    IF (SIZE(SELF%CUM_FRACTION,1) /= KLEV+1) THEN
      WRITE(NULERR,'(A,I0,A,I0,A)') 'Error: 2D aerosol climatology has cumulative fraction CUM_FRACTION defined on ', &
           &  SIZE(SELF%CUM_FRACTION,1), ' half levels but CALC_SURFACE called with ', KLEV+1, &
           &  ' half levels'
      CALL ABOR1('YOEAERC error')
    ENDIF
    
    ! Factor to convert kg/m2 in a layer to mean mixing ratio in kg/kg
    ZFACTOR = RG / (PPRESSUREH(KSTART:KEND,KLEV+1) - PPRESSUREH(KSTART:KEND,KLEV))
    
    DO JA = 1,SELF%NAER
      DO JL = KSTART,KEND
        ! Bilinear interpolation in lon/lat
        ZAERA = SELF%COLUMN_MASS(ILON1(JL),ILAT1(JL),JA) &
             & +ZLATWEIGHT(JL)*( SELF%COLUMN_MASS(ILON1(JL),ILAT2(JL),JA) &
             &                  -SELF%COLUMN_MASS(ILON1(JL),ILAT1(JL),JA))
        ZAERB = SELF%COLUMN_MASS(ILON2(JL),ILAT1(JL),JA) &
             & +ZLATWEIGHT(JL)*( SELF%COLUMN_MASS(ILON2(JL),ILAT2(JL),JA) &
             &                  -SELF%COLUMN_MASS(ILON2(JL),ILAT1(JL),JA))
        ZAER2D(JL) = ZAERA + ZLONWEIGHT(JL) * (ZAERB-ZAERA) 
      ENDDO
      DO JK = 1,KLEV
        PAERO_MMR(KSTART:KEND,IINDEX(JA)) = ZAER2D(KSTART:KEND) * ZFACTOR(KSTART:KEND) &
             &  * (SELF%CUM_FRACTION(KLEV+1,JA) - SELF%CUM_FRACTION(KLEV,JA))
      ENDDO
    ENDDO

  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOEAERC:CALC_SURFACE',1,ZHOOK_HANDLE)

END SUBROUTINE CALC_SURFACE


! ------------------------------------------------------------------
! Calculate interpolation points and weights in latitude and longitude
SUBROUTINE GET_INTERP_POINTS(SELF, KSTART, KEND, PSINLAT, PLON, &
     &                       PLATWEIGHT, PLONWEIGHT, &
     &                       KLAT1, KLAT2, KLON1, KLON2)

  USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

  ! Arguments
  CLASS(TEAERC), INTENT(IN) :: SELF
  ! Indices of start and end columns
  INTEGER(KIND=JPIM), INTENT(IN) :: KSTART, KEND
  ! Sine of latitude of input columns
  REAL(KIND=JPRB), INTENT(IN) :: PSINLAT(KSTART:KEND)
  ! Longitude of input columns in radians
  REAL(KIND=JPRB), INTENT(IN) :: PLON(KSTART:KEND)
  ! Weighting (0-1) between second rather than first point of interpolation
  REAL(KIND=JPRB), INTENT(OUT), DIMENSION(KSTART:KEND) :: PLATWEIGHT, PLONWEIGHT
  ! First and second points of (linear) latitude and longitude
  ! interpolation
  INTEGER(KIND=JPIM), INTENT(OUT), DIMENSION(KSTART:KEND) :: KLAT1, KLAT2, KLON1, KLON2

  ! Local variables
  INTEGER(KIND=JPIM) :: JL
  REAL(KIND=JPRB) :: ZLAT(KSTART:KEND)

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOEAERC:GET_INTERP_POINTS',0,ZHOOK_HANDLE)

  ! Longitude interpolation points and weights
  PLONWEIGHT = 1.0_JPRB + (PLON - SELF%LON1) / SELF%DLON
  KLON1 = MAX(1, MIN(INT(PLONWEIGHT), SELF%NLON))
  KLON2 = MOD(KLON1,SELF%NLON)+1 ! Correct treatment of dateline
  PLONWEIGHT = PLONWEIGHT - KLON1

  ZLAT = ASIN(PSINLAT)
  PLATWEIGHT = 1.0_JPRB + (ZLAT - SELF%LAT1) / SELF%DLAT
  KLAT1 = MAX(0, MIN(INT(PLATWEIGHT), SELF%NLAT-1))
  KLAT2 = KLAT1+1
  ! Linear interpolation in latitude
  !PLATWEIGHT = PLATWEIGHT - KLAT1
  ! Interpolation in sin(latitude)
  PLATWEIGHT = (PSINLAT - SIN(SELF%LAT1+(KLAT1-1)*SELF%DLAT)) &
       &     / (SIN(SELF%LAT1+(KLAT2-1)*SELF%DLAT) - SIN(SELF%LAT1+(KLAT1-1)*SELF%DLAT))

  IF (LHOOK) CALL DR_HOOK('YOEAERC:GET_INTERP_POINTS',1,ZHOOK_HANDLE)

END SUBROUTINE GET_INTERP_POINTS


! Locate maxima in the organic matter and black carbon fields, assumed
! to be forest fires erroneously in the climatology, and compute a
! scaling to remove them in the visibility calculation. If the field
! is 4D (i.e. has a long time variation) then we assume we have the
! newer climatology in which erroneous fires are not present.
SUBROUTINE LOCATE_FIRES(SELF)

  USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOMLUN,  ONLY : NULOUT

  ! What is the maximum ratio of a surface aerosol load and the
  ! maximum of its four immediate neighbours?
  REAL(KIND=JPRB), PARAMETER :: ZMAX_RATIO = 1.25_JPRB

  ! Maximum mass mixing ratio (kg/kg) of organic matter or black
  ! carbon to prevent visibility of less than around 10 km from an
  ! individual aerosol type (although note that only OM ever gets to
  ! this level.
  REAL(KIND=JPRB), PARAMETER :: ZMAX_OM_BC = 0.5E-7_JPRB
  
  ! Arguments
  CLASS(TEAERC), INTENT(INOUT) :: SELF

  INTEGER(KIND=JPIM) :: JAER       ! Loop index for aerosol type
  INTEGER(KIND=JPIM) :: JLON, JLAT ! Loop indices for location

  ! Max neighbour value
  REAL(KIND=JPRB) :: ZMAX_NEIGHBOUR
  ! Limiting value
  REAL(KIND=JPRB) :: ZLIMIT
  
  ! Number of candidate fires found per aerosol species
  INTEGER(KIND=JPIM) :: NFIRES(SELF%NAER)

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOEAERC:LOCATE_FIRES',0,ZHOOK_HANDLE)

  NFIRES = 0

  DO JAER = 1,SELF%NAER
    ! Only organic matter and black carbon that does not have a long
    ! time variation
    IF ((SELF%AER_TYPE(JAER) == 'OM' .AND. .NOT. SELF%IS_SPECIES_4D(JAER)) &
         .OR. (SELF%AER_TYPE(JAER) == 'BC' .AND. .NOT. SELF%IS_SPECIES_4D(JAER))) THEN
      SELF%SURFACE_SCALING(:,:,JAER) = 1.0_JPRB
      DO JLAT = 2,SELF%NLAT-1
        DO JLON = 2,SELF%NLON-1
          ZMAX_NEIGHBOUR = MAX(MAX(SELF%MMR(JLON-1,JLAT,  SELF%NLEV,JAER), &
               &                   SELF%MMR(JLON+1,JLAT,  SELF%NLEV,JAER)), &
               &               MAX(SELF%MMR(JLON,  JLAT-1,SELF%NLEV,JAER), &
               &                   SELF%MMR(JLON,  JLAT+1,SELF%NLEV,JAER)))
          ZLIMIT = MIN(ZMAX_NEIGHBOUR*ZMAX_RATIO, ZMAX_OM_BC)
          IF (SELF%MMR(JLON,JLAT,SELF%NLEV,JAER) > ZLIMIT) THEN
            SELF%SURFACE_SCALING(JLON,JLAT,JAER) = ZLIMIT/SELF%MMR(JLON,JLAT,SELF%NLEV,JAER)
            NFIRES(JAER) = NFIRES(JAER) + 1
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  IF (MAXVAL(NFIRES) > 0) THEN
    WRITE(NULOUT,'(A,I0,A,I0,A)') 'YOEAERC: Visibility correction for possible fires in aerosol climatology: ', &
         MAXVAL(NFIRES), ' OM/BC pixels out of ', SELF%NLON*SELF%NLAT, ' scaled down'
  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOEAERC:LOCATE_FIRES',1,ZHOOK_HANDLE)

END SUBROUTINE LOCATE_FIRES

END MODULE YOEAERC

