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

!****  YOMCOMPO
!
!     PURPOSE.
!     --------
!       MODULE CONTAINING VARIABLES FOR COMPOSITION VALID for GHG, AER and CHEMISTRY 
!
!     PARAMETER        DESCRIPTION                                  
!     ---------         -----------                                  !      
!
!     REFERENCE.
!     ----------
!                   
!     AUTHOR.
!     -------
!     2016-09-20   J. Flemming

MODULE YOMCOMPO

USE PARKIND1  ,ONLY : JPIM, JPRB

IMPLICIT NONE
SAVE

INTEGER(KIND=JPIM), PARAMETER :: NPEMIS2D = 500 ! sync with JPOSEMIS2D in parfpos.F90
INTEGER(KIND=JPIM), PARAMETER :: NPEMIS2DAUX = 10 ! sync with JPOSEMIS2DAUX in parfpos.F90
INTEGER(KIND=JPIM), PARAMETER :: NPEMIS3D = 10 ! sync with JPOSEMIS3D in parfpos.F90

! Type for a single emission specification
! Matches NAMCOMPO_EMIS namelist
TYPE :: TCOMPO_EMIS
  CHARACTER(LEN=16)  :: SPECIES = ''                       ! This should match the CNAME of the tracer's GFL descriptor
  CHARACTER(LEN=20)  :: SECTOR = ''                        ! Purely for information; the model code should NOT special-case this
  INTEGER(KIND=JPIM) :: PARAMID = -9999                    ! The ecCodes paramId of the input field
  INTEGER(KIND=JPIM) :: PARAM_INDEX = -9999                ! The index of the actual field in the YEMISxD array
  INTEGER(KIND=JPIM) :: NREQIN = 1                         ! Is emission required on input (for testing)
  REAL(KIND=JPRB)    :: SCALING = 1.0_JPRB                 ! Global scaling factor
  CHARACTER(LEN=16)  :: DIURNAL_CYCLE_TYPE = 'Uniform'     ! Which type of diurnal cycle to apply
  REAL(KIND=JPRB)    :: DIURNAL_AMPLITUDE = 0.0_JPRB       ! Amplitude of diurnal cycle (mean-to-peak, as fraction of mean)
  REAL(KIND=JPRB)    :: DIURNAL_PEAK_HOUR = 12.0_JPRB      ! Hour at which peak is reached (local solar time)
  REAL(KIND=JPRB)    :: DIURNAL_BASELINE = 0.0_JPRB        ! Baseline value for diurnal cycle
  ! ... optional parameters for other diurnal cycle types
  CHARACTER(LEN=16)  :: VERTICAL_PROFILE_TYPE = 'Surface'  ! Which type of vertical profile to apply
  REAL(KIND=JPRB)    :: VERTICAL_SURFACE_FRACTION = 0.0_JPRB ! Fraction to emit at surface ignoring profile
  REAL(KIND=JPRB)    :: VERTICAL_BASE_HEIGHT = 0.0_JPRB    ! Height (m) above surface or reference of plume base
  REAL(KIND=JPRB)    :: VERTICAL_TOP_HEIGHT = 0.0_JPRB     ! Height (m) above surface or reference of plume top
  REAL(KIND=JPRB)    :: VERTICAL_THRESHOLD = 0.0_JPRB      ! Threshold (m) below which to ignore map and treat as surface flux
  INTEGER(KIND=JPIM) :: VERTICAL_BASE_LEVEL = 0            ! Level number or offset of plume base
  INTEGER(KIND=JPIM) :: VERTICAL_TOP_LEVEL = 0             ! Level number or offset of plume top
  INTEGER(KIND=JPIM) :: VERTICAL_PARAMID = -9999           ! The ecCodes paramId of the input field of injection heights
  INTEGER(KIND=JPIM) :: VERTICAL_PARAM_INDEX = -9999       ! The index of the injection height field in the YEMIS2DAUX array
  ! ... optional parameters for other vertical profile types
  CHARACTER(LEN=8)   :: TEMPORALITY = 'Default'            ! Which type of online temporal processing to apply
  INTEGER(KIND=JPIM) :: LEGACY_CHEM_OVERRIDE = 0           ! Bitmask for legacy overrides:
                                                           !   1=DCBB; 2=INJ_CHEM; 4=VOLC_ALTI; 8=ANT_HIGH;
                                                           !   16=use-for-fire-inj; 32=use-for-volc-alti
  LOGICAL            :: NON_SIMPLE_TRACER = .FALSE.        ! Set if this is a non-simple-tracer species
                                                           !   (e.g. a GLOMAP component represented by multiple tracers)
END TYPE TCOMPO_EMIS

! Type for a single emission auxiliary fields specification
TYPE :: TCOMPO_EMIS_AUX
  INTEGER(KIND=JPIM) :: PARAMID = 0
  CHARACTER(LEN=8)   :: TEMPORALITY = 'Default'
  CHARACTER(LEN=16)  :: CNAME = ''
END TYPE TCOMPO_EMIS_AUX

TYPE :: TCOMPO
!  mass diagnostics switch for global budget   
     LOGICAL    :: LCHEM_DIA
!  period (in hours) on which to output budget
     REAL(KIND=JPRB) :: RCHEM_DIA_PERIOD
!  output of photolysis rates in extra ) 
     LOGICAL    :: LCHEM_DDFLX
!  if depostion as part of diffusion flux, do not add before , direct input to tracer diffusion    
     LOGICAL    :: LCOMPO_DDFLX_DIR
! variable tropopause (humidity)
    LOGICAL :: LCHEM_TROPO
! Use EQSAM4CLIM
    LOGICAL :: LAEREQSAM4CLIM
! aerosol surface fluxes  
    LOGICAL :: LAEROSFC
! Climatological diurnal Cycle (cosine) for climatological dry deposition velocities (AER + CHEM)    
    LOGICAL :: LCOMPO_DCDD
! Use nitrate aerosol scheme
    LOGICAL :: LAERNITRATE
! Use aer resuspension
    LOGICAL :: LAERRESUSPENSION
! Use distinct SOA species in aerosol scheme
    LOGICAL :: LAERSOA
! SOA coupled with chemistry
    LOGICAL :: LAERSOA_COUPLED
! Output CH4 loss tendency from chemistry, and optionally apply online
    INTEGER(KIND=JPIM) :: KGHG_CHEMTEND_CH4
! Use nucleation in IFS-GLOMAP
    LOGICAL :: LAERNUCL
! Aersols scheme (aer/glomap)
    CHARACTER(LEN=10)    :: AERO_SCHEME
! 2D emissions
    INTEGER(KIND=JPIM) :: NEMIS2D_DESC
    TYPE(TCOMPO_EMIS), DIMENSION(:), ALLOCATABLE :: YEMIS2D_DESC
! 3D emissions
    INTEGER(KIND=JPIM) :: NEMIS3D_DESC
    TYPE(TCOMPO_EMIS), DIMENSION(:), ALLOCATABLE :: YEMIS3D_DESC
! 2D emission auxiliary fields (injection heights etc)
    TYPE(TCOMPO_EMIS_AUX), DIMENSION(:), ALLOCATABLE :: YEMIS2DAUX_DESC
END TYPE TCOMPO 

!!TYPE(TCOMPO), POINTER :: YRCOMPO => NULL()
CHARACTER(LEN=10)    :: AERO_SCHEME_OP_OBS

CONTAINS

FUNCTION COMPO_EMIS_READ_NAMELIST(THIS_OUT) RESULT(LLOK)
  USE YOMLUN, ONLY : NULNAM
  TYPE(TCOMPO_EMIS), INTENT(OUT) :: THIS_OUT
  TYPE(TCOMPO_EMIS) :: THIS
  INTEGER :: IOSTAT
  CHARACTER(LEN=8) :: CLSTAT
  LOGICAL :: LLOK
  NAMELIST / NAMCOMPO_EMIS / THIS
  READ(NULNAM, NAMCOMPO_EMIS, IOSTAT=IOSTAT)
  SELECT CASE(IOSTAT)
    CASE(0)
      LLOK = (THIS%SPECIES /= '' .OR. THIS%SECTOR /= '' .OR. THIS%PARAMID>0) 
      IF (LLOK) THIS_OUT=THIS
    CASE DEFAULT
      WRITE(CLSTAT,'(I5)') IOSTAT
      CALL ABOR1 ('YOMCOMPO:COMPO_EMIS_READ_NAMELIST:READ ERROR ' // TRIM(CLSTAT) // ' ON NAMELIST FILE')
  END SELECT

  IF (THIS%VERTICAL_BASE_HEIGHT > THIS%VERTICAL_TOP_HEIGHT) THEN
    CALL ABOR1 ('YOMCOMPO:COMPO_EMIS_READ_NAMELIST:VERTICAL_TOP_HEIGHT SHOULD BE NO LESS THAN VERTICAL_BASE_HEIGHT')
  ENDIF

  IF (THIS%VERTICAL_BASE_LEVEL < THIS%VERTICAL_TOP_LEVEL) THEN
    CALL ABOR1 ('YOMCOMPO:COMPO_EMIS_READ_NAMELIST:VERTICAL_TOP_LEVEL SHOULD BE NO GREATER THAN VERTICAL_BASE_LEVEL')
  ENDIF

  IF (THIS%PARAM_INDEX /= -9999) THEN
    CALL ABOR1 ('YOMCOMPO:COMPO_EMIS_READ_NAMELIST:PARAM_INDEX SHOULD NOT BE EXPLICITLY SPECIFIED')
  ENDIF

  IF (THIS%VERTICAL_PARAM_INDEX /= -9999) THEN
    CALL ABOR1 ('YOMCOMPO:COMPO_EMIS_READ_NAMELIST:VERTICAL_PARAM_INDEX SHOULD NOT BE EXPLICITLY SPECIFIED')
  ENDIF

  IF (THIS%LEGACY_CHEM_OVERRIDE - IAND(THIS%LEGACY_CHEM_OVERRIDE, 63) /= 0) THEN
    CALL ABOR1 ('YOMCOMPO:COMPO_EMIS_READ_NAMELIST:LEGACY_CHEM_OVERRIDE SHOULD BE 0-63 (USING BITS 0-5 ONLY)')
  ENDIF
END FUNCTION COMPO_EMIS_READ_NAMELIST

FUNCTION COMPO_EMIS_AUX_READ_NAMELIST(THIS_OUT) RESULT(LLOK)
  USE YOMLUN, ONLY : NULNAM
  TYPE(TCOMPO_EMIS_AUX), INTENT(OUT) :: THIS_OUT
  TYPE(TCOMPO_EMIS_AUX) :: THIS
  INTEGER :: IOSTAT
  CHARACTER(LEN=8) :: CLSTAT
  LOGICAL :: LLOK
  NAMELIST / NAMCOMPO_EMIS_AUX / THIS
  READ(NULNAM, NAMCOMPO_EMIS_AUX, IOSTAT=IOSTAT)
  SELECT CASE(IOSTAT)
    CASE(0)
      LLOK = (THIS%PARAMID /= 0)
      IF (LLOK) THIS_OUT=THIS
    CASE DEFAULT
      WRITE(CLSTAT,'(I5)') IOSTAT
      CALL ABOR1 ('YOMCOMPO:COMPO_EMIS_AUX_READ_NAMELIST:READ ERROR ' // TRIM(CLSTAT) // ' ON NAMELIST FILE')
  END SELECT

  IF (THIS%CNAME /= '') THEN
    CALL ABOR1 ('YOMCOMPO:COMPO_EMIS_AUX_READ_NAMELIST:CNAME SHOULD NOT BE EXPLICITLY SPECIFIED')
  ENDIF
END FUNCTION COMPO_EMIS_AUX_READ_NAMELIST

END MODULE YOMCOMPO 
