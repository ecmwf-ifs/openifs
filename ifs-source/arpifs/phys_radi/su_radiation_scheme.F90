! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SU_RADIATION_SCHEME(YDMODEL,YDERDI,YDEAERATM,YDEPHY,YDERAD,PRADIATION,LDOUTPUT)
  
  ! SU_RADIATION_SCHEME - Setting up modular radiation scheme ecRad
  !
  ! PURPOSE
  ! -------
  !   The ecRad radiation scheme is contained in a separate
  !   library. SU_RADIATION_SCHEME sets up a small derived type that
  !   contains the configuration object for the radiation scheme, plus a
  !   small number of additional variables needed for its implemenation
  !   in the IFS.
  !
  ! INTERFACE
  ! ---------
  !   SU_RADIATION_SCHEME is called from SUECRAD.  The radiation scheme
  !   is actually run using the RADIATION_SCHEME routine.
  !
  ! AUTHOR
  ! ------
  !   Robin Hogan, ECMWF
  !   Original: 2015-09-16
  !
  ! MODIFICATIONS
  ! -------------
  !   2017-03-03  R. Hogan   Put global variables in TRADIATION derived type
  !   2017-11-17  S. Remy    Add Nitrates and SOA if NAERMACC=0
  !   2017-11-28  R. Hogan   Delta scaling applied to particles only
  !   2018-01-11  R. Hogan   Capability to scale solar spectrum in each band
  !   2018-04-20  A. Bozzo   Added capability to read in aerosol optical properties
  !                          at selected wavelengths
  !   2019-01-21  R. Hogan   Explicit albedo and emissivity spectral definitions
  !                          leading to smarter weighting in ecRad
  !   2021-01-24  R. Hogan   Added ecCKD via NSWGASOPTICS and NLWGASOPTICS
  !   2022-04-13  R. Hogan   Use YAER_RAD_DESC for aerosol optical properties
  !   2022-04-13  R. Hogan   Converted from module to free subroutine
  
  !-----------------------------------------------------------------------
  
  USE PARKIND1,         ONLY : JPRB,JPIM
  USE radiation_config, ONLY : ISolverMcICA, ISolverSpartacus, &
       &                       ISolverTripleclouds, &
       &                       ILiquidModelSlingo, ILiquidModelSOCRATES, &
       &                       IIceModelFu, IIceModelBaran, IIceModelYi, &
       &                       IOverlapExponential, IOverlapMaximumRandom, &
       &                       IOverlapExponentialRandom, IGasModelECCKD, IGasModelIFSRRTMG

  ! This routine copies information between the IFS radiation
  ! configuration (stored mostly in YDERAD) and the radiation
  ! configuration of the modular radiation scheme (stored in
  ! PRADIATION%rad_config).  The optional input logical LDOUTPUT
  ! controls whether to print lots of information during the setup
  ! stage (default is no).

  USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOMLUN,   ONLY : NULNAM, NULOUT, NULERR
  USE TYPE_MODEL,ONLY: MODEL
  USE YOERAD,   ONLY : TERAD
  USE YOEPHY,   ONLY : TEPHY
  USE YOEAERATM,ONLY : TEAERATM
  USE YOERDI  , ONLY : TERDI
  USE YOERADIATION,ONLY:TRADIATION

  USE RADIATION_INTERFACE,      ONLY : SETUP_RADIATION
  USE RADIATION_AEROSOL_OPTICS, ONLY : DRY_AEROSOL_MASS_EXTINCTION
  USE RADIATION_AEROSOL_OPTICS_DESCRIPTION, ONLY : AEROSOL_OPTICS_DESCRIPTION_TYPE

  IMPLICIT NONE


  TYPE(MODEL)       ,INTENT(IN)             :: YDMODEL
  ! Radiation configuration information
  TYPE(TERDI)       ,INTENT(INOUT)          :: YDERDI
  TYPE(TEAERATM)    ,INTENT(INOUT)          :: YDEAERATM
  TYPE(TEPHY)       ,INTENT(IN)             :: YDEPHY
  TYPE(TERAD)       ,INTENT(INOUT)          :: YDERAD
  TYPE(TRADIATION)  ,INTENT(INOUT), TARGET  :: PRADIATION

  ! Whether or not to print out information on the radiation scheme
  ! configuration
  LOGICAL, INTENT(IN), OPTIONAL :: LDOUTPUT

  ! Verbosity of configuration information 0=none, 1=warning,
  ! 2=info, 3=progress, 4=detailed, 5=debug
  INTEGER(KIND=JPIM) :: IVERBOSESETUP
  INTEGER(KIND=JPIM) :: ISTAT

  ! Data directory name
  CHARACTER(LEN=512) :: CL_DATA_DIR

  ! Type for accessing aerosol optics metadata
  TYPE(AEROSOL_OPTICS_DESCRIPTION_TYPE) :: AER_DESC

  ! Arrays to avoid temporaries
  REAL(KIND=JPRB)    :: ZWAVBOUND(15)
  INTEGER(KIND=JPIM) :: IBAND(16)

  ! Do we use the nearest ecRad band to the albedo/emissivity
  ! intervals, or a more intelligent weighting?
  LOGICAL :: LL_DO_NEAREST_SW_ALBEDO, LL_DO_NEAREST_LW_EMISS

  ! Old Tegen aerosol scheme diagnosed by other logicals being false:
  ! the logic setting this needs to be the same as in
  ! su_radiation_scheme.F90
  LOGICAL :: LL_USE_TEGEN_AEROSOLS
  
  ! Aerosol-type loop index
  INTEGER(KIND=JPIM) :: JA,JTAB

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "posname.intfb.h"
#include "abor1.intfb.h"

  IF (LHOOK) CALL DR_HOOK('SU_RADIATION_SCHEME',0,ZHOOK_HANDLE)

  ! *** GENERAL SETUP ***
  ASSOCIATE(RAD_CONFIG=>PRADIATION%RAD_CONFIG,&
       & LAERVOL=>YDEAERATM%LAERVOL, &
       & YSPECTPLANCK=>YDERAD%YSPECTPLANCK, &
       & YAER_RAD_DESC=>YDERAD%YAER_RAD_DESC, &
       & NACTAERO=>YDMODEL%YRML_GCONF%YGFL%NACTAERO, &
       & YAERO_DESC=>YDEAERATM%YAERO_DESC, &
       & ITYPE_TROP_BG_AER=>PRADIATION%ITYPE_TROP_BG_AER, &
       & ITYPE_STRAT_BG_AER=>PRADIATION%ITYPE_STRAT_BG_AER)

  ! Configure verbosity of setup of radiation scheme
  IVERBOSESETUP = 4 ! Provide plenty of information
  IF (PRESENT(LDOUTPUT)) THEN
    IF (.NOT. LDOUTPUT) THEN
      IVERBOSESETUP = 1 ! Warnings and errors only
    ENDIF
  ENDIF
  RAD_CONFIG%IVERBOSESETUP = IVERBOSESETUP

  IF (IVERBOSESETUP > 1) THEN
    WRITE(NULOUT,'(a)') '-------------------------------------------------------------------------------'
    WRITE(NULOUT,'(a)') 'RADIATION_SETUP: ecRad 1.6'
  ENDIF

  ! Normal operation of the radiation scheme displays only errors
  ! and warnings
  RAD_CONFIG%IVERBOSE = 1

  ! Read data directory name from the DATA environment variable
  CALL GETENV("DATA", CL_DATA_DIR)
  IF (CL_DATA_DIR /= " ") THEN
    RAD_CONFIG%DIRECTORY_NAME = TRIM(CL_DATA_DIR) // "/ifsdata"
  ELSE
    ! If DATA not present, use the current directory
    RAD_CONFIG%DIRECTORY_NAME = "."
  ENDIF

  ! Do we do Hogan and Bozzo (2015) approximate longwave updates?
  RAD_CONFIG%DO_LW_DERIVATIVES = YDERAD%LAPPROXLWUPDATE

  ! If we are to perform Hogan and Bozzo (2015) approximate
  ! shortwave updates then we need the downwelling direct and
  ! diffuse shortwave fluxes at the surface in each albedo spectral
  ! interval
  RAD_CONFIG%DO_CANOPY_FLUXES_SW = YDERAD%LAPPROXSWUPDATE

  ! If we are to perform approximate longwave updates and we are
  ! using the new 6-interval longwave emissivity scheme then we need
  ! ecRad to compute the downwelling surface longwave fluxes in each
  ! emissivity spectral interval
  IF (YDERAD%NLWOUT > 1) THEN
    RAD_CONFIG%DO_CANOPY_FLUXES_LW = .TRUE.
  ENDIF

  ! Surface spectral fluxes are needed for UV and PAR calculations
  RAD_CONFIG%DO_SURFACE_SW_SPECTRAL_FLUX = .TRUE.

  ! *** SETUP GAS OPTICS ***

  IF (YDERAD%NSWGASOPTICS == 0 .AND. YDERAD%NLWGASOPTICS == 0) THEN

    ! Classic RRTMG gas optics scheme.  We assume the IFS has
    ! already set-up RRTMG, so the setup_radiation routine below
    ! does not have to.
    RAD_CONFIG%I_GAS_MODEL                = IGasModelIFSRRTMG
    RAD_CONFIG%USE_GENERAL_CLOUD_OPTICS   = .FALSE.
    RAD_CONFIG%USE_GENERAL_AEROSOL_OPTICS = .TRUE.
    RAD_CONFIG%DO_SETUP_IFSRRTM           = .FALSE.
    ! Don't weight by Planck function (of sun or representative
    ! terrestrial temperaure) when mapping from surface
    ! albedo/emissivity bands to gas-optics bands
    RAD_CONFIG%DO_WEIGHTED_SURFACE_MAPPING= .FALSE.

  ELSEIF (YDERAD%NSWGASOPTICS == 0 .OR. YDERAD%NLWGASOPTICS == 0) THEN

    WRITE(NULERR,'(a)') '*** Error: NSWGASOPTICS and NLWGASOPTICS must either both be zero or both nonzero'
    CALL ABOR1('RADIATION_SETUP: error interpreting NSWGASOPTICS and NLWGASOPTICS')

  ELSE

    ! ecCKD gas optics scheme
    RAD_CONFIG%I_GAS_MODEL                = IGasModelECCKD
    RAD_CONFIG%USE_GENERAL_CLOUD_OPTICS   = .TRUE.
    RAD_CONFIG%USE_GENERAL_AEROSOL_OPTICS = .TRUE.

    ! If NSWGASOPTICS or NLWGASOPTICS are negative then use the
    ! default models; if positive then interpret as the number of g
    ! points
    IF (YDERAD%NSWGASOPTICS == 64) THEN
      ! 19-band model with near-IR bands for each window to
      ! accurately capture cloud absorption, and fine bands in the
      ! UV for UV index
      RAD_CONFIG%GAS_OPTICS_SW_OVERRIDE_FILE_NAME &
           & = 'ecckd-1.2_sw_climate_window-64b_ckd-definition.nc'
    ELSEIF (YDERAD%NSWGASOPTICS > 0) THEN
      ! Usually the 16- and 32-point models are available
      WRITE(RAD_CONFIG%GAS_OPTICS_SW_OVERRIDE_FILE_NAME, '(a,i0,a)') &
           &  'ecckd-1.0_sw_climate_rgb-', YDERAD%NSWGASOPTICS, &
           &  'b_ckd-definition.nc'
    ELSE
      RAD_CONFIG%GAS_OPTICS_SW_OVERRIDE_FILE_NAME &
           & = 'ecckd-1.0_sw_climate_rgb-32b_ckd-definition.nc'
    ENDIF

    IF (YDERAD%NLWGASOPTICS == 64) THEN
      ! 13-band model
      RAD_CONFIG%GAS_OPTICS_LW_OVERRIDE_FILE_NAME &
           & = 'ecckd-1.2_lw_climate_narrow-64b_ckd-definition.nc'
    ELSEIF (YDERAD%NLWGASOPTICS > 0) THEN
      ! Usually the 16- and 32-point models are available
      WRITE(RAD_CONFIG%GAS_OPTICS_LW_OVERRIDE_FILE_NAME, '(a,i0,a)') &
           &  'ecckd-1.0_lw_climate_fsck-', YDERAD%NLWGASOPTICS, &
           &  'b_ckd-definition.nc'
    ELSE
      RAD_CONFIG%GAS_OPTICS_LW_OVERRIDE_FILE_NAME &
           &  = 'ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc'
    ENDIF

  ENDIF


  ! *** SETUP CLOUD OPTICS ***

  ! Setup liquid optics for RRTMG configuration
  IF (YDERAD%NLIQOPT == 2) THEN
    RAD_CONFIG%I_LIQ_MODEL = ILIQUIDMODELSLINGO
  ELSEIF (YDERAD%NLIQOPT == 4) THEN
    RAD_CONFIG%I_LIQ_MODEL = ILIQUIDMODELSOCRATES
  ELSE
    WRITE(NULERR,'(a,i0)') '*** Error: Unavailable liquid optics model in modular radiation scheme: NLIQOPT=', &
         &  YDERAD%NLIQOPT
    CALL ABOR1('RADIATION_SETUP: error interpreting NLIQOPT')   
  ENDIF
  ! Setup liquid optics for generalized cloud configuration
  RAD_CONFIG%CLOUD_TYPE_NAME(1) = "mie_droplet"

  ! Setup ice optics
  IF (YDERAD%NICEOPT == 3) THEN
    RAD_CONFIG%I_ICE_MODEL = IICEMODELFU ! ecRad-RRTMG configuration
    RAD_CONFIG%CLOUD_TYPE_NAME(2) = "fu-muskatel_ice" ! ecRad generalized configuration
    IF (YDERAD%LFU_LW_ICE_OPTICS_BUG) THEN
      RAD_CONFIG%DO_FU_LW_ICE_OPTICS_BUG = .TRUE.
    ENDIF
  ELSEIF (YDERAD%NICEOPT == 4) THEN
    RAD_CONFIG%I_ICE_MODEL = IICEMODELBARAN
    IF (RAD_CONFIG%I_GAS_MODEL == IGasModelECCKD) THEN
      WRITE(NULERR,'(a,i0)') '*** Error: Baran ice optics unavailable with generalized cloud optics'
      CALL ABOR1('RADIATION_SETUP: error interpreting NICEOPT')
    ENDIF
  ELSEIF (YDERAD%NICEOPT == 5) THEN
    RAD_CONFIG%I_ICE_MODEL = IICEMODELYI
    RAD_CONFIG%CLOUD_TYPE_NAME(2) = "baum-general-habit-mixture_ice"
  ELSEIF (YDERAD%NICEOPT == 6) THEN
    RAD_CONFIG%I_ICE_MODEL = IICEMODELFU ! ecRad-RRTMG configuration
    RAD_CONFIG%CLOUD_TYPE_NAME(2) = "fu-muskatel-rough_ice" ! ecRad generalized configuration
  ELSE
    WRITE(NULERR,'(a,i0)') '*** Error: Unavailable ice optics model in modular radiation scheme: NICEOPT=', &
         &  YDERAD%NICEOPT
    CALL ABOR1('RADIATION_SETUP: error interpreting NICEOPT')
  ENDIF

  ! For consistency with earlier versions of the IFS radiation
  ! scheme, until 45R1 we performed shortwave delta-Eddington
  ! scaling after the merge of the cloud, aerosol and gas optical
  ! properties.  Setting this to "false" does the scaling on the
  ! cloud and aerosol properties separately before merging with
  ! gases, which is more physically appropriate. The impact is very
  ! small (see item 6 of table 2 of Technical Memo 787).
  RAD_CONFIG%DO_SW_DELTA_SCALING_WITH_GASES = .FALSE.


  ! *** SETUP AEROSOLS ***

  ! Configure aerosol properties in the RAD_CONFIG structure
  RAD_CONFIG%USE_AEROSOLS = .TRUE.

  IF (YDEAERATM%LAERCCN .OR. YDEAERATM%LAERRRTM .OR. YDERAD%NAERMACC == 1 &
       &  .OR. NACTAERO > 0) THEN
    ! Using MACC climatology or prognostic aerosol variables
    LL_USE_TEGEN_AEROSOLS = .FALSE.
  ELSE
    ! Using Tegen climatology
    LL_USE_TEGEN_AEROSOLS = .TRUE.
    IF (RAD_CONFIG%I_GAS_MODEL == IGasModelECCKD) THEN
      WRITE(NULERR,'(a,i0)') '*** Error: Tegen aerosol climatology incompatible with ecCKD gas optics scheme'
      CALL ABOR1('RADIATION_SETUP: Tegen/ecCKD incompatibility')
    ELSE
      RAD_CONFIG%USE_GENERAL_AEROSOL_OPTICS = .FALSE.
    ENDIF
  ENDIF

  ! Flexible configuration of aerosol species seen by radiation: the
  ! YAER_RAD_DESC describes the aerosols seen by radiation; this can
  ! be configured on the namelist, otherwise it will either be
  ! entirely diagnostic or entirely prognostic
  CALL YAER_RAD_DESC%READ_NAMELIST(NULNAM)
  IF (YAER_RAD_DESC%NAEROSOL == 0) THEN
    ! No aerosols have been described on the namelist
    IF (NACTAERO > 0 .AND. YDEAERATM%LAERRRTM) THEN
      ! We have prognostic aerosols (NACTAERO>0) and they are
      ! interactive with radiation (LAERRRTM=true): copy over
      ! information from the aerosol GFL fields described in
      ! YAERO_DESC
      CALL YAER_RAD_DESC%RESERVE(32)
      CALL YAER_RAD_DESC%ADD_ALL_PROGNOSTICS(YAERO_DESC)
      IF (YDEAERATM%LAERFILLCLIM) THEN
        ! Add climatological only if the prognostic for that species
        ! is not already in YAER_RAD_DESC
        CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Sea_Salt_bin1','SS',.TRUE.,BIN=1)
        CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Sea_Salt_bin2','SS',.TRUE.,BIN=2)
        CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Sea_Salt_bin3','SS',.TRUE.,BIN=3)
        ! It is possible that the prognostics contain hydrophilic
        ! dust, in which case we don't want to add hydrophobic
        ! climatological dust, therefore we check there is no
        ! hydrophilic dust for each bin.  ADD_CLIMATOLOGICAL already
        ! checks there is no hydrophobic dust present.
        IF (YAER_RAD_DESC%GET_INDEX('DD',.TRUE.,1) == 0) THEN
          CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Mineral_Dust_bin1','DD',.FALSE.,BIN=1)
        ENDIF
        IF (YAER_RAD_DESC%GET_INDEX('DD',.TRUE.,2) == 0) THEN
          CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Mineral_Dust_bin2','DD',.FALSE.,BIN=2)
        ENDIF
        IF (YAER_RAD_DESC%GET_INDEX('DD',.TRUE.,3) == 0) THEN
          CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Mineral_Dust_bin3','DD',.FALSE.,BIN=3)
        ENDIF
        CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Organic_Matter_hydrophilic','OM',.TRUE.)
        CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Organic_Matter_hydrophobic','OM',.FALSE.)
        ! Note that "hydrophilic" black carbon is treated as hydrophobic
        ! in radiation, but we need to keep the types distinguishable
        ! here since adding is conditional on the type not being present. 
        CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Black_Carbon_hydrophilic','BC',.TRUE.)
        CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Black_Carbon_hydrophobic','BC',.FALSE.)
        CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Sulfates','SU',.TRUE.)
      ENDIF
    ELSEIF (.NOT. LL_USE_TEGEN_AEROSOLS) THEN
      ! Classic CAMS climatology, storing two of the indices to use
      ! for tropospheric and stratospheric "background" aerosol. Note
      ! that the preferred optics models are by default different for
      ! prognostic and climatological aerosol, but can be overridden
      ! with the NAERAD namelist; therefore they set by the
      ! CAEROPTICSMODEL_?? strings in suecrad.F90
      CALL YAER_RAD_DESC%RESERVE(13)
      CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Sea_Salt_bin1','SS',.TRUE.,BIN=1)
      CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Sea_Salt_bin2','SS',.TRUE.,BIN=2)
      CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Sea_Salt_bin3','SS',.TRUE.,BIN=3)
      CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Mineral_Dust_bin1','DD',.FALSE.,BIN=1)
      CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Mineral_Dust_bin2','DD',.FALSE.,BIN=2)
      CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Mineral_Dust_bin3','DD',.FALSE.,BIN=3)
      CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Organic_Matter_hydrophilic','OM',.TRUE.)
      CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Organic_Matter_hydrophobic','OM',.FALSE.)
      ! Note that "hydrophilic" black carbon is treated as hydrophobic
      ! in radiation, but we need to keep the types distinguishable
      ! here since adding is conditional on the type not being present. 
      CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Black_Carbon_hydrophilic','BC',.TRUE.)
      CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Black_Carbon_hydrophobic','BC',.FALSE.)
      CALL YAER_RAD_DESC%ADD_CLIMATOLOGICAL('Sulfates','SU',.TRUE.)
      CALL YAER_RAD_DESC%ADD_PARAMETRIC('Tropospheric_background_organic',0,'OM',.FALSE., &
           &  KINDEX=ITYPE_TROP_BG_AER)
      CALL YAER_RAD_DESC%ADD_PARAMETRIC('Stratospheric_background_sulfates',1,'SU',.FALSE., &
           &  KINDEX=ITYPE_STRAT_BG_AER)
    ELSE
      WRITE(NULOUT,'(A)') 'Warning: using older Tegen aerosol climatology'
    ENDIF
    IF (YAER_RAD_DESC%NAEROSOL > 0) THEN
      CALL YAER_RAD_DESC%SET_HYDROPHOBIC('BC')
    ENDIF
  ELSE
    ! Aerosols have been described manually via the namelist
    ! NAM_AEROSOL_RADIATION: this is incompatible with the Tegen
    ! climatology
    IF (LL_USE_TEGEN_AEROSOLS) THEN
      !
      WRITE(NULOUT,'(A)') 'SU_RADIATION_SCHEME: Warning: aerosol specified manually, so setting NAERMACC=1'
      YDERAD%NAERMACC = 1
      LL_USE_TEGEN_AEROSOLS = .FALSE.
    ENDIF
    ! Locate any prognostic aerosols amongst the GFL fields
    CALL YAER_RAD_DESC%FIND_PROGNOSTICS(YAERO_DESC)
    ! Also locate any tropospheric or stratospheric "background"
    ! aerosols, which have PARAMETRIC_SCHEME parameters set to 0 and
    ! 1, respectively
    ITYPE_TROP_BG_AER  = YAER_RAD_DESC%PARAMETRIC_INDEX(0)
    ITYPE_STRAT_BG_AER = YAER_RAD_DESC%PARAMETRIC_INDEX(1)
  ENDIF

  ! Print aerosol-radiation description, either in the form of a
  ! namelist or a more concise format
  !CALL YAER_RAD_DESC%PRINT(NULOUT)
  CALL YAER_RAD_DESC%PRINT_NAMELIST(NULOUT)

  ! If monochromatic aerosol properties are available they will be
  ! read in automatically so the following is not needed
  !IF (YDEAERATM%LAERRAD) RAD_CONFIG%AEROSOL_OPTICS%READ_MONOCHROMATIC_OPTICS=.TRUE.

  ! But we need to specify the wavelengths for the general aerosol
  ! optics
  CALL RAD_CONFIG%SET_AEROSOL_WAVELENGTH_MONO( &
       &  [3.4e-07_JPRB, 3.55e-07_JPRB, 3.8e-07_JPRB, 4.0e-07_JPRB, 4.4e-07_JPRB, &
       &   4.69e-07_JPRB, 5.0e-07_JPRB, 5.32e-07_JPRB, 5.5e-07_JPRB, 6.45e-07_JPRB, &
       &   6.7e-07_JPRB, 8.0e-07_JPRB, 8.58e-07_JPRB, 8.65e-07_JPRB, 1.02e-06_JPRB, &
       &   1.064e-06_JPRB, 1.24e-06_JPRB, 1.64e-06_JPRB, 2.13e-06_JPRB, 1.0e-05_JPRB])

  IF (.NOT. LL_USE_TEGEN_AEROSOLS) THEN
    ! Using MACC climatology or prognostic aerosol variables
    
    ! Perhaps the following is in the wrong place because it is before
    ! the RADIATION namelist group has been read so the user cannot
    ! override the aerosol optics file
    IF (.NOT. RAD_CONFIG%USE_GENERAL_AEROSOL_OPTICS) THEN
      ! If not using "general" aerosol optics (i.e. compatible with
      ! all gas optics models, then the default aerosol optics file is
      ! the following - please update here, not in
      ! radiation/module/radiation_config.F90.  This file is only
      ! compatible with the RRTM bands.
      RAD_CONFIG%AEROSOL_OPTICS_OVERRIDE_FILE_NAME = 'aerosol_ifs_rrtm_49R1_20230119.nc'
    ELSE
      ! This file is compatible with all spectral discretizations
      RAD_CONFIG%AEROSOL_OPTICS_OVERRIDE_FILE_NAME = 'aerosol_ifs_49R1_20230725.nc'
    ENDIF

    ! Load aerosol optics metadata; ecRad will later read the actual
    ! optical properties
    IF (RAD_CONFIG%AEROSOL_OPTICS_OVERRIDE_FILE_NAME(1:1) == "/") THEN
      CALL AER_DESC%READ(TRIM(RAD_CONFIG%AEROSOL_OPTICS_OVERRIDE_FILE_NAME), IVERBOSE=IVERBOSESETUP)
    ELSE
      CALL AER_DESC%READ(TRIM(RAD_CONFIG%DIRECTORY_NAME) // '/' &
           &  // TRIM(RAD_CONFIG%AEROSOL_OPTICS_OVERRIDE_FILE_NAME), IVERBOSE=IVERBOSESETUP)
    ENDIF

    ! Set preferred optics model for desert dust, black carbon,
    ! organic matter and sulfate; these are set in suecrad.F90 and may
    ! be difference between prognostic and climatological aerosols
    CALL AER_DESC%PREFERRED_OPTICAL_MODEL('DD',TRIM(YDERAD%CAEROPTICSMODEL_DD))
    CALL AER_DESC%PREFERRED_OPTICAL_MODEL('BC',TRIM(YDERAD%CAEROPTICSMODEL_BC))
    CALL AER_DESC%PREFERRED_OPTICAL_MODEL('OM',TRIM(YDERAD%CAEROPTICSMODEL_OM))
    CALL AER_DESC%PREFERRED_OPTICAL_MODEL('SU',TRIM(YDERAD%CAEROPTICSMODEL_SU))

    ! Read aersol-radiation settings from YAER_RAD_DESC structure
    RAD_CONFIG%N_AEROSOL_TYPES = YAER_RAD_DESC%NAEROSOL

    ! Indices to the aerosol optical properties in
    ! aerosol_ifs_rrtm_*.nc, for each class, where negative numbers
    ! index hydrophilic aerosol types and positive numbers index
    ! hydrophobic aerosol types
    RAD_CONFIG%I_AEROSOL_TYPE_MAP = 0 ! There can be up to 256 types

    ! Use the AER_DESC%GET_INDEX function to get indices to the
    ! optical properties of each aerosol species based on properties
    ! in YAER_RAD_DESC, and put them in I_AEROSOL_TYPE_MAP. If a zero
    ! is returned then the aerosol was not found, which is an error.
    DO JA = 1,RAD_CONFIG%N_AEROSOL_TYPES
      ASSOCIATE(DATA=>YAER_RAD_DESC%DATA(JA))
      IF (TRIM(DATA%OPTICS_MODEL) == ' ') THEN
        ! Optics model not specified for this type
        RAD_CONFIG%I_AEROSOL_TYPE_MAP(JA) &
             &  = AER_DESC%GET_INDEX(DATA%AER_TYPE, DATA%LHYDROPHILIC, &
             &      IBIN=DATA%SIZE_BIN)
        IF (RAD_CONFIG%I_AEROSOL_TYPE_MAP(JA) == 0) THEN
          WRITE(NULERR,'(A,A2,A,L,A,I0)') '*** Error finding optical properties of aerosol ', &
               &  DATA%AER_TYPE, ', hydrophilic=', DATA%LHYDROPHILIC, ', bin=', DATA%SIZE_BIN
          CALL ABOR1('RADIATION_SETUP: Error finding aerosol optical properties')
        ENDIF
      ELSE
        ! Optics model is provided
        RAD_CONFIG%I_AEROSOL_TYPE_MAP(JA) &
             &  = AER_DESC%GET_INDEX(DATA%AER_TYPE, DATA%LHYDROPHILIC, &
             &      IBIN=DATA%SIZE_BIN, OPTICAL_MODEL_STR=TRIM(DATA%OPTICS_MODEL))
        IF (RAD_CONFIG%I_AEROSOL_TYPE_MAP(JA) == 0) THEN
          WRITE(NULERR,'(A,A2,A,L,A,I0)') '*** Error finding optical properties of aerosol ', &
               &  DATA%AER_TYPE, ', hydrophilic=', DATA%LHYDROPHILIC, ', bin=', DATA%SIZE_BIN, &
               &  ', optical_model=', TRIM(DATA%OPTICS_MODEL)
          CALL ABOR1('RADIATION_SETUP: Error finding aerosol optical properties')
        ENDIF
      ENDIF
      END ASSOCIATE
    ENDDO

    ! Background aerosol mass-extinction coefficients are obtained
    ! after the configuration files have been read - see later in this
    ! routine.
    
  ELSE
    ! Using Tegen climatology
    IF (RAD_CONFIG%USE_GENERAL_AEROSOL_OPTICS) THEN
      WRITE(NULERR,'(a,i0)') '*** Error: Tegen aerosol climatology incompatible with general aerosol optics'
      CALL ABOR1('RADIATION_SETUP: Tegen incompatibility')
    ENDIF

    RAD_CONFIG%N_AEROSOL_TYPES = 6
    RAD_CONFIG%I_AEROSOL_TYPE_MAP = 0 ! There can be up to 256 types
    RAD_CONFIG%I_AEROSOL_TYPE_MAP(1:6) = (/&
         &  1,&! Continental background
         &  2,&! Maritime
         &  3,&! Desert
         &  4,&! Urban
         &  5,&! Volcanic active
         &  6 /)  ! Stratospheric background

    ! Manually set the aerosol optics file name (the directory will
    ! be added automatically)
    RAD_CONFIG%AEROSOL_OPTICS_OVERRIDE_FILE_NAME = 'aerosol_ifs_rrtm_tegen.nc'
  ENDIF
  
  ! *** SETUP SOLVER ***
  
  ! 3D effects are off by default
  RAD_CONFIG%DO_3D_EFFECTS = .FALSE.

  ! Select longwave solver
  SELECT CASE (YDERAD%NLWSOLVER)
  CASE(0)
    RAD_CONFIG%I_SOLVER_LW = ISOLVERMCICA
  CASE(1)
    RAD_CONFIG%I_SOLVER_LW = ISOLVERSPARTACUS
  CASE(2)
    RAD_CONFIG%I_SOLVER_LW = ISOLVERSPARTACUS
    RAD_CONFIG%DO_3D_EFFECTS = .TRUE.
  CASE(3)
    RAD_CONFIG%I_SOLVER_LW = ISOLVERTRIPLECLOUDS
  CASE DEFAULT
    WRITE(NULERR,'(a,i0)') '*** Error: Unknown value for NLWSOLVER: ', YDERAD%NLWSOLVER
    CALL ABOR1('RADIATION_SETUP: error interpreting NLWSOLVER')
  END SELECT

  ! Select shortwave solver
  SELECT CASE (YDERAD%NSWSOLVER)
  CASE(0)
    RAD_CONFIG%I_SOLVER_SW = ISOLVERMCICA
  CASE(1)
    RAD_CONFIG%I_SOLVER_SW = ISOLVERSPARTACUS
    RAD_CONFIG%DO_3D_EFFECTS = .FALSE.
    IF (YDERAD%NLWSOLVER == 2) THEN
      CALL ABOR1('RADIATION_SETUP: cannot represent 3D effects in LW but not SW')
    ENDIF
  CASE(2)
    RAD_CONFIG%I_SOLVER_SW = ISOLVERSPARTACUS
    RAD_CONFIG%DO_3D_EFFECTS = .TRUE.
    IF (YDERAD%NLWSOLVER == 1) THEN
      CALL ABOR1('RADIATION_SETUP: cannot represent 3D effects in SW but not LW')
    ENDIF
  CASE(3)
    RAD_CONFIG%I_SOLVER_SW = ISOLVERTRIPLECLOUDS
  CASE DEFAULT
    WRITE(NULERR,'(a,i0)') '*** Error: Unknown value for NSWSOLVER: ', YDERAD%NSWSOLVER
    CALL ABOR1('RADIATION_SETUP: error interpreting NSWSOLVER')
  END SELECT

  ! For stability the cloud effective size can't be too small in
  ! SPARTACUS
  RAD_CONFIG%MIN_CLOUD_EFFECTIVE_SIZE = 500.0_JPRB

  ! SPARTACUS solver requires delta scaling to be done separately
  ! for clouds & aerosols
  IF (RAD_CONFIG%I_SOLVER_SW == ISOLVERSPARTACUS) THEN
    RAD_CONFIG%DO_SW_DELTA_SCALING_WITH_GASES = .FALSE.
  ENDIF

  ! Do we represent longwave scattering?
  RAD_CONFIG%DO_LW_CLOUD_SCATTERING = .FALSE.
  RAD_CONFIG%DO_LW_AEROSOL_SCATTERING = .FALSE.
  SELECT CASE (YDERAD%NLWSCATTERING)
  CASE(1)
    RAD_CONFIG%DO_LW_CLOUD_SCATTERING = .TRUE.
  CASE(2)
    RAD_CONFIG%DO_LW_CLOUD_SCATTERING = .TRUE.
    IF (.NOT. LL_USE_TEGEN_AEROSOLS) THEN
      ! Tegen climatology omits data required to do longwave
      ! scattering by aerosols, so only turn this on with a more
      ! recent scattering database
      RAD_CONFIG%DO_LW_AEROSOL_SCATTERING = .TRUE.
    ENDIF
  END SELECT

  SELECT CASE (YDERAD%NCLOUDOVERLAP)
  CASE (1)
    RAD_CONFIG%I_OVERLAP_SCHEME = IOVERLAPMAXIMUMRANDOM
  CASE (2)
    ! Use Exponential-Exponential cloud overlap to match original IFS
    ! implementation of Raisanen cloud generator
    RAD_CONFIG%I_OVERLAP_SCHEME = IOVERLAPEXPONENTIAL
  CASE (3)
    RAD_CONFIG%I_OVERLAP_SCHEME = IOVERLAPEXPONENTIALRANDOM
  CASE DEFAULT
    WRITE(NULERR,'(a,i0)') '*** Error: Unknown value for NCLOUDOVERLAP: ', YDERAD%NCLOUDOVERLAP
    CALL ABOR1('RADIATION_SETUP: error interpreting NCLOUDOVERLAP')
  END SELECT

  ! Change cloud overlap to exponential-random if Tripleclouds or
  ! SPARTACUS selected as both the shortwave and longwave solvers
  IF (RAD_CONFIG%I_OVERLAP_SCHEME /= IOVERLAPEXPONENTIALRANDOM &
       & .AND. (     RAD_CONFIG%I_SOLVER_SW == ISOLVERTRIPLECLOUDS &
       &        .OR. RAD_CONFIG%I_SOLVER_LW == ISOLVERTRIPLECLOUDS &
       &        .OR. RAD_CONFIG%I_SOLVER_SW == ISOLVERSPARTACUS &
       &        .OR. RAD_CONFIG%I_SOLVER_LW == ISOLVERSPARTACUS)) THEN
    IF (RAD_CONFIG%I_SOLVER_SW == RAD_CONFIG%I_SOLVER_LW) THEN
      WRITE(NULOUT,'(a)') 'Warning: Tripleclouds/SPARTACUS solver selected so changing cloud overlap to Exp-Ran'
      RAD_CONFIG%I_OVERLAP_SCHEME = IOVERLAPEXPONENTIALRANDOM
    ELSE
      ! If the solvers are not the same and exponential-random has
      ! not been selected then abort
      WRITE(NULERR,'(a)') '*** Error: Tripleclouds and SPARTACUS solvers can only simulate exponential-random overlap'
      CALL ABOR1('RADIATION_SETUP: Cloud overlap incompatible with solver')
    ENDIF
    
    ! For additional stability in SPARTACUS solver it helps if the
    ! cloud fraction threshold is higher than the default of 1.0e-6
    ! used for McICA; this is done for Tripleclouds too so that it
    ! is a good control for SPARTACUS.
    RAD_CONFIG%CLOUD_FRACTION_THRESHOLD = 2.5E-5_JPRB
  ENDIF
  
  ! Populate the mapping between the 14 RRTM shortwave bands and the
  ! 6 albedo inputs. 
  ZWAVBOUND(1:5) = [ 0.25e-6_jprb, 0.44e-6_jprb, 0.69e-6_jprb, &
       &             1.19e-6_jprb, 2.38e-6_jprb ]
  IBAND(1:6)  = [ 1,2,3,4,5,6 ]
  ! If NALBEDOSCHEME==2 then we are using the 6-component MODIS
  ! albedo climatology, and a weighted average is used to compute
  ! the albedos in each ecRad spectral band. If NALBEDOSCHEME==3
  ! then we use the diffuse part of the 4 components but still with
  ! a weighted average. Otherwise the older behaviour is followed:
  ! the nearest albedo interval to each band is selected, resulting
  ! in a discrete mapping that matches the one in YOESRTWN:NMPSRTM.
  ! Note that this tends to bias albedo high because there is a lot
  ! of energy around the interface between the UV-Vis and Near-IR
  ! channels, so this should be close to the 0.7 microns intended by
  ! the MODIS dataset, not shifted to the nearest RRTM band boundary
  ! at 0.625 microns. Note that if we use the ecCKD gas optics
  ! scheme then we use the weighted average.
  LL_DO_NEAREST_SW_ALBEDO = (YDEPHY%NALBEDOSCHEME<2 .AND. YDERAD%NSWGASOPTICS==0)
  CALL RAD_CONFIG%DEFINE_SW_ALBEDO_INTERVALS(6, ZWAVBOUND, IBAND, &
       &  DO_NEAREST=LL_DO_NEAREST_SW_ALBEDO)

  ! Likewise between the 16 RRTM longwave bands and the NLWEMISS
  ! emissivity inputs - these are defined in suecrad.F90. A weighted
  ! average is used for the ecCKD gas optics scheme due to the
  ! likelihood of the full-spectrum correlated-k method being used,
  ! for which nearest neighbour mapping is very poor.
  LL_DO_NEAREST_LW_EMISS = (YDERAD%NLWGASOPTICS==0)
  CALL RAD_CONFIG%DEFINE_LW_EMISS_INTERVALS(UBOUND(YSPECTPLANCK%INTERVAL_MAP,1), &
       &  YSPECTPLANCK%WAVLEN_BOUND, YSPECTPLANCK%INTERVAL_MAP, &
       &  DO_NEAREST=LL_DO_NEAREST_LW_EMISS)

  ! *** IMPLEMENT SETTINGS ***

  ! For advanced configuration, the configuration data for the
  ! "radiation" project can specified directly in the namelist.
  ! However, the variable naming convention is not consistent with
  ! the rest of the IFS.  For basic configuration there are specific
  ! variables in the NAERAD namelist available in the YDERAD
  ! structure.
  CALL POSNAME(NULNAM, 'RADIATION', ISTAT)
  SELECT CASE (ISTAT)
  CASE(0)
    CALL RAD_CONFIG%READ(UNIT=NULNAM)
  CASE(1)
    WRITE(NULOUT,'(a)') 'Namelist RADIATION not found, using settings from NAERAD only'
  CASE DEFAULT
    CALL ABOR1('RADIATION_SETUP: error reading RADIATION section of namelist file')
  END SELECT

  ! Print configuration
  IF (IVERBOSESETUP > 1) THEN
    WRITE(NULOUT,'(a)') 'Radiation scheme settings:'
    CALL RAD_CONFIG%PRINT(IVERBOSE=IVERBOSESETUP)
  ENDIF

  ! Do we use the solar cycle number to modify the solar spectrum,
  ! only available with shortwave ecCKD gas optics?
  IF (YDERAD%LSPECTRALSOLARCYCLE .AND. YDERAD%NSWGASOPTICS /= 0) THEN
    !RAD_CONFIG%SSI_OVERRIDE_FILE_NAME = '/home/parr/radiation_data/ssi_nrl2.nc'
    RAD_CONFIG%USE_SPECTRAL_SOLAR_CYCLE = .TRUE.
  ELSE
    RAD_CONFIG%USE_SPECTRAL_SOLAR_CYCLE = .FALSE.
  ENDIF
  
  ! Use configuration data to set-up radiation scheme, including
  ! reading scattering datafiles
  CALL SETUP_RADIATION(RAD_CONFIG)

  ! Do we scale the incoming solar radiation in each band?
  IF (YDERAD%NSOLARSPECTRUM > 0 &
       &  .AND. RAD_CONFIG%I_GAS_MODEL == IGasModelIFSRRTMG) THEN
    IF (RAD_CONFIG%N_BANDS_SW /= 14) THEN
      WRITE(NULERR,'(a,i0,a)') '*** Error: ', RAD_CONFIG%N_BANDS_SW, &
           &  ' shortwave bands but need 14 to apply spectral scaling'
      CALL ABOR1('RADIATION_SETUP: Shortwave must have 14 bands to apply spectral scaling')
    ELSE
      RAD_CONFIG%USE_SPECTRAL_SOLAR_SCALING = .TRUE.
    ENDIF
  ENDIF

  ! Get spectral weightings for UV and PAR
  CALL RAD_CONFIG%GET_SW_WEIGHTS(0.2E-6_JPRB, 0.4415E-6_JPRB,&
       &  PRADIATION%NWEIGHT_UV, PRADIATION%IBAND_UV, PRADIATION%WEIGHT_UV,&
       &  'ultraviolet')
  CALL RAD_CONFIG%GET_SW_WEIGHTS(0.4E-6_JPRB, 0.7E-6_JPRB,&
       &  PRADIATION%NWEIGHT_PAR, PRADIATION%IBAND_PAR, PRADIATION%WEIGHT_PAR,&
       &  'photosynthetically active radiation, PAR')

  PRADIATION%TROP_BG_AER_MASS_EXT  = 0.0_JPRB
  PRADIATION%STRAT_BG_AER_MASS_EXT = 0.0_JPRB
  IF (.NOT. LL_USE_TEGEN_AEROSOLS) THEN
    ! With the MACC aerosol climatology we can add in the background
    ! aerosol afterwards using the Tegen arrays.  In this case we
    ! first configure the background aerosol mass-extinction
    ! coefficient at 550 nm
    IF (ITYPE_TROP_BG_AER > 0) THEN
      PRADIATION%TROP_BG_AER_MASS_EXT  = DRY_AEROSOL_MASS_EXTINCTION(RAD_CONFIG,&
           &                                   ITYPE_TROP_BG_AER, 550.0e-9_JPRB)
      WRITE(NULOUT,'(a,i2,a,e12.4,a)') 'Tropospheric background:  aerosol type ',&
           &  ITYPE_TROP_BG_AER, ', 550-nm mass-extinction coefficient ', &
           &  PRADIATION%TROP_BG_AER_MASS_EXT, ' m2 kg-1'
    ELSE
      WRITE(NULOUT,'(a)') 'No tropospheric background aerosol'
    ENDIF

    IF (ITYPE_STRAT_BG_AER > 0) THEN
      PRADIATION%STRAT_BG_AER_MASS_EXT = DRY_AEROSOL_MASS_EXTINCTION(RAD_CONFIG,&
           &                                   ITYPE_STRAT_BG_AER, 550.0e-9_JPRB)
      WRITE(NULOUT,'(a,i2,a,e12.4,a)') 'Stratospheric background: aerosol type ',&
           &  ITYPE_STRAT_BG_AER, ', 550-nm mass-extinction coefficient ', &
           &  PRADIATION%STRAT_BG_AER_MASS_EXT, ' m2 kg-1'
    ELSE
      WRITE(NULOUT,'(a)') 'No stratospheric background aerosol'
    ENDIF
  ENDIF

  IF(YDEAERATM%LAERRAD) THEN
    CALL SETUP_MONO_AER_OPTICS(YAER_RAD_DESC, RAD_CONFIG%AEROSOL_OPTICS, &
         &                     RAD_CONFIG%I_AEROSOL_TYPE_MAP)
  ENDIF

  IF (IVERBOSESETUP > 1) THEN
    WRITE(NULOUT,'(a)') '-------------------------------------------------------------------------------'
  ENDIF

  END ASSOCIATE

  IF (LHOOK) CALL DR_HOOK('SU_RADIATION_SCHEME',1,ZHOOK_HANDLE)

CONTAINS


SUBROUTINE SETUP_MONO_AER_OPTICS(YDESC, AO, KMAP)
  ! Copy the monochromatic aerosol optical properties to the variables
  ! needed by CAMS for aerosol diagnostics and data assimilation
  ! (replaces the old SU_AEROP).
  ! Note that ALF_* are NOT in SI base units [m2 kg-1] as the NetCDF files
  ! are, but [m2 g-1], hence the 1.0e-3 conversion factor.
  
  USE PARKIND1,                      ONLY : JPRB, JPIM
  USE YOMHOOK,                       ONLY : LHOOK, DR_HOOK, JPHOOK
  USE RADIATION_AEROSOL_OPTICS_DATA, ONLY : AEROSOL_OPTICS_TYPE
  USE YOEAERRADDESC,                 ONLY : TAER_RAD_DESC
  USE YOMLUN,                        ONLY : NULOUT, NULERR
  ! For each of the various aerosol species, ALF=mass extinction
  ! coefficient (m2/g), OMG=single-scattering albedo, ASY=asymmetry
  ! factor and RALI=lidar ratio.  The large number of hard-coded types
  ! is rather inflexible and a future modification could be to use a
  ! single set of arrays for the hydrophobic and a single set for the
  ! hydrophilic arrays (with an extra dimension for aerosol type), but
  ! it would require changes to all routines using these arrays, for
  ! which the partitioning here may make sense.
  USE YOEAEROP    ,ONLY : ALF_SU, ALF_OM, ALF_DD, ALF_SS, ALF_BC, ALF_NI, ALF_AM, ALF_SOA, &
                       &  ASY_SU, ASY_OM, ASY_DD, ASY_SS, ASY_BC, ASY_NI, ASY_AM, ASY_SOA, &
                       &  OMG_SU, OMG_OM, OMG_DD, OMG_SS, OMG_BC, OMG_NI, OMG_AM, OMG_SOA, &
                &  RALI_BC, RALI_DD, RALI_OM, RALI_SU, RALI_SS, RALI_NI, RALI_AM, RALI_SOA

  ! Description of each aerosol type passed to the radiation scheme
  TYPE(TAER_RAD_DESC),       INTENT(IN) :: YDESC
  ! Structure containing the aerosol optical properties as read from
  ! the optics file
  TYPE(AEROSOL_OPTICS_TYPE), INTENT(IN) :: AO
  ! Mapping from each aerosol type to be passed to the radiation
  ! scheme, to the aerosol properties from the optics file, where
  ! negative numbers indicate to look-up the hydrophilic rather than
  ! hydrophobic aerosol values
  INTEGER(KIND=JPIM),        INTENT(IN) :: KMAP(:)

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
  
  IF (LHOOK) CALL DR_HOOK('SU_RADIATION_SCHEME:SETUP_MONO_AER_OPTICS',0,ZHOOK_HANDLE)
 
  WRITE(NULOUT,'(a)') 'Setting up monochromatic aerosol optics for diagnostic/DA'

  ! Check the array sizing matches
  IF (SIZE(AO%WAVELENGTH_MONO) /= SIZE(ALF_BC)) THEN
    WRITE(NULERR,'(A,I0,A,I0)') 'SU_RADIATION_SCHEME:SETUP_MONO_AER_OPTICS error: ', SIZE(AO%WAVELENGTH_MONO), &
         &  ' wavelengths requested, arrays hold ', SIZE(ALF_BC)
    CALL ABOR1('SU_RADIATION_SCHEME:SETUP_MONO_AER_OPTICS error')
  ENDIF
  
  ! Desert dust
  IF (YDESC%GET_INDEX('DD',.TRUE.,1) == 0) THEN
    CALL SETUP_MONO_PHOBIC_AER(YDESC,AO,KMAP,'DD',1,ALF_DD(1,:,1),OMG_DD(1,:,1),ASY_DD(1,:,1),RALI_DD(1,:,1))
    CALL SETUP_MONO_PHOBIC_AER(YDESC,AO,KMAP,'DD',2,ALF_DD(1,:,2),OMG_DD(1,:,2),ASY_DD(1,:,2),RALI_DD(1,:,2))
    CALL SETUP_MONO_PHOBIC_AER(YDESC,AO,KMAP,'DD',3,ALF_DD(1,:,3),OMG_DD(1,:,3),ASY_DD(1,:,3),RALI_DD(1,:,3))
    ! Hydrophobic dust: copy values over the rest of the array
    DO JTAB=2,12
      ALF_DD(JTAB,:,1)=ALF_DD(1,:,1)
      ALF_DD(JTAB,:,2)=ALF_DD(1,:,2)
      ALF_DD(JTAB,:,3)=ALF_DD(1,:,3)
      OMG_DD(JTAB,:,1)=OMG_DD(1,:,1)
      OMG_DD(JTAB,:,2)=OMG_DD(1,:,2)
      OMG_DD(JTAB,:,3)=OMG_DD(1,:,3)
      ASY_DD(JTAB,:,1)=ASY_DD(1,:,1)
      ASY_DD(JTAB,:,2)=ASY_DD(1,:,2)
      ASY_DD(JTAB,:,3)=ASY_DD(1,:,3)
      RALI_DD(JTAB,:,1)=RALI_DD(1,:,1)
      RALI_DD(JTAB,:,2)=RALI_DD(1,:,2)
      RALI_DD(JTAB,:,3)=RALI_DD(1,:,3)
    ENDDO
  ELSE
    CALL SETUP_MONO_PHILIC_AER(YDESC,AO,KMAP,'DD',1,ALF_DD(:,:,1),OMG_DD(:,:,1),ASY_DD(:,:,1),RALI_DD(:,:,1))
    CALL SETUP_MONO_PHILIC_AER(YDESC,AO,KMAP,'DD',2,ALF_DD(:,:,2),OMG_DD(:,:,2),ASY_DD(:,:,2),RALI_DD(:,:,2))
    CALL SETUP_MONO_PHILIC_AER(YDESC,AO,KMAP,'DD',3,ALF_DD(:,:,3),OMG_DD(:,:,3),ASY_DD(:,:,3),RALI_DD(:,:,3))
  ENDIF
  ! Black carbon
  CALL SETUP_MONO_PHOBIC_AER(YDESC,AO,KMAP,'BC',0,ALF_BC,OMG_BC,ASY_BC,RALI_BC)
  ! Organic matter
  CALL SETUP_MONO_PHILIC_AER(YDESC,AO,KMAP,'OM',0,ALF_OM,OMG_OM,ASY_OM,RALI_OM)
  ! Sulfates
  CALL SETUP_MONO_PHILIC_AER(YDESC,AO,KMAP,'SU',0,ALF_SU,OMG_SU,ASY_SU,RALI_SU)
  ! Sea salt
  CALL SETUP_MONO_PHILIC_AER(YDESC,AO,KMAP,'SS',1,ALF_SS(:,:,1),OMG_SS(:,:,1),ASY_SS(:,:,1),RALI_SS(:,:,1))
  CALL SETUP_MONO_PHILIC_AER(YDESC,AO,KMAP,'SS',2,ALF_SS(:,:,2),OMG_SS(:,:,2),ASY_SS(:,:,2),RALI_SS(:,:,2))
  CALL SETUP_MONO_PHILIC_AER(YDESC,AO,KMAP,'SS',3,ALF_SS(:,:,3),OMG_SS(:,:,3),ASY_SS(:,:,3),RALI_SS(:,:,3))
  ! Nitrates
  CALL SETUP_MONO_PHILIC_AER(YDESC,AO,KMAP,'NI',1,ALF_NI(:,:,1),OMG_NI(:,:,1),ASY_NI(:,:,1),RALI_NI(:,:,1))
  CALL SETUP_MONO_PHILIC_AER(YDESC,AO,KMAP,'NI',2,ALF_NI(:,:,2),OMG_NI(:,:,2),ASY_NI(:,:,2),RALI_NI(:,:,2))
  ! Ammonium
  CALL SETUP_MONO_PHILIC_AER(YDESC,AO,KMAP,'AM',0,ALF_AM,OMG_AM,ASY_AM,RALI_AM)
  ! Secondary organics (biogenic and anthropogenic)
  CALL SETUP_MONO_PHILIC_AER(YDESC,AO,KMAP,'OB',0,ALF_SOA(:,:,1),OMG_SOA(:,:,1),ASY_SOA(:,:,1),RALI_SOA(:,:,1))
  CALL SETUP_MONO_PHILIC_AER(YDESC,AO,KMAP,'OA',0,ALF_SOA(:,:,2),OMG_SOA(:,:,2),ASY_SOA(:,:,2),RALI_SOA(:,:,2))
  
  IF (LHOOK) CALL DR_HOOK('SU_RADIATION_SCHEME:SETUP_MONO_AER_OPTICS',1,ZHOOK_HANDLE)

END SUBROUTINE SETUP_MONO_AER_OPTICS

SUBROUTINE SETUP_MONO_PHOBIC_AER( YAER_RAD_DESC, AO, KMAP, CDAER_TYPE, KBIN &
                              & , PMASS_EXT, PSSA, PASYMMETRY, PRALI)
  ! Setup the monochromatic optical properties for one particular
  ! hydrophobic aerosol species

  USE PARKIND1,                      ONLY : JPRB,JPIM
  USE YOMHOOK,                       ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOEAERRADDESC,                 ONLY : TAER_RAD_DESC
  USE RADIATION_AEROSOL_OPTICS_DATA, ONLY : AEROSOL_OPTICS_TYPE
  USE YOMLUN,                        ONLY : NULOUT, NULERR

  TYPE(TAER_RAD_DESC),       INTENT(IN) :: YAER_RAD_DESC
  TYPE(AEROSOL_OPTICS_TYPE), INTENT(IN) :: AO
  INTEGER(KIND=JPIM),        INTENT(IN) :: KMAP(:)
  CHARACTER(2),              INTENT(IN) :: CDAER_TYPE
  INTEGER(KIND=JPIM),        INTENT(IN) :: KBIN
  ! Mass-extinction coefficient (m2/g), single-scattering albedo,
  ! asymmetry factor and lidar ratio
  REAL(KIND=JPRB), INTENT(INOUT), DIMENSION(:) :: PMASS_EXT, PSSA, PASYMMETRY, PRALI
  
  INTEGER(KIND=JPIM) :: IAER

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

  IF (LHOOK) CALL DR_HOOK('SU_RADIATION_SCHEME:SETUP_MONO_PHOBIC_AER',0,ZHOOK_HANDLE)

  IAER = YAER_RAD_DESC%GET_INDEX(CDAER_TYPE, .FALSE., KBIN)
  IF (IAER > 0) THEN
    IAER = KMAP(IAER)
    WRITE(NULOUT,'(A,A2,A,I0,A,I0)') '  Hydrophobic ', CDAER_TYPE, ' bin ', KBIN, ' -> ', IAER
    IF (IAER <= 0) THEN
      WRITE(NULERR,'(A)') 'SU_RADIATION_SCHEME:SETUP_MONO_PHOBIC_AER error: aerosol optics bin is negative'
      CALL ABOR1('SU_RADIATION_SCHEME:SETUP_MONO_PHOBIC_AER error')
    ENDIF
   ! TODO check unit consistency and add conversion unit based on metadata of netcdf file(?)    
    PMASS_EXT = AO%MASS_EXT_MONO_PHOBIC(:,IAER) * 1.0e-3_JPRB
    PSSA      = AO%SSA_MONO_PHOBIC(:,IAER)
    PASYMMETRY= AO%G_MONO_PHOBIC(:,IAER)
    PRALI     = AO%LIDAR_RATIO_MONO_PHOBIC(:,IAER)
  ELSE
    WRITE(NULOUT,'(A,A2,A,I0,A)') '  Hydrophobic ', CDAER_TYPE, ' bin ', KBIN, ' NOT FOUND'
    ! Set mass-extinction coefficient to zero and others to plausible
    ! values
    PMASS_EXT = 0.0_JPRB
    PSSA      = 0.9_JPRB
    PASYMMETRY= 0.6_JPRB
    PRALI     = 1.0_JPRB
  ENDIF

  IF (LHOOK) CALL DR_HOOK('SU_RADIATION_SCHEME:SETUP_MONO_PHOBIC_AER',1,ZHOOK_HANDLE)

END SUBROUTINE SETUP_MONO_PHOBIC_AER
  
SUBROUTINE SETUP_MONO_PHILIC_AER( YAER_RAD_DESC, AO, KMAP, CDAER_TYPE, KBIN &
                               &, PMASS_EXT, PSSA, PASYMMETRY, PRALI)
  ! Setup the monochromatic optical properties for one particular
  ! hydrophilic aerosol species
  
  USE PARKIND1,                      ONLY : JPRB,JPIM
  USE YOMHOOK,                       ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOEAERRADDESC,                 ONLY : TAER_RAD_DESC
  USE RADIATION_AEROSOL_OPTICS_DATA, ONLY : AEROSOL_OPTICS_TYPE
  USE YOMLUN,                        ONLY : NULOUT, NULERR

  TYPE(TAER_RAD_DESC),       INTENT(IN) :: YAER_RAD_DESC
  TYPE(AEROSOL_OPTICS_TYPE), INTENT(IN) :: AO
  INTEGER(KIND=JPIM),        INTENT(IN) :: KMAP(:)
  CHARACTER(2),              INTENT(IN) :: CDAER_TYPE
  INTEGER(KIND=JPIM),        INTENT(IN) :: KBIN
  ! Mass-extinction coefficient (m2/g), single-scattering albedo,
  ! asymmetry factor and lidar ratio
  REAL(KIND=JPRB), INTENT(INOUT), DIMENSION(:,:) :: PMASS_EXT, PSSA, PASYMMETRY, PRALI
  
  INTEGER(KIND=JPIM) :: IAER

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

  IF (LHOOK) CALL DR_HOOK('SU_RADIATION_SCHEME:SETUP_MONO_PHILIC_AER',0,ZHOOK_HANDLE)

  IAER = YAER_RAD_DESC%GET_INDEX(CDAER_TYPE, .TRUE., KBIN)
  IF (IAER > 0) THEN
    ! Mapping is negative for hydrophilic types, so needs negating to
    ! create an index
    IAER = -KMAP(IAER)
    WRITE(NULOUT,'(A,A2,A,I0,A,I0)') '  Hydrophilic ', CDAER_TYPE, ' bin ', KBIN, ' -> ', IAER
    IF (IAER <= 0) THEN
      WRITE(NULERR,'(A)') 'SU_RADIATION_SCHEME:SETUP_MONO_PHILIC_AER error: aerosol optics bin is negative'
      CALL ABOR1('SU_RADIATION_SCHEME:SETUP_MONO_PHILIC_AER error')
    ENDIF
   ! TODO check unit consistency and add conversion unit based on metadata of netcdf file(?)    
    PMASS_EXT = TRANSPOSE(AO%MASS_EXT_MONO_PHILIC(:,:,IAER)) * 1.0e-3_JPRB
    PSSA      = TRANSPOSE(AO%SSA_MONO_PHILIC(:,:,IAER))
    PASYMMETRY= TRANSPOSE(AO%G_MONO_PHILIC(:,:,IAER))
    PRALI     = TRANSPOSE(AO%LIDAR_RATIO_MONO_PHILIC(:,:,IAER))
  ELSE
    WRITE(NULOUT,'(A,A2,A,I0,A)') '  Hydrophilic ', CDAER_TYPE, ' bin ', KBIN, ' NOT FOUND'
    ! Set mass-extinction coefficient to zero and others to plausible
    ! values
    PMASS_EXT = 0.0_JPRB
    PSSA      = 0.9_JPRB
    PASYMMETRY= 0.6_JPRB
    PRALI     = 1.0_JPRB
  ENDIF

  IF (LHOOK) CALL DR_HOOK('SU_RADIATION_SCHEME:SETUP_MONO_PHILIC_AER',1,ZHOOK_HANDLE)

END SUBROUTINE SETUP_MONO_PHILIC_AER

END SUBROUTINE SU_RADIATION_SCHEME



