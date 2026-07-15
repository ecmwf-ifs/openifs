! radiation_config.F90 - Derived type to configure the radiation scheme
!
! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! Modifications
!   2017-07-22  R. Hogan  Added Yi et al. ice optics model
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-03-15  R. Hogan  Added logicals controlling surface spectral treatment
!   2018-08-29  R. Hogan  Added monochromatic single-scattering albedo / asymmetry factor
!   2018-09-03  R. Hogan  Added min_cloud_effective_size
!   2018-09-04  R. Hogan  Added encroachment_scaling
!   2018-09-13  R. Hogan  Added IEncroachmentFractal
!   2019-01-02  R. Hogan  Added Cloudless solvers
!   2019-01-14  R. Hogan  Added out_of_bounds_[1,2,3]d for checker routines
!   2019-01-18  R. Hogan  Added albedo weighting
!   2019-02-03  R. Hogan  Added ability to fix out-of-physical-bounds inputs
!   2019-02-10  R. Hogan  Renamed "encroachment" to "entrapment"
!   2020-05-18  R. Hogan  Moved out_of_bounds_* to radiation_check.F90
!   2021-07-04  R. Hogan  Numerous changes for ecCKD and general cloud/aerosol optics
!
! Note: The aim is for ecRad in the IFS to be as similar as possible
! to the offline version, so if you make any changes to this or any
! files in this directory, please inform Robin Hogan.
!

module radiation_config

  use parkind1,                      only : jprb

  use radiation_cloud_optics_data,   only : cloud_optics_type
  use radiation_general_cloud_optics_data,   only : general_cloud_optics_type
  use radiation_aerosol_optics_data, only : aerosol_optics_type
  use radiation_pdf_sampler,         only : pdf_sampler_type
  use radiation_cloud_cover,         only : OverlapName, &
       & IOverlapMaximumRandom, IOverlapExponentialRandom, IOverlapExponential
  use radiation_ecckd,               only : ckd_model_type

  implicit none
  public

  ! Configuration codes: use C-style enumerators to avoid having to
  ! remember the numbers

  ! Solvers: can be specified for longwave and shortwave
  ! independently, except for "Homogeneous", which must be the same
  ! for both
  enum, bind(c) 
     enumerator ISolverCloudless, ISolverHomogeneous, ISolverMcICA, &
          &     ISolverSpartacus, ISolverTripleclouds 
  end enum
  character(len=*), parameter :: SolverName(0:4) = (/ 'Cloudless   ', &
       &                                              'Homogeneous ', &
       &                                              'McICA       ', &
       &                                              'SPARTACUS   ', &
       &                                              'Tripleclouds' /)

  ! SPARTACUS shortwave solver can treat the reflection of radiation
  ! back up into different regions in various ways
  enum, bind(c) 
     enumerator &
       & IEntrapmentZero, &     ! No entrapment, as Tripleclouds
       & IEntrapmentEdgeOnly, & ! Only radiation passed through cloud edge is horizontally homogenized
       & IEntrapmentExplicit, & ! Estimate horiz migration dist, account for fractal clouds
       & IEntrapmentExplicitNonFractal, & ! As above but ignore fractal nature of clouds
       & IEntrapmentMaximum ! Complete horizontal homogenization within regions (old SPARTACUS assumption)
  end enum
  
  ! Names available in the radiation namelist for variable
  ! sw_entrapment_name
  character(len=*), parameter :: EntrapmentName(0:4)   = [ 'Zero       ', &
       &                                                   'Edge-only  ', &
       &                                                   'Explicit   ', &
       &                                                   'Non-fractal', &
       &                                                   'Maximum    ' ]
  ! For backwards compatibility, the radiation namelist also supports
  ! the equivalent variable sw_encroachment_name with the following
  ! names
  character(len=*), parameter :: EncroachmentName(0:4) = [ 'Zero    ', &
       &                                                   'Minimum ', &
       &                                                   'Fractal ', &
       &                                                   'Computed', &
       &                                                   'Maximum ' ]

  ! Two-stream models
  ! This is not configurable at run-time

  ! Gas models
  enum, bind(c) 
     enumerator IGasModelMonochromatic, IGasModelIFSRRTMG, IGasModelECCKD
  end enum
  character(len=*), parameter :: GasModelName(0:2) = (/ 'Monochromatic', &
       &                                                'RRTMG-IFS    ', &
       &                                                'ECCKD        '/)

  ! Hydrometeor scattering models
  enum, bind(c) 
     enumerator ILiquidModelMonochromatic, &
          &     ILiquidModelSOCRATES, ILiquidModelSlingo
  end enum
  character(len=*), parameter :: LiquidModelName(0:2) = (/ 'Monochromatic', &
       &                                                   'SOCRATES     ', &
       &                                                   'Slingo       ' /)

  enum, bind(c) 
     enumerator IIceModelMonochromatic, IIceModelFu, &
          &  IIceModelBaran, IIceModelBaran2016, IIceModelBaran2017,   &
          &  IIceModelYi
  end enum
  character(len=*), parameter :: IceModelName(0:5) = (/ 'Monochromatic', &
       &                                                'Fu-IFS       ', &
       &                                                'Baran        ', &
       &                                                'Baran2016    ', &
       &                                                'Baran2017    ', &
       &                                                'Yi           ' /)

  ! Cloud PDF distribution shapes
  enum, bind(c)
    enumerator IPdfShapeLognormal, IPdfShapeGamma
  end enum
  character(len=*), parameter :: PdfShapeName(0:1) = (/ 'Lognormal', &
       &                                                'Gamma    ' /)

  ! Maximum number of different aerosol types that can be provided
  integer, parameter :: NMaxAerosolTypes = 256

  ! Maximum number of different cloud types that can be provided
  integer, parameter :: NMaxCloudTypes = 12

  ! Maximum number of shortwave albedo and longwave emissivity
  ! intervals
  integer, parameter :: NMaxAlbedoIntervals = 256

  ! Length of string buffer for printing config information
  integer, parameter :: NPrintStringLen = 60

  !---------------------------------------------------------------------
  ! Derived type containing all the configuration information needed
  ! to run the radiation scheme.  The intention is that this is fixed
  ! for a given model run.  The parameters are to list first those
  ! quantities that can be set directly by the user, for example using a
  ! namelist, and second those quantities that are computed afterwards
  ! from the user-supplied numbers, especially the details of the gas
  ! optics model.
  type config_type
    ! USER-CONFIGURABLE PARAMETERS

    ! Scale the solar spectrum per band (or g-point if
    ! do_cloud_aerosol_per_sw_g_point=true) via vector
    ! single_level%spectral_solar_scaling
    logical :: use_spectral_solar_scaling = .false.

    ! Modify the solar spectrum per g-point to account for the current
    ! phase of the solar cycle, via scalar
    ! single_level%spectral_solar_cycle_multiplier
    logical :: use_spectral_solar_cycle = .false.
    
    ! Directory in which gas, cloud and aerosol data files are to be
    ! found
    character(len=511) :: directory_name = '.'

    ! If this is true then support arbitrary hydrometeor types (not
    ! just ice and liquid) and arbitrary spectral discretization (not
    ! just RRTMG). It is required that this is true if the ecCKD gas
    ! optics model is selected. General cloud optics has only been
    ! available from ecRad version 1.5.
    logical :: use_general_cloud_optics = .true.

    ! If this is true then support aerosol properties at an arbitrary
    ! spectral discretization (not just RRTMG). It is required that
    ! this is true if the ecCKD gas optics model is selected.
    logical :: use_general_aerosol_optics = .true.

    ! Cloud is deemed to be present in a layer if cloud fraction
    ! exceeds this value
    real(jprb) :: cloud_fraction_threshold = 1.0e-6_jprb
    ! ...and total cloud water mixing ratio exceeds this value
    real(jprb) :: cloud_mixing_ratio_threshold = 1.0e-9_jprb

    ! Overlap scheme
    integer :: i_overlap_scheme = IOverlapExponentialRandom

    ! Use the Shonk et al. (2010) "beta" overlap parameter, rather
    ! than the "alpha" overlap parameter of Hogan and Illingworth
    ! (2000)?
    logical :: use_beta_overlap = .false.

    ! Use a more vectorizable McICA cloud generator, at the expense of
    ! more random numbers being generated?  This is the default on NEC
    ! SX.
#ifdef __SX__
    logical :: use_vectorizable_generator = .true.
#else
    logical :: use_vectorizable_generator = .false.
#endif

    ! Shape of sub-grid cloud water PDF
    integer :: i_cloud_pdf_shape = IPdfShapeGamma

    ! The ratio of the overlap decorrelation length for cloud
    ! inhomogeneities to the overlap decorrelation length for cloud
    ! boundaries.  Observations suggest this has a value of 0.5
    ! (e.g. from the decorrelation lengths of Hogan and Illingworth
    ! 2003 and Hogan and Illingworth 2000).
    real(jprb) :: cloud_inhom_decorr_scaling = 0.5_jprb

    ! Factor controlling how much of the cloud edge length interfaces
    ! directly between the clear-sky region (region a) and the
    ! optically thick cloudy region (region c).  If Lxy is the length
    ! of the interfaces between regions x and y, and Lab and Lbc have
    ! been computed already, then
    !   Lac=clear_to_thick_fraction*min(Lab,Lbc).
    real(jprb) :: clear_to_thick_fraction = 0.0_jprb

    ! Factor allowing lateral transport when the sun is close to
    ! overhead; consider atand(overhead_sun_factor) to be the number
    ! of degrees that the sun angle is perturbed from zenith for the
    ! purposes of computing lateral transport.  A value of up to 0.1
    ! seems to be necessary to account for the fact that some forward
    ! scattered radiation is treated as unscattered by delta-Eddington
    ! scaling; therefore it ought to have the chance to escape.
    real(jprb) :: overhead_sun_factor = 0.0_jprb

    ! Minimum gas optical depth in a single layer at any wavelength,
    ! for stability
    real(jprb) :: min_gas_od_lw = 1.0e-15_jprb
    real(jprb) :: min_gas_od_sw = 0.0_jprb

    ! Maximum gas optical depth in a layer before that g-point will
    ! not be considered for 3D treatment: a limit is required to avoid
    ! expensive computation of matrix exponentials on matrices with
    ! large elements
    real(jprb) :: max_gas_od_3d = 8.0_jprb

    ! Maximum total optical depth of a cloudy region for stability:
    ! optical depth will be capped at this value in the SPARTACUS
    ! solvers
    real(jprb) :: max_cloud_od = 16.0_jprb

    ! How much longwave scattering is included?
    logical :: do_lw_cloud_scattering = .true.
    logical :: do_lw_aerosol_scattering = .true.

    ! Number of regions used to describe clouds and clear skies. A
    ! value of 2 means one clear and one cloudy region, so clouds are
    ! horizontally homogeneous, while a value of 3 means two cloudy
    ! regions with different optical depth, thereby representing
    ! inhomogeneity via the Shonk & Hogan (2008) "Tripleclouds"
    ! method.
    integer :: nregions = 3

    ! Code specifying the solver to be used: use the enumerations
    ! defined above
    integer :: i_solver_sw = ISolverMcICA
    integer :: i_solver_lw = ISolverMcICA

    ! Do shortwave delta-Eddington scaling on the cloud-aerosol-gas
    ! mixture (as in the original IFS scheme), rather than the more
    ! correct approach of separately scaling the cloud and aerosol
    ! scattering properties before merging with gases.  Note that
    ! .true. is not compatible with the SPARTACUS solver.
    logical :: do_sw_delta_scaling_with_gases = .false.

    ! Codes describing the gas model
    integer :: i_gas_model = IGasModelIFSRRTMG

    ! Optics if i_gas_model==IGasModelMonochromatic.
    ! The wavelength to use for the Planck function in metres. If this
    ! is positive then the output longwave fluxes will be in units of
    ! W m-2 um-1.  If this is zero or negative (the default) then
    ! sigma*T^4 will be used and the output longwave fluxes will be in
    ! W m-2.
    real(jprb) :: mono_lw_wavelength = -1.0_jprb
    ! Total zenith optical depth of the atmosphere in the longwave and
    ! shortwave, distributed vertically according to the pressure.
    ! Default is zero.
    real(jprb) :: mono_lw_total_od = 0.0_jprb
    real(jprb) :: mono_sw_total_od = 0.0_jprb
    ! Single-scattering albedo and asymmetry factor: values typical
    ! for liquid clouds with effective radius of 10 microns, at (SW)
    ! 0.55 micron wavelength and (LW) 10.7 microns wavelength
    real(jprb) :: mono_sw_single_scattering_albedo = 0.999999_jprb
    real(jprb) :: mono_sw_asymmetry_factor = 0.86_jprb
    real(jprb) :: mono_lw_single_scattering_albedo = 0.538_jprb
    real(jprb) :: mono_lw_asymmetry_factor = 0.925_jprb

    ! Codes describing particle scattering models
    integer :: i_liq_model = ILiquidModelSOCRATES
    integer :: i_ice_model = IIceModelBaran
    
    ! The mapping from albedo/emissivity intervals to SW/LW bands can
    ! either be done by finding the interval containing the central
    ! wavenumber of the band (nearest neighbour), or by a weighting
    ! according to the spectral overlap of each interval with each
    ! band
    logical :: do_nearest_spectral_sw_albedo = .false.
    logical :: do_nearest_spectral_lw_emiss  = .false.

    ! User-defined monotonically increasing wavelength bounds (m)
    ! between input surface albedo/emissivity intervals. Implicitly
    ! the first interval starts at zero and the last ends at
    ! infinity. These must be set with define_sw_albedo_intervals and
    ! define_lw_emiss_intervals.
    real(jprb) :: sw_albedo_wavelength_bound(NMaxAlbedoIntervals-1) = -1.0_jprb
    real(jprb) :: lw_emiss_wavelength_bound( NMaxAlbedoIntervals-1) = -1.0_jprb

    ! The index to the surface albedo/emissivity intervals for each of
    ! the wavelength bounds specified in sw_albedo_wavelength_bound
    ! and lw_emiss_wavelength_bound
    integer :: i_sw_albedo_index(NMaxAlbedoIntervals) = 0
    integer :: i_lw_emiss_index (NMaxAlbedoIntervals) = 0

    ! Do we compute longwave and/or shortwave radiation?
    logical :: do_lw = .true.
    logical :: do_sw = .true.

    ! Do we compute clear-sky fluxes and/or solar direct fluxes?
    logical :: do_clear = .true.
    logical :: do_sw_direct = .true.

    ! Do we include 3D effects?
    logical :: do_3d_effects = .true.
    
    character(len=511) :: cloud_type_name(NMaxCloudTypes) = ["","","","","","","","","","","",""]
! &
!         &   = ["mie_droplet                   ", &
!         &      "baum-general-habit-mixture_ice"]

    ! Spectral averaging method to use with generalized cloud optics;
    ! see Edwards & Slingo (1996) for definition.  Experimentation
    ! with ecRad suggests that "thick" averaging is more accurate for
    ! both liquid and ice clouds.
    logical :: use_thick_cloud_spectral_averaging(NMaxCloudTypes) &
         &  = [.true.,.true.,.true.,.true.,.true.,.true., &
         &     .true.,.true.,.true.,.true.,.true.,.true.]

    ! To what extent do we include "entrapment" effects in the
    ! SPARTACUS solver? This essentially means that in a situation
    ! like this
    !
    ! 000111
    ! 222222
    !
    ! Radiation downwelling from region 1 may be reflected back into
    ! region 0 due to some degree of homogenization of the radiation
    ! in region 2.  Hogan and Shonk (2013) referred to this as
    ! "anomalous horizontal transport" for a 1D model, although for 3D
    ! calculations it is desirable to include at least some of it. The
    ! options are described by the IEntrapment* parameters above.
    integer :: i_3d_sw_entrapment = IEntrapmentExplicit

    ! In the longwave, the equivalent process it either "on" (like
    ! maximum entrapment) or "off" (like zero entrapment):
    logical :: do_3d_lw_multilayer_effects = .false.

    ! Do we account for the effective emissivity of the side of
    ! clouds?
    logical :: do_lw_side_emissivity = .true.

    ! The 3D transfer rate "X" is such that if transport out of a
    ! region was the only process occurring then by the base of a
    ! layer only exp(-X) of the original flux would remain in that
    ! region. The transfer rate computed geometrically can be very
    ! high for the clear-sky regions in layers with high cloud
    ! fraction.  For stability reasons it is necessary to provide a
    ! maximum possible 3D transfer rate.
    real(jprb) :: max_3d_transfer_rate = 10.0_jprb

    ! It has also sometimes been found necessary to set a minimum
    ! cloud effective size for stability (metres)
    real(jprb) :: min_cloud_effective_size = 100.0_jprb

    ! Given a horizontal migration distance, there is still
    ! uncertainty about how much entrapment occurs associated with how
    ! one assumes cloud boundaries line up in adjacent layers. This
    ! factor can be varied between 0.0 (the boundaries line up to the
    ! greatest extent possible given the overlap parameter) and 1.0
    ! (the boundaries line up to the minimum extent possible). In the
    ! Hogan et al. entrapment paper it is referred to as the overhang
    ! factor zeta, and a value of 0 matches the Monte Carlo
    ! calculations best.
    real(jprb) :: overhang_factor = 0.0_jprb

    ! By default, the Meador & Weaver (1980) expressions are used
    ! instead of the matrix exponential whenever 3D effects can be
    ! neglected (e.g. cloud-free layers or clouds with infinitely
    ! large effective cloud size), but setting the following to true
    ! uses the matrix exponential everywhere, enabling the two
    ! methods to be compared. Note that Meador & Weaver will still be
    ! used for very optically thick g points where the matrix
    ! exponential can produce incorrect results.
    logical :: use_expm_everywhere = .false.

    ! Aerosol descriptors: aerosol_type_mapping must be of length
    ! n_aerosol_types, and contains 0 if that type is to be ignored,
    ! positive numbers to map on to the indices of hydrophobic
    ! aerosols in the aerosol optics configuration file, and negative
    ! numbers to map on to (the negative of) the indices of
    ! hydrophilic aerosols in the configuration file.
    logical :: use_aerosols = .false.
    integer :: n_aerosol_types = 0
    integer :: i_aerosol_type_map(NMaxAerosolTypes)

    ! Save the gas and cloud optical properties for each g point in
    ! "radiative_properties.nc"?
    logical :: do_save_radiative_properties = .false.

    ! Save the flux profiles in each band?
    logical :: do_save_spectral_flux = .false.

    ! Save the surface downwelling shortwave fluxes in each band?
    logical :: do_surface_sw_spectral_flux = .true.

    ! Save the TOA fluxes in each band?
    logical :: do_toa_spectral_flux = .false.

    ! Compute the longwave derivatives needed to apply the approximate
    ! radiation updates of Hogan and Bozzo (2015)
    logical :: do_lw_derivatives = .false.

    ! Save the flux profiles in each g-point (overrides
    ! do_save_spectral_flux if TRUE)?
    logical :: do_save_gpoint_flux = .false.

    ! In the IFS environment, setting up RRTM has already been done
    ! so not needed to do it again
    logical :: do_setup_ifsrrtm = .true.

    ! In the IFS environment the old scheme has a bug in the Fu
    ! longwave ice optics whereby the single scattering albedo is one
    ! minus what it should be.  Unfortunately fixing it makes
    ! forecasts worse. Setting the following to true reproduces the
    ! bug.
    logical :: do_fu_lw_ice_optics_bug = .false.

    ! Control verbosity: 0=none (no output to standard output; write
    ! to standard error only if an error occurs), 1=warning, 2=info,
    ! 3=progress, 4=detailed, 5=debug.  Separate settings for the
    ! setup of the scheme and the execution of it.
    integer :: iverbosesetup = 2
    integer :: iverbose = 1

    ! Are we doing radiative transfer in complex surface canopies
    ! (streets/vegetation), in which case tailored downward fluxes are
    ! needed at the top of the canopy?
    logical :: do_canopy_fluxes_sw = .false.
    logical :: do_canopy_fluxes_lw = .false.
    ! If so, do we use the full spectrum as in the atmosphere, or just
    ! the reduced spectrum in which the shortwave albedo and longwave
    ! emissivity are provided?
    logical :: use_canopy_full_spectrum_sw = .false.
    logical :: use_canopy_full_spectrum_lw = .false.
    ! Do we treat gas radiative transfer in streets/vegetation?
    logical :: do_canopy_gases_sw = .false.
    logical :: do_canopy_gases_lw = .false.

    ! Optics file names for overriding the ones generated from the
    ! other options. If these remain empty then the generated names
    ! will be used (see the "consolidate_config" routine below). If
    ! the user assigns one of these and it starts with a '/' character
    ! then that will be used instead. If the user assigns one and it
    ! doesn't start with a '/' character then it will be prepended by
    ! the contents of directory_name.
    character(len=511) :: ice_optics_override_file_name     = ''
    character(len=511) :: liq_optics_override_file_name     = ''
    character(len=511) :: aerosol_optics_override_file_name = ''
    character(len=511) :: gas_optics_sw_override_file_name  = ''
    character(len=511) :: gas_optics_lw_override_file_name  = ''

    ! Optionally override the default file describing variations in
    ! the spectral solar irradiance associated with the solar cycle
    character(len=511) :: ssi_override_file_name = ''

    ! Do we use the solar spectral irradiance file to update the solar
    ! irradiance in each g point? Only possible if
    ! use_spectral_solar_cycle==true.
    logical :: use_updated_solar_spectrum = .false.
    
    ! Optionally override the look-up table file for the cloud-water
    ! PDF used by the McICA solver
    character(len=511) :: cloud_pdf_override_file_name = ''

    ! Do we compute cloud, aerosol and surface optical properties per
    ! g point?  Not available with RRTMG gas optics model.
    logical :: do_cloud_aerosol_per_sw_g_point = .true.
    logical :: do_cloud_aerosol_per_lw_g_point = .true.

    ! Do we weight the mapping from surface emissivity/albedo to
    ! g-point/band weighting by a reference Planck function (more
    ! accurate) or weighting each wavenumber equally (less accurate
    ! but consistent with IFS Cycle 48r1 and earlier)?
    logical :: do_weighted_surface_mapping = .true.

    ! COMPUTED PARAMETERS

    ! Users of this library should not edit these parameters directly;
    ! they are set by the "consolidate" routine

    ! Has "consolidate" been called?  
    logical :: is_consolidated = .false.

    ! Fraction of each g point in each wavenumber interval,
    ! dimensioned (n_wav_frac_[l|s]w, n_g_[l|s]w)
    real(jprb), allocatable, dimension(:,:) :: g_frac_sw, g_frac_lw

    ! If the nearest surface albedo/emissivity interval is to be used
    ! for each SW/LW band then the following arrays will be allocated
    ! to the length of the number of bands and contain the index to
    ! the relevant interval
    integer, allocatable, dimension(:) :: i_albedo_from_band_sw
    integer, allocatable, dimension(:) :: i_emiss_from_band_lw

    ! ...alternatively, this matrix dimensioned
    ! (n_albedo_intervals,n_bands_sw) providing the weights needed for
    ! computing the albedo in each ecRad band from the albedo in each
    ! native albedo band - see radiation_single_level.F90
    real(jprb), allocatable, dimension(:,:) :: sw_albedo_weights
    ! ...and similarly in the longwave, dimensioned
    ! (n_emiss_intervals,n_bands_lw)
    real(jprb), allocatable, dimension(:,:) :: lw_emiss_weights

    ! Arrays of length the number of g-points that convert from
    ! g-point to the band index
    integer, allocatable, dimension(:) :: i_band_from_g_lw
    integer, allocatable, dimension(:) :: i_band_from_g_sw

    ! We allow for the possibility for g-points to be ordered in terms
    ! of likely absorption (weakest to strongest) across the shortwave
    ! or longwave spectrum, in order that in SPARTACUS we select only
    ! the first n g-points that will not have too large an absorption,
    ! and therefore matrix exponentials that are both finite and not
    ! too expensive to compute.  The following two arrays map the
    ! reordered g-points to the original ones.
    integer, allocatable, dimension(:) :: i_g_from_reordered_g_lw
    integer, allocatable, dimension(:) :: i_g_from_reordered_g_sw

    ! The following map the reordered g-points to the bands
    integer, allocatable, dimension(:) :: i_band_from_reordered_g_lw
    integer, allocatable, dimension(:) :: i_band_from_reordered_g_sw

    ! The following map the reordered g-points to the spectral
    ! information being saved: if do_save_gpoint_flux==TRUE then this
    ! will map on to the original g points, but if only
    ! do_save_spectral_flux==TRUE then this will map on to the bands
    integer, pointer, dimension(:) :: i_spec_from_reordered_g_lw
    integer, pointer, dimension(:) :: i_spec_from_reordered_g_sw

    ! Number of spectral intervals used for the canopy radiative
    ! transfer calculation; they are either equal to
    ! n_albedo_intervals/n_emiss_intervals or n_g_sw/n_g_lw
    integer :: n_canopy_bands_sw = 1
    integer :: n_canopy_bands_lw = 1

    ! Data structures containing gas optics description in the case of
    ! ecCKD
    type(ckd_model_type)         :: gas_optics_sw, gas_optics_lw

    ! Data structure containing cloud scattering data
    type(cloud_optics_type)      :: cloud_optics

    ! Number of general cloud types, default liquid and ice
    integer :: n_cloud_types = 2

    ! List of data structures (one per cloud type) containing cloud
    ! scattering data
    type(general_cloud_optics_type), allocatable :: cloud_optics_sw(:)
    type(general_cloud_optics_type), allocatable :: cloud_optics_lw(:)

    ! Data structure containing aerosol scattering data
    type(aerosol_optics_type)    :: aerosol_optics

    ! Object for sampling from a gamma or lognormal distribution
    type(pdf_sampler_type)       :: pdf_sampler

    ! Optics file names
    character(len=511) :: ice_optics_file_name, &
         &                liq_optics_file_name, &
         &                aerosol_optics_file_name, &
         &                gas_optics_sw_file_name, &
         &                gas_optics_lw_file_name

    ! Solar spectral irradiance file name
    character(len=511) :: ssi_file_name
    
    ! McICA PDF look-up table file name
    character(len=511) :: cloud_pdf_file_name

    ! Number of gpoints and bands in the shortwave and longwave - set
    ! to zero as will be set properly later
    integer :: n_g_sw = 0, n_g_lw = 0
    integer :: n_bands_sw = 0, n_bands_lw = 0

    ! Number of spectral points to save (equal either to the number of
    ! g points or the number of bands
    integer :: n_spec_sw = 0, n_spec_lw = 0

    ! Number of wavenumber intervals used to describe the mapping from
    ! g-points to wavenumber space
    integer :: n_wav_frac_sw = 0, n_wav_frac_lw = 0

    ! Dimensions to store variables that are only needed if longwave
    ! scattering is included. "n_g_lw_if_scattering" is equal to
    ! "n_g_lw" if aerosols are allowed to scatter in the longwave,
    ! and zero otherwise. "n_bands_lw_if_scattering" is equal to
    ! "n_bands_lw" if clouds are allowed to scatter in the longwave,
    ! and zero otherwise.
    integer :: n_g_lw_if_scattering = 0, n_bands_lw_if_scattering = 0

    ! Treat clouds as horizontally homogeneous within the gribox
    logical :: is_homogeneous = .false.

    ! If the solvers are both "Cloudless" then we don't need to do any
    ! cloud processing
    logical :: do_clouds = .true.

   contains
     procedure :: read => read_config_from_namelist
     procedure :: consolidate => consolidate_config
     procedure :: set  => set_config
     procedure :: print => print_config
     procedure :: get_sw_weights
     procedure :: get_sw_mapping
     procedure :: define_sw_albedo_intervals
     procedure :: define_lw_emiss_intervals
     procedure :: set_aerosol_wavelength_mono
     procedure :: consolidate_sw_albedo_intervals
     procedure :: consolidate_lw_emiss_intervals

  end type config_type

!  procedure, private :: print_logical, print_real, print_int

contains


  !---------------------------------------------------------------------
  ! This subroutine reads configuration data from a namelist file, and
  ! anything that is not in the namelists will be set to default
  ! values. If optional output argument "is_success" is present, then
  ! on error (e.g. missing file) it will be set to .false.; if this
  ! argument is missing then on error the program will be aborted. You
  ! may either specify the file_name or the unit of an open file to
  ! read, but not both.
  subroutine read_config_from_namelist(this, file_name, unit, is_success)

    use yomhook,      only : lhook, dr_hook, jphook
    use radiation_io, only : nulout, nulerr, nulrad, radiation_abort

    class(config_type), intent(inout)         :: this
    character(*),       intent(in),  optional :: file_name
    integer,            intent(in),  optional :: unit
    logical,            intent(out), optional :: is_success

    integer :: iosopen, iosread ! Status after calling open and read

    ! The following variables are read from the namelists and map
    ! directly onto members of the config_type derived type

    ! To be read from the radiation_config namelist 
    logical :: do_sw, do_lw, do_clear, do_sw_direct
    logical :: do_3d_effects, use_expm_everywhere, use_aerosols
    logical :: use_general_cloud_optics, use_general_aerosol_optics
    logical :: do_lw_side_emissivity
    logical :: do_3d_lw_multilayer_effects, do_fu_lw_ice_optics_bug
    logical :: do_lw_aerosol_scattering, do_lw_cloud_scattering
    logical :: do_save_radiative_properties, do_save_spectral_flux
    logical :: do_save_gpoint_flux, do_surface_sw_spectral_flux, do_toa_spectral_flux
    logical :: use_beta_overlap, do_lw_derivatives, use_vectorizable_generator
    logical :: do_sw_delta_scaling_with_gases
    logical :: do_canopy_fluxes_sw, do_canopy_fluxes_lw
    logical :: use_canopy_full_spectrum_sw, use_canopy_full_spectrum_lw
    logical :: do_canopy_gases_sw, do_canopy_gases_lw
    logical :: do_cloud_aerosol_per_sw_g_point, do_cloud_aerosol_per_lw_g_point
    logical :: do_weighted_surface_mapping
    logical :: use_spectral_solar_scaling, use_spectral_solar_cycle, use_updated_solar_spectrum
    integer :: n_regions, iverbose, iverbosesetup, n_aerosol_types
    real(jprb):: mono_lw_wavelength, mono_lw_total_od, mono_sw_total_od
    real(jprb):: mono_lw_single_scattering_albedo, mono_sw_single_scattering_albedo
    real(jprb):: mono_lw_asymmetry_factor, mono_sw_asymmetry_factor
    real(jprb):: cloud_inhom_decorr_scaling, cloud_fraction_threshold
    real(jprb):: clear_to_thick_fraction, max_gas_od_3d, max_cloud_od
    real(jprb):: cloud_mixing_ratio_threshold, overhead_sun_factor
    real(jprb):: max_3d_transfer_rate, min_cloud_effective_size
    real(jprb):: overhang_factor, encroachment_scaling
    character(511) :: directory_name, aerosol_optics_override_file_name
    character(511) :: liq_optics_override_file_name, ice_optics_override_file_name
    character(511) :: cloud_pdf_override_file_name
    character(511) :: gas_optics_sw_override_file_name, gas_optics_lw_override_file_name
    character(511) :: ssi_override_file_name
    character(63)  :: liquid_model_name, ice_model_name, gas_model_name
    character(63)  :: sw_solver_name, lw_solver_name, overlap_scheme_name
    character(63)  :: sw_entrapment_name, sw_encroachment_name, cloud_pdf_shape_name
    character(len=511) :: cloud_type_name(NMaxCloudTypes) = ["","","","","","","","","","","",""]
    logical :: use_thick_cloud_spectral_averaging(NMaxCloudTypes) &
         &  = [.false.,.false.,.false.,.false.,.false.,.false., &
         &     .false.,.false.,.false.,.false.,.false.,.false.]
    integer :: i_aerosol_type_map(NMaxAerosolTypes) ! More than 256 is an error

    logical :: do_nearest_spectral_sw_albedo
    logical :: do_nearest_spectral_lw_emiss
    real(jprb) :: sw_albedo_wavelength_bound(NMaxAlbedoIntervals-1)
    real(jprb) :: lw_emiss_wavelength_bound( NMaxAlbedoIntervals-1)
    integer :: i_sw_albedo_index(NMaxAlbedoIntervals)
    integer :: i_lw_emiss_index (NMaxAlbedoIntervals)

    integer :: iunit ! Unit number of namelist file

    namelist /radiation/ do_sw, do_lw, do_sw_direct, &
         &  do_3d_effects, do_lw_side_emissivity, do_clear, &
         &  do_save_radiative_properties, sw_entrapment_name, sw_encroachment_name, &
         &  do_3d_lw_multilayer_effects, do_fu_lw_ice_optics_bug, &
         &  do_save_spectral_flux, do_save_gpoint_flux, &
         &  do_surface_sw_spectral_flux, do_lw_derivatives, do_toa_spectral_flux, &
         &  do_lw_aerosol_scattering, do_lw_cloud_scattering, &
         &  n_regions, directory_name, gas_model_name, &
         &  ice_optics_override_file_name, liq_optics_override_file_name, &
         &  aerosol_optics_override_file_name, cloud_pdf_override_file_name, &
         &  gas_optics_sw_override_file_name, gas_optics_lw_override_file_name, &
         &  ssi_override_file_name, &
         &  liquid_model_name, ice_model_name, max_3d_transfer_rate, &
         &  min_cloud_effective_size, overhang_factor, encroachment_scaling, &
         &  use_canopy_full_spectrum_sw, use_canopy_full_spectrum_lw, &
         &  do_canopy_fluxes_sw, do_canopy_fluxes_lw, &
         &  do_canopy_gases_sw, do_canopy_gases_lw, &
         &  use_general_cloud_optics, use_general_aerosol_optics, &
         &  do_sw_delta_scaling_with_gases, overlap_scheme_name, &
         &  sw_solver_name, lw_solver_name, use_beta_overlap, use_vectorizable_generator, &
         &  use_expm_everywhere, iverbose, iverbosesetup, &
         &  cloud_inhom_decorr_scaling, cloud_fraction_threshold, &
         &  clear_to_thick_fraction, max_gas_od_3d, max_cloud_od, &
         &  cloud_mixing_ratio_threshold, overhead_sun_factor, &
         &  n_aerosol_types, i_aerosol_type_map, use_aerosols, &
         &  mono_lw_wavelength, mono_lw_total_od, mono_sw_total_od, &
         &  mono_lw_single_scattering_albedo, mono_sw_single_scattering_albedo, &
         &  mono_lw_asymmetry_factor, mono_sw_asymmetry_factor, &
         &  cloud_pdf_shape_name, cloud_type_name, use_thick_cloud_spectral_averaging, &
         &  do_nearest_spectral_sw_albedo, do_nearest_spectral_lw_emiss, &
         &  sw_albedo_wavelength_bound, lw_emiss_wavelength_bound, &
         &  i_sw_albedo_index, i_lw_emiss_index, &
         &  do_cloud_aerosol_per_lw_g_point, &
         &  do_cloud_aerosol_per_sw_g_point, do_weighted_surface_mapping, &
         &  use_spectral_solar_scaling, use_spectral_solar_cycle, use_updated_solar_spectrum
         
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_config:read',0,hook_handle)

    ! Copy default values from the original structure 
    do_sw = this%do_sw
    do_lw = this%do_lw
    do_sw_direct = this%do_sw_direct
    do_3d_effects = this%do_3d_effects
    do_3d_lw_multilayer_effects = this%do_3d_lw_multilayer_effects
    do_lw_side_emissivity = this%do_lw_side_emissivity
    do_clear = this%do_clear
    do_lw_aerosol_scattering = this%do_lw_aerosol_scattering
    do_lw_cloud_scattering = this%do_lw_cloud_scattering
    do_sw_delta_scaling_with_gases = this%do_sw_delta_scaling_with_gases
    do_fu_lw_ice_optics_bug = this%do_fu_lw_ice_optics_bug
    do_canopy_fluxes_sw = this%do_canopy_fluxes_sw
    do_canopy_fluxes_lw = this%do_canopy_fluxes_lw
    use_canopy_full_spectrum_sw = this%use_canopy_full_spectrum_sw
    use_canopy_full_spectrum_lw = this%use_canopy_full_spectrum_lw
    do_canopy_gases_sw = this%do_canopy_gases_sw
    do_canopy_gases_lw = this%do_canopy_gases_lw
    n_regions = this%nregions
    directory_name = this%directory_name
    cloud_pdf_override_file_name = this%cloud_pdf_override_file_name
    liq_optics_override_file_name = this%liq_optics_override_file_name
    ice_optics_override_file_name = this%ice_optics_override_file_name
    aerosol_optics_override_file_name = this%aerosol_optics_override_file_name
    gas_optics_sw_override_file_name = this%gas_optics_sw_override_file_name
    gas_optics_lw_override_file_name = this%gas_optics_lw_override_file_name
    ssi_override_file_name = this%ssi_override_file_name
    use_expm_everywhere = this%use_expm_everywhere
    use_aerosols = this%use_aerosols
    do_save_radiative_properties = this%do_save_radiative_properties
    do_save_spectral_flux = this%do_save_spectral_flux
    do_save_gpoint_flux = this%do_save_gpoint_flux
    do_lw_derivatives = this%do_lw_derivatives
    do_surface_sw_spectral_flux = this%do_surface_sw_spectral_flux
    do_toa_spectral_flux = this%do_toa_spectral_flux
    iverbose = this%iverbose
    iverbosesetup = this%iverbosesetup
    use_general_cloud_optics = this%use_general_cloud_optics
    use_general_aerosol_optics = this%use_general_aerosol_optics
    cloud_fraction_threshold = this%cloud_fraction_threshold
    cloud_mixing_ratio_threshold = this%cloud_mixing_ratio_threshold
    use_beta_overlap = this%use_beta_overlap
    use_vectorizable_generator = this%use_vectorizable_generator
    cloud_inhom_decorr_scaling = this%cloud_inhom_decorr_scaling
    clear_to_thick_fraction = this%clear_to_thick_fraction
    overhead_sun_factor = this%overhead_sun_factor
    max_gas_od_3d = this%max_gas_od_3d
    max_cloud_od = this%max_cloud_od
    max_3d_transfer_rate = this%max_3d_transfer_rate
    min_cloud_effective_size = this%min_cloud_effective_size
    cloud_type_name = this%cloud_type_name
    use_thick_cloud_spectral_averaging = this%use_thick_cloud_spectral_averaging

    overhang_factor = this%overhang_factor
    encroachment_scaling = -1.0_jprb
    gas_model_name = '' !DefaultGasModelName
    liquid_model_name = '' !DefaultLiquidModelName
    ice_model_name = '' !DefaultIceModelName
    sw_solver_name = '' !DefaultSwSolverName
    lw_solver_name = '' !DefaultLwSolverName
    sw_entrapment_name = ''
    sw_encroachment_name = ''
    overlap_scheme_name = ''
    cloud_pdf_shape_name = ''
    n_aerosol_types = this%n_aerosol_types
    mono_lw_wavelength = this%mono_lw_wavelength
    mono_lw_total_od = this%mono_lw_total_od
    mono_sw_total_od = this%mono_sw_total_od
    mono_lw_single_scattering_albedo = this%mono_lw_single_scattering_albedo
    mono_sw_single_scattering_albedo = this%mono_sw_single_scattering_albedo
    mono_lw_asymmetry_factor = this%mono_lw_asymmetry_factor
    mono_sw_asymmetry_factor = this%mono_sw_asymmetry_factor
    i_aerosol_type_map = this%i_aerosol_type_map
    do_nearest_spectral_sw_albedo = this%do_nearest_spectral_sw_albedo
    do_nearest_spectral_lw_emiss  = this%do_nearest_spectral_lw_emiss
    sw_albedo_wavelength_bound    = this%sw_albedo_wavelength_bound
    lw_emiss_wavelength_bound     = this%lw_emiss_wavelength_bound
    i_sw_albedo_index             = this%i_sw_albedo_index
    i_lw_emiss_index              = this%i_lw_emiss_index
    do_cloud_aerosol_per_lw_g_point = this%do_cloud_aerosol_per_lw_g_point
    do_cloud_aerosol_per_sw_g_point = this%do_cloud_aerosol_per_sw_g_point
    do_weighted_surface_mapping   = this%do_weighted_surface_mapping
    use_spectral_solar_scaling    = this%use_spectral_solar_scaling
    use_spectral_solar_cycle      = this%use_spectral_solar_cycle
    use_updated_solar_spectrum    = this%use_updated_solar_spectrum

    if (present(file_name) .and. present(unit)) then
      write(nulerr,'(a)') '*** Error: cannot specify both file_name and unit in call to config_type%read'
      call radiation_abort('Radiation configuration error')
    else if (.not. present(file_name) .and. .not. present(unit)) then
      write(nulerr,'(a)') '*** Error: neither file_name nor unit specified in call to config_type%read'
      call radiation_abort('Radiation configuration error')
    end if

    if (present(file_name)) then
      ! Open the namelist file
      iunit = nulrad
      open(unit=iunit, iostat=iosopen, file=trim(file_name))
    else
      ! Assume that iunit represents and open file
      iosopen = 0
      iunit = unit
    end if

    if (iosopen /= 0) then
      ! An error occurred opening the file
      if (present(is_success)) then
        is_success = .false.
        ! We now continue the subroutine so that the default values
        ! are placed in the config structure
      else
        write(nulerr,'(a,a,a)') '*** Error: namelist file "', &
             &                trim(file_name), '" not found'
        call radiation_abort('Radiation configuration error')
      end if
    else

      ! This version exits correctly, but provides less information
      ! about how the namelist was incorrect
      read(unit=iunit, iostat=iosread, nml=radiation)

      ! Depending on compiler this version provides more information
      ! about the error in the namelist
      !read(unit=iunit, nml=radiation)

      if (iosread /= 0) then
        ! An error occurred reading the file
        if (present(is_success)) then
          is_success = .false.
          ! We now continue the subroutine so that the default values
          ! are placed in the config structure
        else if (present(file_name)) then
          write(nulerr,'(a,a,a)') '*** Error reading namelist "radiation" from file "', &
               &      trim(file_name), '"'
          close(unit=iunit)
          call radiation_abort('Radiation configuration error')
        else
          write(nulerr,'(a,i0)') '*** Error reading namelist "radiation" from unit ', &
               &      iunit
          call radiation_abort('Radiation configuration error')
        end if
      end if

      if (present(file_name)) then
        close(unit=iunit)
      end if
    end if

    ! Copy namelist data into configuration object

    ! Start with verbosity levels, which should be within limits
    if (iverbose < 0) then
      iverbose = 0
    end if
    this%iverbose = iverbose

    if (iverbosesetup < 0) then
      iverbosesetup = 0
    end if
    this%iverbosesetup = iverbosesetup

    this%do_lw = do_lw
    this%do_sw = do_sw
    this%do_clear = do_clear
    this%do_sw_direct = do_sw_direct
    this%do_3d_effects = do_3d_effects
    this%do_3d_lw_multilayer_effects = do_3d_lw_multilayer_effects
    this%do_lw_side_emissivity = do_lw_side_emissivity
    this%use_expm_everywhere = use_expm_everywhere
    this%use_aerosols = use_aerosols
    this%do_lw_cloud_scattering = do_lw_cloud_scattering
    this%do_lw_aerosol_scattering = do_lw_aerosol_scattering
    this%nregions = n_regions
    this%do_surface_sw_spectral_flux = do_surface_sw_spectral_flux
    this%do_toa_spectral_flux = do_toa_spectral_flux
    this%do_sw_delta_scaling_with_gases = do_sw_delta_scaling_with_gases
    this%do_fu_lw_ice_optics_bug = do_fu_lw_ice_optics_bug
    this%do_canopy_fluxes_sw = do_canopy_fluxes_sw
    this%do_canopy_fluxes_lw = do_canopy_fluxes_lw
    this%use_canopy_full_spectrum_sw = use_canopy_full_spectrum_sw
    this%use_canopy_full_spectrum_lw = use_canopy_full_spectrum_lw
    this%do_canopy_gases_sw = do_canopy_gases_sw
    this%do_canopy_gases_lw = do_canopy_gases_lw
    this%mono_lw_wavelength = mono_lw_wavelength
    this%mono_lw_total_od = mono_lw_total_od
    this%mono_sw_total_od = mono_sw_total_od
    this%mono_lw_single_scattering_albedo = mono_lw_single_scattering_albedo
    this%mono_sw_single_scattering_albedo = mono_sw_single_scattering_albedo
    this%mono_lw_asymmetry_factor = mono_lw_asymmetry_factor
    this%mono_sw_asymmetry_factor = mono_sw_asymmetry_factor
    this%use_beta_overlap = use_beta_overlap
    this%use_vectorizable_generator = use_vectorizable_generator
    this%cloud_inhom_decorr_scaling = cloud_inhom_decorr_scaling
    this%clear_to_thick_fraction = clear_to_thick_fraction
    this%overhead_sun_factor = overhead_sun_factor
    this%max_gas_od_3d = max_gas_od_3d
    this%max_cloud_od = max_cloud_od
    this%max_3d_transfer_rate = max_3d_transfer_rate
    this%min_cloud_effective_size = max(1.0e-6_jprb, min_cloud_effective_size)
    this%cloud_type_name = cloud_type_name
    this%use_thick_cloud_spectral_averaging = use_thick_cloud_spectral_averaging
    if (encroachment_scaling >= 0.0_jprb) then
      this%overhang_factor = encroachment_scaling
      if (iverbose >= 1) then
        write(nulout, '(a)') 'Warning: radiation namelist parameter "encroachment_scaling" is deprecated: use "overhang_factor"'
      end if
    else
      this%overhang_factor = overhang_factor
    end if
    this%directory_name = directory_name
    this%cloud_pdf_override_file_name = cloud_pdf_override_file_name
    this%liq_optics_override_file_name = liq_optics_override_file_name
    this%ice_optics_override_file_name = ice_optics_override_file_name
    this%aerosol_optics_override_file_name = aerosol_optics_override_file_name
    this%gas_optics_sw_override_file_name = gas_optics_sw_override_file_name
    this%gas_optics_lw_override_file_name = gas_optics_lw_override_file_name
    this%ssi_override_file_name = ssi_override_file_name
    this%use_general_cloud_optics      = use_general_cloud_optics
    this%use_general_aerosol_optics    = use_general_aerosol_optics
    this%cloud_fraction_threshold = cloud_fraction_threshold
    this%cloud_mixing_ratio_threshold = cloud_mixing_ratio_threshold
    this%n_aerosol_types = n_aerosol_types
    this%do_save_radiative_properties = do_save_radiative_properties
    this%do_lw_derivatives = do_lw_derivatives
    this%do_save_spectral_flux = do_save_spectral_flux
    this%do_save_gpoint_flux = do_save_gpoint_flux
    this%do_nearest_spectral_sw_albedo = do_nearest_spectral_sw_albedo
    this%do_nearest_spectral_lw_emiss  = do_nearest_spectral_lw_emiss
    this%sw_albedo_wavelength_bound    = sw_albedo_wavelength_bound
    this%lw_emiss_wavelength_bound     = lw_emiss_wavelength_bound
    this%i_sw_albedo_index             = i_sw_albedo_index
    this%i_lw_emiss_index              = i_lw_emiss_index
    this%do_cloud_aerosol_per_lw_g_point = do_cloud_aerosol_per_lw_g_point
    this%do_cloud_aerosol_per_sw_g_point = do_cloud_aerosol_per_sw_g_point
    this%do_weighted_surface_mapping   = do_weighted_surface_mapping
    this%use_spectral_solar_scaling    = use_spectral_solar_scaling
    this%use_spectral_solar_cycle      = use_spectral_solar_cycle
    this%use_updated_solar_spectrum    = use_updated_solar_spectrum

    if (do_save_gpoint_flux) then
      ! Saving the fluxes every g-point overrides saving as averaged
      ! in a band, but this%do_save_spectral_flux needs to be TRUE as
      ! it is tested inside the solver routines to decide whether to
      ! save anything
      this%do_save_spectral_flux = .true.
    end if

    ! Determine liquid optics model
    call get_enum_code(liquid_model_name, LiquidModelName, &
         &            'liquid_model_name', this%i_liq_model)

    ! Determine ice optics model
    call get_enum_code(ice_model_name, IceModelName, &
         &            'ice_model_name', this%i_ice_model)

    ! Determine gas optics model
    call get_enum_code(gas_model_name, GasModelName, &
         &            'gas_model_name', this%i_gas_model)

    ! Determine solvers
    call get_enum_code(sw_solver_name, SolverName, &
         &            'sw_solver_name', this%i_solver_sw)
    call get_enum_code(lw_solver_name, SolverName, &
         &            'lw_solver_name', this%i_solver_lw)

    if (len_trim(sw_encroachment_name) > 1) then
      call get_enum_code(sw_encroachment_name, EncroachmentName, &
           &             'sw_encroachment_name', this%i_3d_sw_entrapment)
      write(nulout, '(a)') 'Warning: radiation namelist string "sw_encroachment_name" is deprecated: use "sw_entrapment_name"'
    else
      call get_enum_code(sw_entrapment_name, EntrapmentName, &
           &             'sw_entrapment_name', this%i_3d_sw_entrapment)
    end if

    ! Determine overlap scheme
    call get_enum_code(overlap_scheme_name, OverlapName, &
         &             'overlap_scheme_name', this%i_overlap_scheme)
    
    ! Determine cloud PDF shape 
    call get_enum_code(cloud_pdf_shape_name, PdfShapeName, &
         &             'cloud_pdf_shape_name', this%i_cloud_pdf_shape)

    this%i_aerosol_type_map = 0
    if (this%use_aerosols) then
      this%i_aerosol_type_map(1:n_aerosol_types) &
           &  = i_aerosol_type_map(1:n_aerosol_types)
    end if

    ! Will clouds be used at all?
    if ((this%do_sw .and. this%i_solver_sw /= ISolverCloudless) &
         &  .or. (this%do_lw .and. this%i_solver_lw /= ISolverCloudless)) then
      this%do_clouds = .true.
    else
      this%do_clouds = .false.
    end if

    if (this%i_gas_model == IGasModelIFSRRTMG &
         & .and. (this%use_general_cloud_optics &
         &        .or. this%use_general_aerosol_optics)) then
      if (this%do_sw .and. this%do_cloud_aerosol_per_sw_g_point) then
        write(nulout,'(a)') 'Warning: RRTMG SW only supports cloud/aerosol/surface optical properties per band, not per g-point'
        this%do_cloud_aerosol_per_sw_g_point = .false.
      end if
      if (this%do_lw .and. this%do_cloud_aerosol_per_lw_g_point) then
        write(nulout,'(a)') 'Warning: RRTMG LW only supports cloud/aerosol/surface optical properties per band, not per g-point'
        this%do_cloud_aerosol_per_lw_g_point = .false.
      end if
    end if


    ! Normal subroutine exit
    if (present(is_success)) then
      is_success = .true.
    end if

    if (lhook) call dr_hook('radiation_config:read',1,hook_handle)

  end subroutine read_config_from_namelist


  !---------------------------------------------------------------------
  ! This routine is called by radiation_interface:setup_radiation and
  ! it converts the user specified options into some more specific
  ! data such as data file names
  subroutine consolidate_config(this)

    use parkind1,     only : jprd
    use yomhook,      only : lhook, dr_hook, jphook
    use radiation_io, only : nulout, nulerr, radiation_abort

    class(config_type), intent(inout)         :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_config:consolidate',0,hook_handle)

    ! Check consistency of models
    if (this%do_canopy_fluxes_sw .and. .not. this%do_surface_sw_spectral_flux) then
      if (this%iverbosesetup >= 1) then
        write(nulout,'(a)') 'Warning: turning on do_surface_sw_spectral_flux as required by do_canopy_fluxes_sw'
      end if
      this%do_surface_sw_spectral_flux = .true.
    end if

    ! Will clouds be used at all?
    if ((this%do_sw .and. this%i_solver_sw /= ISolverCloudless) &
         &  .or. (this%do_lw .and. this%i_solver_lw /= ISolverCloudless)) then
      this%do_clouds = .true.
    else
      this%do_clouds = .false.
    end if

    ! SPARTACUS only works with Exp-Ran overlap scheme
    if ((       this%i_solver_sw == ISolverSPARTACUS &
         & .or. this%i_solver_lw == ISolverSPARTACUS &
         & .or. this%i_solver_sw == ISolverTripleclouds &
         & .or. this%i_solver_lw == ISolverTripleclouds) &
         & .and. this%i_overlap_scheme /= IOverlapExponentialRandom) then
      write(nulerr,'(a)') '*** Error: SPARTACUS/Tripleclouds solvers can only do Exponential-Random overlap'
      call radiation_abort('Radiation configuration error')
    end if

    if (jprb < jprd .and. this%iverbosesetup >= 1 &
         &  .and. (this%i_solver_sw == ISolverSPARTACUS &
         &    .or. this%i_solver_lw == ISolverSPARTACUS)) then
      write(nulout,'(a)') 'Warning: the SPARTACUS solver may be unstable in single precision'
    end if

    ! If ecCKD gas optics model is being used set relevant file names
    if (this%i_gas_model == IGasModelECCKD) then

      ! This gas optics model requires the general cloud and
      ! aerosol optics settings
      if (.not. this%use_general_cloud_optics) then
        write(nulerr,'(a)') '*** Error: ecCKD gas optics model requires general cloud optics'
        call radiation_abort('Radiation configuration error')
      end if
      if (.not. this%use_general_aerosol_optics) then
        write(nulerr,'(a)') '*** Error: ecCKD gas optics model requires general aerosol optics'
        call radiation_abort('Radiation configuration error')
      end if

      if (len_trim(this%gas_optics_sw_override_file_name) > 0) then
        if (this%gas_optics_sw_override_file_name(1:1) == '/') then
          this%gas_optics_sw_file_name = trim(this%gas_optics_sw_override_file_name)
        else
          this%gas_optics_sw_file_name = trim(this%directory_name) &
               &  // '/' // trim(this%gas_optics_sw_override_file_name)
        end if
      else
        ! In the IFS, the gas optics files should be specified in
        ! ifs/module/radiation_setup.F90, not here
        this%gas_optics_sw_file_name = trim(this%directory_name) &
             &  // "/ecckd-1.0_sw_climate_rgb-32b_ckd-definition.nc"
      end if

      if (len_trim(this%gas_optics_lw_override_file_name) > 0) then
        if (this%gas_optics_lw_override_file_name(1:1) == '/') then
          this%gas_optics_lw_file_name = trim(this%gas_optics_lw_override_file_name)
        else
          this%gas_optics_lw_file_name = trim(this%directory_name) &
               &  // '/' // trim(this%gas_optics_lw_override_file_name)
        end if
      else
        ! In the IFS, the gas optics files should be specified in
        ! ifs/module/radiation_setup.F90, not here
        this%gas_optics_lw_file_name = trim(this%directory_name) &
             &  // "/ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc"
      end if

    end if

    if (this%use_spectral_solar_cycle) then
      if (this%i_gas_model /= IGasModelECCKD) then
        write(nulerr,'(a)') '*** Error: solar cycle only available with ecCKD gas optics model'
        call radiation_abort('Radiation configuration error')
      else
        ! Add directory name to solar spectral irradiance file, if
        ! provided and does not start with '/'
        if (len_trim(this%ssi_override_file_name) > 0) then
          if (this%ssi_override_file_name(1:1) /= '/') then
            this%ssi_file_name = trim(this%directory_name) &
                 &  // '/' // trim(this%ssi_override_file_name)
          else
            this%ssi_file_name = trim(this%ssi_override_file_name)
          end if
        else
          this%ssi_file_name = 'ssi_nrl2.nc'
        end if
      end if
    end if
    
    ! Set aerosol optics file name
    if (len_trim(this%aerosol_optics_override_file_name) > 0) then
      if (this%aerosol_optics_override_file_name(1:1) == '/') then
        this%aerosol_optics_file_name = trim(this%aerosol_optics_override_file_name)
      else
        this%aerosol_optics_file_name = trim(this%directory_name) &
             &  // '/' // trim(this%aerosol_optics_override_file_name)
      end if
    else
      ! In the IFS, the aerosol optics file should be specified in
      ! ifs/module/radiation_setup.F90, not here
      if (this%use_general_aerosol_optics) then
         this%aerosol_optics_file_name &
             &   = trim(this%directory_name) // "/aerosol_ifs_48R1.nc"       
      else
        this%aerosol_optics_file_name &
             &   = trim(this%directory_name) // "/aerosol_ifs_rrtm_46R1_with_NI_AM.nc"
      end if
    end if

    ! Set liquid optics file name
    if (len_trim(this%liq_optics_override_file_name) > 0) then
      if (this%liq_optics_override_file_name(1:1) == '/') then
        this%liq_optics_file_name = trim(this%liq_optics_override_file_name)
      else
        this%liq_optics_file_name = trim(this%directory_name) &
             &  // '/' // trim(this%liq_optics_override_file_name)
      end if
    else if (this%i_liq_model == ILiquidModelSOCRATES) then
      this%liq_optics_file_name &
           &  = trim(this%directory_name) // "/socrates_droplet_scattering_rrtm.nc"
    else if (this%i_liq_model == ILiquidModelSlingo) then
      this%liq_optics_file_name &
           &  = trim(this%directory_name) // "/slingo_droplet_scattering_rrtm.nc"
    end if

    ! Set ice optics file name
    if (len_trim(this%ice_optics_override_file_name) > 0) then
      if (this%ice_optics_override_file_name(1:1) == '/') then
        this%ice_optics_file_name = trim(this%ice_optics_override_file_name)
      else
        this%ice_optics_file_name = trim(this%directory_name) &
             &  // '/' // trim(this%ice_optics_override_file_name)
      end if
    else if (this%i_ice_model == IIceModelFu) then
      this%ice_optics_file_name &
           &   = trim(this%directory_name) // "/fu_ice_scattering_rrtm.nc"
    else if (this%i_ice_model == IIceModelBaran) then
      this%ice_optics_file_name &
           &   = trim(this%directory_name) // "/baran_ice_scattering_rrtm.nc"
    else if (this%i_ice_model == IIceModelBaran2016) then
      this%ice_optics_file_name &
           &   = trim(this%directory_name) // "/baran2016_ice_scattering_rrtm.nc"
    else if (this%i_ice_model == IIceModelBaran2017) then
      this%ice_optics_file_name &
           &   = trim(this%directory_name) // "/baran2017_ice_scattering_rrtm.nc"
    else if (this%i_ice_model == IIceModelYi) then
      this%ice_optics_file_name &
           &   = trim(this%directory_name) // "/yi_ice_scattering_rrtm.nc"
    end if

    ! Set cloud-water PDF look-up table file name
    if (len_trim(this%cloud_pdf_override_file_name) > 0) then
      if (this%cloud_pdf_override_file_name(1:1) == '/') then
        this%cloud_pdf_file_name = trim(this%cloud_pdf_override_file_name)
      else
        this%cloud_pdf_file_name = trim(this%directory_name) &
             &  // '/' // trim(this%cloud_pdf_override_file_name)
      end if
    elseif (this%i_cloud_pdf_shape == IPdfShapeLognormal) then
      this%cloud_pdf_file_name = trim(this%directory_name) // "/mcica_lognormal.nc"
    else
      this%cloud_pdf_file_name = trim(this%directory_name) // "/mcica_gamma.nc"
    end if

    ! Aerosol data
    if (this%n_aerosol_types < 0 &
         &  .or. this%n_aerosol_types > NMaxAerosolTypes) then
      write(nulerr,'(a,i0)') '*** Error: number of aerosol types must be between 0 and ', &
           &  NMaxAerosolTypes
      call radiation_abort('Radiation configuration error')
    end if

    if (this%use_aerosols .and. this%n_aerosol_types == 0) then
      if (this%iverbosesetup >= 2) then
        write(nulout, '(a)') 'Aerosols on but n_aerosol_types=0: optical properties to be computed outside ecRad'
      end if
    end if

    ! In the monochromatic case we need to override the liquid, ice
    ! and aerosol models to ensure compatibility
    if (this%i_gas_model == IGasModelMonochromatic) then
      this%i_liq_model = ILiquidModelMonochromatic
      this%i_ice_model = IIceModelMonochromatic
      this%use_aerosols = .false.
    end if

    ! McICA solver currently can't store full profiles of spectral fluxes
    if (this%i_solver_sw == ISolverMcICA) then
      this%do_save_spectral_flux = .false.
    end if

    if (this%i_solver_sw == ISolverSPARTACUS .and. this%do_sw_delta_scaling_with_gases) then
      write(nulerr,'(a)') '*** Error: SW delta-Eddington scaling with gases not possible with SPARTACUS solver'
      call radiation_abort('Radiation configuration error')
    end if

    if ((this%do_lw .and. this%do_sw) .and. &
         & (     (      this%i_solver_sw == ISolverHomogeneous  &
         &        .and. this%i_solver_lw /= ISolverHomogeneous) &
         &  .or. (      this%i_solver_sw /= ISolverHomogeneous  &
         &        .and. this%i_solver_lw == ISolverHomogeneous) &
         & ) ) then
      write(nulerr,'(a)') '*** Error: if one solver is "Homogeneous" then the other must be'
      call radiation_abort('Radiation configuration error')
    end if

    ! Set is_homogeneous if the active solvers are homogeneous, since
    ! this affects how "in-cloud" water contents are computed
    if (        (this%do_sw .and. this%i_solver_sw == ISolverHomogeneous) &
         & .or. (this%do_lw .and. this%i_solver_lw == ISolverHomogeneous)) then
      this%is_homogeneous = .true.
    end if

    this%is_consolidated = .true.

    if (lhook) call dr_hook('radiation_config:consolidate',1,hook_handle)

  end subroutine consolidate_config


  !---------------------------------------------------------------------
  ! This subroutine sets members of the configuration object via
  ! optional arguments, and any member not specified is left
  ! untouched. Therefore, this should be called after taking data from
  ! the namelist.
  subroutine set_config(config, directory_name, &
       &  do_lw, do_sw, &
       &  do_lw_aerosol_scattering, do_lw_cloud_scattering, &
       &  do_sw_direct)

    class(config_type), intent(inout):: config
    character(len=*), intent(in), optional  :: directory_name
    logical, intent(in), optional           :: do_lw, do_sw
    logical, intent(in), optional           :: do_lw_aerosol_scattering
    logical, intent(in), optional           :: do_lw_cloud_scattering
    logical, intent(in), optional           :: do_sw_direct

    if (present(do_lw)) then
       config%do_lw = do_lw
    end if

    if(present(do_sw)) then
       config%do_sw = do_sw
    end if

    if (present(do_sw_direct)) then
       config%do_sw_direct = do_sw_direct
    end if

    if (present(directory_name)) then
       config%directory_name = trim(directory_name)
    end if

    if (present(do_lw_aerosol_scattering)) then
       config%do_lw_aerosol_scattering = .true.
    end if

    if (present(do_lw_cloud_scattering)) then
       config%do_lw_cloud_scattering = .true.
    end if

  end subroutine set_config


  !---------------------------------------------------------------------
  ! Print configuration information to standard output
  subroutine print_config(this, iverbose)

    use radiation_io, only : nulout

    class(config_type), intent(in) :: this

    integer, optional,  intent(in) :: iverbose
    integer                        :: i_local_verbose

    if (present(iverbose)) then
      i_local_verbose = iverbose
    else
      i_local_verbose = this%iverbose
    end if

    if (i_local_verbose >= 2) then
      !---------------------------------------------------------------------
      write(nulout, '(a)') 'General settings:'
      write(nulout, '(a,a,a)') '  Data files expected in "', &
           &                   trim(this%directory_name), '"'
      call print_logical('  Clear-sky calculations are', 'do_clear', this%do_clear)
      call print_logical('  Saving intermediate radiative properties', &
           &   'do_save_radiative_properties', this%do_save_radiative_properties)
      call print_logical('  Saving spectral flux profiles', &
           &   'do_save_spectral_flux', this%do_save_spectral_flux)
      call print_enum('  Gas model is', GasModelName, 'i_gas_model', &
           &          this%i_gas_model)
      call print_logical('  Aerosols are', 'use_aerosols', this%use_aerosols)
      if (this%use_aerosols) then
        call print_logical('  General aerosol optics', &
             &             'use_general_aerosol_optics', this%use_general_aerosol_optics)
      end if
      if (this%do_clouds) then
        write(nulout,'(a)') '  Clouds are ON'
      else
        write(nulout,'(a)') '  Clouds are OFF'
      end if
      if (this%do_sw) then
        call print_logical('  Do cloud/aerosol/surface SW properties per g-point', &
             &  'do_cloud_aerosol_per_sw_g_point', this%do_cloud_aerosol_per_sw_g_point)
      end if
      if (this%do_lw) then
        call print_logical('  Do cloud/aerosol/surface LW properties per g-point', &
             &  'do_cloud_aerosol_per_lw_g_point', this%do_cloud_aerosol_per_lw_g_point)
      end if
      if (this%do_sw) then
        call print_logical('  Represent solar cycle in spectral irradiance', &
             &  'use_spectral_solar_cycle', this%use_spectral_solar_cycle)
        call print_logical('  Scale spectral solar irradiance', &
             &  'use_spectral_solar_scaling', this%use_spectral_solar_scaling)
      end if
      
      !---------------------------------------------------------------------
      write(nulout, '(a)') 'Surface and top-of-atmosphere settings:'
      call print_logical('  Saving top-of-atmosphere spectral fluxes', &
           &   'do_toa_spectral_flux', this%do_toa_spectral_flux)
      if (this%do_sw) then
        call print_logical('  Saving surface shortwave spectral fluxes', &
             &   'do_surface_sw_spectral_flux', this%do_surface_sw_spectral_flux)
        call print_logical('  Saving surface shortwave fluxes in abledo bands', &
             &   'do_canopy_fluxes_sw', this%do_canopy_fluxes_sw)
      end if
      if (this%do_lw) then
        call print_logical('  Saving surface longwave fluxes in emissivity bands', &
             &   'do_canopy_fluxes_lw', this%do_canopy_fluxes_lw)
        call print_logical('  Longwave derivative calculation is', &
             &   'do_lw_derivatives',this%do_lw_derivatives)
      end if
      if (this%do_sw) then
        call print_logical('  Nearest-neighbour spectral albedo mapping', &
             &   'do_nearest_spectral_sw_albedo', this%do_nearest_spectral_sw_albedo)
      end if
      if (this%do_lw) then
        call print_logical('  Nearest-neighbour spectral emissivity mapping', &
             &   'do_nearest_spectral_lw_emiss', this%do_nearest_spectral_lw_emiss)
      end if
      call print_logical('  Planck-weighted surface albedo/emiss mapping', &
           &   'do_weighted_surface_mapping', this%do_weighted_surface_mapping)

      !---------------------------------------------------------------------
      if (this%do_clouds) then
        write(nulout, '(a)') 'Cloud settings:'
        call print_real('  Cloud fraction threshold', &
             &   'cloud_fraction_threshold', this%cloud_fraction_threshold)
        call print_real('  Cloud mixing-ratio threshold', &
             &   'cloud_mixing_ratio_threshold', this%cloud_mixing_ratio_threshold)
        call print_logical('  General cloud optics', &
             &             'use_general_cloud_optics', this%use_general_cloud_optics)
        if (.not. this%use_general_cloud_optics) then
          call print_enum('  Liquid optics scheme is', LiquidModelName, &
               &          'i_liq_model',this%i_liq_model)
          call print_enum('  Ice optics scheme is', IceModelName, &
               &          'i_ice_model',this%i_ice_model)
          if (this%i_ice_model == IIceModelFu) then
            call print_logical('  Longwave ice optics bug in Fu scheme is', &
                 &   'do_fu_lw_ice_optics_bug',this%do_fu_lw_ice_optics_bug)
          end if
        end if
        call print_enum('  Cloud overlap scheme is', OverlapName, &
             &          'i_overlap_scheme',this%i_overlap_scheme)
        call print_logical('  Use "beta" overlap parameter is', &
             &   'use_beta_overlap', this%use_beta_overlap)
        call print_enum('  Cloud PDF shape is', PdfShapeName, &
             &          'i_cloud_pdf_shape',this%i_cloud_pdf_shape)
        call print_real('  Cloud inhom decorrelation scaling', &
             &   'cloud_inhom_decorr_scaling', this%cloud_inhom_decorr_scaling)
      end if

      !---------------------------------------------------------------------
      write(nulout, '(a)') 'Solver settings:'
      if (this%do_sw) then
        call print_enum('  Shortwave solver is', SolverName, &
             &          'i_solver_sw', this%i_solver_sw)
        
        if (this%i_gas_model == IGasModelMonochromatic) then
          call print_real('  Shortwave atmospheric optical depth', &
               &   'mono_sw_total_od', this%mono_sw_total_od)
          call print_real('  Shortwave particulate single-scattering albedo', &
               &   'mono_sw_single_scattering_albedo', &
               &   this%mono_sw_single_scattering_albedo)
          call print_real('  Shortwave particulate asymmetry factor', &
               &   'mono_sw_asymmetry_factor', &
               &   this%mono_sw_asymmetry_factor)
        end if
        call print_logical('  Shortwave delta scaling after merge with gases', &
             &   'do_sw_delta_scaling_with_gases', &
             &   this%do_sw_delta_scaling_with_gases)
      else
        call print_logical('  Shortwave calculations are','do_sw',this%do_sw)
      end if

      if (this%do_lw) then
        call print_enum('  Longwave solver is', SolverName, 'i_solver_lw', &
             &          this%i_solver_lw)

        if (this%i_gas_model == IGasModelMonochromatic) then
          if (this%mono_lw_wavelength > 0.0_jprb) then
            call print_real('  Longwave effective wavelength (m)', &
                 &   'mono_lw_wavelength', this%mono_lw_wavelength)
          else
            write(nulout,'(a)') '  Longwave fluxes are broadband                              (mono_lw_wavelength<=0)'
          end if
          call print_real('  Longwave atmospheric optical depth', &
               &   'mono_lw_total_od', this%mono_lw_total_od)  
          call print_real('  Longwave particulate single-scattering albedo', &
               &   'mono_lw_single_scattering_albedo', &
               &   this%mono_lw_single_scattering_albedo)
          call print_real('  Longwave particulate asymmetry factor', &
               &   'mono_lw_asymmetry_factor', &
               &   this%mono_lw_asymmetry_factor)
        end if
        call print_logical('  Longwave cloud scattering is', &
             &   'do_lw_cloud_scattering',this%do_lw_cloud_scattering)
        call print_logical('  Longwave aerosol scattering is', &
             &   'do_lw_aerosol_scattering',this%do_lw_aerosol_scattering)
      else
        call print_logical('  Longwave calculations are','do_lw', this%do_lw)
      end if

      if (this%i_solver_sw == ISolverSpartacus &
           &  .or. this%i_solver_lw == ISolverSpartacus) then
        write(nulout, '(a)') '  SPARTACUS options:'
        call print_integer('    Number of regions', 'n_regions', this%nregions)
        call print_real('    Max cloud optical depth per layer', &
             &   'max_cloud_od', this%max_cloud_od)
        call print_enum('    Shortwave entrapment is', EntrapmentName, &
             &          'i_3d_sw_entrapment', this%i_3d_sw_entrapment)
        call print_logical('    Multilayer longwave horizontal transport is', &
             'do_3d_lw_multilayer_effects', this%do_3d_lw_multilayer_effects)
        call print_logical('    Use matrix exponential everywhere is', &
             &   'use_expm_everywhere', this%use_expm_everywhere)
        call print_logical('    3D effects are', 'do_3d_effects', &
             &             this%do_3d_effects)

        if (this%do_3d_effects) then
          call print_logical('    Longwave side emissivity parameterization is', &
               &  'do_lw_side_emissivity', this%do_lw_side_emissivity)
          call print_real('    Clear-to-thick edge fraction is', &
               &  'clear_to_thick_fraction', this%clear_to_thick_fraction)
          call print_real('    Overhead sun factor is', &
               &  'overhead_sun_factor', this%overhead_sun_factor)
          call print_real('    Max gas optical depth for 3D effects', &
               &   'max_gas_od_3d', this%max_gas_od_3d)
          call print_real('    Max 3D transfer rate', &
               &   'max_3d_transfer_rate', this%max_3d_transfer_rate)
          call print_real('    Min cloud effective size (m)', &
               &   'min_cloud_effective_size', this%min_cloud_effective_size)
          call print_real('    Overhang factor', &
               &   'overhang_factor', this%overhang_factor)
        end if

      else if (this%i_solver_sw == ISolverMcICA &
           &  .or. this%i_solver_lw == ISolverMcICA) then
        call print_logical('  Use vectorizable McICA cloud generator', &
             &   'use_vectorizable_generator', this%use_vectorizable_generator)
      end if
            
    end if
    
  end subroutine print_config


  !---------------------------------------------------------------------
  ! In order to estimate UV and photosynthetically active radiation,
  ! we need weighted sum of fluxes considering wavelength range
  ! required.  This routine returns information for how to correctly
  ! weight output spectral fluxes for a range of input wavelengths.
  ! Note that this is approximate; internally it may be assumed that
  ! the energy is uniformly distributed in wavenumber space, for
  ! example.  If the character string "weighting_name" is present, and
  ! iverbose>=2, then information on the weighting will be provided on
  ! nulout.
  subroutine get_sw_weights(this, wavelength1, wavelength2, &
       &                    nweights, iband, weight, weighting_name)

    use parkind1, only : jprb
    use radiation_io, only : nulout, nulerr, radiation_abort
    use radiation_spectral_definition, only : SolarReferenceTemperature

    class(config_type), intent(in) :: this
    ! Range of wavelengths to get weights for (m)
    real(jprb), intent(in) :: wavelength1, wavelength2
    ! Output number of weights needed
    integer,    intent(out)   :: nweights
    ! Only write to the first nweights of these arrays: they contain
    ! the indices to the non-zero bands, and the weight in each of
    ! those bands
    integer,    intent(out)   :: iband(:)
    real(jprb), intent(out)   :: weight(:)
    character(len=*), optional, intent(in) :: weighting_name

    real(jprb), allocatable   :: mapping(:,:)

    ! Internally we deal with wavenumber
    real(jprb) :: wavenumber1, wavenumber2 ! cm-1

    real(jprb) :: wavenumber1_band, wavenumber2_band ! cm-1

    integer :: jband ! Loop index for spectral band

    if (this%n_bands_sw <= 0) then
      write(nulerr,'(a)') '*** Error: get_sw_weights called before number of shortwave bands set'
      call radiation_abort('Radiation configuration error')
    end if

    ! Convert wavelength range (m) to wavenumber (cm-1)
    wavenumber1 = 0.01_jprb / wavelength2
    wavenumber2 = 0.01_jprb / wavelength1

    call this%gas_optics_sw%spectral_def%calc_mapping_from_bands( &
         &  [wavelength1, wavelength2], [1, 2, 3], mapping, &
         &  use_bands=(.not. this%do_cloud_aerosol_per_sw_g_point), use_fluxes=.true.)

    ! "mapping" now contains a 3*nband matrix, where mapping(2,:)
    ! contains the weights of interest.  We now find the non-zero weights
    nweights = 0
    do jband = 1,size(mapping,2)
      if (mapping(2,jband) > 0.0_jprb) then
        nweights = nweights+1
        iband(nweights) = jband;
        weight(nweights) = mapping(2,jband)
      end if
    end do

    if (nweights == 0) then
      write(nulerr,'(a,e8.4,a,e8.4,a)') '*** Error: wavelength range ', &
           &  wavelength1, ' to ', wavelength2, ' m is outside shortwave band'
      call radiation_abort('Radiation configuration error')
    else if (this%iverbosesetup >= 2 .and. present(weighting_name)) then
      write(nulout,'(a,a,a,f6.0,a,f6.0,a)') 'Spectral weights for ', &
           &  weighting_name, ' (', wavenumber1, ' to ', &
           &  wavenumber2, ' cm-1):'
      if (this%do_cloud_aerosol_per_sw_g_point) then
        do jband = 1, nweights
          write(nulout, '(a,i0,a,f8.4)') '  Shortwave g point ', iband(jband), ': ', weight(jband)
        end do
      else
        do jband = 1, nweights
          wavenumber1_band = this%gas_optics_sw%spectral_def%wavenumber1_band(iband(jband))
          wavenumber2_band = this%gas_optics_sw%spectral_def%wavenumber2_band(iband(jband))
          write(nulout, '(a,i0,a,f6.0,a,f6.0,a,f8.4)') '  Shortwave band ', &
               &  iband(jband), ' (', wavenumber1_band, ' to ', &
               &  wavenumber2_band, ' cm-1): ', weight(jband)
        end do
      end if
    end if

  end subroutine get_sw_weights

  
  !---------------------------------------------------------------------
  ! As get_sw_weights but suitable for a larger number of spectral
  ! diagnostics at once: a set of monotonically increasing wavelength
  ! bounds are provided (m), and a mapping matrix is allocated and
  ! returned such that y=matmul(mapping,x), where x is a set of
  ! band-wise fluxes after calling ecRad, e.g. flux%sw_dn_surf_band,
  ! and y is the resulting fluxes in each of the wavenumber
  ! intervals. If the character string "weighting_name" is present,
  ! and iverbose>=2, then information on the weighting will be
  ! provided on nulout.
  subroutine get_sw_mapping(this, wavelength_bound, mapping, weighting_name)

    use parkind1, only : jprb
    use radiation_io, only : nulout, nulerr, radiation_abort
    use radiation_spectral_definition, only : SolarReferenceTemperature

    class(config_type), intent(in) :: this
    ! Range of wavelengths to get weights for (m)
    real(jprb), intent(in)  :: wavelength_bound(:)
    real(jprb), intent(out), allocatable :: mapping(:,:)
    character(len=*), optional, intent(in) :: weighting_name

    real(jprb), allocatable :: mapping_local(:,:)
    integer,    allocatable :: diag_ind(:)

    integer :: ninterval
    
    integer :: jint  ! Loop for interval
    
    if (this%n_bands_sw <= 0) then
      write(nulerr,'(a)') '*** Error: get_sw_mapping called before number of shortwave bands set'
      call radiation_abort('Radiation configuration error')
    end if

    ninterval = size(wavelength_bound)-1
    allocate(diag_ind(ninterval+2))
    diag_ind = 0
    do jint = 1,ninterval+2
      diag_ind(jint) = jint
    end do
    
    call this%gas_optics_sw%spectral_def%calc_mapping_from_bands( &
         &  wavelength_bound, diag_ind, mapping_local, &
         &  use_bands=(.not. this%do_cloud_aerosol_per_sw_g_point), use_fluxes=.false.)

    ! "mapping" now contains a (ninterval+2)*nband matrix, where the
    ! first and last rows correspond to wavelengths smaller than the
    ! first and larger than the last, which we discard
    mapping = mapping_local(2:ninterval+1,:)

    if (this%iverbosesetup >= 2 .and. present(weighting_name)) then
      write(nulout,'(a,a)') 'Spectral mapping generated for ', &
           &  weighting_name
        if (this%do_cloud_aerosol_per_sw_g_point) then
          write(nulout,'(a,i0,a,i0,a,f9.3,a,f9.3,a)') '  from ', size(mapping,2), ' g-points to ', &
             &  size(mapping,1), ' wavelength intervals between ', &
             &  wavelength_bound(1)*1.0e6_jprb, ' um and ', wavelength_bound(ninterval+1)*1.0e6_jprb, ' um'
        else
          write(nulout,'(a,i0,a,i0,a,f9.3,a,f9.3,a)') '  from ', size(mapping,2), ' bands to ', &
               &  size(mapping,1), ' wavelength intervals between ', &
               &  wavelength_bound(1)*1.0e6_jprb, ' um and ', wavelength_bound(ninterval+1)*1.0e6_jprb, ' um'
      end if
    end if
    
  end subroutine get_sw_mapping


  !---------------------------------------------------------------------
  ! The input shortwave surface albedo coming in is likely to be in
  ! different spectral intervals to the gas model in the radiation
  ! scheme. We assume that the input albedo is defined within
  ! "ninterval" spectral intervals covering the wavelength range 0 to
  ! infinity, but allow for the possibility that two intervals may be
  ! indexed back to the same albedo band.  
  subroutine define_sw_albedo_intervals(this, ninterval, wavelength_bound, &
       &                                i_intervals, do_nearest)

    use radiation_io, only : nulerr, radiation_abort
    use radiation_spectral_definition, only : SolarReferenceTemperature

    class(config_type),   intent(inout) :: this
    ! Number of spectral intervals in which albedo is defined
    integer,              intent(in)    :: ninterval
    ! Monotonically increasing wavelength bounds between intervals,
    ! not including the outer bounds (which are assumed to be zero and
    ! infinity)
    real(jprb),           intent(in)    :: wavelength_bound(ninterval-1)
    ! The albedo indices corresponding to each interval
    integer,              intent(in)    :: i_intervals(ninterval)
    logical,    optional, intent(in)    :: do_nearest
    
    if (ninterval > NMaxAlbedoIntervals) then
      write(nulerr,'(a,i0,a,i0)') '*** Error: ', ninterval, &
           &  ' albedo intervals exceeds maximum of ', NMaxAlbedoIntervals
      call radiation_abort('Radiation configuration error')
    end if

    if (present(do_nearest)) then
      this%do_nearest_spectral_sw_albedo = do_nearest
    else
      this%do_nearest_spectral_sw_albedo = .false.
    end if
    if (ninterval > 1) then
      this%sw_albedo_wavelength_bound(1:ninterval-1) = wavelength_bound(1:ninterval-1)
    end if
    this%sw_albedo_wavelength_bound(ninterval:)    = -1.0_jprb
    this%i_sw_albedo_index(1:ninterval)            = i_intervals(1:ninterval)
    this%i_sw_albedo_index(ninterval+1:)           = 0

    ! If this routine is called before setup_radiation then the
    ! spectral intervals are not yet known
    ! consolidate_sw_albedo_intervals is called later.  Otherwise it
    ! is called immediately and overwrites any existing mapping.
    if (this%is_consolidated) then
      call this%consolidate_sw_albedo_intervals
    end if

  end subroutine define_sw_albedo_intervals


  !---------------------------------------------------------------------
  ! As define_sw_albedo_intervals but for longwave emissivity
  subroutine define_lw_emiss_intervals(this, ninterval, wavelength_bound, &
       &                                i_intervals, do_nearest)

    use radiation_io, only : nulerr, radiation_abort
    use radiation_spectral_definition, only : TerrestrialReferenceTemperature

    class(config_type),   intent(inout) :: this
    ! Number of spectral intervals in which emissivity is defined
    integer,              intent(in)    :: ninterval
    ! Monotonically increasing wavelength bounds between intervals,
    ! not including the outer bounds (which are assumed to be zero and
    ! infinity)
    real(jprb),           intent(in)    :: wavelength_bound(ninterval-1)
    ! The emissivity indices corresponding to each interval
    integer,              intent(in)    :: i_intervals(ninterval)
    logical,    optional, intent(in)    :: do_nearest
    
    if (ninterval > NMaxAlbedoIntervals) then
      write(nulerr,'(a,i0,a,i0)') '*** Error: ', ninterval, &
           &  ' emissivity intervals exceeds maximum of ', NMaxAlbedoIntervals
      call radiation_abort('Radiation configuration error')
    end if

    if (present(do_nearest)) then
      this%do_nearest_spectral_lw_emiss = do_nearest
    else
      this%do_nearest_spectral_lw_emiss = .false.
    end if
    if (ninterval > 1) then
      this%lw_emiss_wavelength_bound(1:ninterval-1) = wavelength_bound(1:ninterval-1)
    end if
    this%lw_emiss_wavelength_bound(ninterval:)    = -1.0_jprb
    this%i_lw_emiss_index(1:ninterval)            = i_intervals(1:ninterval)
    this%i_lw_emiss_index(ninterval+1:)           = 0

    if (this%is_consolidated) then
      call this%consolidate_lw_emiss_intervals
    end if

  end subroutine define_lw_emiss_intervals


  !---------------------------------------------------------------------
  ! Set the wavelengths (m) at which monochromatic aerosol properties
  ! are required. This routine must be called before consolidation of
  ! settings.
  subroutine set_aerosol_wavelength_mono(this, wavelength_mono)

    use radiation_io, only : nulerr, radiation_abort
    
    class(config_type), intent(inout) :: this
    real(jprb),         intent(in)    :: wavelength_mono(:)

    if (this%is_consolidated) then
      write(nulerr,'(a)') '*** Errror: set_aerosol_wavelength_mono must be called before setup_radiation'
      call radiation_abort('Radiation configuration error')
    end if
   
    if (allocated(this%aerosol_optics%wavelength_mono)) then
      deallocate(this%aerosol_optics%wavelength_mono)
    end if
    allocate(this%aerosol_optics%wavelength_mono(size(wavelength_mono)))
    this%aerosol_optics%wavelength_mono = wavelength_mono

  end subroutine set_aerosol_wavelength_mono


  !---------------------------------------------------------------------
  ! Consolidate the surface shortwave albedo intervals with the
  ! band/g-point intervals
  subroutine consolidate_sw_albedo_intervals(this)

    use radiation_io, only : nulout
    use radiation_spectral_definition, only : SolarReferenceTemperature

    class(config_type),   intent(inout) :: this

    integer :: ninterval, jint, jband

    ! Count the number of albedo/emissivity intervals
    ninterval = 0
    do jint = 1,NMaxAlbedoIntervals
      if (this%i_sw_albedo_index(jint) > 0) then
        ninterval = jint
      else
        exit
      end if
    end do

    if (ninterval < 1) then
      ! The user has not specified shortwave albedo bands - assume
      ! only one
      ninterval = 1
      this%i_sw_albedo_index(1) = 1
      this%i_sw_albedo_index(2:) = 0
      if (this%use_canopy_full_spectrum_sw) then
        this%n_canopy_bands_sw = this%n_g_sw
      else 
        this%n_canopy_bands_sw = 1
      end if
    else
      if (this%use_canopy_full_spectrum_sw) then
        this%n_canopy_bands_sw = this%n_g_sw
      else 
        this%n_canopy_bands_sw = maxval(this%i_sw_albedo_index(1:ninterval))
      end if
    end if
    
    if (this%do_weighted_surface_mapping) then
      call this%gas_optics_sw%spectral_def%calc_mapping_from_bands( &
           &  this%sw_albedo_wavelength_bound(1:ninterval-1), this%i_sw_albedo_index(1:ninterval), &
           &  this%sw_albedo_weights, use_bands=(.not. this%do_cloud_aerosol_per_sw_g_point))
    else
      ! Weight each wavenumber equally as in IFS Cycles 48 and earlier
      call this%gas_optics_sw%spectral_def%calc_mapping_from_bands( &
           &  this%sw_albedo_wavelength_bound(1:ninterval-1), this%i_sw_albedo_index(1:ninterval), &
           &  this%sw_albedo_weights, use_bands=(.not. this%do_cloud_aerosol_per_sw_g_point))
    end if

    ! Legacy method uses input band with largest weight
    if (this%do_nearest_spectral_sw_albedo) then
      allocate(this%i_albedo_from_band_sw(this%n_bands_sw))
      this%i_albedo_from_band_sw = maxloc(this%sw_albedo_weights, dim=1)
    end if

    if (this%iverbosesetup >= 2) then
      write(nulout, '(a)') 'Surface shortwave albedo'
      if (.not. this%do_nearest_spectral_sw_albedo) then
        call this%gas_optics_sw%spectral_def%print_mapping_from_bands(this%sw_albedo_weights, &
             &       use_bands=(.not. this%do_cloud_aerosol_per_sw_g_point))
      else if (ninterval <= 1) then
        write(nulout, '(a)') 'All shortwave bands will use the same albedo'
      else
        write(nulout, '(a,i0,a)',advance='no') 'Mapping from ', size(this%i_albedo_from_band_sw), &
             &  ' shortwave intervals to albedo intervals:'
        do jband = 1,size(this%i_albedo_from_band_sw)
          write(nulout,'(a,i0)',advance='no') ' ', this%i_albedo_from_band_sw(jband)
        end do
        write(nulout, '()')
      end if
    end if
    
  end subroutine consolidate_sw_albedo_intervals


  !---------------------------------------------------------------------
  ! Consolidate the surface longwave emissivity intervals with the
  ! band/g-point intervals
  subroutine consolidate_lw_emiss_intervals(this)

    use radiation_io, only : nulout
    use radiation_spectral_definition, only : TerrestrialReferenceTemperature

    class(config_type),   intent(inout) :: this

    integer :: ninterval, jint, jband

    ! Count the number of albedo/emissivity intervals
    ninterval = 0
    do jint = 1,NMaxAlbedoIntervals
      if (this%i_lw_emiss_index(jint) > 0) then
        ninterval = jint
      else
        exit
      end if
    end do

    if (ninterval < 1) then
      ! The user has not specified longwave emissivity bands - assume
      ! only one
      ninterval = 1
      this%i_lw_emiss_index(1) = 1
      this%i_lw_emiss_index(2:) = 0
      if (this%use_canopy_full_spectrum_sw) then
        this%n_canopy_bands_lw = this%n_g_lw
      else 
        this%n_canopy_bands_lw = 1
      end if
    else
      if (this%use_canopy_full_spectrum_lw) then
        this%n_canopy_bands_lw = this%n_g_lw
      else 
        this%n_canopy_bands_lw = maxval(this%i_lw_emiss_index(1:ninterval))
      end if
    end if

    if (this%do_weighted_surface_mapping) then
      call this%gas_optics_lw%spectral_def%calc_mapping_from_bands( &
           &  this%lw_emiss_wavelength_bound(1:ninterval-1), this%i_lw_emiss_index(1:ninterval), &
           &  this%lw_emiss_weights, use_bands=(.not. this%do_cloud_aerosol_per_lw_g_point))
    else
      ! Weight each wavenumber equally as in IFS Cycles 48 and earlier
      call this%gas_optics_lw%spectral_def%calc_mapping_from_bands( &
           &  this%lw_emiss_wavelength_bound(1:ninterval-1), this%i_lw_emiss_index(1:ninterval), &
           &  this%lw_emiss_weights, use_bands=(.not. this%do_cloud_aerosol_per_lw_g_point))
    end if

    ! Legacy method uses input band with largest weight
    if (this%do_nearest_spectral_lw_emiss) then
      allocate(this%i_emiss_from_band_lw(this%n_bands_lw))
      this%i_emiss_from_band_lw = maxloc(this%lw_emiss_weights, dim=1)
    end if

    if (this%iverbosesetup >= 2) then
      write(nulout, '(a)') 'Surface longwave emissivity'
      if (.not. this%do_nearest_spectral_lw_emiss) then
        call this%gas_optics_lw%spectral_def%print_mapping_from_bands(this%lw_emiss_weights, &
             &                          use_bands=(.not. this%do_cloud_aerosol_per_lw_g_point))
      else if (ninterval <= 1) then
        write(nulout, '(a)') 'All longwave bands will use the same emissivty'
      else
        write(nulout, '(a,i0,a)',advance='no') 'Mapping from ', size(this%i_emiss_from_band_lw), &
             &  ' longwave intervals to emissivity intervals:'
        do jband = 1,size(this%i_emiss_from_band_lw)
          write(nulout,'(a,i0)',advance='no') ' ', this%i_emiss_from_band_lw(jband)
        end do
        write(nulout, '()')
      end if
    end if

  end subroutine consolidate_lw_emiss_intervals


  !---------------------------------------------------------------------
  ! Return the 0-based index for str in enum_str, or abort if it is
  ! not found
  subroutine get_enum_code(str, enum_str, var_name, icode)

    use radiation_io, only : nulerr, radiation_abort

    character(len=*), intent(in)  :: str
    character(len=*), intent(in)  :: enum_str(0:)
    character(len=*), intent(in)  :: var_name
    integer,          intent(out) :: icode

    integer :: jc
    logical :: is_not_found

    ! If string is empty then we don't modify icode but assume it has
    ! a sensible default value
    if (len_trim(str) > 1) then
      is_not_found = .true.

      do jc = 0,size(enum_str)-1
        if (trim(str) == trim(enum_str(jc))) then
          icode = jc
          is_not_found = .false.
          exit
        end if
      end do
      if (is_not_found) then
        write(nulerr,'(a,a,a,a,a)',advance='no') '*** Error: ', trim(var_name), &
             &  ' must be one of: "', enum_str(0), '"'
        do jc = 1,size(enum_str)-1
          write(nulerr,'(a,a,a)',advance='no') ', "', trim(enum_str(jc)), '"'
        end do
        write(nulerr,'(a)') ''
        call radiation_abort('Radiation configuration error')
      end if
    end if

  end subroutine get_enum_code


  !---------------------------------------------------------------------
  ! Print one line of information: logical
  subroutine print_logical(message_str, name, val)
    use radiation_io, only : nulout
    character(len=*),   intent(in) :: message_str
    character(len=*),   intent(in) :: name
    logical,            intent(in) :: val
    character(4)                   :: on_or_off
    character(NPrintStringLen)     :: str
    if (val) then
      on_or_off = ' ON '
    else
      on_or_off = ' OFF'
    end if
    write(str, '(a,a4)') message_str, on_or_off
    write(nulout,'(a,a,a,a,l1,a)') str, ' (', name, '=', val,')'
  end subroutine print_logical


  !---------------------------------------------------------------------
  ! Print one line of information: integer
  subroutine print_integer(message_str, name, val)
    use radiation_io, only : nulout
    character(len=*),   intent(in) :: message_str
    character(len=*),   intent(in) :: name
    integer,            intent(in) :: val
    character(NPrintStringLen)     :: str
    write(str, '(a,a,i0)') message_str, ' = ', val
    write(nulout,'(a,a,a,a)') str, ' (', name, ')'
  end subroutine print_integer


  !---------------------------------------------------------------------
  ! Print one line of information: real
  subroutine print_real(message_str, name, val)
    use parkind1,     only : jprb
    use radiation_io, only : nulout
    character(len=*),   intent(in) :: message_str
    character(len=*),   intent(in) :: name
    real(jprb),         intent(in) :: val
    character(NPrintStringLen)     :: str
    write(str, '(a,a,g8.3)') message_str, ' = ', val
    write(nulout,'(a,a,a,a)') str, ' (', name, ')'
  end subroutine print_real


  !---------------------------------------------------------------------
  ! Print one line of information: enum
  subroutine print_enum(message_str, enum_str, name, val)
    use radiation_io, only : nulout
    character(len=*),   intent(in) :: message_str
    character(len=*),   intent(in) :: enum_str(0:)
    character(len=*),   intent(in) :: name
    integer,            intent(in) :: val
    character(NPrintStringLen)     :: str
    write(str, '(a,a,a,a)') message_str, ' "', trim(enum_str(val)), '"'
    write(nulout,'(a,a,a,a,i0,a)') str, ' (', name, '=', val,')'
  end subroutine print_enum

end module radiation_config
