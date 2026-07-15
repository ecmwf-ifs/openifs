! radiation_interface.F90 - Public interface to radiation scheme
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
!   2017-04-11  R. Hogan  Changes to enable generalized surface description
!   2017-09-08  R. Hogan  Reverted some changes
!
! To use the radiation scheme, create a configuration_type object,
! call "setup_radiation" on it once to load the look-up-tables and
! data describing how gas and hydrometeor absorption/scattering are to
! be represented, and call "radiation" multiple times on different
! input profiles.

module radiation_interface

  implicit none

  public  :: setup_radiation, set_gas_units, radiation
  private :: radiation_reverse

contains

  !---------------------------------------------------------------------
  ! Load the look-up-tables and data describing how gas and
  ! hydrometeor absorption/scattering are to be represented
  subroutine setup_radiation(config)

    use parkind1,         only : jprb
    use yomhook,          only : lhook, dr_hook, jphook
    use radiation_io,     only : nulerr, radiation_abort
    use radiation_config, only : config_type, ISolverMcICA, &
         &   IGasModelMonochromatic, IGasModelIFSRRTMG, IGasModelECCKD
    use radiation_spectral_definition, only &
         &  : SolarReferenceTemperature, TerrestrialReferenceTemperature
    ! Currently there are two gas absorption models: RRTMG (default)
    ! and monochromatic
    use radiation_monochromatic,  only : &
         &   setup_gas_optics_mono     => setup_gas_optics, &
         &   setup_cloud_optics_mono   => setup_cloud_optics, &
         &   setup_aerosol_optics_mono => setup_aerosol_optics
    use radiation_ifs_rrtm,       only :  setup_gas_optics_rrtmg => setup_gas_optics
    use radiation_ecckd_interface,only :  setup_gas_optics_ecckd => setup_gas_optics
    use radiation_cloud_optics,   only :  setup_cloud_optics
    use radiation_general_cloud_optics, only :  setup_general_cloud_optics
    use radiation_aerosol_optics, only :  setup_aerosol_optics

    
    type(config_type), intent(inout) :: config

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_interface:setup_radiation',0,hook_handle)

    ! Consolidate configuration data, including setting data file
    ! names
    call config%consolidate()

    ! Load the look-up tables from files in the specified directory
    if (config%i_gas_model == IGasModelMonochromatic) then
      call setup_gas_optics_mono(config, trim(config%directory_name))
    else if (config%i_gas_model == IGasModelIFSRRTMG) then
      call setup_gas_optics_rrtmg(config, trim(config%directory_name))
    else if (config%i_gas_model == IGasModelECCKD) then
      call setup_gas_optics_ecckd(config)
    end if

    if (config%do_lw_aerosol_scattering &
         & .and. .not. config%do_lw_cloud_scattering) then
      write(nulerr, '(a)') '*** Error: longwave aerosol scattering requires longwave cloud scattering'
      call radiation_abort('Radiation configuration error')
    end if

    
    ! Whether or not the "radiation" subroutine needs ssa_lw and g_lw
    ! arrays depends on whether longwave scattering by aerosols is to
    ! be included.  If not, one of the array dimensions will be set to
    ! zero.
    if (config%do_lw_aerosol_scattering) then
      config%n_g_lw_if_scattering = config%n_g_lw
    else
      config%n_g_lw_if_scattering = 0
    end if

    ! Whether or not the "radiation" subroutine needs ssa_lw_cloud and
    ! g_lw_cloud arrays depends on whether longwave scattering by
    ! hydrometeors is to be included.  If not, one of the array
    ! dimensions will be set to zero.
    if (config%do_lw_cloud_scattering) then
      config%n_bands_lw_if_scattering = config%n_bands_lw
    else
      config%n_bands_lw_if_scattering = 0
    end if

    ! If we have longwave scattering and McICA then even if there is
    ! no aerosol, it is convenient if single scattering albedo and
    ! g factor arrays are allocated before the call to
    ! solver_lw as they will be needed.
    if (config%do_lw_cloud_scattering &
         &  .and. config%i_solver_lw == ISolverMcICA) then
      config%n_g_lw_if_scattering = config%n_g_lw
    end if

    ! Consolidate the albedo/emissivity intervals with the shortwave
    ! and longwave spectral bands
    if (config%do_sw) then
      call config%consolidate_sw_albedo_intervals
    end if
    if (config%do_lw) then
      call config%consolidate_lw_emiss_intervals
    end if

    if (config%do_clouds) then
      if (config%i_gas_model == IGasModelMonochromatic) then
        !      call setup_cloud_optics_mono(config)
      elseif (config%use_general_cloud_optics) then
        call setup_general_cloud_optics(config)
      else
        call setup_cloud_optics(config)
      end if
    end if

    if (config%use_aerosols) then
      if (config%i_gas_model == IGasModelMonochromatic) then
!        call setup_aerosol_optics_mono(config)
      else 
        call setup_aerosol_optics(config)
      end if
    end if

    ! Load cloud water PDF look-up table for McICA
    if (         config%i_solver_sw == ISolverMcICA &
         &  .or. config%i_solver_lw == ISolverMcICA) then
      call config%pdf_sampler%setup(config%cloud_pdf_file_name, &
           &                        iverbose=config%iverbosesetup)
    end if

    if (lhook) call dr_hook('radiation_interface:setup_radiation',1,hook_handle)

  end subroutine setup_radiation


  !---------------------------------------------------------------------
  ! Scale the gas mixing ratios so that they have the units (and
  ! possibly scale factors) required by the specific gas absorption
  ! model.  This subroutine simply passes the gas object on to the
  ! module of the currently active gas model.
  subroutine set_gas_units(config, gas)
    
    use radiation_config
    use radiation_gas,             only : gas_type
    use radiation_monochromatic,   only : set_gas_units_mono  => set_gas_units
    use radiation_ifs_rrtm,        only : set_gas_units_ifs   => set_gas_units
    use radiation_ecckd_interface, only : set_gas_units_ecckd => set_gas_units

    type(config_type), intent(in)    :: config
    type(gas_type),    intent(inout) :: gas

    if (config%i_gas_model == IGasModelMonochromatic) then
      call set_gas_units_mono(gas)
    elseif (config%i_gas_model == IGasModelECCKD) then
      call set_gas_units_ecckd(gas)
    else
      call set_gas_units_ifs(gas)
    end if

  end subroutine set_gas_units


  !---------------------------------------------------------------------
  ! Run the radiation scheme according to the configuration in the
  ! config object. There are ncol profiles of which only istartcol to
  ! iendcol are to be processed, and there are nlev model levels.  The
  ! output fluxes are written to the flux object, and all other
  ! objects contain the input variables.  The variables may be defined
  ! either in order of increasing or decreasing pressure, but if in
  ! order of decreasing pressure then radiation_reverse will be called
  ! to reverse the order for the computation and then reverse the
  ! order of the output fluxes to match the inputs.
  subroutine radiation(ncol, nlev, istartcol, iendcol, config, &
       &  single_level, thermodynamics, gas, cloud, aerosol, flux)

    use parkind1,                 only : jprb
    use yomhook,                  only : lhook, dr_hook, jphook

    use radiation_io,             only : nulout
    use radiation_config,         only : config_type, &
         &   IGasModelMonochromatic, IGasModelIFSRRTMG, &
         &   ISolverMcICA, ISolverSpartacus, ISolverHomogeneous, &
         &   ISolverTripleclouds
    use radiation_single_level,   only : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_gas,            only : gas_type
    use radiation_cloud,          only : cloud_type
    use radiation_aerosol,        only : aerosol_type
    use radiation_flux,           only : flux_type
    use radiation_spartacus_sw,   only : solver_spartacus_sw
    use radiation_spartacus_lw,   only : solver_spartacus_lw
    use radiation_tripleclouds_sw,only : solver_tripleclouds_sw
    use radiation_tripleclouds_lw,only : solver_tripleclouds_lw
    use radiation_mcica_sw,       only : solver_mcica_sw
    use radiation_mcica_lw,       only : solver_mcica_lw
    use radiation_cloudless_sw,   only : solver_cloudless_sw
    use radiation_cloudless_lw,   only : solver_cloudless_lw
    use radiation_homogeneous_sw, only : solver_homogeneous_sw
    use radiation_homogeneous_lw, only : solver_homogeneous_lw
    use radiation_save,           only : save_radiative_properties

    ! Treatment of gas and hydrometeor optics 
    use radiation_monochromatic,  only : &
         &   gas_optics_mono         => gas_optics, &
         &   cloud_optics_mono       => cloud_optics, &
         &   add_aerosol_optics_mono => add_aerosol_optics
    use radiation_ifs_rrtm,       only : gas_optics_rrtmg => gas_optics
    use radiation_ecckd_interface,only : gas_optics_ecckd => gas_optics
    use radiation_cloud_optics,   only : cloud_optics
    use radiation_general_cloud_optics, only : general_cloud_optics
    use radiation_aerosol_optics, only : add_aerosol_optics

    ! Inputs
    integer, intent(in) :: ncol               ! number of columns
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in)   :: config
    type(single_level_type),  intent(in)   :: single_level
    type(thermodynamics_type),intent(in)   :: thermodynamics
    type(gas_type),           intent(in)   :: gas
    type(cloud_type),         intent(inout):: cloud
    type(aerosol_type),       intent(in)   :: aerosol
    ! Output
    type(flux_type),          intent(inout):: flux


    ! Local variables

    ! Layer optical depth, single scattering albedo and asymmetry factor of
    ! gases and aerosols at each longwave g-point, where the latter
    ! two variables are only defined if aerosol longwave scattering is
    ! enabled (otherwise both are treated as zero).
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: od_lw
    real(jprb), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: &
         &  ssa_lw, g_lw

    ! Layer in-cloud optical depth, single scattering albedo and
    ! asymmetry factor of hydrometeors in each longwave band, where
    ! the latter two variables are only defined if hydrometeor
    ! longwave scattering is enabled (otherwise both are treated as
    ! zero).
    real(jprb), dimension(config%n_bands_lw,nlev,istartcol:iendcol) :: od_lw_cloud
    real(jprb), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol) :: &
         &  ssa_lw_cloud, g_lw_cloud

    ! Layer optical depth, single scattering albedo and asymmetry factor of
    ! gases and aerosols at each shortwave g-point
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol) :: od_sw, ssa_sw, g_sw

    ! Layer in-cloud optical depth, single scattering albedo and
    ! asymmetry factor of hydrometeors in each shortwave band
    real(jprb), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
         &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud

    ! The Planck function (emitted flux from a black body) at half
    ! levels
    real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: planck_hl

    ! The longwave emission from and albedo of the surface in each
    ! longwave g-point; note that these are weighted averages of the
    ! values from individual tiles
    real(jprb), dimension(config%n_g_lw, istartcol:iendcol) :: lw_emission
    real(jprb), dimension(config%n_g_lw, istartcol:iendcol) :: lw_albedo

    ! Direct and diffuse shortwave surface albedo in each shortwave
    ! g-point; note that these are weighted averages of the values
    ! from individual tiles
    real(jprb), dimension(config%n_g_sw, istartcol:iendcol) :: sw_albedo_direct
    real(jprb), dimension(config%n_g_sw, istartcol:iendcol) :: sw_albedo_diffuse

    ! The incoming shortwave flux into a plane perpendicular to the
    ! incoming radiation at top-of-atmosphere in each of the shortwave
    ! g-points
    real(jprb), dimension(config%n_g_sw,istartcol:iendcol) :: incoming_sw

    character(len=100) :: rad_prop_file_name
    character(*), parameter :: rad_prop_base_file_name = "radiative_properties"

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_interface:radiation',0,hook_handle)

    if (thermodynamics%pressure_hl(istartcol,2) &
         &  < thermodynamics%pressure_hl(istartcol,1)) then
      ! Input arrays are arranged in order of decreasing pressure /
      ! increasing height: the following subroutine reverses them,
      ! call the radiation scheme and then reverses the returned
      ! fluxes
      call radiation_reverse(ncol, nlev, istartcol, iendcol, config, &
           &  single_level, thermodynamics, gas, cloud, aerosol, flux)
    else

      ! Input arrays arranged in order of increasing pressure /
      ! decreasing height: progress normally

      ! Extract surface albedos at each gridpoint
      call single_level%get_albedos(istartcol, iendcol, config, &
           &                        sw_albedo_direct, sw_albedo_diffuse, &
           &                        lw_albedo)

      ! Compute gas absorption optical depth in shortwave and
      ! longwave, shortwave single scattering albedo (i.e. fraction of
      ! extinction due to Rayleigh scattering), Planck functions and
      ! incoming shortwave flux at each g-point, for the specified
      ! range of atmospheric columns
      if (config%i_gas_model == IGasModelMonochromatic) then
        call gas_optics_mono(ncol,nlev,istartcol,iendcol, config, &
             &  single_level, thermodynamics, gas, lw_albedo, &
             &  od_lw, od_sw, ssa_sw, &
             &  planck_hl, lw_emission, incoming_sw)
      else if (config%i_gas_model == IGasModelIFSRRTMG) then
        call gas_optics_rrtmg(ncol,nlev,istartcol,iendcol, config, &
             &  single_level, thermodynamics, gas, &
             &  od_lw, od_sw, ssa_sw, lw_albedo=lw_albedo, &
             &  planck_hl=planck_hl, lw_emission=lw_emission, &
             &  incoming_sw=incoming_sw)
      else
        call gas_optics_ecckd(ncol,nlev,istartcol,iendcol, config, &
             &  single_level, thermodynamics, gas, &
             &  od_lw, od_sw, ssa_sw, lw_albedo=lw_albedo, &
             &  planck_hl=planck_hl, lw_emission=lw_emission, &
             &  incoming_sw=incoming_sw)
      end if

      if (config%do_clouds) then
        ! Crop the cloud fraction to remove clouds that have too small
        ! a fraction or water content; after this, we can safely
        ! assume that a cloud is present if cloud%fraction > 0.0.
        call cloud%crop_cloud_fraction(istartcol, iendcol, &
             &            config%cloud_fraction_threshold, &
             &            config%cloud_mixing_ratio_threshold)

        ! Compute hydrometeor absorption/scattering properties in each
        ! shortwave and longwave band
        if (config%i_gas_model == IGasModelMonochromatic) then
          call cloud_optics_mono(nlev, istartcol, iendcol, &
               &  config, thermodynamics, cloud, &
               &  od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
               &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud)
        elseif (config%use_general_cloud_optics) then
          call general_cloud_optics(nlev, istartcol, iendcol, &
               &  config, thermodynamics, cloud, & 
               &  od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
               &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud)
        else
          call cloud_optics(nlev, istartcol, iendcol, &
               &  config, thermodynamics, cloud, & 
               &  od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
               &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud)
        end if
      end if ! do_clouds

      if (config%use_aerosols) then
        if (config%i_gas_model == IGasModelMonochromatic) then
!          call add_aerosol_optics_mono(nlev,istartcol,iendcol, &
!               &  config, thermodynamics, gas, aerosol, & 
!               &  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)
        else
          call add_aerosol_optics(nlev,istartcol,iendcol, &
               &  config, thermodynamics, gas, aerosol, & 
               &  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)
        end if
      else
        g_sw(:,:,istartcol:iendcol) = 0.0_jprb
        if (config%do_lw_aerosol_scattering) then
          ssa_lw(:,:,istartcol:iendcol) = 0.0_jprb
          g_lw(:,:,istartcol:iendcol)   = 0.0_jprb
        end if
      end if

      ! For diagnostic purposes, save these intermediate variables to
      ! a NetCDF file
      if (config%do_save_radiative_properties) then
        if (istartcol == 1 .and. iendcol == ncol) then
          rad_prop_file_name = rad_prop_base_file_name // ".nc"
        else
          write(rad_prop_file_name,'(a,a,i4.4,a,i4.4,a)') &
               &  rad_prop_base_file_name, '_', istartcol, '-',iendcol,'.nc'
        end if
        call save_radiative_properties(trim(rad_prop_file_name), &
             &  nlev, istartcol, iendcol, &
             &  config, single_level, thermodynamics, cloud, &
             &  planck_hl, lw_emission, lw_albedo, &
             &  sw_albedo_direct, sw_albedo_diffuse, incoming_sw, &
             &  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw, &
             &  od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
             &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud)
      end if

      if (config%do_lw) then
        if (config%iverbose >= 2) then
          write(nulout,'(a)') 'Computing longwave fluxes'
        end if

        if (config%i_solver_lw == ISolverMcICA) then
          ! Compute fluxes using the McICA longwave solver
          call solver_mcica_lw(nlev,istartcol,iendcol, &
               &  config, single_level, cloud, & 
               &  od_lw, ssa_lw, g_lw, od_lw_cloud, ssa_lw_cloud, &
               &  g_lw_cloud, planck_hl, lw_emission, lw_albedo, flux)
        else if (config%i_solver_lw == ISolverSPARTACUS) then
          ! Compute fluxes using the SPARTACUS longwave solver
          call solver_spartacus_lw(nlev,istartcol,iendcol, &
               &  config, thermodynamics, cloud, & 
               &  od_lw, ssa_lw, g_lw, od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
               &  planck_hl, lw_emission, lw_albedo, flux)
        else if (config%i_solver_lw == ISolverTripleclouds) then
          ! Compute fluxes using the Tripleclouds longwave solver
          call solver_tripleclouds_lw(nlev,istartcol,iendcol, &
               &  config, cloud, & 
               &  od_lw, ssa_lw, g_lw, od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
               &  planck_hl, lw_emission, lw_albedo, flux)
        elseif (config%i_solver_lw == ISolverHomogeneous) then
          ! Compute fluxes using the homogeneous solver
          call solver_homogeneous_lw(nlev,istartcol,iendcol, &
               &  config, cloud, & 
               &  od_lw, ssa_lw, g_lw, od_lw_cloud, ssa_lw_cloud, &
               &  g_lw_cloud, planck_hl, lw_emission, lw_albedo, flux)
        else
          ! Compute fluxes using the cloudless solver
          call solver_cloudless_lw(nlev,istartcol,iendcol, &
               &  config, od_lw, ssa_lw, g_lw, &
               &  planck_hl, lw_emission, lw_albedo, flux)
        end if
      end if

      if (config%do_sw) then
        if (config%iverbose >= 2) then
          write(nulout,'(a)') 'Computing shortwave fluxes'
        end if

        if (config%i_solver_sw == ISolverMcICA) then
          ! Compute fluxes using the McICA shortwave solver
          call solver_mcica_sw(nlev,istartcol,iendcol, &
               &  config, single_level, cloud, & 
               &  od_sw, ssa_sw, g_sw, od_sw_cloud, ssa_sw_cloud, &
               &  g_sw_cloud, sw_albedo_direct, sw_albedo_diffuse, &
               &  incoming_sw, flux)
        else if (config%i_solver_sw == ISolverSPARTACUS) then
          ! Compute fluxes using the SPARTACUS shortwave solver
          call solver_spartacus_sw(nlev,istartcol,iendcol, &
               &  config, single_level, thermodynamics, cloud, & 
               &  od_sw, ssa_sw, g_sw, od_sw_cloud, ssa_sw_cloud, &
               &  g_sw_cloud, sw_albedo_direct, sw_albedo_diffuse, &
               &  incoming_sw, flux)
        else if (config%i_solver_sw == ISolverTripleclouds) then
          ! Compute fluxes using the Tripleclouds shortwave solver
          call solver_tripleclouds_sw(nlev,istartcol,iendcol, &
               &  config, single_level, cloud, & 
               &  od_sw, ssa_sw, g_sw, od_sw_cloud, ssa_sw_cloud, &
               &  g_sw_cloud, sw_albedo_direct, sw_albedo_diffuse, &
               &  incoming_sw, flux)
        elseif (config%i_solver_sw == ISolverHomogeneous) then
          ! Compute fluxes using the homogeneous solver
          call solver_homogeneous_sw(nlev,istartcol,iendcol, &
               &  config, single_level, cloud, & 
               &  od_sw, ssa_sw, g_sw, od_sw_cloud, ssa_sw_cloud, &
               &  g_sw_cloud, sw_albedo_direct, sw_albedo_diffuse, &
               &  incoming_sw, flux)
        else
          ! Compute fluxes using the cloudless solver
          call solver_cloudless_sw(nlev,istartcol,iendcol, &
               &  config, single_level, od_sw, ssa_sw, g_sw, &
               &  sw_albedo_direct, sw_albedo_diffuse, &
               &  incoming_sw, flux)
        end if
      end if

      ! Store surface downwelling, and TOA, fluxes in bands from
      ! fluxes in g points
      call flux%calc_surface_spectral(config, istartcol, iendcol)
      call flux%calc_toa_spectral    (config, istartcol, iendcol)

    end if
    
    if (lhook) call dr_hook('radiation_interface:radiation',1,hook_handle)

  end subroutine radiation


  !---------------------------------------------------------------------
  ! If the input arrays are arranged in order of decreasing pressure /
  ! increasing height then this subroutine reverses them, calls the
  ! radiation scheme and then reverses the returned fluxes. Since this
  ! subroutine calls, and is called by "radiation", it must be in this
  ! module to avoid circular dependencies.
  subroutine radiation_reverse(ncol, nlev, istartcol, iendcol, config, &
       &  single_level, thermodynamics, gas, cloud, aerosol, flux)
 
    use parkind1, only : jprb

    use radiation_io,             only : nulout
    use radiation_config,         only : config_type
    use radiation_single_level,   only : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_gas,            only : gas_type
    use radiation_cloud,          only : cloud_type
    use radiation_aerosol,        only : aerosol_type
    use radiation_flux,           only : flux_type

    ! Inputs
    integer, intent(in) :: ncol               ! number of columns
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(gas_type),           intent(in) :: gas
    type(cloud_type),         intent(in) :: cloud
    type(aerosol_type),       intent(in) :: aerosol
    ! Output
    type(flux_type),          intent(inout):: flux

    ! Reversed data structures
    type(thermodynamics_type) :: thermodynamics_rev
    type(gas_type)            :: gas_rev
    type(cloud_type)          :: cloud_rev
    type(aerosol_type)        :: aerosol_rev
    type(flux_type)           :: flux_rev

    ! Start and end levels for aerosol data
    integer :: istartlev, iendlev

    if (config%iverbose >= 2) then
      write(nulout,'(a)') 'Reversing arrays to be in order of increasing pressure'
    end if

    ! Allocate reversed arrays
    call thermodynamics_rev%allocate(ncol, nlev)
    call cloud_rev%allocate(ncol, nlev)
    call flux_rev%allocate(config, istartcol, iendcol, nlev)
    if (allocated(aerosol%mixing_ratio)) then
      istartlev = nlev + 1 - aerosol%iendlev
      iendlev   = nlev + 1 - aerosol%istartlev
      call aerosol_rev%allocate(ncol, istartlev, iendlev, &
           &                    config%n_aerosol_types)
    end if

    ! Fill reversed thermodynamic arrays
    thermodynamics_rev%pressure_hl(istartcol:iendcol,:) &
         &  = thermodynamics%pressure_hl(istartcol:iendcol, nlev+1:1:-1)
    thermodynamics_rev%temperature_hl(istartcol:iendcol,:) &
         &  = thermodynamics%temperature_hl(istartcol:iendcol, nlev+1:1:-1)

    ! Fill reversed gas arrays
    call gas%reverse(istartcol, iendcol, gas_rev)

    if (config%do_clouds) then
      ! Fill reversed cloud arrays
      cloud_rev%q_liq(istartcol:iendcol,:) &
           &  = cloud%q_liq(istartcol:iendcol,nlev:1:-1)
      cloud_rev%re_liq(istartcol:iendcol,:) &
           &  = cloud%re_liq(istartcol:iendcol,nlev:1:-1)
      cloud_rev%q_ice(istartcol:iendcol,:) &
           &  = cloud%q_ice(istartcol:iendcol,nlev:1:-1)
      cloud_rev%re_ice(istartcol:iendcol,:) &
           &  = cloud%re_ice(istartcol:iendcol,nlev:1:-1)
      cloud_rev%fraction(istartcol:iendcol,:) &
           &  = cloud%fraction(istartcol:iendcol,nlev:1:-1)
      cloud_rev%overlap_param(istartcol:iendcol,:) &
           &  = cloud%overlap_param(istartcol:iendcol,nlev-1:1:-1)
      if (allocated(cloud%fractional_std)) then
        cloud_rev%fractional_std(istartcol:iendcol,:) &
             &  = cloud%fractional_std(istartcol:iendcol,nlev:1:-1)
      else
        cloud_rev%fractional_std(istartcol:iendcol,:) = 0.0_jprb       
      end if
      if (allocated(cloud%inv_cloud_effective_size)) then
        cloud_rev%inv_cloud_effective_size(istartcol:iendcol,:) &
             &  = cloud%inv_cloud_effective_size(istartcol:iendcol,nlev:1:-1)
      else
        cloud_rev%inv_cloud_effective_size(istartcol:iendcol,:) = 0.0_jprb
      end if
    end if

    if (allocated(aerosol%mixing_ratio)) then
      aerosol_rev%mixing_ratio(:,istartlev:iendlev,:) &
           &  = aerosol%mixing_ratio(:,aerosol%iendlev:aerosol%istartlev:-1,:)
    end if

    ! Run radiation scheme on reversed profiles
    call radiation(ncol, nlev,istartcol,iendcol, &
         &  config, single_level, thermodynamics_rev, gas_rev, &
         &  cloud_rev, aerosol_rev, flux_rev)

    ! Reorder fluxes
    if (allocated(flux%lw_up)) then
      flux%lw_up(istartcol:iendcol,:) &
           &  = flux_rev%lw_up(istartcol:iendcol,nlev+1:1:-1)
      flux%lw_dn(istartcol:iendcol,:) &
           &  = flux_rev%lw_dn(istartcol:iendcol,nlev+1:1:-1)
      if (allocated(flux%lw_up_clear)) then
        flux%lw_up_clear(istartcol:iendcol,:) &
             &  = flux_rev%lw_up_clear(istartcol:iendcol,nlev+1:1:-1)
        flux%lw_dn_clear(istartcol:iendcol,:) &
             &  = flux_rev%lw_dn_clear(istartcol:iendcol,nlev+1:1:-1)
      end if
    end if
    if (allocated(flux%sw_up)) then
      flux%sw_up(istartcol:iendcol,:) &
           &  = flux_rev%sw_up(istartcol:iendcol,nlev+1:1:-1)
      flux%sw_dn(istartcol:iendcol,:) &
           &  = flux_rev%sw_dn(istartcol:iendcol,nlev+1:1:-1)
      if (allocated(flux%sw_dn_direct)) then
        flux%sw_dn_direct(istartcol:iendcol,:) &
             &  = flux_rev%sw_dn_direct(istartcol:iendcol,nlev+1:1:-1)
      end if
      if (allocated(flux%sw_up_clear)) then
        flux%sw_up_clear(istartcol:iendcol,:) &
             &  = flux_rev%sw_up_clear(istartcol:iendcol,nlev+1:1:-1)
        flux%sw_dn_clear(istartcol:iendcol,:) &
             &  = flux_rev%sw_dn_clear(istartcol:iendcol,nlev+1:1:-1)
        if (allocated(flux%sw_dn_direct_clear)) then
          flux%sw_dn_direct_clear(istartcol:iendcol,:) &
               &  = flux_rev%sw_dn_direct_clear(istartcol:iendcol,nlev+1:1:-1)
        end if
      end if
    end if

    ! Deallocate reversed arrays
    call thermodynamics_rev%deallocate
    call gas_rev%deallocate
    call cloud_rev%deallocate
    call flux_rev%deallocate
    if (allocated(aerosol%mixing_ratio)) then
      call aerosol_rev%deallocate
    end if

  end subroutine radiation_reverse

end module radiation_interface
