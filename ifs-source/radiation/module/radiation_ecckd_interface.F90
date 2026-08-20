! radiation_ecckd_interface.F90 - Interface to ecCKD gas optics model
!
! (C) Copyright 2020- ECMWF.
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
! License: see the COPYING file for details
!

module radiation_ecckd_interface

  implicit none

  public  :: setup_gas_optics, set_gas_units, gas_optics !, planck_function

contains

  !---------------------------------------------------------------------
  ! Setup the ecCKD generalized gas optics model
  subroutine setup_gas_optics(config)

    use parkind1, only : jprb
    use radiation_config
    use yomhook,  only : lhook, dr_hook, jphook

    type(config_type), intent(inout), target :: config

    integer :: jj
    
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ecckd_interface:setup_gas_optics',0,hook_handle)

    if (config%do_sw) then

      ! Read shortwave ecCKD gas optics NetCDF file
      call config%gas_optics_sw%read(trim(config%gas_optics_sw_file_name), &
           &                         config%iverbosesetup)

      ! Copy over relevant properties
      config%n_g_sw     = config%gas_optics_sw%ng

      if (config%do_cloud_aerosol_per_sw_g_point) then
        ! Bands and g points are the same
        config%n_bands_sw = config%n_g_sw
      else
        ! Bands are groups of g points and span a continuous region of
        ! wavenumber space
        config%n_bands_sw = config%gas_optics_sw%spectral_def%nband
      end if

      allocate(config%i_band_from_g_sw          (config%n_g_sw))
      allocate(config%i_band_from_reordered_g_sw(config%n_g_sw))
      allocate(config%i_g_from_reordered_g_sw   (config%n_g_sw))
        
      if (config%do_cloud_aerosol_per_sw_g_point) then
        config%i_band_from_g_sw           = [ (jj, jj = 1,config%n_g_sw) ]
        config%i_band_from_reordered_g_sw = [ (jj, jj = 1,config%n_g_sw) ]
      else
        config%i_band_from_g_sw &
             &  = config%gas_optics_sw%spectral_def%i_band_number
        config%i_band_from_reordered_g_sw &
             &  = config%gas_optics_sw%spectral_def%i_band_number
      end if
      config%i_g_from_reordered_g_sw      = [ (jj, jj = 1,config%n_g_sw) ]

      if (config%iverbosesetup >= 2) then
        call config%gas_optics_sw%print()
      end if

      ! Override solar spectral irradiance, if filename provided
      if (config%use_spectral_solar_cycle) then
        call config%gas_optics_sw%read_spectral_solar_cycle(config%ssi_file_name, &
             &           config%iverbosesetup, config%use_updated_solar_spectrum)
      end if

    end if

    if (config%do_lw) then

      ! Read longwave ecCKD gas optics NetCDF file
      call config%gas_optics_lw%read(trim(config%gas_optics_lw_file_name), &
           &                         config%iverbosesetup)

      ! Copy over relevant properties
      config%n_g_lw     = config%gas_optics_lw%ng

      if (config%do_cloud_aerosol_per_lw_g_point) then
        ! Bands and g points are the same
        config%n_bands_lw = config%n_g_lw
      else
        ! Bands are groups of g points and span a continuous region of
        ! wavenumber space
        config%n_bands_lw = config%gas_optics_lw%spectral_def%nband
      end if

      allocate(config%i_band_from_g_lw          (config%n_g_lw))
      allocate(config%i_band_from_reordered_g_lw(config%n_g_lw))
      allocate(config%i_g_from_reordered_g_lw   (config%n_g_lw))

      if (config%do_cloud_aerosol_per_lw_g_point) then
        config%i_band_from_g_lw           = [ (jj, jj = 1,config%n_g_lw) ]
        config%i_band_from_reordered_g_lw = [ (jj, jj = 1,config%n_g_lw) ]
      else
        config%i_band_from_g_lw &
             &  = config%gas_optics_lw%spectral_def%i_band_number
        config%i_band_from_reordered_g_lw &
             &  = config%gas_optics_lw%spectral_def%i_band_number
      end if
      config%i_g_from_reordered_g_lw      = [ (jj, jj = 1,config%n_g_lw) ]

      if (config%iverbosesetup >= 2) then
        call config%gas_optics_lw%print()
      end if

    end if

    ! The i_spec_* variables are used solely for storing spectral
    ! data, and this can either be by band or by g-point
    if (config%do_save_spectral_flux) then
      if (config%do_save_gpoint_flux) then
        config%n_spec_sw = config%n_g_sw
        config%n_spec_lw = config%n_g_lw
        config%i_spec_from_reordered_g_sw => config%i_g_from_reordered_g_sw
        config%i_spec_from_reordered_g_lw => config%i_g_from_reordered_g_lw
      else
        config%n_spec_sw = config%n_bands_sw
        config%n_spec_lw = config%n_bands_lw
        config%i_spec_from_reordered_g_sw => config%i_band_from_reordered_g_sw
        config%i_spec_from_reordered_g_lw => config%i_band_from_reordered_g_lw
      end if
    else
      config%n_spec_sw = 0
      config%n_spec_lw = 0
      nullify(config%i_spec_from_reordered_g_sw)
      nullify(config%i_spec_from_reordered_g_lw)
    end if

    if (lhook) call dr_hook('radiation_ecckd_interface:setup_gas_optics',1,hook_handle)
    
  end subroutine setup_gas_optics


  !---------------------------------------------------------------------
  ! Scale gas mixing ratios according to required units
  subroutine set_gas_units(gas)

    use radiation_gas, only : gas_type, IVolumeMixingRatio
    use yomhook,       only : lhook, dr_hook, jphook
    
    type(gas_type),    intent(inout) :: gas

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ecckd_interface:set_gas_units',0,hook_handle)

    call gas%set_units(IVolumeMixingRatio)

    if (lhook) call dr_hook('radiation_ecckd_interface:set_gas_units',1,hook_handle)

  end subroutine set_gas_units


  !---------------------------------------------------------------------
  ! Compute gas optical depths, shortwave scattering, Planck function
  ! and incoming shortwave radiation at top-of-atmosphere
  subroutine gas_optics(ncol,nlev,istartcol,iendcol, &
       &  config, single_level, thermodynamics, gas, & 
       &  od_lw, od_sw, ssa_sw, lw_albedo, planck_hl, lw_emission, &
       &  incoming_sw)

    use parkind1, only : jprb
    use yomhook,  only : lhook, dr_hook, jphook

    use radiation_config,         only : config_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_single_level,   only : single_level_type
    use radiation_gas_constants,  only : NMaxGases
    use radiation_gas

    integer,                  intent(in) :: ncol               ! number of columns
    integer,                  intent(in) :: nlev               ! number of levels
    integer,                  intent(in) :: istartcol, iendcol ! range of cols to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(gas_type),           intent(in) :: gas

    ! Longwave albedo of the surface
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
         &  intent(in), optional :: lw_albedo

    ! Gaseous layer optical depth in longwave and shortwave, and
    ! shortwave single scattering albedo (i.e. fraction of extinction
    ! due to Rayleigh scattering) at each g-point
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_lw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_sw, ssa_sw

    ! The Planck function (emitted flux from a black body) at half
    ! levels at each longwave g-point
    real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), &
         &   intent(out), optional :: planck_hl
    ! Planck function for the surface (W m-2)
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
         &   intent(out), optional :: lw_emission

    ! The incoming shortwave flux into a plane perpendicular to the
    ! incoming radiation at top-of-atmosphere in each of the shortwave
    ! g-points
    real(jprb), dimension(config%n_g_sw,istartcol:iendcol), &
         &   intent(out), optional :: incoming_sw

    ! Temperature at full levels (K)
    real(jprb) :: temperature_fl(istartcol:iendcol,nlev)

    integer :: jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ecckd_interface:gas_optics',0,hook_handle)

    !temperature_fl(istartcol:iendcol,:) &
    !     &  = 0.5_jprb * (thermodynamics%temperature_hl(istartcol:iendcol,1:nlev) &
    !     &               +thermodynamics%temperature_hl(istartcol:iendcol,2:nlev+1))
 
    temperature_fl(istartcol:iendcol,:) &
         &  = (thermodynamics%temperature_hl(istartcol:iendcol,1:nlev) &
         &     *thermodynamics%pressure_hl(istartcol:iendcol,1:nlev) &
         &    +thermodynamics%temperature_hl(istartcol:iendcol,2:nlev+1) &
         &     *thermodynamics%pressure_hl(istartcol:iendcol,2:nlev+1)) &
         &  / (thermodynamics%pressure_hl(istartcol:iendcol,1:nlev) &
         &    +thermodynamics%pressure_hl(istartcol:iendcol,2:nlev+1))
 
    if (config%do_sw) then

      call config%gas_optics_sw%calc_optical_depth(ncol,nlev,istartcol,iendcol, &
           &  NMaxGases, thermodynamics%pressure_hl, &
           &  temperature_fl, &
           &  gas%mixing_ratio, &
!           &  reshape(gas%mixing_ratio(istartcol:iendcol,:,:), &
!           &          [nlev,iendcol-istartcol+1,NMaxGases],order=[2,1,3]), &
           &  od_sw, rayleigh_od_fl=ssa_sw)
      ! At this point od_sw = absorption optical depth and ssa_sw =
      ! rayleigh optical depth: convert to total optical depth and
      ! single-scattering albedo
      od_sw = od_sw + ssa_sw
      ssa_sw = ssa_sw / od_sw

      if (present(incoming_sw)) then
        if (single_level%spectral_solar_cycle_multiplier == 0.0_jprb) then
          call config%gas_optics_sw%calc_incoming_sw(single_level%solar_irradiance, &
               &        incoming_sw)
        else
          call config%gas_optics_sw%calc_incoming_sw(single_level%solar_irradiance, &
               &        incoming_sw, single_level%spectral_solar_cycle_multiplier)
        end if
      end if

    end if

    if (config%do_lw) then

      call config%gas_optics_lw%calc_optical_depth(ncol,nlev,istartcol,iendcol, &
           &  NMaxGases, thermodynamics%pressure_hl, &
           &  temperature_fl, &
           &  gas%mixing_ratio, &
!           &  reshape(gas%mixing_ratio(istartcol:iendcol,:,:), &
!           &          [nlev,iendcol-istartcol+1,NMaxGases],order=[2,1,3]), &
           &  od_lw)

      ! Calculate the Planck function for each g point
      do jcol = istartcol,iendcol
        call config%gas_optics_lw%calc_planck_function(nlev+1, &
             &  thermodynamics%temperature_hl(jcol,:), planck_hl(:,:,jcol))
      end do
      call config%gas_optics_lw%calc_planck_function(iendcol+1-istartcol, &
           &  single_level%skin_temperature(istartcol:iendcol), &
           &  lw_emission(:,:))
      lw_emission = lw_emission * (1.0_jprb - lw_albedo)

    end if

    if (lhook) call dr_hook('radiation_ecckd_interface:gas_optics',1,hook_handle)
    
  end subroutine gas_optics

  ! !---------------------------------------------------------------------
  ! ! Externally facing function for computing the Planck function
  ! ! without reference to any gas profile; typically this would be used
  ! ! for computing the emission by a surface.
  ! subroutine planck_function(config, temperature, planck_surf)

  !   use parkind1,                 only : jprb
  !   use radiation_config,         only : config_type

  !   type(config_type), intent(in) :: config
  !   real(jprb),        intent(in) :: temperature

  !   ! Planck function of the surface (W m-2)
  !   real(jprb), dimension(config%n_g_lw), intent(out) :: planck_surf

  ! end subroutine planck_function

end module radiation_ecckd_interface
