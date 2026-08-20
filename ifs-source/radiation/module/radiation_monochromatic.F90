! radiation_interface.F90 - Monochromatic gas/cloud optics for testing
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
!   2017-04-11  R. Hogan  Receive "surface" dummy argument
!   2017-09-13  R. Hogan  Revert
!   2018-08-29  R. Hogan  Particulate single-scattering albedo / asymmetry from namelist

module radiation_monochromatic

  implicit none

  public  :: setup_gas_optics, gas_optics, set_gas_units, &
       &     setup_cloud_optics, cloud_optics,            &
       &     setup_aerosol_optics, add_aerosol_optics

contains

  ! Provides elemental function "delta_eddington"
#include "radiation_delta_eddington.h"

  !---------------------------------------------------------------------
  ! Setup the arrays in the config object corresponding to the
  ! monochromatic gas optics model.  The directory argument is not
  ! used, since no look-up tables need to be loaded.
  subroutine setup_gas_optics(config, directory)

    use radiation_config, only : config_type
    
    type(config_type), intent(inout) :: config
    character(len=*),  intent(in)    :: directory

    ! In the monochromatic model we have simply one band and g-point
    ! in both the longwave and shortwave parts of the spectrum
    config%n_g_sw     = 1
    config%n_g_lw     = 1
    config%n_bands_sw = 1
    config%n_bands_lw = 1

    ! Allocate arrays
    allocate(config%i_band_from_g_sw          (config%n_g_sw))
    allocate(config%i_band_from_g_lw          (config%n_g_lw))
    allocate(config%i_band_from_reordered_g_sw(config%n_g_sw))
    allocate(config%i_band_from_reordered_g_lw(config%n_g_lw))
    allocate(config%i_g_from_reordered_g_sw(config%n_g_sw))
    allocate(config%i_g_from_reordered_g_lw(config%n_g_lw))

    ! Indices are trivial...
    config%i_band_from_g_sw           = 1
    config%i_band_from_g_lw           = 1
    config%i_g_from_reordered_g_sw    = 1
    config%i_g_from_reordered_g_lw    = 1
    config%i_band_from_reordered_g_sw = 1
    config%i_band_from_reordered_g_lw = 1

  end subroutine setup_gas_optics


  !---------------------------------------------------------------------
  ! Dummy routine for scaling gas mixing ratios
  subroutine set_gas_units(gas)

    use radiation_gas,           only : gas_type
    type(gas_type),    intent(inout) :: gas

  end subroutine set_gas_units


  !---------------------------------------------------------------------
  ! Dummy setup routine for cloud optics: in fact, no setup is
  ! required for monochromatic case
  subroutine setup_cloud_optics(config)

    use radiation_config, only : config_type
    type(config_type), intent(inout) :: config

  end subroutine setup_cloud_optics


  !---------------------------------------------------------------------
  ! Dummy subroutine since no aerosols are represented in
  ! monochromatic case
  subroutine setup_aerosol_optics(config)

    use radiation_config,              only : config_type
    type(config_type), intent(inout) :: config

  end subroutine setup_aerosol_optics


  !---------------------------------------------------------------------
  ! Compute gas optical depths, shortwave scattering, Planck function
  ! and incoming shortwave radiation at top-of-atmosphere
  subroutine gas_optics(ncol,nlev,istartcol,iendcol, &
       config, single_level, thermodynamics, gas, lw_albedo, & 
       od_lw, od_sw, ssa_sw, planck_hl, lw_emission, &
       incoming_sw)

    use parkind1,                 only : jprb
    use radiation_config,         only : config_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_single_level,   only : single_level_type
    use radiation_gas,            only : gas_type
    use radiation_constants,      only : Pi, StefanBoltzmann

    ! Inputs
    integer, intent(in) :: ncol               ! number of columns
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(gas_type),           intent(in) :: gas

    ! Longwave albedo of the surface
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
         &  intent(in) :: lw_albedo

    ! Outputs

    ! Gaseous layer optical depth in longwave and shortwave, and
    ! shortwave single scattering albedo (i.e. fraction of extinction
    ! due to Rayleigh scattering) at each g-point
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_lw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_sw, ssa_sw

    ! The Planck function (emitted flux from a black body) at half
    ! levels and at the surface at each longwave g-point
    real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), intent(out) :: &
         &   planck_hl
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), intent(out) :: &
         &   lw_emission

    ! The incoming shortwave flux into a plane perpendicular to the
    ! incoming radiation at top-of-atmosphere in each of the shortwave
    ! g-points
    real(jprb), dimension(config%n_g_sw,istartcol:iendcol), intent(out) :: &
         &   incoming_sw
    
    ! Ratio of the optical depth of the entire atmospheric column that
    ! is in the current layer
    real(jprb), dimension(istartcol:iendcol) :: extinction_fraction

    ! In the monochromatic model, the absorption by the atmosphere is
    ! assumed proportional to the mass in each layer, so is defined in
    ! terms of a total zenith optical depth and then distributed with
    ! height according to the pressure.
    !real(jprb), parameter :: total_od_sw = 0.10536_jprb ! Transmittance 0.9
    !real(jprb), parameter :: total_od_lw = 1.6094_jprb  ! Transmittance 0.2

    integer :: jlev

    do jlev = 1,nlev
      ! The fraction of the total optical depth in the current layer
      ! is proportional to the fraction of the mass of the atmosphere
      ! in the current layer, computed from pressure assuming
      ! hydrostatic balance
      extinction_fraction = &
           &   (thermodynamics%pressure_hl(istartcol:iendcol,jlev+1) &
           &   -thermodynamics%pressure_hl(istartcol:iendcol,jlev)) &
           &   /thermodynamics%pressure_hl(istartcol:iendcol,nlev)
      
      ! Assign longwave and shortwave optical depths
      od_lw(1,jlev,:) = config%mono_lw_total_od*extinction_fraction
      od_sw(1,jlev,:) = config%mono_sw_total_od*extinction_fraction
    end do

    ! Shortwave single-scattering albedo is almost entirely Rayleigh
    ! scattering
    ssa_sw = 0.999999_jprb

    ! Entire shortwave spectrum represented in one band
    incoming_sw(1,:) = single_level%solar_irradiance

    if (single_level%is_simple_surface) then
      if (config%mono_lw_wavelength <= 0.0_jprb) then
        ! Entire longwave spectrum represented in one band
        lw_emission(1,istartcol:iendcol) &
             &  = StefanBoltzmann * single_level%skin_temperature(istartcol:iendcol)**4 &
             &  * single_level%lw_emissivity(istartcol:iendcol,1)
        do jlev = 1,nlev+1
          planck_hl(1,jlev,istartcol:iendcol) = StefanBoltzmann * thermodynamics%temperature_hl(istartcol:iendcol,jlev)**4
        end do
      else
        ! Single wavelength: multiply by pi to convert W sr-1 m-3 to W m-3
        lw_emission(1,istartcol:iendcol) = Pi*planck_function(config%mono_lw_wavelength, &
             &             single_level%skin_temperature(istartcol:iendcol)) &
             &  * single_level%lw_emissivity(istartcol:iendcol,1)
        do jlev = 1,nlev+1
          planck_hl(1,jlev,istartcol:iendcol) = Pi*planck_function(config%mono_lw_wavelength, &
               &             thermodynamics%temperature_hl(istartcol:iendcol,jlev))
        end do
      end if
    else
      lw_emission = transpose(single_level%lw_emission)
    end if

  end subroutine gas_optics


  !---------------------------------------------------------------------
  ! Compute cloud optical depth, single-scattering albedo and
  ! g factor in the longwave and shortwave
  subroutine cloud_optics(nlev,istartcol,iendcol, &
       &   config, thermodynamics, cloud, & 
       &   od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
       &   od_sw_cloud, ssa_sw_cloud, g_sw_cloud)

    use parkind1,                 only : jprb
    use radiation_config,         only : config_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_cloud,          only : cloud_type
    use radiation_constants,      only : AccelDueToGravity, &
         &   DensityLiquidWater, DensitySolidIce

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in) :: config
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(cloud_type),   intent(in) :: cloud

    ! Outputs

    ! Layer optical depth, single scattering albedo and g factor of
    ! clouds in each longwave band, where the latter two
    ! variables are only defined if cloud longwave scattering is
    ! enabled (otherwise both are treated as zero).
    real(jprb), dimension(config%n_bands_lw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_lw_cloud
    real(jprb), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol), &
         &   intent(out) :: ssa_lw_cloud, g_lw_cloud

    ! Layer optical depth, single scattering albedo and g factor of
    ! clouds in each shortwave band
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_sw_cloud, ssa_sw_cloud, g_sw_cloud

    ! In-cloud liquid and ice water path in a layer, in kg m-2
    real(jprb), dimension(nlev,istartcol:iendcol) :: lwp_kg_m2, iwp_kg_m2

    integer  :: jlev, jcol

    ! Factor to convert from gridbox-mean mass mixing ratio to
    ! in-cloud water path
    real(jprb) :: factor

    ! Convert cloud mixing ratio into liquid and ice water path
    ! in each layer
    do jlev = 1, nlev
      do jcol = istartcol, iendcol
        ! Factor to convert from gridbox-mean mass mixing ratio to
        ! in-cloud water path involves the pressure difference in
        ! Pa, acceleration due to gravity and cloud fraction
        ! adjusted to avoid division by zero.
        factor = ( thermodynamics%pressure_hl(jcol,jlev+1)    &
             &    -thermodynamics%pressure_hl(jcol,jlev  )  ) &
             &   / (AccelDueToGravity &
             &   * max(epsilon(1.0_jprb), cloud%fraction(jcol,jlev)))
        lwp_kg_m2(jlev,jcol) = factor * cloud%q_liq(jcol,jlev)
        iwp_kg_m2(jlev,jcol) = factor * cloud%q_ice(jcol,jlev)
      end do
    end do

    ! Geometric optics approximation: particles treated as much larger
    ! than the wavelength in both shortwave and longwave
    od_sw_cloud(1,:,:) &
         &   = (3.0_jprb/(2.0_jprb*DensityLiquidWater)) &
         &   * lwp_kg_m2 / transpose(cloud%re_liq(istartcol:iendcol,:)) &
         &   + (3.0_jprb / (2.0_jprb * DensitySolidIce)) &
         &   * iwp_kg_m2 / transpose(cloud%re_ice(istartcol:iendcol,:))
    od_lw_cloud(1,:,:) = lwp_kg_m2 * 137.22_jprb &
         &   + (3.0_jprb / (2.0_jprb * DensitySolidIce)) &
         &   * iwp_kg_m2 / transpose(cloud%re_ice(istartcol:iendcol,:))

    if (config%iverbose >= 4) then
      do jcol = istartcol,iendcol
        write(*,'(a,i0,a,f7.3,a,f7.3)') 'Profile ', jcol, ': shortwave optical depth = ', &
             &  sum(od_sw_cloud(1,:,jcol)*cloud%fraction(jcol,:)), &
             &  ', longwave optical depth = ', &
             &  sum(od_lw_cloud(1,:,jcol)*cloud%fraction(jcol,:))
        !    print *, 'LWP = ', sum(lwp_kg_m2(:,istartcol)*cloud%fraction(istartcol,:))
      end do
    end if

    ssa_sw_cloud = config%mono_sw_single_scattering_albedo
    g_sw_cloud   = config%mono_sw_asymmetry_factor

    ! In-place delta-Eddington scaling
    call delta_eddington(od_sw_cloud, ssa_sw_cloud, g_sw_cloud)

    if (config%do_lw_cloud_scattering) then
      ssa_lw_cloud = config%mono_lw_single_scattering_albedo
      g_lw_cloud   = config%mono_lw_asymmetry_factor
      ! In-place delta-Eddington scaling
      call delta_eddington(od_lw_cloud, ssa_lw_cloud, g_lw_cloud)
    end if

  end subroutine cloud_optics


  !---------------------------------------------------------------------
  ! Dummy subroutine since no aerosols are represented in
  ! monochromatic case
  subroutine add_aerosol_optics(nlev,istartcol,iendcol, &
       &  config, thermodynamics, gas, aerosol, & 
       &  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)

    use parkind1,                      only : jprb

    use radiation_config,              only : config_type
    use radiation_thermodynamics,      only : thermodynamics_type
    use radiation_gas,                 only : gas_type
    use radiation_aerosol,             only : aerosol_type

    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in), target :: config
    type(thermodynamics_type),intent(in)  :: thermodynamics
    type(gas_type),           intent(in)  :: gas
    type(aerosol_type),       intent(in)  :: aerosol
    ! Optical depth, single scattering albedo and asymmetry factor of
    ! the atmosphere (gases on input, gases and aerosols on output)
    ! for each g point. Note that longwave ssa and asymmetry and
    ! shortwave asymmetry are all zero for gases, so are not yet
    ! defined on input and are therefore intent(out).
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(inout) :: od_lw
    real(jprb), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol), &
         &  intent(out) :: ssa_lw, g_lw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(inout) &
         &  :: od_sw, ssa_sw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: g_sw

    g_sw(:,:,istartcol:iendcol) = 0.0_jprb

    if (config%do_lw_aerosol_scattering) then
      ssa_lw(:,:,istartcol:iendcol) = 0.0_jprb
      g_lw(:,:,istartcol:iendcol)   = 0.0_jprb
    end if

  end subroutine add_aerosol_optics

  !---------------------------------------------------------------------
  ! Planck function in terms of wavelength
  elemental function planck_function(wavelength, temperature)

    use parkind1,            only : jprb

    use radiation_constants, only : BoltzmannConstant, PlanckConstant, &
         &                          SpeedOfLight

    real(jprb), intent(in) :: wavelength  ! metres
    real(jprb), intent(in) :: temperature ! Kelvin

    ! Output in W sr-1 m-3
    real(jprb)             :: planck_function

    if (temperature > 0.0_jprb) then
      planck_function = 2.0_jprb * PlanckConstant * SpeedOfLight**2 &
           &   / (wavelength**5 &
           &   * (exp(PlanckConstant*SpeedOfLight &
           &         /(wavelength*BoltzmannConstant*temperature)) - 1.0_jprb))
    else
      planck_function = 0.0_jprb
    end if

  end function planck_function

end module radiation_monochromatic
