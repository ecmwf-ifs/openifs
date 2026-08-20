! radiation_aerosol_optics.F90 - Computing aerosol optical properties
!
! (C) Copyright 2015- ECMWF.
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
!   2018-04-15  R. Hogan  Add "direct" option
!   2020-11-14  R. Hogan  Add setup_general_aerosol_optics for ecCKD compatibility
!   2022-03-27  R. Hogan  Add setup_general_aerosol_optics_legacy to use RRTM aerosol files with ecCKD
!   2022-11-22  P. Ukkonen / R. Hogan  Optimizations to enhance vectorization

module radiation_aerosol_optics

  implicit none
  public

contains

  ! Provides the elemental function "delta_eddington_extensive"
#include "radiation_delta_eddington.h"

  !---------------------------------------------------------------------
  ! Load aerosol scattering data; this subroutine delegates to one
  ! in radiation_aerosol_optics_data.F90
  subroutine setup_aerosol_optics(config)

    use parkind1,                      only : jprb
    use yomhook,                       only : lhook, dr_hook, jphook
    use radiation_config,              only : config_type
    use radiation_aerosol_optics_data, only : aerosol_optics_type
    use radiation_io,                  only : nulerr, radiation_abort

    type(config_type), intent(inout) :: config

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics:setup_aerosol_optics',0,hook_handle)

    if (config%n_aerosol_types > 0) then
      ! Load data from file and prepare to map config%n_aerosol_types
      ! aerosol types
      if (config%use_general_aerosol_optics) then
        ! Read file containing high spectral resolution optical
        ! properties and average to the spectral intervals of the
        ! current gas-optics scheme
        call setup_general_aerosol_optics(config)
      else
        ! Read file containing optical properties already in the bands
        ! of the gas-optics scheme
        call config%aerosol_optics%setup(trim(config%aerosol_optics_file_name), &
             &                           iverbose=config%iverbosesetup)
      end if

      call config%aerosol_optics%initialize_types(config%n_aerosol_types)

      ! Check agreement in number of bands
      if (config%n_bands_lw /= config%aerosol_optics%n_bands_lw) then
        write(nulerr,'(a,i0,a,i0,a)') '*** Error: number of longwave bands (', &
             &  config%n_bands_lw, ') does not match aerosol optics look-up table (', &
             &  config%aerosol_optics%n_bands_lw, ')'
        call radiation_abort()
      end if
      if (config%n_bands_sw /= config%aerosol_optics%n_bands_sw) then
        write(nulerr,'(a)') '*** Error: number of shortwave bands does not match aerosol optics look-up table'
        call radiation_abort()
      end if

      ! Map aerosol types to those loaded from the data file
      call config%aerosol_optics%set_types(config%i_aerosol_type_map(1:config%n_aerosol_types))
    end if

    if (config%iverbosesetup >= 1) then
      call config%aerosol_optics%print_description(config%i_aerosol_type_map(1:config%n_aerosol_types))
    end if

    if (lhook) call dr_hook('radiation_aerosol_optics:setup_aerosol_optics',1,hook_handle)

  end subroutine setup_aerosol_optics


  !---------------------------------------------------------------------
  ! Read file containing high spectral resolution optical properties
  ! and average to the spectral intervals of the current gas-optics
  ! scheme
  subroutine setup_general_aerosol_optics(config)

    use parkind1,                      only : jprb
    use yomhook,                       only : lhook, dr_hook, jphook
    use easy_netcdf_read_mpi,          only : netcdf_file
    use radiation_config,              only : config_type
    use radiation_aerosol_optics_data, only : aerosol_optics_type
    use radiation_spectral_definition, only : SolarReferenceTemperature, &
         &                                    TerrestrialReferenceTemperature
    use radiation_io,                  only : nulout

    type(config_type), intent(inout), target :: config

    ! The NetCDF file containing the aerosol optics data
    type(netcdf_file)  :: file

    ! Wavenumber points in NetCDF file
    real(jprb), allocatable :: wavenumber(:) ! cm-1

    ! Hydrophilic aerosol properties
    real(jprb), allocatable :: mass_ext_philic(:,:,:)    ! Mass-ext coefficient (m2 kg-1)
    real(jprb), allocatable :: ssa_philic(:,:,:)         ! Single-scattering albedo
    real(jprb), allocatable :: g_philic(:,:,:)           ! Asymmetry factor
    real(jprb), allocatable :: lidar_ratio_philic(:,:,:) ! Lidar ratio (sr)

    ! Hydrophobic aerosol properties
    real(jprb), allocatable :: mass_ext_phobic(:,:)      ! Mass-ext coefficient (m2 kg-1)
    real(jprb), allocatable :: ssa_phobic(:,:)           ! Single-scattering albedo
    real(jprb), allocatable :: g_phobic(:,:)             ! Asymmetry factor
    real(jprb), allocatable :: lidar_ratio_phobic(:,:)   ! Lidar ratio (sr)

    ! Mapping matrix between optical properties at the wavenumbers in
    ! the file, and spectral intervals used by the gas-optics scheme
    real(jprb), allocatable :: mapping(:,:)

    ! Pointer to the aerosol optics coefficients for brevity of access
    type(aerosol_optics_type), pointer :: ao

    ! Target monochromatic wavenumber for interpolation (cm-1)
    real(jprb) :: wavenumber_target

    ! Number of spectral points describing aerosol properties in the
    ! shortwave and longwave
    integer    :: nspecsw, nspeclw

    ! Number of monochromatic wavelengths required
    integer    :: nmono

    integer    :: n_type_philic, n_type_phobic, nrh, nwn
    integer    :: jtype, jwl, iwn

    ! Weight of first point in interpolation
    real(jprb) :: weight1

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics:setup_general_aerosol_optics',0,hook_handle)

    ao => config%aerosol_optics

    call file%open(trim(config%aerosol_optics_file_name), iverbose=config%iverbosesetup)

    if (.not. file%exists('wavenumber')) then
      ! Assume we have an old-style aerosol optics file with optical
      ! properties provided per pre-defined band
      call file%close()
      if (config%iverbosesetup >= 2) then
        write(nulout,'(a)') 'Legacy aerosol optics file: mapping between bands'
      end if
      call setup_general_aerosol_optics_legacy(config, trim(config%aerosol_optics_file_name))
      if (lhook) call dr_hook('radiation_aerosol_optics:setup_general_aerosol_optics',1,hook_handle)
      return
    end if

    if (file%exists('mass_ext_hydrophilic')) then
      ao%use_hydrophilic = .true.
    else
      ao%use_hydrophilic = .false.
    end if
 
    call file%get('wavenumber', wavenumber)
    nwn = size(wavenumber)

    ! Read the raw scattering data
    call file%get('mass_ext_hydrophobic',    mass_ext_phobic)
    call file%get('ssa_hydrophobic',         ssa_phobic)
    call file%get('asymmetry_hydrophobic',   g_phobic)
    call file%get('lidar_ratio_hydrophobic', lidar_ratio_phobic)

    call file%get_global_attribute('description_hydrophobic', &
         &                         ao%description_phobic_str)

    if (ao%use_hydrophilic) then
      call file%get('mass_ext_hydrophilic',    mass_ext_philic)
      call file%get('ssa_hydrophilic',         ssa_philic)
      call file%get('asymmetry_hydrophilic',   g_philic)
      call file%get('lidar_ratio_hydrophilic', lidar_ratio_philic)

      call file%get('relative_humidity1',      ao%rh_lower)

      call file%get_global_attribute('description_hydrophilic', &
           &                         ao%description_philic_str)
    end if

    ! Close aerosol scattering file
    call file%close()

    n_type_phobic = size(mass_ext_phobic, 2)
    if (ao%use_hydrophilic) then
      n_type_philic = size(mass_ext_philic, 3)
      nrh = size(ao%rh_lower)
    else
      n_type_philic = 0
      nrh = 0
    end if

    if (config%do_cloud_aerosol_per_sw_g_point) then
      nspecsw = config%gas_optics_sw%spectral_def%ng
    else
      nspecsw = config%gas_optics_sw%spectral_def%nband
    end if

    if (config%do_cloud_aerosol_per_lw_g_point) then
      nspeclw = config%gas_optics_lw%spectral_def%ng
    else
      nspeclw = config%gas_optics_lw%spectral_def%nband
    end if

    if (allocated(ao%wavelength_mono)) then
      ! Monochromatic wavelengths also required
      nmono = size(ao%wavelength_mono)
    else
      nmono = 0
    end if

    call ao%allocate(n_type_phobic, n_type_philic, nrh, nspeclw, nspecsw, nmono)

    if (config%do_sw) then
      call config%gas_optics_sw%spectral_def%calc_mapping(wavenumber, mapping, &
           &          use_bands=(.not. config%do_cloud_aerosol_per_sw_g_point))

      ao%mass_ext_sw_phobic = matmul(mapping, mass_ext_phobic)
      ao%ssa_sw_phobic = matmul(mapping, mass_ext_phobic*ssa_phobic) &
           &           / ao%mass_ext_sw_phobic
      ao%g_sw_phobic = matmul(mapping, mass_ext_phobic*ssa_phobic*g_phobic) &
           &         / (ao%mass_ext_sw_phobic*ao%ssa_sw_phobic)

      if (ao%use_hydrophilic) then
        do jtype = 1,n_type_philic
          ao%mass_ext_sw_philic(:,:,jtype) = matmul(mapping, mass_ext_philic(:,:,jtype))
          ao%ssa_sw_philic(:,:,jtype) = matmul(mapping, mass_ext_philic(:,:,jtype) &
               &                                        *ssa_philic(:,:,jtype)) &
               &           / ao%mass_ext_sw_philic(:,:,jtype)
          ao%g_sw_philic(:,:,jtype) = matmul(mapping, mass_ext_philic(:,:,jtype) &
               &                       *ssa_philic(:,:,jtype)*g_philic(:,:,jtype)) &
               &         / (ao%mass_ext_sw_philic(:,:,jtype)*ao%ssa_sw_philic(:,:,jtype))
        end do
      end if
    end if

    if (config%do_lw) then
      call config%gas_optics_lw%spectral_def%calc_mapping(wavenumber, mapping, &
           &          use_bands=(.not. config%do_cloud_aerosol_per_lw_g_point))

      ao%mass_ext_lw_phobic = matmul(mapping, mass_ext_phobic)
      ao%ssa_lw_phobic = matmul(mapping, mass_ext_phobic*ssa_phobic) &
           &           / ao%mass_ext_lw_phobic
      ao%g_lw_phobic = matmul(mapping, mass_ext_phobic*ssa_phobic*g_phobic) &
           &         / (ao%mass_ext_lw_phobic*ao%ssa_lw_phobic)

      if (ao%use_hydrophilic) then
        do jtype = 1,n_type_philic
          ao%mass_ext_lw_philic(:,:,jtype) = matmul(mapping, mass_ext_philic(:,:,jtype))
          ao%ssa_lw_philic(:,:,jtype) = matmul(mapping, mass_ext_philic(:,:,jtype) &
               &                                        *ssa_philic(:,:,jtype)) &
               &           / ao%mass_ext_lw_philic(:,:,jtype)
          ao%g_lw_philic(:,:,jtype) = matmul(mapping, mass_ext_philic(:,:,jtype) &
               &                       *ssa_philic(:,:,jtype)*g_philic(:,:,jtype)) &
               &         / (ao%mass_ext_lw_philic(:,:,jtype)*ao%ssa_lw_philic(:,:,jtype))
        end do
      end if
    end if

    if (allocated(ao%wavelength_mono)) then
      ! Monochromatic wavelengths also required
      do jwl = 1,nmono
        ! Wavelength (m) to wavenumber (cm-1)
        wavenumber_target = 0.01_jprb / ao%wavelength_mono(jwl)
        ! Find index to first interpolation point, and its weight
        if (wavenumber_target <= wavenumber(1)) then
          weight1 = 1.0_jprb
          iwn = 1
        else if (wavenumber_target >= wavenumber(nwn)) then
          iwn = nwn-1
          weight1 = 0.0_jprb
        else
          iwn = 1
          do while (wavenumber(iwn+1) < wavenumber_target .and. iwn < nwn-1)
            iwn = iwn + 1
          end do
          weight1 = (wavenumber(iwn+1)-wavenumber_target) &
               &  / (wavenumber(iwn+1)-wavenumber(iwn))
        end if
        ! Linear interpolation
        ao%mass_ext_mono_phobic(jwl,:) = weight1 * mass_ext_phobic(iwn,:) &
             &             + (1.0_jprb - weight1)* mass_ext_phobic(iwn+1,:)
        ao%ssa_mono_phobic(jwl,:)      = weight1 * ssa_phobic(iwn,:) &
             &             + (1.0_jprb - weight1)* ssa_phobic(iwn+1,:)
        ao%g_mono_phobic(jwl,:)        = weight1 * g_phobic(iwn,:) &
             &             + (1.0_jprb - weight1)* g_phobic(iwn+1,:)
        ao%lidar_ratio_mono_phobic(jwl,:) = weight1 * lidar_ratio_phobic(iwn,:) &
             &                + (1.0_jprb - weight1)* lidar_ratio_phobic(iwn+1,:)
        if (ao%use_hydrophilic) then
          ao%mass_ext_mono_philic(jwl,:,:) = weight1 * mass_ext_philic(iwn,:,:) &
               &               + (1.0_jprb - weight1)* mass_ext_philic(iwn+1,:,:)
          ao%ssa_mono_philic(jwl,:,:)      = weight1 * ssa_philic(iwn,:,:) &
               &               + (1.0_jprb - weight1)* ssa_philic(iwn+1,:,:)
          ao%g_mono_philic(jwl,:,:)        = weight1 * g_philic(iwn,:,:) &
               &               + (1.0_jprb - weight1)* g_philic(iwn+1,:,:)
          ao%lidar_ratio_mono_philic(jwl,:,:) = weight1 * lidar_ratio_philic(iwn,:,:) &
               &                  + (1.0_jprb - weight1)* lidar_ratio_philic(iwn+1,:,:)
        end if
      end do
    end if

    ! Deallocate memory local to this routine
    deallocate(mass_ext_phobic)
    deallocate(ssa_phobic)
    deallocate(g_phobic)
    deallocate(lidar_ratio_phobic)
    if (ao%use_hydrophilic) then
      deallocate(mass_ext_philic)
      deallocate(ssa_philic)
      deallocate(g_philic)
      deallocate(lidar_ratio_philic)
    end if

    if (lhook) call dr_hook('radiation_aerosol_optics:setup_general_aerosol_optics',1,hook_handle)

  end subroutine setup_general_aerosol_optics


  !---------------------------------------------------------------------
  ! Read file containing legacy-style band-wise aerosol optical
  ! properties and average to the spectral intervals of the current
  ! gas-optics scheme
  subroutine setup_general_aerosol_optics_legacy(config, file_name)

    use parkind1,                      only : jprb
    use yomhook,                       only : lhook, dr_hook, jphook
    use easy_netcdf_read_mpi,          only : netcdf_file
    use radiation_config,              only : config_type
    use radiation_aerosol_optics_data, only : aerosol_optics_type
    use radiation_spectral_definition, only : SolarReferenceTemperature, &
         &                                    TerrestrialReferenceTemperature

    type(config_type), intent(inout), target :: config

    ! The NetCDF file containing the aerosol optics data
    character(len=*), intent(in) :: file_name

    ! Mapping matrix between optical properties at the wavenumbers in
    ! the file, and spectral intervals used by the gas-optics scheme
    real(jprb), allocatable :: mapping(:,:), mapping_transp(:,:)

    ! Pointer to the aerosol optics coefficients for brevity of access
    type(aerosol_optics_type), pointer :: ao

    ! Local copy of aerosol optical properties in the spectral
    ! intervals of the file, which is deallocated when it goes out of
    ! scope
    type(aerosol_optics_type) :: ao_legacy

    integer :: jtype

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics:setup_general_aerosol_optics_legacy',0,hook_handle)
    ao => config%aerosol_optics

    ! Load file into a local structure
    call ao_legacy%setup(file_name, iverbose=config%iverbosesetup)

    ! Copy over scalars and coordinate variables
    call ao%allocate(ao_legacy%n_type_phobic, ao_legacy%n_type_philic, ao_legacy%nrh, &
         &           config%n_bands_lw, config%n_bands_sw, ao_legacy%n_mono_wl)
    ao%description_phobic_str = ao_legacy%description_phobic_str
    if (ao_legacy%use_hydrophilic) then
      ao%description_philic_str = ao_legacy%description_philic_str
      ao%rh_lower = ao_legacy%rh_lower
    end if

    ! use_hydrophilic = ao_legacy%use_hydrophilic
    ! ao%iclass = ao_legacy%iclass
    ! ao%itype = ao_legacy%itype
    ! ao%ntype = ao_legacy%ntype
    ! ao%n_type_phobic = ao_legacy%n_type_phobic
    ! ao%n_type_philic = ao_legacy%n_type_philic
    ! ao%n_mono_wl = ao_legacy%n_mono_wl
    ! ao%use_monochromatic = ao_legacy%use_monochromatic

    if (config%do_sw) then
      call config%gas_optics_sw%spectral_def%calc_mapping_from_wavenumber_bands( &
           &  ao_legacy%wavenumber1_sw, ao_legacy%wavenumber2_sw, mapping_transp, &
           &  use_bands=(.not. config%do_cloud_aerosol_per_sw_g_point))
      if (allocated(mapping)) then
        deallocate(mapping)
      end if
      allocate(mapping(config%n_bands_sw,ao_legacy%n_bands_sw))
      mapping = transpose(mapping_transp)
      ao%mass_ext_sw_phobic = matmul(mapping, ao_legacy%mass_ext_sw_phobic)
      ao%ssa_sw_phobic = matmul(mapping, ao_legacy%mass_ext_sw_phobic*ao_legacy%ssa_sw_phobic) &
           &           / ao%mass_ext_sw_phobic
      ao%g_sw_phobic = matmul(mapping, ao_legacy%mass_ext_sw_phobic*ao_legacy%ssa_sw_phobic &
           &                           *ao_legacy%g_sw_phobic) &
           &         / (ao%mass_ext_sw_phobic*ao%ssa_sw_phobic)

      if (ao%use_hydrophilic) then
        do jtype = 1,ao%n_type_philic
          ao%mass_ext_sw_philic(:,:,jtype) = matmul(mapping, ao_legacy%mass_ext_sw_philic(:,:,jtype))
          ao%ssa_sw_philic(:,:,jtype) = matmul(mapping, ao_legacy%mass_ext_sw_philic(:,:,jtype) &
               &                                        *ao_legacy%ssa_sw_philic(:,:,jtype)) &
               &           / ao%mass_ext_sw_philic(:,:,jtype)
          ao%g_sw_philic(:,:,jtype) = matmul(mapping, ao_legacy%mass_ext_sw_philic(:,:,jtype) &
               &               *ao_legacy%ssa_sw_philic(:,:,jtype)*ao_legacy%g_sw_philic(:,:,jtype)) &
               &         / (ao%mass_ext_sw_philic(:,:,jtype)*ao%ssa_sw_philic(:,:,jtype))
        end do
      end if
    end if

    if (config%do_lw) then
      if (allocated(mapping_transp)) then
        deallocate(mapping_transp)
      end if
      call config%gas_optics_lw%spectral_def%calc_mapping_from_wavenumber_bands( &
           &  ao_legacy%wavenumber1_lw, ao_legacy%wavenumber2_lw, mapping_transp, &
           &  use_bands=(.not. config%do_cloud_aerosol_per_lw_g_point))
      if (allocated(mapping)) then
        deallocate(mapping)
      end if
      allocate(mapping(config%n_bands_lw,ao_legacy%n_bands_lw))
      mapping = transpose(mapping_transp)
      ao%mass_ext_lw_phobic = matmul(mapping, ao_legacy%mass_ext_lw_phobic)
      ao%ssa_lw_phobic = matmul(mapping, ao_legacy%mass_ext_lw_phobic*ao_legacy%ssa_lw_phobic) &
           &           / ao%mass_ext_lw_phobic
      ao%g_lw_phobic = matmul(mapping, ao_legacy%mass_ext_lw_phobic*ao_legacy%ssa_lw_phobic &
           &                           *ao_legacy%g_lw_phobic) &
           &         / (ao%mass_ext_lw_phobic*ao%ssa_lw_phobic)

      if (ao%use_hydrophilic) then
        do jtype = 1,ao%n_type_philic
          ao%mass_ext_lw_philic(:,:,jtype) = matmul(mapping, ao_legacy%mass_ext_lw_philic(:,:,jtype))
          ao%ssa_lw_philic(:,:,jtype) = matmul(mapping, ao_legacy%mass_ext_lw_philic(:,:,jtype) &
               &                                        *ao_legacy%ssa_lw_philic(:,:,jtype)) &
               &           / ao%mass_ext_lw_philic(:,:,jtype)
          ao%g_lw_philic(:,:,jtype) = matmul(mapping, ao_legacy%mass_ext_lw_philic(:,:,jtype) &
               &               *ao_legacy%ssa_lw_philic(:,:,jtype)*ao_legacy%g_lw_philic(:,:,jtype)) &
               &         / (ao%mass_ext_lw_philic(:,:,jtype)*ao%ssa_lw_philic(:,:,jtype))
        end do
      end if
    end if

    if (allocated(ao_legacy%wavelength_mono)) then
      ao%wavelength_mono = ao_legacy%wavelength_mono
      ao%mass_ext_mono_phobic = ao_legacy%mass_ext_mono_phobic
      ao%ssa_mono_phobic = ao_legacy%ssa_mono_phobic
      ao%g_mono_phobic = ao_legacy%g_mono_phobic
      ao%lidar_ratio_mono_phobic = ao_legacy%lidar_ratio_mono_phobic
      if (ao%use_hydrophilic) then
        ao%mass_ext_mono_philic = ao_legacy%mass_ext_mono_philic
        ao%ssa_mono_philic = ao_legacy%ssa_mono_philic
        ao%g_mono_philic = ao_legacy%g_mono_philic
        ao%lidar_ratio_mono_philic = ao_legacy%lidar_ratio_mono_philic
      end if
    end if

    if (lhook) call dr_hook('radiation_aerosol_optics:setup_general_aerosol_optics_legacy',1,hook_handle)

  end subroutine setup_general_aerosol_optics_legacy


  !---------------------------------------------------------------------
  ! Compute aerosol optical properties and add to existing gas optical
  ! depth and scattering properties
  subroutine add_aerosol_optics(nlev,istartcol,iendcol, &
       &  config, thermodynamics, gas, aerosol, & 
       &  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)

    use parkind1,                      only : jprb
    use radiation_io,                  only : nulout, nulerr, radiation_abort
    use yomhook,                       only : lhook, dr_hook, jphook
    use radiation_config,              only : config_type
    use radiation_thermodynamics,      only : thermodynamics_type
    use radiation_gas,                 only : gas_type, IH2O, IMassMixingRatio
    use radiation_aerosol,             only : aerosol_type
    use radiation_constants,           only : AccelDueToGravity
    use radiation_aerosol_optics_data, only : aerosol_optics_type, &
         &  IAerosolClassUndefined,   IAerosolClassIgnored, &
         &  IAerosolClassHydrophobic, IAerosolClassHydrophilic

    real(jprb), parameter :: OneOverAccelDueToGravity = 1.0_jprb / AccelDueToGravity
    
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
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), &
         &   intent(inout) :: od_lw
    real(jprb), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol), &
         &   intent(out)   :: ssa_lw, g_lw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), &
         &   intent(inout) :: od_sw, ssa_sw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), &
         &   intent(out)   :: g_sw

    ! Extinction optical depth, scattering optical depth and
    ! asymmetry-times-scattering-optical-depth for all the aerosols in
    ! a column for each spectral band of the shortwave and longwave
    ! spectrum
    real(jprb), dimension(config%n_bands_sw,nlev) &
         &  :: od_sw_aerosol, scat_sw_aerosol, scat_g_sw_aerosol
    real(jprb), dimension(config%n_bands_lw,nlev) &
         &  :: od_lw_aerosol
    real(jprb), dimension(config%n_bands_lw_if_scattering,nlev) &
         &  :: scat_lw_aerosol, scat_g_lw_aerosol

    real(jprb) :: local_od_sw, local_od_lw
    
    real(jprb) :: h2o_mmr(istartcol:iendcol,nlev)

    real(jprb) :: rh ! Relative humidity with respect to liquid water

    ! Factor (kg m-2) to convert mixing ratio (kg kg-1) to mass in
    ! path (kg m-2)
    real(jprb) :: factor(nlev)

    ! Temporary extinction and scattering optical depths of aerosol
    ! plus gas
    real(jprb) :: local_od, local_scat

    ! Aerosol mixing ratio as a scalar
    real(jprb) :: mixing_ratio
    
    ! Loop indices for column, level, g point, band and aerosol type
    integer :: jcol, jlev, jg, jtype, jband

    ! Range of levels over which aerosols are present
    integer :: istartlev, iendlev

    ! Indices to spectral band and relative humidity look-up table
    integer :: iband, irh, irhs(nlev)

    ! Short cut for ao%itype(jtype)
    integer :: itype

    ! Pointer to the aerosol optics coefficients for brevity of access
    type(aerosol_optics_type), pointer :: ao

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics:add_aerosol_optics',0,hook_handle)

    if (aerosol%is_direct) then
      ! Aerosol optical properties have been provided in each band
      ! directly by the user
      call add_aerosol_optics_direct(nlev,istartcol,iendcol, &
           &  config, aerosol, & 
           &  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)
    else
      ! Aerosol mixing ratios have been provided

      do jtype = 1,config%n_aerosol_types
        if (config%aerosol_optics%iclass(jtype) == IAerosolClassUndefined) then
          write(nulerr,'(a)') '*** Error: not all aerosol types are defined'
          call radiation_abort()
        end if
      end do

      if (config%iverbose >= 2) then
        write(nulout,'(a)') 'Computing aerosol absorption/scattering properties'
      end if

      ao => config%aerosol_optics

      istartlev = lbound(aerosol%mixing_ratio,2)
      iendlev   = ubound(aerosol%mixing_ratio,2)

      if (ubound(aerosol%mixing_ratio,3) /= config%n_aerosol_types) then
        write(nulerr,'(a,i0,a,i0)') '*** Error: aerosol%mixing_ratio contains ', &
             &  ubound(aerosol%mixing_ratio,3), ' aerosol types, expected ', &
             &  config%n_aerosol_types
        call radiation_abort()
      end if

      ! Set variables to zero that may not have been previously
      g_sw(:,:,istartcol:iendcol) = 0.0_jprb
      if (config%do_lw_aerosol_scattering) then
        ssa_lw(:,:,istartcol:iendcol) = 0.0_jprb
        g_lw(:,:,istartcol:iendcol)   = 0.0_jprb
      end if

      call gas%get(IH2O, IMassMixingRatio, h2o_mmr, istartcol=istartcol)

      ! Loop over column
      do jcol = istartcol,iendcol

        ! Reset temporary arrays
        od_sw_aerosol     = 0.0_jprb
        scat_sw_aerosol   = 0.0_jprb
        scat_g_sw_aerosol = 0.0_jprb
        od_lw_aerosol     = 0.0_jprb
        scat_lw_aerosol   = 0.0_jprb
        scat_g_lw_aerosol = 0.0_jprb

        do jlev = istartlev,iendlev
          ! Compute relative humidity with respect to liquid
          ! saturation and the index to the relative-humidity index of
          ! hydrophilic-aerosol data
          rh  = h2o_mmr(jcol,jlev) / thermodynamics%h2o_sat_liq(jcol,jlev)
          irhs(jlev) = ao%calc_rh_index(rh)
          
          factor(jlev) = ( thermodynamics%pressure_hl(jcol,jlev+1) &
               &    -thermodynamics%pressure_hl(jcol,jlev  )  ) &
               &   * OneOverAccelDueToGravity  
        end do
        
        do jtype = 1,config%n_aerosol_types
          itype = ao%itype(jtype)

          ! Add the optical depth, scattering optical depth and
          ! scattering optical depth-weighted asymmetry factor for
          ! this aerosol type to the total for all aerosols.  Note
          ! that the following expressions are array-wise, the
          ! dimension being spectral band.
          if (ao%iclass(jtype) == IAerosolClassHydrophobic) then
            do jlev = istartlev,iendlev
              mixing_ratio = aerosol%mixing_ratio(jcol,jlev,jtype)
              do jband = 1,config%n_bands_sw
                local_od_sw = factor(jlev) * mixing_ratio &
                     &  * ao%mass_ext_sw_phobic(jband,itype)
                od_sw_aerosol(jband,jlev) = od_sw_aerosol(jband,jlev) + local_od_sw
                scat_sw_aerosol(jband,jlev) = scat_sw_aerosol(jband,jlev) &
                     &  + local_od_sw * ao%ssa_sw_phobic(jband,itype)
                scat_g_sw_aerosol(jband,jlev) = scat_g_sw_aerosol(jband,jlev) &
                     &  + local_od_sw * ao%ssa_sw_phobic(jband,itype) &
                     &  * ao%g_sw_phobic(jband,itype)
              end do
              if (config%do_lw_aerosol_scattering) then
                do jband = 1,config%n_bands_lw
                  local_od_lw = factor(jlev) * mixing_ratio &
                       &  * ao%mass_ext_lw_phobic(jband,itype)
                  od_lw_aerosol(jband,jlev) = od_lw_aerosol(jband,jlev) + local_od_lw
                  scat_lw_aerosol(jband,jlev) = scat_lw_aerosol(jband,jlev) &
                       &  + local_od_lw * ao%ssa_lw_phobic(jband,itype)
                  scat_g_lw_aerosol(jband,jlev) = scat_g_lw_aerosol(jband,jlev) &
                       &  + local_od_lw * ao%ssa_lw_phobic(jband,itype) &
                       &  * ao%g_lw_phobic(jband,itype)
                end do
              else
                ! If aerosol longwave scattering is not included then we
                ! weight the optical depth by the single scattering
                ! co-albedo
                do jband = 1,config%n_bands_lw
                  od_lw_aerosol(jband,jlev) = od_lw_aerosol(jband,jlev) &
                       &  + factor(jlev) * mixing_ratio &
                       &  * ao%mass_ext_lw_phobic(jband,itype) &
                       &  * (1.0_jprb - ao%ssa_lw_phobic(jband,itype))
                end do
              end if
            end do
            
          else if (ao%iclass(jtype) == IAerosolClassHydrophilic) then
            ! Hydrophilic aerosols require the look-up tables to
            ! be indexed with irh
            do jlev = istartlev,iendlev  
              mixing_ratio = aerosol%mixing_ratio(jcol,jlev,jtype)
              irh = irhs(jlev)
              do jband = 1,config%n_bands_sw
                local_od_sw = factor(jlev) * mixing_ratio &
                     &  * ao%mass_ext_sw_philic(jband,irh,itype)
                od_sw_aerosol(jband,jlev) = od_sw_aerosol(jband,jlev) + local_od_sw
                scat_sw_aerosol(jband,jlev) = scat_sw_aerosol(jband,jlev) &
                     &  + local_od_sw * ao%ssa_sw_philic(jband,irh,itype)
                scat_g_sw_aerosol(jband,jlev) = scat_g_sw_aerosol(jband,jlev) &
                     &  + local_od_sw * ao%ssa_sw_philic(jband,irh,itype) &
                     &  * ao%g_sw_philic(jband,irh,itype)
              end do
              if (config%do_lw_aerosol_scattering) then
                do jband = 1,config%n_bands_lw
                  local_od_lw = factor(jlev) * mixing_ratio &
                       &  * ao%mass_ext_lw_philic(jband,irh,itype)
                  od_lw_aerosol(jband,jlev) = od_lw_aerosol(jband,jlev) + local_od_lw
                  scat_lw_aerosol(jband,jlev) = scat_lw_aerosol(jband,jlev) &
                       &  + local_od_lw * ao%ssa_lw_philic(jband,irh,itype)
                  scat_g_lw_aerosol(jband,jlev) = scat_g_lw_aerosol(jband,jlev) &
                       &  + local_od_lw * ao%ssa_lw_philic(jband,irh,itype) &
                       &  * ao%g_lw_philic(jband,irh,itype)
                end do
              else
                ! If aerosol longwave scattering is not included then we
                ! weight the optical depth by the single scattering
                ! co-albedo
                do jband = 1,config%n_bands_lw
                  od_lw_aerosol(jband,jlev) = od_lw_aerosol(jband,jlev) &
                       &  + factor(jlev) * mixing_ratio &
                       &  * ao%mass_ext_lw_philic(jband,irh,itype) &
                       &  * (1.0_jprb - ao%ssa_lw_philic(jband,irh,itype))
                end do
              end if
            end do
            
            ! Implicitly, if ao%iclass(jtype) == IAerosolClassNone, then
            ! no aerosol scattering properties are added
          end if

        end do ! Loop over aerosol type
        
        if (.not. config%do_sw_delta_scaling_with_gases) then
          ! Delta-Eddington scaling on aerosol only.  Note that if
          ! do_sw_delta_scaling_with_gases==.true. then the delta
          ! scaling is done to the cloud-aerosol-gas mixture inside
          ! the solver
          call delta_eddington_extensive_vec(config%n_bands_sw*nlev, od_sw_aerosol, &
               &                             scat_sw_aerosol, scat_g_sw_aerosol)
        end if

        ! Combine aerosol shortwave scattering properties with gas
        ! properties (noting that any gas scattering will have an
        ! asymmetry factor of zero)
        if (config%do_cloud_aerosol_per_sw_g_point) then

          ! We can assume the band and g-point indices are the same
          do jlev = 1,nlev
            do jg = 1,config%n_g_sw
              local_scat = ssa_sw(jg,jlev,jcol)*od_sw(jg,jlev,jcol) + scat_sw_aerosol(jg,jlev)
              od_sw(jg,jlev,jcol) = od_sw(jg,jlev,jcol) + od_sw_aerosol(jg,jlev)
              g_sw(jg,jlev,jcol) = scat_g_sw_aerosol(jg,jlev) / max(local_scat, 1.0e-24_jprb)
              ssa_sw(jg,jlev,jcol) = min(local_scat / max(od_sw(jg,jlev,jcol), 1.0e-24_jprb), 1.0_jprb)
            end do
          end do

        else

          do jlev = 1,nlev
            do jg = 1,config%n_g_sw
              ! Need to map between bands and g-points
              iband = config%i_band_from_reordered_g_sw(jg)
              local_od = od_sw(jg,jlev,jcol) + od_sw_aerosol(iband,jlev)
              if (local_od > 0.0_jprb .and. od_sw_aerosol(iband,jlev) > 0.0_jprb) then
                local_scat = ssa_sw(jg,jlev,jcol) * od_sw(jg,jlev,jcol) &
                     &  + scat_sw_aerosol(iband,jlev)
                ! Note that asymmetry_sw of gases is zero so the following
                ! simply weights the aerosol asymmetry by the scattering
                ! optical depth
                if (local_scat > 0.0_jprb) then
                  g_sw(jg,jlev,jcol) = scat_g_sw_aerosol(iband,jlev) / local_scat
                end if
                ssa_sw(jg,jlev,jcol) = local_scat / local_od
                od_sw (jg,jlev,jcol) = local_od
              end if
            end do
          end do

        end if

        ! Combine aerosol longwave scattering properties with gas
        ! properties, noting that in the longwave, gases do not
        ! scatter at all
        if (config%do_lw_aerosol_scattering) then

          call delta_eddington_extensive_vec(config%n_bands_lw*nlev, od_lw_aerosol, &
               &                             scat_lw_aerosol, scat_g_lw_aerosol)  

          do jlev = istartlev,iendlev
            do jg = 1,config%n_g_lw
              iband = config%i_band_from_reordered_g_lw(jg)
              local_od = od_lw(jg,jlev,jcol) + od_lw_aerosol(iband,jlev)
              if (local_od > 0.0_jprb .and. od_lw_aerosol(iband,jlev) > 0.0_jprb) then
                ! All scattering is due to aerosols, therefore the
                ! asymmetry factor is equal to the value for aerosols
                if (scat_lw_aerosol(iband,jlev) > 0.0_jprb) then
                  g_lw(jg,jlev,jcol) = scat_g_lw_aerosol(iband,jlev) &
                       &  / scat_lw_aerosol(iband,jlev)
                end if
                ssa_lw(jg,jlev,jcol) = scat_lw_aerosol(iband,jlev) / local_od
                od_lw (jg,jlev,jcol) = local_od
              end if
            end do
          end do
          
        else

          if (config%do_cloud_aerosol_per_lw_g_point) then
            ! We can assume band and g-point indices are the same
            do jlev = istartlev,iendlev
              do jg = 1,config%n_g_lw
                od_lw(jg,jlev,jcol) = od_lw(jg,jlev,jcol) + od_lw_aerosol(jg,jlev)
              end do
            end do
          else
            do jlev = istartlev,iendlev
              do jg = 1,config%n_g_lw
                od_lw(jg,jlev,jcol) = od_lw(jg,jlev,jcol) &
                     &  + od_lw_aerosol(config%i_band_from_reordered_g_lw(jg),jlev)
              end do
            end do
          end if

        end if
        
      end do ! Loop over column

    end if

    if (lhook) call dr_hook('radiation_aerosol_optics:add_aerosol_optics',1,hook_handle)

  end subroutine add_aerosol_optics


  !---------------------------------------------------------------------
  ! Add precomputed optical properties to gas optical depth and
  ! scattering properties
  subroutine add_aerosol_optics_direct(nlev,istartcol,iendcol, &
       &  config, aerosol, & 
       &  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)

    use parkind1,                      only : jprb
    use radiation_io,                  only : nulerr, radiation_abort
    use yomhook,                       only : lhook, dr_hook, jphook
    use radiation_config,              only : config_type
    use radiation_aerosol,             only : aerosol_type

    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in), target :: config
    type(aerosol_type),       intent(in)  :: aerosol
    ! Optical depth, single scattering albedo and asymmetry factor of
    ! the atmosphere (gases on input, gases and aerosols on output)
    ! for each g point. Note that longwave ssa and asymmetry and
    ! shortwave asymmetry are all zero for gases, so are not yet
    ! defined on input and are therefore intent(out).
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), &
         &   intent(inout) :: od_lw
    real(jprb), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol), &
         &   intent(out)   :: ssa_lw, g_lw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), &
         &   intent(inout) :: od_sw, ssa_sw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), &
         &   intent(out)   :: g_sw

    ! Temporary extinction and scattering optical depths of aerosol
    ! plus gas
    real(jprb) :: local_od, local_scat

    ! Extinction optical depth, scattering optical depth and
    ! asymmetry-times-scattering-optical-depth for all the aerosols at
    ! a point in space for each spectral band of the shortwave and
    ! longwave spectrum
    real(jprb), dimension(config%n_bands_sw,nlev) &
         & :: od_sw_aerosol, scat_sw_aerosol, scat_g_sw_aerosol
    real(jprb), dimension(config%n_bands_lw,nlev) :: od_lw_aerosol
    real(jprb), dimension(config%n_bands_lw_if_scattering,nlev) &
         & :: scat_lw_aerosol, scat_g_lw_aerosol

    ! Loop indices for column, level, g point and band
    integer :: jcol, jlev, jg, jb

    ! Range of levels over which aerosols are present
    integer :: istartlev, iendlev

    ! Indices to spectral band
    integer :: iband

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics:add_aerosol_optics_direct',0,hook_handle)

    if (config%do_sw) then
      ! Check array dimensions
      if (ubound(aerosol%od_sw,1) /= config%n_bands_sw) then
        write(nulerr,'(a,i0,a,i0)') '*** Error: aerosol%od_sw contains ', &
           &  ubound(aerosol%od_sw,1), ' band, expected ', &
           &  config%n_bands_sw
        call radiation_abort()
      end if

      istartlev = lbound(aerosol%od_sw,2)
      iendlev   = ubound(aerosol%od_sw,2)

      ! Set variables to zero that may not have been previously
      g_sw(:,:,istartcol:iendcol) = 0.0_jprb

      ! Loop over position
      do jcol = istartcol,iendcol
! Added for DWD (2020)
!NEC$ forced_collapse
        do jlev = istartlev,iendlev
          do jb = 1,config%n_bands_sw
            od_sw_aerosol(jb,jlev) = aerosol%od_sw(jb,jlev,jcol)
            scat_sw_aerosol(jb,jlev) = aerosol%ssa_sw(jb,jlev,jcol) * od_sw_aerosol(jb,jlev)
            scat_g_sw_aerosol(jb,jlev) = aerosol%g_sw(jb,jlev,jcol) * scat_sw_aerosol(jb,jlev)

            if (.not. config%do_sw_delta_scaling_with_gases) then
              ! Delta-Eddington scaling on aerosol only.  Note that if
              ! do_sw_delta_scaling_with_gases==.true. then the delta
              ! scaling is done to the cloud-aerosol-gas mixture
              ! inside the solver
              call delta_eddington_extensive(od_sw_aerosol(jb,jlev), scat_sw_aerosol(jb,jlev), &
                   &                         scat_g_sw_aerosol(jb,jlev))
            end if
          end do
        end do
        ! Combine aerosol shortwave scattering properties with gas
        ! properties (noting that any gas scattering will have an
        ! asymmetry factor of zero)
        do jlev = istartlev,iendlev
          if (od_sw_aerosol(1,jlev) > 0.0_jprb) then
            do jg = 1,config%n_g_sw
              iband = config%i_band_from_reordered_g_sw(jg)
              local_od = od_sw(jg,jlev,jcol) + od_sw_aerosol(iband,jlev)
              local_scat = ssa_sw(jg,jlev,jcol) * od_sw(jg,jlev,jcol) &
                   &  + scat_sw_aerosol(iband,jlev)
              ! Note that asymmetry_sw of gases is zero so the following
              ! simply weights the aerosol asymmetry by the scattering
              ! optical depth
              g_sw(jg,jlev,jcol) = scat_g_sw_aerosol(iband,jlev) / local_scat
              local_od = od_sw(jg,jlev,jcol) + od_sw_aerosol(iband,jlev)
              ssa_sw(jg,jlev,jcol) = local_scat / local_od
              od_sw (jg,jlev,jcol) = local_od
            end do
          end if
        end do
      end do

    end if

    if (config%do_lw) then

      if (ubound(aerosol%od_lw,1) /= config%n_bands_lw) then
        write(nulerr,'(a,i0,a,i0)') '*** Error: aerosol%od_lw contains ', &
           &  ubound(aerosol%od_lw,1), ' band, expected ', &
           &  config%n_bands_lw
        call radiation_abort()
      end if

      istartlev = lbound(aerosol%od_lw,2)
      iendlev   = ubound(aerosol%od_lw,2)

      if (config%do_lw_aerosol_scattering) then
        ssa_lw(:,:,istartcol:iendcol) = 0.0_jprb
        g_lw(:,:,istartcol:iendcol)   = 0.0_jprb
 
        ! Loop over position
        do jcol = istartcol,iendcol
! Added for DWD (2020)
!NEC$ forced_collapse
          do jlev = istartlev,iendlev
            do jb = 1,config%n_bands_lw
              od_lw_aerosol(jb,jlev) = aerosol%od_lw(jb,jlev,jcol)
              scat_lw_aerosol(jb,jlev) = aerosol%ssa_lw(jb,jlev,jcol) * od_lw_aerosol(jb,jlev)
              scat_g_lw_aerosol(jb,jlev) = aerosol%g_lw(jb,jlev,jcol) * scat_lw_aerosol(jb,jlev)
            
              call delta_eddington_extensive(od_lw_aerosol(jb,jlev), scat_lw_aerosol(jb,jlev), &
                   &                         scat_g_lw_aerosol(jb,jlev))
            end do
          end do
          do jlev = istartlev,iendlev
            do jg = 1,config%n_g_lw
              iband = config%i_band_from_reordered_g_lw(jg)
              if (od_lw_aerosol(iband,jlev) > 0.0_jprb) then
                ! All scattering is due to aerosols, therefore the
                ! asymmetry factor is equal to the value for aerosols
                if (scat_lw_aerosol(iband,jlev) > 0.0_jprb) then
                  g_lw(jg,jlev,jcol) = scat_g_lw_aerosol(iband,jlev) &
                       &  / scat_lw_aerosol(iband,jlev)
                end if
                local_od = od_lw(jg,jlev,jcol) + od_lw_aerosol(iband,jlev)
                ssa_lw(jg,jlev,jcol) = scat_lw_aerosol(iband,jlev) / local_od
                od_lw (jg,jlev,jcol) = local_od
              end if
            end do
          end do
        end do

      else ! No longwave scattering

        ! Loop over position
        do jcol = istartcol,iendcol
! Added for DWD (2020)
!NEC$ forced_collapse
          do jlev = istartlev,iendlev
            ! If aerosol longwave scattering is not included then we
            ! weight the optical depth by the single scattering
            ! co-albedo
            do jb = 1, config%n_bands_lw
              od_lw_aerosol(jb,jlev) = aerosol%od_lw(jb,jlev,jcol) &
                 &  * (1.0_jprb - aerosol%ssa_lw(jb,jlev,jcol))
            end do
          end do
          do jlev = istartlev,iendlev
            do jg = 1,config%n_g_lw
              od_lw(jg,jlev,jcol) = od_lw(jg,jlev,jcol) &
                   &  + od_lw_aerosol(config%i_band_from_reordered_g_lw(jg),jlev)
            end do
          end do
        end do

      end if
    end if


    if (lhook) call dr_hook('radiation_aerosol_optics:add_aerosol_optics_direct',1,hook_handle)

  end subroutine add_aerosol_optics_direct
 

  !---------------------------------------------------------------------
  ! Sometimes it is useful to specify aerosol in terms of its optical
  ! depth at a particular wavelength.  This function returns the dry
  ! mass-extinction coefficient, i.e. the extinction cross section per
  ! unit mass, for aerosol of type "itype" at the specified wavelength
  ! (m). For hydrophilic types, the value at the first relative
  ! humidity bin is taken.
  function dry_aerosol_mass_extinction(config, itype, wavelength)

    use parkind1,                      only : jprb
    use radiation_io,                  only : nulerr, radiation_abort
    use radiation_config,              only : config_type
    use radiation_aerosol_optics_data, only : aerosol_optics_type, &
         &  IAerosolClassUndefined,   IAerosolClassIgnored, &
         &  IAerosolClassHydrophobic, IAerosolClassHydrophilic

    type(config_type), intent(in), target :: config

    ! Aerosol type
    integer, intent(in) :: itype

    ! Wavelength (m)
    real(jprb), intent(in) :: wavelength
    
    real(jprb) :: dry_aerosol_mass_extinction

    ! Index to the monochromatic wavelength requested
    integer :: imono

    ! Pointer to the aerosol optics coefficients for brevity of access
    type(aerosol_optics_type), pointer :: ao

    ao => config%aerosol_optics

    imono = minloc(abs(wavelength - ao%wavelength_mono), 1)

    if (abs(wavelength - ao%wavelength_mono(imono))/wavelength > 0.02_jprb) then
      write(nulerr,'(a,e11.4,a)') '*** Error: requested wavelength ', &
           &  wavelength, ' not within 2% of stored wavelengths'
      call radiation_abort()
     end if

    if (ao%iclass(itype) == IAerosolClassHydrophobic) then
      dry_aerosol_mass_extinction = ao%mass_ext_mono_phobic(imono,ao%itype(itype))
    else if (ao%iclass(itype) == IAerosolClassHydrophilic) then
      ! Take the value at the first relative-humidity bin for the
      ! "dry" aerosol value
      dry_aerosol_mass_extinction = ao%mass_ext_mono_philic(imono,1,ao%itype(itype))
    else
      dry_aerosol_mass_extinction = 0.0_jprb
    end if

  end function dry_aerosol_mass_extinction


  !---------------------------------------------------------------------
  ! Compute aerosol extinction coefficient at a particular wavelength
  ! and a single height - this is useful for visibility diagnostics
  subroutine aerosol_extinction(ncol,istartcol,iendcol, &
       &  config, wavelength, mixing_ratio, relative_humidity, extinction)

    use parkind1,                      only : jprb
    use yomhook,                       only : lhook, dr_hook, jphook
    use radiation_io,                  only : nulerr, radiation_abort
    use radiation_config,              only : config_type
    use radiation_aerosol_optics_data, only : aerosol_optics_type, &
         &  IAerosolClassUndefined,   IAerosolClassIgnored, &
         &  IAerosolClassHydrophobic, IAerosolClassHydrophilic

    integer, intent(in) :: ncol               ! number of columns
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in), target :: config
    real(jprb), intent(in)  :: wavelength ! Requested wavelength (m)
    real(jprb), intent(in)  :: mixing_ratio(ncol,config%n_aerosol_types)
    real(jprb), intent(in)  :: relative_humidity(ncol)
    real(jprb), intent(out) :: extinction(ncol)

    ! Local aerosol extinction
    real(jprb) :: ext

    ! Index to the monochromatic wavelength requested
    integer :: imono

    ! Pointer to the aerosol optics coefficients for brevity of access
    type(aerosol_optics_type), pointer :: ao
    
    ! Loop indices for column and aerosol type
    integer :: jcol, jtype

    ! Relative humidity index
    integer :: irh

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics:aerosol_extinction',0,hook_handle)

    do jtype = 1,config%n_aerosol_types
      if (config%aerosol_optics%iclass(jtype) == IAerosolClassUndefined) then
        write(nulerr,'(a)') '*** Error: not all aerosol types are defined'
        call radiation_abort()
      end if
    end do

    ao => config%aerosol_optics

    imono = minloc(abs(wavelength - ao%wavelength_mono), 1)

    if (abs(wavelength - ao%wavelength_mono(imono))/wavelength > 0.02_jprb) then
      write(nulerr,'(a,e11.4,a)') '*** Error: requested wavelength ', &
           &  wavelength, ' not within 2% of stored wavelengths'
      call radiation_abort()
     end if

    ! Loop over position
    do jcol = istartcol,iendcol
      ext = 0.0_jprb
      ! Get relative-humidity index
      irh = ao%calc_rh_index(relative_humidity(jcol))
      ! Add extinction coefficients from each aerosol type
      do jtype = 1,config%n_aerosol_types
        if (ao%iclass(jtype) == IAerosolClassHydrophobic) then
          ext = ext + mixing_ratio(jcol,jtype) &
               &    * ao%mass_ext_mono_phobic(imono,ao%itype(jtype))
        else if (ao%iclass(jtype) == IAerosolClassHydrophilic) then
          ext = ext + mixing_ratio(jcol,jtype) &
               &    * ao%mass_ext_mono_philic(imono,irh,ao%itype(jtype))
        end if
      end do

      extinction(jcol) = ext
    end do

    if (lhook) call dr_hook('radiation_aerosol_optics:aerosol_extinction',1,hook_handle)

  end subroutine aerosol_extinction

end module radiation_aerosol_optics
