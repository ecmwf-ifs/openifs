! radiation_general_cloud_optics.F90 - Computing generalized cloud optical properties
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

module radiation_general_cloud_optics

  implicit none

  public
  
contains

  ! Provides elemental function "delta_eddington_scat_od"
#include "radiation_delta_eddington.h"


  !---------------------------------------------------------------------
  ! Load cloud scattering data; this subroutine delegates to one
  ! in radiation_general_cloud_optics_data.F90
  subroutine setup_general_cloud_optics(config)

    use parkind1,         only : jprb
    use yomhook,          only : lhook, dr_hook, jphook

    use radiation_io,     only : nulout
    use radiation_config, only : config_type, NMaxCloudTypes
    use radiation_spectral_definition, only : SolarReferenceTemperature, &
         &                                    TerrestrialReferenceTemperature

    type(config_type), intent(inout) :: config

    character(len=511) :: file_name

    integer :: jtype ! loop index
    integer :: strlen

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_general_cloud_optics:setup_general_cloud_optics',0,hook_handle)

    ! Count number of cloud types
    config%n_cloud_types = 0
    do jtype = 1,NMaxCloudTypes
      if (len_trim(config%cloud_type_name(jtype)) > 0) then
        config%n_cloud_types = jtype
      else
        exit
      end if
    end do

    ! If cloud_type_name has not been provided then assume liquid,ice
    ! noting that the default spectral averaging (defined in
    ! radiation_config.F90) is "thick"
    if (config%n_cloud_types == 0) then
      config%cloud_type_name(1) = "mie_droplet"
      config%cloud_type_name(2) = "baum-general-habit-mixture_ice"
      ! Optionally override spectral averaging method
      !config%use_thick_cloud_spectral_averaging(1) = .false.
      !config%use_thick_cloud_spectral_averaging(2) = .false.
      config%n_cloud_types = 2
    end if

    ! Allocate structures
    if (config%do_sw) then
      allocate(config%cloud_optics_sw(config%n_cloud_types))
    end if

    if (config%do_lw) then
      allocate(config%cloud_optics_lw(config%n_cloud_types))
    end if

    ! Load cloud optics data
    do jtype = 1,config%n_cloud_types
      if (config%cloud_type_name(jtype)(1:1) == '/') then
        file_name = trim(config%cloud_type_name(jtype))
      else
        strlen = len_trim(config%cloud_type_name(jtype))
        if (config%cloud_type_name(jtype)(strlen-2:strlen) == ".nc") then
          file_name = trim(config%directory_name) &
               &  // '/' // trim(config%cloud_type_name(jtype))
        else
          file_name = trim(config%directory_name) &
               &  // '/' // trim(config%cloud_type_name(jtype)) &
               &  // '_scattering.nc'
        end if
      end if

      if (config%do_sw) then
        if (config%iverbosesetup >= 2) then
          write(nulout,'(a,i0,a)') 'Shortwave cloud type ', jtype, ':'
        end if
        call config%cloud_optics_sw(jtype)%setup(file_name, &
             &  config%gas_optics_sw%spectral_def, &
             &  use_bands=(.not. config%do_cloud_aerosol_per_sw_g_point), &
             &  use_thick_averaging=config%use_thick_cloud_spectral_averaging(jtype), &
             &  weighting_temperature=SolarReferenceTemperature, &
             &  iverbose=config%iverbosesetup)
      end if

      if (config%do_lw) then
        if (config%iverbosesetup >= 2) then
          write(nulout,'(a,i0,a)') 'Longwave cloud type ', jtype, ':'
        end if
        call config%cloud_optics_lw(jtype)%setup(file_name, &
             &  config%gas_optics_lw%spectral_def, &
             &  use_bands=(.not. config%do_cloud_aerosol_per_lw_g_point), &
             &  use_thick_averaging=config%use_thick_cloud_spectral_averaging(jtype), &
             &  weighting_temperature=TerrestrialReferenceTemperature, &
             &  iverbose=config%iverbosesetup)
      end if

    end do

    if (lhook) call dr_hook('radiation_general_cloud_optics:setup_general_cloud_optics',1,hook_handle)

  end subroutine setup_general_cloud_optics

  !---------------------------------------------------------------------
  ! Compute cloud optical properties
  subroutine general_cloud_optics(nlev,istartcol,iendcol, &
       &  config, thermodynamics, cloud, & 
       &  od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
       &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    use radiation_io,     only : nulout
    use radiation_config, only : config_type
    use radiation_thermodynamics, only    : thermodynamics_type
    use radiation_cloud, only             : cloud_type
    use radiation_constants, only         : AccelDueToGravity
    !use radiation_general_cloud_optics_data, only : general_cloud_optics_type

    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in), target :: config
    type(thermodynamics_type),intent(in)  :: thermodynamics
    type(cloud_type),   intent(in)        :: cloud

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
    real(jprb), dimension(config%n_bands_sw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_sw_cloud, ssa_sw_cloud, g_sw_cloud

    ! In-cloud water path of one cloud type (kg m-2)
    real(jprb), dimension(istartcol:iendcol,nlev) :: water_path

    integer :: jtype, jcol, jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_general_cloud_optics:general_cloud_optics',0,hook_handle)

    if (config%iverbose >= 2) then
      write(nulout,'(a)') 'Computing cloud absorption/scattering properties'
    end if

    ! Array-wise assignment
    od_lw_cloud  = 0.0_jprb
    od_sw_cloud  = 0.0_jprb
    ssa_sw_cloud = 0.0_jprb
    g_sw_cloud   = 0.0_jprb
    if (config%do_lw_cloud_scattering) then
      ssa_lw_cloud = 0.0_jprb
      g_lw_cloud   = 0.0_jprb
    end if

    ! Loop over cloud types
    do jtype = 1,config%n_cloud_types
      ! Compute in-cloud water path
      if (config%is_homogeneous) then
        water_path = cloud%mixing_ratio(istartcol:iendcol,:,jtype) &
             &  *  (thermodynamics%pressure_hl(istartcol:iendcol, 2:nlev+1) &
             &     -thermodynamics%pressure_hl(istartcol:iendcol, 1:nlev)) &
             &  * (1.0_jprb / AccelDueToGravity)
      else
        water_path = cloud%mixing_ratio(istartcol:iendcol,:,jtype) &
             &  *  (thermodynamics%pressure_hl(istartcol:iendcol, 2:nlev+1) &
             &     -thermodynamics%pressure_hl(istartcol:iendcol, 1:nlev)) &
             &  * (1.0_jprb / (AccelDueToGravity &
             &                 * max(config%cloud_fraction_threshold, &
             &                       cloud%fraction(istartcol:iendcol,:))))
      end if

      ! Add optical properties to the cumulative total for the
      ! longwave and shortwave
      if (config%do_lw) then
        ! For the moment, we use ssa_lw_cloud and g_lw_cloud as
        ! containers for scattering optical depth and scattering
        ! coefficient x asymmetry factor, then scale after
        if (config%do_lw_cloud_scattering) then
          call config%cloud_optics_lw(jtype)%add_optical_properties(config%n_bands_lw, nlev, &
               &  iendcol+1-istartcol, cloud%fraction(istartcol:iendcol,:), &
               &  water_path, cloud%effective_radius(istartcol:iendcol,:,jtype), &
               &  od_lw_cloud, ssa_lw_cloud, g_lw_cloud)
        else
          call config%cloud_optics_lw(jtype)%add_optical_properties(config%n_bands_lw, nlev, &
               &  iendcol+1-istartcol, cloud%fraction(istartcol:iendcol,:), &
               &  water_path, cloud%effective_radius(istartcol:iendcol,:,jtype), od_lw_cloud)
        end if
      end if
      
      if (config%do_sw) then
        ! For the moment, we use ssa_sw_cloud and g_sw_cloud as
        ! containers for scattering optical depth and scattering
        ! coefficient x asymmetry factor, then scale after
        call config%cloud_optics_sw(jtype)%add_optical_properties(config%n_bands_sw, nlev, &
             &  iendcol+1-istartcol, cloud%fraction(istartcol:iendcol,:), &
             &  water_path, cloud%effective_radius(istartcol:iendcol,:,jtype), &
             &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud)
      end if
    end do

    ! Scale the combined longwave optical properties
    if (config%do_lw_cloud_scattering) then
      do jcol = istartcol, iendcol
        do jlev = 1,nlev
          if (cloud%fraction(jcol,jlev) > 0.0_jprb) then
            ! Note that original cloud optics does not do
            ! delta-Eddington scaling for liquid clouds in longwave
            call delta_eddington_extensive(od_lw_cloud(:,jlev,jcol), &
                 &  ssa_lw_cloud(:,jlev,jcol), g_lw_cloud(:,jlev,jcol))
            
            ! Scale to get asymmetry factor and single scattering albedo
            g_lw_cloud(:,jlev,jcol) = g_lw_cloud(:,jlev,jcol) &
                 &  / max(ssa_lw_cloud(:,jlev,jcol), 1.0e-15_jprb)
            ssa_lw_cloud(:,jlev,jcol) = ssa_lw_cloud(:,jlev,jcol) &
                 &  / max(od_lw_cloud(:,jlev,jcol),  1.0e-15_jprb)
          end if
        end do
      end do
    end if
    
    ! Scale the combined shortwave optical properties
    if (config%do_sw) then
      if (.not. config%do_sw_delta_scaling_with_gases) then
        do jcol = istartcol, iendcol
          do jlev = 1,nlev
            if (cloud%fraction(jcol,jlev) > 0.0_jprb) then
              call delta_eddington_extensive(od_sw_cloud(:,jlev,jcol), &
                   &  ssa_sw_cloud(:,jlev,jcol), g_sw_cloud(:,jlev,jcol))
            end if
          end do
        end do
      end if

      do jcol = istartcol, iendcol
        do jlev = 1,nlev
          if (cloud%fraction(jcol,jlev) > 0.0_jprb) then
            ! Scale to get asymmetry factor and single scattering albedo
            g_sw_cloud(:,jlev,jcol) = g_sw_cloud(:,jlev,jcol) &
                 &  / max(ssa_sw_cloud(:,jlev,jcol), 1.0e-15_jprb)
            ssa_sw_cloud(:,jlev,jcol) = ssa_sw_cloud(:,jlev,jcol) &
                 &  / max(od_sw_cloud(:,jlev,jcol),  1.0e-15_jprb)
          end if
        end do
      end do
    end if

    if (lhook) call dr_hook('radiation_general_cloud_optics:general_cloud_optics',1,hook_handle)

  end subroutine general_cloud_optics

end module radiation_general_cloud_optics
