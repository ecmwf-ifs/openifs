! radiation_save.F90 - Save data to NetCDF files
!
! Copyright (C) 2014-2019 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! Modifications
!   2017-04-22  R. Hogan  Adapt for new way of describing longwave properties
!   2019-01-02  R. Hogan  Only save cloud properties if do_clouds==.true.

module radiation_save

  use parkind1, only : jprb

  implicit none

  ! Two available subroutines: save final fluxes and save intermediate
  ! radiative properties
  public :: save_fluxes, save_radiative_properties, save_inputs

contains

  !---------------------------------------------------------------------
  ! Save fluxes in "flux" to NetCDF file_name, plus pressure from the
  ! thermodynamics object
  subroutine save_fluxes(file_name, config, thermodynamics, flux, &
       &                 iverbose)

    use yomhook,  only           : lhook, dr_hook, jphook
    use easy_netcdf

    use radiation_io,             only : nulout
    use radiation_config,         only : config_type, IGasModelMonochromatic
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_flux,           only : flux_type

    character(len=*), intent(in)         :: file_name
    type(config_type), intent(in)        :: config
    type(thermodynamics_type), intent(in):: thermodynamics
    type(flux_type), intent(in)          :: flux

    type(netcdf_file)                    :: out_file
    integer                              :: ncol, n_lev_plus1

    character(5), parameter :: default_lw_units_str = 'W m-2'
    character(5)            :: lw_units_str

    integer, optional, intent(in) :: iverbose
    integer                       :: i_local_verbose

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_save:save_fluxes',0,hook_handle)
    
    if (present(iverbose)) then
      i_local_verbose = iverbose
    else
      i_local_verbose = config%iverbose
    end if

    ! Work out array dimensions
    if (config%do_sw) then
      ncol = size(flux%sw_up,1)
      n_lev_plus1 = size(flux%sw_up,2)
    elseif (config%do_lw) then
      ncol = size(flux%lw_up,1)
      n_lev_plus1 = size(flux%lw_up,2)
    else
      if (i_local_verbose >= 1) then
        write(nulout,'(a,a,a)') 'Warning: neither longwave nor shortwave computed so ', &
             &                  file_name,' not written'
      end if
      return
    end if

    if (config%i_gas_model == IGasModelMonochromatic &
         .and. config%mono_lw_wavelength > 0.0_jprb) then
      lw_units_str = 'W m-3'
    else
      lw_units_str = default_lw_units_str
    end if

    ! Open the file
    call out_file%create(trim(file_name), iverbose=i_local_verbose)

    ! Variables stored internally with column varying fastest, but in
    ! output file column varies most slowly so need to transpose
    call out_file%transpose_matrices(.true.)

    ! Spectral fluxes in memory are dimensioned (nband,ncol,nlev), but
    ! are reoriented in the output file to be (nband,nlev,ncol), where
    ! the convention here is first dimension varying fastest
    call out_file%permute_3d_arrays( (/ 1, 3, 2 /) )

    ! Define dimensions
    call out_file%define_dimension("column", ncol)
    call out_file%define_dimension("half_level", n_lev_plus1)

    if (config%do_save_spectral_flux) then
      if (config%do_lw) then
        call out_file%define_dimension("band_lw", config%n_spec_lw)
      end if
      if (config%do_sw) then
        call out_file%define_dimension("band_sw", config%n_spec_sw)
      end if
    else if (config%do_surface_sw_spectral_flux) then
      if (config%do_sw) then
        call out_file%define_dimension("band_sw", config%n_bands_sw)
      end if
    end if

    if (config%do_lw .and. config%do_canopy_fluxes_lw) then
      call out_file%define_dimension("canopy_band_lw", &
           &  size(flux%lw_dn_surf_canopy, 1))
    end if
    if (config%do_sw .and. config%do_canopy_fluxes_sw) then
      call out_file%define_dimension("canopy_band_sw", &
           &  size(flux%sw_dn_diffuse_surf_canopy, 1))
    end if

    ! Put global attributes
    call out_file%put_global_attributes(title_str &
         & = "Radiative flux profiles from SPARTACUS offline radiation model")

    ! Define variables
    call out_file%define_variable("pressure_hl", &
         &   dim2_name="column", dim1_name="half_level", &
         &   units_str="Pa", long_name="Pressure")

    if (config%do_lw) then
      call out_file%define_variable("flux_up_lw", &
           &   dim2_name="column", dim1_name="half_level", &
           &   units_str=lw_units_str, long_name="Upwelling longwave flux")
      call out_file%define_variable("flux_dn_lw", &
           &   dim2_name="column", dim1_name="half_level", &
           &   units_str=lw_units_str, long_name="Downwelling longwave flux")
      if (config%do_clear) then
        call out_file%define_variable("flux_up_lw_clear", &
             &   dim2_name="column", dim1_name="half_level", &
             &   units_str=lw_units_str, &
             &   long_name="Upwelling clear-sky longwave flux")
        call out_file%define_variable("flux_dn_lw_clear", &
             &   dim2_name="column", dim1_name="half_level", &
             &   units_str=lw_units_str, &
             &   long_name="Downwelling clear-sky longwave flux")
      end if

      if (config%do_lw_derivatives) then
        call out_file%define_variable("lw_derivative", &
             &  dim2_name="column", dim1_name="half_level", &
             &  units_str="1", &
             &  long_name="Derivative of upwelling LW flux w.r.t. surface value")
      end if

      if (config%do_save_spectral_flux) then
        call out_file%define_variable("spectral_flux_up_lw", &
             &   dim3_name="column", dim2_name="half_level", &
             &   dim1_name="band_lw", units_str=lw_units_str, &
             &   long_name="Spectral upwelling longwave flux")
        call out_file%define_variable("spectral_flux_dn_lw", &
             &   dim3_name="column", dim2_name="half_level", &
             &   dim1_name="band_lw", units_str=lw_units_str, &
             &   long_name="Spectral downwelling longwave flux")
        if (config%do_clear) then
          call out_file%define_variable("spectral_flux_up_lw_clear", &
               &   dim3_name="column", dim2_name="half_level", &
               &   dim1_name="band_lw", units_str=lw_units_str, &
               &   long_name="Spectral upwelling clear-sky longwave flux")
          call out_file%define_variable("spectral_flux_dn_lw_clear", &
               &   dim3_name="column", dim2_name="half_level", &
               &   dim1_name="band_lw", units_str=lw_units_str, &
               &   long_name="Spectral downwelling clear-sky longwave flux")
        end if
      end if
   
      if (config%do_canopy_fluxes_lw) then
        call out_file%define_variable("canopy_flux_dn_lw_surf", &
             &   dim2_name="column", dim1_name="canopy_band_lw", units_str=lw_units_str, &
             &   long_name="Surface downwelling longwave flux in canopy bands")
      end if

    end if

    if (config%do_sw) then
      call out_file%define_variable("flux_up_sw", &
           &   dim2_name="column", dim1_name="half_level", &
           &   units_str="W m-2", long_name="Upwelling shortwave flux")
      call out_file%define_variable("flux_dn_sw", &
           &   dim2_name="column", dim1_name="half_level", &
           &   units_str="W m-2", long_name="Downwelling shortwave flux")
      if (config%do_sw_direct) then
        call out_file%define_variable("flux_dn_direct_sw", &
             &   dim2_name="column", dim1_name="half_level", &
             &   units_str="W m-2", &
             &   long_name="Downwelling direct shortwave flux")
      end if
      if (config%do_clear) then
        call out_file%define_variable("flux_up_sw_clear", &
             &   dim2_name="column", dim1_name="half_level", &
             &   units_str="W m-2", &
             &   long_name="Upwelling clear-sky shortwave flux")
        call out_file%define_variable("flux_dn_sw_clear", &
             &   dim2_name="column", dim1_name="half_level", &
             &   units_str="W m-2", &
             &   long_name="Downwelling clear-sky shortwave flux")
        if (config%do_sw_direct) then
          call out_file%define_variable("flux_dn_direct_sw_clear", &
               &   dim2_name="column", dim1_name="half_level", &
               &   units_str="W m-2", &
               &   long_name="Downwelling clear-sky direct shortwave flux")
        end if
      end if

      if (config%do_save_spectral_flux) then
        call out_file%define_variable("spectral_flux_up_sw", &
             &   dim3_name="column", dim2_name="half_level", &
             &   dim1_name="band_sw", units_str="W m-2", &
             &   long_name="Spectral upwelling shortwave flux")
        call out_file%define_variable("spectral_flux_dn_sw", &
             &   dim3_name="column", dim2_name="half_level", &
             &   dim1_name="band_sw", units_str="W m-2", &
             &   long_name="Spectral downwelling shortwave flux")
        if (config%do_sw_direct) then
          call out_file%define_variable("spectral_flux_dn_direct_sw", &
               &   dim3_name="column", dim2_name="half_level", &
               &   dim1_name="band_sw", units_str="W m-2", &
               &   long_name="Spectral downwelling direct shortwave flux")
        end if
        if (config%do_clear) then
          call out_file%define_variable("spectral_flux_up_sw_clear", &
               &   dim3_name="column", dim2_name="half_level", &
               &   dim1_name="band_sw", units_str="W m-2", &
               &   long_name="Spectral upwelling clear-sky shortwave flux")
          call out_file%define_variable("spectral_flux_dn_sw_clear", &
               &   dim3_name="column", dim2_name="half_level", &
               &   dim1_name="band_sw", units_str="W m-2", &
               &   long_name="Spectral downwelling clear-sky shortwave flux")
          if (config%do_sw_direct) then
            call out_file%define_variable("spectral_flux_dn_direct_sw_clear", &
                 &   dim3_name="column", dim2_name="half_level", &
                 &   dim1_name="band_sw", units_str="W m-2", &
                 &   long_name="Spectral downwelling clear-sky direct shortwave flux")
          end if
        end if
      else if (config%do_surface_sw_spectral_flux) then
        call out_file%define_variable("spectral_flux_dn_sw_surf", &
             &   dim2_name="column", dim1_name="band_sw", units_str="W m-2", &
             &   long_name="Spectral downwelling shortwave flux at surface")
        call out_file%define_variable("spectral_flux_dn_direct_sw_surf", &
             &   dim2_name="column", dim1_name="band_sw", units_str="W m-2", &
             &   long_name="Spectral downwelling direct shortwave flux at surface")
        if (config%do_clear) then
          call out_file%define_variable("spectral_flux_dn_sw_surf_clear", &
               &   dim2_name="column", dim1_name="band_sw", units_str="W m-2", &
               &   long_name="Spectral downwelling clear-sky shortwave flux at surface")
          call out_file%define_variable("spectral_flux_dn_direct_sw_surf_clear", &
               &   dim2_name="column", dim1_name="band_sw", units_str="W m-2", &
               &   long_name="Spectral downwelling clear-sky direct shortwave flux at surface")
        end if
      end if
   
      if (config%do_canopy_fluxes_sw) then
        call out_file%define_variable("canopy_flux_dn_diffuse_sw_surf", &
             &   dim2_name="column", dim1_name="canopy_band_sw", units_str="W m-2", &
             &   long_name="Surface downwelling diffuse shortwave flux in canopy bands")
        call out_file%define_variable("canopy_flux_dn_direct_sw_surf", &
             &   dim2_name="column", dim1_name="canopy_band_sw", units_str="W m-2", &
             &   long_name="Surface downwelling direct shortwave flux in canopy bands")
      end if

    end if
   
    if (config%do_lw) then
      call out_file%define_variable("cloud_cover_lw", &
           &  dim1_name="column", units_str="1", &
           &  long_name="Total cloud cover diagnosed by longwave solver")
    end if
    if (config%do_sw) then
      call out_file%define_variable("cloud_cover_sw", &
           &  dim1_name="column", units_str="1", &
           &  long_name="Total cloud cover diagnosed by shortwave solver")
    end if

    ! Write variables

    call out_file%put("pressure_hl", thermodynamics%pressure_hl)

    if (config%do_lw) then
      call out_file%put("flux_up_lw", flux%lw_up)
      call out_file%put("flux_dn_lw", flux%lw_dn)
      if (config%do_clear) then
        call out_file%put("flux_up_lw_clear", flux%lw_up_clear)
        call out_file%put("flux_dn_lw_clear", flux%lw_dn_clear)
      end if

      if (config%do_lw_derivatives) then
        call out_file%put("lw_derivative", flux%lw_derivatives)
      end if

      if (config%do_save_spectral_flux) then
        call out_file%put("spectral_flux_up_lw", flux%lw_up_band)
        call out_file%put("spectral_flux_dn_lw", flux%lw_dn_band)
        if (config%do_clear) then
          call out_file%put("spectral_flux_up_lw_clear", flux%lw_up_clear_band)
          call out_file%put("spectral_flux_dn_lw_clear", flux%lw_dn_clear_band)
        end if
      end if

      if (config%do_canopy_fluxes_lw) then
        call out_file%put("canopy_flux_dn_lw_surf", flux%lw_dn_surf_canopy, &
             &            do_transp = .false.)
      end if

    end if

    if (config%do_sw) then
      call out_file%put("flux_up_sw", flux%sw_up)
      call out_file%put("flux_dn_sw", flux%sw_dn)
      if (config%do_sw_direct) then
        call out_file%put("flux_dn_direct_sw", flux%sw_dn_direct)
      end if
      if (config%do_clear) then
        call out_file%put("flux_up_sw_clear", flux%sw_up_clear)
        call out_file%put("flux_dn_sw_clear", flux%sw_dn_clear)
        if (config%do_sw_direct) then
          call out_file%put("flux_dn_direct_sw_clear", flux%sw_dn_direct_clear)
        end if
      end if

      if (config%do_save_spectral_flux) then
        call out_file%put("spectral_flux_up_sw", flux%sw_up_band)
        call out_file%put("spectral_flux_dn_sw", flux%sw_dn_band)
        if (config%do_sw_direct) then
          call out_file%put("spectral_flux_dn_direct_sw", &
               &   flux%sw_dn_direct_band)
        end if
        if (config%do_clear) then
          call out_file%put("spectral_flux_up_sw_clear", flux%sw_up_clear_band)
          call out_file%put("spectral_flux_dn_sw_clear", flux%sw_dn_clear_band)
          if (config%do_sw_direct) then
            call out_file%put("spectral_flux_dn_direct_sw_clear", &
                 &   flux%sw_dn_direct_clear_band)
          end if
        end if
      else if (config%do_surface_sw_spectral_flux) then
        call out_file%put("spectral_flux_dn_sw_surf", flux%sw_dn_surf_band, &
               &   do_transp=.false.)
        call out_file%put("spectral_flux_dn_direct_sw_surf", flux%sw_dn_direct_surf_band, &
               &   do_transp=.false.)
        if (config%do_clear) then
          call out_file%put("spectral_flux_dn_sw_surf_clear", flux%sw_dn_surf_clear_band, &
               &   do_transp=.false.)
          call out_file%put("spectral_flux_dn_direct_sw_surf_clear", &
               &            flux%sw_dn_direct_surf_clear_band, do_transp=.false.)
        end if
      end if

      if (config%do_canopy_fluxes_sw) then
        call out_file%put("canopy_flux_dn_diffuse_sw_surf", flux%sw_dn_diffuse_surf_canopy, &
             &            do_transp = .false.)
        call out_file%put("canopy_flux_dn_direct_sw_surf",  flux%sw_dn_direct_surf_canopy, &
             &            do_transp = .false.)
      end if

    end if

    if (config%do_lw) then
      call out_file%put("cloud_cover_lw", flux%cloud_cover_lw)
    end if
    if (config%do_sw) then
      call out_file%put("cloud_cover_sw", flux%cloud_cover_sw)
    end if

    ! Close file
    call out_file%close()

    if (lhook) call dr_hook('radiation_save:save_fluxes',1,hook_handle)

  end subroutine save_fluxes
  

  !---------------------------------------------------------------------
  ! Save intermediate radiative properties, specifically the
  ! scattering and absorption properties at each g-point/band
  subroutine save_radiative_properties(file_name, nlev, &
       &  istartcol, iendcol, &
       &  config, single_level, thermodynamics, cloud, &
       &  planck_hl, lw_emission, lw_albedo, &
       &  sw_albedo_direct, sw_albedo_diffuse, &
       &  incoming_sw, &
       &  od_lw, ssa_lw, g_lw, &
       &  od_sw, ssa_sw, g_sw, &
       &  od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
       &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud)

    use radiation_config,        only : config_type
    use radiation_single_level,  only : single_level_type
    use radiation_thermodynamics,only : thermodynamics_type
    use radiation_cloud,         only : cloud_type
    use easy_netcdf

    character(len=*),         intent(in) :: file_name
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(cloud_type),         intent(in) :: cloud

    integer, intent(in) :: nlev, istartcol, iendcol

    ! Input variables, as defined in radiation_interface.F90

    ! Layer optical depth, single scattering albedo and asymmetry factor of
    ! gases and aerosols at each shortwave g-point
    real(jprb), intent(in), dimension(config%n_g_sw,nlev,istartcol:iendcol) :: od_sw, ssa_sw, g_sw

   ! Layer optical depth, single scattering albedo and asymmetry factor of
    ! hydrometeors in each shortwave band
    real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
         &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud

    ! Direct and diffuse surface albedo, and the incoming shortwave
    ! flux into a plane perpendicular to the incoming radiation at
    ! top-of-atmosphere in each of the shortwave g-points
    real(jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol) &
         &  :: sw_albedo_direct, sw_albedo_diffuse, incoming_sw

    ! Layer optical depth, single scattering albedo and asymmetry factor of
    ! gases and aerosols at each longwave g-point, where the latter
    ! two variables are only defined if aerosol longwave scattering is
    ! enabled (otherwise both are treated as zero).
    real(jprb), intent(in), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: od_lw
    real(jprb), intent(in), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: &
         &  ssa_lw, g_lw

    ! Layer optical depth, single scattering albedo and asymmetry factor of
    ! hydrometeors in each longwave band, where the latter two
    ! variables are only defined if hydrometeor longwave scattering is
    ! enabled (otherwise both are treated as zero).
    real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol) :: od_lw_cloud
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol) :: &
         &  ssa_lw_cloud, g_lw_cloud

    ! The Planck function (emitted flux from a black body) at half
    ! levels and at the surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: planck_hl

    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), dimension(config%n_g_lw, istartcol:iendcol) :: lw_emission, lw_albedo

    integer :: n_col_local ! Number of columns from istartcol to iendcol

    ! Object for output NetCDF file
    type(netcdf_file) :: out_file

    n_col_local = iendcol + 1 - istartcol

    ! Alas the NetCDF library is not thread-safe for writing, so we
    ! must write radiative-property files serially

    !$OMP CRITICAL

    ! Open the file
    call out_file%create(trim(file_name), iverbose=config%iverbose)

    ! Configure matrix and 3D-array orientation
    call out_file%transpose_matrices(.true.)

    ! Sometimes the Planck function values are very large or small
    ! even if the fluxes are within a manageable range
    call out_file%double_precision(.true.)

    ! Define dimensions
    !    call out_file%define_dimension("column", n_col_local)
    call out_file%define_dimension("column", 0) ! "Unlimited" dimension
    call out_file%define_dimension("level", nlev)
    call out_file%define_dimension("half_level", nlev+1)
    call out_file%define_dimension("level_interface", nlev-1)

    if (config%do_lw) then
      call out_file%define_dimension("gpoint_lw", config%n_g_lw) 
      call out_file%define_dimension("band_lw", config%n_bands_lw) 
    end if
    if (config%do_sw) then
      call out_file%define_dimension("gpoint_sw", config%n_g_sw) 
      call out_file%define_dimension("band_sw", config%n_bands_sw) 
    end if

    ! Define variables
    call out_file%define_variable("pressure_hl", &
         &  dim2_name="column", dim1_name="half_level", &
         &  units_str="Pa", long_name="Pressure on half-levels")

    if (allocated(thermodynamics%h2o_sat_liq)) then
      call out_file%define_variable("q_sat_liquid", &
           &  dim2_name="column", dim1_name="level", &
           &  units_str="kg kg-1", long_name="Specific humidity at liquid saturation")
    end if

    if (config%do_sw) then
      call out_file%define_variable("cos_solar_zenith_angle", &
           &  dim1_name="column", units_str="1", &
           &  long_name="Cosine of the solar zenith angle")
    end if

    if (config%do_clouds) then
      call out_file%define_variable("cloud_fraction", &
           &  dim2_name="column", dim1_name="level", &
           &  units_str="1", long_name="Cloud fraction")
      call out_file%define_variable("overlap_param", &
           &  dim2_name="column", dim1_name="level_interface", &
           &  units_str="1", long_name="Cloud overlap parameter")
    end if

    if (config%do_lw) then
      call out_file%define_variable("planck_hl", &
           &  dim3_name="column", dim2_name="half_level", dim1_name="gpoint_lw", &
           &  units_str="W m-2", long_name="Planck function on half-levels")
      call out_file%define_variable("lw_emission", &
           &  dim2_name="column", dim1_name="gpoint_lw", &
           &  units_str="W m-2", long_name="Longwave surface emission")
      call out_file%define_variable("lw_emissivity", &
           &  dim2_name="column", dim1_name="gpoint_lw", &
           &  units_str="1", long_name="Surface longwave emissivity")

      call out_file%define_variable("od_lw", &
           &  dim3_name="column", dim2_name="level", dim1_name="gpoint_lw", &
           &  units_str="1", long_name="Clear-sky longwave optical depth")
      if (config%do_lw_aerosol_scattering) then
        call out_file%define_variable("ssa_lw", &
           &  dim3_name="column", dim2_name="level", dim1_name="gpoint_lw", &
           &  units_str="1", long_name="Clear-sky longwave single scattering albedo")
        call out_file%define_variable("asymmetry_lw", &
           &  dim3_name="column", dim2_name="level", dim1_name="gpoint_lw", &
           &  units_str="1", long_name="Clear-sky longwave asymmetry factor")
      end if

      if (config%do_clouds) then
        call out_file%define_variable("od_lw_cloud", &
             &  dim3_name="column", dim2_name="level", dim1_name="band_lw", &
             &  units_str="1", long_name="In-cloud longwave optical depth")
        if (config%do_lw_cloud_scattering) then
          call out_file%define_variable("ssa_lw_cloud", &
               &  dim3_name="column", dim2_name="level", dim1_name="band_lw", &
               &  units_str="1", long_name="Cloud longwave single scattering albedo")
          call out_file%define_variable("asymmetry_lw_cloud", &
               &  dim3_name="column", dim2_name="level", dim1_name="band_lw", &
               &  units_str="1", long_name="Cloud longwave asymmetry factor")
        end if
      end if ! do_clouds
    end if ! do_lw
    
    if (config%do_sw) then
      call out_file%define_variable("incoming_sw", &
           &  dim2_name="column", dim1_name="gpoint_sw", &
           &  units_str="W m-2", long_name="Incoming shortwave flux at top-of-atmosphere in direction of sun")

      call out_file%define_variable("sw_albedo", &
           &  dim2_name="column", dim1_name="gpoint_sw", &
           &  units_str="1", long_name="Surface shortwave albedo to diffuse radiation")
      call out_file%define_variable("sw_albedo_direct", &
           &  dim2_name="column", dim1_name="gpoint_sw", &
           &  units_str="1", long_name="Surface shortwave albedo to direct radiation")

      call out_file%define_variable("od_sw", &
           &  dim3_name="column", dim2_name="level", dim1_name="gpoint_sw", &
           &  units_str="1", long_name="Clear-sky shortwave optical depth")
      call out_file%define_variable("ssa_sw", &
           &  dim3_name="column", dim2_name="level", dim1_name="gpoint_sw", &
           &  units_str="1", long_name="Clear-sky shortwave single scattering albedo")
      call out_file%define_variable("asymmetry_sw", &
           &  dim3_name="column", dim2_name="level", dim1_name="gpoint_sw", &
           &  units_str="1", long_name="Clear-sky shortwave asymmetry factor")

      if (config%do_clouds) then
        call out_file%define_variable("od_sw_cloud", &
             &  dim3_name="column", dim2_name="level", dim1_name="band_sw", &
             &  units_str="1", long_name="In-cloud shortwave optical depth")
        call out_file%define_variable("ssa_sw_cloud", &
             &  dim3_name="column", dim2_name="level", dim1_name="band_sw", &
             &  units_str="1", long_name="Cloud shortwave single scattering albedo")
        call out_file%define_variable("asymmetry_sw_cloud", &
             &  dim3_name="column", dim2_name="level", dim1_name="band_sw", &
             &  units_str="1", long_name="Cloud shortwave asymmetry factor")
      end if
    end if
   
    if (config%do_clouds) then
      if (allocated(cloud%fractional_std)) then
        call out_file%define_variable("fractional_std", &
             &  dim2_name="column", dim1_name="level", units_str="1", &
             &  long_name="Fractional standard deviation of cloud optical depth")
      end if
      if (allocated(cloud%inv_cloud_effective_size)) then
        call out_file%define_variable("inv_cloud_effective_size", &
             &  dim2_name="column", dim1_name="level", units_str="m-1", &
             &  long_name="Inverse of cloud effective horizontal size")
      end if
      if (allocated(cloud%inv_inhom_effective_size)) then
        call out_file%define_variable("inv_inhom_effective_size", &
             &  dim2_name="column", dim1_name="level", units_str="m-1", &
             &  long_name="Inverse of cloud inhomogeneity effective horizontal size")
      end if
   end if

    ! Write variables
    call out_file%put("pressure_hl", thermodynamics%pressure_hl(istartcol:iendcol,:))

    if (allocated(thermodynamics%h2o_sat_liq)) then
      call out_file%put("q_sat_liquid", thermodynamics%h2o_sat_liq(istartcol:iendcol,:))
    end if

    if (config%do_clouds) then
      call out_file%put("cloud_fraction", cloud%fraction(istartcol:iendcol,:))
      call out_file%put("overlap_param", cloud%overlap_param(istartcol:iendcol,:))
    end if

    if (config%do_sw) then
      call out_file%put("cos_solar_zenith_angle", single_level%cos_sza(istartcol:iendcol))
      call out_file%put("sw_albedo", sw_albedo_diffuse, do_transp=.false.)
      call out_file%put("sw_albedo_direct", sw_albedo_direct, do_transp=.false.)
    end if

    if (config%do_lw) then
      call out_file%put("lw_emissivity", 1.0_jprb - lw_albedo, do_transp=.false.)
      call out_file%put("planck_hl", planck_hl)
      call out_file%put("lw_emission", lw_emission, do_transp=.false.)
      
      call out_file%put("od_lw", od_lw)
      if (config%do_lw_aerosol_scattering) then
        call out_file%put("ssa_lw", ssa_lw)
        call out_file%put("asymmetry_lw", g_lw)
      end if
      
      if (config%do_clouds) then
        call out_file%put("od_lw_cloud", od_lw_cloud)
        if (config%do_lw_cloud_scattering) then
          call out_file%put("ssa_lw_cloud", ssa_lw_cloud)
          call out_file%put("asymmetry_lw_cloud", g_lw_cloud)
        end if
      end if
    end if
    
    if (config%do_sw) then
      call out_file%put("incoming_sw", incoming_sw, do_transp=.false.)
      
      call out_file%put("od_sw", od_sw)
      call out_file%put("ssa_sw", ssa_sw)
      call out_file%put("asymmetry_sw", g_sw)
      
      if (config%do_clouds) then
        call out_file%put("od_sw_cloud", od_sw_cloud)
        call out_file%put("ssa_sw_cloud", ssa_sw_cloud)
        call out_file%put("asymmetry_sw_cloud", g_sw_cloud)
      end if
    end if
    
    if (config%do_clouds) then
      if (allocated(cloud%fractional_std)) then
        call out_file%put("fractional_std", cloud%fractional_std(istartcol:iendcol,:))
      end if
      if (allocated(cloud%inv_cloud_effective_size)) then
        call out_file%put("inv_cloud_effective_size", cloud%inv_cloud_effective_size(istartcol:iendcol,:))
      end if
      if (allocated(cloud%inv_inhom_effective_size)) then
        call out_file%put("inv_inhom_effective_size", cloud%inv_inhom_effective_size(istartcol:iendcol,:))
      end if
    end if

    ! Close the file
    call out_file%close()

    !$OMP END CRITICAL
    
  end subroutine save_radiative_properties
  

  !---------------------------------------------------------------------
  ! Save inputs to the radiation scheme
  subroutine save_inputs(file_name, config, single_level, thermodynamics, &
       &                 gas, cloud, aerosol, lat, lon, iverbose)
    use yomhook,  only           : lhook, dr_hook, jphook
    use radiation_config,         only : config_type
    use radiation_single_level,   only : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_gas
    use radiation_cloud,          only : cloud_type
    use radiation_aerosol,        only : aerosol_type
    use easy_netcdf

    character(len=*),             intent(in)   :: file_name
    type(config_type),            intent(in)   :: config
    type(single_level_type),      intent(in)   :: single_level
    type(thermodynamics_type),    intent(in)   :: thermodynamics
    type(gas_type),               intent(inout):: gas
    type(cloud_type),             intent(in)   :: cloud
    type(aerosol_type), optional, intent(in)   :: aerosol
    real(jprb),         optional, intent(in)   :: lat(:), lon(:)
    integer,            optional, intent(in)   :: iverbose

    real(jprb), allocatable :: mixing_ratio(:,:)
    real(jprb), allocatable :: aerosol_mmr(:,:,:)
    real(jprb), allocatable :: seed(:)
    integer       :: i_local_verbose
    integer       :: ncol, nlev
    integer       :: jgas
    character(32) :: var_name, long_name

    ! Object for output NetCDF file
    type(netcdf_file) :: out_file

    logical :: do_aerosol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_save:save_inputs',0,hook_handle)

    if (present(iverbose)) then
      i_local_verbose = iverbose
    else
      i_local_verbose = config%iverbose
    end if

    ! Work out array dimensions
    ncol = size(thermodynamics%pressure_hl,1)
    nlev = size(thermodynamics%pressure_hl,2)
    nlev = nlev - 1
    
    do_aerosol = config%use_aerosols .and. present(aerosol)

    ! Open the file
    call out_file%create(trim(file_name), iverbose=i_local_verbose)

    ! Variables stored internally with column varying fastest, but in
    ! output file column varies most slowly so need to transpose
    call out_file%transpose_matrices(.true.)

    ! In the case of aerosols we convert dimensions (ncol,nlev,ntype)
    ! in memory to (nlev,ntype,ncol) in file (in both cases the first
    ! dimension varying fastest).
    call out_file%permute_3d_arrays( (/ 2, 3, 1 /) ) ! For aerosols

    ! Define dimensions
    call out_file%define_dimension("column",     ncol)
    call out_file%define_dimension("level",      nlev)
    call out_file%define_dimension("half_level", nlev+1)
    if (allocated(cloud%overlap_param)) then
      call out_file%define_dimension("level_interface", nlev-1)
    end if
    call out_file%define_dimension("sw_albedo_band", &
         &                         size(single_level%sw_albedo,2))
    call out_file%define_dimension("lw_emissivity_band", &
         &                         size(single_level%lw_emissivity,2))

    if (do_aerosol) then
      call out_file%define_dimension("aerosol_type", size(aerosol%mixing_ratio,3))
    end if

    ! Put global attributes
    call out_file%put_global_attributes( &
         &   title_str="Input profiles to radiation scheme")

    ! Define single-level variables
    call out_file%define_variable("solar_irradiance", &
         &   units_str="W m-2", long_name="Solar irradiance at Earth's orbit")
    if (present(lat)) then
      call out_file%define_variable("lat", &
           &   dim1_name="column", units_str="degrees_north", long_name="Latitude")
    end if
    if (present(lon)) then
      call out_file%define_variable("lon", &
           &   dim1_name="column", units_str="degrees_east", long_name="Longitude")
    end if
    call out_file%define_variable("skin_temperature", &
         &   dim1_name="column", units_str="K", long_name="Skin_temperature")
    call out_file%define_variable("cos_solar_zenith_angle", &
         &   dim1_name="column", units_str="1", &
         &   long_name="Cosine of the solar zenith angle")

    if (allocated(single_level%sw_albedo_direct)) then
      call out_file%define_variable("sw_albedo", &
           &   dim2_name="column", dim1_name="sw_albedo_band", &
           &   units_str="1", long_name="Shortwave surface albedo to diffuse radiation")
            call out_file%define_variable("sw_albedo_direct", &
           &   dim2_name="column", dim1_name="sw_albedo_band", &
           &   units_str="1", long_name="Shortwave surface albedo to direct radiation")
    else
      call out_file%define_variable("sw_albedo", &
           &   dim2_name="column", dim1_name="sw_albedo_band", &
           &   units_str="1", long_name="Shortwave surface albedo")
      
    end if
    call out_file%define_variable("lw_emissivity", &
         &   dim2_name="column", dim1_name="lw_emissivity_band", &
         &   units_str="1", long_name="Longwave surface emissivity")

    if (allocated(single_level%iseed)) then
      call out_file%define_variable("iseed", &
           &   dim1_name="column", units_str="1", is_double=.true., &
           &   long_name="Seed for random-number generator")
    end if

    ! Define thermodynamic variables on half levels
    call out_file%define_variable("pressure_hl", &
         &   dim2_name="column", dim1_name="half_level", &
         &   units_str="Pa", long_name="Pressure")
    call out_file%define_variable("temperature_hl", &
         &   dim2_name="column", dim1_name="half_level", &
         &   units_str="K", long_name="Temperature")

    ! Define gas mixing ratios on full levels
    call out_file%define_variable("q", &
         &   dim2_name="column", dim1_name="level", &
         &   units_str="1", long_name="Specific humidity")
    call out_file%define_variable("o3_mmr", &
         &   dim2_name="column", dim1_name="level", &
         &   units_str="1", long_name="Ozone mass mixing ratio")
    do jgas = 1,NMaxGases
      if (gas%is_present(jgas) .and. jgas /= IH2O .and. jgas /= IO3) then
        write(var_name,'(a,a)') trim(GasLowerCaseName(jgas)), '_vmr'
        write(long_name,'(a,a)') trim(GasName(jgas)), ' volume mixing ratio'
        call out_file%define_variable(trim(var_name), &
             &   dim2_name="column", dim1_name="level", &
             &   units_str="1", long_name=trim(long_name))
      end if
    end do

    if (config%do_clouds) then
      ! Define cloud variables on full levels
      call out_file%define_variable("cloud_fraction", &
           &   dim2_name="column", dim1_name="level", &
           &   units_str="1", long_name="Cloud fraction")
      call out_file%define_variable("q_liquid", &
           &   dim2_name="column", dim1_name="level", &
           &   units_str="1", long_name="Gridbox-mean liquid water mixing ratio")
      call out_file%define_variable("q_ice", &
           &   dim2_name="column", dim1_name="level", &
           &   units_str="1", long_name="Gridbox-mean ice water mixing ratio")
      call out_file%define_variable("re_liquid", &
           &   dim2_name="column", dim1_name="level", &
           &   units_str="m", long_name="Ice effective radius")
      if (allocated(cloud%re_ice)) then
        call out_file%define_variable("re_ice", &
             &   dim2_name="column", dim1_name="level", &
             &   units_str="m", long_name="Ice effective radius")
      end if
      if (allocated(cloud%overlap_param)) then
        call out_file%define_variable("overlap_param", &
             &  dim2_name="column", dim1_name="level_interface", &
             &  units_str="1", long_name="Cloud overlap parameter")
      end if
      if (allocated(cloud%fractional_std)) then
        call out_file%define_variable("fractional_std", &
             &  dim2_name="column", dim1_name="level", units_str="1", &
             &  long_name="Fractional standard deviation of cloud optical depth")
      end if
      if (allocated(cloud%inv_cloud_effective_size)) then
        call out_file%define_variable("inv_cloud_effective_size", &
             &   dim2_name="column", dim1_name="level", units_str="m-1", &
             &   long_name="Inverse of cloud effective horizontal size")
      end if
      if (allocated(cloud%inv_inhom_effective_size)) then
        call out_file%define_variable("inv_inhom_effective_size", &
             &  dim2_name="column", dim1_name="level", units_str="m-1", &
             &  long_name="Inverse of cloud inhomogeneity effective horizontal size")
      end if
    end if ! do_clouds

    ! Define aerosol mass mixing ratio
    if (do_aerosol) then
      call out_file%define_variable("aerosol_mmr", &
           &   dim3_name="column", dim2_name="aerosol_type", &
           &   dim1_name="level", units_str="kg kg-1", &
           &   long_name="Aerosol mass mixing ratio")
    end if

    ! Write variables
    call out_file%put("solar_irradiance", single_level%solar_irradiance)
    if (present(lat)) then
      call out_file%put("lat", lat)
    end if
    if (present(lon)) then
      call out_file%put("lon", lon)
    end if
    call out_file%put("skin_temperature", single_level%skin_temperature)
    call out_file%put("cos_solar_zenith_angle", single_level%cos_sza)
    call out_file%put("sw_albedo", single_level%sw_albedo)
    if (allocated(single_level%sw_albedo_direct)) then
      call out_file%put("sw_albedo_direct", single_level%sw_albedo_direct)
    end if
    call out_file%put("lw_emissivity", single_level%lw_emissivity)
    if (config%do_clouds .and. allocated(single_level%iseed)) then
      allocate(seed(ncol))
      seed = single_level%iseed
      call out_file%put("iseed", seed)
      deallocate(seed)
    end if

    call out_file%put("pressure_hl", thermodynamics%pressure_hl)
    call out_file%put("temperature_hl", thermodynamics%temperature_hl)

    allocate(mixing_ratio(ncol,nlev))
    call gas%get(IH2O, IMassMixingRatio, mixing_ratio)
    call out_file%put("q", mixing_ratio)
    call gas%get(IO3, IMassMixingRatio, mixing_ratio)
    call out_file%put("o3_mmr", mixing_ratio)
    do jgas = 1,NMaxGases
      if (gas%is_present(jgas) .and. jgas /= IH2O .and. jgas /= IO3) then
        write(var_name,'(a,a)') trim(GasLowerCaseName(jgas)), '_vmr'
        call gas%get(jgas, IVolumeMixingRatio, mixing_ratio)
        call out_file%put(trim(var_name), mixing_ratio)
      end if
    end do
    deallocate(mixing_ratio)

    if (config%do_clouds) then
      call out_file%put("cloud_fraction", cloud%fraction)
      call out_file%put("q_liquid", cloud%q_liq)
      call out_file%put("q_ice", cloud%q_ice)
      call out_file%put("re_liquid", cloud%re_liq)
      if (allocated(cloud%re_ice)) then
        call out_file%put("re_ice", cloud%re_ice)
      end if
      if (allocated(cloud%overlap_param)) then
        call out_file%put("overlap_param", cloud%overlap_param)
      end if
      if (allocated(cloud%fractional_std)) then
        call out_file%put("fractional_std", cloud%fractional_std)
      end if
      if (allocated(cloud%inv_cloud_effective_size)) then
        call out_file%put("inv_cloud_effective_size", cloud%inv_cloud_effective_size)
      end if
      if (allocated(cloud%inv_inhom_effective_size)) then
        call out_file%put("inv_inhom_effective_size", cloud%inv_inhom_effective_size)
      end if
    end if

    if (do_aerosol) then
      allocate(aerosol_mmr(ncol, nlev, size(aerosol%mixing_ratio,3)))
      aerosol_mmr = 0.0_jprb
      aerosol_mmr(:,aerosol%istartlev:aerosol%iendlev,:) = aerosol%mixing_ratio
      call out_file%put("aerosol_mmr", aerosol_mmr)
      deallocate(aerosol_mmr)
    end if

    ! Close the file
    call out_file%close()

    if (lhook) call dr_hook('radiation_save:save_inputs',1,hook_handle)
    
  end subroutine save_inputs

end module radiation_save
