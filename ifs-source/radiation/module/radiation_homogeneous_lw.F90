! radiation_homogeneous_lw.F90 - Longwave homogeneous-column (no cloud fraction) solver
!
! (C) Copyright 2016- ECMWF.
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
!   2017-04-11  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2019-01-14  R. Hogan  Save spectral flux profile if required

module radiation_homogeneous_lw

  public

contains

  !---------------------------------------------------------------------
  ! Longwave homogeneous solver, in which clouds are assumed to fill
  ! the gridbox horizontally
  subroutine solver_homogeneous_lw(nlev,istartcol,iendcol, &
       &  config, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    use radiation_config, only         : config_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type, indexed_sum_profile
    use radiation_two_stream, only     : calc_two_stream_gammas_lw, &
         &                               calc_reflectance_transmittance_lw, &
         &                               calc_no_scattering_transmittance_lw
    use radiation_adding_ica_lw, only  : adding_ica_lw, calc_fluxes_no_scattering_lw
    use radiation_lw_derivatives, only : calc_lw_derivatives_ica

    implicit none

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
         &  od
    real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
         &  ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each longwave band
    real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol)   :: &
         &  od_cloud
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &  nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

    ! Planck function at each half-level and the surface
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
         &  planck_hl
  
    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) &
         &  :: emission, albedo

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Diffuse reflectance and transmittance for each layer in clear
    ! and all skies
    real(jprb), dimension(config%n_g_lw, nlev) :: reflectance, transmittance

    ! Emission by a layer into the upwelling or downwelling diffuse
    ! streams, in clear and all skies
    real(jprb), dimension(config%n_g_lw, nlev) :: source_up, source_dn

    ! Fluxes per g point
    real(jprb), dimension(config%n_g_lw, nlev+1) :: flux_up, flux_dn

    ! Combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    real(jprb), dimension(config%n_g_lw) :: od_total, ssa_total, g_total

    ! Two-stream coefficients
    real(jprb), dimension(config%n_g_lw) :: gamma1, gamma2

    ! Optical depth of cloud in g-point space
    real(jprb), dimension(config%n_g_lw) :: od_cloud_g

    ! Is there any cloud in the profile?
    logical :: is_cloudy_profile

    ! Number of g points
    integer :: ng

    ! Loop indices for level and column
    integer :: jlev, jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_homogeneous_lw:solver_homogeneous_lw',0,hook_handle)

    ng = config%n_g_lw

    ! Loop through columns
    do jcol = istartcol,iendcol

      ! Is there any cloud in the profile?
      is_cloudy_profile = .false.
      do jlev = 1,nlev
        if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
          is_cloudy_profile = .true.
          exit
        end if
      end do

      ! If clear-sky fluxes need to be computed then we first compute
      ! the reflectance and transmittance of all layers, neglecting
      ! clouds. If clear-sky fluxes are not required then we only do
      ! the clear-sky layers since these will be needed when we come
      ! to do the total-sky fluxes.
      do jlev = 1,nlev
        if (config%do_clear .or. cloud%fraction(jcol,jlev) &
             &                 < config%cloud_fraction_threshold) then
          if (config%do_lw_aerosol_scattering) then
            ! Scattering case: first compute clear-sky reflectance,
            ! transmittance etc at each model level
            ssa_total = ssa(:,jlev,jcol)
            g_total   = g(:,jlev,jcol)
            call calc_two_stream_gammas_lw(ng, ssa_total, g_total, &
                 &  gamma1, gamma2)
            call calc_reflectance_transmittance_lw(ng, &
                 &  od(:,jlev,jcol), gamma1, gamma2, &
                 &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1,jcol), &
                 &  reflectance(:,jlev), transmittance(:,jlev), &
                 &  source_up(:,jlev), source_dn(:,jlev))
          else
            ! Non-scattering case: use simpler functions for
            ! transmission and emission
            call calc_no_scattering_transmittance_lw(ng, od(:,jlev,jcol), &
                 &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1, jcol), &
                 &  transmittance(:,jlev), source_up(:,jlev), source_dn(:,jlev))          
            ! Ensure that clear-sky reflectance is zero since it may be
            ! used in cloudy-sky case
            reflectance(:,jlev) = 0.0_jprb
          end if
         
        end if
      end do

      if (config%do_clear) then
        if (config%do_lw_aerosol_scattering) then
          ! Then use adding method to compute fluxes
          call adding_ica_lw(ng, nlev, &
               &  reflectance, transmittance, source_up, source_dn, &
               &  emission(:,jcol), albedo(:,jcol), &
               &  flux_up, flux_dn)
        else
          ! Simpler down-then-up method to compute fluxes
          call calc_fluxes_no_scattering_lw(ng, nlev, &
               &  transmittance, source_up, source_dn, &
               &  emission(:,jcol), albedo(:,jcol), &
               &  flux_up, flux_dn)
          
        end if

        ! Sum over g-points to compute broadband fluxes
        flux%lw_up_clear(jcol,:) = sum(flux_up,1)
        flux%lw_dn_clear(jcol,:) = sum(flux_dn,1)
        ! Store surface spectral downwelling fluxes
        flux%lw_dn_surf_clear_g(:,jcol) = flux_dn(:,nlev+1)

        ! Save the spectral fluxes if required
        if (config%do_save_spectral_flux) then
          call indexed_sum_profile(flux_up, config%i_spec_from_reordered_g_lw, &
               &                   flux%lw_up_clear_band(:,jcol,:))
          call indexed_sum_profile(flux_dn, config%i_spec_from_reordered_g_lw, &
               &                   flux%lw_dn_clear_band(:,jcol,:))
        end if

      end if ! Do clear-sky calculations

      ! Now the total-sky calculation.  If this is a clear profile and
      ! clear-sky fluxes have been calculated then we can simply copy
      ! over the clear-sky fluxes, otherwise we need to compute fluxes
      ! now.
      if (is_cloudy_profile .or. .not. config%do_clear) then
        do jlev = 1,nlev
          ! Compute combined gas+aerosol+cloud optical properties;
          ! note that for clear layers, the reflectance and
          ! transmittance have already been calculated
          if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
            od_cloud_g = od_cloud(config%i_band_from_reordered_g_lw,jlev,jcol)
            od_total = od(:,jlev,jcol) + od_cloud_g
            ssa_total = 0.0_jprb
            g_total   = 0.0_jprb

            if (config%do_lw_cloud_scattering) then
              ! Scattering case: calculate reflectance and
              ! transmittance at each model level
              if (config%do_lw_aerosol_scattering) then
                where (od_total > 0.0_jprb)
                  ssa_total = (ssa(:,jlev,jcol)*od(:,jlev,jcol) &
                       &     + ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     *  od_cloud_g) & 
                       &     / od_total
                end where
                where (ssa_total > 0.0_jprb .and. od_total > 0.0_jprb)
                  g_total = (g(:,jlev,jcol)*ssa(:,jlev,jcol)*od(:,jlev,jcol) &
                       &     +   g_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     * ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     *  od_cloud_g) &
                       &     / (ssa_total*od_total)
                end where
              else
                where (od_total > 0.0_jprb)
                  ssa_total = ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     * od_cloud_g / od_total
                end where
                where (ssa_total > 0.0_jprb .and. od_total > 0.0_jprb)
                  g_total = g_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     * ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     *  od_cloud_g / (ssa_total*od_total)
                end where
              end if
            
              ! Compute cloudy-sky reflectance, transmittance etc at
              ! each model level
              call calc_two_stream_gammas_lw(ng, ssa_total, g_total, &
                   &  gamma1, gamma2)
              call calc_reflectance_transmittance_lw(ng, &
                   &  od_total, gamma1, gamma2, &
                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1,jcol), &
                   &  reflectance(:,jlev), transmittance(:,jlev), &
                   &  source_up(:,jlev), source_dn(:,jlev))
            else
              ! No-scattering case: use simpler functions for
              ! transmission and emission
              call calc_no_scattering_transmittance_lw(ng, od_total, &
                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1, jcol), &
                   &  transmittance(:,jlev), source_up(:,jlev), source_dn(:,jlev))
            end if
          end if ! is cloudy layer
        end do
        
        if (config%do_lw_cloud_scattering) then
          ! Use adding method to compute fluxes for an overcast sky
          call adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
               &  emission(:,jcol), albedo(:,jcol), &
               &  flux_up, flux_dn)
        else
          ! Simpler down-then-up method to compute fluxes
          call calc_fluxes_no_scattering_lw(ng, nlev, &
               &  transmittance, source_up, source_dn, emission(:,jcol), albedo(:,jcol), &
               &  flux_up, flux_dn)
        end if
        
        ! Store overcast broadband fluxes
        flux%lw_up(jcol,:) = sum(flux_up,1)
        flux%lw_dn(jcol,:) = sum(flux_dn,1)
        ! Store surface spectral downwelling fluxes
        flux%lw_dn_surf_g(:,jcol) = flux_dn(:,nlev+1)

        ! Save the spectral fluxes if required
        if (config%do_save_spectral_flux) then
          call indexed_sum_profile(flux_up, config%i_spec_from_reordered_g_lw, &
               &                   flux%lw_up_band(:,jcol,:))
          call indexed_sum_profile(flux_dn, config%i_spec_from_reordered_g_lw, &
               &                   flux%lw_dn_band(:,jcol,:))
        end if

      else
        ! No cloud in profile and clear-sky fluxes already
        ! calculated: copy them over
        flux%lw_up(jcol,:) = flux%lw_up_clear(jcol,:)
        flux%lw_dn(jcol,:) = flux%lw_dn_clear(jcol,:)
        flux%lw_dn_surf_g(:,jcol) = flux%lw_dn_surf_clear_g(:,jcol)

        if (config%do_save_spectral_flux) then
          flux%lw_up_band(:,jcol,:) = flux%lw_up_clear_band(:,jcol,:)
          flux%lw_dn_band(:,jcol,:) = flux%lw_dn_clear_band(:,jcol,:)
        end if

     end if

      ! Compute the longwave derivatives needed by Hogan and Bozzo
      ! (2015) approximate radiation update scheme, using clear-sky
      ! transmittance if no clouds were present in the profile,
      ! all-sky transmittance otherwise
      if (config%do_lw_derivatives) then
        call calc_lw_derivatives_ica(ng, nlev, jcol, transmittance, flux_up(:,nlev+1), &
             &                       flux%lw_derivatives)
 
      end if

    end do

    if (lhook) call dr_hook('radiation_homogeneous_lw:solver_homogeneous_lw',1,hook_handle)
    
  end subroutine solver_homogeneous_lw

end module radiation_homogeneous_lw
