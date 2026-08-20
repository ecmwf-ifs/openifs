! (C) Copyright 2019- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! radiation_cloudless_lw.F90 - Longwave homogeneous cloudless solver
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!

module radiation_cloudless_lw

public :: solver_cloudless_lw

contains

  !---------------------------------------------------------------------
  ! Longwave homogeneous solver containing no clouds
  subroutine solver_cloudless_lw(nlev,istartcol,iendcol, &
       &  config, od, ssa, g, planck_hl, emission, albedo, flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    use radiation_config, only         : config_type
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

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
         &  od
    real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
         &  ssa, g

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

    ! Combined optical depth, single scattering albedo and asymmetry
    ! factor
    real(jprb), dimension(config%n_g_lw) :: ssa_total, g_total

    ! Two-stream coefficients
    real(jprb), dimension(config%n_g_lw) :: gamma1, gamma2

    ! Number of g points
    integer :: ng

    ! Loop indices for level and column
    integer :: jlev, jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloudless_lw:solver_cloudless_lw',0,hook_handle)

    ng = config%n_g_lw

    ! Loop through columns
    do jcol = istartcol,iendcol

      ! Compute the reflectance and transmittance of all layers,
      ! neglecting clouds
      do jlev = 1,nlev
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
          ! Ensure that clear-sky reflectance is zero
          reflectance(:,jlev) = 0.0_jprb
        end if
      end do

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

      if (config%do_clear) then
        ! Clear-sky calculations are equal to all-sky for this solver:
        ! copy fluxes over
        flux%lw_up_clear(jcol,:) = flux%lw_up(jcol,:)
        flux%lw_dn_clear(jcol,:) = flux%lw_dn(jcol,:)
        flux%lw_dn_surf_clear_g(:,jcol) = flux%lw_dn_surf_g(:,jcol)
        if (config%do_save_spectral_flux) then
          flux%lw_up_clear_band(:,jcol,:) = flux%lw_up_band(:,jcol,:)
          flux%lw_dn_clear_band(:,jcol,:) = flux%lw_dn_band(:,jcol,:)
        end if
      end if

      ! Compute the longwave derivatives needed by Hogan and Bozzo
      ! (2015) approximate radiation update scheme
      if (config%do_lw_derivatives) then
        call calc_lw_derivatives_ica(ng, nlev, jcol, transmittance, flux_up(:,nlev+1), &
             &                       flux%lw_derivatives)
       end if

    end do

    if (lhook) call dr_hook('radiation_cloudless_lw:solver_cloudless_lw',1,hook_handle)

  end subroutine solver_cloudless_lw

end module radiation_cloudless_lw
