! radiation_tripleclouds_lw.F90 - Longwave "Tripleclouds" solver
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
!   2017-04-28  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-10-08  R. Hogan  Call calc_region_properties
!   2020-09-18  R. Hogan  Replaced some array expressions with loops
!   2020-09-19  R. Hogan  Implement the cloud-only-scattering optimization

module radiation_tripleclouds_lw

  public

contains
  ! Small routine for scaling cloud optical depth in the cloudy
  ! regions
#include "radiation_optical_depth_scaling.h"

  !---------------------------------------------------------------------
  ! This module contains just one subroutine, the longwave
  ! "Tripleclouds" solver in which cloud inhomogeneity is treated by
  ! dividing each model level into three regions, one clear and two
  ! cloudy (with differing optical depth). This approach was described
  ! by Shonk and Hogan (2008).
  subroutine solver_tripleclouds_lw(nlev,istartcol,iendcol, &
       &  config, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

!    use radiation_io, only             : nulout
    use radiation_config, only         : config_type, IPdfShapeGamma
    use radiation_cloud, only          : cloud_type
    use radiation_regions, only        : calc_region_properties
    use radiation_overlap, only        : calc_overlap_matrices
    use radiation_flux, only           : flux_type, indexed_sum
    use radiation_matrix, only         : singlemat_x_vec
    use radiation_two_stream, only     : calc_ref_trans_lw, &
         &                               calc_no_scattering_transmittance_lw
    use radiation_adding_ica_lw, only  : adding_ica_lw, calc_fluxes_no_scattering_lw
    use radiation_lw_derivatives, only : calc_lw_derivatives_region

    implicit none

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth of each layer at each longwave
    ! g-point
    real(jprb), intent(in), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: od

    ! Gas and aerosol single-scattering albedo and asymmetry factor,
    ! only if longwave scattering by aerosols is to be represented
    real(jprb), intent(in), &
         &  dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: ssa, g

    ! Cloud and precipitation optical depth of each layer in each
    ! longwave band
    real(jprb), intent(in) :: od_cloud(config%n_bands_lw,nlev,istartcol:iendcol)

    ! Cloud and precipitation single-scattering albedo and asymmetry
    ! factor, only if longwave scattering by clouds is to be
    ! represented
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &                            nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

    ! Planck function (emitted flux from a black body) at half levels
    ! and at the surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: planck_hl

    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: emission, albedo

    ! Optical depth, single scattering albedo and asymmetry factor in
    ! each g-point including gas, aerosol and clouds
    real(jprb), dimension(config%n_g_lw) :: od_total, ssa_total, g_total

    ! Modified optical depth after Tripleclouds scaling to represent
    ! cloud inhomogeneity
    real(jprb), dimension(config%n_g_lw) :: od_cloud_new

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local constants
    integer, parameter :: nregions = 3

    ! In a clear-sky layer this will be 1, otherwise equal to nregions
    integer :: nreg

    ! Local variables
    
    ! The area fractions of each region
    real(jprb) :: region_fracs(1:nregions,nlev,istartcol:iendcol)

    ! The scaling used for the optical depth in the cloudy regions
    real(jprb) :: od_scaling(2:nregions,nlev,istartcol:iendcol)

    ! Directional overlap matrices defined at all layer interfaces
    ! including top-of-atmosphere and the surface
    real(jprb), dimension(nregions,nregions,nlev+1, &
         &                istartcol:iendcol) :: u_matrix, v_matrix

    ! Diffuse reflection and transmission matrices of each layer
    real(jprb), dimension(config%n_g_lw, nregions, nlev) :: reflectance, transmittance

    ! Emission by a layer into the upwelling or downwelling diffuse
    ! streams
    real(jprb), dimension(config%n_g_lw, nregions, nlev) &
         &  :: source_up, source_dn

    ! Clear-sky reflectance and transmittance
    real(jprb), dimension(config%n_g_lw, nlev) &
         &  :: ref_clear, trans_clear

    ! ...clear-sky equivalent
    real(jprb), dimension(config%n_g_lw, nlev) &
         &  :: source_up_clear, source_dn_clear

    ! Total albedo of the atmosphere/surface just above a layer
    ! interface with respect to downwelling diffuse radiation at that
    ! interface, where level index = 1 corresponds to the
    ! top-of-atmosphere
    real(jprb), dimension(config%n_g_lw, nregions, nlev+1) :: total_albedo

    ! Upwelling radiation just above a layer interface due to emission
    ! below that interface, where level index = 1 corresponds to the
    ! top-of-atmosphere
    real(jprb), dimension(config%n_g_lw, nregions, nlev+1) :: total_source

    ! Total albedo and source of the atmosphere just below a layer interface
    real(jprb), dimension(config%n_g_lw, nregions) &
         &  :: total_albedo_below, total_source_below

    ! Downwelling flux below and above an interface between
    ! layers into a plane perpendicular to the direction of the sun
    real(jprb), dimension(config%n_g_lw, nregions) &
         &  :: flux_dn, flux_dn_below, flux_up

    ! ...clear-sky equivalent (no distinction between "above/below")
    real(jprb), dimension(config%n_g_lw, nlev+1) &
         &  :: flux_dn_clear, flux_up_clear

    ! Clear-sky equivalent, but actually its reciprocal to replace
    ! some divisions by multiplications
    real(jprb), dimension(config%n_g_lw, nregions) :: inv_denom

    ! Identify clear-sky layers, with pseudo layers for outer space
    ! and below the ground, both treated as single-region clear skies
    logical :: is_clear_sky_layer(0:nlev+1)

    ! Index of the highest cloudy layer
    integer :: i_cloud_top

    integer :: jcol, jlev, jg, jreg, jreg2, ng

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_tripleclouds_lw:solver_tripleclouds_lw',0,hook_handle)

    ! --------------------------------------------------------
    ! Section 1: Prepare general variables and arrays
    ! --------------------------------------------------------
    ! Copy array dimensions to local variables for convenience
    ng   = config%n_g_lw

    ! Compute the wavelength-independent region fractions and
    ! optical-depth scalings
    call calc_region_properties(nlev,nregions,istartcol,iendcol, &
         &  config%i_cloud_pdf_shape == IPdfShapeGamma, &
         &  cloud%fraction, cloud%fractional_std, region_fracs, &
         &  od_scaling, config%cloud_fraction_threshold)

    ! Compute wavelength-independent overlap matrices u_matrix and v_matrix
    call calc_overlap_matrices(nlev,nregions,istartcol,iendcol, &
         &  region_fracs, cloud%overlap_param, &
         &  u_matrix, v_matrix, &
         &  decorrelation_scaling=config%cloud_inhom_decorr_scaling, &
         &  cloud_fraction_threshold=config%cloud_fraction_threshold, &
         &  use_beta_overlap=config%use_beta_overlap, &
         &  cloud_cover=flux%cloud_cover_lw)

    ! Main loop over columns
    do jcol = istartcol, iendcol
      ! --------------------------------------------------------
      ! Section 2: Prepare column-specific variables and arrays
      ! --------------------------------------------------------

      ! Define which layers contain cloud; assume that
      ! cloud%crop_cloud_fraction has already been called
      is_clear_sky_layer = .true.
      i_cloud_top = nlev+1
      do jlev = 1,nlev
        if (cloud%fraction(jcol,jlev) > 0.0_jprb) then
          is_clear_sky_layer(jlev) = .false.
          ! Get index to the first cloudy layer from the top
          if (i_cloud_top > jlev) then
            i_cloud_top = jlev
          end if
        end if
      end do
      if (config%do_lw_aerosol_scattering) then
        ! This is actually the first layer in which we need to
        ! consider scattering
        i_cloud_top = 1
      end if

      ! --------------------------------------------------------
      ! Section 3: Clear-sky calculation
      ! --------------------------------------------------------

      if (.not. config%do_lw_aerosol_scattering) then
        ! No scattering in clear-sky flux calculation; note that here
        ! the first two dimensions of the input arrays are unpacked
        ! into vectors inside the routine        
        call calc_no_scattering_transmittance_lw(ng*nlev, od(:,:,jcol), &
             &  planck_hl(:,1:nlev,jcol), planck_hl(:,2:nlev+1, jcol), &
             &  trans_clear, source_up_clear, source_dn_clear)
        ! Ensure that clear-sky reflectance is zero since it may be
        ! used in cloudy-sky case
        ref_clear = 0.0_jprb
        ! Simple down-then-up method to compute fluxes
        call calc_fluxes_no_scattering_lw(ng, nlev, &
             &  trans_clear, source_up_clear, source_dn_clear, &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up_clear, flux_dn_clear)
      else
        ! Scattering in clear-sky flux calculation
        call calc_ref_trans_lw(ng*nlev, &
             &  od(:,:,jcol), ssa(:,:,jcol), g(:,:,jcol), &
             &  planck_hl(:,1:jlev,jcol), planck_hl(:,2:jlev+1,jcol), &
             &  ref_clear, trans_clear, &
             &  source_up_clear, source_dn_clear)
        ! Use adding method to compute fluxes
        call adding_ica_lw(ng, nlev, &
             &  ref_clear, trans_clear, source_up_clear, source_dn_clear, &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up_clear, flux_dn_clear)
      end if

      if (config%do_clear) then
        ! Sum over g-points to compute broadband fluxes
        flux%lw_up_clear(jcol,:) = sum(flux_up_clear,1)
        flux%lw_dn_clear(jcol,:) = sum(flux_dn_clear,1)
        ! Store surface spectral downwelling fluxes / TOA upwelling
        flux%lw_dn_surf_clear_g(:,jcol) = flux_dn_clear(:,nlev+1)
        flux%lw_up_toa_clear_g (:,jcol) = flux_up_clear(:,1)
        ! Save the spectral fluxes if required
        if (config%do_save_spectral_flux) then
          do jlev = 1,nlev+1
            call indexed_sum(flux_up_clear(:,jlev), &
                 &           config%i_spec_from_reordered_g_lw, &
                 &           flux%lw_up_clear_band(:,jcol,jlev))
            call indexed_sum(flux_dn_clear(:,jlev), &
                 &           config%i_spec_from_reordered_g_lw, &
                 &           flux%lw_dn_clear_band(:,jcol,jlev))
          end do
        end if
      end if

      ! --------------------------------------------------------
      ! Section 4: Loop over cloudy layers to compute reflectance and transmittance
      ! --------------------------------------------------------
      ! In this section the reflectance, transmittance and sources
      ! are computed for each layer
      
      ! Firstly, ensure clear-sky transmittance is valid for whole
      ! depth of the atmosphere, because even above cloud it is used
      ! by the LW derivatives
      transmittance(:,1,:) = trans_clear(:,:)
      ! Dummy values in cloudy regions above cloud top
      if (i_cloud_top > 0) then
        transmittance(:,2:,1:min(i_cloud_top,nlev)) = 1.0_jprb
      end if

      do jlev = i_cloud_top,nlev ! Start at cloud top and work down

        ! Copy over clear-sky properties
        reflectance(:,1,jlev)    = ref_clear(:,jlev)
        source_up(:,1,jlev)      = source_up_clear(:,jlev) ! Scaled later by region size
        source_dn(:,1,jlev)      = source_dn_clear(:,jlev) ! Scaled later by region size
        nreg = nregions
        if (is_clear_sky_layer(jlev)) then
          nreg = 1
          reflectance(:,2:,jlev)   = 0.0_jprb
          transmittance(:,2:,jlev) = 1.0_jprb
          source_up(:,2:,jlev)     = 0.0_jprb
          source_dn(:,2:,jlev)     = 0.0_jprb
        else
          do jreg = 2,nreg
            ! Cloudy sky
            ! Add scaled cloud optical depth to clear-sky value
            od_cloud_new = od_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                 &       * od_scaling(jreg,jlev,jcol)
            od_total = od(:,jlev,jcol) + od_cloud_new

            if (config%do_lw_cloud_scattering) then
              ssa_total = 0.0_jprb
              g_total   = 0.0_jprb              
              if (config%do_lw_aerosol_scattering) then
                where (od_total > 0.0_jprb)
                  ssa_total = (ssa(:,jlev,jcol)*od(:,jlev,jcol) &
                       &     + ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     *  od_cloud_new) & 
                       &     / od_total
                end where
                where (ssa_total > 0.0_jprb .and. od_total > 0.0_jprb)
                  g_total = (g(:,jlev,jcol)*ssa(:,jlev,jcol)*od(:,jlev,jcol) &
                       &     +   g_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     * ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     *  od_cloud_new) &
                       &     / (ssa_total*od_total)
                end where
              else
                where (od_total > 0.0_jprb)
                  ssa_total = ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     * od_cloud_new / od_total
                end where
                where (ssa_total > 0.0_jprb .and. od_total > 0.0_jprb)
                  g_total = g_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     * ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     *  od_cloud_new / (ssa_total*od_total)
                end where
              end if
              call calc_ref_trans_lw(ng, &
                   &  od_total, ssa_total, g_total, &
                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1,jcol), &
                   &  reflectance(:,jreg,jlev), transmittance(:,jreg,jlev), &
                   &  source_up(:,jreg,jlev), source_dn(:,jreg,jlev))
            else
              ! No-scattering case: use simpler functions for
              ! transmission and emission
              call calc_no_scattering_transmittance_lw(ng, od_total, &
                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1, jcol), &
                   &  transmittance(:,jreg,jlev), source_up(:,jreg,jlev), source_dn(:,jreg,jlev))
              reflectance(:,jreg,jlev) = 0.0_jprb
            end if
          end do
          ! Emission is scaled by the size of each region
          do jreg = 1,nregions
            source_up(:,jreg,jlev) = region_fracs(jreg,jlev,jcol) * source_up(:,jreg,jlev)
            source_dn(:,jreg,jlev) = region_fracs(jreg,jlev,jcol) * source_dn(:,jreg,jlev)
          end do
        end if

      end do ! Loop over levels

      ! --------------------------------------------------------
      ! Section 5: Compute total sources and albedos at each half level
      ! --------------------------------------------------------

      total_albedo(:,:,:) = 0.0_jprb
      total_source(:,:,:) = 0.0_jprb

      ! Calculate the upwelling radiation emitted by the surface, and
      ! copy the surface albedo into total_albedo 
      do jreg = 1,nregions
        do jg = 1,ng
          ! region_fracs(jreg,nlev,jcol) is the fraction of each region in the
          ! lowest model level
          total_source(jg,jreg,nlev+1) = region_fracs(jreg,nlev,jcol)*emission(jg,jcol)
          total_albedo(jg,jreg,nlev+1) = albedo(jg,jcol)
        end do
      end do

      ! Work up from the surface computing the total albedo of the
      ! atmosphere and the total upwelling due to emission below each
      ! level below using the adding method
      do jlev = nlev,i_cloud_top,-1

        total_albedo_below        = 0.0_jprb

        if (is_clear_sky_layer(jlev)) then
          total_albedo_below = 0.0_jprb
          total_source_below = 0.0_jprb
          do jg = 1,ng
            inv_denom(jg,1) = 1.0_jprb &
                 &  / (1.0_jprb - total_albedo(jg,1,jlev+1)*reflectance(jg,1,jlev))
            total_albedo_below(jg,1) = reflectance(jg,1,jlev) &
                 &  + transmittance(jg,1,jlev)*transmittance(jg,1,jlev)*total_albedo(jg,1,jlev+1) &
                 &  * inv_denom(jg,1)
            total_source_below(jg,1) = source_up(jg,1,jlev) &
                 &  + transmittance(jg,1,jlev)*(total_source(jg,1,jlev+1) &
                 &  + total_albedo(jg,1,jlev+1)*source_dn(jg,1,jlev)) &
                 &  * inv_denom(jg,1)
          end do
        else
          inv_denom = 1.0_jprb / (1.0_jprb - total_albedo(:,:,jlev+1)*reflectance(:,:,jlev))
          total_albedo_below = reflectance(:,:,jlev) &
               &  + transmittance(:,:,jlev)*transmittance(:,:,jlev)*total_albedo(:,:,jlev+1) &
               &  * inv_denom
          total_source_below = source_up(:,:,jlev) &
               &  + transmittance(:,:,jlev)*(total_source(:,:,jlev+1) &
               &  + total_albedo(:,:,jlev+1)*source_dn(:,:,jlev)) &
               &  * inv_denom
        end if

        ! Account for cloud overlap when converting albedo below a
        ! layer interface to the equivalent values just above
        if (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev-1)) then
          total_albedo(:,:,jlev) = total_albedo_below(:,:)
          total_source(:,:,jlev) = total_source_below(:,:)
        else
          total_source(:,:,jlev) = singlemat_x_vec(ng,ng,nregions,&
               &  u_matrix(:,:,jlev,jcol), total_source_below)
          ! Use overlap matrix and exclude "anomalous" horizontal
          ! transport described by Shonk & Hogan (2008).  Therefore,
          ! the operation we perform is essentially diag(total_albedo)
          ! = matmul(transpose(v_matrix), diag(total_albedo_below)).
          do jreg = 1,nregions
            do jreg2 = 1,nregions
              total_albedo(:,jreg,jlev) &
                   &  = total_albedo(:,jreg,jlev) &
                   &  + total_albedo_below(:,jreg2) &
                   &  * v_matrix(jreg2,jreg,jlev,jcol)

            end do
          end do

        end if
        
      end do ! Reverse loop over levels

      ! --------------------------------------------------------
      ! Section 6: Copy over downwelling fluxes above cloud top
      ! --------------------------------------------------------
      do jlev = 1,i_cloud_top
        if (config%do_clear) then
          ! Clear-sky fluxes have already been averaged: use these
          flux%lw_dn(jcol,jlev) = flux%lw_dn_clear(jcol,jlev)
          if (config%do_save_spectral_flux) then
            flux%lw_dn_band(:,jcol,jlev) = flux%lw_dn_clear_band(:,jcol,jlev)
          end if
        else
          flux%lw_dn(jcol,:) = sum(flux_dn_clear(:,jlev))
          if (config%do_save_spectral_flux) then
            call indexed_sum(flux_dn_clear(:,jlev), &
                 &           config%i_spec_from_reordered_g_lw, &
                 &           flux%lw_dn_band(:,jcol,jlev))
          end if
        end if
      end do

      ! --------------------------------------------------------
      ! Section 7: Compute fluxes up to top-of-atmosphere
      ! --------------------------------------------------------

      ! Compute the fluxes just above the highest cloud
      flux_up(:,1) = total_source(:,1,i_cloud_top) &
           &  + total_albedo(:,1,i_cloud_top)*flux_dn_clear(:,i_cloud_top)
      flux_up(:,2:) = 0.0_jprb
      flux%lw_up(jcol,i_cloud_top) = sum(flux_up(:,1))
      if (config%do_save_spectral_flux) then
        call indexed_sum(flux_up(:,1), &
             &           config%i_spec_from_reordered_g_lw, &
             &           flux%lw_up_band(:,jcol,i_cloud_top))
      end if
      do jlev = i_cloud_top-1,1,-1
        flux_up(:,1) = trans_clear(:,jlev)*flux_up(:,1) + source_up_clear(:,jlev)
        flux%lw_up(jcol,jlev) = sum(flux_up(:,1))
        if (config%do_save_spectral_flux) then
          call indexed_sum(flux_up(:,1), &
               &           config%i_spec_from_reordered_g_lw, &
               &           flux%lw_up_band(:,jcol,jlev))
        end if
      end do
      flux%lw_up_toa_g(:,jcol) = sum(flux_up,2)

      ! --------------------------------------------------------
      ! Section 8: Compute fluxes down to surface
      ! --------------------------------------------------------

      ! Copy over downwelling spectral fluxes at top of first
      ! scattering layer, using overlap matrix to translate to the
      ! regions of the first layer of cloud
      do jreg = 1,nregions
        flux_dn(:,jreg)  = v_matrix(jreg,1,i_cloud_top,jcol)*flux_dn_clear(:,i_cloud_top)
      end do

      ! Final loop back down through the atmosphere to compute fluxes
      do jlev = i_cloud_top,nlev

        if (is_clear_sky_layer(jlev)) then
          do jg = 1,ng
            flux_dn(jg,1) = (transmittance(jg,1,jlev)*flux_dn(jg,1) &
                 &  + reflectance(jg,1,jlev)*total_source(jg,1,jlev+1) + source_dn(jg,1,jlev) ) &
                 &  / (1.0_jprb - reflectance(jg,1,jlev)*total_albedo(jg,1,jlev+1))
            flux_up(jg,1) = total_source(jg,1,jlev+1) + flux_dn(jg,1)*total_albedo(jg,1,jlev+1)
          end do
          flux_dn(:,2:)  = 0.0_jprb
          flux_up(:,2:)  = 0.0_jprb
        else
          flux_dn = (transmittance(:,:,jlev)*flux_dn &
               &     + reflectance(:,:,jlev)*total_source(:,:,jlev+1) + source_dn(:,:,jlev) ) &
               &  / (1.0_jprb - reflectance(:,:,jlev)*total_albedo(:,:,jlev+1))
          flux_up = total_source(:,:,jlev+1) + flux_dn*total_albedo(:,:,jlev+1)
        end if

        if (.not. (is_clear_sky_layer(jlev) &
             &    .and. is_clear_sky_layer(jlev+1))) then
          ! Account for overlap rules in translating fluxes just above
          ! a layer interface to the values just below
          flux_dn_below = singlemat_x_vec(ng,ng,nregions, &
               &  v_matrix(:,:,jlev+1,jcol), flux_dn)
          flux_dn = flux_dn_below
        end if ! Otherwise the fluxes in each region are the same so
               ! nothing to do

        ! Store the broadband fluxes
        flux%lw_up(jcol,jlev+1) = sum(sum(flux_up,1))
        flux%lw_dn(jcol,jlev+1) = sum(sum(flux_dn,1))

        ! Save the spectral fluxes if required
        if (config%do_save_spectral_flux) then
          call indexed_sum(sum(flux_up,2), &
               &           config%i_spec_from_reordered_g_lw, &
               &           flux%lw_up_band(:,jcol,jlev+1))
          call indexed_sum(sum(flux_dn,2), &
               &           config%i_spec_from_reordered_g_lw, &
               &           flux%lw_dn_band(:,jcol,jlev+1))
         end if

      end do ! Final loop over levels
      
      ! Store surface spectral downwelling fluxes, which at this point
      ! are at the surface
      flux%lw_dn_surf_g(:,jcol) = sum(flux_dn,2)

      ! Compute the longwave derivatives needed by Hogan and Bozzo
      ! (2015) approximate radiation update scheme
      if (config%do_lw_derivatives) then
        ! Note that at this point flux_up contains the spectral
        ! fluxes into the regions of the lowest layer; we sum over
        ! regions first to provide a simple spectral flux upwelling
        ! from the surface
        call calc_lw_derivatives_region(ng, nlev, nregions, jcol, transmittance, &
             &  u_matrix(:,:,:,jcol), sum(flux_up,2), flux%lw_derivatives)
      end if
      
    end do ! Loop over columns

    if (lhook) call dr_hook('radiation_tripleclouds_lw:solver_tripleclouds_lw',1,hook_handle)

  end subroutine solver_tripleclouds_lw

end module radiation_tripleclouds_lw
