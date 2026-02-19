! radiation_tripleclouds_sw.F90 - Shortwave "Tripleclouds" solver
!
! Copyright (C) 2016-2019 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! Modifications
!   2017-04-11  R. Hogan  Receive albedos at g-points
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-10-08  R. Hogan  Call calc_region_properties
!   2019-01-02  R. Hogan  Fixed problem of do_save_spectral_flux .and. .not. do_sw_direct

module radiation_tripleclouds_sw

contains
  ! Provides elemental function "delta_eddington"
#include "radiation_delta_eddington.h"

  ! Small routine for scaling cloud optical depth in the cloudy
  ! regions
#include "radiation_optical_depth_scaling.h"

  ! This module contains just one subroutine, the shortwave
  ! "Tripleclouds" solver in which cloud inhomogeneity is treated by
  ! dividing each model level into three regions, one clear and two
  ! cloudy (with differing optical depth). This approach was described
  ! by Shonk and Hogan (2008).

  subroutine solver_tripleclouds_sw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
       &  albedo_direct, albedo_diffuse, incoming_sw, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook
!    use radiation_io, only             : nulout
    use radiation_config, only         : config_type, IPdfShapeGamma
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_regions, only        : calc_region_properties
    use radiation_overlap, only        : calc_overlap_matrices
    use radiation_flux, only           : flux_type, &
         &                               indexed_sum, add_indexed_sum
    use radiation_matrix, only         : singlemat_x_vec
    use radiation_two_stream, only     : calc_two_stream_gammas_sw, &
         &                       calc_reflectance_transmittance_sw

    implicit none

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each shortwave g-point
!    real(jprb), intent(in), dimension(istartcol:iendcol,nlev,config%n_g_sw) :: &
    real(jprb), intent(in), dimension(config%n_g_sw,nlev,istartcol:iendcol) :: &
         &  od, ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each shortwave band
    real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol) :: &
         &  od_cloud, ssa_cloud, g_cloud

    ! Optical depth, single scattering albedo and asymmetry factor in
    ! each g-point including gas, aerosol and clouds
    real(jprb), dimension(config%n_g_sw) :: od_total, ssa_total, g_total

    ! Direct and diffuse surface albedos, and the incoming shortwave
    ! flux into a plane perpendicular to the incoming radiation at
    ! top-of-atmosphere in each of the shortwave g points
    real(jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol) :: &
         &  albedo_direct, albedo_diffuse, incoming_sw

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

    ! Two-stream variables
    real(jprb), dimension(config%n_g_sw) :: gamma1, gamma2, gamma3

    ! Diffuse reflection and transmission matrices of each layer
    real(jprb), dimension(config%n_g_sw, nregions, nlev) &
         &  :: reflectance, transmittance

    ! Terms translating the direct flux entering the layer from above
    ! to the reflected radiation exiting upwards (ref_dir) and the
    ! scattered radiation exiting downwards (trans_dir_diff), along with the
    ! direct unscattered transmission matrix (trans_dir_dir).
    real(jprb), dimension(config%n_g_sw, nregions, nlev) &
         &  :: ref_dir, trans_dir_diff, trans_dir_dir

    ! Total albedo of the atmosphere/surface just above a layer
    ! interface with respect to downwelling diffuse and direct
    ! (respectively) radiation at that interface, where level index =
    ! 1 corresponds to the top-of-atmosphere
    real(jprb), dimension(config%n_g_sw, nregions, nlev+1) &
         &  :: total_albedo, total_albedo_direct

    ! ...equivalent values for clear-skies
    real(jprb), dimension(config%n_g_sw, nlev+1) &
         &  :: total_albedo_clear, total_albedo_clear_direct

    ! Total albedo of the atmosphere just below a layer interface
    real(jprb), dimension(config%n_g_sw, nregions) &
         &  :: total_albedo_below, total_albedo_below_direct

    ! Direct downwelling flux below and above an interface between
    ! layers into a plane perpendicular to the direction of the sun
    real(jprb), dimension(config%n_g_sw, nregions) &
         &  :: direct_dn, direct_dn_below
    ! Diffuse equivalents
    real(jprb), dimension(config%n_g_sw, nregions) &
         &  :: flux_dn, flux_dn_below, flux_up

    ! ...clear-sky equivalent (no distinction between "above/below")
    real(jprb), dimension(config%n_g_sw) &
         &  :: direct_dn_clear, flux_dn_clear, flux_up_clear

    ! Clear-sky equivalent, but actually its reciprocal to replace
    ! some divisions by multiplications
    real(jprb), dimension(config%n_g_sw, nregions) :: inv_denom

    ! Identify clear-sky layers, with pseudo layers for outer space
    ! and below the ground, both treated as single-region clear skies
    logical :: is_clear_sky_layer(0:nlev+1)

    ! Scattering optical depth of gas+aerosol and of cloud
    real(jprb) :: scat_od, scat_od_cloud

    real(jprb) :: mu0

    integer :: jcol, jlev, jg, jreg, iband, jreg2, ng

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_tripleclouds_sw:solver_tripleclouds_sw',0,hook_handle)

    ! --------------------------------------------------------
    ! Section 1: Prepare general variables and arrays
    ! --------------------------------------------------------
    ! Copy array dimensions to local variables for convenience
    ng   = config%n_g_sw

    ! Compute the wavelength-independent region fractions and
    ! optical-depth scalings
    call calc_region_properties(nlev,nregions,istartcol,iendcol, &
         &  config%i_cloud_pdf_shape == IPdfShapeGamma, &
         &  cloud%fraction, cloud%fractional_std, region_fracs, &
         &  od_scaling, config%cloud_fraction_threshold)

    ! Compute wavelength-independent overlap matrices u_matrix and
    ! v_matrix
    call calc_overlap_matrices(nlev,nregions,istartcol,iendcol, &
         &  region_fracs, cloud%overlap_param, &
         &  u_matrix, v_matrix, &
         &  decorrelation_scaling=config%cloud_inhom_decorr_scaling, &
         &  cloud_fraction_threshold=config%cloud_fraction_threshold, &
         &  use_beta_overlap=config%use_beta_overlap, &
         &  cloud_cover=flux%cloud_cover_sw)

    ! Main loop over columns
    do jcol = istartcol, iendcol
      ! --------------------------------------------------------
      ! Section 2: Prepare column-specific variables and arrays
      ! --------------------------------------------------------

      ! Copy local cosine of the solar zenith angle
      mu0 = single_level%cos_sza(jcol)

      ! Skip profile if sun is too low in the sky
      if (mu0 < 1.0e-10_jprb) then
        flux%sw_dn(jcol,:) = 0.0_jprb
        flux%sw_up(jcol,:) = 0.0_jprb
        if (allocated(flux%sw_dn_direct)) then
          flux%sw_dn_direct(jcol,:) = 0.0_jprb
        end if
        if (config%do_clear) then
          flux%sw_dn_clear(jcol,:) = 0.0_jprb
          flux%sw_up_clear(jcol,:) = 0.0_jprb
          if (allocated(flux%sw_dn_direct_clear)) then
            flux%sw_dn_direct_clear(jcol,:) = 0.0_jprb
          end if
        end if

        if (config%do_save_spectral_flux) then
          flux%sw_dn_band(:,jcol,:) = 0.0_jprb
          flux%sw_up_band(:,jcol,:) = 0.0_jprb
          if (allocated(flux%sw_dn_direct_band)) then
            flux%sw_dn_direct_band(:,jcol,:) = 0.0_jprb
          end if
          if (config%do_clear) then
            flux%sw_dn_clear_band(:,jcol,:) = 0.0_jprb
            flux%sw_up_clear_band(:,jcol,:) = 0.0_jprb
            if (allocated(flux%sw_dn_direct_clear_band)) then
              flux%sw_dn_direct_clear_band(:,jcol,:) = 0.0_jprb
            end if
          end if
        end if

        flux%sw_dn_diffuse_surf_g(:,jcol) = 0.0_jprb
        flux%sw_dn_direct_surf_g(:,jcol)  = 0.0_jprb
        if (config%do_clear) then
          flux%sw_dn_diffuse_surf_clear_g(:,jcol) = 0.0_jprb
          flux%sw_dn_direct_surf_clear_g(:,jcol)  = 0.0_jprb
        end if

        cycle
      end if ! sun is below the horizon

      ! At this point mu0 >= 1.0e-10

      ! Define which layers contain cloud; assume that
      ! cloud%crop_cloud_fraction has already been called
      is_clear_sky_layer = .true.
      do jlev = 1,nlev
        if (cloud%fraction(jcol,jlev) > 0.0_jprb) then
          is_clear_sky_layer(jlev) = .false.
        end if
      end do

      ! --------------------------------------------------------
      ! Section 3: Loop over layers to compute reflectance and transmittance
      ! --------------------------------------------------------
      ! In this section the reflectance, transmittance and sources
      ! are computed for each layer
      do jlev = 1,nlev ! Start at top-of-atmosphere

        ! Array-wise assignments
        gamma1 = 0.0_jprb
        gamma2 = 0.0_jprb
        gamma3 = 0.0_jprb

        nreg = nregions
        if (is_clear_sky_layer(jlev)) then
          nreg = 1
        end if
        do jreg = 1,nreg
          if (jreg == 1) then
            od_total  =  od(:,jlev,jcol)
            ssa_total = ssa(:,jlev,jcol)
            g_total   =   g(:,jlev,jcol)
          else
            do jg = 1,ng
              ! Mapping from g-point to band
              iband = config%i_band_from_reordered_g_sw(jg)
              scat_od = od(jg,jlev,jcol)*ssa(jg,jlev,jcol)
              scat_od_cloud = od_cloud(iband,jlev,jcol) &
                   &  * ssa_cloud(iband,jlev,jcol) * od_scaling(jreg,jlev,jcol)
              ! Add scaled cloud optical depth to clear-sky value
              od_total(jg) = od(jg,jlev,jcol) &
                   &  + od_cloud(iband,jlev,jcol)*od_scaling(jreg,jlev,jcol)
              ! Compute single-scattering albedo and asymmetry
              ! factor of gas-cloud combination
              ssa_total(jg) = (scat_od+scat_od_cloud) &
                   &  / od_total(jg)
              g_total(jg) = (scat_od*g(jg,jlev,jcol) &
                   &         + scat_od_cloud * g_cloud(iband,jlev,jcol)) &
                   &      / (scat_od + scat_od_cloud)
            end do
          end if

          if (config%do_sw_delta_scaling_with_gases) then
            ! Apply delta-Eddington scaling to the aerosol-gas(-cloud)
            ! mixture
            call delta_eddington(od_total, ssa_total, g_total)
          end if
          call calc_two_stream_gammas_sw(ng, &
               &  mu0, ssa_total, g_total, &
               &  gamma1, gamma2, gamma3)
          call calc_reflectance_transmittance_sw(ng, &
               &  mu0, od_total, ssa_total, &
               &  gamma1, gamma2, gamma3, &
               &  reflectance(:,jreg,jlev), transmittance(:,jreg,jlev), &
               &  ref_dir(:,jreg,jlev), trans_dir_diff(:,jreg,jlev), &
               &  trans_dir_dir(:,jreg,jlev) )
        end do
      end do

      ! --------------------------------------------------------
      ! Section 4: Compute total albedos
      ! --------------------------------------------------------

      total_albedo(:,:,:) = 0.0_jprb
      total_albedo_direct(:,:,:) = 0.0_jprb

      ! Copy surface albedo in clear-sky region
      do jg = 1,ng
        total_albedo(jg,1,nlev+1) = albedo_diffuse(jg,jcol)
      end do

      ! If direct albedo is available, use it; otherwise copy from
      ! diffuse albedo
      do jg = 1,ng
        total_albedo_direct(jg,1,nlev+1) &
             &  = mu0 * albedo_direct(jg,jcol)
      end do

      ! If there is cloud in the lowest layer then we need the albedos
      ! underneath
      if (.not. is_clear_sky_layer(nlev)) then
        do jreg = 2,nregions
          total_albedo(:,jreg,nlev+1)        = total_albedo(:,1,nlev+1)
          total_albedo_direct(:,jreg,nlev+1) = total_albedo_direct(:,1,nlev+1)
        end do
      end if

      if (config%do_clear) then
        total_albedo_clear(:,nlev+1)        = total_albedo(:,1,nlev+1)
        total_albedo_clear_direct(:,nlev+1) = total_albedo_direct(:,1,nlev+1)
      end if

      ! Work up from the surface computing the total albedo of the
      ! atmosphere below that point using the adding method
      do jlev = nlev,1,-1

        total_albedo_below        = 0.0_jprb
        total_albedo_below_direct = 0.0_jprb

        if (config%do_clear) then
          ! For clear-skies there is no need to consider "above" and
          ! "below" quantities since with no cloud overlap to worry
          ! about, these are the same
          inv_denom(:,1) = 1.0_jprb &
               &  / (1.0_jprb - total_albedo_clear(:,jlev+1)*reflectance(:,1,jlev))
          total_albedo_clear(:,jlev) = reflectance(:,1,jlev) &
               &  + transmittance(:,1,jlev) * transmittance(:,1,jlev) &
               &  * total_albedo_clear(:,jlev+1) * inv_denom(:,1)
          total_albedo_clear_direct(:,jlev) = ref_dir(:,1,jlev) &
               &  + (trans_dir_dir(:,1,jlev)*total_albedo_clear_direct(:,jlev+1) &
               &     +trans_dir_diff(:,1,jlev)*total_albedo_clear(:,jlev+1)) &
               &  * transmittance(:,1,jlev) * inv_denom(:,1)
        end if

        if (is_clear_sky_layer(jlev)) then
          inv_denom(:,1) = 1.0_jprb &
               &  / (1.0_jprb - total_albedo(:,1,jlev+1)*reflectance(:,1,jlev))
          total_albedo_below(:,1) = reflectance(:,1,jlev) &
               &  + transmittance(:,1,jlev)  * transmittance(:,1,jlev) &
               &  * total_albedo(:,1,jlev+1) * inv_denom(:,1)
          total_albedo_below_direct(:,1) = ref_dir(:,1,jlev) &
               &  + (trans_dir_dir(:,1,jlev)*total_albedo_direct(:,1,jlev+1) &
               &     +trans_dir_diff(:,1,jlev)*total_albedo(:,1,jlev+1)) &
               &  * transmittance(:,1,jlev) * inv_denom(:,1)
        else
          inv_denom = 1.0_jprb / (1.0_jprb - total_albedo(:,:,jlev+1)*reflectance(:,:,jlev))
          total_albedo_below = reflectance(:,:,jlev) &
               &  + transmittance(:,:,jlev)  * transmittance(:,:,jlev) &
               &  * total_albedo(:,:,jlev+1) * inv_denom
          total_albedo_below_direct = ref_dir(:,:,jlev) &
               &  + (trans_dir_dir(:,:,jlev)*total_albedo_direct(:,:,jlev+1) &
               &     +trans_dir_diff(:,:,jlev)*total_albedo(:,:,jlev+1)) &
               &  * transmittance(:,:,jlev) * inv_denom
        end if

        ! Account for cloud overlap when converting albedo below a
        ! layer interface to the equivalent values just above
        if (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev-1)) then
          total_albedo(:,:,jlev)        = total_albedo_below(:,:)
          total_albedo_direct(:,:,jlev) = total_albedo_below_direct(:,:)
        else
          ! Use overlap matrix and exclude "anomalous" horizontal
          ! transport described by Shonk & Hogan (2008).  Therefore,
          ! the operation we perform is essentially diag(total_albedo)
          ! = matmul(transpose(v_matrix)), diag(total_albedo_below)).
          do jreg = 1,nregions
            do jreg2 = 1,nregions
              total_albedo(:,jreg,jlev) &
                   &  = total_albedo(:,jreg,jlev) &
                   &  + total_albedo_below(:,jreg2) &
                   &  * v_matrix(jreg2,jreg,jlev,jcol)
              total_albedo_direct(:,jreg,jlev) &
                   &  = total_albedo_direct(:,jreg,jlev) &
                   &  + total_albedo_below_direct(:,jreg2) &
                   &  * v_matrix(jreg2,jreg,jlev,jcol)
            end do
          end do

        end if

      end do ! Reverse loop over levels

      ! --------------------------------------------------------
      ! Section 5: Compute fluxes
      ! --------------------------------------------------------

      ! Top-of-atmosphere fluxes into the regions of the top-most
      ! layer, zero since we assume no diffuse downwelling
      flux_dn = 0.0_jprb
      ! Direct downwelling flux (into a plane perpendicular to the
      ! sun) entering the top of each region in the top-most layer
      do jreg = 1,nregions
        direct_dn(:,jreg) = incoming_sw(:,jcol)*region_fracs(jreg,1,jcol)
      end do
      flux_up = direct_dn*total_albedo_direct(:,:,1)

      if (config%do_clear) then
        flux_dn_clear = 0.0_jprb
        direct_dn_clear(:) = incoming_sw(:,jcol)
        flux_up_clear = direct_dn_clear*total_albedo_clear_direct(:,1)
      end if

      ! Store the TOA broadband fluxes
      flux%sw_up(jcol,1) = sum(sum(flux_up,1))
      flux%sw_dn(jcol,1) = mu0 * sum(sum(direct_dn,1))
      if (allocated(flux%sw_dn_direct)) then
        flux%sw_dn_direct(jcol,1) = flux%sw_dn(jcol,1)
      end if
      if (config%do_clear) then
        flux%sw_up_clear(jcol,1) = sum(flux_up_clear)
        flux%sw_dn_clear(jcol,1) = mu0 * sum(direct_dn_clear)
        if (allocated(flux%sw_dn_direct_clear)) then
          flux%sw_dn_direct_clear(jcol,1) = flux%sw_dn_clear(jcol,1)
        end if
      end if

      ! Save the spectral fluxes if required; some redundancy here as
      ! the TOA downwelling flux is the same in clear and cloudy skies
      if (config%do_save_spectral_flux) then
        call indexed_sum(sum(flux_up,2), &
             &           config%i_spec_from_reordered_g_sw, &
             &           flux%sw_up_band(:,jcol,1))
        call indexed_sum(sum(direct_dn,2), &
             &           config%i_spec_from_reordered_g_sw, &
             &           flux%sw_dn_band(:,jcol,1))
        flux%sw_dn_band(:,jcol,1) = &
             &  mu0 * flux%sw_dn_band(:,jcol,1)
        if (allocated(flux%sw_dn_direct_band)) then
          flux%sw_dn_direct_band(:,jcol,1) = flux%sw_dn_band(:,jcol,1)
        end if
        call add_indexed_sum(sum(flux_dn,2), &
             &           config%i_spec_from_reordered_g_sw, &
             &           flux%sw_dn_band(:,jcol,1))
        if (config%do_clear) then
          call indexed_sum(flux_up_clear, &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_up_clear_band(:,jcol,1))
          call indexed_sum(direct_dn_clear, &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_dn_clear_band(:,jcol,1))
          flux%sw_dn_clear_band(:,jcol,1) &
               &  = mu0 * flux%sw_dn_clear_band(:,jcol,1)
          if (allocated(flux%sw_dn_direct_clear_band)) then
            flux%sw_dn_direct_clear_band(:,jcol,1) &
                 &  = flux%sw_dn_clear_band(:,jcol,1)
          end if
          call add_indexed_sum(flux_dn_clear, &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_dn_clear_band(:,jcol,1))
        end if
      end if

      ! Final loop back down through the atmosphere to compute fluxes
      do jlev = 1,nlev
        if (config%do_clear) then
          flux_dn_clear = (transmittance(:,1,jlev)*flux_dn_clear + direct_dn_clear &
               &  * (trans_dir_dir(:,1,jlev)*total_albedo_clear_direct(:,jlev+1)*reflectance(:,1,jlev) &
               &     + trans_dir_diff(:,1,jlev) )) &
               &  / (1.0_jprb - reflectance(:,1,jlev)*total_albedo_clear(:,jlev+1))
          direct_dn_clear = trans_dir_dir(:,1,jlev)*direct_dn_clear
          flux_up_clear = direct_dn_clear*total_albedo_clear_direct(:,jlev+1) &
               &        +   flux_dn_clear*total_albedo_clear(:,jlev+1)
        end if

        if (is_clear_sky_layer(jlev)) then
          flux_dn(:,1) = (transmittance(:,1,jlev)*flux_dn(:,1) + direct_dn(:,1) &
               &  * (trans_dir_dir(:,1,jlev)*total_albedo_direct(:,1,jlev+1)*reflectance(:,1,jlev) &
               &     + trans_dir_diff(:,1,jlev) )) &
               &  / (1.0_jprb - reflectance(:,1,jlev)*total_albedo(:,1,jlev+1))
          direct_dn(:,1) = trans_dir_dir(:,1,jlev)*direct_dn(:,1)
          flux_up(:,1) = direct_dn(:,1)*total_albedo_direct(:,1,jlev+1) &
               &  +        flux_dn(:,1)*total_albedo(:,1,jlev+1)
          flux_dn(:,2:)  = 0.0_jprb
          flux_up(:,2:)  = 0.0_jprb
          direct_dn(:,2:)= 0.0_jprb
        else
          flux_dn = (transmittance(:,:,jlev)*flux_dn + direct_dn &
               &  * (trans_dir_dir(:,:,jlev)*total_albedo_direct(:,:,jlev+1)*reflectance(:,:,jlev) &
               &     + trans_dir_diff(:,:,jlev) )) &
               &  / (1.0_jprb - reflectance(:,:,jlev)*total_albedo(:,:,jlev+1))
          direct_dn = trans_dir_dir(:,:,jlev)*direct_dn
          flux_up = direct_dn*total_albedo_direct(:,:,jlev+1) &
               &  +   flux_dn*total_albedo(:,:,jlev+1)
        end if

        if (.not. (is_clear_sky_layer(jlev) &
             &    .and. is_clear_sky_layer(jlev+1))) then
          ! Account for overlap rules in translating fluxes just above
          ! a layer interface to the values just below
          flux_dn_below = singlemat_x_vec(ng,ng,nregions, &
               &  v_matrix(:,:,jlev+1,jcol), flux_dn)
          direct_dn_below = singlemat_x_vec(ng,ng,nregions, &
               &  v_matrix(:,:,jlev+1,jcol), direct_dn)
          flux_dn = flux_dn_below
          direct_dn = direct_dn_below
        end if ! Otherwise the fluxes in each region are the same so
               ! nothing to do

        ! Store the broadband fluxes
        flux%sw_up(jcol,jlev+1) = sum(sum(flux_up,1))
        if (allocated(flux%sw_dn_direct)) then
          flux%sw_dn_direct(jcol,jlev+1) = mu0 * sum(sum(direct_dn,1))
          flux%sw_dn(jcol,jlev+1) &
               &  = flux%sw_dn_direct(jcol,jlev+1) + sum(sum(flux_dn,1))
        else
          flux%sw_dn(jcol,jlev+1) = mu0 * sum(sum(direct_dn,1)) + sum(sum(flux_dn,1))   
        end if
        if (config%do_clear) then
          flux%sw_up_clear(jcol,jlev+1) = sum(flux_up_clear)
          if (allocated(flux%sw_dn_direct_clear)) then
            flux%sw_dn_direct_clear(jcol,jlev+1) = mu0 * sum(direct_dn_clear)
            flux%sw_dn_clear(jcol,jlev+1) &
                 &  = flux%sw_dn_direct_clear(jcol,jlev+1) + sum(flux_dn_clear)
          else
            flux%sw_dn_clear(jcol,jlev+1) = mu0 * sum(direct_dn_clear) &
                 &  + sum(flux_dn_clear)
          end if
        end if

        ! Save the spectral fluxes if required
        if (config%do_save_spectral_flux) then
          call indexed_sum(sum(flux_up,2), &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_up_band(:,jcol,jlev+1))
          call indexed_sum(sum(direct_dn,2), &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_dn_band(:,jcol,jlev+1))
          flux%sw_dn_band(:,jcol,jlev+1) = &
               &  mu0 * flux%sw_dn_band(:,jcol,jlev+1)
          if (allocated(flux%sw_dn_direct_band)) then
            flux%sw_dn_direct_band(:,jcol,jlev+1) &
                 &  = flux%sw_dn_band(:,jcol,jlev+1)
          end if
          call add_indexed_sum(sum(flux_dn,2), &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_dn_band(:,jcol,jlev+1))
          if (config%do_clear) then
            call indexed_sum(flux_up_clear, &
                 &           config%i_spec_from_reordered_g_sw, &
                 &           flux%sw_up_clear_band(:,jcol,jlev+1))
            call indexed_sum(direct_dn_clear, &
                 &           config%i_spec_from_reordered_g_sw, &
                 &           flux%sw_dn_clear_band(:,jcol,jlev+1))
            flux%sw_dn_clear_band(:,jcol,jlev+1) = &
                 &  mu0 * flux%sw_dn_clear_band(:,jcol,jlev+1)
            if (allocated(flux%sw_dn_direct_clear_band)) then
              flux%sw_dn_direct_clear_band(:,jcol,jlev+1) &
                   &  = flux%sw_dn_clear_band(:,jcol,jlev+1)
            end if
            call add_indexed_sum(flux_dn_clear, &
                 &           config%i_spec_from_reordered_g_sw, &
                 &           flux%sw_dn_clear_band(:,jcol,jlev+1))
          end if
        end if

      end do ! Final loop over levels
      
      ! Store surface spectral fluxes, if required (after the end of
      ! the final loop over levels, the current values of these arrays
      ! will be the surface values)
      flux%sw_dn_diffuse_surf_g(:,jcol) = sum(flux_dn,2)
      flux%sw_dn_direct_surf_g(:,jcol)  = mu0 * sum(direct_dn,2)
      if (config%do_clear) then
        flux%sw_dn_diffuse_surf_clear_g(:,jcol) = flux_dn_clear
        flux%sw_dn_direct_surf_clear_g(:,jcol)  = mu0 * direct_dn_clear
      end if

    end do ! Loop over columns

    if (lhook) call dr_hook('radiation_tripleclouds_sw:solver_tripleclouds_sw',1,hook_handle)

  end subroutine solver_tripleclouds_sw

end module radiation_tripleclouds_sw
