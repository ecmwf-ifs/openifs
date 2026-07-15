! radiation_adding_ica_sw.F90 - Shortwave adding method in independent column approximation
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
!   2017-10-23  R. Hogan  Renamed single-character variables

module radiation_adding_ica_sw

  public

contains

  subroutine adding_ica_sw(ncol, nlev, incoming_toa, &
       &  albedo_surf_diffuse, albedo_surf_direct, cos_sza, &
       &  reflectance, transmittance, ref_dir, trans_dir_diff, trans_dir_dir, &
       &  flux_up, flux_dn_diffuse, flux_dn_direct)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! Inputs
    integer, intent(in) :: ncol ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! Incoming downwelling solar radiation at top-of-atmosphere (W m-2)
    real(jprb), intent(in),  dimension(ncol)         :: incoming_toa

    ! Surface albedo to diffuse and direct radiation
    real(jprb), intent(in),  dimension(ncol)         :: albedo_surf_diffuse, &
         &                                              albedo_surf_direct

    ! Cosine of the solar zenith angle
    real(jprb), intent(in),  dimension(ncol)         :: cos_sza

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance

    ! Fraction of direct-beam solar radiation entering the top of a
    ! layer that is reflected back up or scattered forward into the
    ! diffuse stream at the base of the layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: ref_dir, trans_dir_diff

    ! Direct transmittance, i.e. fraction of direct beam that
    ! penetrates a layer without being scattered or absorbed
    real(jprb), intent(in),  dimension(ncol, nlev)   :: trans_dir_dir

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling,
    ! diffuse downwelling and direct downwelling
    real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn_diffuse, &
         &                                              flux_dn_direct
    
    ! Albedo of the entire earth/atmosphere system below each half
    ! level
    real(jprb), dimension(ncol, nlev+1) :: albedo

    ! Upwelling radiation at each half-level due to scattering of the
    ! direct beam below that half-level (W m-2)
    real(jprb), dimension(ncol, nlev+1) :: source

    ! Equal to 1/(1-albedo*reflectance)
    real(jprb), dimension(ncol, nlev)   :: inv_denominator

    ! Loop index for model level and column
    integer :: jlev, jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_adding_ica_sw:adding_ica_sw',0,hook_handle)

    ! Compute profile of direct (unscattered) solar fluxes at each
    ! half-level by working down through the atmosphere
    flux_dn_direct(:,1) = incoming_toa
    do jlev = 1,nlev
      flux_dn_direct(:,jlev+1) = flux_dn_direct(:,jlev)*trans_dir_dir(:,jlev)
    end do

    albedo(:,nlev+1) = albedo_surf_diffuse

    ! At the surface, the direct solar beam is reflected back into the
    ! diffuse stream
    source(:,nlev+1) = albedo_surf_direct * flux_dn_direct(:,nlev+1) * cos_sza

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to direct
    ! radiation that is scattered below that level
! Added for DWD (2020)
!NEC$ outerloop_unroll(8)
    do jlev = nlev,1,-1
      ! Next loop over columns. We could do this by indexing the
      ! entire inner dimension as follows, e.g. for the first line:
      !   inv_denominator(:,jlev) = 1.0_jprb / (1.0_jprb-albedo(:,jlev+1)*reflectance(:,jlev))
      ! and similarly for subsequent lines, but this slows down the
      ! routine by a factor of 2!  Rather, we do it with an explicit
      ! loop.
      do jcol = 1,ncol
        ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
        inv_denominator(jcol,jlev) = 1.0_jprb / (1.0_jprb-albedo(jcol,jlev+1)*reflectance(jcol,jlev))
        ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
        albedo(jcol,jlev) = reflectance(jcol,jlev) + transmittance(jcol,jlev) * transmittance(jcol,jlev) &
             &                                     * albedo(jcol,jlev+1) * inv_denominator(jcol,jlev)
        ! Shonk & Hogan (2008) Eq 11:
        source(jcol,jlev) = ref_dir(jcol,jlev)*flux_dn_direct(jcol,jlev) &
             &  + transmittance(jcol,jlev)*(source(jcol,jlev+1) &
             &        + albedo(jcol,jlev+1)*trans_dir_diff(jcol,jlev)*flux_dn_direct(jcol,jlev)) &
             &  * inv_denominator(jcol,jlev)
      end do
    end do

    ! At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn_diffuse(:,1) = 0.0_jprb

    ! At top-of-atmosphere, all upwelling radiation is due to
    ! scattering by the direct beam below that level
    flux_up(:,1) = source(:,1)

    ! Work back down through the atmosphere computing the fluxes at
    ! each half-level
! Added for DWD (2020)
!NEC$ outerloop_unroll(8)
    do jlev = 1,nlev
      do jcol = 1,ncol
        ! Shonk & Hogan (2008) Eq 14 (after simplification):
        flux_dn_diffuse(jcol,jlev+1) &
             &  = (transmittance(jcol,jlev)*flux_dn_diffuse(jcol,jlev) &
             &     + reflectance(jcol,jlev)*source(jcol,jlev+1) &
             &     + trans_dir_diff(jcol,jlev)*flux_dn_direct(jcol,jlev)) * inv_denominator(jcol,jlev)
        ! Shonk & Hogan (2008) Eq 12:
        flux_up(jcol,jlev+1) = albedo(jcol,jlev+1)*flux_dn_diffuse(jcol,jlev+1) &
             &            + source(jcol,jlev+1)
        flux_dn_direct(jcol,jlev) = flux_dn_direct(jcol,jlev)*cos_sza(jcol)
      end do
    end do
    flux_dn_direct(:,nlev+1) = flux_dn_direct(:,nlev+1)*cos_sza

    if (lhook) call dr_hook('radiation_adding_ica_sw:adding_ica_sw',1,hook_handle)

  end subroutine adding_ica_sw

end module radiation_adding_ica_sw
