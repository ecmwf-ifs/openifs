! radiation_overlap.F90 - Module to compute cloud overlap quantities
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
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-10-05  R. Hogan  Generalized alpha overlap for non-equal regions
!   2018-10-08  R. Hogan  Removed calc_region_fractions

module radiation_overlap

  implicit none

  public :: calc_overlap_matrices

contains


  ! This function now superceded by calc_region_properties in module
  ! radiation_regions
  ! !---------------------------------------------------------------------
  ! ! Return an array of length nreg containing the fraction of the
  ! ! gridbox occupied by each region for the specified cloud fraction.
  ! pure function calc_region_fractions(nreg, cloud_fraction)

  !   use parkind1, only : jprb

  !   integer,    intent(in)      :: nreg
  !   real(jprb), intent(in)      :: cloud_fraction

  !   real(jprb), dimension(nreg) :: calc_region_fractions
  !   integer :: jreg

  !   if (nreg == 1) then
  !     ! Only one region: must occupy all of gridbox
  !     calc_region_fractions(1) = 1.0_jprb
  !   else
  !     ! Two or more regions: the first is the cloud-free region
  !     calc_region_fractions(1) = 1.0_jprb - cloud_fraction

  !     do jreg = 2,nreg
  !       ! The cloudy regions are assumed to each have the same
  !       ! fraction - see Shonk and Hogan (2008) for justification
  !       calc_region_fractions(jreg) = cloud_fraction / (nreg - 1.0_jprb)
  !     end do
  !   end if

  ! end function calc_region_fractions

  !---------------------------------------------------------------------
  ! Calculate a matrix expressing the overlap of regions in adjacent
  ! layers, using the method of Shonk et al. (2010) in terms of their
  ! "beta" overlap parameter
  pure function calc_beta_overlap_matrix(nreg, op, frac_upper, frac_lower, &
       &  frac_threshold) result(overlap_matrix)

    use parkind1, only : jprb
    
    integer, intent(in) :: nreg ! Number of regions

    ! Overlap parameter for each region, and fraction of the gridbox
    ! occupied by each region in the upper and lower layers
    real(jprb), intent(in), dimension(nreg) :: op, frac_upper, frac_lower

    ! Cloud-fraction threshold below which cloud is deemed not to be
    ! present
    real(jprb), intent(in) :: frac_threshold

    ! Output overlap matrix
    real(jprb) :: overlap_matrix(nreg,nreg)

    ! Denominator and its reciprocal in computing the random part of
    ! the overlap matrix
    real(jprb) :: denominator, factor

    ! Beta overlap parameter multiplied by the minimum region fraction
    ! of the upper and lower layers
    real(jprb) :: op_x_frac_min(nreg)

    integer :: jupper, jlower, jreg

    ! In computing the random part of the overlap matrix we need
    ! to divide all elements by "denominator", or for efficiency
    ! multiply by "factor"
    denominator = 1.0_jprb
    do jreg = 1,nreg
      op_x_frac_min(jreg) = op(jreg) &
           &  * min(frac_upper(jreg), frac_lower(jreg))
      denominator = denominator - op_x_frac_min(jreg)
    end do
    ! In principle the denominator can be zero
    if (denominator >= frac_threshold) then
      factor = 1.0_jprb / denominator
      ! Create the random part of the overlap matrix
      do jupper = 1,nreg
        do jlower = 1,nreg
          overlap_matrix(jupper,jlower) = factor &
               &  * (frac_lower(jlower)-op_x_frac_min(jlower)) &
               &  * (frac_upper(jupper)-op_x_frac_min(jupper))
        end do
      end do
    else
      overlap_matrix = 0.0_jprb
    end if
    
    ! Add on the maximum part of the overlap matrix
    do jreg = 1,nreg
      overlap_matrix(jreg,jreg) = overlap_matrix(jreg,jreg) &
           &  + op_x_frac_min(jreg)
    end do

  end function calc_beta_overlap_matrix


  !---------------------------------------------------------------------
  ! Calculate a matrix expressing the overlap of regions in adjacent
  ! layers, using the Hogan and Illingworth (2000) "alpha" overlap
  ! parameter, but allowing for the two cloudy regions in the
  ! Tripleclouds assumption to have different areas
  pure function calc_alpha_overlap_matrix(nreg, op, op_inhom, &
       &  frac_upper, frac_lower) result(overlap_matrix)

    use parkind1, only : jprb
    
    integer, intent(in) :: nreg ! Number of regions

    ! Overlap parameter for cloud boundaries and for internal
    ! inhomogeneities
    real(jprb), intent(in) :: op, op_inhom

    ! Fraction of the gridbox occupied by each region in the upper and
    ! lower layers
    real(jprb), intent(in), dimension(nreg) :: frac_upper, frac_lower

    ! Output overlap matrix
    real(jprb) :: overlap_matrix(nreg,nreg)

    ! Combined cloud cover of pair of layers
    real(jprb) :: pair_cloud_cover

    ! Cloud fraction of upper and lower layers
    real(jprb) :: cf_upper, cf_lower

    ! One divided by cloud fraction
    real(jprb) :: one_over_cf

    ! Fraction of domain with cloud in both layers
    real(jprb) :: frac_both

    cf_upper = sum(frac_upper(2:nreg))
    cf_lower = sum(frac_lower(2:nreg))

    pair_cloud_cover = op*max(cf_upper,cf_lower) &
           &  + (1.0_jprb - op) &
           &  * (cf_upper+cf_lower-cf_upper*cf_lower)

    ! Clear in both layers
    overlap_matrix(1,1) = 1.0_jprb - pair_cloud_cover
    if (nreg == 2) then
      ! Clear in upper layer, cloudy in lower layer
      overlap_matrix(1,2) = pair_cloud_cover - cf_upper
      ! Clear in lower layer, cloudy in upper layer
      overlap_matrix(2,1) = pair_cloud_cover - cf_lower
      ! Cloudy in both layers
      overlap_matrix(2,2) = cf_upper + cf_lower - pair_cloud_cover
    else
       ! Clear in upper layer, cloudy in lower layer
      one_over_cf = 1.0_jprb / max(cf_lower, 1.0e-6_jprb)
      overlap_matrix(1,2) = (pair_cloud_cover - cf_upper) &
           &              * frac_lower(2) * one_over_cf
      overlap_matrix(1,3) = (pair_cloud_cover - cf_upper) &
           &              * frac_lower(3) * one_over_cf
      ! Clear in lower layer, cloudy in upper layer
      one_over_cf = 1.0_jprb / max(cf_upper, 1.0e-6_jprb)
      overlap_matrix(2,1) = (pair_cloud_cover - cf_lower) &
           &              * frac_upper(2) * one_over_cf
      overlap_matrix(3,1) = (pair_cloud_cover - cf_lower) &
           &              * frac_upper(3) * one_over_cf
      ! Cloudy in both layers: frac_both is the fraction of the
      ! gridbox with cloud in both layers
      frac_both = cf_upper + cf_lower - pair_cloud_cover
      ! Treat low and high optical-depth regions within frac_both as
      ! one treats clear and cloudy skies in the whole domain;
      ! redefine the following variables treating the high
      ! optical-depth region as the cloud
      cf_upper = frac_upper(3) / max(cf_upper, 1.0e-6_jprb)
      cf_lower = frac_lower(3) / max(cf_lower, 1.0e-6_jprb)
      pair_cloud_cover = op_inhom*max(cf_upper,cf_lower) &
           &  + (1.0_jprb - op_inhom) &
           &  * (cf_upper+cf_lower-cf_upper*cf_lower)
      ! Assign overlaps for this 2x2 section of the 3x3 matrix as for
      ! the 2-region case above, but multiplied by frac_both
      overlap_matrix(2,2) = frac_both * (1.0_jprb - pair_cloud_cover)
      overlap_matrix(2,3) = frac_both * (pair_cloud_cover - cf_upper)
      overlap_matrix(3,2) = frac_both * (pair_cloud_cover - cf_lower)
      overlap_matrix(3,3) = frac_both * (cf_upper+cf_lower-pair_cloud_cover)
    end if

  end function calc_alpha_overlap_matrix


  !---------------------------------------------------------------------
  ! Calculate a matrix expressing the overlap of regions in adjacent
  ! layers, using the Hogan and Illingworth (2000) "alpha" overlap
  ! parameter, and assuming the two cloudy regions in the Tripleclouds
  ! assumption have the same area
  pure function calc_alpha_overlap_matrix_simple(nreg, op, op_inhom, &
       &  cf_upper, cf_lower) result(overlap_matrix)

    use parkind1, only : jprb
    
    integer, intent(in) :: nreg ! Number of regions

    ! Overlap parameter for cloud boundaries and for internal
    ! inhomogeneities
    real(jprb), intent(in) :: op, op_inhom

    ! Cloud fraction in the upper and lower layers
    real(jprb), intent(in) :: cf_upper, cf_lower

    ! Output overlap matrix
    real(jprb) :: overlap_matrix(nreg,nreg)

    ! Combined cloud cover of pair of layers
    real(jprb) :: pair_cloud_cover

    real(jprb) :: cloud_unit

    pair_cloud_cover = op*max(cf_upper,cf_lower) &
           &  + (1.0_jprb - op) &
           &  * (cf_upper+cf_lower-cf_upper*cf_lower)

    ! Clear in both layers
    overlap_matrix(1,1) = 1.0_jprb - pair_cloud_cover
    if (nreg == 2) then
      ! Clear in upper layer, cloudy in lower layer
      overlap_matrix(1,2) = pair_cloud_cover - cf_upper
      ! Clear in lower layer, cloudy in upper layer
      overlap_matrix(2,1) = pair_cloud_cover - cf_lower
      ! Cloudy in both layers
      overlap_matrix(2,2) = cf_upper + cf_lower - pair_cloud_cover
    else
      ! The following assumes that the two cloudy regions are of equal area.
      ! Clear in upper layer, cloudy in lower layer
      overlap_matrix(1,2) = 0.5_jprb * (pair_cloud_cover - cf_upper)
      overlap_matrix(1,3) = overlap_matrix(1,2)
      ! Clear in lower layer, cloudy in upper layer
      overlap_matrix(2,1) = 0.5_jprb * (pair_cloud_cover - cf_lower)
      overlap_matrix(3,1) = overlap_matrix(2,1)
      ! Cloudy in both layers
      cloud_unit = 0.25_jprb * (cf_upper + cf_lower - pair_cloud_cover)
      overlap_matrix(2,2) = cloud_unit * (1.0_jprb + op_inhom)
      overlap_matrix(2,3) = cloud_unit * (1.0_jprb - op_inhom)
      overlap_matrix(3,3) = overlap_matrix(2,2)
      overlap_matrix(3,2) = overlap_matrix(2,3)
    end if

  end function calc_alpha_overlap_matrix_simple


  !---------------------------------------------------------------------
  ! Compute the upward and downward overlap matrices u_matrix and
  ! v_matrix, respectively, where u_matrix is defined such that
  ! y=u_matrix*x, where x is a vector of upwelling fluxes in each
  ! region just below an interface, and y is a vector of upwelling
  ! fluxes in each region just above that interface. For nlev model
  ! levels there are nlev+1 interfaces including the ground and
  ! top-of-atmosphere, and so that is one of the dimensions of
  ! u_matrix and v_matrix.
  subroutine calc_overlap_matrices(nlev,nreg,istartcol,iendcol, &
       &     region_fracs, overlap_param, u_matrix, v_matrix, decorrelation_scaling, &
       &     cloud_fraction_threshold, cloud_cover, use_beta_overlap)

    use parkind1,     only : jprb
    use yomhook,      only : lhook, dr_hook, jphook

    ! Number of levels and regions
    integer,  intent(in) :: nlev, nreg

    ! Range of columns to process (also outer dimensions of u_matrix
    ! and v_matrix)
    integer, intent(in) :: istartcol, iendcol

    ! Area fraction of each region: region 1 is clear sky, and 2+ are
    ! the cloudy regions (only one or two cloudy regions are
    ! supported)
    real(jprb), intent(in), dimension(1:nreg,nlev,istartcol:iendcol)  :: region_fracs

    ! The overlap parameter: either the "alpha" of Hogan & Illingworth
    ! (2000) or the "beta" of Shonk et al. (2010)
    real(jprb), intent(in), dimension(:,:)  :: overlap_param  ! (ncol,nlev-1)

    ! Output overlap matrices
    real(jprb), intent(out), dimension(nreg,nreg,nlev+1,istartcol:iendcol) &
         &  :: u_matrix, v_matrix

    ! For regions 2 and above, the overlap decorrelation length for
    ! cloud boundaries is scaled by this amount to obtain the overlap
    ! decorrelation length for cloud inhomogeneities. Typically this
    ! number is 0.5, but if omitted it will be assumed to be one (same
    ! decorrelation for cloud boundaries and in-cloud inhomogeneities)
    real(jprb), intent(in), optional :: decorrelation_scaling

    ! Regions smaller than this are ignored
    real(jprb), intent(in), optional :: cloud_fraction_threshold

    ! The diagnosed cloud cover is an optional output
    real(jprb), intent(out), optional :: cloud_cover(:)

    ! Do we use Shonk et al.'s (2010) "beta" overlap parameter?
    logical, intent(in), optional :: use_beta_overlap

    ! Loop indices for column, level, region and the regions in the
    ! upper and lower layers for an interface
    integer  :: jcol, jlev, jupper, jlower

    ! Overlap matrix (non-directional)
    real(jprb) :: overlap_matrix(nreg,nreg)

    ! Fraction of the gridbox occupied by each region in the upper and
    ! lower layers for an interface
    real(jprb) :: frac_upper(nreg), frac_lower(nreg)

    ! Beta overlap parameter for each region
    real(jprb) :: op(nreg)

    ! In case the user doesn't supply cloud_fraction_threshold we use
    ! a default value
    real(jprb) :: frac_threshold

    ! The decorrelation scaling to use, in case decorrelation_scaling
    ! was not provided
    real(jprb) :: used_decorrelation_scaling

    logical :: use_beta_overlap_param

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_overlap:calc_overlap_matrices',0,hook_handle)

    if (present(decorrelation_scaling)) then
       used_decorrelation_scaling = decorrelation_scaling
    else
       used_decorrelation_scaling = 1.0_jprb
    end if

    if (present(cloud_fraction_threshold)) then
      frac_threshold = cloud_fraction_threshold
    else
      frac_threshold = 1.0e-20_jprb
    end if

    if (present(use_beta_overlap)) then
      use_beta_overlap_param = use_beta_overlap
    else
      use_beta_overlap_param = .false.
    end if

    ! Loop through each atmospheric column
    do jcol = istartcol, iendcol
      ! For this column, outer space is treated as one clear-sky
      ! region, so the fractions are assigned as such
      frac_upper(1) = 1.0_jprb
      frac_upper(2:nreg) = 0.0_jprb

      ! Overlap parameter is irrelevant when there is only one region
      ! in the upper layer
      op = 1.0_jprb

      ! Loop down through the atmosphere, where jlev indexes each
      ! half-level starting at 1 for the top-of-atmosphere, as well
      ! as indexing each level starting at 1 for the top-most level.
      do jlev = 1,nlev+1
        ! Fraction of each region just below the interface
        if (jlev > nlev) then
          ! We are at the surface: treat as a single clear-sky
          ! region
          frac_lower(1) = 1.0_jprb
          frac_lower(2:nreg) = 0.0_jprb
        else
          frac_lower = region_fracs(1:nreg,jlev,jcol)
        end if
   
        ! Compute the overlap parameter of the interface just below
        ! the current full level
        if (jlev == 1 .or. jlev > nlev) then
          ! We are at the surface or top-of-atmosphere: overlap
          ! parameter is irrelevant
          op = 1.0_jprb
        else
          ! We are not at the surface
          op(1) = overlap_param(jcol,jlev-1)
          ! For cloudy regions, scale the cloud-boundary overlap
          ! parameter to obtain the cloud-inhomogeneity overlap
          ! parameter as follows
          if (op(1) >= 0.0_jprb) then
            op(2:nreg) = op(1)**(1.0_jprb/used_decorrelation_scaling)
          else
            op(2:nreg) = op(1)
          end if
        end if
     
        if (use_beta_overlap_param) then
          overlap_matrix = calc_beta_overlap_matrix(nreg, op, &
               &  frac_upper, frac_lower, frac_threshold)
        else
          ! Simpler scheme assuming the two cloudy regions have the
          ! same fraction
          !overlap_matrix = calc_alpha_overlap_matrix_simple(nreg, &
          !     &  op(1), op(2), &
          !     &  1.0_jprb - frac_upper(1), 1.0_jprb - frac_lower(1))
          ! More general scheme
          overlap_matrix = calc_alpha_overlap_matrix(nreg, &
               &  op(1), op(2), frac_upper, frac_lower)
        end if

        ! Convert to directional overlap matrices
        do jupper = 1,nreg
          do jlower = 1,nreg
            if (frac_lower(jlower) >= frac_threshold) then
              u_matrix(jupper,jlower,jlev,jcol) = overlap_matrix(jupper,jlower) &
                   &  / frac_lower(jlower)
            else
              u_matrix(jupper,jlower,jlev,jcol) = 0.0_jprb
            end if
            if (frac_upper(jupper) >= frac_threshold) then
              v_matrix(jlower,jupper,jlev,jcol) = overlap_matrix(jupper,jlower) &
                   &  / frac_upper(jupper)
            else
              v_matrix(jlower,jupper,jlev,jcol) = 0.0_jprb
            end if
          end do
        end do
        frac_upper = frac_lower
        
      end do ! levels

      ! Compute cloud cover from one of the directional overlap matrices
      if (present(cloud_cover)) then
        cloud_cover(jcol) = 1.0_jprb - product(v_matrix(1,1,:,jcol))
      end if

    end do ! columns

    if (lhook) call dr_hook('radiation_overlap:calc_overlap_matrices',1,hook_handle)

  end subroutine calc_overlap_matrices
  
end module radiation_overlap
