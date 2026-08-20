! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! radiation_cloud_cover.F90 - Compute cumulative cloud cover for McICA
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! Generate profiles of the cumulative cloud cover as seen from TOA,
! used in the McICA cloud generator.
!
! Modifications
!   2020-10-07  R. Hogan  Ensure iobj1 initialized in case of alpha_obj==0

module radiation_cloud_cover

  use parkind1, only           : jprb

  public

  ! Three overlap schemes.  Note that "Exponential" means that
  ! clear-sky regions have no special significance for computing the
  ! cumulative cloud cover: non-contiguous clouds are exponentially
  ! rather than randomly overlapped. This is the situaition in the
  ! McRad radiation scheme at ECMWF.
  enum, bind(c)
    enumerator IOverlapMaximumRandom, IOverlapExponentialRandom, &
         &     IOverlapExponential
  end enum
  character(len=*), parameter :: OverlapName(0:2) = (/ 'Max-Ran', &
       &                                               'Exp-Ran', &
       &                                               'Exp-Exp' /)

  ! Maximum cloud fraction distinguishable from 1
  real(jprb), parameter :: MaxCloudFrac = 1.0_jprb-epsilon(1.0_jprb)*10.0_jprb


contains

  !---------------------------------------------------------------------
  ! Convert "beta" overlap parameter of Shonk et al. (2010) to "alpha"
  ! overlap parameter of Hogan and Illingworth (2000)
  elemental function beta2alpha(beta, frac1, frac2)

    implicit none

    ! Beta overlap parameter and the cloud fractions in the upper and
    ! lower layers
    real(jprb), intent(in) :: beta, frac1, frac2

    real(jprb)             :: beta2alpha

    ! Absolute difference in cloud fraction
    real(jprb)             :: frac_diff

    if (beta < 1.0_jprb) then
      frac_diff = abs(frac1-frac2)
      beta2alpha = beta &
           &  + (1.0_jprb-beta)*frac_diff &
           &  / (frac_diff + 1.0_jprb/beta - 1.0_jprb)
    else
      beta2alpha = 1.0_jprb
    end if

  end function beta2alpha


  !---------------------------------------------------------------------
  ! Compute total cloud cover according to the specified overlap
  ! rule. This can be used to compute the high, mid and low cloud
  ! cover by passing in subsets of the cloud fraction array
  function cloud_cover(nlev, i_overlap_scheme, frac, overlap_param, &
       &               is_beta_overlap)

    implicit none

    ! Number of levels and the overlap scheme to be applied
    integer, intent(in)    :: nlev, i_overlap_scheme

    ! Cloud fraction and the overlap parameter between adjacent pairs
    ! of levels
    real(jprb), intent(in) :: frac(nlev), overlap_param(nlev-1)

    ! Do we use the "beta" overlap scheme of Shonk et al. (2010)?
    ! Default is false.
    logical, intent(in), optional :: is_beta_overlap

    ! Return cloud cover
    real(jprb)             :: cloud_cover

    ! Cumulative cloud cover from TOA to the base of each layer
    real(jprb) :: cum_cloud_cover(nlev)

    ! Cloud cover of a pair of layers
    real(jprb) :: pair_cloud_cover(nlev-1)

    if (i_overlap_scheme == IOverlapExponentialRandom) then
      call cum_cloud_cover_exp_ran(nlev, frac, overlap_param, &
           &   cum_cloud_cover, pair_cloud_cover, is_beta_overlap)
    else if (i_overlap_scheme == IOverlapExponential) then
      call cum_cloud_cover_exp_exp(nlev, frac, overlap_param, &
           &   cum_cloud_cover, pair_cloud_cover, is_beta_overlap)
    else
      call cum_cloud_cover_max_ran(nlev, frac, cum_cloud_cover, &
           &                       pair_cloud_cover)
    end if

    cloud_cover = cum_cloud_cover(nlev)

  end function cloud_cover


  !---------------------------------------------------------------------
  ! Maximum-random overlap: Geleyn & Hollingsworth formula
  subroutine cum_cloud_cover_max_ran(nlev, frac, &
       & cum_cloud_cover, pair_cloud_cover)

    use yomhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! Inputs
    integer, intent(in)     :: nlev  ! number of model levels

    ! Cloud fraction on full levels
    real(jprb), intent(in)  :: frac(nlev)

    ! Outputs

    ! Cumulative cloud cover from TOA to the base of each layer
    real(jprb), intent(out) :: cum_cloud_cover(nlev)

    ! Cloud cover of a pair of layers
    real(jprb), intent(out) :: pair_cloud_cover(nlev-1)

    ! Local variables

    ! Cumulative product needed in computation of total_cloud_cover
    real(jprb) :: cum_product

    ! Loop index for model level
    integer :: jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud_cover:cum_cloud_cover_max_ran',0,hook_handle)

    ! Loop to compute total cloud cover and the cumulative cloud cover

    ! down to the base of each layer
    cum_product = 1.0_jprb - frac(1)
    cum_cloud_cover(1) = frac(1)
    do jlev = 1,nlev-1
      ! Compute the combined cloud cover of layers jlev and jlev+1
      pair_cloud_cover(jlev) = max(frac(jlev),frac(jlev+1))

      if (frac(jlev) >= MaxCloudFrac) then
        ! Cloud cover has reached one
        cum_product = 0.0_jprb
      else
        cum_product = cum_product * (1.0_jprb - pair_cloud_cover(jlev)) &
             &  / (1.0_jprb - frac(jlev))
      end if
      cum_cloud_cover(jlev+1) = 1.0_jprb - cum_product
    end do

    if (lhook) call dr_hook('radiation_cloud_cover:cum_cloud_cover_max_ran',1,hook_handle)

  end subroutine cum_cloud_cover_max_ran


  !---------------------------------------------------------------------
  ! Exponential-random overlap: exponential overlap for contiguous
  ! clouds, random overlap for non-contiguous clouds
  subroutine cum_cloud_cover_exp_ran(nlev, frac, overlap_param, &
       & cum_cloud_cover, pair_cloud_cover, is_beta_overlap)

    use yomhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! Inputs
    integer, intent(in)     :: nlev  ! number of model levels

    ! Cloud fraction on full levels
    real(jprb), intent(in)  :: frac(nlev)

    ! Cloud overlap parameter for interfaces between model layers,
    ! where 0 indicates random overlap and 1 indicates maximum-random
    ! overlap
    real(jprb), intent(in)  :: overlap_param(nlev-1)

    ! This routine has been coded using the "alpha" overlap parameter
    ! of Hogan and Illingworth (2000). If the following logical is
    ! present and true then the input is interpretted to be the "beta"
    ! overlap parameter of Shonk et al. (2010), and needs to be
    ! converted to alpha.
    logical, intent(in), optional :: is_beta_overlap

    ! Outputs

    ! Cumulative cloud cover from TOA to the base of each layer
    real(jprb), intent(out) :: cum_cloud_cover(nlev)

    ! Cloud cover of a pair of layers
    real(jprb), intent(out) :: pair_cloud_cover(nlev-1)

    ! Local variables

    ! Cumulative product needed in computation of total_cloud_cover
    real(jprb) :: cum_product

    ! "Alpha" overlap parameter
    real(jprb) :: overlap_alpha
    logical    :: do_overlap_conversion

    ! Loop index for model level
    integer :: jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud_cover:cum_cloud_cover_exp_ran',0,hook_handle)


    if (present(is_beta_overlap)) then
      do_overlap_conversion = is_beta_overlap
    else
      do_overlap_conversion = .false.
    end if

    ! Loop to compute total cloud cover and the cumulative cloud cover

    ! down to the base of each layer
    cum_product = 1.0_jprb - frac(1)
    cum_cloud_cover(1) = frac(1)
    do jlev = 1,nlev-1
      ! Convert to "alpha" overlap parameter if necessary
      if (do_overlap_conversion) then
        overlap_alpha = beta2alpha(overlap_param(jlev), &
             &                     frac(jlev), frac(jlev+1))
      else
        overlap_alpha = overlap_param(jlev)
      end if

      ! Compute the combined cloud cover of layers jlev and jlev+1
      pair_cloud_cover(jlev) = overlap_alpha*max(frac(jlev),frac(jlev+1)) &
           &  + (1.0_jprb - overlap_alpha) &
           &  * (frac(jlev)+frac(jlev+1)-frac(jlev)*frac(jlev+1))
! Added for DWD (2020)
#ifdef __SX__
    end do
    do jlev = 1,nlev-1
#endif
      if (frac(jlev) >= MaxCloudFrac) then
        ! Cloud cover has reached one
        cum_product = 0.0_jprb
      else
        cum_product = cum_product * (1.0_jprb - pair_cloud_cover(jlev)) &
             &  / (1.0_jprb - frac(jlev))
      end if
      cum_cloud_cover(jlev+1) = 1.0_jprb - cum_product
    end do


    if (lhook) call dr_hook('radiation_cloud_cover:cum_cloud_cover_exp_ran',1,hook_handle)

  end subroutine cum_cloud_cover_exp_ran



  !---------------------------------------------------------------------
  ! Exponential-exponential overlap: exponential overlap for both
  ! contiguous and non-contiguous clouds. This is the result of the
  ! simple Raisanen cloud generator, but unfortunately it has no
  ! (known) analytic formula for the total cloud cover, or the
  ! cumulative cloud cover.  In partially cloudy columns, The McICA
  ! scheme needs this info in order to devote all the cloudy g-points
  ! to columns containing cloud, which reduces McICA noise. The
  ! following routine provides an approximate estimate of cumulative
  ! cloud cover consistent with the exponential-exponential scheme.
  subroutine cum_cloud_cover_exp_exp(nlev, frac, overlap_param, &
       & cum_cloud_cover, pair_cloud_cover, is_beta_overlap)

    use yomhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! Inputs
    integer, intent(in)     :: nlev  ! number of model levels

    ! Cloud fraction on full levels
    real(jprb), intent(in)  :: frac(nlev)

    ! Cloud overlap parameter for interfaces between model layers,
    ! where 0 indicates random overlap and 1 indicates maximum-random
    ! overlap
    real(jprb), intent(in)  :: overlap_param(nlev-1)

    ! This routine has been coded using the "alpha" overlap parameter
    ! of Hogan and Illingworth (2000). If the following logical is
    ! present and true then the input is interpretted to be the "beta"
    ! overlap parameter of Shonk et al. (2010), and needs to be
    ! converted to alpha.
    logical, intent(in), optional :: is_beta_overlap

    ! Outputs

    ! Cumulative cloud cover from TOA to the base of each layer
    real(jprb), intent(out) :: cum_cloud_cover(nlev)

    ! Cloud cover of a pair of layers
    real(jprb), intent(out) :: pair_cloud_cover(nlev-1)

    ! Local variables

    ! If this routine is called from the radiation_interface tree then
    ! very low cloud fractions have already been set to zero, but if
    ! it is called as a cloud cover diagnostic then this can't be
    ! guaranteed so a small non-zero numbers is required
    real(jprb), parameter :: min_frac = 1.0e-6_jprb

    ! "Alpha" overlap parameter
    real(jprb) :: overlap_alpha(nlev-1)
    logical    :: do_overlap_conversion

    ! Variables describing "concave cloud objects", i.e. contiguous
    ! layers where the cloud fraction monotonically increases then
    ! monotonically decreases

    ! Number of objects
    integer    :: nobj

    ! Indices to the location of the top, maximum cloud fraction, and
    ! base, of each concave cloud object
    integer    :: i_top_obj(nlev)
    integer    :: i_max_obj(nlev)
    integer    :: i_base_obj(nlev)
    ! Poor-man's linked list to allow for deletion of objects: this
    ! variable points to the index of the next active object
    integer    :: i_next_obj(nlev)

    ! Cloud cover of object
    real(jprb) :: cc_obj(nlev)

    ! Overlap parameter between objects
    real(jprb) :: alpha_obj(nlev)

    ! Do (while) loop index for model level
    integer :: jlev

    ! Do loop index for object
    integer :: jobj

    ! Maximum correlation between adjacent objects
    real(jprb) :: alpha_max

    ! Indices to pair of objects
    integer    :: iobj1, iobj2

    ! Combined cloud cover of pair of objects, and scaling to modify
    ! cumulative cloud cover of lower layer
    real(jprb) :: cc_pair, scaling

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud_cover:cum_cloud_cover_exp_exp',0,hook_handle)


    if (present(is_beta_overlap)) then
      do_overlap_conversion = is_beta_overlap
    else
      do_overlap_conversion = .false.
    end if

    ! Loop down through atmosphere to locate objects and compute their
    ! basic properties
    jlev = 1
    nobj = 0
    do while (jlev <= nlev)
      if (frac(jlev) > min_frac) then
        ! Starting a new object: store its top
        nobj = nobj + 1
        i_top_obj(nobj) = jlev;
        ! Find its maximum cloud fraction
        jlev = jlev + 1
        do while (jlev <= nlev)
          if (frac(jlev) < frac(jlev-1)) then
            exit
          end if
          jlev = jlev + 1
        end do
        i_max_obj(nobj) = jlev - 1
        ! Find its base
        do while (jlev <= nlev)
          if (frac(jlev) > frac(jlev-1) .or. frac(jlev) <= min_frac) then
            exit
          end if
          jlev = jlev + 1
        end do
        ! In the case of cloud fraction profile starting from the top
        ! like this: 0.1 0.2 0.1 0.2 0.1, we may want the object grouping
        ! to be (0.1 0.2) (0.1 0.2 0.1), not (0.1 0.2 0.1) (0.2 0.1),
        ! in which case the following should be uncommented
        !        if (jlev < nlev) then
        !          if (frac(jlev) > frac(jlev-1)) then
        !            jlev = jlev - 1
        !          end if
        !        end if
        i_base_obj(nobj) = jlev - 1
        ! Index to the next object
        i_next_obj(nobj) = nobj+1
      else
        jlev = jlev + 1
      end if
    end do

    ! Array assignments
    cum_cloud_cover = 0.0_jprb
    pair_cloud_cover = 0.0_jprb

    if (nobj > 0) then
      ! Only do any more work if there is cloud present

      ! To minimize the potential calls to beta2alpha, we do all the
      ! computations related to overlap parameter here
      if (.not. do_overlap_conversion) then

        ! Compute the combined cloud cover of pairs of layers
        do jlev = 1,nlev-1
          pair_cloud_cover(jlev) &
               &  = overlap_param(jlev)*max(frac(jlev),frac(jlev+1)) &
               &  + (1.0_jprb - overlap_param(jlev)) &
               &  * (frac(jlev)+frac(jlev+1)-frac(jlev)*frac(jlev+1))
        end do
        ! Estimate the effective overlap parameter "alpha_obj" between
        ! adjacent objects as the product of the layerwise overlap
        ! parameters between their layers of maximum cloud fraction
        do jobj = 1,nobj-1
          alpha_obj(jobj) &
               &  = product(overlap_param(i_max_obj(jobj):i_max_obj(jobj+1)-1))
        end do

      else

        ! Convert Shonk et al overlap parameter to Hogan and
        ! Illingworth definition
        overlap_alpha = beta2alpha(overlap_param, &
             &                     frac(1:nlev-1), frac(2:nlev))
        ! Compute the combined cloud cover of pairs of layers
        do jlev = 1,nlev-1
          pair_cloud_cover(jlev) &
               &  = overlap_alpha(jlev)*max(frac(jlev),frac(jlev+1)) &
               &  + (1.0_jprb - overlap_alpha(jlev)) &
               &  * (frac(jlev)+frac(jlev+1)-frac(jlev)*frac(jlev+1))
        end do
        ! Estimate the effective overlap parameter "alpha_obj" between
        ! adjacent objects as the product of the layerwise overlap
        ! parameters between their layers of maximum cloud fraction
        do jobj = 1,nobj-1
          alpha_obj(jobj) &
               &  = product(overlap_alpha(i_max_obj(jobj):i_max_obj(jobj+1)-1))
        end do

      end if


      ! Compute the cumulative cloud cover working down from the top
      ! of each object: this will later be converted to the cumulative
      ! cloud cover working down from TOA
      do jobj = 1,nobj
        cum_cloud_cover(i_top_obj(jobj)) = frac(i_top_obj(jobj))
        do jlev = i_top_obj(jobj), i_base_obj(jobj)-1
          if (frac(jlev) >= MaxCloudFrac) then
            ! Cloud cover has reached one
            cum_cloud_cover(jlev+1) = 1.0_jprb
          else
            cum_cloud_cover(jlev+1) = 1.0_jprb &
                 &  - (1.0_jprb - cum_cloud_cover(jlev)) &
                 &  * (1.0_jprb - pair_cloud_cover(jlev))  &
                 &  / (1.0_jprb - frac(jlev))
          end if
        end do
        cc_obj(jobj) = cum_cloud_cover(i_base_obj(jobj))
      end do

      iobj1 = 1

      ! Sequentially combine objects until there is only one left
      ! covering the full vertical extent of clouds in the profile
      do while (nobj > 1)
        ! Find the most correlated adjacent pair of objects
        alpha_max = 0.0_jprb

        ! Need to re-initialize iobj1 here in case alpha_obj(:)==0.0,
        ! which would mean that the "if" statement in the following
        ! loop would never get triggered
        iobj1 = 1

        jobj = 1
        do while (jobj < nobj)
          if (alpha_obj(jobj) > alpha_max) then
            alpha_max = alpha_obj(jobj)
            iobj1 = jobj
          end if
          jobj = i_next_obj(jobj)
        end do

        ! iobj1 is the index to the first object in the pair, set
        ! iobj2 to the second
        iobj2 = i_next_obj(iobj1)

        ! Set the cumulative cloud cover in the clear-sky gap between
        ! the objects to the value at the base of the upper object
        cum_cloud_cover(i_base_obj(iobj1)+1:i_top_obj(iobj2)-1) &
             &  = cum_cloud_cover(i_base_obj(iobj1))

        ! Calculate the combined cloud cover of the pair of objects
        cc_pair = alpha_obj(iobj1)*max(cc_obj(iobj1), cc_obj(iobj2)) &
             &  + (1.0_jprb - alpha_obj(iobj1)) &
             &  * (cc_obj(iobj1) + cc_obj(iobj2) - cc_obj(iobj1)*cc_obj(iobj2))
        scaling = min(max((cc_pair-cc_obj(iobj1)) / max(min_frac, cc_obj(iobj2)), &
             &            0.0_jprb), &
             &        1.0_jprb)

        ! Scale the combined cloud cover of the lower object to
        ! account for its overlap with the upper object
        do jlev = i_top_obj(iobj2),i_base_obj(iobj2)
          cum_cloud_cover(jlev) = cum_cloud_cover(i_base_obj(iobj1)) &
               +  cum_cloud_cover(jlev) * scaling
        end do

        ! Merge the objects by setting the properties of the upper
        ! object to the combined properties of both.  Note that
        ! i_max_obj is not modified because it is no longer needed.
        cc_obj(iobj1) = cc_pair
        i_base_obj(iobj1) = i_base_obj(iobj2)
        i_next_obj(iobj1) = i_next_obj(iobj2)
        alpha_obj(iobj1)  = alpha_obj(iobj2)
        nobj = nobj - 1
      end do

      ! Finish off the total cloud cover below cloud
      cum_cloud_cover(i_base_obj(iobj1)+1:nlev) &
           &  = cum_cloud_cover(i_base_obj(iobj1))

      ! Ensure that the combined cloud cover of pairs of layers is
      ! consistent with the overhang
      do jlev = 1,nlev-1
        pair_cloud_cover(jlev) = max(pair_cloud_cover(jlev), &
             &     frac(jlev)+cum_cloud_cover(jlev+1)-cum_cloud_cover(jlev))
      end do

      ! Sometimes round-off error can lead to cloud cover just above
      ! one, which in turn can lead to direct shortwave fluxes just
      ! below zero
      cum_cloud_cover = min(cum_cloud_cover, 1.0_jprb)

    end if ! cloud is present in profile

    if (lhook) call dr_hook('radiation_cloud_cover:cum_cloud_cover_exp_exp',1,hook_handle)

  end subroutine cum_cloud_cover_exp_exp

end module radiation_cloud_cover
