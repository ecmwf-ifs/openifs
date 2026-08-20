! radiation_cloud_generator.F90 - Generate water-content or optical-depth scalings for McICA
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
! Generate clouds for McICA using a method modified from Raisanen et
! al. (2002)
!
! Modifications
!   2018-02-22  R. Hogan  Call masked version of PDF sampler for speed
!   2020-03-31  R. Hogan  More vectorizable version of Exp-Ran

module radiation_cloud_generator

  public

contains

  !---------------------------------------------------------------------
  ! Generate scaling factors for the cloud optical depth to represent
  ! cloud overlap, the overlap of internal cloud inhomogeneities and
  ! the fractional standard deviation of these inhomogeneities, for
  ! use in a Monte Carlo Independent Column Approximation radiation
  ! scheme. All returned profiles contain cloud, and the total cloud
  ! cover is also returned, so the calling function can then do a
  ! weighted average of clear and cloudy skies; this is a way to
  ! reduce the Monte Carlo noise in profiles with low cloud cover.
  subroutine cloud_generator(ng, nlev, i_overlap_scheme, &
       &  iseed, frac_threshold, &
       &  frac, overlap_param, decorrelation_scaling, &
       &  fractional_std, pdf_sampler, &
       &  od_scaling, total_cloud_cover, &
       &  use_beta_overlap, use_vectorizable_generator)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook
    use radiation_io,   only     : nulerr, radiation_abort
    use random_numbers_mix, only : randomnumberstream, &
         initialize_random_numbers, uniform_distribution
    use radiation_pdf_sampler, only : pdf_sampler_type
    use radiation_cloud_cover, only : cum_cloud_cover_exp_ran, &
         &       cum_cloud_cover_max_ran, cum_cloud_cover_exp_exp, &
         &       IOverlapMaximumRandom, IOverlapExponentialRandom, &
         &       IOverlapExponential

    implicit none

    ! Inputs
    integer, intent(in)     :: ng    ! number of g points
    integer, intent(in)     :: nlev  ! number of model levels
    integer, intent(in)     :: i_overlap_scheme
    integer, intent(in)     :: iseed ! seed for random number generator

    ! Only cloud fractions above this threshold are considered to be
    ! clouds
    real(jprb), intent(in)  :: frac_threshold

    ! Cloud fraction on full levels
    real(jprb), intent(in)  :: frac(nlev)

    ! Cloud overlap parameter for interfaces between model layers,
    ! where 0 indicates random overlap and 1 indicates maximum-random
    ! overlap
    real(jprb), intent(in)  :: overlap_param(nlev-1)

    ! Overlap parameter for internal inhomogeneities
    real(jprb), intent(in)  :: decorrelation_scaling

    ! Fractional standard deviation at each layer
    real(jprb), intent(in)  :: fractional_std(nlev)

    ! Object for sampling from a lognormal or gamma distribution
    type(pdf_sampler_type), intent(in) :: pdf_sampler

    ! This routine has been coded using the "alpha" overlap parameter
    ! of Hogan and Illingworth (2000). If the following logical is
    ! present and true then the input is interpretted to be the "beta"
    ! overlap parameter of Shonk et al. (2010), and needs to be
    ! converted to alpha.
    logical, intent(in), optional :: use_beta_overlap

    ! Do we use the more vectorizable cloud generator, at the expense
    ! of more random numbers being needed?
    logical, intent(in), optional :: use_vectorizable_generator

    ! Outputs

    ! Cloud optical depth scaling factor, with 0 indicating clear sky
    real(jprb), intent(out) :: od_scaling(ng,nlev)

    ! Total cloud cover using cloud fraction and overlap parameter
    real(jprb), intent(out) :: total_cloud_cover

    ! Local variables

    ! Cumulative cloud cover from TOA to the base of each layer
    real(jprb) :: cum_cloud_cover(nlev)

    ! Scaled random number for finding cloud
    real(jprb) :: trigger

    ! Uniform deviates between 0 and 1
    real(jprb) :: rand_top(ng)

    ! Overlap parameter of inhomogeneities
    real(jprb) :: overlap_param_inhom(nlev-1)

    ! Seed for random number generator and stream for producing random
    ! numbers
    type(randomnumberstream) :: random_stream
    
    ! First and last cloudy layers
    integer :: ibegin, iend

    integer :: itrigger

    ! Loop index for model level and g-point
    integer :: jlev, jg

    ! Cloud cover of a pair of layers, and amount by which cloud at
    ! next level increases total cloud cover as seen from above
    real(jprb), dimension(nlev-1) :: pair_cloud_cover, overhang

    logical :: use_vec_gen

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud_generator:cloud_generator',0,hook_handle)

    if (i_overlap_scheme == IOverlapExponentialRandom) then
      call cum_cloud_cover_exp_ran(nlev, frac, overlap_param, &
           &   cum_cloud_cover, pair_cloud_cover, use_beta_overlap)
    else if (i_overlap_scheme == IOverlapMaximumRandom) then
      call cum_cloud_cover_max_ran(nlev, frac, &
           &   cum_cloud_cover, pair_cloud_cover)
    else if (i_overlap_scheme == IOverlapExponential) then
      call cum_cloud_cover_exp_exp(nlev, frac, overlap_param, &
           &   cum_cloud_cover, pair_cloud_cover, use_beta_overlap)
    else
      write(nulerr,'(a)') '*** Error: cloud overlap scheme not recognised'
      call radiation_abort()
    end if

    total_cloud_cover = cum_cloud_cover(nlev);
    do jlev = 1,nlev-1
      overhang(jlev) = cum_cloud_cover(jlev+1)-cum_cloud_cover(jlev)
    end do

    if (total_cloud_cover < frac_threshold) then
      ! Treat column as clear sky: calling function therefore will not
      ! use od_scaling so we don't need to calculate it
      total_cloud_cover = 0.0_jprb

    else
      ! Cloud is present: need to calculate od_scaling

      ! Find range of cloudy layers
      jlev = 1
      do while (frac(jlev) <= 0.0_jprb) 
        jlev = jlev + 1
      end do
      ibegin = jlev
      iend = jlev
      do jlev = jlev+1,nlev
        if (frac(jlev) > 0.0_jprb) then
          iend = jlev
        end if
      end do

      ! Set overlap parameter of inhomogeneities
      overlap_param_inhom = overlap_param

      do jlev = ibegin,iend-1
        if (overlap_param(jlev) > 0.0_jprb) then
          overlap_param_inhom(jlev) &
               &  = overlap_param(jlev)**(1.0_jprb/decorrelation_scaling)
        end if
      end do

      ! Reset optical depth scaling to clear skies
      od_scaling = 0.0_jprb

      if (present(use_vectorizable_generator)) then
        use_vec_gen = use_vectorizable_generator
      else
        use_vec_gen = .false.
      end if

      if (.not. use_vec_gen) then
        ! Original generator that minimizes the number of random
        ! numbers used, but is not vectorizable

        ! Expensive operation: initialize random number generator for
        ! this column
        call initialize_random_numbers(iseed, random_stream)

        ! Compute ng random numbers to use to locate cloud top
        call uniform_distribution(rand_top, random_stream)
        
        ! Loop over ng columns
        do jg = 1,ng
          ! Find the cloud top height corresponding to the current
          ! random number, and store in itrigger
          trigger = rand_top(jg) * total_cloud_cover
          jlev = ibegin
          do while (trigger > cum_cloud_cover(jlev) .and. jlev < iend)
            jlev = jlev + 1
          end do
          itrigger = jlev
          
          if (i_overlap_scheme /= IOverlapExponential) then
            call generate_column_exp_ran(ng, nlev, jg, random_stream, pdf_sampler, &
                 &  frac, pair_cloud_cover, &
                 &  cum_cloud_cover, overhang, fractional_std, overlap_param_inhom, &
                 &  itrigger, iend, od_scaling)
          else
            call generate_column_exp_exp(ng, nlev, jg, random_stream, pdf_sampler, &
                 &  frac, pair_cloud_cover, &
                 &  cum_cloud_cover, overhang, fractional_std, overlap_param_inhom, &
                 &  itrigger, iend, od_scaling)
          end if
          
        end do

      else
        ! Alternative generator (only for Exp-Ran overlap so far) that
        ! should be vectorizable but generates more random numbers,
        ! some of which are not used

        if (i_overlap_scheme == IOverlapExponential) then
          write(nulerr,'(a)') '*** Error: vectorizable cloud generator is not available with Exp-Exp overlap'
          call radiation_abort()
        end if

        call generate_columns_exp_ran(ng, nlev, iseed, pdf_sampler, &
             &  total_cloud_cover, frac_threshold, frac, pair_cloud_cover, &
             &  cum_cloud_cover, overhang, fractional_std, overlap_param_inhom, &
             &  ibegin, iend, od_scaling)

      end if

    end if

    if (lhook) call dr_hook('radiation_cloud_generator:cloud_generator',1,hook_handle)

  end subroutine cloud_generator


  !---------------------------------------------------------------------
  ! Generate a column of optical depth scalings using
  ! exponential-random overlap (which includes maximum-random overlap
  ! as a limiting case)
  subroutine generate_column_exp_ran(ng, nlev, ig, random_stream, pdf_sampler, &
       &  frac, pair_cloud_cover, &
       &  cum_cloud_cover, overhang, fractional_std, overlap_param_inhom, &
       &  itrigger, iend, od_scaling)

    use parkind1,              only : jprb
    use radiation_pdf_sampler, only : pdf_sampler_type
    use random_numbers_mix,    only : randomnumberstream, &
         initialize_random_numbers, uniform_distribution


    implicit none

    ! Number of g points / columns, and number of current column
    integer, intent(in) :: ng, ig

    ! Number of levels
    integer, intent(in) :: nlev

    ! Stream for producing random numbers
    type(randomnumberstream), intent(inout) :: random_stream

    ! Object for sampling from a lognormal or gamma distribution
    type(pdf_sampler_type), intent(in) :: pdf_sampler

    ! Cloud fraction, cumulative cloud cover and fractional standard
    ! deviation in each layer
    real(jprb), intent(in), dimension(nlev) :: frac, cum_cloud_cover, fractional_std

    ! Cloud cover of a pair of layers, and amount by which cloud at
    ! next level increases total cloud cover as seen from above
    real(jprb), intent(in), dimension(nlev-1) :: pair_cloud_cover, overhang

    ! Overlap parameter of inhomogeneities
    real(jprb), intent(in), dimension(nlev-1) :: overlap_param_inhom

    ! Top of highest cloudy layer (in this subcolumn) and base of
    ! lowest
    integer, intent(in) :: itrigger, iend

    ! Optical depth scaling to output
    real(jprb), intent(inout), dimension(ng,nlev) :: od_scaling

    ! Height indices
    integer :: jlev, jcloud

    ! Number of contiguous cloudy layers for which to compute optical
    ! depth scaling
    integer :: n_layers_to_scale

    integer :: iy

    ! Is it time to fill the od_scaling variable?
    logical :: do_fill_od_scaling

    real(jprb) :: rand_cloud(nlev)
    real(jprb) :: rand_inhom1(nlev), rand_inhom2(nlev)

    ! So far our vertically contiguous cloud contains only one layer
    n_layers_to_scale = 1
    iy = 0

    ! Locate the clouds below this layer: first generate some more
    ! random numbers
    call uniform_distribution(rand_cloud(1:(iend+1-itrigger)),random_stream)

    ! Loop from the layer below the local cloud top down to the
    ! bottom-most cloudy layer
    do jlev = itrigger+1,iend+1
      do_fill_od_scaling = .false.
      if (jlev <= iend) then
        iy = iy+1
        if (n_layers_to_scale > 0) then
          ! There is a cloud above, in which case the probability
          ! of cloud in the layer below is as follows
          if (rand_cloud(iy)*frac(jlev-1) &
               &  < frac(jlev) + frac(jlev-1) - pair_cloud_cover(jlev-1)) then
            ! Add another cloudy layer
            n_layers_to_scale = n_layers_to_scale + 1
          else 
            ! Reached the end of a contiguous set of cloudy layers and
            ! will compute the optical depth scaling immediately.
            do_fill_od_scaling = .true.
          end if
        else
          ! There is clear-sky above, in which case the
          ! probability of cloud in the layer below is as follows
          if (rand_cloud(iy)*(cum_cloud_cover(jlev-1) - frac(jlev-1)) &
               &  < pair_cloud_cover(jlev-1) - overhang(jlev-1) - frac(jlev-1)) then
            ! A new cloud top
            n_layers_to_scale = 1
          end if
        end if
      else
        ! We are at the bottom of the cloudy layers in the model,
        ! so in a moment need to populate the od_scaling array
        do_fill_od_scaling = .true.
      end if

      if (do_fill_od_scaling) then
        ! We have a contiguous range of layers for which we
        ! compute the od_scaling elements using some random
        ! numbers
        call uniform_distribution(rand_inhom1(1:n_layers_to_scale),random_stream)
        call uniform_distribution(rand_inhom2(1:n_layers_to_scale),random_stream)

        ! Loop through the sequence of cloudy layers
        do jcloud = 2,n_layers_to_scale
          ! Use second random number, and inhomogeneity overlap
          ! parameter, to decide whether the first random number
          ! should be repeated (corresponding to maximum overlap)
          ! or not (corresponding to random overlap)
          if (rand_inhom2(jcloud) &
               &  < overlap_param_inhom(jlev-n_layers_to_scale+jcloud-2)) then
            rand_inhom1(jcloud) = rand_inhom1(jcloud-1)
          end if
        end do
        
        ! Sample from a lognormal or gamma distribution to obtain
        ! the optical depth scalings
        call pdf_sampler%sample(fractional_std(jlev-n_layers_to_scale:jlev-1), &
             & rand_inhom1(1:n_layers_to_scale), od_scaling(ig,jlev-n_layers_to_scale:jlev-1))

        n_layers_to_scale = 0
      end if
          
    end do

  end subroutine generate_column_exp_ran


  !---------------------------------------------------------------------
  ! Generate a column of optical depth scalings using
  ! exponential-exponential overlap
  subroutine generate_column_exp_exp(ng, nlev, ig, random_stream, pdf_sampler, &
       &  frac, pair_cloud_cover, &
       &  cum_cloud_cover, overhang, fractional_std, overlap_param_inhom, &
       &  itrigger, iend, od_scaling)

    use parkind1,              only : jprb
    use radiation_pdf_sampler, only : pdf_sampler_type
    use random_numbers_mix,    only : randomnumberstream, &
         initialize_random_numbers, uniform_distribution

    implicit none

    ! Number of g points / columns, and number of current column
    integer, intent(in) :: ng, ig

    ! Number of levels
    integer, intent(in) :: nlev

    ! Stream for producing random numbers
    type(randomnumberstream), intent(inout) :: random_stream

    ! Object for sampling from a lognormal or gamma distribution
    type(pdf_sampler_type), intent(in) :: pdf_sampler

    ! Cloud fraction, cumulative cloud cover and fractional standard
    ! deviation in each layer
    real(jprb), intent(in), dimension(nlev) :: frac, cum_cloud_cover, fractional_std

    ! Cloud cover of a pair of layers, and amount by which cloud at
    ! next level increases total cloud cover as seen from above
    real(jprb), intent(in), dimension(nlev-1) :: pair_cloud_cover, overhang

    ! Overlap parameter of inhomogeneities
    real(jprb), intent(in), dimension(nlev-1) :: overlap_param_inhom

    ! Top of highest cloudy layer (in this subcolumn) and base of
    ! lowest
    integer, intent(in) :: itrigger, iend

    ! Optical depth scaling to output
    real(jprb), intent(inout), dimension(ng,nlev) :: od_scaling

    ! Height indices
    integer :: jlev, jcloud

    integer :: iy

    real(jprb) :: rand_cloud(nlev)
    real(jprb) :: rand_inhom1(nlev), rand_inhom2(nlev)

    ! For each column analysed, this vector locates the clouds. It is
    ! only actually used for Exp-Exp overlap
    logical :: is_cloudy(nlev)

    ! Number of contiguous cloudy layers for which to compute optical
    ! depth scaling
    integer :: n_layers_to_scale

    iy = 0

    is_cloudy = .false.
    is_cloudy(itrigger) = .true.

    ! Locate the clouds below this layer: first generate some more
    ! random numbers
    call uniform_distribution(rand_cloud(1:(iend+1-itrigger)),random_stream)

    ! Loop from the layer below the local cloud top down to the
    ! bottom-most cloudy layer
    do jlev = itrigger+1,iend
      iy = iy+1
      if (is_cloudy(jlev-1)) then
        ! There is a cloud above, in which case the probability
        ! of cloud in the layer below is as follows
        if (rand_cloud(iy)*frac(jlev-1) &
             &  < frac(jlev) + frac(jlev-1) - pair_cloud_cover(jlev-1)) then
          ! Add another cloudy layer
          is_cloudy(jlev) = .true.
        end if
      else
        ! There is clear-sky above, in which case the
        ! probability of cloud in the layer below is as follows
        if (rand_cloud(iy)*(cum_cloud_cover(jlev-1) - frac(jlev-1)) &
             &  < pair_cloud_cover(jlev-1) - overhang(jlev-1) - frac(jlev-1)) then
            ! A new cloud top
          is_cloudy(jlev) = .true.
        end if
      end if
    end do

    ! We have a contiguous range of layers for which we compute the
    ! od_scaling elements using some random numbers

    ! In the Exp-Exp overlap scheme we do all layers at once
    n_layers_to_scale = iend+1 - itrigger
        
    call uniform_distribution(rand_inhom1(1:n_layers_to_scale),random_stream)
    call uniform_distribution(rand_inhom2(1:n_layers_to_scale),random_stream)
        
    ! Loop through the sequence of cloudy layers
    do jcloud = 2,n_layers_to_scale
      ! Use second random number, and inhomogeneity overlap
      ! parameter, to decide whether the first random number
      ! should be repeated (corresponding to maximum overlap)
      ! or not (corresponding to random overlap)
      if (rand_inhom2(jcloud) &
           &  < overlap_param_inhom(iend-n_layers_to_scale+jcloud-1)) then
        rand_inhom1(jcloud) = rand_inhom1(jcloud-1)
      end if
    end do
        
    ! Sample from a lognormal or gamma distribution to obtain the
    ! optical depth scalings

    ! Masked version assuming values outside the range itrigger:iend
    ! are already zero:
    call pdf_sampler%masked_sample(n_layers_to_scale, &
         &  fractional_std(itrigger:iend), &
         &  rand_inhom1(1:n_layers_to_scale), od_scaling(ig,itrigger:iend), &
         &  is_cloudy(itrigger:iend))
        
    ! ! IFS version:
    ! !$omp simd 
    ! do jlev=itrigger,iend
    !    if (.not. is_cloudy(jlev)) then
    !       od_scaling(ig,jlev) = 0.0_jprb
    !    else
    !       call sample_from_pdf_simd(&
    !            pdf_sampler,fractional_std(jlev),&
    !            rand_inhom1(jlev-itrigger+1), &
    !            od_scaling(ig,jlev))
    !    end if
    ! end do

  end subroutine generate_column_exp_exp


  !---------------------------------------------------------------------
  ! Extract the value of a lognormal distribution with fractional
  ! standard deviation "fsd" corresponding to the cumulative
  ! distribution function value "cdf", and return it in x. Since this
  ! is an elemental subroutine, fsd, cdf and x may be arrays. SIMD version.
  subroutine sample_from_pdf_simd(this, fsd, cdf, x)
    use parkind1,              only : jprb
    use radiation_pdf_sampler, only : pdf_sampler_type
    implicit none
#if defined(__GFORTRAN__) || defined(__PGI) || defined(__NEC__)
#else
    !$omp declare simd(sample_from_pdf_simd) uniform(this) &
    !$omp linear(ref(fsd)) linear(ref(cdf))
#endif
    type(pdf_sampler_type), intent(in)  :: this

    ! Fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb),              intent(in)  :: fsd, cdf

    ! Sample from distribution
    real(jprb),              intent(out) :: x

    ! Index to look-up table
    integer    :: ifsd, icdf

    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    ! Bilinear interpolation with bounds
    wcdf = cdf * (this%ncdf-1) + 1.0_jprb
    icdf = max(1, min(int(wcdf), this%ncdf-1))
    wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

    wfsd = (fsd-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
    ifsd = max(1, min(int(wfsd), this%nfsd-1))
    wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

    x =      (1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
         & + (1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
         & +           wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
         & +           wcdf *          wfsd  * this%val(icdf+1,ifsd+1)

  end subroutine sample_from_pdf_simd


  !---------------------------------------------------------------------
  ! Generate columns of optical depth scalings using
  ! exponential-random overlap (which includes maximum-random overlap
  ! as a limiting case).  This version is intended to work better on
  ! hardware with long vector lengths.  As with all calculations in
  ! this file, we zoom into the fraction of the column with cloud at
  ! any height, so that all spectral intervals see a cloud somewhere.
  ! In the McICA solver, this is combined appropriately with the
  ! clear-sky calculation.
  subroutine generate_columns_exp_ran(ng, nlev, iseed, pdf_sampler, &
       &  total_cloud_cover, frac_threshold, frac, pair_cloud_cover, &
       &  cum_cloud_cover, overhang, fractional_std, overlap_param_inhom, &
       &  ibegin, iend, od_scaling)

    use parkind1,              only : jprb
    use radiation_pdf_sampler, only : pdf_sampler_type
    use radiation_random_numbers, only : rng_type, IRngMinstdVector, IRngNative

    implicit none

    ! Number of g points / columns
    integer, intent(in) :: ng

    ! Number of levels
    integer, intent(in) :: nlev

    integer, intent(in) :: iseed ! seed for random number generator

    ! Stream for producing random numbers
    !type(randomnumberstream) :: random_stream
    type(rng_type) :: random_number_generator

    ! Object for sampling from a lognormal or gamma distribution
    type(pdf_sampler_type), intent(in) :: pdf_sampler

    ! Total cloud cover using cloud fraction and overlap parameter
    real(jprb), intent(in) :: total_cloud_cover

    real(jprb), intent(in) :: frac_threshold

    ! Cloud fraction, cumulative cloud cover and fractional standard
    ! deviation in each layer
    real(jprb), intent(in), dimension(nlev) :: frac, cum_cloud_cover, fractional_std

    ! Cloud cover of a pair of layers, and amount by which cloud at
    ! next level increases total cloud cover as seen from above
    real(jprb), intent(in), dimension(nlev-1) :: pair_cloud_cover, overhang

    ! Overlap parameter of inhomogeneities
    real(jprb), intent(in), dimension(nlev-1) :: overlap_param_inhom

    ! Top of highest cloudy layer and base of lowest
    integer, intent(inout) :: ibegin, iend

    ! Optical depth scaling to output
    real(jprb), intent(inout), dimension(ng,nlev) :: od_scaling

    ! Loop indices
    integer :: jlev, jg

    real(jprb) :: rand_cloud(ng,ibegin:iend)
    real(jprb) :: rand_inhom(ng,ibegin-1:iend), rand_inhom2(ng,ibegin:iend)

    ! Is the cloud fraction above the minimum threshold at each level
    logical :: is_any_cloud(ibegin:iend)

    ! Scaled random number for finding cloud
    real(jprb) :: trigger(ng)

    logical :: is_cloud(ng)    ! Is there cloud at this level and spectral interval?
    logical :: prev_cloud(ng)  ! Was there cloud at level above?
    logical :: first_cloud(ng) ! At level of first cloud counting down from top?
    logical :: found_cloud(ng) ! Cloud found in this column counting down from top?

    is_any_cloud = (frac(ibegin:iend) >= frac_threshold)

    ! Initialize random number generator for this column, and state
    ! that random numbers will be requested in blocks of length the
    ! number of spectral intervals ng.
    call random_number_generator%initialize(IRngMinstdVector, iseed=iseed, &
         &                                  nmaxstreams=ng)

    ! Random numbers to use to locate cloud top
    call random_number_generator%uniform_distribution(trigger)

    ! Random numbers to work out whether to transition vertically from
    ! clear to cloudy, cloudy to clear, clear to clear or cloudy to
    ! cloudy
    call random_number_generator%uniform_distribution(rand_cloud, is_any_cloud)

    ! Random numbers to generate sub-grid cloud structure
    call random_number_generator%uniform_distribution(rand_inhom)
    call random_number_generator%uniform_distribution(rand_inhom2, is_any_cloud)

    trigger = trigger * total_cloud_cover

    ! Initialize logicals for clear-sky above first cloudy layer
    found_cloud = .false.
    is_cloud    = .false.
    first_cloud = .false.

    ! Loop down through layers starting at the first cloudy layer
    do jlev = ibegin,iend

      if (is_any_cloud(jlev)) then

! Added for DWD (2020)
!NEC$ shortloop
        do jg = 1,ng
          ! The intention is that all these operations are vectorizable,
          ! since all are vector operations on vectors of length ng...

          ! Copy the cloud mask between levels
          prev_cloud(jg) = is_cloud(jg)

          ! For each spectral interval, has the first cloud appeared at this level?
          first_cloud(jg) = (trigger(jg) <= cum_cloud_cover(jlev) .and. .not. found_cloud(jg))

          ! ...if so, add to found_cloud
          found_cloud(jg) = found_cloud(jg) .or. first_cloud(jg)

          ! There is cloud at this level either if a first cloud has
          ! appeared, or using separate probability calculations
          ! depending on whether there is a cloud above (given by
          ! prev_cloud)
          is_cloud(jg) = first_cloud(jg) &
               &  .or. found_cloud(jg) .and. merge(rand_cloud(jg,jlev)*frac(jlev-1) &
               &               < frac(jlev)+frac(jlev-1)-pair_cloud_cover(jlev-1), &
               &             rand_cloud(jg,jlev)*(cum_cloud_cover(jlev-1) - frac(jlev-1)) &
               &               < pair_cloud_cover(jlev-1) - overhang(jlev-1) - frac(jlev-1), &
               &             prev_cloud(jg))
          ! The random number determining cloud structure decorrelates
          ! with the one above it according to the overlap parameter,
          ! but always decorrelates if there is clear-sky above.  If
          ! there is clear-sky in the present level, the random number
          ! is set to zero to ensure that the optical depth scaling is
          ! also zero.
          rand_inhom(jg,jlev) = merge(merge(rand_inhom(jg,jlev-1), rand_inhom(jg,jlev), &
               &                           rand_inhom2(jg,jlev) < overlap_param_inhom(jlev-1) &
               &                           .and. prev_cloud(jg)), &
               &                     0.0_jprb, is_cloud(jg))
        end do
      else
        ! No cloud at this level
        is_cloud = .false.
      end if
    end do
       
    ! Sample from a lognormal or gamma distribution to obtain the
    ! optical depth scalings, calling the faster masked version and
    ! assuming values outside the range ibegin:iend are already zero
    call pdf_sampler%masked_block_sample(iend-ibegin+1, ng, &
         &  fractional_std(ibegin:iend), &
         &  rand_inhom(:,ibegin:iend), od_scaling(:,ibegin:iend), &
         &  is_any_cloud)

  end subroutine generate_columns_exp_ran

end module radiation_cloud_generator
