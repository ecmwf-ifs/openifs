! radiation_optical_depth_scaling.h - Cloud optical-depth scaling for Tripleclouds 
!
! Copyright (C) 2016-2017 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! Modifications
!   2017-07-14  R. Hogan  Incorporate gamma distribution option
!
! This file is intended to be included inside a module to ensure that
! this simple routine may be inlined

!---------------------------------------------------------------------
! Compute the optical depth scalings for the optically "thick" and
! "thin" regions of a Tripleclouds representation of a sub-grid PDF of
! cloud optical depth. Following Shonk and Hogan (2008), the 16th
! percentile is used for the thin region, and the formulas estimate
! this for both lognormal and gamma distributions.
pure subroutine optical_depth_scaling(nreg, frac_std, do_gamma, od_scaling)

  use parkind1, only : jprb

  ! Number of regions
  integer, intent(in)     :: nreg

  ! Fractional standard deviation of in-cloud water content
  real(jprb), intent(in)  :: frac_std

  ! Do we do a lognormal or gamma distribution?
  logical, intent(in) :: do_gamma

  ! Optical depth scaling for the cloudy regions
  real(jprb), intent(out) :: od_scaling(2:nreg)

  if (nreg == 2) then
    ! Only one clear-sky and one cloudy region: cloudy region is
    ! homogeneous
    od_scaling(2) = 1.0_jprb
  else
    ! Two cloudy regions with optical depth scaled by 1-x and
    ! 1+x.
    ! Simple version which fails when fractional_std >= 1:
    !od_scaling(2) = 1.0_jprb-cloud%fractional_std(jcol,jlev)
    ! According to Shonk and Hogan (2008), 1-x should correspond to
    ! the 16th percentile. 
    if (.not. do_gamma) then
      ! If we treat the distribution as a lognormal such that the
      ! equivalent Normal has a mean mu and standard deviation sigma,
      ! then the 16th percentile of the lognormal is very close to
      ! exp(mu-sigma).
      od_scaling(2) &
           &  = exp(-sqrt(log(frac_std**2+1))) / sqrt(frac_std**2+1)
    else
      ! If we treat the distribution as a gamma then the 16th
      ! percentile is close to the following
      od_scaling(2) = exp(-frac_std*(1.0_jprb + 0.5_jprb*frac_std &
           &                                   *(1.0_jprb+0.5_jprb*frac_std)))
    end if

    ! Ensure mean optical depth is conserved
    od_scaling(3) = 2.0_jprb-od_scaling(2)
  end if

end subroutine optical_depth_scaling
