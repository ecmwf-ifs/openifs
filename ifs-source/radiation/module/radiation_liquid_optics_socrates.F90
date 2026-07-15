! radiation_liquid_optics_socrates.F90 - SOCRATES method for parameterizing liquid droplet optics
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
!   2020-08-10  R. Hogan  Bounded re to be >=1.2um and <=50um

module radiation_liquid_optics_socrates

  use parkind1, only : jprb

  implicit none
  public

  ! SOCRATES (Edwards-Slingo) parameterizes info on the dependence of
  ! the scattering properties in each band on effective radius in
  ! terms of 16 coefficients
  integer, parameter :: NLiqOpticsCoeffsSOCRATES = 16

  ! Range of valid input effective radius, in microns
  real(jprb), parameter :: MinEffectiveRadius = 1.2e-6
  real(jprb), parameter :: MaxEffectiveRadius = 50.0e-6

contains

  !---------------------------------------------------------------------
  ! Compute liquid-droplet scattering properties using a
  ! parameterization consisting of Pade approximants from the
  ! SOCRATES (Edwards-Slingo) code
  subroutine calc_liq_optics_socrates(nb, coeff, lwp, re_in, od, scat_od, g)

    use parkind1, only : jprb
    !use yomhook,  only : lhook, dr_hook, jphook

    ! Number of bands
    integer, intent(in)  :: nb
    ! Coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! Liquid water path (kg m-2) and effective radius (m)
    real(jprb), intent(in) :: lwp, re_in
    ! Total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    integer    :: jb
    ! Local effective radius (m), after applying bounds
    real(jprb) :: re

    !real(jphook) :: hook_handle

    !if (lhook) call dr_hook('radiation_liquid_optics_socrates:calc_liq_optics_socrates',0,hook_handle)

    ! Apply the bounds of validity to effective radius
    re = max(MinEffectiveRadius, min(re_in, MaxEffectiveRadius))

! Added for DWD (2020)
!NEC$ shortloop
    do jb = 1, nb
      od(jb) = lwp * (coeff(jb,1) + re*(coeff(jb,2) + re*coeff(jb,3))) &
         &  / (1.0_jprb + re*(coeff(jb,4) + re*(coeff(jb,5) &
         &  + re*coeff(jb,6))))
      scat_od(jb) = od(jb) * (1.0_jprb &
         &  - (coeff(jb,7) + re*(coeff(jb,8) + re*coeff(jb,9))) &
         &  / (1.0_jprb + re*(coeff(jb,10) + re*coeff(jb,11))))
      g(jb) = (coeff(jb,12) + re*(coeff(jb,13) + re*coeff(jb,14))) &
         &  / (1.0_jprb + re*(coeff(jb,15) + re*coeff(jb,16)))
    end do

    !if (lhook) call dr_hook('radiation_liquid_optics_socrates:calc_liq_optics_socrates',1,hook_handle)

  end subroutine calc_liq_optics_socrates

end module radiation_liquid_optics_socrates
