! (C) Copyright 2017- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! radiation_ice_optics_baran2017.f90 - 2017 parameterization of Baran's ice optical properties
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!

module radiation_ice_optics_baran2017

  implicit none
  public

  ! The number of ice coefficients depends on the parameterization
  integer, parameter :: NIceOpticsCoeffsBaran2017 = 9
  integer, parameter :: NIceOpticsGeneralCoeffsBaran2017 = 5

contains


  !---------------------------------------------------------------------
  ! Compute ice-particle scattering properties using a
  ! parameterization as a function of ice water mixing ratio and
  ! temperature
  subroutine calc_ice_optics_baran2017(nb, coeff_gen, coeff, ice_wp, &
       &  qi, temperature, od, scat_od, g)

    use parkind1, only : jprb
    !use yomhook,  only : lhook, dr_hook, jphook

    ! Number of bands
    integer, intent(in)  :: nb
    ! General coefficients read from a data file
    real(jprb), intent(in) :: coeff_gen(:)
    ! Band-specific coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! Ice water path (kg m-2) and mixing ratio (kg kg-1)
    real(jprb), intent(in) :: ice_wp, qi
    ! Temperature (K)
    real(jprb), intent(in) :: temperature
    ! Total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    ! Modified ice mixing ratio, and the same raised to an appropriate power
    real(jprb) :: qi_mod, qi_mod_od, qi_mod_ssa, qi_mod_g

    !real(jphook) :: hook_handle

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_baran2017',0,hook_handle)

    qi_mod     = qi * exp(coeff_gen(1)*(temperature-coeff_gen(2)))
    qi_mod_od  = qi_mod ** coeff_gen(3)
    qi_mod_ssa = qi_mod ** coeff_gen(4)
    qi_mod_g   = qi_mod ** coeff_gen(5)

    od      = ice_wp * (coeff(1:nb,1) + coeff(1:nb,2)/(1.0_jprb+qi_mod_od *coeff(1:nb,3)))
    scat_od = od     * (coeff(1:nb,4) + coeff(1:nb,5)/(1.0_jprb+qi_mod_ssa*coeff(1:nb,6)))
    g       =           coeff(1:nb,7) + coeff(1:nb,8)/(1.0_jprb+qi_mod_g  *coeff(1:nb,9))

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_baran2017',1,hook_handle)

  end subroutine calc_ice_optics_baran2017

end module radiation_ice_optics_baran2017
