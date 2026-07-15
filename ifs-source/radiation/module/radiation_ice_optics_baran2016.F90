! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! radiation_ice_optics_baran2016.f90 - Baran et al. (2016) scheme for ice optical properties
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!

module radiation_ice_optics_baran2016

  implicit none

  ! The number of ice coefficients depends on the parameterization
  integer, parameter :: NIceOpticsCoeffsBaran2016 = 5

contains


  !---------------------------------------------------------------------
  ! Compute ice-particle scattering properties using a
  ! parameterization as a function of ice water mixing ratio and
  ! temperature
  subroutine calc_ice_optics_baran2016(nb, coeff, ice_wp, &
       &  qi, temperature, od, scat_od, g)

    use parkind1, only : jprb
    !use yomhook,  only : lhook, dr_hook, jphook

    ! Number of bands
    integer, intent(in)  :: nb
    ! Coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! Ice water path (kg m-2) and mixing ratio (kg kg-1)
    real(jprb), intent(in) :: ice_wp, qi
    ! Temperature (K)
    real(jprb), intent(in) :: temperature
    ! Total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    ! Powers of temperature, some multiplied by qi
    real(jprb) :: qi_T, T2, qi_over_T4

    !real(jphook) :: hook_handle

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_baran2016',0,hook_handle)

    T2 = temperature * temperature

    if (qi < 1.0e-3_jprb) then
      qi_T = qi * temperature
      qi_over_T4 = 1.0_jprb / (T2 * T2)
    else
      qi_T = 1.0e-3_jprb * temperature
      qi_over_T4 = 1.0_jprb / (T2 * T2)
    end if

    od      = ice_wp * coeff(1:nb,1) * qi_over_T4
    scat_od = od * (coeff(1:nb,2) + coeff(1:nb,3) * qi_T)
    g       = coeff(1:nb,4) + coeff(1:nb,5) * qi_T

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_baran2016',1,hook_handle)

  end subroutine calc_ice_optics_baran2016

end module radiation_ice_optics_baran2016
