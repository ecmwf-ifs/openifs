! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! radiation_liquid_optics_slingo.F90 - Slingo SW & Lindner-Li LW parameterization of liquid droplet optics
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!

module radiation_liquid_optics_slingo

  implicit none
  public

  integer, parameter :: NLiqOpticsCoeffsSlingoSW = 6
  integer, parameter :: NLiqOpticsCoeffsLindnerLiLW = 13

contains

  !---------------------------------------------------------------------
  ! Compute liquid-droplet scattering properties in the shortwave from
  ! Slingo (1989). WARNING: this parameterization is known not to be
  ! very accurate: see Nielsen et al. (GMD 2014).
  subroutine calc_liq_optics_slingo(nb, coeff, lwp, re, od, scat_od, g)

    use parkind1, only : jprb
    !use yomhook,  only : lhook, dr_hook, jphook

    ! Number of bands
    integer, intent(in)  :: nb
    ! Coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! Liquid water path (kg m-2) and effective radius (m)
    real(jprb), intent(in) :: lwp, re
    ! Total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    ! Liquid water path in g m-2
    real(jprb) :: lwp_gm_2
    ! Effective radius in microns, and its inverse
    real(jprb) :: re_um, inv_re_um

    !real(jphook) :: hook_handle

    !if (lhook) call dr_hook('radiation_liquid_optics_slingo:calc_liq_optics_slingo',0,hook_handle)

    lwp_gm_2 = lwp * 1000.0_jprb
    ! Range of validity reported by Slingo (1989): 4.2-16.6 microns
    re_um = min(max(4.2_jprb, re * 1.0e6_jprb), 16.6_jprb)
    inv_re_um = 1.0_jprb / re_um

    od = lwp_gm_2 * (coeff(1:nb,1) + inv_re_um*coeff(1:nb,2))
    scat_od = od * (1.0_jprb - coeff(1:nb,3) - re_um*coeff(1:nb,4))
    g = coeff(1:nb,5) + re_um*coeff(1:nb,6)

    !if (lhook) call dr_hook('radiation_liquid_optics_slingo:calc_liq_optics_slingo',1,hook_handle)

  end subroutine calc_liq_optics_slingo


  !---------------------------------------------------------------------
  ! Compute liquid-droplet scattering properties in the longwave from
  ! Lindner & Li (2000)
  subroutine calc_liq_optics_lindner_li(nb, coeff, lwp, re, od, scat_od, g)

    use parkind1, only : jprb
    use yomhook,  only : lhook, dr_hook, jphook

    ! Number of bands
    integer, intent(in)  :: nb
    ! Coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! Liquid water path (kg m-2) and effective radius (m)
    real(jprb), intent(in) :: lwp, re
    ! Total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    ! Liquid water path in g m-2
    real(jprb) :: lwp_gm_2
    ! Effective radius in microns, and its inverse
    real(jprb) :: re_um, inv_re_um

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_liquid_optics_slingo:calc_liq_optics_lindner_li',0,hook_handle)

    lwp_gm_2 = lwp * 1000.0_jprb
    ! Range of validity reported by Lindner and Li (2000): 2-40 microns
    re_um = min(max(2.0_jprb, re * 1.0e6_jprb), 40.0_jprb)
    inv_re_um = 1.0_jprb / re_um

    od = lwp_gm_2 * (coeff(1:nb,1) + re_um*coeff(1:nb,2) + inv_re_um*(coeff(1:nb,3) &
         &  + inv_re_um*(coeff(1:nb,4) + inv_re_um*coeff(1:nb,5))))
    scat_od = od * (1.0_jprb - (coeff(1:nb,6) + inv_re_um*coeff(1:nb,7) &
         &                      + re_um*(coeff(1:nb,8) + re_um*coeff(1:nb,9))))
    g = coeff(1:nb,10) + inv_re_um*coeff(1:nb,11) &
         &  + re_um*(coeff(1:nb,12) + re_um*coeff(1:nb,13))

    if (lhook) call dr_hook('radiation_liquid_optics_slingo:calc_liq_optics_lindner_li',1,hook_handle)

  end subroutine calc_liq_optics_lindner_li

end module radiation_liquid_optics_slingo
