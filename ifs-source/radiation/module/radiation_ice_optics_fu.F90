! radiation_ice_optics_fu.F90 - Fu's scheme for ice optical properties
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
!   2020-08-10  R. Hogan  Bounded re to be <= 100um and g to be < 1.0

module radiation_ice_optics_fu

  use parkind1, only : jprb

  implicit none
  public

  ! The number of ice coefficients depends on the parameterization
  integer, parameter :: NIceOpticsCoeffsFuSW  = 10
  integer, parameter :: NIceOpticsCoeffsFuLW  = 11

  ! Limits based on the range of validity of the parameterizations
  real(jprb), parameter :: MaxAsymmetryFactor = 1.0_jprb - 10.0_jprb*epsilon(1.0_jprb)
  real(jprb), parameter :: MaxEffectiveRadius = 100.0e-6_jprb ! metres

contains

  !---------------------------------------------------------------------
  ! Compute shortwave ice-particle scattering properties using Fu
  ! (1996) parameterization.  The asymmetry factor in band 14 goes
  ! larger than one for re > 100.8 um, so we cap re at 100 um.
  ! Asymmetry factor is capped at just less than 1 because if it is
  ! exactly 1 then delta-Eddington scaling leads to a zero scattering
  ! optical depth and then division by zero.
  subroutine calc_ice_optics_fu_sw(nb, coeff, ice_wp, &
       &  re, od, scat_od, g)

    !use yomhook,  only : lhook, dr_hook, jphook

    ! Number of bands
    integer, intent(in)  :: nb
    ! Coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! Ice water path (kg m-2)
    real(jprb), intent(in) :: ice_wp
    ! Effective radius (m)
    real(jprb), intent(in) :: re
    ! Total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    ! Fu's effective diameter (microns) and its inverse
    real(jprb) :: de_um, inv_de_um
    ! Ice water path in g m-2
    real (jprb) :: iwp_gm_2

    integer :: jb
    !real(jphook) :: hook_handle

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_fu_sw',0,hook_handle)

    ! Convert to effective diameter using the relationship in the IFS
    de_um     = min(re, MaxEffectiveRadius) * (1.0e6_jprb / 0.64952_jprb)
    inv_de_um = 1.0_jprb / de_um
    iwp_gm_2  = ice_wp * 1000.0_jprb

! Added for DWD (2020)
!NEC$ shortloop
    do jb = 1, nb
      od(jb) = iwp_gm_2 * (coeff(jb,1) + coeff(jb,2) * inv_de_um)
      scat_od(jb) = od(jb) * (1.0_jprb - (coeff(jb,3) + de_um*(coeff(jb,4) &
         &  + de_um*(coeff(jb,5) + de_um*coeff(jb,6)))))
      g(jb) = min(coeff(jb,7) + de_um*(coeff(jb,8) &
         &  + de_um*(coeff(jb,9) + de_um*coeff(jb,10))), &
         &  MaxAsymmetryFactor)
    end do

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_fu_sw',1,hook_handle)

  end subroutine calc_ice_optics_fu_sw


  !---------------------------------------------------------------------
  ! Compute longwave ice-particle scattering properties using Fu et
  ! al. (1998) parameterization
  subroutine calc_ice_optics_fu_lw(nb, coeff, ice_wp, &
       &  re, od, scat_od, g)

    !use yomhook,  only : lhook, dr_hook, jphook

    ! Number of bands
    integer, intent(in)  :: nb
    ! Coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! Ice water path (kg m-2)
    real(jprb), intent(in) :: ice_wp
    ! Effective radius (m)
    real(jprb), intent(in) :: re
    ! Total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    ! Fu's effective diameter (microns) and its inverse
    real(jprb) :: de_um, inv_de_um
    ! Ice water path in g m-2
    real (jprb) :: iwp_gm_2

    integer :: jb
    !real(jphook) :: hook_handle

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_fu_lw',0,hook_handle)

    ! Convert to effective diameter using the relationship in the IFS
    de_um = min(re, MaxEffectiveRadius) * (1.0e6_jprb / 0.64952_jprb)

    inv_de_um = 1.0_jprb / de_um
    iwp_gm_2  = ice_wp * 1000.0_jprb

! Added for DWD (2020)
!NEC$ shortloop
    do jb = 1, nb
      od(jb) = iwp_gm_2 * (coeff(jb,1) + inv_de_um*(coeff(jb,2) &
         &  + inv_de_um*coeff(jb,3)))
      scat_od(jb) = od(jb) - iwp_gm_2*inv_de_um*(coeff(jb,4) + de_um*(coeff(jb,5) &
         &  + de_um*(coeff(jb,6) + de_um*coeff(jb,7))))
      g(jb) = min(coeff(jb,8) + de_um*(coeff(jb,9) &
         &  + de_um*(coeff(jb,10) + de_um*coeff(jb,11))), &
         &  MaxAsymmetryFactor)
    end do

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_fu_lw',1,hook_handle)

  end subroutine calc_ice_optics_fu_lw

end module radiation_ice_optics_fu
