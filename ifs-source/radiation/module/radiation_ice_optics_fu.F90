! radiation_ice_optics_fu.F90 - Fu's scheme for ice optical properties
!
! Copyright (C) 2014-2016 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radiation_ice_optics_fu

  implicit none

  ! The number of ice coefficients depends on the parameterization
  integer, parameter :: NIceOpticsCoeffsFuSW  = 10
  integer, parameter :: NIceOpticsCoeffsFuLW  = 11


contains

  !---------------------------------------------------------------------
  ! Compute shortwave ice-particle scattering properties using Fu
  ! (1996) parameterization
  subroutine calc_ice_optics_fu_sw(nb, coeff, ice_wp, &
       &  re, od, scat_od, g)

    use parkind1, only : jprb
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

    !real(jprb)  :: hook_handle

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_fu_sw',0,hook_handle)

    ! Convert to effective diameter using the relationship in the IFS
    de_um     = re * (1.0e6_jprb / 0.64952_jprb)
    inv_de_um = 1.0_jprb / de_um
    iwp_gm_2  = ice_wp * 1000.0_jprb
    !    de_um = 20.0_jprb

    od = iwp_gm_2 * (coeff(1:nb,1) + coeff(1:nb,2) * inv_de_um)
    scat_od = od * (1.0_jprb - (coeff(1:nb,3) + de_um*(coeff(1:nb,4) &
         &  + de_um*(coeff(1:nb,5) + de_um*coeff(1:nb,6)))))
    g = min(coeff(1:nb,7) + de_um*(coeff(1:nb,8) &
         &  + de_um*(coeff(1:nb,9) + de_um*coeff(1:nb,10))), 1.0_jprb)

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_fu_sw',1,hook_handle)

  end subroutine calc_ice_optics_fu_sw


  !---------------------------------------------------------------------
  ! Compute longwave ice-particle scattering properties using Fu et
  ! al. (1998) parameterization
  subroutine calc_ice_optics_fu_lw(nb, coeff, ice_wp, &
       &  re, od, scat_od, g)

    use parkind1, only : jprb
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

    !real(jprb)  :: hook_handle

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_fu_lw',0,hook_handle)

    ! Convert to effective diameter using the relationship in the IFS
    de_um = re * (1.0e6_jprb / 0.64952_jprb)
    !    de_um = 20.0_jprb

    inv_de_um = 1.0_jprb / de_um
    iwp_gm_2  = ice_wp * 1000.0_jprb

    od = iwp_gm_2 * (coeff(1:nb,1) + inv_de_um*(coeff(1:nb,2) &
         &  + inv_de_um*coeff(1:nb,3)))
    scat_od = od - iwp_gm_2*inv_de_um*(coeff(1:nb,4) + de_um*(coeff(1:nb,5) &
         &  + de_um*(coeff(1:nb,6) + de_um*coeff(1:nb,7))))
    g = min(coeff(1:nb,8) + de_um*(coeff(1:nb,9) &
         &  + de_um*(coeff(1:nb,10) + de_um*coeff(1:nb,11))), 1.0_jprb)

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_fu_lw',1,hook_handle)

  end subroutine calc_ice_optics_fu_lw

end module radiation_ice_optics_fu
