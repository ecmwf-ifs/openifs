! radiation_cloud_optics.F90 - Computing cloud optical properties
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
!   2017-07-22  R. Hogan  Added Yi et al. ice optics model

module radiation_cloud_optics

  implicit none

  public

contains

  ! Provides elemental function "delta_eddington_scat_od"
#include "radiation_delta_eddington.h"

  !---------------------------------------------------------------------
  ! Load cloud scattering data; this subroutine delegates to one
  ! in radiation_cloud_optics_data.F90, but checks the size of
  ! what is returned
  subroutine setup_cloud_optics(config)

    use parkind1,         only : jprb
    use yomhook,          only : lhook, dr_hook, jphook

    use radiation_io,     only : nulerr, radiation_abort
    use radiation_config, only : config_type, IIceModelFu, IIceModelBaran, &
         &                       IIceModelBaran2016, IIceModelBaran2017, &
         &                       IIceModelYi, &
         &                       ILiquidModelSOCRATES, ILiquidModelSlingo
    use radiation_cloud_optics_data, only  : cloud_optics_type
    use radiation_ice_optics_fu, only    : NIceOpticsCoeffsFuSW, &
         &                                 NIceOpticsCoeffsFuLW
    use radiation_ice_optics_baran, only : NIceOpticsCoeffsBaran, &
         &                                 NIceOpticsCoeffsBaran2016
    use radiation_ice_optics_baran2017, only : NIceOpticsCoeffsBaran2017, &
         &                                 NIceOpticsGeneralCoeffsBaran2017
    use radiation_ice_optics_yi, only    : NIceOpticsCoeffsYiSW, &
         &                                 NIceOpticsCoeffsYiLW
    use radiation_liquid_optics_socrates, only : NLiqOpticsCoeffsSOCRATES
    use radiation_liquid_optics_slingo, only : NLiqOpticsCoeffsSlingoSW, &
         &                                     NLiqOpticsCoeffsLindnerLiLW

    type(config_type), intent(inout) :: config

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud_optics:setup_cloud_optics',0,hook_handle)

    call config%cloud_optics%setup(trim(config%liq_optics_file_name), &
         &     trim(config%ice_optics_file_name), iverbose=config%iverbosesetup)

    ! Check liquid coefficients
    if (size(config%cloud_optics%liq_coeff_lw, 1) /= config%n_bands_lw) then
      write(nulerr,'(a,i0,a,i0,a)') &
           &  '*** Error: number of longwave bands for droplets (', &
           &  size(config%cloud_optics%liq_coeff_lw, 1), &
           &  ') does not match number for gases (', config%n_bands_lw, ')'
      call radiation_abort()
    end if
    if (size(config%cloud_optics%liq_coeff_sw, 1) /= config%n_bands_sw) then
      write(nulerr,'(a,i0,a,i0,a)') &
           &  '*** Error: number of shortwave bands for droplets (', &
           &  size(config%cloud_optics%liq_coeff_sw, 1), &
           &  ') does not match number for gases (', config%n_bands_sw, ')'
      call radiation_abort()        
    end if

    if (config%i_liq_model == ILiquidModelSOCRATES) then
      if (size(config%cloud_optics%liq_coeff_lw, 2) /= NLiqOpticsCoeffsSOCRATES) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** Error: number of liquid cloud optical coefficients (', &
             &  size(config%cloud_optics%liq_coeff_lw, 2), &
             &  ') does not match number expected (', NLiqOpticsCoeffsSOCRATES,')'
        call radiation_abort()
      end if
    else if (config%i_liq_model == ILiquidModelSlingo) then
      if (size(config%cloud_optics%liq_coeff_sw, 2) /= NLiqOpticsCoeffsSlingoSW) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** Error: number of shortwave liquid cloud optical coefficients (', &
             &  size(config%cloud_optics%liq_coeff_sw, 2), &
             &  ') does not match number expected (', NLiqOpticsCoeffsSlingoSW,')'
        call radiation_abort()
      end if
      if (size(config%cloud_optics%liq_coeff_lw, 2) /= NLiqOpticsCoeffsLindnerLiLW) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** Error: number of longwave liquid cloud optical coefficients (', &
             &  size(config%cloud_optics%liq_coeff_lw, 2), &
             &  ') does not match number expected (', NLiqOpticsCoeffsLindnerLiLw,')'
        call radiation_abort()
      end if
    end if

    ! Check ice coefficients
    if (size(config%cloud_optics%ice_coeff_lw, 1) /= config%n_bands_lw) then
      write(nulerr,'(a,i0,a,i0,a)') &
           &  '*** Error: number of longwave bands for ice particles (', &
           &  size(config%cloud_optics%ice_coeff_lw, 1), &
           &  ') does not match number for gases (', config%n_bands_lw, ')'
      call radiation_abort()
    end if
    if (size(config%cloud_optics%ice_coeff_sw, 1) /= config%n_bands_sw) then
      write(nulerr,'(a,i0,a,i0,a)') &
           &  '*** Error: number of shortwave bands for ice particles (', &
           &  size(config%cloud_optics%ice_coeff_sw, 1), &
           &  ') does not match number for gases (', config%n_bands_sw, ')'
      call radiation_abort()
    end if

    if (config%i_ice_model == IIceModelFu) then
      if (size(config%cloud_optics%ice_coeff_lw, 2) &
           &  /= NIceOpticsCoeffsFuLW) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** Error: number of LW ice-particle optical coefficients (', &
             &  size(config%cloud_optics%ice_coeff_lw, 2), &
             &  ') does not match number expected (', NIceOpticsCoeffsFuLW,')'
        call radiation_abort()
      end if
      if (size(config%cloud_optics%ice_coeff_sw, 2) &
           &  /= NIceOpticsCoeffsFuSW) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** Error: number of SW ice-particle optical coefficients (', &
             &  size(config%cloud_optics%ice_coeff_sw, 2), &
             &  ') does not match number expected (', NIceOpticsCoeffsFuSW,')'
        call radiation_abort()
      end if
    else if (config%i_ice_model == IIceModelBaran &
         &  .and. size(config%cloud_optics%ice_coeff_lw, 2) &
         &  /= NIceOpticsCoeffsBaran) then
      write(nulerr,'(a,i0,a,i0,a,i0,a)') &
           &  '*** Error: number of ice-particle optical coefficients (', &
           &  size(config%cloud_optics%ice_coeff_lw, 2), &
           &  ') does not match number expected (', NIceOpticsCoeffsBaran,')'
      call radiation_abort()
    else if (config%i_ice_model == IIceModelBaran2016 &
         &  .and. size(config%cloud_optics%ice_coeff_lw, 2) &
         &  /= NIceOpticsCoeffsBaran2016) then
      write(nulerr,'(a,i0,a,i0,a,i0,a)') &
           &  '*** Error: number of ice-particle optical coefficients (', &
           &  size(config%cloud_optics%ice_coeff_lw, 2), &
           &  ') does not match number expected (', NIceOpticsCoeffsBaran2016,')'
      call radiation_abort()
    else if (config%i_ice_model == IIceModelBaran2017) then
      if (size(config%cloud_optics%ice_coeff_lw, 2) &
           &  /= NIceOpticsCoeffsBaran2017) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** Error: number of ice-particle optical coefficients (', &
             &  size(config%cloud_optics%ice_coeff_lw, 2), &
             &  ') does not match number expected (', NIceOpticsCoeffsBaran2017,')'
        call radiation_abort()
      else if (.not. allocated(config%cloud_optics%ice_coeff_gen)) then
        write(nulerr,'(a)') &
             &  '*** Error: coeff_gen needed for Baran-2017 ice optics parameterization'
        call radiation_abort()
      else if (size(config%cloud_optics%ice_coeff_gen) &
           &  /= NIceOpticsGeneralCoeffsBaran2017) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** Error: number of general ice-particle optical coefficients (', &
             &  size(config%cloud_optics%ice_coeff_gen), &
             &  ') does not match number expected (', NIceOpticsGeneralCoeffsBaran2017,')'
        call radiation_abort()
      end if
    else if (config%i_ice_model == IIceModelYi) then
      if (size(config%cloud_optics%ice_coeff_lw, 2) &
           &  /= NIceOpticsCoeffsYiLW) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** Error: number of LW ice-particle optical coefficients (', &
             &  size(config%cloud_optics%ice_coeff_lw, 2), &
             &  ') does not match number expected (', NIceOpticsCoeffsYiLW,')'
        call radiation_abort()
      end if
      if (size(config%cloud_optics%ice_coeff_sw, 2) &
           &  /= NIceOpticsCoeffsYiSW) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** Error: number of SW ice-particle optical coefficients (', &
             &  size(config%cloud_optics%ice_coeff_sw, 2), &
             &  ') does not match number expected (', NIceOpticsCoeffsYiSW,')'
        call radiation_abort()
      end if
    end if

    if (lhook) call dr_hook('radiation_cloud_optics:setup_cloud_optics',1,hook_handle)

  end subroutine setup_cloud_optics


  !---------------------------------------------------------------------
  ! Compute cloud optical properties
  subroutine cloud_optics(nlev,istartcol,iendcol, &
       &  config, thermodynamics, cloud, & 
       &  od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
       &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    use radiation_io,     only : nulout, nulerr, radiation_abort
    use radiation_config, only : config_type, IIceModelFu, IIceModelBaran, &
         &                       IIceModelBaran2016, IIceModelBaran2017, &
         &                       IIceModelYi, &
         &                       ILiquidModelSOCRATES, ILiquidModelSlingo
    use radiation_thermodynamics, only    : thermodynamics_type
    use radiation_cloud, only             : cloud_type
    use radiation_constants, only         : AccelDueToGravity
    use radiation_cloud_optics_data, only : cloud_optics_type
    use radiation_ice_optics_fu, only     : calc_ice_optics_fu_sw, &
         &                                  calc_ice_optics_fu_lw
    use radiation_ice_optics_baran, only  : calc_ice_optics_baran, &
         &                                  calc_ice_optics_baran2016
    use radiation_ice_optics_baran2017, only  : calc_ice_optics_baran2017
    use radiation_ice_optics_yi, only     : calc_ice_optics_yi_sw, &
         &                                  calc_ice_optics_yi_lw
    use radiation_liquid_optics_socrates, only:calc_liq_optics_socrates
    use radiation_liquid_optics_slingo, only:calc_liq_optics_slingo, &
         &                                   calc_liq_optics_lindner_li

    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in), target :: config
    type(thermodynamics_type),intent(in)  :: thermodynamics
    type(cloud_type),   intent(in)        :: cloud

    ! Layer optical depth, single scattering albedo and g factor of
    ! clouds in each longwave band, where the latter two
    ! variables are only defined if cloud longwave scattering is
    ! enabled (otherwise both are treated as zero).
    real(jprb), dimension(config%n_bands_lw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_lw_cloud
    real(jprb), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol), &
         &   intent(out) :: ssa_lw_cloud, g_lw_cloud

    ! Layer optical depth, single scattering albedo and g factor of
    ! clouds in each shortwave band
    real(jprb), dimension(config%n_bands_sw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_sw_cloud, ssa_sw_cloud, g_sw_cloud

    ! Longwave and shortwave optical depth, scattering optical depth
    ! and asymmetry factor, for liquid and ice in all bands but a
    ! single cloud layer
    real(jprb), dimension(config%n_bands_lw) :: &
         &  od_lw_liq, scat_od_lw_liq, g_lw_liq, &
         &  od_lw_ice, scat_od_lw_ice, g_lw_ice
    real(jprb), dimension(config%n_bands_sw) :: &
         &  od_sw_liq, scat_od_sw_liq, g_sw_liq, &
         &  od_sw_ice, scat_od_sw_ice, g_sw_ice

    ! In-cloud water path of cloud liquid or ice (i.e. liquid or ice
    ! gridbox-mean water path divided by cloud fraction); kg m-2
    real(jprb) :: lwp_in_cloud, iwp_in_cloud

    ! Full-level temperature (K)
    real(jprb) :: temperature

    ! Factor to convert gridbox-mean mixing ratio to in-cloud water
    ! path
    real(jprb) :: factor

    ! Pointer to the cloud optics coefficients for brevity of
    ! access
    type(cloud_optics_type), pointer :: ho

    integer    :: jcol, jlev, jb

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud_optics:cloud_optics',0,hook_handle)

    if (config%iverbose >= 2) then
      write(nulout,'(a)') 'Computing cloud absorption/scattering properties'
    end if

    ho => config%cloud_optics

    ! Array-wise assignment
    od_lw_cloud  = 0.0_jprb
    od_sw_cloud  = 0.0_jprb
    ssa_sw_cloud = 0.0_jprb
    g_sw_cloud   = 0.0_jprb
    if (config%do_lw_cloud_scattering) then
      ssa_lw_cloud = 0.0_jprb
      g_lw_cloud   = 0.0_jprb
    end if

    do jlev = 1,nlev
      do jcol = istartcol,iendcol
        ! Only do anything if cloud is present (assume that
        ! cloud%crop_cloud_fraction has already been called)
        if (cloud%fraction(jcol,jlev) > 0.0_jprb) then

          ! Compute in-cloud liquid and ice water path
          if (config%is_homogeneous) then
            ! Homogeneous solvers assume cloud fills the box
            ! horizontally, so we don't divide by cloud fraction
            factor = ( thermodynamics%pressure_hl(jcol,jlev+1)    &
                 &  -thermodynamics%pressure_hl(jcol,jlev  )  ) &
                 &  / AccelDueToGravity
          else
            factor = ( thermodynamics%pressure_hl(jcol,jlev+1)    &
                 &  -thermodynamics%pressure_hl(jcol,jlev  )  ) &
                 &  / (AccelDueToGravity * cloud%fraction(jcol,jlev))
          end if
          lwp_in_cloud = factor * cloud%q_liq(jcol,jlev)
          iwp_in_cloud = factor * cloud%q_ice(jcol,jlev)

          ! Only compute liquid properties if liquid cloud is
          ! present
          if (lwp_in_cloud > 0.0_jprb) then
            if (config%i_liq_model == ILiquidModelSOCRATES) then
              ! Compute longwave properties
              call calc_liq_optics_socrates(config%n_bands_lw, &
                   &  config%cloud_optics%liq_coeff_lw, &
                   &  lwp_in_cloud, cloud%re_liq(jcol,jlev), &
                   &  od_lw_liq, scat_od_lw_liq, g_lw_liq)
              ! Compute shortwave properties
              call calc_liq_optics_socrates(config%n_bands_sw, &
                   &  config%cloud_optics%liq_coeff_sw, &
                   &  lwp_in_cloud, cloud%re_liq(jcol,jlev), &
                   &  od_sw_liq, scat_od_sw_liq, g_sw_liq)
            else if (config%i_liq_model == ILiquidModelSlingo) then
              ! Compute longwave properties
              call calc_liq_optics_lindner_li(config%n_bands_lw, &
                   &  config%cloud_optics%liq_coeff_lw, &
                   &  lwp_in_cloud, cloud%re_liq(jcol,jlev), &
                   &  od_lw_liq, scat_od_lw_liq, g_lw_liq)
              ! Compute shortwave properties
              call calc_liq_optics_slingo(config%n_bands_sw, &
                   &  config%cloud_optics%liq_coeff_sw, &
                   &  lwp_in_cloud, cloud%re_liq(jcol,jlev), &
                   &  od_sw_liq, scat_od_sw_liq, g_sw_liq)
            else
              write(nulerr,*) '*** Error: Unknown liquid model with code', &
                   &          config%i_liq_model
              call radiation_abort()
            end if

            ! Delta-Eddington scaling in the shortwave only
            if (.not. config%do_sw_delta_scaling_with_gases) then
              call delta_eddington_scat_od(od_sw_liq, scat_od_sw_liq, g_sw_liq)
            end if
            !call delta_eddington_scat_od(od_lw_liq, scat_od_lw_liq, g_lw_liq)

          else
            ! Liquid not present: set properties to zero
            od_lw_liq = 0.0_jprb
            scat_od_lw_liq = 0.0_jprb
            g_lw_liq = 0.0_jprb

            od_sw_liq = 0.0_jprb
            scat_od_sw_liq = 0.0_jprb
            g_sw_liq = 0.0_jprb
          end if ! Liquid present

          ! Only compute ice properties if ice cloud is present
          if (iwp_in_cloud > 0.0_jprb) then
            if (config%i_ice_model == IIceModelBaran) then
              ! Compute longwave properties
              call calc_ice_optics_baran(config%n_bands_lw, &
                   &  config%cloud_optics%ice_coeff_lw, &
                   &  iwp_in_cloud, cloud%q_ice(jcol,jlev), &
                   &  od_lw_ice, scat_od_lw_ice, g_lw_ice)
              ! Compute shortwave properties
              call calc_ice_optics_baran(config%n_bands_sw, &
                   &  config%cloud_optics%ice_coeff_sw, &
                   &  iwp_in_cloud, cloud%q_ice(jcol,jlev), &
                   &  od_sw_ice, scat_od_sw_ice, g_sw_ice)
            else if (config%i_ice_model == IIceModelBaran2016) then
              temperature = 0.5_jprb * (thermodynamics%temperature_hl(jcol,jlev) &
                   &                   +thermodynamics%temperature_hl(jcol,jlev+1))
              ! Compute longwave properties
              call calc_ice_optics_baran2016(config%n_bands_lw, &
                   &  config%cloud_optics%ice_coeff_lw, &
                   &  iwp_in_cloud, cloud%q_ice(jcol,jlev), &
                   &  temperature, &
                   &  od_lw_ice, scat_od_lw_ice, g_lw_ice)
              ! Compute shortwave properties
              call calc_ice_optics_baran2016(config%n_bands_sw, &
                   &  config%cloud_optics%ice_coeff_sw, &
                   &  iwp_in_cloud, cloud%q_ice(jcol,jlev), &
                   &  temperature, &
                   &  od_sw_ice, scat_od_sw_ice, g_sw_ice)
            else if (config%i_ice_model == IIceModelBaran2017) then
              temperature = 0.5_jprb * (thermodynamics%temperature_hl(jcol,jlev) &
                   &                   +thermodynamics%temperature_hl(jcol,jlev+1))
              ! Compute longwave properties
              call calc_ice_optics_baran2017(config%n_bands_lw, &
                   &  config%cloud_optics%ice_coeff_gen, &
                   &  config%cloud_optics%ice_coeff_lw, &
                   &  iwp_in_cloud, cloud%q_ice(jcol,jlev), &
                   &  temperature, &
                   &  od_lw_ice, scat_od_lw_ice, g_lw_ice)
              ! Compute shortwave properties
              call calc_ice_optics_baran2017(config%n_bands_sw, &
                   &  config%cloud_optics%ice_coeff_gen, &
                   &  config%cloud_optics%ice_coeff_sw, &
                   &  iwp_in_cloud, cloud%q_ice(jcol,jlev), &
                   &  temperature, &
                   &  od_sw_ice, scat_od_sw_ice, g_sw_ice)
            else if (config%i_ice_model == IIceModelFu) then
              ! Compute longwave properties
              call calc_ice_optics_fu_lw(config%n_bands_lw, &
                   &  config%cloud_optics%ice_coeff_lw, &
                   &  iwp_in_cloud, cloud%re_ice(jcol,jlev), &
                   &  od_lw_ice, scat_od_lw_ice, g_lw_ice)
              if (config%do_fu_lw_ice_optics_bug) then
                ! Reproduce bug in old IFS scheme
                scat_od_lw_ice = od_lw_ice - scat_od_lw_ice
              end if
              ! Compute shortwave properties
              call calc_ice_optics_fu_sw(config%n_bands_sw, &
                   &  config%cloud_optics%ice_coeff_sw, &
                   &  iwp_in_cloud, cloud%re_ice(jcol,jlev), &
                   &  od_sw_ice, scat_od_sw_ice, g_sw_ice)
            else if (config%i_ice_model == IIceModelYi) then
              ! Compute longwave properties
              call calc_ice_optics_yi_lw(config%n_bands_lw, &
                   &  config%cloud_optics%ice_coeff_lw, &
                   &  iwp_in_cloud, cloud%re_ice(jcol,jlev), &
                   &  od_lw_ice, scat_od_lw_ice, g_lw_ice)
              ! Compute shortwave properties
              call calc_ice_optics_yi_sw(config%n_bands_sw, &
                   &  config%cloud_optics%ice_coeff_sw, &
                   &  iwp_in_cloud, cloud%re_ice(jcol,jlev), &
                   &  od_sw_ice, scat_od_sw_ice, g_sw_ice)
            else
              write(nulerr,*) '*** Error: Unknown ice model with code', &
                   &          config%i_ice_model
              call radiation_abort()
            end if

            ! Delta-Eddington scaling in both longwave and shortwave
            ! (assume that particles are larger than wavelength even
            ! in longwave)
            if (.not. config%do_sw_delta_scaling_with_gases) then
              call delta_eddington_scat_od(od_sw_ice, scat_od_sw_ice, g_sw_ice)
            end if
            call delta_eddington_scat_od(od_lw_ice, scat_od_lw_ice, g_lw_ice)

          else
            ! Ice not present: set properties to zero
            od_lw_ice = 0.0_jprb
            scat_od_lw_ice = 0.0_jprb
            g_lw_ice = 0.0_jprb

            od_sw_ice = 0.0_jprb
            scat_od_sw_ice = 0.0_jprb
            g_sw_ice = 0.0_jprb
          end if ! Ice present

          ! Combine liquid and ice 
          if (config%do_lw_cloud_scattering) then
! Added for DWD (2020)
!NEC$ shortloop
            do jb = 1, config%n_bands_lw
              od_lw_cloud(jb,jlev,jcol) = od_lw_liq(jb) + od_lw_ice(jb)
              if (scat_od_lw_liq(jb)+scat_od_lw_ice(jb) > 0.0_jprb) then
                g_lw_cloud(jb,jlev,jcol) = (g_lw_liq(jb) * scat_od_lw_liq(jb) &
                   &  + g_lw_ice(jb) * scat_od_lw_ice(jb)) &
                   &  / (scat_od_lw_liq(jb)+scat_od_lw_ice(jb))
              else
                g_lw_cloud(jb,jlev,jcol) = 0.0_jprb
              end if
              ssa_lw_cloud(jb,jlev,jcol) = (scat_od_lw_liq(jb) + scat_od_lw_ice(jb)) &
                 &                    / (od_lw_liq(jb) + od_lw_ice(jb))
            end do
          else
            ! If longwave scattering is to be neglected then the
            ! best approximation is to set the optical depth equal
            ! to the absorption optical depth
! Added for DWD (2020)
!NEC$ shortloop
            do jb = 1, config%n_bands_lw
              od_lw_cloud(jb,jlev,jcol) = od_lw_liq(jb) - scat_od_lw_liq(jb) &
                    &                   + od_lw_ice(jb) - scat_od_lw_ice(jb)
            end do
          end if
! Added for DWD (2020)
!NEC$ shortloop
          do jb = 1, config%n_bands_sw
            od_sw_cloud(jb,jlev,jcol) = od_sw_liq(jb) + od_sw_ice(jb)
            g_sw_cloud(jb,jlev,jcol) = (g_sw_liq(jb) * scat_od_sw_liq(jb) &
               &  + g_sw_ice(jb) * scat_od_sw_ice(jb)) &
               &  / (scat_od_sw_liq(jb) + scat_od_sw_ice(jb))
            ssa_sw_cloud(jb,jlev,jcol) &
               &  = (scat_od_sw_liq(jb) + scat_od_sw_ice(jb)) / (od_sw_liq(jb) + od_sw_ice(jb))
          end do
        end if ! Cloud present
      end do ! Loop over column
    end do ! Loop over level

    if (lhook) call dr_hook('radiation_cloud_optics:cloud_optics',1,hook_handle)

  end subroutine cloud_optics

end module radiation_cloud_optics
