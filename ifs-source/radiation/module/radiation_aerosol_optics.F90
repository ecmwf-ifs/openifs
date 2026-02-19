! radiation_aerosol_optics.F90 - Computing aerosol optical properties
!
! Copyright (C) 2015-2018 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! Modifications
!   2018-04-15  R Hogan  Add "direct" option

module radiation_aerosol_optics

  implicit none

contains

  ! Provides the elemental function "delta_eddington_extensive"
#include "radiation_delta_eddington.h"

  !---------------------------------------------------------------------
  ! Load aerosol scattering data; this subroutine delegates to one
  ! in radiation_aerosol_optics_data.F90
  subroutine setup_aerosol_optics(config)

    use parkind1,                      only : jprb
    use yomhook,  only           : lhook, dr_hook, jphook
    use radiation_config,              only : config_type
    use radiation_aerosol_optics_data, only : aerosol_optics_type
    use radiation_io,                  only : nulerr, radiation_abort

    type(config_type), intent(inout) :: config

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics:setup_aerosol_optics',0,hook_handle)

    if (config%n_aerosol_types > 0) then
      ! Load data from file and prepare to map config%n_aerosol_types
      ! aerosol types
      call config%aerosol_optics%setup(trim(config%aerosol_optics_file_name), &
           &                        config%n_aerosol_types, iverbose=config%iverbosesetup)

      ! Check agreement in number of bands
      if (config%n_bands_lw /= config%aerosol_optics%n_bands_lw) then
        write(nulerr,'(a)') '*** Error: number of longwave bands does not match aerosol optics look-up table'
        call radiation_abort()
      end if
      if (config%n_bands_sw /= config%aerosol_optics%n_bands_sw) then
        write(nulerr,'(a)') '*** Error: number of shortwave bands does not match aerosol optics look-up table'
        call radiation_abort()
      end if

      ! Map aerosol types to those loaded from the data file
      call config%aerosol_optics%set_types(config%i_aerosol_type_map(1:config%n_aerosol_types))
    end if

    call config%aerosol_optics%print_description(config%i_aerosol_type_map(1:config%n_aerosol_types))

    if (lhook) call dr_hook('radiation_aerosol_optics:setup_aerosol_optics',1,hook_handle)

  end subroutine setup_aerosol_optics


  !---------------------------------------------------------------------
  ! Compute aerosol optical properties and add to existing gas optical
  ! depth and scattering properties
  subroutine add_aerosol_optics(nlev,istartcol,iendcol, &
       &  config, thermodynamics, gas, aerosol, & 
       &  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)

    use parkind1,                      only : jprb
    use radiation_io,                  only : nulout, nulerr, radiation_abort
    use yomhook,  only           : lhook, dr_hook, jphook
    use radiation_config,              only : config_type
    use radiation_thermodynamics,      only : thermodynamics_type
    use radiation_gas,                 only : gas_type, IH2O, IMassMixingRatio
    use radiation_aerosol,             only : aerosol_type
    use radiation_constants,           only : AccelDueToGravity
    use radiation_aerosol_optics_data, only : aerosol_optics_type, &
         &  IAerosolClassUndefined,   IAerosolClassIgnored, &
         &  IAerosolClassHydrophobic, IAerosolClassHydrophilic

    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in), target :: config
    type(thermodynamics_type),intent(in)  :: thermodynamics
    type(gas_type),           intent(in)  :: gas
    type(aerosol_type),       intent(in)  :: aerosol
    ! Optical depth, single scattering albedo and asymmetry factor of
    ! the atmosphere (gases on input, gases and aerosols on output)
    ! for each g point. Note that longwave ssa and asymmetry and
    ! shortwave asymmetry are all zero for gases, so are not yet
    ! defined on input and are therefore intent(out).
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), &
         &   intent(inout) :: od_lw
    real(jprb), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol), &
         &   intent(out)   :: ssa_lw, g_lw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), &
         &   intent(inout) :: od_sw, ssa_sw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), &
         &   intent(out)   :: g_sw

    ! Extinction optical depth, scattering optical depth and
    ! asymmetry-times-scattering-optical-depth for all the aerosols at
    ! a point in space for each spectral band of the shortwave and
    ! longwave spectrum
    real(jprb), dimension(config%n_bands_sw) &
         & :: od_sw_aerosol, scat_sw_aerosol, scat_g_sw_aerosol, local_od_sw
    real(jprb), dimension(config%n_bands_lw) :: od_lw_aerosol, local_od_lw
    real(jprb), dimension(config%n_bands_lw_if_scattering) &
         & :: scat_lw_aerosol, scat_g_lw_aerosol

    real(jprb) :: h2o_mmr(istartcol:iendcol,nlev)

    real(jprb) :: rh ! Relative humidity with respect to liquid water

    ! Factor (kg m-2) to convert mixing ratio (kg kg-1) to mass in
    ! path (kg m-2)
    real(jprb) :: factor

    ! Temporary extinction and scattering optical depths of aerosol
    ! plus gas
    real(jprb) :: local_od, local_scat

    ! Loop indices for column, level, g point, band and aerosol type
    integer :: jcol, jlev, jg, jtype

    ! Range of levels over which aerosols are present
    integer :: istartlev, iendlev

    ! Indices to spectral band and relative humidity look-up table
    integer :: iband, irh

    ! Pointer to the aerosol optics coefficients for brevity of access
    type(aerosol_optics_type), pointer :: ao

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics:add_aerosol_optics',0,hook_handle)

    if (aerosol%is_direct) then
      ! Aerosol optical properties have been provided in each band
      ! directly by the user
      call add_aerosol_optics_direct(nlev,istartcol,iendcol, &
           &  config, aerosol, & 
           &  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)
    else
      ! Aerosol mixing ratios have been provided

      do jtype = 1,config%n_aerosol_types
        if (config%aerosol_optics%iclass(jtype) == IAerosolClassUndefined) then
          write(nulerr,'(a)') '*** Error: not all aerosol types are defined'
          call radiation_abort()
        end if
      end do

      if (config%iverbose >= 2) then
        write(nulout,'(a)') 'Computing aerosol absorption/scattering properties'
      end if

      ao => config%aerosol_optics

      istartlev = lbound(aerosol%mixing_ratio,2)
      iendlev   = ubound(aerosol%mixing_ratio,2)

      if (ubound(aerosol%mixing_ratio,3) /= config%n_aerosol_types) then
        write(nulerr,'(a,i0,a,i0)') '*** Error: aerosol%mixing_ratio contains ', &
             &  ubound(aerosol%mixing_ratio,3), ' aerosol types, expected ', &
             &  config%n_aerosol_types
        call radiation_abort()
      end if

      ! Set variables to zero that may not have been previously
      g_sw = 0.0_jprb
      if (config%do_lw_aerosol_scattering) then
        ssa_lw = 0.0_jprb
        g_lw   = 0.0_jprb
      end if

      call gas%get(IH2O, IMassMixingRatio, h2o_mmr, istartcol=istartcol)

      ! Loop over position
      do jlev = istartlev,iendlev
        do jcol = istartcol,iendcol
          ! Compute relative humidity with respect to liquid
          ! saturation and the index to the relative-humidity index of
          ! hydrophilic-aerosol data
          rh  = h2o_mmr(jcol,jlev) / thermodynamics%h2o_sat_liq(jcol,jlev)
          irh = ao%calc_rh_index(rh)

          factor = ( thermodynamics%pressure_hl(jcol,jlev+1) &
               &    -thermodynamics%pressure_hl(jcol,jlev  )  ) &
               &   / AccelDueToGravity  

          ! Reset temporary arrays
          od_sw_aerosol     = 0.0_jprb
          scat_sw_aerosol   = 0.0_jprb
          scat_g_sw_aerosol = 0.0_jprb
          od_lw_aerosol     = 0.0_jprb
          scat_lw_aerosol   = 0.0_jprb
          scat_g_lw_aerosol = 0.0_jprb

          do jtype = 1,config%n_aerosol_types
            ! Add the optical depth, scattering optical depth and
            ! scattering optical depth-weighted asymmetry factor for
            ! this aerosol type to the total for all aerosols.  Note
            ! that the following expressions are array-wise, the
            ! dimension being spectral band.
            if (ao%iclass(jtype) == IAerosolClassHydrophobic) then
              local_od_sw = factor * aerosol%mixing_ratio(jcol,jlev,jtype) &
                   &  * ao%mass_ext_sw_phobic(:,ao%itype(jtype))
              od_sw_aerosol = od_sw_aerosol + local_od_sw
              scat_sw_aerosol = scat_sw_aerosol &
                   &  + local_od_sw * ao%ssa_sw_phobic(:,ao%itype(jtype))
              scat_g_sw_aerosol = scat_g_sw_aerosol &
                   &  + local_od_sw * ao%ssa_sw_phobic(:,ao%itype(jtype)) &
                   &  * ao%g_sw_phobic(:,ao%itype(jtype))
              if (config%do_lw_aerosol_scattering) then
                local_od_lw = factor * aerosol%mixing_ratio(jcol,jlev,jtype) &
                     &  * ao%mass_ext_lw_phobic(:,ao%itype(jtype))
                od_lw_aerosol = od_lw_aerosol + local_od_lw
                scat_lw_aerosol = scat_lw_aerosol &
                     &  + local_od_lw * ao%ssa_lw_phobic(:,ao%itype(jtype))
                scat_g_lw_aerosol = scat_g_lw_aerosol &
                     &  + local_od_lw * ao%ssa_lw_phobic(:,ao%itype(jtype)) &
                     &  * ao%g_lw_phobic(:,ao%itype(jtype))
              else
                ! If aerosol longwave scattering is not included then we
                ! weight the optical depth by the single scattering
                ! co-albedo
                od_lw_aerosol = od_lw_aerosol &
                     &  + factor * aerosol%mixing_ratio(jcol,jlev,jtype) &
                     &  * ao%mass_ext_lw_phobic(:,ao%itype(jtype)) &
                     &  * (1.0_jprb - ao%ssa_lw_phobic(:,ao%itype(jtype)))
              end if
            else if (ao%iclass(jtype) == IAerosolClassHydrophilic) then
              ! Hydrophilic aerosols require the look-up tables to
              ! be indexed with irh
              local_od_sw = factor * aerosol%mixing_ratio(jcol,jlev,jtype) &
                   &  * ao%mass_ext_sw_philic(:,irh,ao%itype(jtype))
              od_sw_aerosol = od_sw_aerosol + local_od_sw
              scat_sw_aerosol = scat_sw_aerosol &
                   &  + local_od_sw * ao%ssa_sw_philic(:,irh,ao%itype(jtype))
              scat_g_sw_aerosol = scat_g_sw_aerosol &
                   &  + local_od_sw * ao%ssa_sw_philic(:,irh,ao%itype(jtype)) &
                   &  * ao%g_sw_philic(:,irh,ao%itype(jtype))
              if (config%do_lw_aerosol_scattering) then
                local_od_lw = factor * aerosol%mixing_ratio(jcol,jlev,jtype) &
                     &  * ao%mass_ext_lw_philic(:,irh,ao%itype(jtype))
                od_lw_aerosol = od_lw_aerosol + local_od_lw
                scat_lw_aerosol = scat_lw_aerosol &
                     &  + local_od_lw * ao%ssa_lw_philic(:,irh,ao%itype(jtype))
                scat_g_lw_aerosol = scat_g_lw_aerosol &
                     &  + local_od_lw * ao%ssa_lw_philic(:,irh,ao%itype(jtype)) &
                     &  * ao%g_lw_philic(:,irh,ao%itype(jtype))
              else
                ! If aerosol longwave scattering is not included then we
                ! weight the optical depth by the single scattering
                ! co-albedo
                od_lw_aerosol = od_lw_aerosol &
                     &  + factor * aerosol%mixing_ratio(jcol,jlev,jtype) &
                     &  * ao%mass_ext_lw_philic(:,irh,ao%itype(jtype)) &
                     &  * (1.0_jprb - ao%ssa_lw_philic(:,irh,ao%itype(jtype)))
              end if
            end if
            ! Implicitly, if ao%iclass(jtype) == IAerosolClassNone, then
            ! no aerosol scattering properties are added

          end do ! Loop over aerosol type

          if (.not. config%do_sw_delta_scaling_with_gases) then
            ! Delta-Eddington scaling on aerosol only.  Note that if
            ! do_sw_delta_scaling_with_gases==.true. then the delta
            ! scaling is done to the cloud-aerosol-gas mixture inside
            ! the solver
            call delta_eddington_extensive(od_sw_aerosol, scat_sw_aerosol, &
                 &                         scat_g_sw_aerosol)
          end if

          ! Combine aerosol shortwave scattering properties with gas
          ! properties (noting that any gas scattering will have an
          ! asymmetry factor of zero)
          if (od_sw_aerosol(1) > 0.0_jprb) then
            do jg = 1,config%n_g_sw
              iband = config%i_band_from_reordered_g_sw(jg)
              local_od = od_sw(jg,jlev,jcol) + od_sw_aerosol(iband)
              local_scat = ssa_sw(jg,jlev,jcol) * od_sw(jg,jlev,jcol) &
                   &  + scat_sw_aerosol(iband)
              ! Note that asymmetry_sw of gases is zero so the following
              ! simply weights the aerosol asymmetry by the scattering
              ! optical depth
              g_sw(jg,jlev,jcol) = scat_g_sw_aerosol(iband) / local_scat
              ssa_sw(jg,jlev,jcol) = local_scat / local_od
              od_sw (jg,jlev,jcol) = local_od
            end do
          end if

          ! Combine aerosol longwave scattering properties with gas
          ! properties, noting that in the longwave, gases do not
          ! scatter at all
          if (config%do_lw_aerosol_scattering) then

            call delta_eddington_extensive(od_lw_aerosol, scat_lw_aerosol, &
                 &                         scat_g_lw_aerosol)  

            do jg = 1,config%n_g_lw
              iband = config%i_band_from_reordered_g_lw(jg)
              if (od_lw_aerosol(iband) > 0.0_jprb) then
                ! All scattering is due to aerosols, therefore the
                ! asymmetry factor is equal to the value for aerosols
                if (scat_lw_aerosol(iband) > 0.0_jprb) then
                  g_lw(jg,jlev,jcol) = scat_g_lw_aerosol(iband) &
                       &  / scat_lw_aerosol(iband)
                else
                  g_lw(jg,jlev,jcol) = 0.0_jprb
                end if
                local_od = od_lw(jg,jlev,jcol) + od_lw_aerosol(iband)
                ssa_lw(jg,jlev,jcol) = scat_lw_aerosol(iband) / local_od
                od_lw (jg,jlev,jcol) = local_od
              end if
            end do
          else
            do jg = 1,config%n_g_lw
              od_lw(jg,jlev,jcol) = od_lw(jg,jlev,jcol) &
                   &  + od_lw_aerosol(config%i_band_from_reordered_g_lw(jg))
            end do
          end if

        end do ! Loop over column
      end do ! Loop over level

    end if

    if (lhook) call dr_hook('radiation_aerosol_optics:add_aerosol_optics',1,hook_handle)

  end subroutine add_aerosol_optics


  !---------------------------------------------------------------------
  ! Add precomputed optical properties to gas optical depth and
  ! scattering properties
  subroutine add_aerosol_optics_direct(nlev,istartcol,iendcol, &
       &  config, aerosol, & 
       &  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)

    use parkind1,                      only : jprb
    use radiation_io,                  only : nulerr, radiation_abort
    use yomhook,  only           : lhook, dr_hook, jphook
    use radiation_config,              only : config_type
    use radiation_aerosol,             only : aerosol_type

    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in), target :: config
    type(aerosol_type),       intent(in)  :: aerosol
    ! Optical depth, single scattering albedo and asymmetry factor of
    ! the atmosphere (gases on input, gases and aerosols on output)
    ! for each g point. Note that longwave ssa and asymmetry and
    ! shortwave asymmetry are all zero for gases, so are not yet
    ! defined on input and are therefore intent(out).
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), &
         &   intent(inout) :: od_lw
    real(jprb), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol), &
         &   intent(out)   :: ssa_lw, g_lw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), &
         &   intent(inout) :: od_sw, ssa_sw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), &
         &   intent(out)   :: g_sw

    ! Temporary extinction and scattering optical depths of aerosol
    ! plus gas
    real(jprb) :: local_od, local_scat

    ! Extinction optical depth, scattering optical depth and
    ! asymmetry-times-scattering-optical-depth for all the aerosols at
    ! a point in space for each spectral band of the shortwave and
    ! longwave spectrum
    real(jprb), dimension(config%n_bands_sw) &
         & :: od_sw_aerosol, scat_sw_aerosol, scat_g_sw_aerosol
    real(jprb), dimension(config%n_bands_lw) :: od_lw_aerosol
    real(jprb), dimension(config%n_bands_lw_if_scattering) &
         & :: scat_lw_aerosol, scat_g_lw_aerosol

    ! Loop indices for column, level, g point and band
    integer :: jcol, jlev, jg

    ! Range of levels over which aerosols are present
    integer :: istartlev, iendlev

    ! Indices to spectral band
    integer :: iband

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics:add_aerosol_optics_direct',0,hook_handle)

    if (config%do_sw) then
      ! Check array dimensions
      if (ubound(aerosol%od_sw,1) /= config%n_bands_sw) then
        write(nulerr,'(a,i0,a,i0)') '*** Error: aerosol%od_sw contains ', &
           &  ubound(aerosol%od_sw,1), ' band, expected ', &
           &  config%n_bands_sw
        call radiation_abort()
      end if

      istartlev = lbound(aerosol%od_sw,2)
      iendlev   = ubound(aerosol%od_sw,2)

      ! Set variables to zero that may not have been previously
      g_sw = 0.0_jprb

      ! Loop over position
      do jcol = istartcol,iendcol
        do jlev = istartlev,iendlev
          od_sw_aerosol = aerosol%od_sw(:,jlev,jcol)
          scat_sw_aerosol = aerosol%ssa_sw(:,jlev,jcol) * od_sw_aerosol
          scat_g_sw_aerosol = aerosol%g_sw(:,jlev,jcol) * scat_sw_aerosol

          if (.not. config%do_sw_delta_scaling_with_gases) then
            ! Delta-Eddington scaling on aerosol only.  Note that if
            ! do_sw_delta_scaling_with_gases==.true. then the delta
            ! scaling is done to the cloud-aerosol-gas mixture inside
            ! the solver
            call delta_eddington_extensive(od_sw_aerosol, scat_sw_aerosol, &
                 &                         scat_g_sw_aerosol)
          end if

          ! Combine aerosol shortwave scattering properties with gas
          ! properties (noting that any gas scattering will have an
          ! asymmetry factor of zero)
          if (od_sw_aerosol(1) > 0.0_jprb) then
            do jg = 1,config%n_g_sw
              iband = config%i_band_from_reordered_g_sw(jg)
              local_od = od_sw(jg,jlev,jcol) + od_sw_aerosol(iband)
              local_scat = ssa_sw(jg,jlev,jcol) * od_sw(jg,jlev,jcol) &
                   &  + scat_sw_aerosol(iband)
              ! Note that asymmetry_sw of gases is zero so the following
              ! simply weights the aerosol asymmetry by the scattering
              ! optical depth
              g_sw(jg,jlev,jcol) = scat_g_sw_aerosol(iband) / local_scat
              ssa_sw(jg,jlev,jcol) = local_scat / local_od
              od_sw (jg,jlev,jcol) = local_od
            end do
          end if
        end do
      end do

    end if

    if (config%do_lw) then

      if (ubound(aerosol%od_lw,1) /= config%n_bands_lw) then
        write(nulerr,'(a,i0,a,i0)') '*** Error: aerosol%od_lw contains ', &
           &  ubound(aerosol%od_lw,1), ' band, expected ', &
           &  config%n_bands_lw
        call radiation_abort()
      end if

      istartlev = lbound(aerosol%od_lw,2)
      iendlev   = ubound(aerosol%od_lw,2)

      if (config%do_lw_aerosol_scattering) then
        ssa_lw = 0.0_jprb
        g_lw   = 0.0_jprb
 
        ! Loop over position
        do jcol = istartcol,iendcol
          do jlev = istartlev,iendlev
            od_lw_aerosol = aerosol%od_lw(:,jlev,jcol)
            scat_lw_aerosol = aerosol%ssa_lw(:,jlev,jcol) * od_lw_aerosol
            scat_g_lw_aerosol = aerosol%g_lw(:,jlev,jcol) * scat_lw_aerosol
            
            call delta_eddington_extensive(od_lw_aerosol, scat_lw_aerosol, &
                 &                         scat_g_lw_aerosol)
            
            do jg = 1,config%n_g_lw
              iband = config%i_band_from_reordered_g_lw(jg)
              if (od_lw_aerosol(iband) > 0.0_jprb) then
                ! All scattering is due to aerosols, therefore the
                ! asymmetry factor is equal to the value for aerosols
                if (scat_lw_aerosol(iband) > 0.0_jprb) then
                  g_lw(jg,jlev,jcol) = scat_g_lw_aerosol(iband) &
                       &  / scat_lw_aerosol(iband)
                else
                  g_lw(jg,jlev,jcol) = 0.0_jprb
                end if
                local_od = od_lw(jg,jlev,jcol) + od_lw_aerosol(iband)
                ssa_lw(jg,jlev,jcol) = scat_lw_aerosol(iband) / local_od
                od_lw (jg,jlev,jcol) = local_od
              end if
            end do
          end do
        end do

      else ! No longwave scattering

        ! Loop over position
        do jcol = istartcol,iendcol
          do jlev = istartlev,iendlev
            ! If aerosol longwave scattering is not included then we
            ! weight the optical depth by the single scattering
            ! co-albedo
            od_lw_aerosol = aerosol%od_lw(:,jlev,jcol) &
                 &  * (1.0_jprb - aerosol%ssa_lw(:,jlev,jcol))
            do jg = 1,config%n_g_lw
              od_lw(jg,jlev,jcol) = od_lw(jg,jlev,jcol) &
                   &  + od_lw_aerosol(config%i_band_from_reordered_g_lw(jg))
            end do
          end do
        end do

      end if
    end if


    if (lhook) call dr_hook('radiation_aerosol_optics:add_aerosol_optics_direct',1,hook_handle)

  end subroutine add_aerosol_optics_direct
 

  !---------------------------------------------------------------------
  ! Sometimes it is useful to specify aerosol in terms of its optical
  ! depth at a particular wavelength.  This function returns the dry
  ! shortwave mass-extinction coefficient, i.e. the extinction cross
  ! section per unit mass, for aerosol of type "itype" at shortwave
  ! band "iband". For hydrophilic types, the value at the first
  ! relative humidity bin is taken.
  function dry_aerosol_sw_mass_extinction(config, itype, iband)

    use parkind1,                      only : jprb
    use radiation_config,              only : config_type
    use radiation_aerosol_optics_data, only : aerosol_optics_type, &
         &  IAerosolClassUndefined,   IAerosolClassIgnored, &
         &  IAerosolClassHydrophobic, IAerosolClassHydrophilic

    type(config_type), intent(in), target :: config

    ! Aerosol type and shortwave band as indices to the array
    integer, intent(in) :: itype, iband
    
    real(jprb) dry_aerosol_sw_mass_extinction

    ! Pointer to the aerosol optics coefficients for brevity of access
    type(aerosol_optics_type), pointer :: ao

    ao => config%aerosol_optics

    if (ao%iclass(itype) == IAerosolClassHydrophobic) then
      dry_aerosol_sw_mass_extinction = ao%mass_ext_sw_phobic(iband,ao%itype(itype))
    else if (ao%iclass(itype) == IAerosolClassHydrophilic) then
      ! Take the value at the first relative-humidity bin for the
      ! "dry" aerosol value
      dry_aerosol_sw_mass_extinction = ao%mass_ext_sw_philic(iband,1,ao%itype(itype))
    else
      dry_aerosol_sw_mass_extinction = 0.0_jprb
    end if

  end function dry_aerosol_sw_mass_extinction


  !---------------------------------------------------------------------
  ! Compute aerosol extinction coefficient at a particular shortwave
  ! band and a single height - this is useful for visibility
  ! diagnostics
  subroutine aerosol_sw_extinction(ncol,istartcol,iendcol, &
       &  config, iband, mixing_ratio, relative_humidity, extinction)

    use parkind1,                      only : jprb
    use yomhook,  only           : lhook, dr_hook, jphook
    use radiation_io,                  only : nulerr, radiation_abort
    use radiation_config,              only : config_type
    use radiation_aerosol_optics_data, only : aerosol_optics_type, &
         &  IAerosolClassUndefined,   IAerosolClassIgnored, &
         &  IAerosolClassHydrophobic, IAerosolClassHydrophilic

    integer, intent(in) :: ncol               ! number of columns
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in), target :: config
    integer, intent(in)     :: iband ! Index of required spectral band
    real(jprb), intent(in)  :: mixing_ratio(ncol,config%n_aerosol_types)
    real(jprb), intent(in)  :: relative_humidity(ncol)
    real(jprb), intent(out) :: extinction(ncol)

    ! Local aerosol extinction
    real(jprb) :: ext

    ! Pointer to the aerosol optics coefficients for brevity of access
    type(aerosol_optics_type), pointer :: ao
    
    ! Loop indices for column and aerosol type
    integer :: jcol, jtype

    ! Relative humidity index
    integer :: irh

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics:aerosol_sw_extinction',0,hook_handle)

    do jtype = 1,config%n_aerosol_types
      if (config%aerosol_optics%iclass(jtype) == IAerosolClassUndefined) then
        write(nulerr,'(a)') '*** Error: not all aerosol types are defined'
        call radiation_abort()
      end if
    end do

    ao => config%aerosol_optics

    ! Loop over position
    do jcol = istartcol,iendcol
      ext = 0.0_jprb
      ! Get relative-humidity index
      irh = ao%calc_rh_index(relative_humidity(jcol))
      ! Add extinction coefficients from each aerosol type
      do jtype = 1,config%n_aerosol_types
        if (ao%iclass(jtype) == IAerosolClassHydrophobic) then
          ext = ext + mixing_ratio(jcol,jtype) &
               &    * ao%mass_ext_sw_phobic(iband,ao%itype(jtype))
        else if (ao%iclass(jtype) == IAerosolClassHydrophilic) then
          ext = ext + mixing_ratio(jcol,jtype) &
               &    * ao%mass_ext_sw_philic(iband,irh,ao%itype(jtype))
        end if
      end do

      extinction(jcol) = ext
    end do

    if (lhook) call dr_hook('radiation_aerosol_optics:aerosol_sw_extinction',1,hook_handle)

  end subroutine aerosol_sw_extinction

end module radiation_aerosol_optics
