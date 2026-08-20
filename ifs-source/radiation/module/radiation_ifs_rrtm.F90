! radiation_ifs_rrtm.F90 - Interface to IFS implementation of RRTM-G
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
! Modifications
!   2017-04-11  R. Hogan  Receive "surface" dummy argument
!   2017-09-08  R. Hogan  Reverted some changes
!   2017-10-18  R. Hogan  Added planck_function public function
!   2018-01-11  R. Hogan  Added optional spectral scaling of incoming solar radiation
!   2018-02-22  R. Hogan  Optimized reverse indexing of heights
!   2018-05-05  R. Hogan  gas_optics can be called for reduced number of levels
!   2019-01-02  R. Hogan  Initialize shortwave props to zero in case sun below horizon

module radiation_ifs_rrtm

  implicit none

  public  :: setup_gas_optics, gas_optics, planck_function, set_gas_units

contains

  !---------------------------------------------------------------------
  ! Setup the IFS implementation of RRTM-G gas absorption model
  subroutine setup_gas_optics(config, directory)

    use yoerrtm,   only : jpglw
    use yoesrtm,   only : jpgsw
    use yoerrtftr, only : ngb_lw => ngb
    use yoesrtm,   only : ngb_sw => ngbsw
    use yomhook,   only : lhook, dr_hook, jphook

    use radiation_config

    type(config_type), intent(inout), target :: config
    character(len=*), intent(in)     :: directory

    integer :: irep ! For implied do

    integer, parameter :: RRTM_GPOINT_REORDERING_LW(140) = (/ &
          &   89, 90, 139, 77, 137, 69, 131, 97, 91, 70, 78, 71, 53, 72, 123, 54, 79, 98,  &
          &   92, 55, 80, 132, 124, 81, 73, 56, 99, 82, 57, 23, 125, 100, 24, 74, 93, 58, 25,  &
          &   83, 126, 75, 26, 11, 101, 133, 59, 27, 76, 140, 12, 84, 102, 94, 28, 127, 85,  &
          &   13, 39, 60, 86, 103, 87, 109, 14, 29, 115, 40, 95, 15, 61, 88, 41, 110, 104, 1,  &
          &   116, 42, 30, 134, 128, 138, 96, 62, 16, 43, 117, 63, 111, 44, 2, 64, 31, 65,  &
          &   105, 17, 45, 66, 118, 32, 3, 33, 67, 18, 129, 135, 46, 112, 34, 106, 68, 35, 4,  &
          &   119, 36, 47, 107, 19, 37, 38, 113, 48, 130, 5, 120, 49, 108, 20, 50, 51, 114,  &
          &   21, 121, 52, 136, 122, 6, 22, 7, 8, 9, 10 &
          & /)
    integer, parameter :: RRTM_GPOINT_REORDERING_SW(112) = (/ &
          &   35, 45, 19, 27, 36, 57, 20, 46, 58, 21, 28, 67, 55, 68, 37, 1, 69, 22, 29, 59,  &
          &   78, 101, 79, 77, 70, 76, 47, 75, 30, 81, 60, 102, 80, 82, 23, 2, 83, 84, 85,  &
          &   86, 103, 61, 31, 87, 56, 38, 71, 48, 88, 3, 62, 89, 24, 7, 49, 32, 104, 72, 90,  &
          &   63, 39, 4, 8, 50, 91, 64, 40, 33, 25, 51, 95, 96, 73, 65, 9, 41, 97, 92, 105,  &
          &   52, 5, 98, 10, 42, 99, 100, 66, 11, 74, 34, 53, 26, 6, 106, 12, 43, 13, 54, 93,  &
          &   44, 107, 94, 14, 108, 15, 16, 109, 17, 18, 110, 111, 112 &
          & /)
    
    real(jphook) :: hook_handle

!#include "surdi.intfb.h"
#include "surrtab.intfb.h"
#include "surrtpk.intfb.h"
#include "surrtrf.intfb.h"
#include "rrtm_init_140gp.intfb.h"
#include "srtm_init.intfb.h"

    if (lhook) call dr_hook('radiation_ifs_rrtm:setup_gas_optics',0,hook_handle)

    ! The IFS implementation of RRTMG uses many global variables.  In
    ! the IFS these will have been set up already; otherwise set them
    ! up now.
    if (config%do_setup_ifsrrtm) then
      !call SURDI
      call SURRTAB
      call SURRTPK
      call SURRTRF
      call RRTM_INIT_140GP(directory)
      call SRTM_INIT(directory)
    end if

    ! Cloud and aerosol properties can only be defined per band
    config%do_cloud_aerosol_per_sw_g_point = .false.
    config%do_cloud_aerosol_per_lw_g_point = .false.

    config%n_g_sw = jpgsw
    config%n_g_lw = jpglw
    config%n_bands_sw = 14
    config%n_bands_lw = 16

    ! Wavenumber ranges of each band may be needed so that the user
    ! can compute UV and photosynthetically active radiation for a
    ! particular wavelength range
    call config%gas_optics_sw%spectral_def%allocate_bands_only( &
         &  [2600.0_jprb, 3250.0_jprb, 4000.0_jprb, 4650.0_jprb, 5150.0_jprb, 6150.0_jprb, 7700.0_jprb, &
         &   8050.0_jprb, 12850.0_jprb, 16000.0_jprb, 22650.0_jprb, 29000.0_jprb, 38000.0_jprb, 820.0_jprb], &
         &  [3250.0_jprb, 4000.0_jprb, 4650.0_jprb, 5150.0_jprb, 6150.0_jprb, 7700.0_jprb, 8050.0_jprb, &
         &   12850.0_jprb, 16000.0_jprb, 22650.0_jprb, 29000.0_jprb, 38000.0_jprb, 50000.0_jprb, 2600.0_jprb])
    call config%gas_optics_lw%spectral_def%allocate_bands_only( &
         &  [10.0_jprb, 350.0_jprb, 500.0_jprb, 630.0_jprb, 700.0_jprb, 820.0_jprb, 980.0_jprb, 1080.0_jprb, &
         &   1180.0_jprb, 1390.0_jprb, 1480.0_jprb, 1800.0_jprb, 2080.0_jprb, 2250.0_jprb, 2380.0_jprb, 2600.0_jprb], &
         &  [350.0_jprb, 500.0_jprb, 630.0_jprb, 700.0_jprb, 820.0_jprb, 980.0_jprb, 1080.0_jprb, 1180.0_jprb, &
         &   1390.0_jprb, 1480.0_jprb, 1800.0_jprb, 2080.0_jprb, 2250.0_jprb, 2380.0_jprb, 2600.0_jprb, 3250.0_jprb])

    allocate(config%i_band_from_g_sw          (config%n_g_sw))
    allocate(config%i_band_from_g_lw          (config%n_g_lw))
    allocate(config%i_band_from_reordered_g_sw(config%n_g_sw))
    allocate(config%i_band_from_reordered_g_lw(config%n_g_lw))
    allocate(config%i_g_from_reordered_g_sw(config%n_g_sw))
    allocate(config%i_g_from_reordered_g_lw(config%n_g_lw))

    ! Shortwave starts at 16: need to start at 1
    config%i_band_from_g_sw = ngb_sw - ngb_sw(1)+1
    config%i_band_from_g_lw = ngb_lw

    if (config%i_solver_sw == ISolverSpartacus) then
      ! SPARTACUS requires g points ordered in approximately
      ! increasing order of optical depth
      config%i_g_from_reordered_g_sw = RRTM_GPOINT_REORDERING_SW
    else
      ! Implied-do for no reordering
!      config%i_g_from_reordered_g_sw = RRTM_GPOINT_REORDERING_SW
      config%i_g_from_reordered_g_sw = (/ (irep, irep=1,config%n_g_sw) /)
    end if

    if (config%i_solver_lw == ISolverSpartacus) then
      ! SPARTACUS requires g points ordered in approximately
      ! increasing order of optical depth
      config%i_g_from_reordered_g_lw = RRTM_GPOINT_REORDERING_LW
    else
      ! Implied-do for no reordering
      config%i_g_from_reordered_g_lw = (/ (irep, irep=1,config%n_g_lw) /)
    end if

    config%i_band_from_reordered_g_sw &
         = config%i_band_from_g_sw(config%i_g_from_reordered_g_sw)

    config%i_band_from_reordered_g_lw &
         = config%i_band_from_g_lw(config%i_g_from_reordered_g_lw)

    ! The i_spec_* variables are used solely for storing spectral
    ! data, and this can either be by band or by g-point
    if (config%do_save_spectral_flux .or. config%do_toa_spectral_flux) then
      if (config%do_save_gpoint_flux) then
        config%n_spec_sw = config%n_g_sw
        config%n_spec_lw = config%n_g_lw
        config%i_spec_from_reordered_g_sw => config%i_g_from_reordered_g_sw
        config%i_spec_from_reordered_g_lw => config%i_g_from_reordered_g_lw
      else
        config%n_spec_sw = config%n_bands_sw
        config%n_spec_lw = config%n_bands_lw
        config%i_spec_from_reordered_g_sw => config%i_band_from_reordered_g_sw
        config%i_spec_from_reordered_g_lw => config%i_band_from_reordered_g_lw
      end if
    else
      config%n_spec_sw = 0
      config%n_spec_lw = 0
      nullify(config%i_spec_from_reordered_g_sw)
      nullify(config%i_spec_from_reordered_g_lw)
    end if

    if (lhook) call dr_hook('radiation_ifs_rrtm:setup_gas_optics',1,hook_handle)

  end subroutine setup_gas_optics


  !---------------------------------------------------------------------
  ! Scale gas mixing ratios according to required units
  subroutine set_gas_units(gas)

    use radiation_gas,           only : gas_type, IMassMixingRatio
    type(gas_type),    intent(inout) :: gas

    call gas%set_units(IMassMixingRatio)

  end subroutine set_gas_units


  !---------------------------------------------------------------------
  ! Compute gas optical depths, shortwave scattering, Planck function
  ! and incoming shortwave radiation at top-of-atmosphere
  subroutine gas_optics(ncol,nlev,istartcol,iendcol, &
       &  config, single_level, thermodynamics, gas, & 
       &  od_lw, od_sw, ssa_sw, lw_albedo, planck_hl, lw_emission, &
       &  incoming_sw)

    use parkind1,                 only : jprb, jpim

    USE PARRRTM  , ONLY : JPBAND, JPXSEC, JPINPX 
    USE YOERRTM  , ONLY : JPGPT_LW => JPGPT
    USE YOESRTM  , ONLY : JPGPT_SW => JPGPT  
    use yomhook  , only : lhook, dr_hook, jphook

    use radiation_config,         only : config_type, ISolverSpartacus
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_single_level,   only : single_level_type
    use radiation_gas

    integer, intent(in) :: ncol               ! number of columns
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(gas_type),           intent(in) :: gas

    ! Longwave albedo of the surface
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
         &  intent(in), optional :: lw_albedo

    ! Gaseous layer optical depth in longwave and shortwave, and
    ! shortwave single scattering albedo (i.e. fraction of extinction
    ! due to Rayleigh scattering) at each g-point
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_lw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_sw, ssa_sw

    ! The Planck function (emitted flux from a black body) at half
    ! levels at each longwave g-point
    real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), &
         &   intent(out), optional :: planck_hl
    ! Planck function for the surface (W m-2)
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
         &   intent(out), optional :: lw_emission

    ! The incoming shortwave flux into a plane perpendicular to the
    ! incoming radiation at top-of-atmosphere in each of the shortwave
    ! g-points
    real(jprb), dimension(config%n_g_sw,istartcol:iendcol), &
         &   intent(out), optional :: incoming_sw

    real(jprb) :: incoming_sw_scale(istartcol:iendcol)

    ! The variables in capitals are used in the same way as the
    ! equivalent routine in the IFS

    real(jprb) :: ZOD_LW(JPGPT_LW,nlev,istartcol:iendcol) ! Note ordering of dimensions
    real(jprb) :: ZOD_SW(istartcol:iendcol,nlev,JPGPT_SW)
    real(jprb) :: ZSSA_SW(istartcol:iendcol,nlev,JPGPT_SW)
    real(jprb) :: ZINCSOL(istartcol:iendcol,JPGPT_SW)

    real(jprb) :: ZCOLMOL(istartcol:iendcol,nlev)
    real(jprb) :: ZCOLDRY(istartcol:iendcol,nlev)
    real(jprb) :: ZWBRODL(istartcol:iendcol,nlev) !BROADENING GASES,column density (mol/cm2)
    real(jprb) :: ZCOLBRD(istartcol:iendcol,nlev) !BROADENING GASES, column amount
    real(jprb) :: ZWKL(istartcol:iendcol,JPINPX,nlev)

    real(jprb) :: ZWX(istartcol:iendcol,JPXSEC,nlev) ! Amount of trace gases
    
    real(jprb) :: ZFLUXFAC, ZPI

    ! - from AER
    real(jprb) :: ZTAUAERL(istartcol:iendcol,nlev,JPBAND)

    !- from INTFAC      
    real(jprb) :: ZFAC00(istartcol:iendcol,nlev)
    real(jprb) :: ZFAC01(istartcol:iendcol,nlev)
    real(jprb) :: ZFAC10(istartcol:iendcol,nlev)
    real(jprb) :: ZFAC11(istartcol:iendcol,nlev)
    
    !- from FOR
    real(jprb) :: ZFORFAC(istartcol:iendcol,nlev)
    real(jprb) :: ZFORFRAC(istartcol:iendcol,nlev)
    integer    :: INDFOR(istartcol:iendcol,nlev) 

    !- from MINOR
    integer    :: INDMINOR(istartcol:iendcol,nlev) 
    real(jprb) :: ZSCALEMINOR(istartcol:iendcol,nlev) 
    real(jprb) :: ZSCALEMINORN2(istartcol:iendcol,nlev) 
    real(jprb) :: ZMINORFRAC(istartcol:iendcol,nlev) 
    
    real(jprb)     :: &                 
         &  ZRAT_H2OCO2(istartcol:iendcol,nlev),ZRAT_H2OCO2_1(istartcol:iendcol,nlev), &
         &  ZRAT_H2OO3(istartcol:iendcol,nlev) ,ZRAT_H2OO3_1(istartcol:iendcol,nlev), & 
         &  ZRAT_H2ON2O(istartcol:iendcol,nlev),ZRAT_H2ON2O_1(istartcol:iendcol,nlev), &
         &  ZRAT_H2OCH4(istartcol:iendcol,nlev),ZRAT_H2OCH4_1(istartcol:iendcol,nlev), &
         &  ZRAT_N2OCO2(istartcol:iendcol,nlev),ZRAT_N2OCO2_1(istartcol:iendcol,nlev), &
         &  ZRAT_O3CO2(istartcol:iendcol,nlev) ,ZRAT_O3CO2_1(istartcol:iendcol,nlev)
    
    !- from INTIND
    integer :: JP(istartcol:iendcol,nlev)
    integer :: JT(istartcol:iendcol,nlev)
    integer :: JT1(istartcol:iendcol,nlev)

    !- from PRECISE             
    real(jprb) :: ZONEMINUS, ZONEMINUS_ARRAY(istartcol:iendcol)

    !- from PROFDATA             
    real(jprb) :: ZCOLH2O(istartcol:iendcol,nlev)
    real(jprb) :: ZCOLCO2(istartcol:iendcol,nlev)
    real(jprb) :: ZCOLO3(istartcol:iendcol,nlev)
    real(jprb) :: ZCOLN2O(istartcol:iendcol,nlev)
    real(jprb) :: ZCOLCH4(istartcol:iendcol,nlev)
    real(jprb) :: ZCOLO2(istartcol:iendcol,nlev)
    real(jprb) :: ZCO2MULT(istartcol:iendcol,nlev)
    integer    :: ILAYTROP(istartcol:iendcol)
    integer    :: ILAYSWTCH(istartcol:iendcol)
    integer    :: ILAYLOW(istartcol:iendcol)

    !- from PROFILE             
    real(jprb) :: ZPAVEL(istartcol:iendcol,nlev)
    real(jprb) :: ZTAVEL(istartcol:iendcol,nlev)
    real(jprb) :: ZPZ(istartcol:iendcol,0:nlev)
    real(jprb) :: ZTZ(istartcol:iendcol,0:nlev)
    
    !- from SELF             
    real(jprb) :: ZSELFFAC(istartcol:iendcol,nlev)
    real(jprb) :: ZSELFFRAC(istartcol:iendcol,nlev)
    integer :: INDSELF(istartcol:iendcol,nlev)

    !- from SP             
    real(jprb) :: ZPFRAC(istartcol:iendcol,JPGPT_LW,nlev)
    
    !- from SURFACE             
    integer :: IREFLECT(istartcol:iendcol)

    real(jprb) :: pressure_fl(ncol, nlev), temperature_fl(ncol, nlev)

    ! If nlev is less than the number of heights at which gas mixing
    ! ratios are stored, then we assume that the lower part of the
    ! atmosphere is required. This enables nlev=1 to be passed in to
    ! the routine, in which case the gas properties of the lowest
    ! layer are provided, useful for canopy radiative transfer.
    integer :: istartlev, iendlev

    integer :: jlev, jgreorder, jg, ig, iband, jcol

    real(jphook) :: hook_handle

#include "rrtm_prepare_gases.intfb.h"
#include "rrtm_setcoef_140gp.intfb.h"
#include "rrtm_gas_optical_depth.intfb.h"
#include "srtm_setcoef.intfb.h"
#include "srtm_gas_optical_depth.intfb.h"

    if (lhook) call dr_hook('radiation_ifs_rrtm:gas_optics',0,hook_handle)

    ! Compute start and end levels for indexing the gas mixing ratio
    ! and thermodynamics arrays
    iendlev   = ubound(gas%mixing_ratio,2)
    istartlev = iendlev - nlev + 1

    ZPI = 2.0_jprb*ASIN(1.0_jprb)
    ZFLUXFAC = ZPI * 1.E+4
    ZONEMINUS = 1.0_jprb - 1.0e-6_jprb
    ZONEMINUS_ARRAY = ZONEMINUS

    do jlev=1,nlev
      do jcol= istartcol,iendcol
        pressure_fl(jcol,jlev) &
            &  = 0.5_jprb * (thermodynamics%pressure_hl(jcol,jlev+istartlev-1) &
            &               +thermodynamics%pressure_hl(jcol,jlev+istartlev))
        temperature_fl(jcol,jlev) &
            &  = 0.5_jprb * (thermodynamics%temperature_hl(jcol,jlev+istartlev-1) &
            &               +thermodynamics%temperature_hl(jcol,jlev+istartlev))
      end do
    end do
    
    ! Check we have gas mixing ratios in the right units
    call gas%assert_units(IMassMixingRatio)

    ! Warning: O2 is hard-coded within the following function so the
    ! user-provided concentrations of this gas are ignored for both
    ! the longwave and shortwave
    CALL RRTM_PREPARE_GASES &
         & ( istartcol, iendcol, ncol, nlev, &
         &   thermodynamics%pressure_hl(:,istartlev:iendlev+1), &
         &   pressure_fl, &
         &   thermodynamics%temperature_hl(:,istartlev:iendlev+1), &
         &   temperature_fl, &
         &   gas%mixing_ratio(:,istartlev:iendlev,IH2O), &
         &   gas%mixing_ratio(:,istartlev:iendlev,ICO2), &
         &   gas%mixing_ratio(:,istartlev:iendlev,ICH4), &
         &   gas%mixing_ratio(:,istartlev:iendlev,IN2O), &
         &   gas%mixing_ratio(:,istartlev:iendlev,INO2), &
         &   gas%mixing_ratio(:,istartlev:iendlev,ICFC11), &
         &   gas%mixing_ratio(:,istartlev:iendlev,ICFC12), &
         &   gas%mixing_ratio(:,istartlev:iendlev,IHCFC22), &
         &   gas%mixing_ratio(:,istartlev:iendlev,ICCl4), &
         &   gas%mixing_ratio(:,istartlev:iendlev,IO3), &
         &  ZCOLDRY, ZWBRODL,ZWKL, ZWX, &
         &  ZPAVEL , ZTAVEL , ZPZ , ZTZ, IREFLECT)  

    CALL RRTM_SETCOEF_140GP &
         &( istartcol, iendcol, nlev , ZCOLDRY  , ZWBRODL , ZWKL , &
         &  ZFAC00 , ZFAC01   , ZFAC10 , ZFAC11 , ZFORFAC,ZFORFRAC,INDFOR, JP, JT, JT1 , &
         &  ZCOLH2O, ZCOLCO2  , ZCOLO3 , ZCOLN2O, ZCOLCH4, ZCOLO2,ZCO2MULT , ZCOLBRD, & 
         &  ILAYTROP,ILAYSWTCH, ILAYLOW, ZPAVEL , ZTAVEL , ZSELFFAC, ZSELFFRAC, INDSELF, &
         &  INDMINOR,ZSCALEMINOR,ZSCALEMINORN2,ZMINORFRAC,&
         &  ZRAT_H2OCO2, ZRAT_H2OCO2_1, ZRAT_H2OO3, ZRAT_H2OO3_1, &
         &  ZRAT_H2ON2O, ZRAT_H2ON2O_1, ZRAT_H2OCH4, ZRAT_H2OCH4_1, &
         &  ZRAT_N2OCO2, ZRAT_N2OCO2_1, ZRAT_O3CO2, ZRAT_O3CO2_1)   

    ZTAUAERL(istartcol:iendcol,:,:) = 0.0_jprb

    CALL RRTM_GAS_OPTICAL_DEPTH &
         &( istartcol, iendcol, nlev, ZOD_LW, ZPAVEL, ZCOLDRY, ZCOLBRD, ZWX ,&
         &  ZTAUAERL, ZFAC00 , ZFAC01, ZFAC10 , ZFAC11 , ZFORFAC,ZFORFRAC,INDFOR, &
         &  JP, JT, JT1, ZONEMINUS ,&
         &  ZCOLH2O , ZCOLCO2, ZCOLO3, ZCOLN2O, ZCOLCH4, ZCOLO2,ZCO2MULT ,&
         &  ILAYTROP, ILAYSWTCH,ILAYLOW, ZSELFFAC, ZSELFFRAC, INDSELF, ZPFRAC, &
         &  INDMINOR,ZSCALEMINOR,ZSCALEMINORN2,ZMINORFRAC,&
         &  ZRAT_H2OCO2, ZRAT_H2OCO2_1, ZRAT_H2OO3, ZRAT_H2OO3_1, &
         &  ZRAT_H2ON2O, ZRAT_H2ON2O_1, ZRAT_H2OCH4, ZRAT_H2OCH4_1, &
         &  ZRAT_N2OCO2, ZRAT_N2OCO2_1, ZRAT_O3CO2, ZRAT_O3CO2_1)      

    if (present(lw_albedo)) then
    
      call planck_function_atmos(nlev, istartcol, iendcol, config, &
           &                     thermodynamics, ZPFRAC, planck_hl)

      if (single_level%is_simple_surface) then
        call planck_function_surf(istartcol, iendcol, config, &
             &                    single_level%skin_temperature, ZPFRAC(:,:,1), &
             &                    lw_emission)
        
        ! The following can be used to extract the parameters defined at
        ! the top of the planck_function routine below:
        !write(*,'(a,140(e12.5,","),a)') 'ZPFRAC_surf=[', &
        !&  sum(ZPFRAC(istartcol:iendcol,:,1),1) / (iendcol+1-istartcol), ']'
        
        ! lw_emission at this point is actually the planck function of
        ! the surface
        lw_emission = lw_emission * (1.0_jprb - lw_albedo)
      else
        ! Longwave emission has already been computed
        if (config%use_canopy_full_spectrum_lw) then
          lw_emission = transpose(single_level%lw_emission(istartcol:iendcol,:))
        else
          lw_emission = transpose(single_level%lw_emission(istartcol:iendcol, &
               & config%i_emiss_from_band_lw(config%i_band_from_reordered_g_lw)))
        end if
      end if

    end if

    if (config%i_solver_lw == ISolverSpartacus) then
      ! We need to rearrange the gas optics info in memory: reordering
      ! the g points in order of approximately increasing optical
      ! depth (for efficient 3D processing on only the regions of the
      ! spectrum that are optically thin for gases) and reorder in
      ! pressure since the the functions above treat pressure
      ! decreasing with increasing index.  Note that the output gas
      ! arrays have dimensions in a different order to the inputs,
      ! so there is some inefficiency here.
      do jgreorder = 1,config%n_g_lw
        iband = config%i_band_from_reordered_g_lw(jgreorder)
        ig = config%i_g_from_reordered_g_lw(jgreorder)
        
        ! Top-of-atmosphere half level
        do jlev = 1,nlev
          do jcol = istartcol,iendcol
            ! Some g points can return negative optical depths;
            ! specifically original g points 54-56 which causes
            ! unphysical single-scattering albedo when combined with
            ! aerosol
            od_lw(jgreorder,jlev,jcol) &
                 &   = max(config%min_gas_od_lw, ZOD_LW(ig,nlev+1-jlev,jcol))
          end do
        end do
      end do
    else
      ! G points have not been reordered 
      do jcol = istartcol,iendcol
        do jlev = 1,nlev
          ! Check for negative optical depth
          od_lw(:,jlev,jcol) = max(config%min_gas_od_lw, ZOD_LW(:,nlev+1-jlev,jcol))
        end do
      end do
    end if
    
    CALL SRTM_SETCOEF &
         & ( istartcol, iendcol, nlev,&
         & ZPAVEL  , ZTAVEL,&
         & ZCOLDRY , ZWKL,&
         & ILAYTROP,&
         & ZCOLCH4  , ZCOLCO2 , ZCOLH2O , ZCOLMOL  , ZCOLO2 , ZCOLO3,&
         & ZFORFAC , ZFORFRAC , INDFOR  , ZSELFFAC, ZSELFFRAC, INDSELF, &
         & ZFAC00  , ZFAC01   , ZFAC10  , ZFAC11,&
         & JP      , JT       , JT1     , single_level%cos_sza(istartcol:iendcol)  &
         & )  
    
    ! SRTM_GAS_OPTICAL_DEPTH will not initialize profiles when the sun
    ! is below the horizon, so we do it here
    ZOD_SW(istartcol:iendcol,:,:)  = 0.0_jprb
    ZSSA_SW(istartcol:iendcol,:,:) = 0.0_jprb
    ZINCSOL(istartcol:iendcol,:)   = 0.0_jprb

    CALL SRTM_GAS_OPTICAL_DEPTH &
         &( istartcol, iendcol , nlev  , ZONEMINUS_ARRAY,&
         & single_level%cos_sza(istartcol:iendcol), ILAYTROP,&
         & ZCOLCH4 , ZCOLCO2  , ZCOLH2O, ZCOLMOL , ZCOLO2   , ZCOLO3,&
         & ZFORFAC , ZFORFRAC , INDFOR , ZSELFFAC, ZSELFFRAC, INDSELF,&
         & ZFAC00  , ZFAC01   , ZFAC10 , ZFAC11  ,&
         & JP      , JT       , JT1    ,&
         & ZOD_SW  , ZSSA_SW  , ZINCSOL )
    
    ! Scale the incoming solar per band, if requested
    if (config%use_spectral_solar_scaling) then
      do jg = 1,JPGPT_SW
        do jcol = istartcol,iendcol 
          ZINCSOL(jcol,jg) = ZINCSOL(jcol,jg) * &
            &   single_level%spectral_solar_scaling(config%i_band_from_reordered_g_sw(jg))
        end do
      end do
    end if

    ! Scaling factor to ensure that the total solar irradiance is as
    ! requested.  Note that if the sun is below the horizon then
    ! ZINCSOL will be zero.
    if (present(incoming_sw)) then
      do jcol = istartcol,iendcol
        if (single_level%cos_sza(jcol) > 0.0_jprb) then
! Added for DWD (2020)
!NEC$ nounroll
          incoming_sw_scale(jcol) = single_level%solar_irradiance / sum(ZINCSOL(jcol,:))
        else
          incoming_sw_scale(jcol) = 1.0_jprb
        end if
      end do
    end if

    if (config%i_solver_sw == ISolverSpartacus) then
!    if (.true.) then
      ! Account for reordered g points
      do jgreorder = 1,config%n_g_sw
        ig = config%i_g_from_reordered_g_sw(jgreorder)
        do jlev = 1,nlev
          do jcol = istartcol,iendcol
            ! Check for negative optical depth
            od_sw (jgreorder,nlev+1-jlev,jcol) &
                 &  = max(config%min_gas_od_sw, ZOD_SW (jcol,jlev,ig))
            ssa_sw(jgreorder,nlev+1-jlev,jcol) = ZSSA_SW(jcol,jlev,ig)
           end do
        end do
        if (present(incoming_sw)) then
          incoming_sw(jgreorder,:) &
               &  = incoming_sw_scale(:) * ZINCSOL(:,ig)
        end if
      end do
    else
      ! G points have not been reordered
      do jcol = istartcol,iendcol
        do jlev = 1,nlev
          do jg = 1,config%n_g_sw
            ! Check for negative optical depth
            od_sw (jg,nlev+1-jlev,jcol) = max(config%min_gas_od_sw, ZOD_SW(jcol,jlev,jg))
            ssa_sw(jg,nlev+1-jlev,jcol) = ZSSA_SW(jcol,jlev,jg)
          end do
        end do
        if (present(incoming_sw)) then
          do jg = 1,config%n_g_sw
            incoming_sw(jg,jcol) = incoming_sw_scale(jcol) * ZINCSOL(jcol,jg)
          end do
        end if
      end do
    end if
    
    if (lhook) call dr_hook('radiation_ifs_rrtm:gas_optics',1,hook_handle)
    
  end subroutine gas_optics
  

  !---------------------------------------------------------------------
  ! Compute Planck function of the atmosphere
  subroutine planck_function_atmos(nlev,istartcol,iendcol, &
       config, thermodynamics, PFRAC, &
       planck_hl)

    use parkind1,                 only : jprb, jpim

    USE YOERRTM  , ONLY : JPGPT_LW => JPGPT
    use yoerrtwn, only : totplnk, delwave

    use yomhook, only : lhook, dr_hook, jphook

    use radiation_config,         only : config_type, ISolverSpartacus
    use radiation_thermodynamics, only : thermodynamics_type
    !use radiation_gas

    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in) :: config
    type(thermodynamics_type),intent(in) :: thermodynamics
    real(jprb), intent(in) :: PFRAC(istartcol:iendcol,JPGPT_LW,nlev)

    ! The Planck function (emitted flux from a black body) at half
    ! levels at each longwave g-point
    real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), intent(out) :: &
         &   planck_hl

    ! Planck function values per band
    real(jprb), dimension(istartcol:iendcol, config%n_bands_lw) :: planck_store

    ! Look-up table variables for Planck function
    real(jprb), dimension(istartcol:iendcol) :: frac
    integer,    dimension(istartcol:iendcol) :: ind

    ! Temperature (K) of a half-level
    real(jprb) :: temperature

    real(jprb) :: factor, planck_tmp(istartcol:iendcol,config%n_g_lw)
    real(jprb) :: ZFLUXFAC

    integer :: jlev, jgreorder, jg, ig, iband, jband, jcol, ilevoffset

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ifs_rrtm:planck_function_atmos',0,hook_handle)

    ZFLUXFAC = 2.0_jprb*ASIN(1.0_jprb) * 1.0e4_jprb
    
    ! nlev may be less than the number of original levels, in which
    ! case we assume that the user wants the lower part of the
    ! atmosphere
    ilevoffset = ubound(thermodynamics%temperature_hl,2)-nlev-1

    ! Work out interpolations: for each half level, the index of the
    ! lowest interpolation bound, and the fraction into interpolation
    ! interval
    do jlev = 1,nlev+1
      do jcol = istartcol,iendcol
        temperature = thermodynamics%temperature_hl(jcol,jlev+ilevoffset)
        if (temperature < 339.0_jprb .and. temperature >= 160.0_jprb) then
          ! Linear interpolation between -113 and 66 degC
          ind(jcol)  = int(temperature - 159.0_jprb)
          frac(jcol) = temperature - int(temperature)
        else if(temperature >= 339.0_jprb) then
          ! Extrapolation above 66 degC
          ind(jcol)  = 180
          frac(jcol) = temperature - 339.0_jprb
        else
          ! Cap below -113 degC (to avoid possible negative Planck
          ! function values)
          ind(jcol)  = 1
          frac(jcol) = 0.0_jprb
        end if
      end do

      ! Calculate Planck functions per band
      do jband = 1,config%n_bands_lw
        factor = zfluxfac * delwave(jband)
        do jcol = istartcol,iendcol
          planck_store(jcol,jband) = factor &
               &  * (totplnk(ind(jcol),jband) &
               &  + frac(jcol)*(totplnk(ind(jcol)+1,jband)-totplnk(ind(jcol),jband)))
        end do
      end do

      if (config%i_solver_lw == ISolverSpartacus) then
        ! We need to rearrange the gas optics info in memory:
        ! reordering the g points in order of approximately increasing
        ! optical depth (for efficient 3D processing on only the
        ! regions of the spectrum that are optically thin for gases)
        ! and reorder in pressure since the the functions above treat
        ! pressure decreasing with increasing index.
        if (jlev == 1) then
          ! Top-of-atmosphere half level - note that PFRAC is on model
          ! levels not half levels
          do jgreorder = 1,config%n_g_lw
            iband = config%i_band_from_reordered_g_lw(jgreorder)
            ig = config%i_g_from_reordered_g_lw(jgreorder)
            planck_hl(jgreorder,1,:) = planck_store(:,iband) &
                 &   * PFRAC(:,ig,nlev)
          end do
        else
          do jgreorder = 1,config%n_g_lw
            iband = config%i_band_from_reordered_g_lw(jgreorder)
            ig = config%i_g_from_reordered_g_lw(jgreorder)
            planck_hl(jgreorder,jlev,:) &
                   &   = planck_store(:,iband) &
                   &   * PFRAC(:,ig,nlev+2-jlev)
          end do
        end if
      else
        ! G points have not been reordered 
        if (jlev == 1) then
          ! Top-of-atmosphere half level - note that PFRAC is on model
          ! levels not half levels
          do jg = 1,config%n_g_lw
            iband = config%i_band_from_g_lw(jg)
            planck_hl(jg,1,:) = planck_store(:,iband) * PFRAC(:,jg,nlev)
          end do
        else
          do jg = 1,config%n_g_lw
            iband = config%i_band_from_g_lw(jg)
            planck_tmp(:,jg) = planck_store(:,iband) * PFRAC(:,jg,nlev+2-jlev)
          end do
          do jcol = istartcol,iendcol
            planck_hl(:,jlev,jcol) = planck_tmp(jcol,:)
          end do
        end if
      end if
    end do

    if (lhook) call dr_hook('radiation_ifs_rrtm:planck_function_atmos',1,hook_handle)

  end subroutine planck_function_atmos


  !---------------------------------------------------------------------
  ! Compute Planck function of the surface
  subroutine planck_function_surf(istartcol, iendcol, config, temperature, PFRAC, &
       &  planck_surf)

    use parkind1,                 only : jprb, jpim

    USE YOERRTM  , ONLY : JPGPT_LW => JPGPT
    use yoerrtwn, only : totplnk, delwave

    use yomhook, only : lhook, dr_hook, jphook

    use radiation_config,         only : config_type, ISolverSpartacus
    !    use radiation_gas

    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in) :: config
    real(jprb), intent(in) :: temperature(:)

    real(jprb), intent(in) :: PFRAC(istartcol:iendcol,JPGPT_LW)

    ! Planck function of the surface (W m-2)
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
         &  intent(out) :: planck_surf

    ! Planck function values per band
    real(jprb), dimension(istartcol:iendcol, config%n_bands_lw) :: planck_store

    ! Look-up table variables for Planck function
    real(jprb), dimension(istartcol:iendcol) :: frac
    integer,    dimension(istartcol:iendcol) :: ind

    ! Temperature (K)
    real(jprb) :: Tsurf

    real(jprb) :: factor
    real(jprb) :: ZFLUXFAC

    integer :: jgreorder, jg, ig, iband, jband, jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ifs_rrtm:planck_function_surf',0,hook_handle)

    ZFLUXFAC = 2.0_jprb*ASIN(1.0_jprb) * 1.0e4_jprb

    ! Work out surface interpolations
    do jcol = istartcol,iendcol
      Tsurf = temperature(jcol)
      if (Tsurf < 339.0_jprb .and. Tsurf >= 160.0_jprb) then
        ! Linear interpolation between -113 and 66 degC
        ind(jcol)  = int(Tsurf - 159.0_jprb)
        frac(jcol) = Tsurf - int(Tsurf)
      else if(Tsurf >= 339.0_jprb) then
        ! Extrapolation above 66 degC
        ind(jcol)  = 180
        frac(jcol) = Tsurf - 339.0_jprb
      else
        ! Cap below -113 degC (to avoid possible negative Planck
        ! function values)
        ind(jcol)  = 1
        frac(jcol) = 0.0_jprb
      end if
    end do

    ! Calculate Planck functions per band
    do jband = 1,config%n_bands_lw
      factor = zfluxfac * delwave(jband)
      do jcol = istartcol,iendcol
        planck_store(jcol,jband) = factor &
             &  * (totplnk(ind(jcol),jband) &
             &  + frac(jcol)*(totplnk(ind(jcol)+1,jband)-totplnk(ind(jcol),jband)))
      end do
    end do

    if (config%i_solver_lw == ISolverSpartacus) then
      ! We need to rearrange the gas optics info in memory: reordering
      ! the g points in order of approximately increasing optical
      ! depth (for efficient 3D processing on only the regions of
      ! the spectrum that are optically thin for gases) and reorder
      ! in pressure since the the functions above treat pressure
      ! decreasing with increasing index.
      do jgreorder = 1,config%n_g_lw
        iband = config%i_band_from_reordered_g_lw(jgreorder)
        ig = config%i_g_from_reordered_g_lw(jgreorder)
        planck_surf(jgreorder,:) = planck_store(:,iband) * PFRAC(:,ig)
      end do
    else
      ! G points have not been reordered 
      do jg = 1,config%n_g_lw
        iband = config%i_band_from_g_lw(jg)
        planck_surf(jg,:) = planck_store(:,iband) * PFRAC(:,jg)
      end do
    end if

    if (lhook) call dr_hook('radiation_ifs_rrtm:planck_function_surf',1,hook_handle)
    
  end subroutine planck_function_surf


  !---------------------------------------------------------------------
  ! Externally facing function for computing the Planck function
  ! without reference to any gas profile; typically this would be used
  ! for computing the emission by facets of a complex surface.  Note
  ! that this uses fixed "PFRAC" values, obtained by averaging over
  ! those derived from RRTM-G for near-surface conditions over a line
  ! of meridian from the ECMWF model.
  subroutine planck_function(config, temperature, planck_surf)

    use parkind1,                 only : jprb, jpim

    use radiation_config,         only : config_type

    type(config_type), intent(in) :: config
    real(jprb), intent(in) :: temperature

    ! Planck function of the surface (W m-2)
    real(jprb), dimension(config%n_g_lw), &
         &  intent(out) :: planck_surf

    ! Fraction of each band contributed by each g-point within
    ! it. Since there are 16 bands, this array sums to 16
    real(jprb), parameter, dimension(1,140) :: frac &
         = reshape( (/ 0.21227E+00, 0.18897E+00, 0.25491E+00, 0.17864E+00, 0.11735E+00, 0.38298E-01, 0.57871E-02, &
         &    0.31753E-02, 0.53169E-03, 0.76476E-04, 0.16388E+00, 0.15241E+00, 0.14290E+00, 0.12864E+00, &
         &    0.11615E+00, 0.10047E+00, 0.80013E-01, 0.60445E-01, 0.44918E-01, 0.63395E-02, 0.32942E-02, &
         &    0.54541E-03, 0.15380E+00, 0.15194E+00, 0.14339E+00, 0.13138E+00, 0.11701E+00, 0.10081E+00, &
         &    0.82296E-01, 0.61735E-01, 0.41918E-01, 0.45918E-02, 0.37743E-02, 0.30121E-02, 0.22500E-02, &
         &    0.14490E-02, 0.55410E-03, 0.78364E-04, 0.15938E+00, 0.15146E+00, 0.14213E+00, 0.13079E+00, &
         &    0.11672E+00, 0.10053E+00, 0.81566E-01, 0.61126E-01, 0.41150E-01, 0.44488E-02, 0.36950E-02, &
         &    0.29101E-02, 0.21357E-02, 0.19609E-02, 0.14134E+00, 0.14390E+00, 0.13913E+00, 0.13246E+00, &
         &    0.12185E+00, 0.10596E+00, 0.87518E-01, 0.66164E-01, 0.44862E-01, 0.49402E-02, 0.40857E-02, &
         &    0.32288E-02, 0.23613E-02, 0.15406E-02, 0.58258E-03, 0.82171E-04, 0.29127E+00, 0.28252E+00, &
         &    0.22590E+00, 0.14314E+00, 0.45494E-01, 0.71792E-02, 0.38483E-02, 0.65712E-03, 0.29810E+00, &
         &    0.27559E+00, 0.11997E+00, 0.10351E+00, 0.84515E-01, 0.62253E-01, 0.41050E-01, 0.44217E-02, &
         &    0.36946E-02, 0.29113E-02, 0.34290E-02, 0.55993E-03, 0.31441E+00, 0.27586E+00, 0.21297E+00, &
         &    0.14064E+00, 0.45588E-01, 0.65665E-02, 0.34232E-02, 0.53199E-03, 0.19811E+00, 0.16833E+00, &
         &    0.13536E+00, 0.11549E+00, 0.10649E+00, 0.93264E-01, 0.75720E-01, 0.56405E-01, 0.41865E-01, &
         &    0.59331E-02, 0.26510E-02, 0.40040E-03, 0.32328E+00, 0.26636E+00, 0.21397E+00, 0.14038E+00, &
         &    0.52142E-01, 0.38852E-02, 0.14601E+00, 0.13824E+00, 0.27703E+00, 0.22388E+00, 0.15446E+00, &
         &    0.48687E-01, 0.98054E-02, 0.18870E-02, 0.11961E+00, 0.12106E+00, 0.13215E+00, 0.13516E+00, &
         &    0.25249E+00, 0.16542E+00, 0.68157E-01, 0.59725E-02, 0.49258E+00, 0.33651E+00, 0.16182E+00, &
         &    0.90984E-02, 0.95202E+00, 0.47978E-01, 0.91716E+00, 0.82857E-01, 0.77464E+00, 0.22536E+00 /), (/ 1,140 /) )

    call planck_function_surf(1, 1, config, spread(temperature,1,1), &
         &                    frac, planck_surf)

  end subroutine planck_function

end module radiation_ifs_rrtm
