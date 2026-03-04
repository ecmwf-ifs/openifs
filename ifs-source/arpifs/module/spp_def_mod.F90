! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

module spp_def_mod
  use parkind1, only: jpim, jprb
  use spp_gen_mod

  type spp_pert_pointer
    !
    !   perturbation pointer type
    !
    !   A value of -1 implies that the perturbation is off, otherwise
    !     it points to the spp_pert perturbation in the spp_model
    !
    integer(kind=jpim) :: cfm        = -1
    integer(kind=jpim) :: cfm1       = -1
    integer(kind=jpim) :: cfm2       = -1
    integer(kind=jpim) :: cfm3       = -1
    integer(kind=jpim) :: rkap       = -1
    integer(kind=jpim) :: rkap1      = -1
    integer(kind=jpim) :: rkap2      = -1
    integer(kind=jpim) :: rkap3      = -1
    integer(kind=jpim) :: tofdc      = -1
    integer(kind=jpim) :: hsdt       = -1
    integer(kind=jpim) :: vdexc      = -1
    integer(kind=jpim) :: entrorg    = -1
    integer(kind=jpim) :: entshalp   = -1
    integer(kind=jpim) :: entstpc1   = -1
    integer(kind=jpim) :: detrpen    = -1
    integer(kind=jpim) :: rprcon     = -1
    integer(kind=jpim) :: rtau       = -1
    integer(kind=jpim) :: cududv     = -1
    integer(kind=jpim) :: cududvs    = -1
    integer(kind=jpim) :: ramid      = -1
    integer(kind=jpim) :: rcldiff    = -1
    integer(kind=jpim) :: rclcrit    = -1
    integer(kind=jpim) :: rlcritsnow = -1
    integer(kind=jpim) :: rainevap   = -1
    integer(kind=jpim) :: snowsublim = -1
    integer(kind=jpim) :: cloudinhom = -1
    integer(kind=jpim) :: qsatvervel = -1
    integer(kind=jpim) :: zdecorr    = -1
    integer(kind=jpim) :: zsigqcw    = -1
    integer(kind=jpim) :: zradeffl   = -1
    integer(kind=jpim) :: zradeffi   = -1
    integer(kind=jpim) :: phr        = -1
    integer(kind=jpim) :: zhs_vdaero = -1
    integer(kind=jpim) :: delta_aero = -1
  contains
    procedure :: set => set_spp_pert_pointer

  end type spp_pert_pointer


contains

  subroutine set_spp_pert_pointer( this, cd_pert, kptr )
    !
    !   set this%<cd_pert> = kptr
    !
#include "abor1.intfb.h"
    class(spp_pert_pointer)   :: this
    character(len=jp_lab_len), intent(in) :: cd_pert   ! name of perturbation
    integer(kind=jpim),        intent(in) :: kptr      ! pointer to perturbation

    select case(cd_pert)
      case ('CFM')
        this%cfm         = kptr
      case ('CFM1')
        this%cfm1        = kptr
      case ('CFM2')
        this%cfm2        = kptr
      case ('CFM3')
        this%cfm3        = kptr
      case ('RKAP')
        this%rkap        = kptr
      case ('RKAP1')
        this%rkap1       = kptr
      case ('RKAP2')
        this%rkap2       = kptr
      case ('RKAP3')
        this%rkap3       = kptr
      case ('TOFDC')
        this%tofdc       = kptr
      case ('HSDT')
        this%hsdt        = kptr
      case ('VDEXC')
        this%vdexc       = kptr
      case ('ENTRORG')
        this%entrorg     = kptr
      case ('ENTSHALP')
        this%entshalp    = kptr
      case ('ENTSTPC1')
        this%entstpc1    = kptr
      case ('DETRPEN')
        this%detrpen     = kptr
      case ('RPRCON')
        this%rprcon      = kptr
      case ('RTAU')
        this%rtau        = kptr
      case ('CUDUDV')
        this%cududv      = kptr
      case ('CUDUDVS')
        this%cududvs     = kptr
      case ('RAMID')
        this%ramid       = kptr
      case ('RCLDIFF')
        this%rcldiff     = kptr
      case ('RCLCRIT')
        this%rclcrit     = kptr
      case ('RLCRITSNOW')
        this%rlcritsnow  = kptr
      case ('CLOUDINHOM')
        this%cloudinhom  = kptr
      case ('RAINEVAP')
        this%rainevap    = kptr
      case ('SNOWSUBLIM')
        this%snowsublim  = kptr
      case ('QSATVERVEL')
        this%qsatvervel  = kptr
      case ('ZDECORR')
        this%zdecorr     = kptr
      case ('ZSIGQCW')
        this%zsigqcw     = kptr
      case ('ZRADEFFL')
        this%zradeffl    = kptr
      case ('ZRADEFFI')
        this%zradeffi    = kptr
      case ('PHR')
        this%phr         = kptr
      case ('ZHS_VDAERO')
        this%zhs_vdaero  = kptr
      case ('DELTA_AERO')
        this%delta_aero  = kptr
      case default
        call abor1("spp_pert_pointer%set: the perturbation label is unknown.")
    end select

  end subroutine set_spp_pert_pointer


subroutine define_spp_arpege( sm, nulout )

#include "abor1.intfb.h"
  type(spp_model), intent(inout) :: sm
  integer(kind=jpim),intent(in)  :: nulout
  integer(kind=jpim) :: indef

  sm%version = "SPP for ARPEGE physics, v0.1"
  write(nulout,*) "****************************************"
  write(nulout,*) "define "//sm%version
  write(nulout,*) "****************************************"

  !
  !   set global default values for SPP configuration
  !
  sm%tau=259.2E3_jprb ! 72h;  decorrelation time (s)
  sm%sdev=1.0_jprb
  sm%xlcor = 2000.E3_jprb  ! 2000km; (m)

  sm%kseed_off = 2**8
  ! . . .

  !
  !   allocate arrays
  !
  call allocate_spp_model( sm, 20 )
  sm%defined_perts = cp_undefined

  !
  !  implement perturbations
  !
  !
  call implement_spp_pertn( sm, 'PSIGQSAT'    , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'CLDDPTH'     , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'CLDDPTHDP'   , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'ICE_CLD_WGT' , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'ICENU'       , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'KGN_ACON'    , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'KGN_SBGR'    , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'RADGR'       , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'RADSN'       , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'RFAC_TWOC'   , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'RZC_H'       , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'RZL_INF'     , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'RSWINHF'     , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'RLWINHF'     , [0.1_jprb]  )

  indef = count(sm%defined_perts /= cp_undefined )
  if (indef /= sm%ndef ) then
    call abor1('Inconsistency in terms of the number of defined perturbations.')
  endif
  write(nulout,*) "There are ",sm%ndef," possible perturbations"
  !
  !   inherit global settings
  !
  sm%pndef(1:sm%ndef)%tau    = sm%tau
  sm%pndef(1:sm%ndef)%sdev   = sm%sdev
  sm%pndef(1:sm%ndef)%xlcor  = sm%xlcor


end subroutine define_spp_arpege

subroutine define_spp_ifs( sm, nulout  )
  !
  !   For a science reference see
  !       Ollinaho et al. (2017), https://doi.org/10.1002/qj.3926
  !       Lang et al. (2021) https://doi.org/10.1002/qj.3978

#include "abor1.intfb.h"

  type(spp_model), intent(inout) :: sm
  integer(kind=jpim),intent(in)  :: nulout
  integer(kind=jpim) :: indef

  sm%version = "SPP for IFS physics"
  write(nulout,*) "****************************************"
  write(nulout,*) "define "//sm%version
  write(nulout,*) "****************************************"

  !
  !   set global default values for SPP configuration
  !
  sm%tau=259.2E3_jprb ! 72h;  decorrelation time (s)
  sm%sdev=1.0_jprb
  sm%xlcor = 1000.E3_jprb  ! 2000km; (m)

  sm%kseed_off = 2**8
  !
  !   allocate arrays
  !
  call allocate_spp_model( sm, 36 )
  sm%defined_perts = cp_undefined
  !
  !  implement perturbations
  !
  !
  call implement_spp_pertn( sm, 'CFM'       , [0.78_jprb,0.33_jprb] , knseed_off =[0], ln1=.true.) ! amplitudes for [land, ocean]
  call implement_spp_pertn( sm, 'CFM1'      , [0.65_jprb,0.26_jprb] , cd_rf_pert_label='CFM' , ln1=.true.)
  call implement_spp_pertn( sm, 'CFM2'      , [0.78_jprb,0.33_jprb] , cd_rf_pert_label='CFM' , ln1=.false.)
  call implement_spp_pertn( sm, 'CFM3'      , [0.78_jprb,0.33_jprb] , cd_rf_pert_label='CFM' , ln1=.false.)
  call implement_spp_pertn( sm, 'RKAP'      , [0.26_jprb] , knseed_off =[2**30 + 2**15], ln1=.false.)
  call implement_spp_pertn( sm, 'RKAP1'     , [0.065_jprb], cd_rf_pert_label='RKAP', ln1=.true.)
  call implement_spp_pertn( sm, 'RKAP2'     , [0.26_jprb] , cd_rf_pert_label='RKAP', ln1=.false.)
  call implement_spp_pertn( sm, 'RKAP3'     , [0.26_jprb] , cd_rf_pert_label='RKAP', ln1=.false.)
  call implement_spp_pertn( sm, 'TOFDC'     , [0.78_jprb] , knseed_off = [2** 8], ln1=.true. )
  call implement_spp_pertn( sm, 'HSDT'      , [0.52_jprb] , knseed_off = [2** 9], ln1=.true. )
  call implement_spp_pertn( sm, 'VDEXC'     , [1.04_jprb] , knseed_off = [2**10], ln1=.true. )
  call implement_spp_pertn( sm, 'ENTRORG'   , [0.39_jprb] , knseed_off = [2**11] , ln1=.true.)
  call implement_spp_pertn( sm, 'ENTSHALP'  , [0.39_jprb] , knseed_off = [2**12] , ln1=.true.)
  call implement_spp_pertn( sm, 'ENTSTPC1'  , [0.39_jprb] , knseed_off = [2**30 + 2**12], ln1=.true.)
  call implement_spp_pertn( sm, 'DETRPEN'   , [0.39_jprb] , knseed_off = [2**13] , ln1=.true.)
  call implement_spp_pertn( sm, 'RPRCON'    , [0.52_jprb] , knseed_off = [2**14] , ln1=.true.)
  call implement_spp_pertn( sm, 'RTAU'      , [0.78_jprb] , knseed_off = [2**15] , ln1=.true.)
  call implement_spp_pertn( sm, 'CUDUDV'    , [0.122_jprb, 0.122_jprb] , knseed_off = [2**23] , kidistr=0)
  call implement_spp_pertn( sm, 'CUDUDVS'   , [0.133_jprb, 0.133_jprb] , knseed_off = [2**24] , kidistr=0)
  call implement_spp_pertn( sm, 'RAMID'     , [0.13_jprb] , knseed_off = [2**16] , ln1=.true.)
  call implement_spp_pertn( sm, 'RCLDIFF'   , [1.04_jprb] , knseed_off = [2**17] , ln1=.true.)
  call implement_spp_pertn( sm, 'RCLCRIT'   , [1.04_jprb, 1.04_jprb] , knseed_off = [2**18] , ln1=.true.) ! land, sea
  call implement_spp_pertn( sm, 'RLCRITSNOW', [0.78_jprb] , knseed_off = [2**19] , ln1=.true.)
  call implement_spp_pertn( sm, 'RAINEVAP'  , [0.65_jprb] , knseed_off = [2**30] , ln1=.true.)
  call implement_spp_pertn( sm, 'SNOWSUBLIM', [0.65_jprb] , knseed_off = [2**30 + 2**8] , ln1=.true.)
  call implement_spp_pertn( sm, 'CLOUDINHOM', [0.3_jprb,0.3_jprb] , knseed_off = [2**30 + 2**10, 2**30 + 2**11]&
       , ln1=.true., knrf=2) ! autoconversion, accretion
  call implement_spp_pertn( sm, 'QSATVERVEL', [0.39_jprb] , knseed_off = [2**30 + 2**9] , ln1=.true.)
  call implement_spp_pertn( sm, 'ZDECORR'   , [0.78_jprb] , knseed_off = [2**20] , ln1=.true., radgrid=.true.)
  call implement_spp_pertn( sm, 'ZSIGQCW'   , [0.52_jprb] , knseed_off = [2**21] , ln1=.true., radgrid=.true.)
  call implement_spp_pertn( sm, 'ZRADEFFL'  , [0.78_jprb] , knseed_off = [2**22] , ln1=.true., radgrid=.true.) ! liquid phase
  call implement_spp_pertn( sm, 'ZRADEFFI'  , [0.78_jprb] , ln1=.true., radgrid=.true., cd_rf_pert_label='ZRADEFFL') ! ice phase
  call implement_spp_pertn( sm, 'PHR'       , [0.044_jprb,0.064_jprb, 0.064_jprb, 0.034_jprb], &
       & knseed_off = [2**28,2**29] , knrf=2, ln1=.false.)  !four different values for LW_ST/ LW_TR/  SW_ST/ SW_TR
  call implement_spp_pertn( sm, 'ZHS_VDAERO', [1.04_jprb] , knseed_off = [2**25] , ln1=.true., radgrid=.true.)
  call implement_spp_pertn( sm, 'DELTA_AERO', [0.78_jprb] , knseed_off = [2**26] , ln1=.true., radgrid=.true.)

  indef = count(sm%defined_perts /= cp_undefined )
  if (indef /= sm%ndef ) then
    call abor1('Inconsistency in terms of the number of defined perturbations.')
  endif
  write(nulout,*) "There are ",sm%ndef," possible perturbations"
end subroutine define_spp_ifs

subroutine define_spp_lam( sm, nulout  )
  type(spp_model), intent(inout) :: sm
  integer(kind=jpim),intent(in)  :: nulout
  integer(kind=jpim) :: indef

  sm%version = "SPP for LAM physics, v0.1"
  write(nulout,*) "****************************************"
  write(nulout,*) "define "//sm%version
  write(nulout,*) "****************************************"

  !
  !   set global default values for SPP configuration
  !
  sm%tau=259.2E3_jprb ! 72h;  decorrelation time (s)
  sm%sdev=1.0_jprb
  sm%xlcor = 2000.E3_jprb  ! 2000km; (m)

  sm%kseed_off = 2**8
  ! . . .

  !
  !   allocate arrays
  !
  call allocate_spp_model( sm, 20 )
  sm%defined_perts = cp_undefined

  !
  !  implement perturbations
  !
  !
  call implement_spp_pertn( sm, 'PSIGQSAT'    , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'CLDDPTH'     , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'CLDDPTHDP'   , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'ICE_CLD_WGT' , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'ICENU'       , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'KGN_ACON'    , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'KGN_SBGR'    , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'RADGR'       , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'RADSN'       , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'RFAC_TWOC'   , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'RZC_H'       , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'RZL_INF'     , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'RSWINHF'     , [0.1_jprb]  )
  call implement_spp_pertn( sm, 'RLWINHF'     , [0.1_jprb]  )

  indef = count(sm%defined_perts /= cp_undefined )
  if (indef /= sm%ndef ) then
    call abor1('Inconsistency in terms of the number of defined perturbations.')
  endif
  write(nulout,*) "There are ",sm%ndef," possible perturbations"
  !
  !   inherit global settings
  !
  sm%pndef(1:sm%ndef)%tau    = sm%tau
  sm%pndef(1:sm%ndef)%sdev   = sm%sdev
  sm%pndef(1:sm%ndef)%xlcor  = sm%xlcor


end subroutine define_spp_lam


end module spp_def_mod
