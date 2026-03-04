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

module supergom_class
!-----------------------------------------------------------------------
!   Supergom - Manages the interpolation from model to observation space. Uses:
!
!      gom_mod  - low-level management of model input, interpolation to obs location, message passing to obs PE
!      yomglobs, mkglobstab - basic GOM tables for message passing and obs location
!      cobs, cobsall, slint, obshor - higher-level management of horizontal interpolation and model input
!      gom_plus - reconstructs additional parts of the model state (e.g. pressure profile)
!                 and provides a container for accessing model data in the obs operator
!
!   Author.
!   -------
!   Alan Geer 11-Feb-2016 - model <-> observation space interpolation in one object

!   Modifications.
!   --------------
!   F. Duruisseau 10-Jan-2019 : Allow to use BAYRAD (multi-point with SATEM), not ECMWF
!   Y.Hirahara 01-Apr-2019  Add Lake variable
!   A.Geer     21-May-2019  NetCDF gom_plus_dumps; timestep GOM variable
!   N.Bormann  05-Aug-2020  Support for different 2dGOM treatment per obs-set
!
!
!------------------------------------------------------------------------
use parkind1, only : jpim, jprb
use geometry_mod , only : geometry
use surface_fields_mix , only : tsurf
use yomgmv, only : tgmv
use yomectab, only : tecvar


implicit none
type, public :: class_supergom
  contains
  ! When using forecast-only or OpenIFS, the only routine used is model_in
  ! create a dummy
  procedure :: model_in
end type
contains
subroutine model_in(this,&
 & ydgeometry,ydgmv,ydsurf,ydphy,ydml_gconf,ydslint,yddyna,kdimgmv,ydecvar,ygfl,ydphysmwave,pgfl,pgmv,pgmvs,&
 & psp_sb,psp_sg,psp_sl,psp_rr,psp_cl,psp_x2,psp_ci,&
 & psd_vf,psd_vv,psd_vd,psd_ws,psd_vx,psd_vn,ktimestep,ldgp_done,ldupd_done)
  use yomdyn   , only : tdyn
  use yomphy   , only : tphy
  use model_general_conf_mod , only : model_general_conf_type
  use yomslint, only : tslint
  use yomdyna , only  : tdyna
  use yom_ygfl , only : type_gfld
  use yoe_phys_mwave, only : tephysmwave
  class(class_supergom),intent(inout) :: this
  type(geometry)       ,intent(in)    :: ydgeometry
  type(tgmv)           ,intent(inout) :: ydgmv
  type(tsurf)          ,intent(inout) :: ydsurf
  type(tphy)           ,intent(inout) :: ydphy
  type(model_general_conf_type),intent(in):: ydml_gconf
  type(tslint)         ,intent(in)    :: ydslint
  type(tdyna)          , intent(in)     :: yddyna
  integer(kind=jpim)   ,intent(in)    :: kdimgmv
  type(tecvar)         ,intent(in)    :: ydecvar
  type(type_gfld)      ,intent(in)    :: ygfl
  type(TEPHYSMWAVE)    ,intent(in)    :: ydphysmwave
  real(kind=jprb)      ,intent(in)    :: pgfl(ydgeometry%yrdim%nproma,ydgeometry%yrdimv%nflevg,ygfl%ndim,ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: pgmv(ydgeometry%yrdim%nproma,ydgeometry%yrdimv%nflevg,kdimgmv,ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: pgmvs(ydgeometry%yrdim%nproma,ydgmv%ndimgmvs,ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: psp_sb(ydgeometry%yrdim%nproma,ydsurf%ysp_sbd%nlevs,ydsurf%ysp_sbd%ndim,&
                                       & ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: psp_sg(ydgeometry%yrdim%nproma,ydsurf%ysp_sgd%nlevs,ydsurf%ysp_sgd%ndim,ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: psp_sl(ydgeometry%yrdim%nproma,ydsurf%ysp_sld%ndim,ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: psp_rr(ydgeometry%yrdim%nproma,ydsurf%ysp_rrd%ndim,ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: psp_cl(ydgeometry%yrdim%nproma,ydsurf%ysp_cld%ndim,ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: psp_x2(ydgeometry%yrdim%nproma,ydsurf%ysp_x2d%ndim,ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: psp_ci(ydgeometry%yrdim%nproma,ydsurf%ysp_cid%ndim,ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: psd_vf(ydgeometry%yrdim%nproma,ydsurf%ysd_vfd%ndim,ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: psd_vv(ydgeometry%yrdim%nproma,ydsurf%ysd_vvd%ndim,ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: psd_vd(ydgeometry%yrdim%nproma,ydsurf%ysd_vdd%ndim,ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: psd_ws(ydgeometry%yrdim%nproma,ydsurf%ysd_wsd%ndim,ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: psd_vx(ydgeometry%yrdim%nproma,ydsurf%ysd_vxd%ndim,ydgeometry%yrdim%ngpblks)
  real(kind=jprb)      ,intent(in)    :: psd_vn(ydgeometry%yrdim%nproma,ydsurf%ysd_vnd%ndim,ydgeometry%yrdim%ngpblks)
  integer(kind=jpim)   ,intent(in)    :: ktimestep ! Temporarily the best available time information
  logical              ,intent(in)    :: ldgp_done ! to handle the last timestep, when gp model does not run
  logical              ,intent(in)    :: ldupd_done ! used by neutral winds (affects scatterometer only)
  call abort1("oifs/fc-only - supergom_class%model, YDGOM%MODEL_IN or YDGOM5%MODEL_IN should never be called")

end subroutine model_in
end module supergom_class
