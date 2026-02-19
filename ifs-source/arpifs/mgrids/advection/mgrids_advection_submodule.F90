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

#ifdef WITH_MGRIDS
submodule (mgrids_advection_module) mgrids_advection_submodule
#ifdef WITH_MGRIDS
implicit none
contains

#define __FILENAME__ "mgrids_advection_submodule"

!------------------------------------------------------------------------------

module function mgrids_advection_create( &
  ydrip,                                 &
  yddyn,                                 &
  ydgeometry ) result(this)

  use fckit_module, only : fckit_exception, fckit_log
  use yom_atlas_ifs, only : latlas_mesh
  use yom_ygfl, only : ygfl
!Note: shouldn't be needed but required for Intel compiler
  use yomrip, only : trip
  use yomdyn, only : tdyn
  use type_geometry, only : geometry

  type(mgrids_advection) :: this
  type(trip),     intent(in), target :: ydrip
  type(tdyn),     intent(in), target :: yddyn
  type(geometry), intent(in), target :: ydgeometry

  type(mgrids_advection_config) :: config
  character(len=1024) :: string

  if( latlas_mesh ) then
    call config%read_namelist()
    config%ctype = to_lower(config%ctype)
  else
    config%ctype = "none"
  endif

  if( config%ctype /= "none" .and. ygfl%nfmg == 0 ) then
    call fckit_log%info('mgrids_advection type set from "'//trim(config%ctype)//&
      & '" to "none" as there were no tracers declared up for mgrids_advection.')
    config%ctype = "none"
  endif

  if( config%ctype == "none" )   call setup_no_advection()
  if( config%ctype == "sladv" )  call setup_mgrids_sladv()
  if( config%ctype == "mpdata" ) call setup_mgrids_mpdata()

  if( config%ctype /= "none" ) then
    if( .not. associated(this%implementation) ) then
      write(string,'(A,A,A,A)') 'ERROR: nam_mgrids_advection/ctype "',trim(config%ctype),'" unrecognised. Valid types are: ', &
        & '"none" (no advection), "sladv" (semi-Lagrangian advection), "mpdata" (MPDATA advection)'
      call fckit_exception%abort(string, __FILENAME__, __LINE__ )
    endif
  endif
contains

  subroutine setup_no_advection()
    call this%final()
  end subroutine

  subroutine setup_mgrids_sladv()
    use mgrids_sladv_module , only : create_mgrids_sladv
    call assert_latlasmesh()
    call create_mgrids_sladv( this%implementation,         &
              config,                                      &
              ydrip,                                       &
              yddyn,                                       &
              ydgeometry )
  end subroutine

  subroutine setup_mgrids_mpdata()
    use mgrids_mpdata_module  , only : create_mgrids_mpdata
    call assert_latlasmesh()
    call create_mgrids_mpdata( this%implementation,        &
              config,                                      &
              ydrip,                                       &
              yddyn,                                       &
              ydgeometry )
  end subroutine

  subroutine assert_latlasmesh()
    if( .not. latlas_mesh ) then
      call fckit_exception%abort('ERROR: latlas_mesh was set to False. mgrids advection requires latlas_mesh True', &
        & __FILENAME__,__LINE__)
    endif
  end subroutine

  function to_lower(string)
    character (len=:), allocatable :: to_lower
    character (len=*) , intent(in) :: string
    integer :: i,ic,nlen
    nlen = len(string)
    to_lower = string
    do i=1,nlen
      ic = ichar(to_lower(i:i))
      if (ic >= 65 .and. ic < 90) to_lower(i:i) = char(ic+32)
    end do
  end function

end function

!------------------------------------------------------------------------------

module function mgrids_advection_args_dynamic_cast(args) result(this)
  use fckit_module, only : fckit_exception
!Note: shouldn't be needed but required for Intel compiler
  use dwarf_module, only : dwarf_args

  type(mgrids_advection_args) :: this
  class(dwarf_args), intent(inout), target :: args
  select type( ptr => args )
    class is( mgrids_advection_args )
      this = ptr
    class default
      call fckit_exception%abort( "Cannot cast dwarf_args to mgrids_advection_args", &
        & __FILENAME__, __LINE__ )
  end select
end function

!------------------------------------------------------------------------------

module function active( this ) result( is_active )
  logical :: is_active
  class(mgrids_advection), intent(in) :: this
  if( associated( this%implementation ) ) then
    is_active = .true.
  else
    is_active = .false.
  endif
end function

!------------------------------------------------------------------------------

module subroutine execute(  &
    this,                   &
    ydgmv,                  &
    pgmv,                   &
    pgfl,                   &
    pwrl9,                  &
    pgeo0,                  &
    prcp0,                  &
    pre0f,                  &
    pcty0,                  &
    pgflt1                  &
  )
  !-------------------------------------------------------------------
  use parkind1  ,only : jprb
  use yomhook   ,only : lhook, dr_hook, jphook
!Note: shouldn't be needed but required for Intel compiler
  use yomgmv    ,only : tgmv
  !-------------------------------------------------------------------
  class(mgrids_advection), intent(inout)         :: this
  type(tgmv),              intent(in),    target :: ydgmv
  real(jprb),              intent(in),    target :: pgmv(:,:,:,:)
  real(jprb),              intent(in),    target :: pgfl(:,:,:,:)
  real(jprb),              intent(in),    target :: pwrl9(:,:,:)
  real(jprb),              intent(in),    target :: pgeo0(:,:,:)
  real(jprb),              intent(in),    target :: prcp0(:,:,:,:)
  real(jprb),              intent(in),    target :: pre0f(:,:,:)
  real(jprb),              intent(in),    target :: pcty0(:,0:,:,:)
  real(jprb),              intent(inout), target :: pgflt1(:,:,:,:)
  !-------------------------------------------------------------------
  type(mgrids_advection_args) :: args
  real(kind=jphook) :: zhook_handle
  !-------------------------------------------------------------------

  if( lhook) call dr_hook('mgrids_advection % execute',0,zhook_handle)

  if( associated( this%implementation ) ) then
    args%yrgmv  => ydgmv
    args%zgmv   => pgmv
    args%zgfl   => pgfl
    args%zwrl9  => pwrl9
    args%zgeo0  => pgeo0
    args%zrcp0  => prcp0
    args%zpre0f => pre0f
    args%zcty0  => pcty0
    args%zgflt1 => pgflt1
    call this%implementation%execute( args )
  endif

  if( lhook) call dr_hook('mgrids_advection % execute',1,zhook_handle)

end subroutine

!------------------------------------------------------------------------------

module subroutine final(this)
  class(mgrids_advection) :: this
  if( associated( this%implementation ) ) then
    call this%implementation%final()
  endif
  this%implementation => null()
end subroutine

!------------------------------------------------------------------------------

module subroutine dont_optimize_out( this )
  class(mgrids_advection_args), intent(inout)         :: this
end subroutine

!------------------------------------------------------------------------------

module subroutine mgrids_advection_args__final( this )
  !! Avoids that loops editing arguments will be compiled out (observed behaviour on Cray cce/8.6.2)
  type(mgrids_advection_args) :: this
end subroutine

!------------------------------------------------------------------------------

module subroutine read_namelist( this )
  use yomlun, only : nulnam
!Note: shouldn't be needed but required for Intel compiler
  use parkind1, only : jpim, jprb

  class(mgrids_advection_config) :: this

  character(len=32)    :: ctype = "none"
  character(len=32)    :: cgrid = "same"
  logical              :: lremesh = .false.

  character(len=512)   :: sladv_departurepoint_method = "SETTLS"
  character(len=512)   :: sladv_interpolation_method  = "quasicubic"
  logical              :: sladv_interpolation_limiter = .true.
  integer(jpim)        :: sladv_nhalo = 14

  logical              :: mpdata_lusecfl = .true.
  real(jprb)           :: mpdata_zcflmax = 0.8
  integer(jpim)        :: mpdata_nsubsteps = 8

  namelist/NAM_MGRIDS_ADVECTION/ctype,cgrid,lremesh,&
    & sladv_departurepoint_method,sladv_interpolation_method, sladv_interpolation_limiter, sladv_nhalo, &
    & mpdata_lusecfl,mpdata_zcflmax,mpdata_nsubsteps

#include "posnam.intfb.h"

  if( .not. this%initialised ) then
    call posnam(nulnam,'NAM_MGRIDS_ADVECTION')
    read(nulnam,NAM_MGRIDS_ADVECTION)
    this%ctype = ctype
    this%cgrid = cgrid
    this%lremesh = lremesh
    this%sladv_departurepoint_method = sladv_departurepoint_method
    this%sladv_interpolation_method = sladv_interpolation_method
    this%sladv_interpolation_limiter = sladv_interpolation_limiter
    this%sladv_nhalo = sladv_nhalo
    this%mpdata_lusecfl = mpdata_lusecfl
    this%mpdata_zcflmax = mpdata_zcflmax
    this%mpdata_nsubsteps = mpdata_nsubsteps

    if( trim(this%cgrid) /= "same" ) this%lremesh = .true.

  endif

  this%initialised = .true.

end subroutine

!------------------------------------------------------------------------------

#endif
end submodule
#else
module mgrids_advection_submodule
end module mgrids_advection_submodule
#endif
