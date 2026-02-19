! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#include "dwarf_sladv.h"

module dwarf_sladv_departurepoint_method_module

  use dwarf_sladv_use_module

  implicit none

  private

  type, abstract, private :: Implementation
    class(geometry_type), pointer :: geometry
        !! Note gfortran 6.3 requires "class" instead of "type" to avoid internal compiler error
      
    type(atlas_Field) :: wind_ext
    type(atlas_Nabla) :: nabla
    
  contains

    procedure(intf__setup),       public, deferred :: setup
    procedure(intf__execute),     public, deferred :: execute
    procedure(intf__needs_nabla), public, deferred :: needs_nabla
    
    procedure, private :: extrapolate_wind_constant
    procedure, private :: extrapolate_wind_half
    procedure, private :: extrapolate_wind_next
    
    procedure, private :: create_wind_ext

  end type

  type, public :: DeparturePointMethod_type

    class(Implementation), pointer :: impl => null()

  contains

    procedure, public, nopass :: needs_nabla
    procedure, public :: setup
    procedure, public :: execute
    
    procedure, public :: final
    final :: auto_final

  end type

interface

  subroutine intf__setup( this, config, geometry, nabla )
    import Implementation
    import geometry_type
    import atlas_Config
    import atlas_Nabla
    class(Implementation),            intent(inout)    :: this
    type(geometry_type), target,      intent(in)       :: geometry
    type(atlas_Config),               intent(in)       :: config
    type(atlas_Nabla),                intent(in)       :: nabla
  end subroutine

  subroutine intf__execute( this, istep, dt, wind_tm1, wind_t0, departure_points )
    import Implementation
    import ip, wp
    import atlas_Field
    import atlas_FieldSet
    class(Implementation),            intent(inout)    :: this
    integer(ip),                      intent(in)       :: istep
    real(wp),                         intent(in)       :: dt
    type(atlas_Field),                intent(in)       :: wind_tm1
    type(atlas_Field),                intent(in)       :: wind_t0
    type(atlas_FieldSet),             intent(inout)    :: departure_points
  end subroutine

  function intf__needs_nabla( this ) result(needs_nabla)
    import Implementation
    import atlas_Nabla
    logical :: needs_nabla
    class(Implementation), intent(in)    :: this
  end function

end interface


  type, public, extends(Implementation) :: SETTLS
    integer(ip) :: dp_niter
  contains
    procedure, public :: setup => SETTLS_setup
    procedure, public :: execute => SETTLS_execute
    procedure, public :: needs_nabla => SETTLS_needs_nabla
  end type

  type, public, extends(SETTLS) :: SETTLSVF
    ! Like SETTLS, but introduces a limiter (Vertical Filter) 
    ! in stratosphere levels
    real(wp)    :: scale, scale_off
    integer(ip) :: nlev_vf
  contains
    procedure, public :: setup => SETTLSVF_setup
    procedure, public :: execute => SETTLSVF_execute
    procedure, private :: extrapolate_wind_next_VF
  end type

  type, public, extends(Implementation) :: NonIterative
  contains
    procedure, public :: setup => NonIterative_setup
    procedure, public :: execute => NonIterative_execute
    procedure, public :: needs_nabla => NonIterative_needs_nabla
  end type

contains

!--------------------------------------------------------------------------------------

subroutine implementation_factory( config, impl )
  type(atlas_Config), intent(in) :: config
  class(Implementation), pointer :: impl
  character(len=:), allocatable :: dp_meth
  if( .not. config%get("departurepoint_method",dp_meth) ) then
    call ABORT("departurepoint_method not specified")
  endif
  if( dp_meth == "NonIterative" ) then
    allocate( NonIterative :: impl )
  elseif( dp_meth == "SETTLS" ) then
    allocate( SETTLS :: impl )
  elseif( dp_meth == "SETTLSVF" ) then
    allocate( SETTLSVF :: impl )
  else
    call ABORT("departurepoint_method "//dp_meth//" not valid")
  endif
end subroutine

!--------------------------------------------------------------------------------------

function needs_nabla( config )
  logical :: needs_nabla
  type(atlas_Config), intent(in) :: config
  !---------------------------------------------
  class(Implementation), pointer :: impl
  !---------------------------------------------
  call implementation_factory(config,impl)
  needs_nabla = impl%needs_nabla()
  deallocate( impl )
end function

!--------------------------------------------------------------------------------------

subroutine setup(this, config, geometry, nabla )
  class(DeparturePointMethod_type), intent(inout) :: this
  type(atlas_Config) :: config
  type(geometry_type) :: geometry
  type(atlas_Nabla), optional :: nabla
  !---------------------------------------------
  call implementation_factory(config,this%impl)
  if( this%impl%needs_nabla() ) then
    if( .not. present(nabla) ) then
      call ABORT("nabla not present")
    endif
    this%impl%nabla = nabla
  endif
  call this%impl%setup(config,geometry,nabla)
end subroutine

!--------------------------------------------------------------------------------------

subroutine execute(this, istep, dt, wind_tm1, wind_t0, departure_points )
  class(DeparturePointMethod_type), intent(inout)    :: this
  integer(ip),                      intent(in)       :: istep
  real(wp),                         intent(in)       :: dt
  type(atlas_Field),                intent(in)       :: wind_tm1
  type(atlas_Field),                intent(in)       :: wind_t0
  type(atlas_FieldSet),             intent(inout)    :: departure_points
  call this%impl%execute(istep,dt,wind_tm1,wind_t0,departure_points)
end subroutine

!--------------------------------------------------------------------------------------

subroutine extrapolate_wind_constant(this, wind_t0, wind_ext )
  class(Implementation), intent(inout) :: this
  type(atlas_Field),     intent(in)    :: wind_t0
  type(atlas_Field),     intent(inout) :: wind_ext
  !------------------------------------------------------------
  integer(ip) :: jnode
  real(wp), pointer :: uvw_t0(:,:,:)
  real(wp), pointer :: uvw_ext(:,:,:)
  !------------------------------------------------------------

  call wind_t0%halo_exchange()

  call wind_t0%data(uvw_t0)
  call wind_ext%data(uvw_ext)

  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
  do jnode=1,this%geometry%ngp
    uvw_ext(:,:,jnode)=uvw_t0(:,:,jnode)
  enddo
  !$OMP END PARALLEL DO
end subroutine

!--------------------------------------------------------------------------------------

subroutine extrapolate_wind_half(this, wind_tm1, wind_t0, wind_ext )
  !!> U(t+0.5*dt) = 1.5*U(t) - 0.5*U(t-dt)

  class(Implementation), intent(in) :: this
  type(atlas_Field),                intent(in)    :: wind_tm1
  type(atlas_Field),                intent(in)    :: wind_t0
  type(atlas_Field),                intent(inout) :: wind_ext
  !------------------------------------------------------------
  integer(ip) :: jnode
  real(wp), pointer :: uvw_tm1(:,:,:)
  real(wp), pointer :: uvw_t0(:,:,:)
  real(wp), pointer :: uvw_ext(:,:,:)
  !------------------------------------------------------------

  call wind_tm1%data(uvw_tm1)
  call wind_t0%data(uvw_t0)
  call wind_ext%data(uvw_ext)

  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
  do jnode=1,this%geometry%ngp
    uvw_ext(:,:,jnode)=1.5_wp*uvw_t0(:,:,jnode)-0.5_wp*uvw_tm1(:,:,jnode)
  enddo
  !$OMP END PARALLEL DO
  call wind_ext%set_dirty()
end subroutine

!--------------------------------------------------------------------------------------

subroutine extrapolate_wind_next(this, wind_tm1, wind_t0, wind_ext )
  !!> U(t+dt)     = 2.0*U(t) - 1.0*U(t-dt)
  !!> As used by SETTLS

  class(Implementation), intent(in) :: this
  type(atlas_Field),                intent(in)    :: wind_tm1
  type(atlas_Field),                intent(in)    :: wind_t0
  type(atlas_Field),                intent(inout) :: wind_ext
  !------------------------------------------------------------
  integer(ip) :: jnode, jlev, jvar
  real(wp), pointer :: uvw_tm1(:,:,:)
  real(wp), pointer :: uvw_t0(:,:,:)
  real(wp), pointer :: uvw_ext(:,:,:)
  !------------------------------------------------------------

  call wind_tm1%data(uvw_tm1)
  call wind_t0%data(uvw_t0)
  call wind_ext%data(uvw_ext)

  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
  do jnode=1,this%geometry%ngp
    uvw_ext(:,:,jnode)=2.0_wp*uvw_t0(:,:,jnode)-uvw_tm1(:,:,jnode)
  enddo
  !$OMP END PARALLEL DO
  call wind_ext%set_dirty()
end subroutine

!--------------------------------------------------------------------------------------

subroutine extrapolate_wind_next_VF(this, wind_tm1, wind_t0, wind_ext, nlev_vf, scale, scale_off )
  !!> U(t+dt)     = 2.0*U(t) - 1.0*U(t-dt)
  !!> As used by SETTLS

  ! Difference with no VF: detect oscillating gridpoints which should not be extrapolated

  class(SETTLSVF), intent(in) :: this
  type(atlas_Field),                intent(in)    :: wind_tm1
  type(atlas_Field),                intent(in)    :: wind_t0
  type(atlas_Field),                intent(inout) :: wind_ext
  integer(ip),                      intent(in)    :: nlev_vf
  real(wp),                         intent(in)    :: scale
  real(wp),                         intent(in)    :: scale_off
  !------------------------------------------------------------
  integer(ip) :: jnode, jlev, jvar
  real(wp), pointer :: uvw_tm1(:,:,:)
  real(wp), pointer :: uvw_t0(:,:,:)
  real(wp), pointer :: uvw_ext(:,:,:)
  real(wp) :: alpha
  !------------------------------------------------------------

  call wind_tm1%data(uvw_tm1)
  call wind_t0%data(uvw_t0)
  call wind_ext%data(uvw_ext)

  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev,alpha)
  do jnode=1,this%geometry%ngp
    do jlev=1,nlev_vf
      uvw_ext(1,jlev,jnode)=2.0_wp*uvw_t0(1,jlev,jnode)-uvw_tm1(1,jlev,jnode)
      uvw_ext(2,jlev,jnode)=2.0_wp*uvw_t0(2,jlev,jnode)-uvw_tm1(2,jlev,jnode)

      alpha = 0.5_wp*(1._wp-tanh(-scale*(uvw_t0(3,jlev,jnode)*uvw_tm1(3,jlev,jnode)+scale_off)))
      uvw_ext(3,jlev,jnode) = uvw_t0(3,jlev,jnode) + alpha*(uvw_t0(3,jlev,jnode)-uvw_tm1(3,jlev,jnode))
    enddo
    do jlev=nlev_vf+1,this%geometry%nlev
      uvw_ext(:,jlev,jnode) = 2.0_wp*uvw_t0(:,jlev,jnode)-uvw_tm1(:,jlev,jnode)
    enddo
  enddo
  !$OMP END PARALLEL DO
  call wind_ext%set_dirty()
end subroutine

!--------------------------------------------------------------------------------------

subroutine SETTLS_setup( this, config, geometry, nabla )
  class(SETTLS),                    intent(inout)    :: this
  type(atlas_Config),               intent(in)       :: config
  type(geometry_type), target,      intent(in)       :: geometry
  type(atlas_Nabla),                intent(in)       :: nabla
  if( .not. config%get("departurepoint_iterations", this%dp_niter ) ) then
    call ABORT( "departurepoint_iterations not found in config" )
  endif
  this%geometry => geometry
  call this%create_wind_ext()
end subroutine

!--------------------------------------------------------------------------------------

subroutine SETTLSVF_setup( this, config, geometry, nabla )
  class(SETTLSVF),                  intent(inout)    :: this
  type(atlas_Config),               intent(in)       :: config
  type(geometry_type), target,      intent(in)       :: geometry
  type(atlas_Nabla),                intent(in)       :: nabla
  call SETTLS_setup(this,config,geometry,nabla)
  if( .not. config%get("SETTLSVF_levels", this%nlev_vf ) ) then
    call ABORT( "SETTLSVF_levels not found in config" )
  endif
  if( .not. config%get("SETTLSVF_scale", this%scale ) ) then
    call ABORT( "SETTLSVF_scale not found in config" )
  endif
end subroutine

!--------------------------------------------------------------------------------------

subroutine NonIterative_setup( this, config, geometry, nabla )
  class(NonIterative),              intent(inout)    :: this
  type(atlas_Config),               intent(in)       :: config
  type(geometry_type), target,      intent(in)       :: geometry
  type(atlas_Nabla),                intent(in)       :: nabla
  this%geometry => geometry
  this%nabla = nabla
  call this%create_wind_ext()
end subroutine

!--------------------------------------------------------------------------------------

subroutine SETTLS_execute( this, istep, dt, wind_tm1, wind_t0, departure_points )
  use dwarf_sladv_compute_dp_iter_module
  class(SETTLS),                    intent(inout)    :: this
  integer(ip),                      intent(in)       :: istep
  real(wp),                         intent(in)       :: dt
  type(atlas_Field),                intent(in)       :: wind_tm1
  type(atlas_Field),                intent(in)       :: wind_t0
  type(atlas_FieldSet),             intent(inout)    :: departure_points
  if( istep == 0 ) then
    call this%extrapolate_wind_constant(wind_t0, this%wind_ext)
  else
    call this%extrapolate_wind_next(wind_tm1,wind_t0,this%wind_ext)
  endif
  call compute_dp_iter(dt,this%dp_niter,this%geometry,wind_t0,this%wind_ext,departure_points)
end subroutine

!--------------------------------------------------------------------------------------

subroutine SETTLSVF_execute( this, istep, dt, wind_tm1, wind_t0, departure_points )
  use dwarf_sladv_compute_dp_iter_module
  class(SETTLSVF),                  intent(inout)    :: this
  integer(ip),                      intent(in)       :: istep
  real(wp),                         intent(in)       :: dt
  type(atlas_Field),                intent(in)       :: wind_tm1
  type(atlas_Field),                intent(in)       :: wind_t0
  type(atlas_FieldSet),             intent(inout)    :: departure_points
  if( istep == 0 ) then
    call this%extrapolate_wind_constant(wind_t0, this%wind_ext)
  else
    call this%extrapolate_wind_next_VF(wind_tm1,wind_t0,this%wind_ext,this%nlev_vf,this%scale,this%scale_off)
  endif
  call compute_dp_iter(dt,this%dp_niter,this%geometry,wind_t0,this%wind_ext,departure_points)
end subroutine

!--------------------------------------------------------------------------------------

subroutine NonIterative_execute( this, istep,dt, wind_tm1, wind_t0, departure_points )
  use dwarf_sladv_compute_dp_noniter_module
  class(NonIterative),              intent(inout)    :: this
  integer(ip),                      intent(in)       :: istep
  real(wp),                         intent(in)       :: dt
  type(atlas_Field),                intent(in)       :: wind_tm1
  type(atlas_Field),                intent(in)       :: wind_t0
  type(atlas_FieldSet),             intent(inout)    :: departure_points
  if( istep == 0 ) then
    call this%extrapolate_wind_constant(wind_t0, this%wind_ext)
  else
    call this%extrapolate_wind_half(wind_tm1,wind_t0,this%wind_ext)
  endif
  call compute_dp_noniter(dt,this%geometry,wind_t0,this%wind_ext,this%nabla,departure_points)
end subroutine

!--------------------------------------------------------------------------------------

function SETTLS_needs_nabla( this ) result(needs_nabla)
  logical :: needs_nabla
  class(SETTLS), intent(in)    :: this
  needs_nabla = .false.
end function

!--------------------------------------------------------------------------------------

function NonIterative_needs_nabla( this ) result(needs_nabla)
  logical :: needs_nabla
  class(NonIterative), intent(in)    :: this
  needs_nabla = .true.
end function

!--------------------------------------------------------------------------------------

subroutine create_wind_ext( this )
    use atlas_module
    class( Implementation ) :: this
    this%wind_ext = this%geometry%fs_structuredcolumns%create_field( &
      name="wind(t+x)", kind=atlas_real(wp), variables=3, type="vector" )
end subroutine

subroutine final( this )
    class(DeparturePointMethod_type) :: this
    if( associated(this%impl) ) deallocate(this%impl)
end subroutine

subroutine auto_final( this )
    type(DeparturePointMethod_type) :: this
    call this%final()
end subroutine

end module
