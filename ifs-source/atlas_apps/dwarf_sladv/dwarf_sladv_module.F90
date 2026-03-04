! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#include "dwarf_sladv.h"

module dwarf_sladv_module

use dwarf_sladv_use_module

use dwarf_sladv_departurepoint_method_module

implicit none

public  :: dwarf_sladv
public  :: dwarf_args
private :: dwarf_sladv_setup
private :: execute
private :: final

!------------------------------------------------------------------------------

type, abstract, extends( dwarf ) :: dwarf_sladv
  ! Needs setting up in concrete class
  integer(ip) :: ntrac
  integer(ip) :: istep
  real(wp)    :: dt
  
  type(Geometry_type) :: geometry
  logical                       :: interpolation_limiter
  character(len=:), allocatable :: interpolation_method

  type(DeparturePointMethod_type) :: departurepoint_calculation

  ! For computation of gradients, optional for McGregor scheme to departure point calculation
  type(atlas_nabla) :: nabla

  ! No need to set up in concrete class
  type(atlas_FieldSet) :: departure_points
  type(atlas_Field) :: tracer_field, arrival_tracer_field
  type(atlas_Field) :: wind_t9, wind_t0
contains
  procedure, public :: execute
  procedure(execute), deferred, public :: preprocess
  procedure(execute), deferred, public :: postprocess
  procedure, public :: final
  procedure, public :: dwarf_sladv_setup

  procedure, public, nopass :: needs_nabla

end type

!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------

function needs_nabla( config )
  logical :: needs_nabla
  type(atlas_Config), intent(in) :: config
  !---------------------------------------------
  integer(ip) :: dp_meth
  type(DeparturePointMethod_type) :: departurepoint_method
  !---------------------------------------------
  needs_nabla = departurepoint_method%needs_nabla( config )
end function

!------------------------------------------------------------------------------

subroutine dwarf_sladv_setup( this, config, functionspace, nabla )
  class(dwarf_sladv), intent(inout) :: this
  type(atlas_Config), intent(in)           :: config
  type(atlas_functionspace_StructuredColumns) :: functionspace
  type(atlas_Nabla), intent(in), optional :: nabla

  integer :: ngp_halo
  integer :: ngp
  integer :: nlev

  call this%geometry%setup( config, functionspace )

  if( this%needs_nabla(config) ) then
    if( present(nabla) .and. .not. nabla%is_null() ) then
       this%nabla = nabla
    else
       call ABORT("optional nabla argument required for given configuration")
    endif
  endif
  
  call this%departurepoint_calculation%setup( config, this%geometry, this%nabla )

  ngp_halo = this%geometry%ngptot
  ngp = this%geometry%ngp
  nlev = this%geometry%nlev
  
  !-------------------------------------------------------------------------
  ! departure point fields
  !-------------------------------------------------------------------------
  this%departure_points = atlas_FieldSet()
  call this%departure_points%add( atlas_Field( name="dplonfield",  kind=atlas_real(wp), shape=[nlev,ngp] ) )
  call this%departure_points%add( atlas_Field( name="dplatfield",  kind=atlas_real(wp), shape=[nlev,ngp] ) )
  call this%departure_points%add( atlas_Field( name="dpvertfield", kind=atlas_real(wp), shape=[nlev,ngp] ) )
  call set_units( this%departure_points%field("dplonfield"), "radians" )
  call set_units( this%departure_points%field("dplatfield"), "radians" )

  !-------------------------------------------------------------------------
  ! Wind and tracer fields
  !-------------------------------------------------------------------------
  this%wind_t0      = functionspace%create_field( &
      name="wind(t)",    kind=atlas_real(wp), variables=3, type="vector" )
  this%wind_t9      = functionspace%create_field( &
      name="wind(t-dt)", kind=atlas_real(wp), variables=3, type="vector" )
  this%tracer_field = functionspace%create_field( &
      name="tracers",    kind=atlas_real(wp), variables=this%ntrac, type="scalar" )

end subroutine

!------------------------------------------------------------------------------

#if defined _CRAYFTN || defined __INTEL_COMPILER
! Note: workaround no longer required with cce/8.7 onwards, but doesn't harm either
!       Workaround also required for intel/18.0.1
#define CRAY_WORKAROUND recursive
#else
#define CRAY_WORKAROUND
#endif
CRAY_WORKAROUND subroutine execute(this,args)
  !-------------------------------------------------------------------------
  class(dwarf_sladv), intent(inout) :: this
  class(dwarf_args), intent(inout) :: args
  !-------------------------------------------------------------------------
  type(Timer_type)    :: timer
  character(len=1024) :: string
  type(atlas_Trace)   :: trace
  !-------------------------------------------------------------------------
  
  trace = atlas_Trace( __FILENAME__, __LINE__, "Execute SLADV" )
  call timer%start()

  !-------------------------------------------------------------------------
  ! Transfer data from model
  !-------------------------------------------------------------------------
  call this%preprocess(args)

  !-------------------------------------------------------------------------
  ! timestep from t to t+dt
  !-------------------------------------------------------------------------

  call this%departurepoint_calculation%execute( &
    this%istep,                                 &
    this%dt,                                    &
    this%wind_t9,                               &
    this%wind_t0,                               &
    this%departure_points                       &
  )

  call interpolate_tracers_inplace( this%tracer_field, this%departure_points )

  !-------------------------------------------------------------------------
  ! Transfer data to model
  !-------------------------------------------------------------------------
  call this%postprocess(args)

  !-------------------------------------------------------------------------
  call trace%final()

contains


  subroutine interpolate_tracers_inplace(tracer_field,departure_points)
    type(atlas_Config) :: config
    type(atlas_interpolation) :: interpolation
    type(atlas_Trace) :: trace
    type(atlas_Field) :: tracer_field
    type(atlas_FieldSet) :: departure_points
    trace = atlas_Trace( __FILENAME__, __LINE__, "atlas_tracer_interpolation")

    config = atlas_Config()
    call config%set("type",this%interpolation_method//"3D")
    call config%set("limiter",this%interpolation_limiter)
    call config%set("matrix_free",.true.)

    ! uses field "dplonfield", "dplatfield", and "dpvertfield"
    interpolation = atlas_interpolation(config, this%geometry%fs_structuredcolumns,departure_points)

    call this%tracer_field%halo_exchange()
    
    if( this%arrival_tracer_field%is_null() ) then
      this%arrival_tracer_field = this%geometry%fs_structuredcolumns%create_field( &
        name="arrival_tracers", kind=atlas_real(wp), variables=this%ntrac )
    endif

    call interpolation%execute( tracer_field, this%arrival_tracer_field )

    call update_tracers( tracer_field, this%arrival_tracer_field )

    call config%final()
    call interpolation%final()
    call trace%final()
  end subroutine

  subroutine update_tracers( tracer_field, arrival_tracer_field )
    type(atlas_Field) :: tracer_field, arrival_tracer_field

    real(wp), pointer :: arrival_tracers(:,:,:), tracers(:,:,:)
    integer(ip) :: jnode

    call tracer_field%data( tracers )
    call arrival_tracer_field%data( arrival_tracers )

    !$OMP  PARALLEL DO SCHEDULE(STATIC) &
    !$OMP& PRIVATE(jnode)
    do jnode=1,this%geometry%ngp
        tracers(:,:,jnode) = arrival_tracers(:,:,jnode)
    enddo
    !$OMP END PARALLEL DO

    call tracer_field%set_dirty()
  end subroutine

end subroutine

!------------------------------------------------------------------------------

CRAY_WORKAROUND subroutine final(this)
  class(dwarf_sladv), intent(inout) :: this
  call this%final()
  call this%geometry%final()
  call this%departurepoint_calculation%final()
  call this%departure_points%final()
  call this%tracer_field%final()
  call this%wind_t9%final()
  call this%wind_t0%final()
  call this%nabla%final()
  call this%arrival_tracer_field%final()
end subroutine

!------------------------------------------------------------------------------

end module
