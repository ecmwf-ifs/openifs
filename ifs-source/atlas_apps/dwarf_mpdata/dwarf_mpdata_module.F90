! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#define __FILENAME__ "dwarf_mpdata_module.F90"

module dwarf_mpdata_module

use dwarf_module,                      only : dwarf,       & 
                                            & dwarf_args
use advection_MPDATA_MPDATA_module,    only : mpdata_type
use advection_MPDATA_auxiliary_module, only : wp
use atlas_module,                      only : atlas_Field, &
                                            & atlas_Trace

implicit none
private

private :: dwarf
private :: dwarf_args
private :: mpdata_type
private :: atlas_Field
private :: wp


!------------------------------------------------------------------------------

type, abstract, public, extends( dwarf ) :: dwarf_mpdata
  ! Needs setting up in concrete class
  integer(4) :: ntrac
  integer(4) :: istep
  real(wp)   :: dt
  integer(4) :: nsteps
  type(mpdata_type)       :: mpdata
  type(atlas_Field)       :: tracer_field
  type(atlas_Field)       :: rho_field
  type(atlas_Field)       :: rhojacfield
  type(atlas_Field)       :: rhofacfield
  type(atlas_Field)       :: rhodtfield
  type(atlas_Field)       :: velfield
  type(atlas_Field)       :: vhatfield
  type(atlas_Field)       :: wnfield
  type(atlas_Field)       :: vxyfield
  type(atlas_Field)       :: vnfield
  type(atlas_Field)       :: jacobian_determinant
contains
  procedure, public :: execute
  procedure(execute), deferred, public :: preprocess
  procedure(execute), deferred, public :: postprocess
  procedure, public :: final => dwarf_mpdata_final

  procedure, public :: dwarf_mpdata_final
  procedure, public :: dwarf_mpdata_setup
end type

public :: print_checksum
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------

subroutine dwarf_mpdata_setup( this )
  class(dwarf_mpdata) :: this
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
  use fckit_module,                      only : fckit_log
  use advection_MPDATA_auxiliary_module, only : Timer_type
  !-------------------------------------------------------------------------
  class(dwarf_mpdata), intent(inout) :: this
  class(dwarf_args), intent(inout) :: args
  !-------------------------------------------------------------------------
  real(wp), pointer   :: jdet(:,:)
  real(wp), pointer   :: rhojac(:,:)
  real(wp), pointer   :: rhofac(:,:)
  real(wp), pointer   :: rhodt(:,:)
  real(wp), pointer   :: rhophy(:,:)
  real(wp), pointer   :: vel(:,:,:)
  real(wp), pointer   :: vhat(:,:,:)
  real(wp), pointer   :: vxy(:,:,:)
  real(wp), pointer   :: vn(:,:)
  real(wp), pointer   :: wn(:,:)
  real(wp)            :: zrhojac
  integer             :: jnode, jlev, nb_nodes, nb_levels, istep, itrac
  type(Timer_type)    :: timer
  character(len=1024) :: string
  type(atlas_Trace)   :: trace, trace_timesteps
  !-------------------------------------------------------------------------
  trace = atlas_Trace(__FILENAME__,__LINE__,"dwarf_mpdata::execute")
  call fckit_log%info( "Executing MPDATA..." )
  call timer%start()
  !-------------------------------------------------------------------------
  ! Transfer data from model
  !-------------------------------------------------------------------------
  call this%preprocess(args)

  !-------------------------------------------------------------------------
  ! timestep from t to t+dt
  !-------------------------------------------------------------------------

  trace_timesteps = atlas_Trace(__FILENAME__,__LINE__,"dwarf_mpdata::timesteps")
 
  call this%velfield             % data( vel )
  call this%vhatfield            % data( vhat )
  call this%vxyfield             % data( vxy )
  call this%vnfield              % data( vn )
  call this%wnfield              % data( wn )

  call this%jacobian_determinant % data( jdet )
  call this%rhojacfield          % data( rhojac )
  call this%rhofacfield          % data( rhofac )
  call this%rhodtfield           % data( rhodt )
  call this%rho_field            % data( rhophy )

  nb_levels = this%rho_field%shape(1)
  nb_nodes  = this%rho_field%shape(2)

  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev)
  do jnode=1,nb_nodes
    do jlev=1,nb_levels
      vel(1,jlev,jnode) = vhat(1,jlev,jnode)*jdet(jlev,jnode)
      vel(2,jlev,jnode) = vhat(2,jlev,jnode)*jdet(jlev,jnode)
      vel(3,jlev,jnode) = vhat(3,jlev,jnode)*jdet(jlev,jnode)
    enddo
  enddo
  !$OMP END PARALLEL DO

  do istep=1,this%nsteps

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev)
    do jnode=1,nb_nodes
      do jlev=1,nb_levels
        rhojac(jlev,jnode) = rhophy(jlev,jnode)*jdet(jlev,jnode)
        rhofac(jlev,jnode) = 1._wp
        rhodt(jlev,jnode)  = 0._wp
      enddo
    enddo
    !$OMP END PARALLEL DO

    !call print_checksum(__FILE__,__LINE__,this%mpdata,"velfield",this%velfield)

    call this%mpdata%advectors_in_edges( &
        this%velfield,                   &
        this%vxyfield,                   &
        this%wnfield                     &
    )

    call this%mpdata%face_normals(       &
        this%vxyfield,                   &
        this%vnfield                     &
    )

    call this%mpdata%execute_density(    &
        this%dt,                         &
        this%rho_field,                  &
        this%jacobian_determinant,       &
        this%rhofacfield,                &
        this%rhodtfield,                 &
        this%vnfield,                    &
        this%wnfield                     &
    )

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev,zrhojac)
    do jnode=1,nb_nodes
      do jlev=1,nb_levels
        zrhojac = rhophy(jlev,jnode)*jdet(jlev,jnode)
        rhofac(jlev,jnode) = rhojac(jlev,jnode)/zrhojac
        rhodt(jlev,jnode) = (zrhojac-rhojac(jlev,jnode))/this%dt
      enddo
    enddo
    !$OMP END PARALLEL DO

    do itrac=1, this%ntrac
      call this%mpdata%execute(          &
          this%dt,                       &
          itrac,                         &
          this%tracer_field,             &
          this%rhojacfield,              &
          this%rhofacfield,              &
          this%rhodtfield,               &
          this%vnfield,                  &
          this%wnfield                   &
      )
    enddo

!    call print_checksum(__FILE__,__LINE__,"tracer_field",this%tracer_field)

  enddo

  call trace_timesteps%final()

  !-------------------------------------------------------------------------
  ! Transfer data to model
  !-------------------------------------------------------------------------
  call this%postprocess(args)

  !-------------------------------------------------------------------------
  write(string,'(A,F5.2,A)') "Executing MPDATA... done  ( ",timer%elapsed()," s )"
  call fckit_log%info( string )
  call trace%final()

end subroutine

!------------------------------------------------------------------------------

subroutine dwarf_mpdata_final(this)
  class(dwarf_mpdata), intent(inout) :: this
  call this%tracer_field%final()
  call this%rho_field%final()
  call this%rhojacfield%final()
  call this%rhofacfield%final()
  call this%rhodtfield%final()
  call this%velfield%final()
  call this%vhatfield%final()
  call this%wnfield%final()
  call this%vxyfield%final()
  call this%vnfield%final()
  call this%jacobian_determinant%final()
end subroutine

!------------------------------------------------------------------------------

subroutine print_checksum( file, line, prefix, field )
  use fckit_module, only : fckit_log
  use atlas_module
  character(len=*)  , intent(in) :: file
  integer           , intent(in) :: line
  type(atlas_Field) , intent(in) :: field
  type(atlas_functionspace_NodeColumns) :: fs
  character(len=*)    :: prefix
  character(len=1024) :: cat_prefix
  character(len=1024) :: string
  fs = field%functionspace()
  write(cat_prefix,'(A,A,I0,A,A)') file," +",line," ",prefix
  write(string,'(A,A,A)') trim(cat_prefix)," ",fs%checksum(field)
  call fckit_log%info(string)
  call fs%final()
end subroutine


end module
