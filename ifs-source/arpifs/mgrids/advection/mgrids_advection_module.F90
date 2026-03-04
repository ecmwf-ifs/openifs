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

module mgrids_advection_module
#ifdef WITH_MGRIDS

use dwarf_module , only : dwarf, dwarf_args
use yomrip       , only : trip
use yomdyn       , only : tdyn
use type_geometry, only : geometry
use yomgmv       , only : tgmv
use parkind1     , only : jprb, jpim

implicit none
private

private :: dwarf, dwarf_args
private :: trip
private :: tdyn
private :: geometry
private :: tgmv
private :: jprb, jpim

!------------------------------------------------------------------------------

type, public :: mgrids_advection_config
  character(len=32)    :: ctype
  character(len=32)    :: cgrid
  logical              :: lremesh
  character(len=512)   :: sladv_departurepoint_method
  character(len=512)   :: sladv_interpolation_method
  logical              :: sladv_interpolation_limiter
  integer(jpim)        :: sladv_nhalo
  logical              :: mpdata_lusecfl
  real(jprb)           :: mpdata_zcflmax
  integer(jpim)        :: mpdata_nsubsteps

  logical, private :: initialised = .false.
contains
  procedure, public :: read_namelist
end type

!------------------------------------------------------------------------------

type, public :: mgrids_advection
  class(dwarf), pointer :: implementation => null()
contains
  procedure, public :: active
    !! Notify if mgrids_advection is active ( implementation != null ) 
  procedure, public :: execute
    !! Execute mgrids_advection
  procedure, public :: final
    !! Cleanup mgrids_advection
end type

!------------------------------------------------------------------------------

type, public, extends( dwarf_args ) :: mgrids_advection_args
  type(tgmv), pointer :: yrgmv           => null()
  real(jprb), pointer :: zgmv(:,:,:,:)   => null()
  real(jprb), pointer :: zgfl(:,:,:,:)   => null()
  real(jprb), pointer :: zwrl9(:,:,:)    => null()
  real(jprb), pointer :: zgeo0(:,:,:)    => null()
  real(jprb), pointer :: zrcp0(:,:,:,:)  => null()
  real(jprb), pointer :: zpre0f(:,:,:)   => null()
  real(jprb), pointer :: zcty0(:,:,:,:)  => null()
  real(jprb), pointer :: zgflt1(:,:,:,:) => null()
contains
  final :: mgrids_advection_args__final
end type

!------------------------------------------------------------------------------

interface mgrids_advection
    module function mgrids_advection_create( &
      ydrip,                                 &
      yddyn,                                 &
      ydgeometry ) result(this)
      type(mgrids_advection) :: this
      type(trip),     intent(in), target :: ydrip
      type(tdyn),     intent(in), target :: yddyn
      type(geometry), intent(in), target :: ydgeometry
    end function
end interface

!------------------------------------------------------------------------------

interface
   module function active( this ) result( is_active )
     logical :: is_active
     class(mgrids_advection), intent(in) :: this
   end function
end interface

!------------------------------------------------------------------------------

interface
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
     class(mgrids_advection), intent(inout)         :: this
     type(tgmv),              intent(in),    target :: ydgmv
     real(jprb),              intent(in),    target :: pgmv(:,:,:,:)
     real(jprb),              intent(in),    target :: pgfl(:,:,:,:)
     real(jprb),              intent(in),    target :: pwrl9(:,:,:)
     real(jprb),              intent(in),    target :: pgeo0(:,:,:)
     real(jprb),              intent(in),    target :: prcp0(:,:,:,:)
     real(jprb),              intent(in),    target :: pre0f(:,:,:)
     real(jprb),              intent(in),    target :: pcty0(:,:,:,:)
     real(jprb),              intent(inout), target :: pgflt1(:,:,:,:)
   end subroutine
end interface

!------------------------------------------------------------------------------

interface
   module subroutine final(this)
     class(mgrids_advection) :: this
   end subroutine
end interface

!------------------------------------------------------------------------------
!Not sure if this is required by the standard, but is required for intel/18.0.1
interface
   module subroutine dont_optimize_out(this)
     class(mgrids_advection_args), intent(inout) :: this
   end subroutine
end interface

!------------------------------------------------------------------------------

interface mgrids_advection_args
  module function mgrids_advection_args_dynamic_cast(args) result(this)
    type(mgrids_advection_args) :: this
    class(dwarf_args), intent(inout), target :: args
  end function
end interface

!------------------------------------------------------------------------------
! Avoids that loops editing arguments will be compiled out (observed behaviour on Cray cce/8.6.2)
!------------------------------------------------------------------------------
interface
  module subroutine mgrids_advection_args__final( this )
    type(mgrids_advection_args) :: this
  end subroutine
end interface

!------------------------------------------------------------------------------

interface
  module subroutine read_namelist( this )
    class(mgrids_advection_config) :: this
  end subroutine
end interface 

!------------------------------------------------------------------------------

#endif
end module
