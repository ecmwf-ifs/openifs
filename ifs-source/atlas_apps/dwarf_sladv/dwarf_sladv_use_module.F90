! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

module dwarf_sladv_use_module

use fckit_module
use atlas_module, only: atlas_Field, atlas_Metadata, atlas_Trace, atlas_Config, &
&                       atlas_FieldSet, atlas_Config, atlas_Interpolation,      &
&                       atlas_real, atlas_Functionspace, atlas_Nabla,           &
&                       atlas_functionspace_StructuredColumns
use dwarf_module
use dwarf_sladv_kind_module, only: i8, dp, ip, wp
use dwarf_sladv_parameters_module, only: rpi, rndoff
use dwarf_sladv_geometry_module, only: Geometry_type

implicit none
public

! Simple timer facility
type :: Timer_type
private
  integer(i8) :: clck_rate = 1
  integer(i8) :: counted = 0
  logical    :: paused  = .True.
  integer(i8) :: clck_counts_start
  integer(i8) :: clck_counts_stop
contains
  procedure, public :: start   => Timer_start
  procedure, public :: pause   => Timer_pause
  procedure, public :: resume  => Timer_resume
  procedure, public :: elapsed => Timer_elapsed
end type Timer_type

interface print_checksum
  module procedure print_checksum_field
  module procedure print_checksum_fieldset
end interface

contains

! -----------------------------------------------------------------------------
! Compute elapsed time
! -----------------------------------------------------------------------------
function Timer_elapsed(self) result(time)
  class(Timer_type), intent(inout) :: self
  real(dp) :: time

  if (.not. self%paused) then
    call system_clock(self%clck_counts_stop, self%clck_rate)
    time = (self%counted + self%clck_counts_stop - &
          & self%clck_counts_start) / real(self%clck_rate)
  else if (self%counted .ge. 0) then
    time = self%counted / real(self%clck_rate)
  else
    time = 0.
  endif

end function
! -----------------------------------------------------------------------------



! -----------------------------------------------------------------------------
! Starts the clock
! -----------------------------------------------------------------------------
subroutine Timer_start(self)
  class(Timer_type), intent(inout) :: self

  call system_clock(self%clck_counts_start, self%clck_rate)
  self%paused  = .False.
  self%counted = 0
end subroutine
! -----------------------------------------------------------------------------



! -----------------------------------------------------------------------------
! Pauses the clock
! -----------------------------------------------------------------------------
subroutine Timer_pause(self)
  class(Timer_type), intent(inout) :: self

  if( .not. self%paused ) then
    call system_clock(self%clck_counts_stop, self%clck_rate)
    self%counted = self%counted + self%clck_counts_stop - self%clck_counts_start
    self%paused  = .True.
  endif

end subroutine
! -----------------------------------------------------------------------------



! -----------------------------------------------------------------------------
! Resume the clock
! -----------------------------------------------------------------------------
subroutine Timer_resume(self)
  class(Timer_type), intent(inout) :: self

  call system_clock ( self%clck_counts_start, self%clck_rate )
  self%paused = .False.

end subroutine Timer_resume
! -----------------------------------------------------------------------------


subroutine print_checksum_field( file, line, prefix, field )
  use fckit_module, only : fckit_log
  use atlas_module
  character(len=*)  , intent(in) :: file
  integer           , intent(in) :: line
  type(atlas_Field) , intent(in) :: field
  type(atlas_functionspace_StructuredColumns) :: fs
  character(len=*)    :: prefix
  character(len=1024) :: cat_prefix
  character(len=1024) :: string
  fs = field%functionspace()
  write(cat_prefix,'(A,A,I0,A,A)') file," +",line," ",prefix
  write(string,'(A,A,A)') trim(cat_prefix)," ",fs%checksum(field)
  call fckit_log%info(string)
  call fs%final()
end subroutine

subroutine print_checksum_fieldset( file, line, prefix, fields )
  use fckit_module, only : fckit_log
  use atlas_module
  character(len=*)  , intent(in) :: file
  integer           , intent(in) :: line
  type(atlas_FieldSet) , intent(in) :: fields
  type(atlas_functionspace_StructuredColumns) :: fs
  character(len=*)    :: prefix
  character(len=1024) :: cat_prefix
  character(len=1024) :: string
  type(atlas_Field) :: field1
  field1 = fields%field(1)
  fs = field1%functionspace()
  write(cat_prefix,'(A,A,I0,A,A)') file," +",line," ",prefix
  write(string,'(A,A,A)') trim(cat_prefix)," ",fs%checksum(fields)
  call fckit_log%info(string)
  call fs%final()
  call field1%final()
end subroutine

subroutine set_units( field, value )
  type(atlas_Field) :: field
  character(len=*) :: value
  type(atlas_Metadata) :: metadata
  metadata = field%metadata()
  call metadata%set("units",value)
  call metadata%final()
end subroutine

!------------------------------------------------------------------------------

end module
