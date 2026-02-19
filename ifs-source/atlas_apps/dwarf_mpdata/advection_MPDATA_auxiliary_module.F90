! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

! =============================================================================
! auxiliary_module
! This module contains strictly auxiliary subroutines
! =============================================================================

module advection_MPDATA_auxiliary_module

use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double

use atlas_module
use fckit_module, only : fckit_log, fckit_mpi_comm, fckit_main, fckit_resource

#ifdef HAVE_OMP
  use omp_lib
#endif

implicit none

private

character(len=2048), public :: string ! use for logging

type(fckit_mpi_comm), save, public :: mpi

! Single precision
integer, parameter, public :: sp = c_float  ! selected_real_kind(4,2)

! Double precision
integer, parameter, public :: dp = c_double ! selected_real_kind(13,300)

! Working precision
#ifdef HAVE_SINGLE_PRECISION
integer, parameter, public :: wp = sp
#else
integer, parameter, public :: wp = dp
#endif

! pi and commonly used divisions
real(wp), parameter, public :: rpi     = 2.0_wp * asin(1.0_wp)
real(wp), parameter, public :: rpi_2   = 0.5_wp * rpi
real(wp), parameter, public :: rpi_4   = 0.25_wp * rpi
real(wp), parameter, public :: rad2deg = 180._wp / rpi
real(wp), parameter, public :: deg2rad = rpi / 180._wp

! index helpers for vectors
integer, parameter, public :: MXX = 1
integer, parameter, public :: MYY = 2
integer, parameter, public :: MZZ = 3

! Simple timer facility
public :: Timer_type
type :: Timer_type
private
  integer(8) :: clck_counts_start
  integer(8) :: clck_counts_stop
  integer(8) :: clck_rate
  integer(8) :: counted = 0
  logical    :: paused  = .True.
contains
  procedure, public :: start   => Timer_start
  procedure, public :: pause   => Timer_pause
  procedure, public :: resume  => Timer_resume
  procedure, public :: elapsed => Timer_elapsed
end type Timer_type


public :: workdir
public :: dwarf_D_advection_MPDATA_initialize
public :: omp_get_max_threads

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

end function Timer_elapsed
! -----------------------------------------------------------------------------



! -----------------------------------------------------------------------------
! Starts the clock
! -----------------------------------------------------------------------------
subroutine Timer_start(self)
class(Timer_type), intent(inout) :: self

call system_clock(self%clck_counts_start, self%clck_rate)
self%paused  = .False.
self%counted = 0

end subroutine Timer_start
! -----------------------------------------------------------------------------



! -----------------------------------------------------------------------------
! Pauses the clock
! -----------------------------------------------------------------------------
subroutine Timer_pause(self)
class(Timer_type), intent(inout) :: self

call system_clock(self%clck_counts_stop, self%clck_rate)
self%counted = self%counted + self%clck_counts_stop - self%clck_counts_start
self%paused  = .True.

end subroutine Timer_pause
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



! =============================================================================
! Configuration of OMP
! =============================================================================
#ifndef HAVE_OMP
function omp_get_max_threads()
  integer :: omp_get_max_threads
  omp_get_max_threads = 0
end function
#endif
! =============================================================================


! -----------------------------------------------------------------------------
! Defines the working directory
! -----------------------------------------------------------------------------
function workdir()
  character(len=1024) :: workdir
  call get_environment_variable('PWD',workdir)
end function
! -----------------------------------------------------------------------------



! =============================================================================
! Print versions of the libraries used
! =============================================================================
subroutine dwarf_D_advection_MPDATA_initialize

character(len=:), allocatable :: name

mpi = fckit_mpi_comm()
call atlas_library%initialise()

call fckit_main%name(name)

call fckit_log%info("")
call fckit_log%info("+---------------------------------------------------+")
call fckit_log%info("| dwarf_D_advection_MPDATA                          |")
call fckit_log%info("+---------------------------------------------------+")
call fckit_log%info("    Atlas version ("//trim(atlas_version())//") "//&
& "git-sha1 "//trim(atlas_git_sha1_abbrev(13)) )
call fckit_log%info("    eckit version ("//trim(eckit_version())//") "//&
& "git-sha1 "//trim(eckit_git_sha1_abbrev(13)) )
call fckit_log%info("")
write(string,'(A,A)')  "    Program     = ", trim(name); call fckit_log%info(string)
write(string,'(A,A)')  "    workdir     = ", trim(workdir()); call fckit_log%info(string)
write(string,'(A,I0)') "    MPI Tasks   = ", mpi%size();
call fckit_log%info(string)
write(string,'(A,I0)') "    OMP Threads = ", omp_get_max_threads();
call fckit_log%info(string)
call fckit_log%info("+---------------------------------------------------+")
call fckit_log%info("",flush=.true.)


end subroutine dwarf_D_advection_MPDATA_initialize
! =============================================================================

end module advection_MPDATA_auxiliary_module
