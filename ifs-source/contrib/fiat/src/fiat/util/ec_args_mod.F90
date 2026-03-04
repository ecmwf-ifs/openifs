! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

module ec_args_mod

!**** Interface to ec_args command-line handling
!
!     Purpose.
!     --------
!     Fortran 90 Interface to storing and retrieving command line arguments
!     for the C ec_args API.
!  
!     An example C program:
!
!         #include "ec_args.h"  
!         int main( int argc, char* argv[] ) {
!             ec_args(argc,argv);
!             int num_args     = ec_argc();
!             const char* name = ec_argv()[0];
!         }
!
!     An example Fortran program to do the same:
!     
!         program main
!             use ec_args_mod
!             integer :: num_args
!             character(len=:), allocatable :: name
!             call ec_args()
!             num_args = ec_argc()
!             name     = ec_argv(0)
!         end program  
!
!     Author.
!     -------
!        W.Deconinck, ECMWF
!
!     Modifications.
!     --------------
!        Original: 2021-05-18
!
!     ------------------------------------------------------------------
implicit none

private
public :: ec_argc, ec_argv, ec_args

#define MAX_ARG_LEN 1024
#define EC_MAX_ARGS 512
 !! Matches value in ec_args.c

interface
    function ec_argc() bind(C,name="ec_argc") result(argc)
      use, intrinsic :: iso_c_binding, only : c_int
      integer(c_int) :: argc
    end function
    function ec_argv_bindc() bind(C,name="ec_argv") result(argv)
      use, intrinsic :: iso_c_binding, only : c_ptr
      type(c_ptr) :: argv
    end function
    subroutine ec_args_bindc(argc,argv) BIND(C,NAME="ec_args")
      use, intrinsic :: iso_c_binding, only : c_int, c_ptr
      integer(c_int), value :: argc
      type(c_ptr), dimension(*) :: argv
    end subroutine
end interface

contains

function ec_argv(iarg) result(argv)
    use, intrinsic :: iso_c_binding
    implicit none
    character(len=:), allocatable :: argv
    integer(c_int), intent(in) :: iarg
    integer(c_int) :: argc
    type(c_ptr) :: argv_cptr
    type(c_ptr), pointer :: argv_cptrs(:)
    argc = ec_argc()
    argv_cptr = ec_argv_bindc()
    call c_f_pointer ( argv_cptr , argv_cptrs, (/argc/) )
    argv = to_string( argv_cptrs(iarg+1), MAX_ARG_LEN )
end function

subroutine ec_args()
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int) :: argc
  type(c_ptr) :: argv(EC_MAX_ARGS)
  if( ec_argc() == 0 ) then
    call read_command_line(argc,argv)
    call ec_args_bindc(argc,argv)
  endif
end subroutine

function to_string(cptr,maxlen) result(string)
  ! Convert string from C (char*) to Fortran
  ! (copied from fckit)
  use, intrinsic :: iso_c_binding
  implicit none
  character(len=:), allocatable :: string
  type(c_ptr) :: cptr
  integer(c_int) :: maxlen
  character(kind=c_char,len=1), pointer :: s(:)
  integer i, nchars
  call c_f_pointer ( cptr , s, (/maxlen/) )
  i = 1
  do
      if (s(i) == c_null_char) exit
      i = i + 1
  enddo
  nchars = i - 1  ! Exclude null character from Fortran string
  allocate( character(len=(nchars)) :: string )
  do i=1,nchars
      string(i:i) = s(i)
  enddo
end function


subroutine read_command_line(argc,argv)
  ! Read command line arguments into argc and argv as in C
  ! (copied from fckit)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), parameter :: CMD_MAX_LEN = MAX_ARG_LEN * EC_MAX_ARGS
  integer(c_int) :: argc
  type(c_ptr) :: argv(:)
  character(kind=c_char,len=1), save, target :: args(CMD_MAX_LEN)
  character(kind=c_char,len=CMD_MAX_LEN), save, target :: cmd
  character(kind=c_char,len=CMD_MAX_LEN) :: arg
  integer(c_int) :: iarg, arglen, pos, ich, argpos
  call get_command(cmd)
  do ich=1,len(cmd)
    if (cmd(ich:ich) == " ") then
      cmd(ich:ich) = c_null_char
      exit
    endif
  enddo
  argv(1) = c_loc(cmd(1:1))
  argc = command_argument_count()+1
  pos = 1
  do iarg=1,argc
    argpos = pos
    call get_command_argument(iarg, arg )
    arglen = len_trim(arg)
    do ich=1,arglen
      args(pos) = arg(ich:ich)
      pos = pos+1
    end do
    args(pos) = c_null_char;  pos = pos+1
    args(pos) = " ";          pos = pos+1
    argv(iarg+1) = c_loc(args(argpos))
  enddo
end subroutine

end module
