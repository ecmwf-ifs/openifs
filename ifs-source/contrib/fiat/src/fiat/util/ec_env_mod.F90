! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

module ec_env_mod

!**** Interface to ec_env environment handling
!
!     Purpose.
!     --------
!     Fortran 90 Interface to setting and getting environment variables
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
public :: ec_putenv, ec_getenv, ec_setenv, ec_numenv, ec_environ

interface
    !! ISO-C bindings from file ec_env.c
    subroutine ec_putenv_overwrite_bind_c(env,env_len) bind(C)
      use, intrinsic :: iso_c_binding, only : c_int, c_char
      character(kind=c_char),     intent(in) :: env(*)
      integer(kind=c_int), value, intent(in) :: env_len
    end subroutine
    subroutine ec_putenv_nooverwrite_bind_c(env,env_len) bind(C)
      use, intrinsic :: iso_c_binding, only : c_int, c_char
      character(kind=c_char),     intent(in) :: env(*)
      integer(kind=c_int), value, intent(in) :: env_len
    end subroutine
    subroutine ec_getenv_bind_c(key,value,key_len,value_len) bind(C)
      use, intrinsic :: iso_c_binding, only : c_int, c_char
      character(kind=c_char),     intent(in)    :: key(*)
      character(kind=c_char),     intent(inout) :: value(*)
      integer(kind=c_int), value, intent(in)    :: key_len
      integer(kind=c_int), value, intent(in)    :: value_len
    end subroutine
    subroutine ec_numenv_bind_c(value) bind(C)
      use, intrinsic :: iso_c_binding, only : c_int
      integer(kind=c_int), intent(out) :: value
    end subroutine
    subroutine ec_environ_bind_c(i,value,value_len) bind(C)
      use, intrinsic :: iso_c_binding, only : c_int, c_char
      integer(kind=c_int),        intent(in)    :: i
      character(kind=c_char),     intent(inout) :: value(*)
      integer(kind=c_int), value, intent(in)    :: value_len
    end subroutine
end interface

contains

subroutine ec_putenv(env, overwrite)
  !! Set environment. Add optional argument 'OVERWRITE=.false.' to avoid overwriting if already set.
  !! Example:
  !!   call ec_putenv("DR_HOOK=1")                    ! Forces overwrite (default: OVERWRITE=.TRUE.)
  !!   call ec_putenv("DR_HOOK=0",OVERWRITE=.FALSE)   ! Will have no effect as DR_HOOK=1 already exists
  !!
  character(len=*), intent(in) :: env
  logical, optional, intent(in) :: overwrite ! assume true if not present
  if( present(overwrite) ) then
    if( overwrite ) then
      call ec_putenv_overwrite_bind_c(env,len(env))
    else
      call ec_putenv_nooverwrite_bind_c(env,len(env))  
    endif
  else
    call ec_putenv_overwrite_bind_c(env,len(env))
  endif
end subroutine

subroutine ec_setenv(key, value, overwrite)
  !! Set environment, with non-optional 'overwrite' argument which can be used to force or avoid overwriting if already set.
  !! Example:
  !!   call ec_setenv("DR_HOOK","1",OVERWRITE=.TRUE.)   ! Forces overwrite
  !!   call ec_putenv("DR_HOOK","0",OVERWRITE=.FALSE)   ! Will have no effect as DR_HOOK=1 already exists
  !!
  character(len=*), intent(in) :: key
  character(len=*), intent(in) :: value
  logical, intent(in) :: overwrite
  character(len=:), allocatable :: env
  env = trim(key)//'='//trim(value)
  if( overwrite ) then
    call ec_putenv_overwrite_bind_c(env,len(env))
  else
    call ec_putenv_nooverwrite_bind_c(env,len(env))  
  endif
end subroutine

subroutine ec_getenv(key, value)
  !! Get environment variable by key. If key is not available, value argument will be filled with spaces.
  !! Example:
  !!   CHARACTER(len=255) :: CENV_DR_HOOK
  !!   call ec_getenv("DR_HOOK",CENV_DR_HOOK)
  !!
  character(len=*), intent(in) :: key
  character(len=*), intent(inout) :: value
  call ec_getenv_bind_c(key,value,len(key),len(value))
end subroutine

function ec_numenv() result(value)
  !! Return number of environment variables
  !!
  use, intrinsic :: iso_c_binding, only : c_int
  integer(kind=c_int) :: value
  call ec_numenv_bind_c(value)
end function

subroutine ec_environ(i,value)
  !! Set value to entry of the unix "environ" variable with index 1..ec_numenv
  !! The content of value will be of form "KEY=VALUE", as in ec_putenv
  use, intrinsic :: iso_c_binding, only : c_int
  integer(kind=c_int), intent(in) :: i
  character(len=*), intent(inout) :: value
  call ec_environ_bind_c(i,value,len(value))
end subroutine

end module
