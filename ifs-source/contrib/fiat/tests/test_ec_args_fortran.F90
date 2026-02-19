! (C) Copyright 2021- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!
! Testing of ec_args; expected to be launched with 3 arguments:
!    <program> arg1 arg2 arg3
!

#define FAIL(msg) call fail_impl(msg,__LINE__)

program test_ec_args
use iso_c_binding
use ec_args_mod
implicit none
integer(c_int) :: argc, iarg

argc = ec_argc()
if( argc /= 0 ) FAIL("argc should be zero before initialisation")

! Setup ec_args --> computes command line arguments in Fortran and stores in C for C to use as well
call ec_args()

argc = ec_argc()
if( argc == 0 ) FAIL("argc should be non-zero after initialisation")

write(0,*) "program name = ", ec_argv(0)

if( argc /= 4 ) FAIL("4 arguments expected")
if( ec_argv(1) /= "arg1" ) FAIL("unexpected value 'arg1'")
if( ec_argv(2) /= "arg2" ) FAIL("unexpected value 'arg2'")
if( ec_argv(3) /= "arg3" ) FAIL("unexpected value 'arg3'")

contains
subroutine fail_impl(msg,line)
    character(*) :: msg
    integer :: line
    write(0,'(A,I0,A)') "TEST FAILED in test_ec_args.F90 @ line ",line," :"
    write(0,*) msg
    stop 1
end subroutine
end program
