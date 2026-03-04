! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

program abort_test
use oml_mod
use mpl_module
use yomhook, only : lhook, dr_hook_init
use fckit_module
implicit none

#include "abor1.intfb.h"

interface
   subroutine my_namespace__function_1() bind(C)
   end subroutine
end interface

call test_init()
if( fckit_mpi%rank() == fckit_mpi%size()-1 ) then
    call my_namespace__function_1()
endif
call test_end()

contains

subroutine test_init()
    use fckit_module, only : fckit_main, fckit_signal, fckit_exception
    use ec_env_mod, only : ec_setenv
    implicit none


    if( mpl() ) then
        call mpl_init(LDINFO=.TRUE.)
    endif

    call fckit_main%initialise()

    if( do_abor1() ) then
        call set_abor1_exception_handler()
    endif

    call ec_setenv("DR_HOOK","1",overwrite=.false.)
    call ec_setenv("DR_HOOK_SILENT","1",overwrite=.false.)
    call dr_hook_init()
    if (.not.lhook) then
        call fckit_signal%set_handler(fckit_signal%SIGTRAP())
        call fckit_signal%set_handler(fckit_signal%SIGABRT())
    endif
    
end subroutine

subroutine test_end()
    if( mpl() ) then
        call mpl_barrier()
        call mpl_end(LDMEMINFO=.FALSE.)    
    endif
end subroutine

function mpl() result(lmpl)
    logical :: lmpl
    character(len=512) :: env
    CALL get_environment_variable('MPL',env)
    if( env == '0' ) then
        lmpl = .false.
    else
        lmpl = .true.
    endif
end function

function do_abor1() result(labor1)
    logical :: labor1
    character(len=512) :: env = ' '
    CALL get_environment_variable('ABOR1',env)
    if( env == '0' .or. len_trim(env) == 0 ) then
        labor1 = .false.
    else
        labor1 = .true.
    endif
end function


end program
