! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

program test_drhook
use yomhook
use oml_mod
use mpl_module
#ifdef WITH_FCKIT
use fckit_module
#endif
implicit none
real(jphook) :: zhook_handle
integer :: nproc, myproc

call test_init()

if (lhook) call dr_hook('test_drhook',0,zhook_handle)
call function_1()
if (lhook) call dr_hook('test_drhook',1,zhook_handle)

call test_end()

contains

subroutine function_3
implicit none
#include "abor1.intfb.h"
real(jphook) :: zhook_handle
real(jphook) :: automatic_array(10,10)
if (lhook) call dr_hook('function_3',0,zhook_handle)
automatic_array = 2
if( do_abort() .and. oml_my_thread() == oml_max_threads() .and. ( myproc == nproc ) ) then
    call abor1fl("test_drhook_fortran.F90",__LINE__,"problem in function_3")
endif
if (lhook) call dr_hook('function_3',1,zhook_handle)    
end subroutine

subroutine function_2
implicit none
real(jphook) :: zhook_handle
real(jphook) :: automatic_array(10,10)
if (lhook) call dr_hook('function_2',0,zhook_handle)
automatic_array = 1
!$OMP PARALLEL
call function_3()
!$OMP END PARALLEL
call barrier()
if (lhook) call dr_hook('function_2',1,zhook_handle)    
end subroutine

subroutine function_1
implicit none
real(jphook) :: zhook_handle
real(jphook) :: automatic_array(10,10)
automatic_array = 3
if (lhook) call dr_hook('function_1',0,zhook_handle)
call function_2()
if (lhook) call dr_hook('function_1',1,zhook_handle)
end subroutine

! -----------------------------------------------------------------------------------------
! Initialization and Finalization
! -----------------------------------------------------------------------------------------

subroutine test_init
use ec_env_mod, only : ec_setenv

if( mpl() ) then
    call mpl_init(ldinfo=.true.)
    nproc = mpl_numproc
    myproc = mpl_rank
#ifdef WITH_FCKIT
elseif( fckit() ) then
    call fckit_main%init()
    nproc = fckit_mpi%size()
    myproc = fckit_mpi%rank() + 1
#else
else
    nproc = 1
    myproc = 1
#endif
endif

call ec_setenv("DR_HOOK",        "1",   overwrite=.true.)
call ec_setenv("DR_HOOK_OPT",    "PROF",overwrite=.false.)
call ec_setenv("DR_HOOK_SILENT", "1",   overwrite=.false.)

end subroutine

subroutine test_end
call dr_hook_end()
if( mpl() ) then
  call mpl_end(ldmeminfo=.false.)
endif
end subroutine

! -----------------------------------------------------------------------------------------
! Utility functions, checking environment variables etc
! -----------------------------------------------------------------------------------------

function do_abort() result(labort)
    logical :: labort
    character(len=512) :: env
    call get_environment_variable("ABORT",env)
    if( env == '0' ) then
        labort = .false.
    else
        labort = .true.
    endif
end function

function mpl() result(lmpl)
    logical :: lmpl
    character(len=512) :: env
    call get_environment_variable("MPL",env)
    if( env == '0' ) then
        lmpl = .false.
    else
        lmpl = .true.
    endif
end function

function fckit() result(lfckit)
    logical :: lfckit
    character(len=512) :: env
    call get_environment_variable("FCKIT",env)
    if( env == '0' ) then
        lfckit = .false.
    else
        lfckit = .true.
    endif
end function

subroutine barrier()
    if( mpl() ) then
        call mpl_barrier()
#ifdef WITH_FCKIT
    elseif( fckit() ) then
        call fckit_mpi%barrier()
#endif
    endif
end subroutine

! -----------------------------------------------------------------------------------------

end program
