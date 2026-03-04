! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

! This program aborts from an OpenMP parallel region
! It can be called within a MPI-parallel context when "MPL=1" is set in the environment
! To enable DR_HOOK, set "DR_HOOK=1" in environment
program test_abor1
use oml_mod
use mpl_module
use yomhook
implicit none

#include "abor1.intfb.h"

call test_init()

if( MPL_RANK <= 1 ) write(0,'(A,I0)') "OML_MAX_THREADS = ", OML_MAX_THREADS()

!$OMP PARALLEL
if( OML_MY_THREAD() == OML_GET_NUM_THREADS() .AND. (MPL_RANK == MPL_NUMPROC .OR. MPL_NUMPROC < 0)) then
    if( do_abort() ) call abor1fl("test_abort.F90",__LINE__,"aborting from OpenMP parallel region")
endif
!$OMP END PARALLEL

call test_end()

contains

    subroutine test_init()
        use ec_env_mod, only : ec_setenv
        implicit none

        ! Only enables MPL when environment MPL=1
        if( mpl() ) then
            call mpl_init(LDINFO=.FALSE.)
        endif

        ! Only enables DR_HOOK when environment DR_HOOK=1
        call ec_setenv("DR_HOOK_SILENT","1",overwrite=.false.)
        call dr_hook_init() 
    end subroutine

    subroutine test_end()
        ! Should not reach here unless "ABORT=0" in environment
        if( mpl() ) then
            call mpl_barrier()
            call mpl_end(LDMEMINFO=.FALSE.)    
        endif
    end subroutine

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
        if( env == '1' ) then
            lmpl = .true.
        else
            lmpl = .false.
        endif
    end function


end program
