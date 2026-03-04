! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!
! Simple Test program
!
subroutine fail_impl(msg,line)
  character(*) :: msg
  integer :: line
  write(0,'(A,I0,A)') "TEST FAILED in test_mpl.F90 @ line ",line," :"
  write(0,*) msg
  stop 1
end subroutine

#define FAIL(msg) call fail_impl(msg,__LINE__)

program test_mpl
use ec_parkind, only : jpim
use mpl_module, only: mpl_init, mpl_end, mpl_rank, linitmpi_via_mpl

implicit none

integer(jpim) :: nprocs
logical :: verbose = .false.

call mpl_init(KPROCS=nprocs,ldinfo=verbose,ldenv=.true.)

if( nprocs == 0 )            FAIL("nprocs must be > 0")
if( mpl_rank == 0 )          FAIL("mpl_rank must be >= 1")
if( .not. linitmpi_via_mpl ) FAIL("linitmpi_via_mpl must be True")

call mpl_end(ldmeminfo=verbose)
! Note that with mpi_serial meminfo will not be printed regardless of ldmeminfo

end program
