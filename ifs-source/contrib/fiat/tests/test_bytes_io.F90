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
  write(0,'(A,I0,A)') "TEST FAILED in test_bytes_io.F90 @ line ",line," :"
  write(0,*) msg
  stop 1
end subroutine

#define FAIL(msg) call fail_impl(msg,__LINE__)

program bytes_io_test
use EC_PARKIND, only: jpim, jprd
use bytes_io_mod
use yomhook, only: lhook

implicit none

integer(jpim) :: unit, iret, ibuf, nbytes
integer(jpim) :: nbuffer(100)
integer(jpim), allocatable :: nrbuffer(:)
real(jprd) :: zbuffer(50)
real(jprd), allocatable :: zrbuffer(:)
integer(jpim) :: nrsize

integer(jpim), parameter :: sizeof_int  = 4
integer(jpim), parameter :: sizeof_real = 8

! Turn off DRHOOK, to avoid MPI init during testing
lhook = .False.

do ibuf=1,size(nbuffer)
  nbuffer(ibuf)=ibuf
enddo

call bytes_io_open( unit, "testfile_bytes_io", "w", iret )
if( iret < 0 ) FAIL("open for write failed")

call bytes_io_write( unit, size(nbuffer), sizeof_int, iret  )
if( iret < 0 ) FAIL("writing failed")

call bytes_io_write( unit, nbuffer, size(nbuffer)*sizeof_int, iret  )
if( iret < 0 ) FAIL("writing failed")

call bytes_io_close( unit, iret )
if( iret < 0 ) FAIL("close failed")

call bytes_io_open( unit, "testfile_bytes_io", "r", iret)
if( iret < 0 ) FAIL("open for read failed")

call bytes_io_read( unit, nrsize, sizeof_int, iret  )
if( iret < 0 ) FAIL("reading failed")

allocate( nrbuffer(nrsize) )
call bytes_io_read( unit, nrbuffer, nrsize*sizeof_int, iret  )
if( iret < 0 ) FAIL("reading failed")

call bytes_io_close( unit, iret )
if( iret < 0 ) FAIL("close failed")


do ibuf=1,size(nbuffer)
  if( nrbuffer(ibuf) /= nbuffer(ibuf) ) then
    FAIL("rbuffer read not equal to nbuffer written")
  endif
enddo


!=============================


do ibuf=1,size(zbuffer)
  zbuffer(ibuf)=ibuf
enddo

call bytes_io_open( unit, "testfile_bytes_io", "w", iret )
if( iret < 0 ) FAIL("open for write failed")

call bytes_io_write( unit, size(zbuffer), sizeof_int, iret  )
if( iret < 0 ) FAIL("writing failed")

call bytes_io_write( unit, zbuffer, size(zbuffer)*sizeof_real, iret  )
if( iret < 0 ) FAIL("writing failed")

call bytes_io_close( unit, iret )
if( iret < 0 ) FAIL("close failed")

call bytes_io_open( unit, "testfile_bytes_io", "r", iret)
if( iret < 0 ) FAIL("open for read failed")

call bytes_io_read( unit, nrsize, sizeof_int, iret  )
if( iret < 0 ) FAIL("reading failed")

if( nrsize /= size(zbuffer) ) FAIL("size does not match")

allocate( zrbuffer(nrsize) )
call bytes_io_read( unit, zrbuffer, nrsize*sizeof_real, iret  )
if( iret < 0 ) FAIL("reading failed")

call bytes_io_close( unit, iret )
if( iret < 0 ) FAIL("close failed")

do ibuf=1,size(zbuffer)
  if( zrbuffer(ibuf) /= zbuffer(ibuf) ) then
    FAIL("zrbuffer read not equal to zbuffer written")
  endif
enddo

write(0,'(A)') "SUCCESS"
end program
