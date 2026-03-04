! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

program drhook_ex2
use yomhook, only : jphook, dr_hook
use dr_hook_watch_mod, only : dr_hook_watch
implicit none
integer(4) :: n
real(8), allocatable :: a(:)
integer(4) :: j
character(len=256) a_out
character(len=20) cn
real(jphook) :: zhook_handle
call dr_hook('drhook_ex2',0,zhook_handle)
n = 5
allocate(a(n))
do j=1,n
  a(j) = j-1
enddo
! Watch & fail when A gets altered
call dr_hook_watch('a(:)',a,LDABORT=.TRUE.,LDPRINT=.FALSE.,LDTRBK=.FALSE.)
call sub1(a,n)
deallocate(a)
call dr_hook('drhook_ex2',1,zhook_handle)
end program drhook_ex2

subroutine sub1(a,n)
use yomhook, only : jphook, dr_hook
implicit none
integer(4), intent(in) :: n
real(8), intent(inout) :: a(n)
real(jphook) :: zhook_handle
call dr_hook('sub1',0,zhook_handle)
if (n > 0) call sub2(a(1))
call dr_hook('sub1',1,zhook_handle)
end subroutine sub1

subroutine sub2(s)
use yomhook, only : jphook, dr_hook
implicit none
real(8), intent(inout) :: s
real(jphook) :: zhook_handle
call dr_hook('sub2',0,zhook_handle)
!s = 1/s ! divide by zero removed
s = 1 ! a(1) has now been altered --> watch point should detect this
call dr_hook('sub2',1,zhook_handle)
end subroutine sub2
