! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

program drhook_ex1
use yomhook, only : jphook, dr_hook
implicit none
integer(4), parameter :: n = 100
real(8)    :: a(n)
integer(4) :: j
real(jphook) :: zhook_handle
call dr_hook('drhook_ex1',0,zhook_handle)
do j=1,n
  a(j) = j-1
enddo
call sub1(a,n)
call dr_hook('drhook_ex1',1,zhook_handle)
end program drhook_ex1

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
write(0,*) s
s = 1._8/s ! divide by zero
call dr_hook('sub2',1,zhook_handle)
end subroutine sub2
