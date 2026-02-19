! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

program drhook_ex3
use yomhook, only : jphook, dr_hook
implicit none
integer(4) :: n
real(8), allocatable  :: a(:)
integer(4) :: j
real(jphook) :: zhook_handle
call dr_hook('drhook_ex3',0,zhook_handle)
n = 100
allocate(a(n))
do j=1,n
  a(j) = j-1
enddo
write(0,*)'drhook_ex3: sum#1 = ',sum(a)
deallocate(a)

n = 5 * n
allocate(a(n))
do j=1,n
  a(j) = j + n
enddo
do j=1,n
 call sub1(a,j)
enddo
write(0,*)'drhook_ex3: sum#2 = ',sum(a)
deallocate(a)

n = n/10
allocate(a(-n:n))
do j=-n,n
  a(j) = j + 2*n
enddo
do j=1,n
 call sub1(a,j)
enddo
write(0,*)'drhook_ex3: sum#3 = ',sum(a)
deallocate(a)
call dr_hook('drhook_ex3',1,zhook_handle)
end program drhook_ex3

subroutine sub1(a,n)
use yomhook, only : jphook, dr_hook
implicit none
integer(4), intent(in) :: n
real(8), intent(inout) :: a(n)
integer(4) j
real(jphook) :: zhook_handle
call dr_hook('sub1',0,zhook_handle)
do j=1,n
  if (mod(j,2) == 0) call sub2(a(j))
  a(j) = 2*a(j) + 1
enddo
call dr_hook('sub1',1,zhook_handle)
end subroutine sub1

subroutine sub2(s)
use yomhook, only : jphook, dr_hook
implicit none
real(8), intent(inout) :: s
real(jphook) :: zhook_handle
call dr_hook('sub2',0,zhook_handle)
s = 1/(s+1)
call dr_hook('sub2',1,zhook_handle)
end subroutine sub2
