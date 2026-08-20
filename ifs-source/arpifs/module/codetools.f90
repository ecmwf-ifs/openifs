! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

module codetools
   use parkind1,only: jpim,jprb

   implicit none

   interface setdefault
      module procedure :: setdefaulti,setdefaultr,setdefaultl
   end interface

   interface setptrdefault
      module procedure :: setptrdefaulti,setptrdefaultr,setptrdefaultl,setptrdefaultc
   end interface

   interface setptrdefault1
      module procedure :: setptrdefaulti1,setptrdefaultr1,setptrdefaultl1,setptrdefaultc1
   end interface

   interface setptrdefault2
      module procedure :: setptrdefaulti2,setptrdefaultr2,setptrdefaultl2,setptrdefaultc2
   end interface
contains
   elemental integer(kind=jpim) function ifelse(ltest,yes,no)
      logical,intent(in) :: ltest
      integer(kind=jpim),intent(in) :: yes,no

      if (ltest) then
         ifelse = yes
      else
         ifelse = no
      end if
   end function

   elemental real(kind=jprb) function ifelser(ltest,yes,no)
      logical,intent(in) :: ltest
      real(kind=jprb),intent(in) :: yes,no

      if (ltest) then
         ifelser = yes
      else
         ifelser = no
      end if
   end function

   elemental subroutine setdefaultc(s,val,x)
      character(len=*),intent(out) :: s
      character(len=*),intent(in) :: val
      character(len=*),optional,intent(in) :: x

      s = val
      if (present(x)) s = x
   end subroutine

   elemental integer(kind=jpim) function setdefaulti(val,x)
      integer(kind=jpim),intent(in) :: val
      integer(kind=jpim),optional,intent(in) :: x

      setdefaulti = val
      if (present(x)) setdefaulti = x
   end function

   elemental logical function setdefaultl(val,x)
      logical,intent(in) :: val
      logical,optional,intent(in) :: x

      setdefaultl = val
      if (present(x)) setdefaultl = x
   end function

   elemental real(kind=jprb) function setdefaultr(val,x)
      real(kind=jprb),intent(in) :: val
      real(kind=jprb),optional,intent(in) :: x

      setdefaultr = val
      if (present(x)) setdefaultr = x
   end function

   subroutine setptrdefaulti(xp,x,val)
      integer(kind=jpim),pointer,intent(inout) :: xp
      integer(kind=jpim),target,intent(inout) :: x
      integer(kind=jpim),intent(in) :: val

      xp => x
      xp = val
   end subroutine

   subroutine setptrdefaultl(xp,x,val)
      logical,pointer,intent(inout) :: xp
      logical,target,intent(inout) :: x
      logical,intent(in) :: val

      xp => x
      xp = val
   end subroutine

   subroutine setptrdefaultr(xp,x,val)
      real(kind=jprb),pointer,intent(inout) :: xp
      real(kind=jprb),target,intent(inout) :: x
      real(kind=jprb),intent(in) :: val

      xp => x
      xp = val
   end subroutine

   subroutine setptrdefaultc(xp,x,val)
      character(len=*),pointer,intent(inout) :: xp
      character(len=*),target,intent(inout) :: x
      character(len=*),intent(in) :: val

      xp => x
      xp = val
   end subroutine

   subroutine setptrdefaulti1(xp,x,val)
      integer(kind=jpim),pointer,intent(inout) :: xp(:)
      integer(kind=jpim),target,intent(inout) :: x(:)
      integer(kind=jpim),intent(in) :: val

      xp => x
      xp(:) = val
   end subroutine

   subroutine setptrdefaultl1(xp,x,val)
      logical,pointer,intent(inout) :: xp(:)
      logical,target,intent(inout) :: x(:)
      logical,intent(in) :: val

      xp => x
      xp(:) = val
   end subroutine

   subroutine setptrdefaultr1(xp,x,val)
      real(kind=jprb),pointer,intent(inout) :: xp(:)
      real(kind=jprb),target,intent(inout) :: x(:)
      real(kind=jprb),intent(in) :: val

      xp => x
      xp(:) = val
   end subroutine

   subroutine setptrdefaultc1(xp,x,val)
      character(len=*),pointer,intent(inout) :: xp(:)
      character(len=*),target,intent(inout) :: x(:)
      character(len=*),intent(in) :: val

      xp => x
      xp(:) = val
   end subroutine

   subroutine setptrdefaulti2(xp,x,val)
      integer(kind=jpim),pointer,intent(inout) :: xp(:,:)
      integer(kind=jpim),target,intent(inout) :: x(:,:)
      integer(kind=jpim),intent(in) :: val

      xp => x
      xp(:,:) = val
   end subroutine

   subroutine setptrdefaultl2(xp,x,val)
      logical,pointer,intent(inout) :: xp(:,:)
      logical,target,intent(inout) :: x(:,:)
      logical,intent(in) :: val

      xp => x
      xp(:,:) = val
   end subroutine

   subroutine setptrdefaultr2(xp,x,val)
      real(kind=jprb),pointer,intent(inout) :: xp(:,:)
      real(kind=jprb),target,intent(inout) :: x(:,:)
      real(kind=jprb),intent(in) :: val

      xp => x
      xp(:,:) = val
   end subroutine

   subroutine setptrdefaultc2(xp,x,val)
      character(len=*),pointer,intent(inout) :: xp(:,:)
      character(len=*),target,intent(inout) :: x(:,:)
      character(len=*),intent(in) :: val

      xp => x
      xp(:,:) = val
   end subroutine

   integer(kind=jpim) function addindex(i,n)
      integer(kind=jpim),intent(in) :: n
      integer(kind=jpim),intent(inout) :: i

      i = i+n
      addindex = i
   end function

   integer(kind=jpim) function indexadd(i,n)
      integer(kind=jpim),intent(in) :: n
      integer(kind=jpim),intent(inout) :: i

      indexadd = i+1
      i = i+n
   end function

   real(kind=jprb) function shamean(z,n)
      real(kind=jprb),intent(in) :: z(*)
      integer(kind=jpim),intent(in) :: n

      integer(kind=jpim) :: i,ni,off

      ni = sqrt(real(n,jprb))
      ! assert that last off < n
      if (ni*(ni+1) >= n) ni = ni-1

      shamean = 0
      off = 0
      if (ni > 1) then
         do i=1,ni
            shamean = (shamean+(2*mod(i,2)-1)*sum(z(off+1:off+i))/i)/2
            off = off+i
         end do

         do i=ni,1,-1
            shamean = (shamean-(2*mod(i,2)-1)*sum(z(off+1:off+i))/i)/2
            off = off+i
         end do
      end if

      shamean = (shamean+sum(z(off+1:n))/(n-off))/2
   end function

   subroutine stat2d(out,z,np,nf,nb)
      real(kind=jprb) :: z(:,:,:)
      integer(kind=jpim),intent(in) :: out,np,nf,nb

      integer(kind=jpim) :: i,nz
      real(kind=jprb) :: m

      if (nf == 0) then
         print*,"no fields"
         return
      end if

      nz = np*nb

      do i=1,nf
         m = sum(z(:,i,:))/nz
         write(out,"('#',i0,':',4(x,g0))") i,minval(z(:,i,:)),maxval(z(:,i,:)),&
            shamean(z(:,i,:),nz),sum(z(:,i,:)**2)/nz-m**2
      end do
   end subroutine

   subroutine stat3d(out,z,np,nl,nf,nb)
      real(kind=jprb) :: z(:,:,:,:)
      integer(kind=jpim),intent(in) :: out,np,nl,nf,nb

      integer(kind=jpim) :: i,nz
      real(kind=jprb) :: m

      if (nf == 0) then
         print*,"no fields"
         return
      end if

      nz = np*nl*nb

      do i=1,nf
         m = sum(z(:,:,i,:))/nz
         write(out,"('#',i0,':',4(x,g0))") i,minval(z(:,:,i,:)),maxval(z(:,:,i,:)),&
            shamean(z(:,:,i,:),nz),sum(z(:,:,i,:)**2)/nz-m**2
      end do
   end subroutine

   character(len=90) function statb2d(z,np)
      real(kind=jprb) :: z(:)
      integer(kind=jpim),intent(in) :: np

      real(kind=jprb) :: m

      m = sum(z(1:np))/np
      write(statb2d,"(4(x,g0))") minval(z(1:np)),maxval(z(1:np)),shamean(z(1:np),np),&
         sum(z(1:np)**2)/np-m**2
   end function

   character(len=90) function statb3d(z,np,nl)
      real(kind=jprb) :: z(:,:)
      integer(kind=jpim),intent(in) :: np,nl

      integer(kind=jpim) :: nz
      real(kind=jprb) :: m

      nz = np*nl

      m = sum(z(1:np,1:nl))/nz
      write(statb3d,"(4(x,g0))") minval(z(1:np,1:nl)),maxval(z(1:np,1:nl)),&
         shamean(z(1:np,1:nl),nz),sum(z(1:np,1:nl)**2)/nz-m**2
   end function
end module
