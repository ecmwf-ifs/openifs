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

! Hans Hersbach ECMWF, March 2011
module mlink
  implicit none

  type link
    integer :: il(0:128-1),ic
  end type link 

  type tree
    integer                     :: nl,nlmax
    integer                     :: ns
    type (link), pointer        :: lnk(:)   => NULL()
  end type tree

contains

subroutine init_tree(t,nlmax)
  implicit none
  type (tree) :: t
  integer     :: nlmax

  t%nlmax=nlmax
  if(associated(t%lnk)) deallocate(t%lnk) ; allocate(t%lnk(nlmax))

  t%nl=1 ; t%lnk(1)%il(:)=0
  t%ns=0 ; t%lnk(1)%ic   =0
end subroutine init_tree
!-----------------------------------------------------------------------

subroutine find_index(t,sid,index)
  implicit none
  character(len=*) :: sid
  integer          :: index
  type (tree)      :: t

  integer             :: ichar,iascii,ilink,ilen,nl,ns
  type (link),pointer :: l(:)

  nl =  t%nl
  ns =  t%ns
  l  => t%lnk

! Find link by scanning through tree
  ilink=1
  ilen=len_trim(sid)
  do ichar=1,ilen
     iascii=iachar(sid(ichar:ichar))
     if (l(ilink)%il(iascii)==0) then
       nl=nl+1                   ! Add new link
       if (nl>t%nlmax) write(*,*)'Link Error, # of links too small:', t%nlmax
       l(nl   )%il(:     )=0
       l(nl   )%ic        =0
       l(ilink)%il(iascii)=nl    ! Update parent
     endif
     ilink=l(ilink)%il(iascii)
  enddo
  if (l(ilink)%ic==0) then       ! Add new element
     ns=ns+1
     l(ilink)%ic=ns
  endif

  index=l(ilink)%ic

  t%nl=nl
  t%ns=ns
end subroutine find_index
!-----------------------------------------------------------------------

subroutine print_tree_stats(t)
  implicit none
  type (tree) :: t

  write(*,'(A)')"Tree statistics:"
  write(*,'(A,I10)')"Number of entries allocated: ",t%ns
  write(*,'(A,2I10)')"Number of links used/allocated: ",t%nl,t%nlmax
  if(t%ns/=0) &
& write(*,'(A,I4)')"Average number of links per entry: ",t%nl/t%ns
end subroutine print_tree_stats

end module mlink
