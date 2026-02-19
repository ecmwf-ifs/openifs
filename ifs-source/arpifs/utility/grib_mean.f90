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

program grib_mean
! Average input grib file along supplied keys
! Hans Hersbach July 2013, ECMWF

  use grib_api
  use mlink

  implicit none
  character(len=20)  :: name_space
  integer            :: kiter,ifile,ofile,sfile,vfile,igrib,iret
  character(len=256) :: key,cunit
  character(len=256) :: value
  character(len=512) :: ckey,tkey,cset
  character          :: s,c=','
  logical            :: lnot_ok,lstdv,lvari,lmean,l2,lfilter,lscale,lstrict=.true.
  
  integer,parameter  :: ninmax=10000
  character(len=256) :: gribfile(ninmax),meanfile,stdvfile,varifile,keys,filter,carg
  integer            :: ioptval,getopt,nin,kin
  character          :: options*22,copt

  integer,parameter  :: nmaxmessage=1000000 ! one-million !
  integer            :: i,nvals

  type (tree)        :: thead
  integer            :: nhead,ihead,ncheck,nmin

  character(len=1), allocatable :: message(:)
  real*8          , allocatable :: vals(:)
  integer(kind=kindOfSize)  :: mlen

  type grib_message
    integer             :: igrib,nvals,nm
    real*8              :: xmod
    integer,allocatable :: n(:)
    real*8, allocatable :: x(:),v(:)
  end type grib_message

  logical           :: filter_out
  real*8            :: hmod,rscale
  real*8, parameter :: rmiss=-1.38E31

  type(grib_message), pointer :: gm(:) => NULL(), gg => NULL()
  allocate(gm(nmaxmessage))

 !-Crack options
  nin=0       ; ncheck=0  ; nmin=1 ; rscale=1.0
  meanfile="" ; stdvfile="" ; varifile=""; keys="" ; cset=""
  data options/'d:i:s:k:o:n:v:m:f:r:N:'/
  do
    ioptval=getopt(options,carg)
    carg=trim(carg)
    copt=char(ioptval)

    if (ioptval <=0) exit
    if (copt == 'i') nin=nin+1
    if (copt == 'i') gribfile(nin) = trim(carg)
    if (copt == 'k') keys          = c//trim(carg)//c
    if (copt == 'f') filter        = c//trim(carg)//c
    if (copt == 'o') meanfile      = trim(carg)
    if (copt == 'd') stdvfile      = trim(carg)
    if (copt == 'v') varifile      = trim(carg)
    if (copt == 's') cset          = trim(carg)
    if (copt == 'n') read(carg,*)    ncheck
    if (copt == 'n') lstrict=.true.
    if (copt == 'N') read(carg,*)    ncheck
    if (copt == 'N') lstrict=.false.
    if (copt == 'm') read(carg,*)    nmin      ; nmin=max(1,nmin)
    if (copt == 'r') read(carg,*)    rscale
  enddo
  lmean=(meanfile/="")
  lstdv=(stdvfile/="")
  lvari=(varifile/="")
  l2=lstdv.or.lvari
  lscale=(rscale/=1.0)

!-Some QC
  lnot_ok=.false.
  if (nin==0)                                             lnot_ok=.true.
  if (meanfile=="" .and. stdvfile=="" .and. varifile=="") lnot_ok=.true.
  if (lnot_ok) then
     write(*,*)'Usage: grib_mean -i gribfile1 -i gribfile2,..&
     & { -o mean_gribfile -d stdv_gribfile -v vari_gribfile } -k mkey1,mkey2,..&
     & [ -s key1=value1,key2=value2,.. -n n_fields_expected_per_mean -m minimal_in_mean - r rescale_factor ]'
     stop
  endif

  write(*,'(2A)')"Keys to be averaged: ",trim(keys)

!-Open input GRIB file and loop over all messages
  call init_tree(thead,500000) ; nhead=0 ; ihead=0 ; nvals=0
  do kin=1,nin
   call grib_open_file(ifile,trim(gribfile(kin)),'r') 
   do
     call grib_new_from_file(ifile,igrib, iret) ; if (iret==GRIB_END_OF_FILE) exit

   !-Find all mars keys.
     name_space='mars' ; ckey="" ; tkey="" ; s="" ; lfilter=.false.
     call grib_keys_iterator_new(igrib,kiter,name_space)
     do
        call grib_keys_iterator_next(kiter, iret) ; if (iret /= GRIB_SUCCESS) exit
        call grib_keys_iterator_get_name(kiter,key)
        call grib_get_string(igrib,trim(key),value)

      !-Check filter
        if ( filter_out(filter,key,value) ) lfilter=.true.
      !-Filter out the averaging keys in the search string
        if ( index(keys,     trim(key)//c)/=0 ) cycle

        call set_qmarks(keys,key,value)

        tkey=trim(tkey)//"_"//trim(value)
        ckey=trim(ckey)//trim(s)//trim(key)//'='//trim(value) ; s=','
     end do

   !-Filter out?
     if (lfilter) cycle
    
   !-Find header index and point gg to it
     call find_index(thead,trim(tkey),ihead)    
     gg => gm(ihead)

   !-A new set of headers has been found
     if (ihead> nhead) then
         write(*,'("Mean ",I6.6,": ",A)') ihead,trim(ckey)

         call grib_get_message_size(igrib, mlen)
         allocate(message(mlen), stat=iret)
         call grib_copy_message(igrib,message)
         call grib_new_from_message(gg%igrib,message)
         deallocate(message)

       !-Reserve space for output grib field, and initialize
         call grib_get_size(igrib,'values',gg%nvals)
         if (l2) allocate(gg%v(gg%nvals))
         if (l2) gg%v(:)=0.0
                 allocate(gg%n(gg%nvals),gg%x(gg%nvals))
                 gg%nm  =0
                 gg%n(:)=0
                 gg%x(:)=0.0 

       !-Modulo 360?
         gg%xmod=0.0
         call grib_get_string(igrib,'units',cunit)
         if (trim(cunit)=="degrees") gg%xmod=360.0  ! for some ocean wave fields
     endif
     nhead=max(ihead,nhead)

   !-Update statistics
     if (nvals/=gg%nvals) then
         nvals= gg%nvals
         if (allocated(vals)) deallocate(vals) ; allocate(vals(nvals))
     endif
     call grib_set(igrib,'missingValue',rmiss)
     call grib_get(igrib,'values'      ,vals )
     do i=1,nvals
        if (vals(i)==rmiss) cycle
        if(lscale)vals(i)=rscale*vals(i)
        vals(i)=hmod(vals(i)-gg%x(i),gg%xmod,0.5_8)
        gg%n(i)=gg%n(i)+1
        gg%x(i)=gg%x(i)+vals(i)/gg%n(i)
        if (l2)&
      & gg%v(i)=(gg%v(i)+vals(i)*vals(i)/gg%n(i))*(gg%n(i)-1.)/gg%n(i)
     enddo
     gg%nm=gg%nm+1
     
   !-Move on to next GRIB message
     call grib_keys_iterator_delete(kiter)
     call grib_release(igrib)
   enddo
 !-Move on to next input file
   call grib_close_file(ifile)
  enddo 

!-Open output GRIB file(s), post process messages and write out
  if (lmean) call grib_open_file(ofile,trim(meanfile),'w')
  if (lstdv) call grib_open_file(sfile,trim(stdvfile),'w')
  if (lvari) call grib_open_file(vfile,trim(varifile),'w')
  do ihead=1,nhead
      gg => gm(ihead)
      if (ncheck>0 .and. gg%nm/=ncheck) then
         write(*,'("Mean ",I6.6,": Incorrect # of fields. Found/Imposed:",2I4)')ihead,gg%nm,ncheck
         if (lstrict) call abort()
         cycle
      endif
      if(minval(gg%n)<nmin) &
    & call grib_set  (gg%igrib,'bitmapPresent',1    )
      call set_keys  (gg%igrib,trim(cset))
      call grib_set  (gg%igrib,'missingValue' ,rmiss)
      do i=1,gg%nvals
         if (lmean.and.gg%n(i)<nmin)  gg%x(i) = rmiss
         if (l2   .and.gg%n(i)<nmin)  gg%v(i) = rmiss
      enddo
      if (lmean) then
         do i=1,gg%nvals
            if (gg%n(i)>=nmin) gg%x(i) = hmod(gg%x(i),gg%xmod,0.0_8)
         enddo
         call grib_set  (gg%igrib,'values', gg%x)
         call grib_write(gg%igrib         ,ofile)
      endif
      if (lvari) then
         call grib_set  (gg%igrib,'values', gg%v)
         call grib_write(gg%igrib         ,vfile)
      endif
      if (lstdv) then
         do i=1,gg%nvals
            if (gg%n(i)>=nmin) gg%v(i) = sqrt(max(0.,gg%v(i)))
         enddo
         call grib_set  (gg%igrib,'values', gg%v)
         call grib_write(gg%igrib         ,sfile)
      endif
  enddo
  if (lmean) call grib_close_file(ofile)
  if (lstdv) call grib_close_file(sfile)
  if (lvari) call grib_close_file(vfile)

!-Print out some statistics:
! call print_tree_stats(thead)
  
end program grib_mean

!---------------------------------------------------------
function hmod(x,xm,xf)
  implicit none
  real*8 :: hmod
  real*8, intent(in)  :: x,xm,xf
  hmod=x ; if(xm>0.0) hmod=mod(x+(2.+xf)*xm,xm)-xf*xm
end

!---------------------------------------------------------

subroutine set_qmarks(keys,key,value)
  implicit none
  character(len=*), intent(in)  :: keys,key
  character(len=*), intent(out) :: value

  integer   :: i,j

  i=index(keys,','//trim(key)//'=') ; if(i==0) return
  i=i+index(keys(i:),'=')
  do j=1,len_trim(value)
     if (keys(i:i)=='?') value(j:j)='?'
     i=i+1
  enddo
end subroutine set_qmarks

!---------------------------------------------------------

function filter_out(filter,key,value)
  implicit none
  character(len=*), intent(in) :: filter,key,value

  logical          :: filter_out
  character :: c
  integer   :: i,j

  filter_out=.false.
  i=index(filter,','//trim(key)//'=') ; if(i==0) return
  i=i+index(filter(i:),'=')
  do j=1,len_trim(value)
     c=filter(i:i) ; i=i+1
     if (c=='?') cycle
     if (c/=value(j:j)) filter_out=.true.
  enddo
end function filter_out

!---------------------------------------------------------

subroutine set_keys(igrib,keys)
  use grib_api
  implicit none
  integer, intent(in) :: igrib
  character(len=*), intent(in) :: keys
  integer          :: ikey,ierr

  integer   :: i1,i2,i3,l
  character :: c=',',s='='

  l=len_trim(keys) ; if(l==0) return
  i1=0
  do
     i2=i1+index(keys(i1+1:l)//s,s)
     i3=i1+index(keys(i1+1:l)//c,c)

   !-Some keys can only be set as an integer; i.e. string is not allowed
!    read(keys(i2+1:i3-1),'(i)',iostat=ierr)ikey
     read(keys(i2+1:i3-1),*,iostat=ierr)ikey

     if (ierr==0) then 
        call grib_set_int   (igrib,trim(keys(i1+1:i2-1)),ikey)
     else
        call grib_set_string(igrib,trim(keys(i1+1:i2-1)),trim(keys(i2+1:i3-1)))
     endif
     i1=i3 
     if (i1==l+1) exit
  enddo
end subroutine set_keys
!---------------------------------------------------------
