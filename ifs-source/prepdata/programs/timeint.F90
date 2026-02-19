! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

!module parkind1
!integer,parameter :: jpim=4
!integer,parameter :: jprb=8
!end module parkind1

program timeint

!     Input:  -mon_in  - input experiment  monthly climate fields
!     Output: _mon_out - output experiment geographical fields
use parkind1  ,only : jpim, jprb
use grib_api
use ec_datetime_mod, only : hourdiff

implicit none

integer(kind=jpim) :: iyyyymdh, iyyyy, imm, idd, ihh, ihdiff2, ihdiff, idate
integer(kind=jpim) :: iarg, ioptval, icodemo, ierr, imonin, iret, itable
integer(kind=jpim) :: ifoundfi, icode, iyo, imo, imt1, imt2, iyt1, iyt2
integer(kind=jpim) :: imonout, icount, igribin, igribout, ibitmap, ilev
integer(kind=jpim) :: nvalues, nlevs, nv
		     
real(kind=jprb) :: zmiss=-9999._jprb
real(kind=jprb)    :: zwei1, zwei2, zcodemo
real(kind=jprb),allocatable :: clim1(:,:),clim2(:,:),expt(:)

logical, save :: lfirst = .true.
      
integer(kind=jpim) :: getopt
      
character*256 clarg
character clopt
character*256 yinpfi, youtfi
character*10  ydate

!          1.        Reads data

iarg=0
get_loop:do
  ioptval=GETOPT('d:c:i:m:o:?',clarg)
  clopt=char(ioptval)
  if (ioptval.eq.-1) then
    if (iarg.eq.0) then
      print *,'NO OPTION SUPPLIED ...'
      print *,'USAGE: timeint -d basetime (YYYYMMDDHH) ', &
          & '-c grib_code -i input_mon -o output'
      print *,'more help with "timeint -\?"'
      stop
    else
      exit get_loop
    endif
  endif
  if (clopt.eq.'d') then
    iarg=iarg+1
    read(unit=clarg,fmt='(i10)') iyyyymdh
  else if (clopt.eq.'c') then
    iarg=iarg+1
    read(unit=clarg,fmt='(f7.0)') zcodemo
  else if (clopt.eq.'i') then
    iarg=iarg+1
    read(unit=clarg,fmt='(a256)') yinpfi
  else if (clopt.eq.'o') then
    iarg=iarg+1
    read(unit=clarg,fmt='(a256)') youtfi
  else if (clopt.eq.'?') then
    print *,'timeint interpolates a monthly climate to a given date '
    print *,'basetime   == current date/time (YYYYMMDDHH)'
    print *,'grib_code  == GRIB code number'
    print *,'input_mon  == monthly climate input file'
    print *,'output     == time interpolated output file '
    print *,'USAGE: timeint -d basetime -c grib_code -i input_mon -o output'
    stop
  endif
enddo get_loop

!          2.        Read monthly climate file
!                    -------------------------

!            2.1     Open monthly climate file

ierr=0

call grib_open_file(imonin,YINPFI,'r',iret)
if (iret /= GRIB_SUCCESS) then
   write(ierr,'(''Error opening file _mon_in, iret = '',i10)') iret
   call abort
   stop
endif
 
!            2.2     Read first record

call grib_new_from_file(imonin,igribin,iret)
if (iret /= GRIB_SUCCESS) then
   write(ierr,'(''Error reading monthly climate file, iret = '',i10)') iret
   call abort
   stop
endif

call grib_clone(igribin,igribout)

!            2.3     Unpack header records to determine data lengths


call grib_get(igribin,'bitmapPresent',ibitmap)
if (ibitmap == 1) then
  call grib_set(igribin,'missingValue',zmiss)
endif
call grib_get(igribin,'getNumberOfValues',nvalues)
call grib_get(igribin,'NV',nv)
if (nv > 0) then
  nlevs = nv/2 - 1
else
  nlevs = 1
endif

!            2.4     Allocate arrays dynamically

allocate(clim1(nvalues,nlevs))
allocate(clim2(nvalues,nlevs))
allocate(expt(nvalues))

write(ydate,'(i10)') iyyyymdh
read(ydate,'(i4,3i2)') iyyyy,imm,idd,ihh

if (idd.ge.15) then
  imt1=imm
  imt2=1+mod(imm,12)
  iyt1=iyyyy
  if(imt2.eq.1) then
    iyt2=iyt1+1
  else
    iyt2=iyt1
  endif
else
  imt1=1+mod(imm+10,12)
  imt2=imm
  if(imt1.eq.12) then
    iyt1=iyyyy-1
  else
    iyt1=iyyyy
  endif
  iyt2=iyyyy
endif

ifoundfi=0
 
!            2.5     Loop over climate file, reading GRIB records

scan:do
  if ( .not.lfirst ) then
    call grib_new_from_file(imonin,igribin,iret)
    if (iret /= GRIB_SUCCESS) then
      if (iret == GRIB_END_OF_FILE) then
        if (ifoundfi == 2*nlevs) exit scan
        write(ierr,'(''EOF detected on monthly climate file'')')
        write(ierr,'(''Number of fields located = '',i10)') ifoundfi
        write(ierr,'(''Number of fields required = '',i10)') 2*nlevs
        call abort
        stop
      else
        write(ierr,'(''Error reading monthly climate file'')')
        write(ierr,'(''iret = '',i10)') iret
        call abort
        stop
      endif
    endif
  endif
  lfirst=.false.

  call grib_get(igribin,'paramId',icode)

  ! Handle possible parameter.table notation
  icodemo = floor(zcodemo)
  itable = nint(1000*(zcodemo-icodemo))
  if (itable /= 128) icodemo = 1000*itable + icodemo

  if (icodemo.eq.icode) then

!            2.6     Decode GRIB

    call grib_get(igribin,'bitmapPresent',ibitmap)
    if (ibitmap == 1) call grib_set(igribin,'missingValue',zmiss)

    call grib_get(igribin,'dataDate',idate)
    write(ydate,'(i8)') idate
    read(ydate,'(i4,i2)') iyo,imo

    if (imo == imt1 .or. imo == imt2) then 
      ifoundfi=ifoundfi+1

      if (nv > 0) then
        call grib_get(igribin,'level',ilev)
      else
        ilev = 1
      endif
      
      if (imo == imt1) call grib_get(igribin,'values',clim1(:,ilev),iret)
      if (imo == imt2) call grib_get(igribin,'values',clim2(:,ilev),iret)
      if (iret /= GRIB_SUCCESS) then
        write(ierr,'(''Error decoding monthly climate GRIB record'')')
        write(ierr,'(''iret = '',i10)') iret
        call abort
        stop
      endif
    endif
 
  endif
   
enddo scan

!=======================================================================

!          3.        Recode and write
!                    ----------------

!            3.1     Open output file

call grib_open_file(imonout,YOUTFI,'w',iret)
if (iret /= GRIB_SUCCESS) then
   write(ierr,'(''Error opening file _mon_out, iret = '',i10)') iret
   call abort
   stop
endif
icount=0

!            3.2       Interpolate in time to date iyyyymdh

call hourdiff(iyt2,imt2,15,0,iyyyy,imm,idd,ihh,ihdiff2,iret)
call hourdiff(iyt2,imt2,15,0,iyt1,imt1,15,0,ihdiff,iret)

zwei1=float(ihdiff2)/float(ihdiff)
zwei2=1.-zwei1

do ilev = 1, nlevs
  if(ibitmap == 0) then
    expt(:)=zwei1*clim1(:,ilev)+zwei2*clim2(:,ilev)
  else
    where (clim1(:,ilev) == zmiss .or. clim2(:,ilev) == zmiss)
        expt(:)=zmiss
    elsewhere
      expt(:)=zwei1*clim1(:,ilev)+zwei2*clim2(:,ilev)
    end where
  endif

  !            3.3     Recode GRIB

  call grib_set(igribout,'paramId',icodemo)

  call grib_set(igribout,'dataDate',iyyyy*10000+imm*100+idd)
  call grib_set(igribout,'dataTime',ihh*100)

  if (nv > 0) call grib_set(igribout,'level',ilev)

  call grib_set(igribout,'values',expt)

  !            3.4     Write modified GRIB record to file _mon_out

  call grib_write(igribout,imonout,iret)
  if (iret /= GRIB_SUCCESS) then
    write(ierr,'(''Error writing GRIB record'')')
    write(ierr,'(''iret = '',i10)') iret
    write(ierr,'(''icode = '',i10)') icode
    call abort
    stop
  endif

end do

write(*,*) ' Code ',icodemo,' Date/time= ',iyyyymdh
write(*,"(' Interpolated between month ',i2,' with weight ',f5.3, &
     & ' and month ',i2,' with weight ',f5.3)")  &
     &  imt1,zwei1,imt2,zwei2


!=======================================================================

!          4.        Termination
!                    -----------

!            4.1     Close files

call grib_close_file(imonout,iret)
if (iret /= GRIB_SUCCESS) then
   write(ierr,'(''Error closing file _mon_out, iret = '',i10)') iret
   call abort
   stop
endif

call grib_release(igribin)
call grib_close_file(imonin,iret)
if (iret /= GRIB_SUCCESS) then
   write(ierr,'(''Error closing file _mon_in, iret = '',i10)') iret
   call abort
   stop
endif

stop
end program timeint
