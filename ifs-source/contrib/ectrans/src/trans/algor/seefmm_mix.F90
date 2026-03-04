! (C) Copyright 2009- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

module seefmm_mix
!**** *SEEFMM_MIX*  - Implementation of Simple Exponential Expansion FMM

!     Purpose.
!     --------
!     Implementation of Simple Exponential Expansion FMM

!**   Interface.
!     ----------

!     Method.
!     -------
!     Based on Algorithm described in Section 4 of the article 
!     "An improved fast multipole algorithm for potential fields on the line "


!     Reference.
!     ----------
!     "An improved fast multipole algorithm for potential fields on the line "
!     by Norman Yarvin and Vladimir Rohklin, SIAM J. Numer. Anal. Vol. 36,No. 2,629-666.  [1]
!
!     Author.
!     -------
!     Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!     Original : 2009-06-04
!     ------------------------------------------------------------------


use parkind1,only : jpim     ,jprb, jprd
use ecsort_mix, only : keysort
use wts500_mod, only: wts500

private

integer(kind=jpim) :: nfmm_lim=200 ! Appr. break-even limit for FMM
integer(kind=jpim),parameter :: nquadEm14=28 ! Quadrature size for eps~=1.e-14
integer(kind=jpim),parameter :: nquadEm10=20!  Quadrature size for eps~=1.e-10
integer(kind=jpim),parameter :: nquadEm07=14!  Quadrature size for eps~=1.e-07

type fmm_type
integer(kind=jpim)          :: nxy    ! Total number of point "nx+ny"
integer(kind=jpim)          :: nx     ! Number of 'x' points
integer(kind=jpim)          :: nquad  ! Quadrature N
integer(kind=jpim)          :: ncik   ! Number of elem. in cik
real(kind=jprb)             :: rw(56) ! Quadrature weights
real(kind=jprb)   , pointer :: rdexp(:,:)  ! exp(xy(i)-xy(i-1))
integer(kind=jpim), pointer :: index(:)    ! index for sorted xy
integer(kind=jpim), pointer :: nclose(:)   ! No of "close" points
real(kind=jprb)   , pointer :: cik(:)      ! Correction term (142 in [1])

end type fmm_type

public :: fmm_type, setup_seefmm, free_seefmm, seefmm_mulm

contains
recursive subroutine setup_seefmm(kx,px,ky,py,ydfmm,pdiff)

implicit none

!**** *SETUP_SEEFMM*  - Setup seefmm

! Purpose - Pre-computations for applying SEEFMM

! Explicit arguments :
! --------------------
! kx - Number of x points
! px - x points
! ky - Number of y points 
! py - y points
! ydfmm - result of pre-computations
! pdiff - difference matrix (optional)

integer(kind=jpim),intent(in)  :: kx    
real(kind=jprd)   ,intent(in)  :: px(:)
integer(kind=jpim),intent(in)  :: ky
real(kind=jprb)   ,intent(in)  :: py(:)
type(fmm_type)    ,intent(out) :: ydfmm
real(kind=jprb),optional,intent(in)  :: pdiff(:,:)

real(kind=jprb) :: zxy(kx+ky),zrt(56),zcik((kx+ky)*(kx+ky))
real(kind=jprb) :: zr
integer(kind=jpim) :: ixy
!---------------------------------------------------------------------------
ydfmm%nx=kx
ixy=kx+ky
ydfmm%nxy=ixy
allocate(ydfmm%index(ixy))
!ydfmm%nquad=nquadEm14 !Set precicion to 1.E-14
ydfmm%nquad=nquadEm07 !Set precicion to 1.E-07
! Combine px and py to form xxy, compute ascending index for xxy
call comb_xy(kx,px,ky,py,ixy,zxy,ydfmm%index)
! Setup quadrature, scale (see 3.1.1 in [1])
call suquad(ixy,zxy(ydfmm%index(1))-zxy(ydfmm%index(ixy)),&
 & ydfmm%nquad,ydfmm%rw,zrt,zr)
allocate(ydfmm%rdexp(ydfmm%nquad,ixy))
allocate(ydfmm%nclose(ixy))
! Main pre-computation
call prepotf(kx,ixy,ydfmm%nquad,ydfmm%rw,zrt,zr,zxy,ydfmm%index,&
 & ydfmm%rdexp,ydfmm%nclose,zcik,ydfmm%ncik,pdiff)
! Needed as size of cik unknown beforehand
allocate(ydfmm%cik(ydfmm%ncik))
ydfmm%cik(:)=zcik(1:ydfmm%ncik)

end subroutine setup_seefmm
!==========================================================================

subroutine free_seefmm(ydfmm)
implicit none

!**** *FREE_SEEFMM*  - Release memory

! Purpose - Release memory used by ydfmm

! Explicit arguments :
! --------------------
! ydfmm - result of pre-computations
type(fmm_type)    ,intent(inout) :: ydfmm

deallocate(ydfmm%index)
deallocate(ydfmm%rdexp)
deallocate(ydfmm%nclose)
deallocate(ydfmm%cik)

end subroutine free_seefmm

!==========================================================================
recursive subroutine potf(kn,kx,ldxout,kquad,prw,pq,prdexp,kindex,kclosel,kcik,pcik,ptheta)
implicit none

integer(kind=jpim),intent(in)  :: kn
integer(kind=jpim),intent(in)  :: kx
logical           ,intent(in)  :: ldxout
integer(kind=jpim),intent(in)  :: kquad
real(kind=jprb)   ,intent(in)  :: prw(:)
real(kind=jprb)   ,intent(in)  :: pq(:)
real(kind=jprb)   ,intent(in)  :: prdexp(:,:)
integer(kind=jpim),intent(in)  :: kindex(:)
integer(kind=jpim),intent(in)  :: kclosel(:)
integer(kind=jpim),intent(in)  :: kcik
real(kind=jprb)   ,intent(in)  :: pcik(:)
real(kind=jprb)   ,intent(out) :: ptheta(:)

real(kind=jprb)    :: zalpha(kquad),zq(kn),ztheta(kn)
integer(kind=jpim) :: j1,j2,jm,inumc,idist,iquad
integer(kind=jpim) :: iout,iq,i1,i1p1,i1pd,ik1,ix,iy
logical :: lxy,llxy(kn)
 
lxy(ik1) = (ik1 <= kx .eqv. ldxout)
!-------------------------------------------------------------------------

ztheta(:)=0.0_JPRB
if(ldxout) then
  ix=0
  iy=-kx
else
  ix=-kx
  iy=0
endif

do j1=1,kn
  i1=kindex(j1)
  llxy(j1)=lxy(i1)
  if(llxy(j1)) then
    zq(j1)=pq(kindex(j1)+ix)
  else
    zq(j1)=0.0_jprb
  endif
enddo

zalpha(:)=zq(1)
do j1=2,kn
  if(llxy(j1)) then
    do jm=1,kquad
      zalpha(jm)=zalpha(jm)*prdexp(jm,j1)+zq(j1)
    enddo
  else
    do jm=1,kquad
      zalpha(jm)=zalpha(jm)*prdexp(jm,j1)
      ztheta(j1)=ztheta(j1)+prw(jm)*zalpha(jm)
    enddo
  endif
enddo

zalpha(1:kquad)=zq(kn)
do j1=kn-1,1,-1
  if(llxy(j1)) then
    do jm=1,kquad
      zalpha(jm)=zalpha(jm)*prdexp(jm,j1+1)+zq(j1)
    enddo
  else
    do jm=1,kquad
      zalpha(jm)=zalpha(jm)*prdexp(jm,j1+1)
      ztheta(j1)=ztheta(j1)-prw(jm)*zalpha(jm)
    enddo
  endif
enddo


IF(kcik > 0) then
  inumc=0
  do j1=1,kn-1
    do j2=1,kclosel(j1)
      idist=j2
      if(.not.llxy(j1) .and. llxy(j1+idist)) then
        inumc=inumc+1
        ztheta(j1)=ztheta(j1)-pcik(inumc)*zq(j1+idist)
      elseif(llxy(j1) .and. .not.llxy(j1+idist)) then
        inumc=inumc+1
        ztheta(j1+idist)=ztheta(j1+idist)+pcik(inumc)*zq(j1)
      endif
    enddo
  enddo
endif

do j1=1,kn
  if(.not. llxy(j1)) then
    i1=kindex(j1)
    ptheta(i1+iy)=ztheta(j1)
  endif
enddo

end subroutine potf
!==========================================================================
recursive subroutine seefmm_mulv(ydfmm,ldxout,pq,ptheta)
implicit none

type(fmm_type)    ,intent(in)  :: ydfmm
logical           ,intent(in)  :: ldxout
real(kind=jprb)   ,intent(in)  :: pq(:)
real(kind=jprb)   ,intent(out) :: ptheta(:)

!-------------------------------------------------------------------------
call potf(ydfmm%nxy,ydfmm%nx,ldxout,ydfmm%nquad,&
 & ydfmm%rw,pq,ydfmm%rdexp,ydfmm%index,&
 & ydfmm%nclose,ydfmm%ncik,ydfmm%cik,ptheta)

end subroutine seefmm_mulv
!==========================================================================
recursive subroutine seefmm_mulm(ydfmm,km,kskip,ldxout,pq,ptheta)
implicit none

type(fmm_type)    ,intent(in)  :: ydfmm
integer(kind=jpim),intent(in)  :: km
integer(kind=jpim),intent(in)  :: kskip
logical           ,intent(in)  :: ldxout
real(kind=jprb)   ,intent(in)  :: pq(:,:)
real(kind=jprb)   ,intent(out) :: ptheta(:,:)

!-------------------------------------------------------------------------
call potfm(ydfmm%nxy,km,kskip,ydfmm%nx,ldxout,ydfmm%nquad,&
 & ydfmm%rw,pq,ydfmm%rdexp,ydfmm%index,&
 & ydfmm%nclose,ydfmm%ncik,ydfmm%cik,ptheta)
end subroutine seefmm_mulm
!==========================================================================

recursive subroutine potfm(kn,km,kskip,kx,ldxout,kquad,prw,pq,prdexp,kindex,kclosel,kcik,pcik,ptheta)
implicit none

integer(kind=jpim),intent(in)  :: kn
integer(kind=jpim),intent(in)  :: km
integer(kind=jpim),intent(in)  :: kskip
integer(kind=jpim),intent(in)  :: kx
logical           ,intent(in)  :: ldxout
integer(kind=jpim),intent(in)  :: kquad
real(kind=jprb)   ,intent(in)  :: prw(:)
real(kind=jprb)   ,intent(in)  :: pq(:,:)
real(kind=jprb)   ,intent(in)  :: prdexp(:,:)
integer(kind=jpim),intent(in)  :: kindex(:)
integer(kind=jpim),intent(in)  :: kclosel(:)
integer(kind=jpim),intent(in)  :: kcik
real(kind=jprb)   ,intent(in)  :: pcik(:)
real(kind=jprb)   ,intent(out) :: ptheta(:,:)

real(kind=jprb) :: zalpha(kquad,km)
integer(kind=jpim) :: j1,j2,jm,jq,inumc,idist,iquad
integer(kind=jpim) :: iout,iq,i1,i1p1,i1pd,ik1,ix,iy
logical :: lxy,llxy(kn)
 
lxy(ik1) = (ik1 <= kx .eqv. ldxout)
!-------------------------------------------------------------------------

!CALL GSTATS(209,0)
ptheta(:,:)=0.0_JPRB
if(ldxout) then
  ix=0
  iy=-kx
else
  ix=-kx
  iy=0
endif
do j1=1,kn
  i1=kindex(j1)
  llxy(j1)=lxy(i1)
enddo

if(llxy(1)) then
  do jm=1,km,kskip
    zalpha(:,jm)=pq(jm,kindex(1)+ix)
  enddo
else
  zalpha(:,:)=0.0_jprb
endif
!CALL GSTATS(209,1)
!CALL GSTATS(210,0)
do j1=2,kn
  i1=kindex(j1)
  if(llxy(j1) ) then
    if( kskip==1 )then
      do jq=1,kquad
        do jm=1,km
          zalpha(jq,jm)=zalpha(jq,jm)*prdexp(jq,j1)
          zalpha(jq,jm)=zalpha(jq,jm)+pq(jm,i1+ix)
        enddo
      enddo
   else
      do jq=1,kquad
        do jm=1,km,kskip
          zalpha(jq,jm)=zalpha(jq,jm)*prdexp(jq,j1)
          zalpha(jq,jm)=zalpha(jq,jm)+pq(jm,i1+ix)
        enddo
      enddo
   endif
  else
    if( kskip==1 )then
      do jq=1,kquad
        do jm=1,km
          zalpha(jq,jm)=zalpha(jq,jm)*prdexp(jq,j1)
          ptheta(jm,i1+iy)=ptheta(jm,i1+iy)+prw(jq)*zalpha(jq,jm)
        enddo
      enddo
    else
      do jq=1,kquad
        do jm=1,km,kskip
          zalpha(jq,jm)=zalpha(jq,jm)*prdexp(jq,j1)
          ptheta(jm,i1+iy)=ptheta(jm,i1+iy)+prw(jq)*zalpha(jq,jm)
        enddo
      enddo
    endif
  endif
enddo
!CALL GSTATS(210,1)

!CALL GSTATS(211,0)
if(llxy(kn)) then
  do jm=1,km,kskip
    zalpha(:,jm)=pq(jm,kindex(kn)+ix)
  enddo
else
  zalpha(:,:)=0.0
endif
!CALL GSTATS(211,1)
!CALL GSTATS(212,0)
do j1=kn-1,1,-1
  i1=kindex(j1)
  i1p1=kindex(j1+1)
  if(llxy(j1)) then
    if( kskip==1 )then
      do jq=1,kquad
        do jm=1,km
          zalpha(jq,jm)=zalpha(jq,jm)*prdexp(jq,j1+1)
          zalpha(jq,jm)=zalpha(jq,jm)+pq(jm,i1+ix)
        enddo
      enddo
    else
      do jq=1,kquad
        do jm=1,km,kskip
          zalpha(jq,jm)=zalpha(jq,jm)*prdexp(jq,j1+1)
          zalpha(jq,jm)=zalpha(jq,jm)+pq(jm,i1+ix)
        enddo
      enddo
    endif
  else
    if( kskip==1 )then
      do jq=1,kquad
        do jm=1,km
          zalpha(jq,jm)=zalpha(jq,jm)*prdexp(jq,j1+1)
          ptheta(jm,i1+iy)=ptheta(jm,i1+iy)-prw(jq)*zalpha(jq,jm)
        enddo
      enddo
    else
      do jq=1,kquad
        do jm=1,km,kskip
          zalpha(jq,jm)=zalpha(jq,jm)*prdexp(jq,j1+1)
          ptheta(jm,i1+iy)=ptheta(jm,i1+iy)-prw(jq)*zalpha(jq,jm)
        enddo
      enddo
    endif
  endif
enddo
!CALL GSTATS(212,1)


IF(kcik > 0) then
!  CALL GSTATS(213,0)
  inumc=0
  do j1=1,kn-1
    do j2=1,kclosel(j1)
      idist=j2
      i1=kindex(j1)
      i1pd=kindex(j1+idist)
      if(.not.llxy(j1) .and. llxy(j1+idist)) then
        inumc=inumc+1
        do jm=1,km,kskip
          ptheta(jm,i1+iy)=ptheta(jm,i1+iy)-pcik(inumc)*pq(jm,i1pd+ix)
        enddo
      elseif(llxy(j1) .and. .not.llxy(j1+idist)) then
        inumc=inumc+1
        do jm=1,km,kskip
          ptheta(jm,i1pd+iy)=ptheta(jm,i1pd+iy)+pcik(inumc)*pq(jm,i1+ix)
        enddo
      endif
    enddo
  enddo
!  CALL GSTATS(213,1)
endif

end subroutine potfm
!=========================================================================
recursive subroutine suquad(kn,prange,kquad,prw,prt,pr)
implicit none

integer(kind=jpim)        ,intent(in)  :: kn
real(kind=jprb),intent(in)  :: prange
integer(kind=jpim)        ,intent(in) :: kquad
real(kind=jprb),intent(out) :: prw(:)
real(kind=jprb),intent(out) :: prt(:)
real(kind=jprb),intent(out) :: pr

real(kind=jprb) :: za,zb,zs
integer(kind=jpim) :: jm
!-------------------------------------------------------------------------

za=1.0
zb=500.0
zs=zb/prange
pr=za/zs
call wts500(prt,prw,kquad)
do jm=1,kquad
  prw(jm)=prw(jm)*zs
  prt(jm)=prt(jm)*zs
enddo
end subroutine suquad
!==========================================================================

recursive subroutine comb_xy(kx,px,ky,py,kxy,pxy,kindex)

implicit none

integer(kind=jpim), intent(in)  :: kx,ky
real(kind=jprd),    intent(in)  :: px(:)
real(kind=jprb),    intent(in)  :: py(:)
integer(kind=jpim), intent(in)  :: kxy
real(kind=jprb),    intent(out) :: pxy(:)
integer(kind=jpim), intent(out) :: kindex(:)
integer(kind=jpim) :: jxy,ix,iy,iret

!-------------------------------------------------------------------------

pxy(1:kx)=px(1:kx)
pxy(kx+1:kx+ky)=py(1:ky)
!call m01daf(pxy,1,kxy,'D',irank,ifail)
call keysort(iret,pxy,kxy,descending=.true.,index=kindex,init=.true.)
!!$do jxy=1,kxy
!!$  kindex(irank(jxy))=jxy
!!$enddo

end subroutine comb_xy
!==========================================================================
recursive subroutine prepotf(kx,kxy,kquad,prw,prt,pr,pxy,kindex,prdexp,&
 & kclosel,pcik,knocik,pdiff)

implicit none

integer(kind=jpim), intent(in)  :: kx
integer(kind=jpim), intent(in)  :: kxy
integer(kind=jpim), intent(in)  :: kquad
real(kind=jprb),    intent(in)  :: pxy(:)
real(kind=jprb),    intent(in)  :: prw(:)
real(kind=jprb),    intent(in)  :: pr
real(kind=jprb),    intent(in)  :: prt(:)
integer(kind=jpim), intent(in)  :: kindex(:)
real(kind=jprb),    intent(out) :: prdexp(:,:)
integer(kind=jpim), intent(out) :: kclosel(:)
real(kind=jprb),    intent(out) :: pcik(:)
integer(kind=jpim), intent(out) :: knocik
real(kind=jprb),optional, intent(in)  :: pdiff(:,:)

real(kind=jprb) :: zdx
real(kind=jprb) :: zsum
real(kind=jprb) :: zdiff(kxy,kxy)
integer(kind=jpim)  :: jxy,jq,isize,jdist,ixy,ixym1,i1,i1pd,j1,j2
logical :: llexit
!-------------------------------------------------------------------------
if(present(pdiff)) then
  zdiff(:,:)=pdiff(:,:)
else
  do j1=1,kxy
    do j2=1,kxy
      zdiff(j1,j2)=pxy(j1)-pxy(j2)
    enddo
  enddo
endif
do jxy=2,kxy
  ixy=kindex(jxy)
  ixym1=kindex(jxy-1)
  do jq=1,kquad
    prdexp(jq,jxy)=exp(zdiff(ixy,ixym1)*prt(jq))
  enddo
enddo
kclosel(:)=0
knocik=0
isize=size(pcik)
llexit=.true.
do jxy=1,kxy-1
  do jdist=1,kxy-jxy
    i1=kindex(jxy)
    i1pd=kindex(jxy+jdist)
    zdx=zdiff(i1,i1pd)
    if(zdx < pr) then
      llexit=.false.
      kclosel(jxy)=kclosel(jxy)+1
      if((i1 > kx .and. i1pd <= kx) .or. (i1pd > kx .and.  i1 <= kx)) then
        knocik=knocik+1
        zsum=0.0_jprb
        do jq=1,kquad
          zsum=zsum+prw(jq)*exp(-zdx*prt(jq))
        enddo
        pcik(knocik)=1.0_jprb/zdx-zsum
      endif
    else
      exit
    endif
  enddo
  if(knocik > isize) stop ' precompfint : pcik tto small'
enddo

end subroutine prepotf
!==========================================================================

end module seefmm_mix
