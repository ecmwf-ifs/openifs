! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!*** ecsort_shared.h ***

#if INT_VERSION == 4

#define DATA_TYPE            INTEGER(KIND=JPIM)
#define SIZEOF_ME            sizeof_int4
#define KEYSORT_1D           INT4_KEYSORT_1D
#define KEYSORT_1D_DRHOOKSTR 'ECSORT_MIX:INT4_KEYSORT_1D'
#define KEYSORT_2D           INT4_KEYSORT_2D
#define KEYSORT_2D_DRHOOKSTR 'ECSORT_MIX:INT4_KEYSORT_2D'
#define KEYSORT_NUMBER       11
#define RSORT_DRHOOKSTR      'ECSORT_MIX:RSORT32_FUNC_11'
#define USE_RSORT64          0
#define HEAPSORT             INT4_HEAPSORT
#define HEAPSORT_DRHOOKSTR   'ECSORT_MIX:INT4_HEAPSORT'
#define DBGPRINT             INT4_DBGPRINT
#define DBGFMTNUM            1011
#define ECQSORT_DRHOOKSTR    'ECSORT_MIX:INT4_ECQSORT'
#define COUNT_DRHOOKSTR      'ECSORT_MIX:INT4_COUNT'
#define GNOME_DRHOOKSTR      'ECSORT_MIX:INT4_GNOME'

#elif REAL_VERSION == 8

#define DATA_TYPE            REAL(KIND=JPRD)
#define SIZEOF_ME            sizeof_real8
#define KEYSORT_1D           REAL8_KEYSORT_1D
#define KEYSORT_1D_DRHOOKSTR 'ECSORT_MIX:REAL8_KEYSORT_1D'
#define KEYSORT_2D           REAL8_KEYSORT_2D
#define KEYSORT_2D_DRHOOKSTR 'ECSORT_MIX:REAL8_KEYSORT_2D'
#define KEYSORT_NUMBER       12
#define RSORT_DRHOOKSTR      'ECSORT_MIX:RSORT64_12'
#define USE_RSORT64          1
#define HEAPSORT             REAL8_HEAPSORT
#define HEAPSORT_DRHOOKSTR   'ECSORT_MIX:REAL8_HEAPSORT'
#define DBGPRINT             REAL8_DBGPRINT
#define DBGFMTNUM            1012
#define ECQSORT_DRHOOKSTR    'ECSORT_MIX:REAL8_ECQSORT'
#define COUNT_DRHOOKSTR      'ECSORT_MIX:REAL8_COUNT'
#define GNOME_DRHOOKSTR      'ECSORT_MIX:REAL8_GNOME'

#elif REAL_VERSION == 4

#define DATA_TYPE            REAL(KIND=JPRM)
#define SIZEOF_ME            sizeof_real4
#define KEYSORT_1D           REAL4_KEYSORT_1D
#define KEYSORT_1D_DRHOOKSTR 'ECSORT_MIX:REAL4_KEYSORT_1D'
#define KEYSORT_2D           REAL4_KEYSORT_2D
#define KEYSORT_2D_DRHOOKSTR 'ECSORT_MIX:REAL4_KEYSORT_2D'
#define KEYSORT_NUMBER       13
#define RSORT_DRHOOKSTR      'ECSORT_MIX:RSORT32_FUNC_13'
#define USE_RSORT64          0
#define HEAPSORT             REAL4_HEAPSORT
#define HEAPSORT_DRHOOKSTR   'ECSORT_MIX:REAL4_HEAPSORT'
#define DBGPRINT             REAL4_DBGPRINT
#define DBGFMTNUM            1013
#define ECQSORT_DRHOOKSTR    'ECSORT_MIX:REAL4_ECQSORT'
#define COUNT_DRHOOKSTR      'ECSORT_MIX:REAL4_COUNT'
#define GNOME_DRHOOKSTR      'ECSORT_MIX:REAL4_GNOME'

#elif INT_VERSION == 8

#define DATA_TYPE            INTEGER(KIND=JPIB)
#define SIZEOF_ME            sizeof_int8
#define KEYSORT_1D           INT8_KEYSORT_1D
#define KEYSORT_1D_DRHOOKSTR 'ECSORT_MIX:INT8_KEYSORT_1D'
#define KEYSORT_2D           INT8_KEYSORT_2D
#define KEYSORT_2D_DRHOOKSTR 'ECSORT_MIX:INT8_KEYSORT_2D'
#define KEYSORT_NUMBER       14
#define RSORT_DRHOOKSTR      'ECSORT_MIX:RSORT64_14'
#define USE_RSORT64          1
#define HEAPSORT             INT8_HEAPSORT
#define HEAPSORT_DRHOOKSTR   'ECSORT_MIX:INT8_HEAPSORT'
#define DBGPRINT             INT8_DBGPRINT
#define DBGFMTNUM            1014
#define ECQSORT_DRHOOKSTR    'ECSORT_MIX:INT8_ECQSORT'
#define COUNT_DRHOOKSTR      'ECSORT_MIX:INT8_COUNT'
#define GNOME_DRHOOKSTR      'ECSORT_MIX:INT8_GNOME'

#else

  ERROR in programming : No datatype given (should never have ended up here)

#endif

!----------------------------
!--   Public subroutines   --
!----------------------------


SUBROUTINE KEYSORT_1D(rc, a, n, method, descending, index, init)
!-- Please note that we assume that a(:) occupies consecutive memory locations
INTEGER(KIND=JPIM), intent(out)           :: rc
DATA_TYPE         , intent(inout)         :: a(:)
INTEGER(KIND=JPIM), intent(in)            :: n
INTEGER(KIND=JPIM), intent(in), OPTIONAL  :: method
logical, intent(in), OPTIONAL  :: descending
INTEGER(KIND=JPIM), intent(inout), TARGET, OPTIONAL :: index(:)
logical, intent(in), OPTIONAL  :: init
! === END OF INTERFACE BLOCK ===
DATA_TYPE          , allocatable :: aa(:,:)
INTEGER(KIND=JPIM) :: imethod, irev, idummy, index_adj
logical :: LLfast, LLdescending, LLomp_okay, LLinit
INTEGER(KIND=JPIM) :: ITID, ichunk, iret, inumt
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK(KEYSORT_1D_DRHOOKSTR,0,ZHOOK_HANDLE)

rc = 0
if (n <= 0 .or. size(a) <= 0) goto 99

if (present(descending)) then
  LLdescending = descending
else
  LLdescending = .FALSE.
endif

irev = 0
if (LLdescending) irev = 1

ITID = OML_MY_THREAD()
imethod = current_method(ITID)
if (present(method)) then
  imethod = min(max(min_method,method),max_method)
endif

if (imethod /= quicksort_method .and. &
   &imethod /= countingsort_method) then
  LLfast = .FALSE.
else if (imethod == quicksort_method) then
  !-- hasn't been implemented if index is present ;-(
  LLfast = (&
       & .not.present(index) .and. &
       & .not.present(init))
else if (imethod == countingsort_method) then
  !-- index-presence is ok
  LLfast = .TRUE.
endif

if (LLfast) then
  !- Only Quick-sort & CountingSort covered

  if (imethod == quicksort_method) then
    inumt = OML_MAX_THREADS()
    LLomp_okay = (inumt > 1 .and. nomp >= inumt .and. n >= nomp)
    LLomp_okay = (LLomp_okay .and. .not. OML_IN_PARALLEL()) ! Prevents nested OpenMP
    if (LLomp_okay) then
      !-- Max 2-way OpenMP parallelism for now ...
      ichunk = n/2
!$OMP PARALLEL PRIVATE(iret)
!$OMP SECTIONS
!$OMP SECTION
      CALL ecqsortfast(KEYSORT_NUMBER, ichunk, a(1), irev, iret)
!$OMP SECTION
      CALL ecqsortfast(KEYSORT_NUMBER, n-ichunk, a(ichunk+1), irev, iret)
!$OMP END SECTIONS
!$OMP END PARALLEL
      CALL ecmerge2(KEYSORT_NUMBER, 1, ichunk, n-ichunk, a(1), &
           & idummy, 0, 1, irev, idummy, rc)
    else
      CALL ecqsortfast(KEYSORT_NUMBER, n, a(1), irev, rc)
    endif
    GOTO 99

  else if (imethod == countingsort_method) then
    if (.not.present(index)) then
      CALL ec_countingsort(KEYSORT_NUMBER, n, 1, 1, a(1), idummy, 0, 1, irev, rc)
    else
      LLinit = .FALSE.
      if (present(init)) LLinit = init
      if (LLinit) then
        CALL init_index(index, index_adj=-1)
        index_adj = 0
      else
        index_adj = 1
      endif
      CALL ec_countingsort(KEYSORT_NUMBER, n, 1, 1, a(1), index(1), size(index), index_adj, irev, rc)
      if (index_adj == 0) CALL adjust_index(index, +1)
    endif
    GOTO 99

  else
    LLfast = .false.
  endif
endif

!-- LLfast == .FALSE. :

allocate(aa(n,1))

if (LLdescending) then
  aa(1:n,1) = -a(1:n)
else
  aa(1:n,1) = a(1:n)
endif

CALL keysort(rc, aa, n, method=method, index=index, init=init)

if (LLdescending) then
  a(1:n) = -aa(1:n,1)
else
  a(1:n) = aa(1:n,1)
endif

deallocate(aa)

99 continue
IF (LHOOK) CALL DR_HOOK(KEYSORT_1D_DRHOOKSTR,1,ZHOOK_HANDLE,n)
END SUBROUTINE


SUBROUTINE KEYSORT_2D(&
     &rc, a, n,&
     &key, multikey, method,&
     &index, init, transposed)

INTEGER(KIND=JPIM), intent(out)           :: rc
DATA_TYPE         , intent(inout)         :: a(:,:)
INTEGER(KIND=JPIM), intent(in)            :: n
INTEGER(KIND=JPIM), intent(in), OPTIONAL  :: key, method
INTEGER(KIND=JPIM), intent(in), OPTIONAL  :: multikey(:)
logical, intent(in), OPTIONAL  :: transposed
INTEGER(KIND=JPIM), intent(inout), TARGET, OPTIONAL :: index(:)
logical, intent(in), OPTIONAL  :: init
! === END OF INTERFACE BLOCK ===
INTEGER(KIND=JPIM), POINTER :: iindex(:)
INTEGER(KIND=JPIM) :: ikey, istride, imethod, inumkeys, imethod_1st, imethod_rest
INTEGER(KIND=JPIM) :: lda, iptr, i, j, sda, idiff, irev, inumt, jkey, jj, ilastkey
INTEGER(KIND=JPIM) :: j1, j2, jmid, inum, imax, iadd, imod, iret, inc, iamax, ibmax
DATA_TYPE         , allocatable :: data(:)
INTEGER(KIND=JPIM), allocatable :: ikeys(:), ista(:), ichunk(:), irank(:)
logical LLinit, LLdescending, LLtrans, LLomp_okay, LLadjusted, LLdebug, LLomp_prefix
character(len=1) clenv
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_SUBHANDLE
REAL(KIND=JPHOOK) :: ZHOOK_SUBHANDLE0
REAL(KIND=JPHOOK) :: ZHOOK_SUBHANDLE1
REAL(KIND=JPHOOK) :: ZHOOK_SUBHANDLE2
REAL(KIND=JPHOOK) :: ZHOOK_SUBHANDLE3
INTEGER(KIND=JPIM) :: ITID

IF (LHOOK) CALL DR_HOOK(KEYSORT_2D_DRHOOKSTR,0,ZHOOK_HANDLE)

rc = 0
lda = size(a, dim=1)
sda = size(a, dim=2)
if (n <= 0 .or. lda <= 0 .or. sda <= 0) goto 99

inumt = OML_MAX_THREADS()
ITID = OML_MY_THREAD()
imethod = current_method(ITID)
if (present(method)) then
  imethod = min(max(min_method,method),max_method)
endif
imethod_1st = imethod
imethod_rest = imethod

ikey = 1
if (present(key)) ikey = key

if (present(multikey)) then
  allocate(ikeys(size(multikey)))
  ikeys(:) = multikey(:)
else
  allocate(ikeys(1))
  ikeys(1) = ikey
endif
inumkeys = size(ikeys)

!-- Only the RADIX-sort & now also QUICK-sort & CountingSort give the result we want with multiple keys
if (inumkeys > 1 .and. &
    imethod /= radixsort_method .and. &
    imethod /= quicksort_method .and. &
    imethod /= countingsort_method) then
   imethod = default_method
   imethod_1st = imethod
   imethod_rest = imethod
   !-- Since "default_method" may now be [overridden as] HEAP-sort, make sure its then "radixsort_method"
   !   Note: The first sweep may still be e.g. HEAP-sort
   if (imethod /= radixsort_method .and. &
       imethod /= quicksort_method .and. &
       imethod /= countingsort_method) then
     imethod = radixsort_method
     imethod_rest = imethod
   endif
endif

LLinit = .FALSE.
if (present(init)) LLinit = init

if (present(index)) then
  iindex => index(1:n)
else
  allocate(iindex(n))
  LLinit = .TRUE.
endif

if (LLinit) CALL init_index(iindex)

istride = 1
LLtrans = .FALSE.
if (present(transposed)) LLtrans = transposed
if (LLtrans) then
  istride = lda
else if (sda >= 2 .and. lda >= 1) then
!-- Check for presence of sub-array and adjust lda automatically
  call addrdiff(a(1,1),a(1,2),idiff)
  ! lda below: The true leading dimension; overrides sub-arrays one
  lda = idiff/SIZEOF_ME
endif

ilastkey = 0
LLadjusted = .FALSE.
LLomp_prefix = .FALSE.
!$ LLomp_prefix = (istride == 1 .and. nomp >= inumt .and. n >= nomp)
if (LLomp_prefix) then
  call get_environment_variable('EC_SORTING_DEBUG',clenv)
  LLdebug = (clenv == '1' .and. n < 10000)
  if (LLdebug) write(0,*)'>> EC_SORTING_DEBUG=1'
else
  LLdebug = .FALSE.
endif

1000 format(1x,a,2i12,:,/,(10i5))
1001 format(1x,'[#',i2,']:',a,(10i5))
1002 format(1x,'[#',i2,']:',a,:,/,(10i5))
1003 format(1x,'[#',i2,']:',a,2i12,:,/,(10i5))
1004 format(1x,a,:,(10i5))
1005 format(1x,a,i2,1x,a)

imethod = imethod_1st
KEYLOOP: do jkey=inumkeys,1,-1
!--   Sort by the least significant key first
  ikey = abs(ikeys(jkey))
  if (ikey == 0) cycle KEYLOOP

  if (istride == 1) then
    iptr = lda * (ikey - 1) + 1
  else
    iptr = ikey
  endif

  if (LLdebug) then
    write(0,1000) '<BEGIN>iindex(1:n)=',n,sum(iindex(1:n)),iindex(1:n)
    if (LLadjusted) then
      CALL DBGPRINT(-jkey,'<BEGIN>',a,iindex,n,ikey,1,n,1)
    else
      CALL DBGPRINT(-jkey,'<BEGIN>',a,iindex,n,ikey,1,n,0)
    endif
    ilastkey = ikey
  endif

  LLdescending = (ikeys(jkey) < 0)
  irev = 0
  if (LLdescending) irev = 1

  !-- Since "irev" is passed into the ecqsort, no explicit reversing is needed --> savings
  if (imethod == quicksort_method .or. &
      imethod == countingsort_method) LLdescending = .FALSE.

  if (LLdescending) then
    if (istride == 1) then
      a(1:n,ikey) = -a(1:n,ikey)
    else
      a(ikey,1:n) = -a(ikey,1:n)
    endif
    irev = 0 ! prevents use of "reverse" algorithm in ecmerge2 for radix-sort
  endif

  LLomp_okay = LLomp_prefix .and. (inumt > 1) .and. (&
       & imethod == radixsort_method .or. &
       & imethod == quicksort_method .or. &
       & imethod == countingsort_method)
  LLomp_okay = LLomp_okay .and. (.not. OML_IN_PARALLEL()) ! Prevents nested OpenMP

  if (.not.LLomp_okay) then
    select case (imethod)
    case (radixsort_method)
      IF (LHOOK) CALL DR_HOOK(RSORT_DRHOOKSTR,0,ZHOOK_SUBHANDLE0)
#if USE_RSORT64 == 1
      CALL rsort64(KEYSORT_NUMBER, n, istride, iptr, a(1,1), iindex(1), 1, rc)
#else
      CALL rsort32_func(KEYSORT_NUMBER, n, istride, iptr, a(1,1), iindex(1), 1, rc)
#endif
      IF (LHOOK) CALL DR_HOOK(RSORT_DRHOOKSTR,1,ZHOOK_SUBHANDLE0, n)
    case (heapsort_method)
      if (istride == 1) then
        CALL HEAPSORT(KEYSORT_NUMBER, n, a(1:n, ikey), rc, irev, istride, iindex)
      else
        CALL HEAPSORT(KEYSORT_NUMBER, n, a(ikey, 1:n), rc, irev, istride, iindex)
      endif
    case (quicksort_method)
      IF (LHOOK) CALL DR_HOOK(ECQSORT_DRHOOKSTR,0,ZHOOK_SUBHANDLE0)
      CALL ecqsort(KEYSORT_NUMBER, n, istride, iptr, a(1,1), iindex(1), 1, irev, rc)
      IF (LHOOK) CALL DR_HOOK(ECQSORT_DRHOOKSTR,1,ZHOOK_SUBHANDLE0,n)
    case (countingsort_method)
      IF (LHOOK) CALL DR_HOOK(COUNT_DRHOOKSTR,0,ZHOOK_SUBHANDLE0)
      CALL ec_countingsort(KEYSORT_NUMBER, n, istride, iptr, a(1,1), iindex(1), n, 1, irev, rc)
      IF (LHOOK) CALL DR_HOOK(COUNT_DRHOOKSTR,1,ZHOOK_SUBHANDLE0,n)
    case (gnomesort_method)
      IF (LHOOK) CALL DR_HOOK(GNOME_DRHOOKSTR,0,ZHOOK_SUBHANDLE0)
      CALL ecgnomesort(KEYSORT_NUMBER, n, istride, iptr, a(1,1), iindex(1), n, 1, rc)
      IF (LHOOK) CALL DR_HOOK(GNOME_DRHOOKSTR,1,ZHOOK_SUBHANDLE0,n)
    end select

  else ! i.e. LLomp_okay ; radix, quick & counting -sorts only
    if (.not.allocated(ista)) then
      allocate(ista(inumt+1),ichunk(inumt))
      inc = n/inumt
      iadd = 1
      imod = mod(n,inumt)
      if (imod == 0) iadd = 0
      ista(1) = 1
      do j=2,inumt
        ista(j) = ista(j-1) + inc + iadd
        if (iadd > 0 .and. j > imod) iadd = 0
      enddo
      ista(inumt+1) = n + 1
      do j=1,inumt
        ichunk(j) = ista(j+1) - ista(j)
      enddo
      if (LLdebug) then
        write(0,1005) '>> imethod,name=',imethod,method_name(imethod)
        write(0,1004) '>> inumt,n,nomp=',inumt,n,nomp
        write(0,1004) '>> ista(1:inumt+1)=',ista(1:inumt+1)
        write(0,1004) '>> ichunk(1:inumt)=',ichunk(1:inumt)
      endif
      allocate(irank(n))
    endif

    if (LLdebug) write(0,1004) '>>KEYLOOP: jkey,ikey,irev,iptr=',jkey,ikey,irev,iptr

    if (.not.LLadjusted) then ! only once
      if (LLdebug) write(0,1000) '<1>iindex(1:n)=',n,sum(iindex(1:n)),iindex(1:n)
      call adjust_index(iindex, -1) ! Fortran -> C
      if (LLdebug) write(0,1000) '<2>iindex(1:n)=',n,sum(iindex(1:n))+n,iindex(1:n)
      LLadjusted = .TRUE.
    endif

    if (LLdebug) write(0,*)'>> Sorting inumt-chunks in parallel'
!$OMP PARALLEL PRIVATE(j,j1,j2,inum,iret,inc,ITID,ZHOOK_SUBHANDLE1,ZHOOK_SUBHANDLE2)
    IF (LHOOK) CALL DR_HOOK('ECSORT_MIX:KEYSORT_2D>OMPSORT',0,ZHOOK_SUBHANDLE1)
    ITID = OML_MY_THREAD()
!$OMP DO SCHEDULE(DYNAMIC,1)
    do j=1,inumt
      j1 = ista(j)
      inum = ichunk(j)
      j2 = j1 + inum - 1
      inc = j1
      if (LLdebug) write(0,1001) ITID,'j,j1,j2,inum,inc=',j,j1,j2,inum,inc
      if (LLdebug) write(0,1002) ITID,'iindex(j1:j2) > ',iindex(j1:j2)
      select case (imethod)
      case (radixsort_method)
        IF (LHOOK) CALL DR_HOOK(RSORT_DRHOOKSTR,0,ZHOOK_SUBHANDLE2)
#if USE_RSORT64 == 1
        CALL rsort64(KEYSORT_NUMBER, inum, istride, iptr, a(1,1), iindex(j1), 0, iret)
#else
        CALL rsort32_func(KEYSORT_NUMBER, inum, istride, iptr, a(1,1), iindex(j1), 0, iret)
#endif
        IF (LHOOK) CALL DR_HOOK(RSORT_DRHOOKSTR,1,ZHOOK_SUBHANDLE2, inum)
      case (quicksort_method)
        IF (LHOOK) CALL DR_HOOK(ECQSORT_DRHOOKSTR,0,ZHOOK_SUBHANDLE2)
        CALL ecqsort(KEYSORT_NUMBER, inum, istride, iptr, a(1,1), iindex(j1), 0, irev, iret)
        IF (LHOOK) CALL DR_HOOK(ECQSORT_DRHOOKSTR,1,ZHOOK_SUBHANDLE2,inum)
      case (countingsort_method)
        IF (LHOOK) CALL DR_HOOK(COUNT_DRHOOKSTR,0,ZHOOK_SUBHANDLE2)
        CALL ec_countingsort(KEYSORT_NUMBER, inum, istride, iptr, a(1,1), iindex(j1), inum, 0, irev, iret)
        IF (LHOOK) CALL DR_HOOK(COUNT_DRHOOKSTR,1,ZHOOK_SUBHANDLE2,inum)
      end select
      if (LLdebug) write(0,1002) ITID,'iindex(j1:j2) < ',iindex(j1:j2)
    enddo
!$OMP END DO
    IF (LHOOK) CALL DR_HOOK('ECSORT_MIX:KEYSORT_2D>OMPSORT',1,ZHOOK_SUBHANDLE1)
!$OMP END PARALLEL

    if (LLdebug) write(0,1000) '<after_sort>iindex(1:n)=',n,sum(iindex(1:n))+n,iindex(1:n)
    if (LLdebug) CALL DBGPRINT(0,'<after_sort>',a,iindex,n,ikey,1,n,1)
    CALL get_rank(iindex, irank, index_adj=+1)

    if (LLdebug) write(0,*) '>> Merge neighbouring chunks in parallel as much as possible'
    inc = 2
    imax = (inumt+inc-1)/inc
    do jj=1,imax
      if (LLdebug) write(0,1001) jj,'<before_merge> jj,inc,imax,inumt=',jj,inc,imax,inumt
!$OMP PARALLEL PRIVATE(j,j1,j2,inum,iamax,ibmax,jmid,iret,ZHOOK_SUBHANDLE3,ITID)
      IF (LHOOK) CALL DR_HOOK('ECSORT_MIX:KEYSORT_2D>OMPMERGE',0,ZHOOK_SUBHANDLE3)
      ITID = OML_MY_THREAD()
!$OMP DO SCHEDULE(DYNAMIC,1)
      do j=1,inumt,inc
        j1 = j
        j2 = j + inc - 1
        jmid = (j1 + j2)/2 + 1
        j2 = min(j2,inumt)
        jmid = min(jmid,inumt)
        if (LLdebug) write(0,1001) ITID,'j,j1,j2,jmid=',j,j1,j2,jmid
        iamax = ista(jmid) - ista(j1)
        inum = sum(ichunk(j1:j2))
        ibmax = inum - iamax
        if (LLdebug) write(0,1001) ITID,'j,iamax,ibmax,inum=',j,iamax,ibmax,inum
        if (iamax == 0 .or. ibmax == 0 .or. inum == 0) cycle
        j1 = ista(j1)
        j2 = ista(j2+1) - 1
        if (LLdebug) write(0,1001) ITID,'j,j1,j2,inum=',j,j1,j2,inum
        if (LLdebug) write(0,1002) ITID,'iindex(j1:j2) > ',iindex(j1:j2)
        call ecmerge2(KEYSORT_NUMBER, iptr, iamax, ibmax, a(1,1), &
             & iindex(j1), inum, 0, irev, irank(1), iret)
        if (LLdebug) write(0,1002) ITID,'iindex(j1:j2) < ',iindex(j1:j2)
      enddo ! do j=1,inumt,inc
!$OMP END DO
      IF (LHOOK) CALL DR_HOOK('ECSORT_MIX:KEYSORT_2D>OMPMERGE',1,ZHOOK_SUBHANDLE3)
!$OMP END PARALLEL

      if (LLdebug) write(0,1003) jj,'<after_merge>iindex(1:n)=',n,sum(iindex(1:n))+n,iindex(1:n)
      if (LLdebug) CALL DBGPRINT(jj,'<after_merge>',a,iindex,n,ikey,1,n,1)

      inc = inc * 2
    enddo ! do jj=1,imax
    rc = n
  endif ! if (LLomp_okay)

  if (LLdescending) then
    if (istride == 1) then
      a(1:n,ikey) = -a(1:n,ikey)
    else
      a(ikey,1:n) = -a(ikey,1:n)
    endif
  endif

  if (LLadjusted .and. imethod /= imethod_rest) then ! Restore back immediately
    if (LLdebug) write(0,1000) '<3a>iindex(1:n)=',n,sum(iindex(1:n))+n,iindex(1:n)
    call adjust_index(iindex, +1) ! C -> Fortran
    if (LLdebug) write(0,1000) '<4a>iindex(1:n)=',n,sum(iindex(1:n)),iindex(1:n)
    LLadjusted = .FALSE.
  endif

  imethod = imethod_rest
enddo KEYLOOP

deallocate(ikeys)
if (allocated(ista)) deallocate(ista)
if (allocated(ichunk)) deallocate(ichunk)
if (allocated(irank)) deallocate(irank)

if (LLadjusted) then ! Restore back
  if (LLdebug) write(0,1000) '<3b>iindex(1:n)=',n,sum(iindex(1:n))+n,iindex(1:n)
  call adjust_index(iindex, +1) ! C -> Fortran
  if (LLdebug) write(0,1000) '<4b>iindex(1:n)=',n,sum(iindex(1:n)),iindex(1:n)
  LLadjusted = .FALSE.
endif

if (LLdebug) write(0,1000) '<END>iindex(1:n)=',n,sum(iindex(1:n)),iindex(1:n)
if (LLdebug) CALL DBGPRINT(0,'<END>',a,iindex,n,ilastkey,1,n,0)

if (.not.present(index)) then
  LLomp_okay = (nomp >= inumt .and. n >= nomp)
  if (istride == 1) then
    LLomp_okay = (LLomp_okay .and. sda >= inumt .and. .not. OML_IN_PARALLEL()) ! Prevents nested OpenMP
!$OMP PARALLEL PRIVATE(j,data) IF (LLomp_okay)
    allocate(data(n))
!$OMP DO SCHEDULE(DYNAMIC,1)
    do j=1,sda
      data(1:n) = a(iindex(1:n),j)
      a(1:n,j) = data(1:n)
    enddo
!$OMP END DO
    deallocate(data)
!$OMP END PARALLEL
  else
    LLomp_okay = (LLomp_okay .and. lda >= inumt .and. .not. OML_IN_PARALLEL()) ! Prevents nested OpenMP
!$OMP PARALLEL PRIVATE(i,data) IF (LLomp_okay)
    allocate(data(n))
!$OMP DO SCHEDULE(DYNAMIC,1)
    do i=1,lda
      data(1:n) = a(i,iindex(1:n))
      a(i,1:n) = data(1:n)
    enddo
!$OMP END DO
    deallocate(data)
!$OMP END PARALLEL
  endif

  deallocate(iindex)
endif

99 continue
IF (LHOOK) CALL DR_HOOK(KEYSORT_2D_DRHOOKSTR,1,ZHOOK_HANDLE,n*inumkeys)
END SUBROUTINE

!-----------------------------
!--   Private subroutines   --
!-----------------------------

SUBROUTINE DBGPRINT(jj, cdstr, a, index, n, key, k1, k2, kadd)
character(len=*), intent(in) :: cdstr
INTEGER(KIND=JPIM), intent(in) :: jj, n, key, k1, k2, kadd
INTEGER(KIND=JPIM), intent(in) :: index(:)
DATA_TYPE, intent(in) :: a(:,:)
INTEGER(KIND=JPIM) :: i,j
1000 FORMAT(i3,a,5i5)
1011 FORMAT((5i12)) ! integer*4
1012 FORMAT(1p,(5g20.12)) ! real*8
1013 FORMAT(1p,(5g20.12)) ! real*4
1014 FORMAT((5i12)) ! integer*8
WRITE(0,1000) jj,cdstr//': n,key,k1,k2,kadd,a(index(k1:k2)+kadd,:)=',&
     &                     n,key,k1,k2,kadd
do j=k1,k2
  i = index(j)+kadd
  WRITE(0,'(2i6)',advance='no') j,i-kadd
  WRITE(0,DBGFMTNUM) a(i,:)
enddo
END SUBROUTINE

SUBROUTINE HEAPSORT(id, n, a, rc, irev, istride, index)
INTEGER(KIND=JPIM), intent(in) :: id, n, irev, istride
DATA_TYPE, intent(in) :: a(:)
INTEGER(KIND=JPIM), intent(out) :: rc
INTEGER(KIND=JPIM), intent(inout) :: index(:)
INTEGER(KIND=JPIM) :: i,j,right,left,idx
DATA_TYPE :: tmp
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK(HEAPSORT_DRHOOKSTR,0,ZHOOK_HANDLE)
rc = 0
if (n <= 0 .or. size(a) <= 0) goto 99
left  = n/2+1
right = n
LOOP: do
  if (left > 1) then
    left = left - 1
    idx  = index(left)
  else
    idx = index(right)
    index(right) = index(1)
    right = right - 1
    if (right == 1) then
      index(1) = idx
      exit LOOP
    endif
  endif
  tmp = a(idx)
  i = left
  j = 2*left
  do while (j <= right)
    if (j < right) then
      if (a(index(j)) < a(index(j+1))) j = j + 1
    endif
    if (tmp < a(index(j))) then
      index(i) = index(j)
      i = j
      j = 2*j
    else
      j = right + 1
    endif
  enddo
  index(i) = idx
enddo LOOP
rc = n
99 continue
IF (LHOOK) CALL DR_HOOK(HEAPSORT_DRHOOKSTR,1,ZHOOK_HANDLE)
END SUBROUTINE



#ifndef NO_UNDEF

#undef DATA_TYPE
#undef SIZEOF_ME
#undef KEYSORT_1D
#undef KEYSORT_1D_DRHOOKSTR
#undef KEYSORT_2D
#undef KEYSORT_2D_DRHOOKSTR
#undef KEYSORT_NUMBER
#undef RSORT_DRHOOKSTR
#undef USE_RSORT64
#undef QSORTFAST_DRHOOKSTR
#undef HEAPSORT
#undef HEAPSORT_DRHOOKSTR
#undef DBGPRINT
#undef DBGFMTNUM
#undef ECQSORT_DRHOOKSTR
#undef COUNT_DRHOOKSTR
#undef GNOME_DRHOOKSTR

#endif
