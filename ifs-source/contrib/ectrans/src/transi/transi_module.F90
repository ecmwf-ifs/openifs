! (C) Copyright 2014- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!> @file trans_module.F90
!! @brief Fortran layer to trans
!!
!! This file contains the trans_module,
!! which bridges the IFS trans library to
!! the C-API
!!
!! @author Willem Deconinck
!! @date Jul 2014

#define BOOLEAN integer(c_int)

module trans_module

use, intrinsic :: iso_c_binding, only: &
  c_ptr, &
  c_char, &
  c_int, &
  c_size_t, &
  c_float, &
  c_double, &
  c_null_ptr

use, intrinsic :: iso_fortran_env, only: &
  output_unit, &
  error_unit

use MPL_module, only: &
  MPL_INIT, &
  MPL_END, &
  MPL_NPROC, &
  MPL_MYRANK
implicit none

private :: c_ptr
private :: c_char
private :: c_int
private :: c_size_t
private :: c_float
private :: c_double
private :: c_null_ptr

private :: output_unit
private :: error_unit

#if ECTRANS_HAVE_MPI
private :: MPL_INIT
private :: MPL_END
private :: MPL_NPROC
private :: MPL_MYRANK
#endif

public :: Trans_t
public :: DirTrans_t
public :: InvTrans_t
public :: InvTransAdj_t
public :: GathGrid_t
public :: DistGrid_t
public :: GathSpec_t
public :: DistSpec_t
public :: VorDivToUV_t
public :: trans_use_mpi
public :: trans_set_handles_limit
public :: trans_set_radius
public :: trans_set_nprtrv
public :: trans_init
public :: trans_setup
public :: trans_inquire
public :: trans_dirtrans
public :: trans_invtrans
public :: trans_invtrans_adj
public :: trans_distgrid
public :: trans_gathgrid
public :: trans_distspec
public :: trans_gathspec
public :: trans_vordiv_to_UV
public :: trans_delete
public :: trans_finalize
public :: allocate_ptr
public :: access_ptr
public :: free_ptr

private

#include "setup_trans0.h"
#include "trans_end.h"
#include "trans_release.h"
#include "setup_trans.h"
#include "dir_trans.h"
#include "inv_trans.h"
#include "inv_transad.h"
#include "dist_grid.h"
#include "gath_grid.h"
#include "dist_spec.h"
#include "gath_spec.h"
#include "trans_inq.h"
#include "specnorm.h"
#include "vordiv_to_uv.h"

integer, SAVE                      :: TRANS_MAX_HANDLES = 100
integer, SAVE                      :: N_REGIONS_EW
integer, SAVE                      :: N_REGIONS_NS
integer, SAVE, allocatable, target :: N_REGIONS(:)
integer(c_int), SAVE       :: NPRTRV = 1
real(c_double), SAVE       :: RRAD = 6371.22e+03
logical, SAVE              :: is_init = .False.

integer, SAVE              :: trans_out
logical, SAVE              :: close_devnull

#if ECTRANS_HAVE_MPI
logical, SAVE :: USE_MPI   = .True.
#else
logical, SAVE :: USE_MPI   = .False.
#endif


integer, private, parameter :: MAX_STR_LEN = 1024

integer, parameter :: FFT992 = 1
integer, parameter :: FFTW   = 2
integer, parameter :: TRANS_SUCCESS          =  0
integer, parameter :: TRANS_ERROR            = -1
integer, parameter :: TRANS_NOTIMPL          = -2
integer, parameter :: TRANS_MISSING_ARG      = -3
integer, parameter :: TRANS_UNRECOGNIZED_ARG = -4
integer, parameter :: TRANS_STALE_ARG        = -5

!> @brief Interface to the Trans_t struct in transi/trans.h
type, bind(C) :: Trans_t

  ! FILL IN THESE 4 VALUES YOURSELF BEFORE callING trans_setup() */
  integer(c_int) :: ndgl        ! -- Number of lattitudes
  type(c_ptr)    :: nloen       ! -- Number of longitude points for each latitude
                                !    TYPE: INTEGER(1:NDGL)
  integer(c_int) :: nlon        ! -- Number of longitude points for all latitudes
  integer(c_int) :: nsmax       ! -- Spectral truncation wave number
  BOOLEAN        :: lsplit
  integer(c_int) :: llatlon
  integer(c_int) :: flt
  integer(c_int) :: fft
  type(c_ptr)    :: readfp;
  type(c_ptr)    :: writefp;
  type(c_ptr)    :: cache;
  integer(c_size_t) :: cachesize;


  ! PARALLELISATION
  integer(c_int) :: myproc      ! -- Current MPI task (numbering starting at 1)
  integer(c_int) :: nproc       ! -- Number of parallel MPI tasks

  ! MULTI-TRANSFORMS MANAGEMENT
  integer(c_int) :: handle       ! --  Resolution tag for which info is required ,default is the
                                 !     first defined resulution (input)
  BOOLEAN        :: ldlam        ! --  True if the corresponding resolution is LAM, false if it is global

  ! SPECTRAL SPACE
  integer(c_int) :: nspec        ! --  Number of complex spectral coefficients on this PE
  integer(c_int) :: nspec2       ! --  2*nspec
  integer(c_int) :: nspec2g      ! --  global KSPEC2
  integer(c_int) :: nspec2mx     ! --  Maximun KSPEC2 among all PEs
  integer(c_int) :: nump         ! --  Number of spectral waves handled by this PE
  integer(c_int) :: ngptot       ! --  Total number of grid columns on this PE
  integer(c_int) :: ngptotg      ! --  Total number of grid columns on the Globe
  integer(c_int) :: ngptotmx     ! --  Maximum number of grid columns on any of the PEs
  type(c_ptr)    :: ngptotl      ! --  Number of grid columns one each PE
                                 !     TYPE: INTEGER(1:N_REGIONS_NS,1:N_REGIONS_EW)
  type(c_ptr)    :: nmyms        ! --  This PEs spectral zonal wavenumbers
                                 !     TYPE: INTEGER(1:NUMP)
  type(c_ptr)    :: nasm0        ! --  Address in a spectral array of (m, n=m)
                                 !     TYPE: INTEGER(0:NSMAX)
  integer(c_int) :: nprtrw       ! --  Number of processors in A-direction (input)
  type(c_ptr)    :: numpp        ! --  No. of wave numbers each wave set is responsible for.
                                 !     TYPE: INTEGER(1:NPRTRW)
  type(c_ptr)    :: npossp       ! --  Defines partitioning of global spectral fields among PEs
                                 !     TYPE: INTEGER(1:NPRTRW+1)
  type(c_ptr)    :: nptrms       ! --  Pointer to the first wave number of a given a-set
                                 !     TYPE: INTEGER(1:NPRTRW)
  type(c_ptr)    :: nallms       ! --  Wave numbers for all wave-set concatenated together
                                 !     to give all wave numbers in wave-set order
                                 !     TYPE: INTEGER(1:NSMAX+1)
  type(c_ptr)    :: ndim0g       ! --  Defines partitioning of global spectral fields among PEs
                                 !     TYPE: INTEGER(0:NSMAX)
  type(c_ptr)    :: nvalue       ! --  n value for each KSPEC2 spectral coeffient
                                 !     TYPE: INTEGER(1:NSPEC2)


  ! GRIDPOINT SPACE
  integer(c_int) :: n_regions_NS !
  integer(c_int) :: n_regions_EW !
  integer(c_int) :: my_region_NS !
  integer(c_int) :: my_region_EW !
  type(c_ptr)    :: n_regions    ! --  Number of East-West Regions per band of North-South Regions
  type(c_ptr)    :: nfrstlat     ! --  First latitude of each a-set in grid-point space
                                 !     TYPE: INTEGER(1:N_REGIONS_NS)
  type(c_ptr)    :: nlstlat      ! --  Last latitude of each a-set in grid-point space
                                 !     TYPE: INTEGER(1:N_REGIONS_NS)
  integer(c_int) :: nfrstloff    ! --  Offset for first lat of own a-set in grid-point space
  type(c_ptr)    :: nptrlat      ! --  Pointer to the start of each latitude
                                 !     TYPE: INTEGER(1:NDGL)
  type(c_ptr)    :: nptrfrstlat  ! --  Pointer to the first latitude of each a-set in
                                 !     NSTA and NONL arrays
                                 !     TYPE: INTEGER(1:N_REGIONS_NS)
  type(c_ptr)    :: nptrlstlat   ! --  Pointer to the last latitude of each a-set in
                                 !     NSTA and NONL arrays
                                 !     TYPE: INTEGER(1:N_REGIONS_NS)
  integer(c_int) :: nptrfloff    ! --  Offset for pointer to the first latitude of own a-set
                                 !     NSTA and NONL arrays, i.e. nptrfrstlat(myseta)-1
  type(c_ptr)    :: nsta         ! --  Position of first grid column for the latitudes on a
                                 !     processor. The information is available for all processors.
                                 !     The b-sets are distinguished by the last dimension of
                                 !     nsta(). The latitude band for each a-set is addressed by
                                 !     nptrfrstlat(jaset),nptrlstlat(jaset), and
                                 !     nptrfloff=nptrfrstlat(myseta) on this processors a-set.
                                 !     Each split latitude has two entries in nsta(,:) which
                                 !     necessitates the rather complex addressing of nsta(,:)
                                 !     and the overdimensioning of nsta by N_REGIONS_NS.
                                 !     TYPE: INTEGER(1:NDGL+N_REGIONS_NS-1,1:N_REGIONS_EW)
  type(c_ptr)    :: nonl         ! --  Number of grid columns for the latitudes on a processor.
                                 !     Similar to nsta() in data structure.
                                 !     TYPE: INTEGER(1:NDGL+N_REGIONS_NS-1,1:N_REGIONS_EW)
  type(c_ptr)    :: ldsplitlat   ! --  True if latitude is split in grid point space over
                                 !     two a-sets.
                                 !     TYPE: BOOLEAN(1:NDGL)  (BOOLEAN=c_int)

  ! FOURIER SPACE
  integer(c_int) :: nprtrns      ! --  No. of sets in N-S direction (Fourier space)
                                 !     (always equal to NPRTRW)
  type(c_ptr)    :: nultpp       ! --  Number of latitudes for which each a-set is calculating
                                 !     the FFT's.
                                 !     TYPE: INTEGER(1:NPRTRNS)
  type(c_ptr)    :: nptrls       ! --  Pointer to first global latitude of each a-set for which
                                 !     it performs the Fourier calculations
                                 !     TYPE: INTEGER(1:NPRTRNS)
  type(c_ptr)    :: nnmeng       ! --  associated (with NLOENG) cut-off zonal wavenumber
                                 !     TYPE: INTEGER(1:NDGL)

  ! LEGENDRE
  type(c_ptr)    :: rmu          ! --  sin(Gaussian latitudes)
                                 !     TYPE: REAL(1:NDGL)
  type(c_ptr)    :: rgw          ! --  Gaussian weights
                                 !     TYPE: REAL(1:NDGL)
  type(c_ptr)    :: rpnm         ! --  Legendre polynomials
                                 !     TYPE: REAL(1:NLEI3,1:NSPOLEGL)
  integer(c_int) :: nlei3        ! --  First dimension of Legendre polynomials
  integer(c_int) :: nspolegl     ! --  Second dimension of Legendre polynomials
  type(c_ptr)    :: npms         ! --  Adress for legendre polynomial for given M (NSMAX)
                                 !     TYPE: INTEGER(0:NSMAX)
  type(c_ptr)    :: rlapin       ! --  Eigen-values of the inverse Laplace operator
                                 !     TYPE: REAL(-1:NSMAX+2)
  type(c_ptr)    :: ndglu        ! --  Number of active points in an hemisphere for a given wavenumber "m"
                                 !     TYPE: INTEGER(0:NSMAX)
end type Trans_t

type, bind(C) :: DirTrans_t
  type(c_ptr)    :: rgp
  type(c_ptr)    :: rspscalar
  type(c_ptr)    :: rspvor
  type(c_ptr)    :: rspdiv
  integer(c_int) :: nproma
  integer(c_int) :: nscalar
  integer(c_int) :: nvordiv
  integer(c_int) :: ngpblks
  integer(c_int) :: lglobal
  type(c_ptr)    :: trans
  integer(c_int) :: count
end type DirTrans_t

type, bind(C) :: InvTrans_t
  type(c_ptr)    :: rspscalar
  type(c_ptr)    :: rspvor
  type(c_ptr)    :: rspdiv
  type(c_ptr)    :: rgp
  integer(c_int) :: nproma
  integer(c_int) :: nscalar
  integer(c_int) :: nvordiv
  integer(c_int) :: lscalarders
  integer(c_int) :: luvder_EW
  integer(c_int) :: lvordivgp
  integer(c_int) :: ngpblks
  integer(c_int) :: lglobal
  type(c_ptr)    :: trans
  integer(c_int) :: count
end type InvTrans_t

type, bind(C) :: InvTransAdj_t
  type(c_ptr)    :: rspscalar
  type(c_ptr)    :: rspvor
  type(c_ptr)    :: rspdiv
  type(c_ptr)    :: rgp
  integer(c_int) :: nproma
  integer(c_int) :: nscalar
  integer(c_int) :: nvordiv
  integer(c_int) :: lscalarders
  integer(c_int) :: luvder_EW
  integer(c_int) :: lvordivgp
  integer(c_int) :: ngpblks
  integer(c_int) :: lglobal
  type(c_ptr)    :: trans
  integer(c_int) :: count
end type InvTransAdj_t

type, bind(C) :: DistGrid_t
  type(c_ptr)    :: rgpg
  type(c_ptr)    :: rgp
  type(c_ptr)    :: nfrom
  integer(c_int) :: nproma
  integer(c_int) :: nfld
  integer(c_int) :: ngpblks
  type(c_ptr)    :: trans
  integer(c_int) :: count
end type DistGrid_t

type, bind(C) :: GathGrid_t
  type(c_ptr)    :: rgpg
  type(c_ptr)    :: rgp
  type(c_ptr)    :: nto
  integer(c_int) :: nproma
  integer(c_int) :: nfld
  integer(c_int) :: ngpblks
  type(c_ptr)    :: trans
  integer(c_int) :: count
end type GathGrid_t

type, bind(C) :: DistSpec_t
  type(c_ptr)    :: rspecg
  type(c_ptr)    :: rspec
  type(c_ptr)    :: nfrom
  integer(c_int) :: nfld
  type(c_ptr)    :: trans
  integer(c_int) :: count
end type DistSpec_t

type, bind(C) :: GathSpec_t
  type(c_ptr)    :: rspecg
  type(c_ptr)    :: rspec
  type(c_ptr)    :: nto
  integer(c_int) :: nfld
  type(c_ptr)    :: trans
  integer(c_int) :: count
end type GathSpec_t

type, bind(C) :: VorDivToUV_t
  type(c_ptr)    :: rspvor
  type(c_ptr)    :: rspdiv
  type(c_ptr)    :: rspu
  type(c_ptr)    :: rspv
  integer(c_int) :: nfld
  integer(c_int) :: nsmax
  integer(c_int) :: ncoeff
  integer(c_int) :: count
end type VorDivToUV_t

type, bind(C) :: SpecNorm_t
  type(c_ptr)    :: rspec
  integer(c_int) :: nmaster
  type(c_ptr)    :: rmet
  type(c_ptr)    :: rnorm
  integer(c_int) :: nfld
  type(c_ptr)    :: trans
  integer(c_int) :: count
end type SpecNorm_t

interface
  subroutine transi_malloc_bool(ptr,len) bind(C,name="transi_malloc_bool")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int
    type(c_ptr) :: ptr
    integer(c_int), value :: len
  end subroutine transi_malloc_bool
  subroutine transi_malloc_int(ptr,len) bind(C,name="transi_malloc_int")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int
    type(c_ptr) :: ptr
    integer(c_int), value :: len
  end subroutine transi_malloc_int
  subroutine transi_malloc_float(ptr,len) bind(C,name="transi_malloc_float")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int
    type(c_ptr) :: ptr
    integer(c_int), value :: len
  end subroutine transi_malloc_float
  subroutine transi_malloc_double(ptr,len) bind(C,name="transi_malloc_double")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int
    type(c_ptr) :: ptr
    integer(c_int), value :: len
  end subroutine transi_malloc_double
  subroutine transi_free(ptr) bind(C,name="transi_free")
    use, intrinsic :: iso_c_binding, only: c_ptr
    type(c_ptr), intent(in) :: ptr
  end subroutine transi_free
end interface


interface trans_inquire
  module procedure trans_inquire_cstr
  module procedure trans_inquire_fstr
end interface trans_inquire


interface allocate_ptr
!  module procedure allocate_bool1_ptr
  module procedure allocate_int1_ptr
  module procedure allocate_int2_ptr
  module procedure allocate_double1_ptr
  module procedure allocate_double2_ptr
end interface allocate_ptr

interface access_ptr
!  module procedure allocate_bool1_ptr
  module procedure allocate_int1_ptr
  module procedure allocate_int2_ptr
  module procedure allocate_double1_ptr
  module procedure allocate_double2_ptr
end interface access_ptr

contains

! =============================================================================
! From fckit_c_interop module

function c_str_to_string(s) result(string)
  use, intrinsic :: iso_c_binding
  character(kind=c_char,len=1), intent(in) :: s(*)
  character(len=:), allocatable :: string
  integer i, nchars
  i = 1
  do
     if (s(i) == c_null_char) exit
     i = i + 1
  enddo
  nchars = i - 1  ! Exclude null character from Fortran string
  allocate(character(len=nchars) :: string)
  do i=1,nchars
    string(i:i) = s(i)
  enddo
end function

function c_ptr_to_string(cptr) result(string)
  use, intrinsic :: iso_c_binding
  type(c_ptr), intent(in) :: cptr
  character(kind=c_char,len=:), allocatable :: string
  character, pointer  :: s(:)
  call c_f_pointer ( cptr , s, (/MAX_STR_LEN/) )
  string = c_str_to_string(s)
end function

! =============================================================================

subroutine to_lower(str)
  character(*), intent(in out) :: str
  integer :: i

  do i = 1, len(str)
    select case(str(i:i))
      case("A":"Z")
        str(i:i) = achar(iachar(str(i:i))+32)
    end select
  end do
end subroutine to_lower

subroutine transi_error(err_msg)
  character(len=*), intent(in) :: err_msg
  write(error_unit,'(A)') err_msg
end subroutine

function trans_set_handles_limit(limit) bind(C,name="trans_set_handles_limit")
  integer(c_int) :: trans_set_handles_limit
  integer(c_int), value, intent(in) :: limit
  TRANS_MAX_HANDLES = limit
  trans_set_handles_limit = TRANS_SUCCESS
end function

function trans_set_radius(radius) bind(C,name="trans_set_radius")
  integer(c_int) :: trans_set_radius
  real(c_double), value, intent(in) :: radius
  RRAD = radius
  trans_set_radius = TRANS_SUCCESS
end function

function trans_set_nprtrv(kprtrv) bind(C,name="trans_set_nprtrv")
  integer(c_int) :: trans_set_nprtrv
  integer(c_int), value, intent(in) :: kprtrv
  NPRTRV = kprtrv
  trans_set_nprtrv = TRANS_SUCCESS
end function

function trans_use_mpi(lmpi) bind(C,name="trans_use_mpi")
  integer(c_int) :: trans_use_mpi
  integer(c_int), value, intent(in) :: lmpi
#if ECTRANS_HAVE_MPI
  if( lmpi == 0 ) then
    USE_MPI = .False.
  else
    USE_MPI = .True.
  endif
#endif
  trans_use_mpi = TRANS_SUCCESS
end function

function devnull(opened)
  integer :: devnull
  logical, intent(out), optional :: opened
  integer :: devnull_unit
  inquire(file="/dev/null", number=devnull_unit)
  if( devnull_unit == 5 ) devnull_unit = -1 ! Willem D: Bug in gfortran and openmpi
  if( devnull_unit == -1 ) then
    devnull_unit = 777
    open( unit=devnull_unit, file="/dev/null" )
    if( present(opened) ) opened = .true.
  else
    if( present(opened) ) opened = .false.
  endif
  devnull = devnull_unit
end function

function trans_init() bind(C,name="trans_init") result(iret)
  integer(c_int) :: iret
  integer :: NPRTRW, NPRGPNS
  integer, allocatable :: I_REGIONS(:)
  logical :: LMPOFF
  integer :: devnull_unit

  LMPOFF = .not. USE_MPI

  trans_out = devnull( opened=close_devnull )

  if( USE_MPI ) then
    call MPL_INIT(KOUTPUT=0,KUNIT=trans_out,LDINFO=.False.)
    allocate( I_REGIONS(MPL_NPROC()) )
    NPRGPNS = MPL_NPROC()
    NPRTRW = MPL_NPROC()/NPRTRV;
  else
    allocate( I_REGIONS(1) )
    NPRGPNS = 1
    NPRTRW = 1;
  endif

  call SETUP_TRANS0(KOUT=trans_out,KERR=error_unit,KPRINTLEV=0,KMAX_RESOL=TRANS_MAX_HANDLES,&
    &               KPRTRW=NPRTRW,  LDEQ_REGIONS=.True.,KPRGPNS=NPRGPNS,KPRGPEW=1,&
    &               PRAD=RRAD, K_REGIONS_NS=N_REGIONS_NS,K_REGIONS_EW=N_REGIONS_EW,K_REGIONS=I_REGIONS,&
    &               LDMPOFF=LMPOFF )
  allocate(N_REGIONS(1:N_REGIONS_NS))
  N_REGIONS(1:N_REGIONS_NS)=I_REGIONS(1:N_REGIONS_NS)
  is_init = .True.

  iret = TRANS_SUCCESS

end function trans_init


function trans_setup(trans) bind(C,name="trans_setup") result(iret)
  use, intrinsic :: iso_c_binding
  integer(c_int) :: iret
  type(Trans_t), intent(inout) :: trans
  integer(c_int), pointer :: nloen(:)
  integer(c_int), pointer :: n_regions_fptr(:)
  logical, parameter :: lkeeprpnm =.False.
  logical, parameter :: luserpnm  =.False. ! Don't use Belusov algorithm (uses twice the memory)
  logical :: ldlam ! output
  logical :: lgridonly, lsplit, lusefftw !input
  logical :: lspeconly ! only
  logical :: llatlon ! input
  logical :: llatlonshift ! input
  integer(c_int) :: nlon
  integer(c_int) :: err
  character(len=MAX_STR_LEN) :: readfp, writefp
  logical :: luseflt

  iret = TRANS_SUCCESS

  lsplit = .False.
  if( trans%lsplit /= 0 ) lsplit = .True.

  if( trans%fft == FFTW )   lusefftw   = .True.

  llatlon = .False.
  llatlonshift = .False.
  if( trans%llatlon /= 0 ) llatlon = .True.
  if( trans%llatlon == 2 ) llatlonshift = .True.

  if ( .not. is_init ) then
    err = trans_init()
  endif

  lspeconly = .False.
  if( trans%ndgl < 0 ) then
    lspeconly = .true.
    trans%ndgl = 2
  endif

  nlon = trans%nlon
  if( nlon < 0 .and. trans%ndgl >= 0 ) then
    if( c_associated( trans%nloen ) ) then
      call C_F_POINTER( trans%nloen, nloen, (/trans%ndgl/) )
      nlon = nloen(1)
    else
      nlon = 2*trans%ndgl
    endif
  endif

  lgridonly = .False.
  if( trans%nsmax < 0 ) then
    lgridonly = .true.
  endif

  if( lgridonly .and. lspeconly ) then
    write(error_unit,'(A)') "trans_setup: ERROR: Cannot setup with both lgridonly and lspeconly. Make up your mind."
    iret = TRANS_ERROR
    return
  endif


  writefp=""
  if( c_associated(trans%writefp) ) then
    writefp = c_ptr_to_string(trans%writefp)
    !call cptr_to_f_string(trans%writefp,writefp)
  endif

  readfp=""
  if( c_associated(trans%readfp) ) then
    readfp = c_ptr_to_string(trans%readfp)
    !call cptr_to_f_string(trans%readfp,readfp)
  endif

  if ( trans%cachesize > 0 ) then
    if( .not. c_associated( trans%cache ) ) then
      write(error_unit,'(A)') "Cache memory was not allocated"
      iret = TRANS_MISSING_ARG
      return
    endif
  endif

#define LATLON_FLAGS LDLL=llatlon, LDSHIFTLL=llatlonshift,

! if( trans%flt > 0 .and. trans%nsmax+1 > trans%ndgl ) then
!   write(error_unit,'(A)') "trans_setup: WARNING: A bug in trans doesn't allow to use FLT with "&
!     & //                  "truncation (nsmax+1) > nb_latitudes (ndgl). Continuing with FLT=OFF."
! endif
! if( trans%nsmax+1 > trans%ndgl ) then
!   trans%flt = 0
! endif

 luseflt = .False.
 if( trans%flt > 0 ) luseflt = .True.

 if( .not. c_associated( trans%nloen ) ) then
   ! Setup that involves latlon requires no nloen

   if( len_trim(readfp) > 0 ) then
     if( trans%flt >= 0 ) then

       ! LONLAT; Impose FLT; READ coeffs from file
       call SETUP_TRANS( LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KDLON=nlon, &
         & CDIO_LEGPOL="readf", &
         & CDLEGPOLFNAME=readfp, &
         & LDUSEFLT=luseflt )

     else

       ! LONLAT; Default FLT; READ coeffs from file
       call SETUP_TRANS( LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KDLON=nlon, &
         & CDIO_LEGPOL="readf", &
         & CDLEGPOLFNAME=readfp )

     endif
   elseif( len_trim(writefp) > 0 ) then
     if( trans%flt >= 0 ) then

       ! LONLAT; Impose FLT; WRITE coeffs to file
       call SETUP_TRANS( LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KDLON=nlon, &
         & CDIO_LEGPOL="writef", &
         & CDLEGPOLFNAME=writefp, &
         & LDUSEFLT=luseflt )

     else

       ! LONLAT; Impose FLT; WRITE coeffs to file
       call SETUP_TRANS( LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KDLON=nlon, &
         & CDIO_LEGPOL="writef", &
         & CDLEGPOLFNAME=writefp )

     endif
   elseif( trans%cachesize > 0 ) then
     if( trans%flt >= 0 ) then

       ! LONLAT; Impose FLT; read CACHED coefficients
       call SETUP_TRANS( LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KDLON=nlon, &
         & CDIO_LEGPOL="membuf", &
         & KLEGPOLPTR=trans%cache, &
         & KLEGPOLPTR_LEN=trans%cachesize, &
         & LDUSEFLT=luseflt )

     else

       ! LONLAT; Default FLT; read CACHED coefficients
       call SETUP_TRANS( LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KDLON=nlon, &
         & CDIO_LEGPOL="membuf", &
         & KLEGPOLPTR=trans%cache, &
         & KLEGPOLPTR_LEN=trans%cachesize )

     endif

   else

     if( trans%flt >= 0 ) then

       ! LONLAT; Impose FLT
       call SETUP_TRANS( LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KDLON=nlon, &
         & LDUSEFLT=luseflt )

     else

       ! LONLAT; Default FLT
       call SETUP_TRANS( LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KDLON=nlon )

     endif

   endif

 else ! we have nloen

   call C_F_POINTER( trans%nloen, nloen, (/trans%ndgl/) )

   if( len_trim(readfp) > 0 ) then

     if( trans%flt >= 0 ) then

       ! REDUCEDGAUSSIANGRID; Impose FLT; READ coefficients
       call SETUP_TRANS( LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KLOEN=nloen, CDIO_LEGPOL="readf", &
         & CDLEGPOLFNAME=trim(readfp),&
         & LDUSEFLT=luseflt )

     else

       ! REDUCEDGAUSSIANGRID; Default FLT; READ coefficients
       call SETUP_TRANS( LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KLOEN=nloen, CDIO_LEGPOL="readf", &
         & CDLEGPOLFNAME=trim(readfp) )

     endif
   elseif( len_trim(writefp) > 0 ) then

     if( trans%flt >= 0 ) then

       ! REDUCEDGAUSSIANGRID; Impose FLT; WRITE coefficients
       call SETUP_TRANS(  LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KLOEN=nloen, &
         & CDIO_LEGPOL="writef", &
         & CDLEGPOLFNAME=trim(writefp), &
         & LDUSEFLT=luseflt )

     else

       ! REDUCEDGAUSSIANGRID; Default FLT; READ coefficients
       call SETUP_TRANS(  LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KLOEN=nloen, &
         & CDIO_LEGPOL="writef", &
         & CDLEGPOLFNAME=trim(writefp) )

     endif

   elseif( trans%cachesize > 0 ) then

     if( trans%flt >= 0 ) then

       ! REDUCEDGAUSSIANGRID; Default FLT; read CACHED coefficients
       call SETUP_TRANS(  LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KLOEN=nloen, &
         & CDIO_LEGPOL="membuf", &
         & KLEGPOLPTR=trans%cache, &
         & KLEGPOLPTR_LEN=trans%cachesize, &
         & LDUSEFLT=luseflt )

     else

       ! REDUCEDGAUSSIANGRID; Impose FLT; read CACHED coefficients
       call SETUP_TRANS(  LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KLOEN=nloen, &
         & CDIO_LEGPOL="membuf", &
         & KLEGPOLPTR=trans%cache, &
         & KLEGPOLPTR_LEN=trans%cachesize )

     endif

   else
     if( trans%flt >= 0 ) then

       ! REDUCEDGAUSSIANGRID; Impose FLT
       call SETUP_TRANS(  LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KLOEN=nloen, &
         & LDUSEFLT=luseflt )

     else

       ! REDUCEDGAUSSIANGRID; Default FLT
       call SETUP_TRANS(  LATLON_FLAGS &
         & KSMAX=trans%nsmax, &
         & KRESOL=trans%handle, &
         & KDGL=trans%ndgl, &
         & LDGRIDONLY=LGRIDONLY, &
         & LDSPSETUPONLY=LSPECONLY, &
         & LDSPLIT=LSPLIT, &
         & LDKEEPRPNM=LKEEPRPNM, &
         & LDUSERPNM=LUSERPNM, &
         & LDUSEFFTW=LUSEFFTW, &
         & KLOEN=nloen )
     endif
   endif
 endif

  if( USE_MPI ) then
    trans%myproc = MPL_MYRANK()
    trans%nproc = MPL_NPROC()
  else
    trans%myproc = 1
    trans%nproc = 1
  endif
  trans%n_regions_NS = N_REGIONS_NS
  trans%n_regions_EW = N_REGIONS_EW

  call TRANS_INQ( KRESOL    = trans%handle, &
                  KMY_REGION_NS = trans%my_region_NS, &
                  KMY_REGION_EW = trans%my_region_EW, &
                  KSPEC     = trans%nspec, &
                  KSPEC2    = trans%nspec2, &
                  KSPEC2G   = trans%nspec2g, &
                  KSPEC2MX  = trans%nspec2mx, &
                  KNUMP     = trans%nump, &
                  KGPTOT    = trans%ngptot, &
                  KGPTOTG   = trans%ngptotg, &
                  KGPTOTMX  = trans%ngptotmx, &
                  KFRSTLOFF = trans%nfrstloff, &
                  KPTRFLOFF = trans%nptrfloff, &
                  KPRTRW    = trans%nprtrw, &
                  KLEI3     = trans%nlei3, &
                  KSPOLEGL  = trans%nspolegl, &
                  LDLAM     = ldlam &
                )


  if( trans%llatlon == 1 ) trans%ngptotg = trans%ngptotg-nlon

  trans%ldlam = 0
  if( ldlam ) trans%ldlam = 1

  trans%nprtrns = trans%nprtrw

  trans%ngptotl     = C_NULL_PTR
  trans%nmyms       = C_NULL_PTR
  trans%nasm0       = C_NULL_PTR
  trans%numpp       = C_NULL_PTR
  trans%npossp      = C_NULL_PTR
  trans%nptrms      = C_NULL_PTR
  trans%nallms      = C_NULL_PTR
  trans%ndim0g      = C_NULL_PTR
  trans%n_regions   = C_NULL_PTR
  trans%nfrstlat    = C_NULL_PTR
  trans%nlstlat     = C_NULL_PTR
  trans%nptrlat     = C_NULL_PTR
  trans%nptrfrstlat = C_NULL_PTR
  trans%nptrlstlat  = C_NULL_PTR
  trans%nsta        = C_NULL_PTR
  trans%nonl        = C_NULL_PTR
  trans%nultpp      = C_NULL_PTR
  trans%nptrls      = C_NULL_PTR
  trans%nnmeng      = C_NULL_PTR
  trans%rmu         = C_NULL_PTR
  trans%rgw         = C_NULL_PTR
  trans%rpnm        = C_NULL_PTR
  trans%npms        = C_NULL_PTR
  trans%ndglu       = C_NULL_PTR
  trans%rlapin      = C_NULL_PTR
  trans%nvalue      = C_NULL_PTR
  trans%ldsplitlat  = C_NULL_PTR

  call allocate_ptr( trans%n_regions,N_REGIONS_NS, n_regions_fptr )
  n_regions_fptr(:) = N_REGIONS(:)

end function trans_setup

function trans_inquire_cstr(trans,vars) bind(C,name="trans_inquire") result(iret)
  integer(c_int) :: iret
  type(Trans_t), intent(inout) :: trans
  character(len=1,kind=c_char), dimension(*), intent(in) :: vars
  character(len=MAX_STR_LEN,kind=c_char) :: vars_fstr
  vars_fstr = c_str_to_string(vars)
  iret = trans_inquire_fstr(trans,vars_fstr)
end function trans_inquire_cstr

function trans_inquire_fstr(trans,vars_fstr) result(iret)
  integer(c_int) :: iret
  type(Trans_t), intent(inout) :: trans
  character(len=*) :: vars_fstr
  character(20) :: var_arr(30), var
  integer :: nvars, jvar
  !logical(c_bool), pointer :: bool1(:)
  integer(c_int), pointer :: int1(:), int2(:,:)
  real(c_double), pointer :: double1(:), double2(:,:)
  logical, allocatable :: booltmp(:)

  nvars = count(transfer(vars_fstr, 'a', len(vars_fstr)) == ",") + 1
  read(vars_fstr, *) var_arr(1:nvars)

  do jvar=1,nvars
    var = trim(var_arr(jvar))
    call to_lower(var)

    if    ( var == "numpp" ) then
      call allocate_ptr( trans%numpp, trans%nprtrw, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KUMPP=int1 )

    elseif( var == "ngptotl" ) then
      call allocate_ptr( trans%ngptotl, trans%n_regions_NS, trans%n_regions_EW, int2 )
      call TRANS_INQ( KRESOL=trans%handle,  KGPTOTL=int2 )

    elseif( var == "nmyms" ) then
      call allocate_ptr( trans%nmyms, trans%nump, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KMYMS=int1 )

    elseif( var == "nasm0" ) then
      call allocate_ptr( trans%nasm0, trans%nsmax+1, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KASM0=int1 )

    elseif( var == "npossp" ) then
      call allocate_ptr( trans%npossp, trans%nprtrw+1, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KPOSSP=int1 )

    elseif( var == "nptrms" ) then
      call allocate_ptr( trans%nptrms, trans%nprtrw, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KPTRMS=int1 )

    elseif( var == "nallms" ) then
      call allocate_ptr( trans%nallms, trans%nsmax+1, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KALLMS=int1 )

    elseif( var == "ndim0g" ) then
      call allocate_ptr( trans%ndim0g, trans%nsmax+1, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KDIM0G=int1 )

    elseif( var == "nvalue" ) then
      call allocate_ptr( trans%nvalue, trans%nspec2, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KNVALUE=int1 )

    elseif( var == "nfrstlat" ) then
      call allocate_ptr( trans%nfrstlat, trans%n_regions_NS, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KFRSTLAT=int1 )

    elseif( var == "nlstlat" ) then
      call allocate_ptr( trans%nlstlat, trans%n_regions_NS, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KLSTLAT=int1 )

    elseif( var == "nptrlat" ) then
      call allocate_ptr( trans%nptrlat, trans%ndgl, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KPTRLAT=int1 )

    elseif( var == "nptrfrstlat" ) then
      call allocate_ptr( trans%nptrfrstlat ,trans%n_regions_ns, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KPTRFRSTLAT=int1 )

    elseif( var == "nptrlstlat" ) then
      call allocate_ptr( trans%nptrlstlat, trans%n_regions_ns, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KPTRLSTLAT=int1 )

    elseif( var == "nsta" ) then
      call allocate_ptr( trans%nsta, trans%ndgl+trans%n_regions_NS-1, trans%n_regions_EW, int2 )
      call TRANS_INQ( KRESOL=trans%handle,  KSTA=int2 )

    elseif( var == "nonl" ) then
      call allocate_ptr( trans%nonl, trans%ndgl+trans%n_regions_NS-1, trans%n_regions_EW, int2 )
      call TRANS_INQ( KRESOL=trans%handle,  KONL=int2 )

    elseif( var == "ldsplitlat" ) then
      !!call allocate_ptr( trans%nonl, trans%ndgl, bool1 )
      !allocate( booltmp(trans%ndgl) )
      !call TRANS_INQ( KRESOL=trans%handle,  LDSPLITLAT=booltmp )
      !bool1(:) = booltmp(:)
      !deallocate( booltmp )
      iret = TRANS_NOTIMPL
      return

    elseif( var == "nultpp" ) then
      call allocate_ptr( trans%nultpp, trans%nprtrns, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KULTPP=int1 )

    elseif( var == "nptrls" ) then
      call allocate_ptr( trans%nptrls, trans%nprtrns, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KPTRLS=int1 )

    elseif( var == "nnmeng" ) then
      call allocate_ptr( trans%nnmeng, trans%ndgl, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KNMENG=int1 )

    elseif( var == "rmu" ) then
      call allocate_ptr( trans%rmu, trans%ndgl, double1 )
      call TRANS_INQ( KRESOL=trans%handle,  PMU=double1 )

    elseif( var == "rgw" ) then
      call allocate_ptr( trans%rgw, trans%ndgl, double1 )
      call TRANS_INQ( KRESOL=trans%handle,  PGW=double1 )

    elseif( var == "rpnm" ) then
      call allocate_ptr( trans%rpnm, trans%nlei3, trans%nspolegl, double2 )
      call TRANS_INQ( KRESOL=trans%handle,  PRPNM=double2 )

    elseif( var == "npms" ) then
      call allocate_ptr( trans%npms, trans%nsmax+1, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KPMS=int1 )

    elseif( var == "rlapin" ) then
      call allocate_ptr( trans%rlapin, trans%nsmax+4, double1 )
      call TRANS_INQ( KRESOL=trans%handle,  PLAPIN=double1 )

    elseif( var == "ndglu" ) then
      call allocate_ptr( trans%ndglu, trans%nsmax+1, int1 )
      call TRANS_INQ( KRESOL=trans%handle,  KDGLU=int1 )

    elseif(       var /= "ndgl"         &
        &   .and. var /= "nsmax"        &
        &   .and. var /= "myproc"       &
        &   .and. var /= "nproc"        &
        &   .and. var /= "ldlam"        &
        &   .and. var /= "nspec"        &
        &   .and. var /= "nspec2"       &
        &   .and. var /= "nspec2g"      &
        &   .and. var /= "nspec2mx"     &
        &   .and. var /= "nump"         &
        &   .and. var /= "ngptot"       &
        &   .and. var /= "ngptotg"      &
        &   .and. var /= "ngptotmx"     &
        &   .and. var /= "n_regions_ns" &
        &   .and. var /= "n_regions_ew" &
        &   .and. var /= "my_region_ns" &
        &   .and. var /= "my_region_ew" &
        &   .and. var /= "nfrstloff"    &
        &   .and. var /= "nptrfloff"    &
        &   .and. var /= "nprtrns"      &
        &   .and. var /= "nlei3"        &
        &   .and. var /= "nspolegl"     ) then
      write(error_unit,*) "trans_inqure: ERROR: unrecognized variable ", var
      iret = TRANS_UNRECOGNIZED_ARG
      return
    endif

  enddo

  iret = TRANS_SUCCESS
end function trans_inquire_fstr

subroutine free_ptr(ptr)
  use, intrinsic :: iso_c_binding
  type(c_ptr) :: ptr
  if( c_associated( ptr ) ) then
    call transi_free(ptr)
    ptr = c_null_ptr
  endif
end subroutine free_ptr

!subroutine allocate_bool1_ptr(ptr,len,bool1)
!  use, intrinsic :: iso_c_binding
!  type(c_ptr) :: ptr
!  integer(c_int) :: len
!  logical(c_bool), pointer :: bool1(:)
!  if( .not. c_associated( ptr ) ) call transi_malloc_bool(ptr,len)
!  call c_f_pointer( ptr, bool1, (/len/) )
!end subroutine allocate_bool1_ptr

subroutine allocate_int1_ptr(ptr,len,int1)
  use, intrinsic :: iso_c_binding
  type(c_ptr) :: ptr
  integer(c_int) :: len
  integer(c_int), pointer :: int1(:)
  if( .not. c_associated( ptr ) ) call transi_malloc_int(ptr,len)
  call c_f_pointer( ptr, int1, (/len/) )
end subroutine allocate_int1_ptr

subroutine allocate_int2_ptr(ptr,len1,len2,int2)
  use, intrinsic :: iso_c_binding
  type(c_ptr) :: ptr
  integer(c_int) :: len1, len2
  integer(c_int), pointer :: int2(:,:)
  if( .not. c_associated( ptr ) ) call transi_malloc_int(ptr,len1*len2)
  call c_f_pointer( ptr, int2, (/len1,len2/) )
end subroutine allocate_int2_ptr

subroutine allocate_double1_ptr(ptr,len,double1)
  use, intrinsic :: iso_c_binding
  type(c_ptr) :: ptr
  integer(c_int) :: len
  real(c_double), pointer :: double1(:)
  if( .not. c_associated( ptr ) ) call transi_malloc_double(ptr,len)
  call c_f_pointer( ptr, double1, (/len/) )
end subroutine allocate_double1_ptr

subroutine allocate_double2_ptr(ptr,len1,len2,double2)
  use, intrinsic :: iso_c_binding
  type(c_ptr) :: ptr
  integer(c_int) :: len1, len2
  real(c_double), pointer :: double2(:,:)
  if( .not. c_associated( ptr ) ) call transi_malloc_double(ptr,len1*len2)
  call c_f_pointer( ptr, double2, (/len1,len2/) )
end subroutine allocate_double2_ptr

function trans_delete(trans) bind(C,name="trans_delete")
  use, intrinsic :: iso_c_binding
  integer(c_int) :: trans_delete
  type(Trans_t), intent(inout) :: trans
  call free_ptr( trans%nloen       )
  call free_ptr( trans%readfp      )
  call free_ptr( trans%writefp     )
  call free_ptr( trans%ngptotl     )
  call free_ptr( trans%nmyms       )
  call free_ptr( trans%nasm0       )
  call free_ptr( trans%numpp       )
  call free_ptr( trans%npossp      )
  call free_ptr( trans%nptrms      )
  call free_ptr( trans%nallms      )
  call free_ptr( trans%ndim0g      )
  call free_ptr( trans%nvalue      )
  call free_ptr( trans%n_regions   )
  call free_ptr( trans%nfrstlat    )
  call free_ptr( trans%nlstlat     )
  call free_ptr( trans%nptrlat     )
  call free_ptr( trans%nptrfrstlat )
  call free_ptr( trans%nptrlstlat  )
  call free_ptr( trans%nsta        )
  call free_ptr( trans%nonl        )
  call free_ptr( trans%ldsplitlat  )
  call free_ptr( trans%nultpp      )
  call free_ptr( trans%nptrls      )
  call free_ptr( trans%nnmeng      )
  call free_ptr( trans%rmu         )
  call free_ptr( trans%rgw         )
  call free_ptr( trans%rpnm        )
  call free_ptr( trans%npms        )
  call free_ptr( trans%rlapin      )
  call free_ptr( trans%ndglu       )
  call trans_release( trans%handle )
  trans_delete = TRANS_SUCCESS
end function trans_delete

function trans_finalize() bind(C,name="trans_finalize")
  use, intrinsic :: iso_c_binding
  integer(c_int) :: trans_finalize
  call TRANS_END()
  if( USE_MPI ) call MPL_END(LDMEMINFO=.FALSE.)
  if( close_devnull ) then
    ! Don't close devnull in case other code is also using this unit
    ! close (devnull())
  endif

  is_init = .False.

  trans_finalize = TRANS_SUCCESS
end function trans_finalize


function get_nlon( trans ) result(nlon)
  use, intrinsic :: iso_c_binding, only : c_associated, c_f_pointer
  integer ::nlon
  type(Trans_t) :: trans
  integer, pointer :: nloen(:)
  nlon = trans%nlon
  if( nlon < 0 ) then
    if( c_associated( trans%nloen ) ) then
      call c_f_pointer( trans%nloen, nloen, (/trans%ndgl/) )
      nlon = nloen(1)
    else
      nlon = 2*trans%ndgl
    endif
  endif
end function

function assert_global(trans,RGP) result(iret)
  integer :: iret
  type(Trans_t), intent(in) :: trans
  real(c_double), intent(in) :: RGP(:,:,:)  !(NPROMA==ngptotg,NFLD,NGPBLKS==1)
  integer :: nproma, ngpblks, nlon
  iret = TRANS_SUCCESS
  nproma   = size(RGP,1)
  ngpblks  = size(RGP,3)

  if( trans%nproc /= 1 ) then
    call transi_error("trans_invtrans: ERROR: Configuration only valid for nproc == 1")
    iret = TRANS_ERROR
    return
  endif

  if( trans%llatlon == 1 ) then
    nlon     = get_nlon(trans)
    if( trans%ngptot /= trans%ngptotg + nlon ) then
      call transi_error("trans: Assertion failed for lonlat grids: (ngptot == ngptotg+nlon)")
      iret = TRANS_ERROR
      return
    endif
  endif

  if( nproma  /= trans%ngptotg ) then
    call transi_error("trans_invtrans: ERROR: Configuration only valid for nproma == ngpgot")
    iret = TRANS_ERROR
    return
  endif

  if( ngpblks /= 1 ) then
    call transi_error("trans: ERROR: Configuration only valid for ngpblks == 1")
    iret = TRANS_ERROR
    return
  endif
end function

function prepare_global_invtrans(trans,RGP,RGPM) result(iret)
  integer :: iret
  type(Trans_t), intent(in) :: trans
  real(c_double), target, intent(in) :: RGP(:,:,:)      !(NPROMA==ngptotg,NFLD,NGPBLKS==1)
  real(c_double), pointer, intent(out) :: RGPM(:,:,:)   !(NPROMA==ngptot, NFLD,NGPBLKS==1)
    !! Modified RGP to add one duplicate latitude at equator
  integer :: ilat, ilon, nfld

  iret = assert_global(trans,RGP)
  if( iret /= TRANS_SUCCESS ) return

  if( trans%llatlon == 1 ) then
    nfld     = size(RGP,2)
    allocate( RGPM(trans%ngptot,nfld,1) )
  else
    RGPM => RGP
  endif
end function

function finish_global_invtrans(trans,RGP,RGPM) result(iret)
  integer :: iret
  type(Trans_t), intent(in) :: trans
  real(c_double), intent(inout) :: RGP(:,:,:)               !(NPROMA==ngptotg,FIELD,NGPBLKS==1)
  real(c_double), pointer, intent(inout) :: RGPM(:,:,:)     !(NPROMA==ngptotg,FIELD,NGPBLKS==1)
    !! Modified RGP with an added duplicate latitude at equator
  integer :: nlon, ilat, ilon, icount

  iret = assert_global(trans,RGP)
  if( iret /= TRANS_SUCCESS ) return

  if( trans%llatlon == 1 ) then
    nlon = get_nlon(trans)
    icount = 0
    do ilat=1,trans%ndgl+2
      if( ilat <= trans%ndgl/2 .or. ilat >= trans%ndgl/2+2) then
        do ilon=1,nlon
          icount=icount+1
          RGP(icount,:,1) = RGPM(ilon+(ilat-1)*nlon,:,1)
        enddo ! ilon
      endif
    enddo ! ilat
    deallocate(RGPM)
    nullify(RGPM)
  else
    nullify(RGPM)
  endif
end function

function prepare_global_gptosp_trans(trans,RGP,RGPM) result(iret)
  integer :: iret
  type(Trans_t), intent(in) :: trans
  real(c_double), target, intent(in) :: RGP(:,:,:)      !(NPROMA==ngptotg,NFLD,NGPBLKS==1)
  real(c_double), pointer, intent(out) :: RGPM(:,:,:)   !(NPROMA==ngptot, NFLD,NGPBLKS==1)
    !! Modified RGP to add one duplicate latitude at equator
  integer :: nlon, ilat, ilon, icount, nfld

  iret = assert_global(trans,RGP)
  if( iret /= TRANS_SUCCESS ) return

  if( trans%llatlon == 1 ) then
    nfld     = size(RGP,2)
    nlon     = get_nlon(trans)
    icount = 0
    allocate( RGPM(trans%ngptot,nfld,1) )
    do ilat=1,trans%ndgl+2 ! There is 1 too little latitude in RGPM
      do ilon=1,nlon
        icount = icount+1
        RGPM(ilon+(ilat-1)*nlon,:,1) = RGP(icount,:,1)
      enddo ! ilon
      if( ilat == trans%ndgl/2+1) then
        icount = icount-nlon
      endif
    enddo ! ilat
  else
    RGPM => RGP
  endif
end function

function finish_global_gptosp_trans(trans,RGP,RGPM) result(iret)
  integer :: iret
  type(Trans_t), intent(in) :: trans
  real(c_double), intent(inout) :: RGP(:,:,:)               !(NPROMA==ngptotg,FIELD,NGPBLKS==1)
  real(c_double), pointer, intent(inout) :: RGPM(:,:,:)     !(NPROMA==ngptotg,FIELD,NGPBLKS==1)
    !! Modified RGP with an added duplicate latitude at equator

  iret = assert_global(trans,RGP)
  if( iret /= TRANS_SUCCESS ) return

  if( trans%llatlon == 1 ) then
    deallocate(RGPM)
    nullify(RGPM)
  else
    nullify(RGPM)
  endif
end function

function trans_dirtrans(args) bind(C,name="trans_dirtrans") result(iret)
  use, intrinsic :: iso_c_binding
  integer(c_int) :: iret
  type(DirTrans_t), intent(inout) :: args
  real(c_double), pointer :: RSPVOR(:,:)    !(FIELD,WAVE)
  real(c_double), pointer :: RSPDIV(:,:)    !(FIELD,WAVE)
  real(c_double), pointer :: RSPSCALAR(:,:) !(FIELD,WAVE)
  real(c_double), pointer :: RGP(:,:,:)     !(NPROMA,IF_GP,NGPBLKS)
  real(c_double), pointer :: RGPM(:,:,:)    !(NPROMA,FIELD,NGPBLKS)
  type(Trans_t), pointer :: trans
  logical :: llatlon

  if( args%count > 0 ) then
    call transi_error("trans_dirtrans: ERROR: arguments are not new")
    iret = TRANS_STALE_ARG
    return
  endif
  args%count = 1

  if( .not. c_associated(args%trans) ) then
    call transi_error( "trans_dirtrans: ERROR: trans was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%trans, trans )

  if( .not. c_associated(args%rgp) ) then
    call transi_error( "trans_dirtrans: ERROR: Array RGP was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif

  if( args%nvordiv > 0 ) then
    if( .not. c_associated(args%rspvor)    ) then
      call transi_error( "Array RSPVOR was not allocated" )
      iret = TRANS_MISSING_ARG
      return
    endif
    if( .not. c_associated(args%rspdiv)    ) then
      call transi_error( "Array RSPDIV was not allocated" )
      iret = TRANS_MISSING_ARG
      return
    endif
    call C_F_POINTER(args%rspvor,    RSPVOR,    (/args%nvordiv, trans%nspec2/)  )
    call C_F_POINTER(args%rspdiv,    RSPDIV,    (/args%nvordiv, trans%nspec2/)  )
  endif
  if( args%nscalar > 0 ) then
    if( .not. c_associated(args%rspscalar) ) then
      call transi_error( "Array RSPSCALAR was not allocated" )
      iret = TRANS_MISSING_ARG
      return
    endif
    call C_F_POINTER(args%rspscalar, RSPSCALAR, (/args%nscalar, trans%nspec2/)  )
  endif

  llatlon = .false.
  if( trans%llatlon /= 0 ) llatlon = .true.

  if( args%lglobal == 1 ) then
    call C_F_POINTER( args%rgp, RGP, (/trans%ngptotg,args%nscalar+2*args%nvordiv,1/) )
    iret = prepare_global_gptosp_trans(trans,RGP,RGPM)
    if( iret /= TRANS_SUCCESS ) return
  else
    call C_F_POINTER( args%rgp, RGP, (/args%nproma,args%nscalar+2*args%nvordiv,args%ngpblks/) )
    RGPM => RGP
  endif

  if( args%nvordiv > 0 .and. args%nscalar > 0 ) then
    call DIR_TRANS( KRESOL=trans%handle, &
      &             KPROMA=args%nproma, &
      &             LDLATLON=llatlon, &
      &             PGP=RGPM, &
      &             PSPVOR=RSPVOR,PSPDIV=RSPDIV,PSPSCALAR=RSPSCALAR ) ! unused args: KVSETUV,KVSETSC
  elseif( args%nscalar > 0 ) then
    call DIR_TRANS( KRESOL=trans%handle, &
      &             KPROMA=args%nproma, &
      &             LDLATLON=llatlon, &
      &             PGP=RGPM, &
      &             PSPSCALAR=RSPSCALAR ) ! unused args: KVSETUV,KVSETSC
  elseif( args%nvordiv > 0 ) then
    call DIR_TRANS( KRESOL=trans%handle, &
      &             KPROMA=args%nproma, &
      &             LDLATLON=llatlon, &
      &             PGP=RGPM, &
      &             PSPVOR=RSPVOR,PSPDIV=RSPDIV ) ! unused args: KVSETUV,KVSETSC
  endif

  if( args%lglobal == 1 ) then
    iret = finish_global_gptosp_trans(trans,RGP,RGPM)
    if( iret /= TRANS_SUCCESS ) return
  else
    nullify(RGPM)
  endif

  iret = TRANS_SUCCESS
end function trans_dirtrans

function trans_invtrans(args) bind(C,name="trans_invtrans") result(iret)
  use, intrinsic :: iso_c_binding
  integer(c_int) :: iret
  type(InvTrans_t), intent(inout) :: args
  real(c_double), pointer :: RSPVOR(:,:)    !(FIELD,WAVE)
  real(c_double), pointer :: RSPDIV(:,:)    !(FIELD,WAVE)
  real(c_double), pointer :: RSPSCALAR(:,:) !(FIELD,WAVE)
  real(c_double), pointer :: RGP(:,:,:)     !(NPROMA,FIELD,NGPBLKS)
  real(c_double), pointer :: RGPM(:,:,:)    !(NPROMA,FIELD,NGPBLKS)

  logical :: lscalarders
  logical :: luvder_EW
  logical :: lvordivgp
  logical :: llatlon
  type(Trans_t), pointer :: trans
  integer :: nfld_gp

  if( args%count > 0 ) then
    call transi_error( "trans_invtrans: ERROR: arguments are not new" )
    iret = TRANS_STALE_ARG
    return
  endif
  args%count = 1

  if( .not. c_associated(args%trans) ) then
    call transi_error( "trans was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%trans, trans )

  if( .not. c_associated(args%rgp) ) then
    call transi_error( "Array RGP was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif

  lscalarders = .false.; if( args%lscalarders == 1 ) lscalarders = .true.
  luvder_EW   = .false.; if( args%luvder_EW   == 1 ) luvder_EW   = .true.
  lvordivgp   = .false.; if( args%lvordivgp   == 1 ) lvordivgp   = .true.

  if( args%nvordiv > 0 ) then
    if( .not. c_associated(args%rspvor)    ) then
      call transi_error( "Array RSPVOR was not allocated" )
      iret = TRANS_MISSING_ARG
      return
    endif
    if( .not. c_associated(args%rspdiv)    ) then
      call transi_error( "Array RSPDIV was not allocated" )
      iret = TRANS_MISSING_ARG
      return
    endif
    call C_F_POINTER(args%rspvor,    RSPVOR,    (/args%nvordiv, trans%nspec2/)  )
    call C_F_POINTER(args%rspdiv,    RSPDIV,    (/args%nvordiv, trans%nspec2/)  )
  endif
  if( args%nscalar > 0 ) then
    if( .not. c_associated(args%rspscalar) ) then
      call transi_error( "Array RSPSCALAR was not allocated" )
      iret = TRANS_MISSING_ARG
      return
    endif
    call C_F_POINTER(args%rspscalar, RSPSCALAR, (/args%nscalar, trans%nspec2/)  )
  endif

  llatlon = .false.
  if( trans%llatlon /= 0 ) llatlon = .true.

  nfld_gp = 0
  if( lvordivgp )    nfld_gp = nfld_gp + 2*args%nvordiv  ! voriticty + divergence
  if( .true. )       nfld_gp = nfld_gp + args%nscalar    ! scalars
  if( .true. )       nfld_gp = nfld_gp + 2*args%nvordiv  ! u + v
  if( lscalarders )  nfld_gp = nfld_gp + args%nscalar    ! scalars N-S derivatives
  if( luvder_EW )    nfld_gp = nfld_gp + 2*args%nvordiv  ! u + v   E-W derivatives
  if( lscalarders )  nfld_gp = nfld_gp + args%nscalar    ! scalars E-W derivatives

  if( args%lglobal == 1 ) then
    call C_F_POINTER( args%rgp, RGP, (/trans%ngptotg,nfld_gp,1/) )
    iret = prepare_global_invtrans(trans,RGP,RGPM)
    if( iret /= TRANS_SUCCESS ) return
  else
    call C_F_POINTER( args%rgp, RGP, (/args%nproma,nfld_gp,args%ngpblks/) )
    RGPM => RGP
  endif

  if( args%nvordiv > 0 .and. args%nscalar > 0 ) then
    call INV_TRANS( KRESOL=trans%handle, &
      &             KPROMA=args%nproma, &
      &             LDLATLON=llatlon, &
      &             LDSCDERS=lscalarders, &
      &             LDVORGP=lvordivgp, &
      &             LDDIVGP=lvordivgp, &
      &             LDUVDER=luvder_EW, &
      &             PSPVOR=RSPVOR,PSPDIV=RSPDIV,PSPSCALAR=RSPSCALAR, &
      &             PGP=RGPM ) ! unused args: KVSETUV,KVSETSC
  elseif( args%nscalar > 0 ) then
    call INV_TRANS( KRESOL=trans%handle, &
      &             KPROMA=args%nproma, &
      &             LDLATLON=llatlon, &
      &             LDSCDERS=lscalarders, &
      &             PSPSCALAR=RSPSCALAR, &
      &             PGP=RGPM ) ! unused args: KVSETUV,KVSETSC
  elseif( args%nvordiv > 0 ) then
    call INV_TRANS( KRESOL=trans%handle, &
      &             KPROMA=args%nproma, &
      &             LDLATLON=llatlon, &
      &             LDVORGP=lvordivgp, &
      &             LDDIVGP=lvordivgp, &
      &             LDUVDER=luvder_EW, &
      &             PSPVOR=RSPVOR,PSPDIV=RSPDIV, &
      &             PGP=RGPM ) ! unused args: KVSETUV,KVSETSC
  endif

  if( args%lglobal == 1 ) then
    iret = finish_global_invtrans(trans,RGP,RGPM)
    if( iret /= TRANS_SUCCESS ) return
  else
    nullify(RGPM)
  endif


  iret = TRANS_SUCCESS

end function trans_invtrans

function trans_invtrans_adj(args) bind(C,name="trans_invtrans_adj") result(iret)
  use, intrinsic :: iso_c_binding
  integer(c_int) :: iret
  type(InvTransAdj_t), intent(inout) :: args
  real(c_double), pointer :: RSPVOR(:,:)    !(FIELD,WAVE)
  real(c_double), pointer :: RSPDIV(:,:)    !(FIELD,WAVE)
  real(c_double), pointer :: RSPSCALAR(:,:) !(FIELD,WAVE)
  real(c_double), pointer :: RGP(:,:,:)     !(NPROMA,FIELD,NGPBLKS)
  real(c_double), pointer :: RGPM(:,:,:)    !(NPROMA,FIELD,NGPBLKS)

  logical :: lscalarders
  logical :: luvder_EW
  logical :: lvordivgp
  logical :: llatlon
  type(Trans_t), pointer :: trans
  integer :: nfld_gp

  if( args%count > 0 ) then
    call transi_error( "trans_invtrans_adj: ERROR: arguments are not new" )
    iret = TRANS_STALE_ARG
    return
  endif
  args%count = 1

  if( .not. c_associated(args%trans) ) then
    call transi_error( "trans_invtrans_adj:trans was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%trans, trans )

  if( .not. c_associated(args%rgp) ) then
    call transi_error( "trans_invtrans_adj:Array RGP was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif

  lscalarders = .false.; if( args%lscalarders == 1 ) lscalarders = .true.
  luvder_EW   = .false.; if( args%luvder_EW   == 1 ) luvder_EW   = .true.
  lvordivgp   = .false.; if( args%lvordivgp   == 1 ) lvordivgp   = .true.

  if( args%nvordiv > 0 ) then
    if( .not. c_associated(args%rspvor)    ) then
      call transi_error( "trans_invtrans_adj::Array RSPVOR was not allocated" )
      iret = TRANS_MISSING_ARG
      return
    endif
    if( .not. c_associated(args%rspdiv)    ) then
      call transi_error( "trans_invtrans_adj::Array RSPDIV was not allocated" )
      iret = TRANS_MISSING_ARG
      return
    endif
    call C_F_POINTER(args%rspvor,    RSPVOR,    (/args%nvordiv, trans%nspec2/)  )
    call C_F_POINTER(args%rspdiv,    RSPDIV,    (/args%nvordiv, trans%nspec2/)  )
  endif
  if( args%nscalar > 0 ) then
    if( .not. c_associated(args%rspscalar) ) then
      call transi_error( "trans_invtrans_adj::Array RSPSCALAR was not allocated" )
      iret = TRANS_MISSING_ARG
      return
    endif
    call C_F_POINTER(args%rspscalar, RSPSCALAR, (/args%nscalar, trans%nspec2/)  )
  endif

  llatlon = .false.
  if( trans%llatlon /= 0 ) llatlon = .true.

  nfld_gp = 0
  if( lvordivgp )    nfld_gp = nfld_gp + 2*args%nvordiv  ! voriticty + divergence
  if( .true. )       nfld_gp = nfld_gp + args%nscalar    ! scalars
  if( .true. )       nfld_gp = nfld_gp + 2*args%nvordiv  ! u + v
  if( lscalarders )  nfld_gp = nfld_gp + args%nscalar    ! scalars N-S derivatives
  if( luvder_EW )    nfld_gp = nfld_gp + 2*args%nvordiv  ! u + v   E-W derivatives
  if( lscalarders )  nfld_gp = nfld_gp + args%nscalar    ! scalars E-W derivatives

  if( args%lglobal == 1 ) then
    call C_F_POINTER( args%rgp, RGP, (/trans%ngptotg,nfld_gp,1/) )
    iret = prepare_global_gptosp_trans(trans,RGP,RGPM)
    if( iret /= TRANS_SUCCESS ) return
  else
    call C_F_POINTER( args%rgp, RGP, (/args%nproma,nfld_gp,args%ngpblks/) )
    RGPM => RGP
  endif

  ! Note that llatlon is not an option in INV_TRANSAD unlile INV_TRANS and DIR_TRANS
  if( args%nvordiv > 0 .and. args%nscalar > 0 ) then
    call INV_TRANSAD( KRESOL=trans%handle, &
      &               KPROMA=args%nproma, &
!      &               LDLATLON=llatlon, &
      &               LDSCDERS=lscalarders, &
      &               LDVORGP=lvordivgp, &
      &               LDDIVGP=lvordivgp, &
      &               LDUVDER=luvder_EW, &
      &               PSPVOR=RSPVOR,PSPDIV=RSPDIV,PSPSCALAR=RSPSCALAR, &
      &               PGP=RGPM ) ! unused args: KVSETUV,KVSETSC
  elseif( args%nscalar > 0 ) then
    call INV_TRANSAD( KRESOL=trans%handle, &
      &               KPROMA=args%nproma, &
!      &               LDLATLON=llatlon, &
      &               LDSCDERS=lscalarders, &
      &               PSPSCALAR=RSPSCALAR, &
      &               PGP=RGPM ) ! unused args: KVSETUV,KVSETSC
  elseif( args%nvordiv > 0 ) then
    call INV_TRANSAD( KRESOL=trans%handle, &
      &               KPROMA=args%nproma, &
!      &               LDLATLON=llatlon, &
      &               LDVORGP=lvordivgp, &
      &               LDDIVGP=lvordivgp, &
      &               LDUVDER=luvder_EW, &
      &               PSPVOR=RSPVOR,PSPDIV=RSPDIV, &
      &               PGP=RGPM ) ! unused args: KVSETUV,KVSETSC
  endif

  if( args%lglobal == 1 ) then
    iret = finish_global_gptosp_trans(trans,RGP,RGPM)
    if( iret /= TRANS_SUCCESS ) return
  else
    nullify(RGPM)
  endif

  iret = TRANS_SUCCESS
end function trans_invtrans_adj

function trans_distgrid(args) bind(C,name="trans_distgrid") result(iret)
  use, intrinsic :: iso_c_binding
  integer(c_int) :: iret
  type(DistGrid_t), intent(inout) :: args
  real(c_double), pointer :: RGPG(:,:)     ! (NFLD_from,NGPTOTG  (+nlon) )  (+nlon in case of LonLat grid)
  real(c_double), pointer :: LL_RGPG(:,:)  ! (NFLD_from,NGPTOTG          )
  real(c_double), pointer :: RGP (:,:,:)   ! (NPROMA,IF_GP,NGPBLKS)
  integer(c_int), pointer :: NFROM(:)
  type(Trans_t), pointer  :: trans
  integer :: jfld, isend, jsend
  integer :: icount, ilat, ilon, nlon
  integer :: check
  integer(c_int), pointer :: nloen(:)

  if( args%count > 0 ) then
    call transi_error( "trans_distgrid: ERROR: arguments are not new" )
    iret = TRANS_STALE_ARG
    return
  endif
  args%count = 1

  if( .not. c_associated(args%trans) )   then
    call transi_error( "trans was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%trans, trans )

  if( .not. c_associated(args%nfrom) ) then
    call transi_error( "Array NFROM was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%nfrom, NFROM, (/args%nfld/) )


  isend = 0
  do jfld = 1, args%nfld
    if ( NFROM(jfld) == trans%myproc ) isend = isend + 1
  enddo

  if( .not. c_associated(args%rgp) )  then
    call transi_error( "Array RGP was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%rgp, RGP, (/args%nproma,args%nfld,args%ngpblks/) )

  if( isend > 0 ) then
    if( .not. c_associated(args%rgpg) ) then
      call transi_error( "Array RGPG was not allocated" )
      iret = TRANS_MISSING_ARG
      return
    endif


    if( trans%llatlon == 1 ) then

      nlon = trans%nlon
      if( nlon < 0 ) then
        if( c_associated( trans%nloen ) ) then
          call C_F_POINTER( trans%nloen, nloen, (/trans%ndgl/) )
          nlon = nloen(1)
        else
          nlon = 2*trans%ndgl
        endif
      endif

      call C_F_POINTER( args%rgpg, LL_RGPG, (/trans%ngptotg,isend/) )
      allocate( RGPG(trans%ngptotg+nlon,isend) ) ! 1 extra latitudes should be allocated
      do jsend=1,isend
        check = 0
        icount = 0
        do ilat=1,trans%ndgl+2 ! There is 1 too little latitude in LL_RGPG
          do ilon=1,nlon
            ICOUNT=ICOUNT+1
            RGPG(ILON+(ILAT-1)*nlon,jsend) = LL_RGPG(ICOUNT,jsend)
            check = check+1
          enddo ! ilon
          if( ilat == trans%ndgl/2+1) then
            ICOUNT = ICOUNT-nlon
          endif
        enddo ! ilat
        if( check /= trans%ngptotg+nlon ) then
          call transi_error( "ERROR: not all values are assigned" )
          iret = TRANS_ERROR
          deallocate( RGPG )
          return
        endif
      enddo ! jsend
    else
      call C_F_POINTER( args%rgpg, RGPG, (/trans%ngptotg,isend/) )
    endif ! llatlon

    call DIST_GRID(PGPG=RGPG,KFDISTG=args%nfld,KFROM=NFROM,KPROMA=args%nproma,KRESOL=trans%handle,PGP=RGP)

    if( trans%llatlon == 1 ) then
      deallocate( RGPG )
    endif
  else
    call DIST_GRID(          KFDISTG=args%nfld,KFROM=NFROM,KPROMA=args%nproma,KRESOL=trans%handle,PGP=RGP)
  endif

  iret = TRANS_SUCCESS

end function trans_distgrid


function trans_gathgrid(args) bind(C,name="trans_gathgrid") result(iret)
  use, intrinsic :: iso_c_binding
  integer(c_int) :: iret
  type(GathGrid_t), intent(inout) :: args
  real(c_double), pointer :: RGPG(:,:)     !(NFLD_to,NGPTOTG)
  real(c_double), pointer :: RGP (:,:,:)   !(NPROMA,NFLD,NGPBLKS)
  integer(c_int), pointer :: NTO(:)
  type(Trans_t), pointer  :: trans
  real(c_double), pointer :: LL_RGPG (:,:)
  integer :: jfld, irecv
  integer :: icount, ilat, ilon, jrecv, nlon
  integer(c_int), pointer :: nloen(:)

  if( args%count > 0 ) then
    call transi_error( "trans_gathgrid: ERROR: arguments are not new" )
    iret = TRANS_STALE_ARG
    return
  endif
  args%count = 1

  if( .not. c_associated(args%trans) ) then
    call transi_error( "trans was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%trans, trans )

  if( .not. c_associated(args%nto) ) then
    call transi_error( "trans_gath_grid: Array NTO was not allocated")
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%nto, NTO, (/args%nfld/) )


  irecv = 0
  do jfld = 1, args%nfld
    if ( NTO(jfld) == trans%myproc ) irecv = irecv + 1
  enddo

  if( .not. c_associated(args%rgp) ) then
    call transi_error( "trans_gath_grid: Array RGP was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%rgp, RGP, (/args%nproma,args%nfld,args%ngpblks/) )


  if( irecv > 0 ) then
    if( .not. c_associated(args%rgpg) ) then
      call transi_error( "Array RGPG was not allocated" )
      iret = TRANS_MISSING_ARG
      return
    endif
    if( trans%llatlon == 1 ) then
      nlon = trans%nlon
      if( nlon < 0 ) then
        if( c_associated( trans%nloen ) ) then
          call C_F_POINTER( trans%nloen, nloen, (/trans%ndgl/) )
          nlon = nloen(1)
        else
          nlon = 2*trans%ndgl
        endif
      endif
      allocate( RGPG(trans%ngptotg+nlon,irecv) ) ! 1 extra latitudes
    else
      call C_F_POINTER( args%rgpg, RGPG, (/trans%ngptotg,irecv/) )
    endif
    call GATH_GRID(KRESOL=trans%handle,KFGATHG=args%nfld,KTO=NTO,KPROMA=args%nproma,PGP=RGP,PGPG=RGPG)

    if( trans%llatlon == 1 ) then

      ! There is 1 too many latitude in RGPG
      call C_F_POINTER( args%rgpg, LL_RGPG, (/trans%ngptotg,irecv/) )

      do jrecv=1,irecv
        icount = 0
        do ilat=1,trans%ndgl+2
          do ilon=1,nlon
            if( ilat <= trans%ndgl/2 .or. ilat >= trans%ndgl/2+2) then
              ICOUNT=ICOUNT+1
              LL_RGPG(ICOUNT,jrecv) = RGPG(ILON+(ILAT-1)*nlon,jrecv)
            endif
          enddo ! ilon
        enddo ! ilat
        if( ICOUNT /= trans%ngptotg) then
          call transi_error( "CHECK failed" )
          iret = TRANS_ERROR
          deallocate( RGPG )
          return
        endif
      enddo ! jrecv
      deallocate( RGPG )
    endif
  else
    call GATH_GRID(KRESOL=trans%handle,KFGATHG=args%nfld,KTO=NTO,KPROMA=args%nproma,PGP=RGP)
  endif

  iret = TRANS_SUCCESS

end function trans_gathgrid


function trans_distspec(args) bind(C,name="trans_distspec") result(iret)
  use, intrinsic :: iso_c_binding
  integer(c_int) :: iret
  type(DistSpec_t), intent(inout) :: args
  real(c_double), pointer :: RSPEC (:,:)   ! (NFLD,NSPEC2)
  real(c_double), pointer :: RSPECG(:,:)   ! (NFLD_from,NSPEC2G)
  integer(c_int), pointer :: NFROM(:)
  type(Trans_t), pointer  :: trans
  integer :: jfld, isend

  if( args%count > 0 ) then
    call transi_error( "trans_distspec: ERROR: arguments are not new" )
    iret = TRANS_STALE_ARG
    return
  endif
  args%count = 1

  if( .not. c_associated(args%trans) ) then
    call transi_error( "trans was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%trans, trans )

  if( .not. c_associated(args%nfrom) ) then
    call transi_error( "Array NFROM was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%nfrom, NFROM, (/args%nfld/) )


  isend = 0
  do jfld = 1, args%nfld
    if ( NFROM(jfld) == trans%myproc ) isend = isend + 1
  enddo

  if( .not. c_associated(args%rspec) ) then
    call transi_error( "Array RSPEC was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%rspec, RSPEC, (/args%nfld,trans%nspec2/) )

  if( isend > 0 ) then
    if( .not. c_associated(args%rspecg) ) then
      call transi_error(  "Array RSPECG was not allocated" )
    endif
    call C_F_POINTER( args%rspecg, RSPECG, (/isend,trans%nspec2g/) )
    call DIST_SPEC(KRESOL=trans%handle,KFDISTG=args%nfld,KFROM=NFROM,PSPEC=RSPEC,PSPECG=RSPECG)
  else
    call DIST_SPEC(KRESOL=trans%handle,KFDISTG=args%nfld,KFROM=NFROM,PSPEC=RSPEC)
  endif

  iret = TRANS_SUCCESS
end function trans_distspec



function trans_gathspec(args) bind(C,name="trans_gathspec") result(iret)
  use, intrinsic :: iso_c_binding
  integer(c_int) :: iret
  type(GathSpec_t), intent(inout) :: args
  real(c_double), pointer :: RSPEC(:,:)    ! (NFLD,NSPEC2)
  real(c_double), pointer :: RSPECG(:,:)   ! (NFLD_to,NSPEC2G)
  integer(c_int), pointer :: NTO(:)
  type(Trans_t), pointer  :: trans
  integer :: jfld, irecv

  if( args%count > 0 ) then
    call transi_error( "trans_gathspec: ERROR: arguments are not new" )
    iret = TRANS_STALE_ARG
    return
  endif
  args%count = 1

  if( .not. c_associated(args%trans) ) then
    call transi_error( "trans was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%trans, trans )

  if( .not. c_associated(args%nto) ) then
    call transi_error( "Array NTO was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%nto, NTO, (/args%nfld/) )


  irecv = 0
  do jfld = 1, args%nfld
    if ( NTO(jfld) == trans%myproc ) irecv = irecv + 1
  enddo

  if( .not. c_associated(args%rspec) ) then
    call transi_error( "Array RSPEC was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%rspec, RSPEC, (/args%nfld,trans%nspec2/) )

  if( irecv > 0 ) then
    if( .not. c_associated(args%rspecg) ) then
      call transi_error( "Array RSPECG was not allocated" )
      iret = TRANS_MISSING_ARG
      return
    endif
    call C_F_POINTER( args%rspecg, RSPECG, (/irecv,trans%nspec2g/) )
    call GATH_SPEC(KRESOL=trans%handle,KFGATHG=args%nfld,KTO=NTO,PSPEC=RSPEC,PSPECG=RSPECG)
  else
    call GATH_SPEC(KRESOL=trans%handle,KFGATHG=args%nfld,KTO=NTO,PSPEC=RSPEC)
  endif

  iret = TRANS_SUCCESS

end function trans_gathspec


function trans_vordiv_to_UV(args) bind(C,name="trans_vordiv_to_UV") result(iret)
  use, intrinsic :: iso_c_binding
  integer(c_int) :: iret
  type(VorDivToUV_t), intent(inout) :: args
  real(c_double), pointer :: RSPVOR(:,:)
  real(c_double), pointer :: RSPDIV(:,:)
  real(c_double), pointer :: RSPU(:,:)
  real(c_double), pointer :: RSPV(:,:)
  integer(c_int) :: err

  if( args%count > 0 ) then
    call transi_error( "trans_vordiv_to_UV: ERROR: arguments are not new" )
    iret = TRANS_STALE_ARG
    return
  endif
  args%count = 1

  if( args%ncoeff == 0 ) then
    call transi_error( "trans_vordiv_to_UV: ERROR: missing argument nspec2")
    iret = TRANS_MISSING_ARG
    return
  endif

  if( args%nsmax == 0 ) then
    call transi_error( "trans_vordiv_to_UV: ERROR: missing argument nsmax")
    iret = TRANS_MISSING_ARG
    return
  endif


  ! Set vorticity
  if( .not. c_associated(args%rspvor) ) then
    call transi_error( "Array RSPVOR was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%rspvor, RSPVOR, (/args%nfld,args%ncoeff/) )

  ! Set divergence
  if( .not. c_associated(args%rspdiv) ) then
    call transi_error( "Array RSPDIV was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%rspdiv, RSPDIV, (/args%nfld,args%ncoeff/) )

  ! Set U
  if( .not. c_associated(args%rspu) ) then
    call transi_error( "Array RSPU was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%rspu, RSPU, (/args%nfld,args%ncoeff/) )

  ! Set V
  if( .not. c_associated(args%rspv) ) then
    call transi_error( "Array RSPV was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%rspv, RSPV, (/args%nfld,args%ncoeff/) )


  if ( .not. is_init ) then
    err = trans_init()
  endif

  ! PSPVOR(:,:) - spectral vorticity (input)
  ! PSPDIV(:,:) - spectral divergence (input)
  ! PSPU(:,:)   - spectral U (u*cos(theta) (output)
  ! PSPV(:,:)   - spectral V (v*cos(theta) (output)
  ! KSMAX       - spectral resolution (input)
  call VORDIV_TO_UV(PSPVOR=RSPVOR,PSPDIV=RSPDIV,PSPU=RSPU,PSPV=RSPV,KSMAX=args%nsmax)
  iret = TRANS_SUCCESS

end function trans_vordiv_to_UV

function trans_specnorm(args) bind(C,name="trans_specnorm") result(iret)
  use, intrinsic :: iso_c_binding
  integer(c_int) :: iret
  type(SpecNorm_t), intent(inout) :: args
  real(c_double), pointer :: RSPEC(:,:)     !(IF_GP,NGPTOTG)
  real(c_double), pointer :: RNORM(:)
  real(c_double), pointer :: RMET(:)
  type(Trans_t), pointer  :: trans
  integer :: jfld, irecv

  if( args%count > 0 ) then
    call transi_error( "trans_specnorm: ERROR: arguments are not new" )
    iret = TRANS_STALE_ARG
    return
  endif
  args%count = 1

  if( .not. c_associated(args%trans) ) then
    call transi_error( "trans was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%trans, trans )

  if( args%nfld == 0 ) then
    iret = TRANS_SUCCESS
    return
  endif

  if( .not. c_associated(args%rspec) ) then
    call transi_error( "Array RSPEC was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%rspec, RSPEC, (/args%nfld,trans%nspec2/) )

  if( .not. c_associated(args%rnorm) ) then
    call transi_error( "Array RNORM was not allocated" )
    iret = TRANS_MISSING_ARG
    return
  endif
  call C_F_POINTER( args%rnorm, RNORM, (/args%nfld/) )

  if( c_associated(args%rmet) ) then
    call C_F_POINTER( args%rmet, RMET, (/trans%nsmax+1/) )
    RMET(0:) => RMET(:)
  endif

  if( .not. c_associated(args%rmet) ) then
    call SPECNORM(KRESOL=trans%handle,PSPEC=RSPEC,KMASTER=args%nmaster,PNORM=RNORM)
  else
    call SPECNORM(KRESOL=trans%handle,PSPEC=RSPEC,KMASTER=args%nmaster,PNORM=RNORM,PMET=RMET)
  endif

  iret = TRANS_SUCCESS

end function trans_specnorm


end module trans_module
