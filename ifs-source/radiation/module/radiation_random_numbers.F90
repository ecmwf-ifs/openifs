! radiation_random_numbers.F90 - Generate random numbers for McICA solver
!
! (C) Copyright 2020- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! The derived type "rng_type" is a random number generator that uses
! either (1) Fortran's built-in random_number function, or (2) a
! vectorized version of the MINSTD linear congruential generator.  In
! the case of (2), an rng_type object is initialized with a seed that
! is used to fill up a state of "nmaxstreams" elements using the C++
! minstd_rand0 version of the MINSTD linear congruential generator
! (LNG), which has the form istate[i+1] = mod(istate[i]*A0, M) from
! i=1 to i=nmaxstreams. Subsequent requests for blocks of nmaxstreams
! of random numbers use the C++ minstd_ran algorithm in a vectorizable
! form, which modifies the state elements via istate[i] <-
! mod(istate[i]*A, M). Uniform deviates are returned that normalize
! the state elements to the range 0-1.
!
! The MINSTD generator was coded because the random_numbers_mix
! generator in the IFS was found not to vectorize well on some
! hardware.  I am no expert on random number generators, so my
! implementation should really be looked at and improved by someone
! who knows what they are doing.
!
! Reference for MINSTD: Park, Stephen K.; Miller, Keith
! W. (1988). "Random Number Generators: Good Ones Are Hard To Find"
! (PDF). Communications of the ACM. 31 (10):
! 1192-1201. doi:10.1145/63039.63042
!
! Modifications
!   2022-12-01  R. Hogan  Fixed zeroed state in single precision

module radiation_random_numbers

  use parkind1, only : jprb, jprd, jpim, jpib

  implicit none

  public :: rng_type, IRngMinstdVector, IRngNative

  enum, bind(c) 
    enumerator IRngNative, &    ! Built-in Fortran-90 RNG
         &     IRngMinstdVector ! Vector MINSTD algorithm
  end enum
  
  ! Maximum number of random numbers that can be computed in one call
  ! - this can be increased
  integer(kind=jpim), parameter :: NMaxStreams = 512
  
  ! A requirement of the generator is that the operation mod(A*X,M) is
  ! performed with no loss of precision, so type used for A and X must
  ! be able to hold the largest possible value of A*X without
  ! overflowing, going negative or losing precision. The largest
  ! possible value is 48271*2147483647 = 103661183124337. This number
  ! can be held in either a double-precision real number, or an 8-byte
  ! integer. Either may be used, but on some hardwares it has been
  ! found that operations on double-precision reals are faster. Select
  ! which you prefer by defining USE_REAL_RNG_STATE for double
  ! precision, or undefining it for an 8-byte integer.
#define USE_REAL_RNG_STATE 1

  ! Define RNG_STATE_TYPE based on USE_REAL_RNG_STATE, where jprd
  ! refers to a double-precision number regardless of the working
  ! precision described by jprb, while jpib describes an 8-byte
  ! integer
#ifdef USE_REAL_RNG_STATE
#define RNG_STATE_TYPE real(kind=jprd)
#else
#define RNG_STATE_TYPE integer(kind=jpib)
#endif

  ! The constants used in the main random number generator
  RNG_STATE_TYPE , parameter :: IMinstdA  = 48271
  RNG_STATE_TYPE , parameter :: IMinstdM  = 2147483647

  ! An alternative value of A that can be used to initialize the
  ! members of the state from a single seed
  RNG_STATE_TYPE , parameter :: IMinstdA0 = 16807
  
  ! Scaling to convert the state to a uniform deviate in the range 0
  ! to 1 in working precision
  real(kind=jprb), parameter :: IMinstdScale = 1.0_jprb / real(IMinstdM,jprb)

  !---------------------------------------------------------------------
  ! A random number generator type: after being initialized with a
  ! seed, type and optionally a number of vector streams, subsequent
  ! calls to "uniform_distribution" are used to fill 1D or 2D arrays
  ! with random numbers in a way that ought to be fast.
  type rng_type

    integer(kind=jpim) :: itype = IRngNative
    RNG_STATE_TYPE     :: istate(NMaxStreams)
    integer(kind=jpim) :: nmaxstreams = NMaxStreams
    integer(kind=jpim) :: iseed

  contains
    procedure :: initialize
    procedure :: uniform_distribution_1d, &
         &       uniform_distribution_2d, &
         &       uniform_distribution_2d_masked
    generic   :: uniform_distribution => uniform_distribution_1d, &
         &                               uniform_distribution_2d, &
         &                               uniform_distribution_2d_masked

  end type rng_type

contains

  !---------------------------------------------------------------------
  ! Initialize a random number generator, where "itype" may be either
  ! IRngNative, indicating to use Fortran's native random_number
  ! subroutine, or IRngMinstdVector, indicating to use the MINSTD
  ! linear congruential generator (LCG).  In the latter case
  ! "nmaxstreams" should be provided indicating that random numbers
  ! will be requested in blocks of this length. The generator is
  ! seeded with "iseed".
  subroutine initialize(this, itype, iseed, nmaxstreams)

    class(rng_type), intent(inout) :: this
    integer(kind=jpim), intent(in), optional :: itype
    integer(kind=jpim), intent(in), optional :: iseed
    integer(kind=jpim), intent(in), optional :: nmaxstreams

    integer, allocatable :: iseednative(:)
    integer :: nseed, jseed, jstr
    real(jprd) :: rseed ! Note this must be in double precision

    if (present(itype)) then
      this%itype = itype
    else
      this%itype = IRngNative
    end if
    
    if (present(iseed)) then
      this%iseed = iseed
    else
      this%iseed = 1
    end if

    if (present(nmaxstreams)) then
      this%nmaxstreams = nmaxstreams
    else
      this%nmaxstreams = NMaxStreams
    end if
    
    if (this%itype == IRngMinstdVector) then
      ! ! OPTION 1: Use the C++ minstd_rand0 algorithm to populate the
      ! ! state: this loop is not vectorizable because the state in
      ! ! one stream depends on the one in the previous stream.
      ! this%istate(1) = this%iseed
      ! do jseed = 2,this%nmaxstreams
      !   this%istate(jseed) = mod(IMinstdA0 * this%istate(jseed-1), IMinstdM)
      ! end do

      ! OPTION 2: Use a modified (and vectorized) C++ minstd_rand0 algorithm to
      ! populate the state
      rseed = real(abs(this%iseed),jprd)
      do jstr = 1,this%nmaxstreams
        ! Note that nint returns an integer of type jpib (8-byte)
        ! which may be converted to double if that is the type of
        ! istate
        this%istate(jstr) = nint(mod(rseed*jstr*(1.0_jprd-0.05_jprd*jstr &
             &      +0.005_jprd*jstr**2)*IMinstdA0, real(IMinstdM,jprd)),kind=jpib)
      end do

      ! One warmup of the C++ minstd_rand algorithm
      do jstr = 1,this%nmaxstreams
        this%istate(jstr) = mod(IMinstdA * this%istate(jstr), IMinstdM)
      end do
      
    else
      ! Native generator by default
      call random_seed(size=nseed)
      allocate(iseednative(nseed))
      do jseed = 1,nseed
        iseednative(jseed) = this%iseed + jseed - 1
      end do
      call random_seed(put=iseednative)
      deallocate(iseednative)
    end if

  end subroutine initialize

  !---------------------------------------------------------------------
  ! Populate vector "randnum" with pseudo-random numbers; if rannum is
  ! of length greater than nmaxstreams (specified when the generator
  ! was initialized) then only the first nmaxstreams elements will be
  ! assigned.
  subroutine uniform_distribution_1d(this, randnum)

    class(rng_type), intent(inout) :: this
    real(kind=jprb), intent(out)   :: randnum(:)

    integer :: imax, i

    if (this%itype == IRngMinstdVector) then
      
      imax = min(this%nmaxstreams, size(randnum))

      ! C++ minstd_rand algorithm
      do i = 1, imax
        ! The following calculation is computed entirely with 8-byte
        ! numbers (whether real or integer)
        this%istate(i) = mod(IMinstdA * this%istate(i), IMinstdM)
        ! Scale the current state to a number in working precision
        ! (jprb) between 0 and 1
        randnum(i) = IMinstdScale * this%istate(i)
      end do

    else

      call random_number(randnum)

    end if

  end subroutine uniform_distribution_1d


  !---------------------------------------------------------------------
  ! Populate matrix "randnum" with pseudo-random numbers; if the inner
  ! dimension of rannum is of length greater than nmaxstreams
  ! (specified when the generator was initialized) then only the first
  ! nmaxstreams elements along this dimension will be assigned.
  subroutine uniform_distribution_2d(this, randnum)

    class(rng_type), intent(inout) :: this
    real(kind=jprb), intent(out)   :: randnum(:,:)

    integer :: imax, jblock, i

    if (this%itype == IRngMinstdVector) then
      
      imax = min(this%nmaxstreams, size(randnum,1))

      ! C++ minstd_ran algorithm
      do jblock = 1,size(randnum,2)
        ! These lines should be vectorizable
        do i = 1, imax
          this%istate(i) = mod(IMinstdA * this%istate(i), IMinstdM)
          randnum(i,jblock) = IMinstdScale * this%istate(i)
        end do
      end do

    else

      call random_number(randnum)

    end if

  end subroutine uniform_distribution_2d

  !---------------------------------------------------------------------
  ! Populate matrix "randnum" with pseudo-random numbers; if the inner
  ! dimension of rannum is of length greater than nmaxstreams
  ! (specified when the generator was initialized) then only the first
  ! nmaxstreams elements along this dimension will be assigned. This
  ! version only operates on outer dimensions for which "mask" is true.
  subroutine uniform_distribution_2d_masked(this, randnum, mask)

    class(rng_type), intent(inout) :: this
    real(kind=jprb), intent(inout) :: randnum(:,:)
    logical,         intent(in)    :: mask(:)

    integer :: imax, jblock, i

    if (this%itype == IRngMinstdVector) then
      
      imax = min(this%nmaxstreams, size(randnum,1))

      ! C++ minstd_ran algorithm
      do jblock = 1,size(randnum,2)
        if (mask(jblock)) then
          ! These lines should be vectorizable
          do i = 1, imax
            this%istate(i) = mod(IMinstdA * this%istate(i), IMinstdM)
            randnum(i,jblock) = IMinstdScale * this%istate(i)
          end do
        end if
      end do

    else

      do jblock = 1,size(randnum,2)
        if (mask(jblock)) then
          call random_number(randnum(:,jblock))
        end if
      end do

    end if

  end subroutine uniform_distribution_2d_masked


end module radiation_random_numbers

