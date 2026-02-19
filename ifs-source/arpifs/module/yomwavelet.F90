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

MODULE YOMWAVELET

!     ------------------------------------------------------------------

!*     Wavelet transform

!         LJBWAVELET        =  Global switch to turn on wavelet Jb
!         N_WAVELET_CUTOFFS =  Largest wavenumber for each scale
!         N_WAVELET_SCALES  =  Number of scales in the wavelet transform
!         WAVELET_FILTER    =  Lowpass filter coefficients
!         NSMIN_WAVELET     =  Minimum spherical transform resolution
!                              (The transform package insists on having
!                               more wavenumbers than processors.)
!         JPNSCALES      =  Max. number of allowed spactral bands in wavelet decomp. (default is 20)
!         SKYLINE_TOL    =  To reduce computational cost, small elements of the symmetric
!                           square roots of the vertical covariance matrices are replaced
!                           by zero. This is done by constructing a "skyline" for each matrix.
!                           For each row, M*_SKYLINE_LEFT and M*_SKYLINE_RIGHT give the first
!                           and last elements of each row that must be used, such that the
!                           absolute values of all the unused elements are less than
!                           SKYLINE_TOL times the diagonal element for the row.

!         WAVELET_VCORS  =  vertical covariance matrices

!         LJBWSTATS      =  .T. => calculate and output covariance
!                           matrices, rather than reading them. 
!         LHYBRID_JB: .T. => hybrid computation of B from a pre-existing wavelet file and a set of 
!                        new ensemble background forecasts. Relative weight of the two is controlled by ALPHA_HYBRID_JB
!         ALPHA_HYBRID_JB: Relative weight of new ensemble perturbations to input wavelet B in the
!                         computation of a new B [0.,1.]
!         L_TAPER_ALPHA_HYBRID_JB: Reduce relative weight gradually towards zero between 1hPa and 0.2hPa
!         LONLINEWEIGHTED_HYBRID_JB: - Change weight (=ALPHA) with number of samples/matrix per state (=N),
!                                        ALPHA0 = ALPHA_HYBRID_JB
!                                        ALPHA  = N*ALPHA0/((1-ALPHA0)+N*ALPHA0)
!                                        Example: For ALPHA0=0.3, ALPHA=0.3<T63 increasing to ALPHA=0.94 at T399.
!                                    - Use only samples of the day in powerspectrum/variance (=>ALPHA=1)
!         L_IN_WAVELET_STATS_CALC = .true. when under SUJBWAVSTATS, .false. otherwise
!         L_CALCULATE_JBW_INV = active only when SUJBWAVSTATS = .T. ; if .T. both JB and JB-1
!                           are computed and written out in wavelet_out.cv, otherwise only JB
!         LJBWFDMEM      =  .T. => read forecast differences once and store
!                           in memory. .F. => reduce memory at the expense
!                           of more I/O.
!         LJBWHYB        =  .T. => calculate an hybrid Jb term
!         LDISTREAD      =  .T. => Distribute reading of background fields over processors
!         N_FirstScale   =1  => First scale for diagnostics of wavelet variance
!         N_LastScale    =-1 => Last scale for diagnostics of wavelet variance
!         CINBGSTATES    =  file name prefix for input background states
!         N_BGDATES      =  number of states per member for the input analysis ensemble
!         N_BGSTEPBY     =  number of hours between successive BG states in JB computation
!         N_BGMEMBERS    =  number of members in the input analysis ensemble

!      TYPE_WAVELETJB_VCOR_STRUCT        = derived type for correlation matrices
!      ---------"---------%NLATSG = number of latitudes (all processors)
!      ---------"---------%NLATS  = number of latitudes (this processor)
!      ---------"---------%NLONSG = longitudes for each row (all PEs)
!      ---------"---------%NLONS  = longitudes for each row (this PE)
!      ---------"---------%NPTOTG = number of horizontal points (all PEs)
!      ---------"---------%NPTOT  = number of horizontal points (this PE)
!      ---------"---------%NLEVS  = first dimension of a matrix
!      ---------"---------%NLON_OFFSET offset to this PE's zeroth longitude
!      ---------"---------%NLAT_OFFSET = offset to this PE's zeroth latitude
!      ---------"---------%LON0G  = first longitude of each row (all PEs)
!      ---------"---------%DLON   = longitudinal gridpacing each row (this PE)
!      ---------"---------%LATSG  = latitudes of rows (all processors)
!      ---------"---------%LPELIST= flags showing which PEs update which
!                                   matrices
!      ---------"---------%MASTERPE = master PE for each matrix.
!      ---------"---------%INDEXG2L = converts global matrix number to
!                                     local matrix number
!      ---------"---------%MATS   = square roots of correlation matrices
!      ---------"---------%MATINVS  = inverse square roots of matrices
!      ---------"---------%NSAMPS  = sample size for covariance calc'n
!      ---------"---------%SAMPAV  = sample mean for covariance calc'n
!      ---------"---------%EVALS   = eigenvalues of the covariance matrices
!      ---------"---------%MATS_SKYLINE_LEFT  = matrix skyline for MATS (see SKYLINE_TOL above)
!      ---------"---------%MATS_SKYLINE_RIGHT = matrix skyline for MATS
!      ---------"---------%MATINVS_SKYLINE_LEFT  = matrix skyline for MATINVS
!      ---------"---------%MATINVS_SKYLINE_RIGHT = matrix skyline for MATINVS

!      TYPE_WAVELETJB_GRID_STRUCT        = derived type for grid definitions
!      ---------"---------%NDGLG       = Number of Gaussian latitudes
!      ---------"---------%NLOENG      = Number of points in each row
!      ---------"---------%NRESOL_ID   = Resolution id for transform package
!      ---------"---------%NGPTOT      = Number of gridpoints on this PE
!      ---------"---------%NGPTOTG     = Global number of gridpoints
!      ---------"---------%NGPTOTL     = number of gridpoints on this and other PEs
!      ---------"---------%NGPTOTMX    = max(NGPTOT) over all PEs
!      ---------"---------%NSMAX       = Spectral truncation
!      ---------"---------%NSPEC2      = Number of spectral coeffs on this PE
!      ---------"---------%GELAM       = Longitudes of gridpoints on this PE
!      ---------"---------%GELAT       = Latitudes  of gridpoints on this PE
!      ---------"---------%GW          = Gaussian quadrature weights
!      ---------"---------%NSTA        = KSTA      array returned by TRANS_INQ
!      ---------"---------%NONL        = KONL      array returned by TRANS_INQ
!      ---------"---------%NPTRFRSTLAT = KPTRFRSTLAT returned by TRANS_INQ
!      ---------"---------%NFRSTLAT    = KFRSTLAT  array returned by TRANS_INQ
!      ---------"---------%NLSTLAT     = KLSTLAT   array returned by TRANS_INQ
!      ---------"---------%NPTRFLOFF   = KPTRFLOFF array returned by TRANS_INQ
!      ---------"---------%NUMP        = KUMP      Number of m's on this PE
!      ---------"---------%MYMS        = MYMS      List of m's on this PE
!      ---------"---------%NASM0       = KASM0     array returned by TRANS_INQ
!      ---------"---------%NBLOCKSL    = Number of NPROMA blocks on this PE
!      ---------"---------%NBLOCKSG    = Number of NPROMA blocks on all PEs

!  Modified:
!    M.Fisher   17-Mar-2004 : Wavelet Jb control vector in grid space(s)
!    H.Varella  15-Nov-2011 : Option LJBWFDMEM included
!    h.Varella  04-Jan-2012 : Option LJBWHYB for an hybrid Jb term
!    M. Fisher   7-Mar-2012 : Moved DEALLOCATES out of DEALGES into module
!    M.Bonavita 11-May-2012 : L_CALCULATE_JBW_INV, LJBWFDMEM
!    E.Holm     20-Apr-2013 : LEDA_READ_AV
!    V. Chabot 30-Mar-2016  : Possibility to specify some scales for wavelet diagnostic.
!    E.Holm     10-Oct-2021 : Remove LEDA_READ_AV
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

PRIVATE
PUBLIC :: TYPE_WAVELETJB_CONFIG, TYPE_WAVELETJB_DATA, &
        & TYPE_WAVELETJB_VCOR_STRUCT, TYPE_WAVELETJB_GRID_STRUCT, &
        & JPNSCALES

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER   :: JPNSCALES=20

TYPE TYPE_WAVELETJB_CONFIG
  LOGICAL              :: LJBWAVELET
  LOGICAL              :: LJBWSTATS
  LOGICAL              :: LHYBRID_JB
  LOGICAL              :: L_TAPER_ALPHA_HYBRID_JB
  LOGICAL              :: LONLINEWEIGHTED_HYBRID_JB
  LOGICAL              :: LJBWFDMEM
  LOGICAL              :: L_CALCULATE_JBW_INV
  LOGICAL              :: LJBWHYB
  LOGICAL              :: LDISTREAD
  INTEGER(KIND=JPIM)   :: N_FirstScale 
  INTEGER(KIND=JPIM)   :: N_LastScale 
  INTEGER(KIND=JPIM)   :: N_WAVELET_SCALES
  INTEGER(KIND=JPIM)   :: N_BGMEMBERS
  INTEGER(KIND=JPIM)   :: N_BGDATES
  INTEGER(KIND=JPIM)   :: N_BGSTEPBY
  INTEGER(KIND=JPIM)   :: NSMIN_WAVELET  
  INTEGER(KIND=JPIM)   :: JB_WAVELET_SCALES(JPNSCALES)
  REAL(KIND=JPRB)      :: SKYLINE_TOL
  REAL(KIND=JPRB)      :: ALPHA_HYBRID_JB
  CHARACTER (LEN = 16) ::  CINBGSTATES
END TYPE TYPE_WAVELETJB_CONFIG

!     ------------------------------------------------------------------

TYPE TYPE_WAVELETJB_DATA
  LOGICAL                     :: L_IN_WAVELET_STATS_CALC
  INTEGER(KIND=JPIM), ALLOCATABLE :: N_WAVELET_CUTOFFS(:)
  REAL(KIND=JPRB),    ALLOCATABLE :: WAVELET_FILTERS(:,:)
END TYPE TYPE_WAVELETJB_DATA

!     ------------------------------------------------------------------

TYPE TYPE_WAVELETJB_VCOR_STRUCT
  INTEGER(KIND=JPIM) :: NLATSG
  INTEGER(KIND=JPIM) :: NLATS
  INTEGER(KIND=JPIM) :: NPTOTG
  INTEGER(KIND=JPIM) :: NPTOT
  INTEGER(KIND=JPIM) :: NLEVS
  INTEGER(KIND=JPIM) :: NLAT_OFFSET  
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:)     :: NLONSG
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:)     :: NLONS
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:)     :: NSAMPS
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:)     :: INDEXG2L
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:)     :: MASTERPE
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:)     :: NLON_OFFSET
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:,:)   :: MATS_SKYLINE_LEFT
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:,:)   :: MATS_SKYLINE_RIGHT
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:,:)   :: MATINVS_SKYLINE_LEFT
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:,:)   :: MATINVS_SKYLINE_RIGHT
  LOGICAL,            ALLOCATABLE, DIMENSION(:,:)   :: LPELIST
  REAL(KIND=JPRB),    ALLOCATABLE, DIMENSION(:)     :: LON0G
  REAL(KIND=JPRB),    ALLOCATABLE, DIMENSION(:)     :: DLON
  REAL(KIND=JPRB),    ALLOCATABLE, DIMENSION(:)     :: LATSG
  REAL(KIND=JPRB),    ALLOCATABLE, DIMENSION(:,:)   :: SAMPAV
  REAL(KIND=JPRB),    ALLOCATABLE, DIMENSION(:,:)   :: EVALS
  REAL(KIND=JPRB),    ALLOCATABLE, DIMENSION(:,:,:) :: MATS
  REAL(KIND=JPRB),    ALLOCATABLE, DIMENSION(:,:,:) :: MATINVS
END TYPE TYPE_WAVELETJB_VCOR_STRUCT

!     ------------------------------------------------------------------

TYPE TYPE_WAVELETJB_GRID_STRUCT
  INTEGER(KIND=JPIM)                            :: NDGLG
  INTEGER(KIND=JPIM)                            :: NRESOL_ID
  INTEGER(KIND=JPIM)                            :: NGPTOT
  INTEGER(KIND=JPIM)                            :: NGPTOTG
  INTEGER(KIND=JPIM)                            :: NGPTOTMX
  INTEGER(KIND=JPIM)                            :: NSPEC2
  INTEGER(KIND=JPIM)                            :: NSMAX
  INTEGER(KIND=JPIM)                            :: NPTRFLOFF
  INTEGER(KIND=JPIM)                            :: NUMP
  INTEGER(KIND=JPIM)                            :: NBLOCKSL
  INTEGER(KIND=JPIM)                            :: NBLOCKSG  
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:)     :: NLOENG
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:)     :: NPTRFRSTLAT
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:)     :: NFRSTLAT
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:)     :: NLSTLAT
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:)     :: MYMS
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:)     :: NASM0
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:,:)   :: NSTA
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:,:)   :: NONL
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:,:)   :: NGPTOTL
  REAL(KIND=JPRB),    ALLOCATABLE, DIMENSION(:)     :: GELAM
  REAL(KIND=JPRB),    ALLOCATABLE, DIMENSION(:)     :: GELAT
  REAL(KIND=JPRB),    ALLOCATABLE, DIMENSION(:)     :: GW
END TYPE TYPE_WAVELETJB_GRID_STRUCT

!     ------------------------------------------------------------------

END MODULE YOMWAVELET
