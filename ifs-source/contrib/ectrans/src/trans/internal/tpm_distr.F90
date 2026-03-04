! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_DISTR

! Module for distributed memory environment.

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!*    Variables describing distributed memory parallelization

INTEGER(KIND=JPIM) :: NPROC     ! Number of processors (NPRGPNS*NPRGPEW) 
INTEGER(KIND=JPIM) :: NPRGPNS   ! No. of sets in N-S direction (grid-point space)
INTEGER(KIND=JPIM) :: NPRGPEW   ! No. of sets in E-W direction (grid-point space)
INTEGER(KIND=JPIM) :: NPRTRW    ! No. of sets in wave direction (spectral space)
INTEGER(KIND=JPIM) :: NPRTRV    ! NPROC/NPRTRW
INTEGER(KIND=JPIM) :: NPRTRNS   ! No. of sets in N-S direction (Fourier space)
                                ! (always equal to NPRTRW)
LOGICAL            :: LEQ_REGIONS ! TRUE - Use new eq_regions partitioning
                                  ! FALSE- Use old NPRGPNS x NPRGPEW partitioning
INTEGER(KIND=JPIM) :: MYPROC    ! My processor number
INTEGER(KIND=JPIM) :: MYSETW    ! My set number in wave direction (spectral space) 
INTEGER(KIND=JPIM) :: MYSETV    ! My set number in field direction(S.S and F.S)
INTEGER(KIND=JPIM) :: NCOMBFLEN ! Size of communication buffer

INTEGER(KIND=JPIM) :: MTAGLETR   ! Tag
INTEGER(KIND=JPIM) :: MTAGML     ! Tag
INTEGER(KIND=JPIM) :: MTAGLG     ! Tag
INTEGER(KIND=JPIM) :: MTAGGL     ! Tag
INTEGER(KIND=JPIM) :: MTAGPART   ! Tag
INTEGER(KIND=JPIM) :: MTAGDISTSP ! Tag
INTEGER(KIND=JPIM) :: MTAGLM     ! Tag
INTEGER(KIND=JPIM) :: MTAGDISTGP ! Tag
 
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPRCIDS(:) ! Array containing the process ids

TYPE DISTR_TYPE
LOGICAL   :: LGRIDONLY          ! TRUE - only grid space structures are available
LOGICAL   :: LWEIGHTED_DISTR    ! TRUE - weighted distribution
LOGICAL   :: LSPLIT             ! TRUE - latitudes are shared between a-sets
LOGICAL   :: LCPNMONLY          ! TRUE - Compute Legendre polynomials only, not FFTs

! SPECTRAL SPACE

INTEGER(KIND=JPIM) :: NUMP      ! No. of spectral waves handled by this processor
INTEGER(KIND=JPIM) :: NSPEC     ! No. of complex spectral coefficients (on this PE)
INTEGER(KIND=JPIM) :: NSPEC2    ! 2*NSPEC
INTEGER(KIND=JPIM) :: NSPEC2MX  ! maximun NSPEC2 among all PEs
INTEGER(KIND=JPIM) :: NTPEC2    ! cf. NSPEC2 but for truncation NTMAX
INTEGER(KIND=JPIM) :: NUMTP     ! cf. NUMP but for truncation NTMAX

INTEGER(KIND=JPIM) :: NSPOLEGL  ! No. of legendre polynomials on this PE
INTEGER(KIND=JPIM) :: NLEI3D    ! (NLEI3-1)/NPRTRW+1

INTEGER(KIND=JPIM) ,ALLOCATABLE :: MYMS(:)    ! Wave numbers handled by this PE
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NUMPP(:)   ! No. of wave numbers each wave set is
                                 ! responsible for
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPOSSP(:)  ! Not needed in transform?
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPROCM(:)  ! Process that does the calc. for certain 
                                 ! wavenumber M
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NDIM0G(:)  ! Defines partitioning of global spectral
                                 ! fields among PEs

INTEGER(KIND=JPIM) ,ALLOCATABLE :: NASM0(:)  ! Address in a spectral array of (m, n=m)
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NATM0(:)  ! Same as NASM0 but for NTMAX
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NALLMS(:) ! Wave numbers for all a-set concatenated
                                ! together to give all wave numbers in a-set
                                ! order. Used when global spectral norms 
                                ! have to be gathered.
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPTRMS(:) ! Pointer to the first wave number of a given
                                ! a-set in nallms array.


! Legendre polynomials

INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLATLS(:,:) ! First latitude for which each a-set,bset calcul.
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLATLE(:,:) ! Last latitude for which each a-set,bset calcul.

INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPMT(:) ! Adress for legendre polynomial for
                              ! given M (NTMAX)
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPMS(:) ! Adress for legendre polynomial for
                              ! given M (NSMAX)
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPMG(:) ! Global version of NPMS

! FOURIER SPACE

INTEGER(KIND=JPIM) :: NDGL_FS ! Number of rows of latitudes for which this process is
                     ! performing Fourier Space calculations

INTEGER(KIND=JPIM) ,ALLOCATABLE  :: NSTAGTF(:) ! Offset for specific latitude in 
                                  ! Fourier/gridpoint buffer
INTEGER(KIND=JPIM) :: NLENGTF ! Second dimension of Fourier/gridpoint buffer 
                     ! (sum of (NLOEN+3) over local latitudes)

INTEGER(KIND=JPIM) ,ALLOCATABLE :: NULTPP(:) ! No of lats. for each wave_set  (F.S)
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPROCL(:) ! Process responsible for each lat. (F.S)
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPTRLS(:) ! Pointer to first lat. (F.S)

! NSTAGT0B to NLENGT1B: help arrays for spectral to fourier space transposition
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NSTAGT0B(:) ! Start adresses for segments within buffer
                                  ! (according to processors to whom data 
                                  ! is going to be sent) 
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NSTAGT1B(:) 
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPNTGTB0(:,:)
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPNTGTB1(:,:)
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLTSFTB(:)  
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLTSGTB(:)
INTEGER(KIND=JPIM) ,ALLOCATABLE :: MSTABF(:)
INTEGER(KIND=JPIM) :: NLENGT0B  ! dimension
INTEGER(KIND=JPIM) :: NLENGT1B  ! dimension

! GRIDPOINT SPACE

INTEGER(KIND=JPIM) :: NDGL_GP ! D%NLSTLAT(MY_REGION_NS)-D%NFRSTLOFF
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NFRSTLAT(:) ! First lat of each a-set 
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLSTLAT(:)  ! Last lat of each a-set 
INTEGER(KIND=JPIM) :: NFRSTLOFF ! Offset for first lat of own a-set 
                       ! i.e. NFRSTLOFF=NFRSTLAT(MYSETA)-1
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPTRLAT(:) ! Pointer to start of latitude 
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPTRFRSTLAT(:) ! Pointer to the first latitude of each 
                                     ! a-set in NSTA and NONL arrays
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPTRLSTLAT(:) ! Pointer to the last latitude of each
                                    ! a-set in NSTA and NONL arrays
INTEGER(KIND=JPIM) :: NPTRFLOFF ! Offset for pointer to the first latitude of own a-set
                       ! NSTA and NONL arrays, i.e. NPTRFRSTLAT(MYSETA)-1
LOGICAL   ,ALLOCATABLE :: LSPLITLAT(:) ! True if latitude is split over 2 a-sets

!  NSTA(R%NDGL+NPRGPNS-1,NPRGPEW) :  Position of first grid column
!             for the latitudes on a processor. The information is
!             available for all processors. The b-sets are distinguished
!             by the last dimension of NSTA(). The latitude band for
!             each a-set is addressed by NPTRFRSTLAT(JASET),
!             NPTRLSTLAT(JASET), and NPTRFLOFF=NPTRFRSTLAT(MYSETA) on
!             this processors a-set. Each split latitude has two entries
!             in NSTA(,:) which necessitates the rather complex
!             addressing of NSTA(,:) and the overdimensioning of NSTA by
!             NPRGPNS.
!  NONL(R%NDGL+NPRGPNS-1,NPRGPEW)  :  Number of grid columns for
!             the latitudes on a processor. Similar to NSTA() in data
!             structure.
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NSTA(:,:)
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NONL(:,:)

INTEGER(KIND=JPIM) :: NGPTOT   ! Total number of grid columns on this PE
INTEGER(KIND=JPIM) :: NGPTOTG  ! Total number of grid columns on the Globe
INTEGER(KIND=JPIM) :: NGPTOTMX ! Maximum number of grid columns on any of the PEs
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NGPTOTL(:,:) ! Number of grid columns on each PE.

REAL(KIND=JPRB) ,ALLOCATABLE :: RWEIGHT(:) ! Weight per grid-point (if weighted distribution)
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NPROCA_GP(:) ! Number of grid-points per a-set

END TYPE DISTR_TYPE

TYPE(DISTR_TYPE),ALLOCATABLE,TARGET :: DISTR_RESOL(:)
TYPE(DISTR_TYPE),POINTER     :: D

END MODULE TPM_DISTR

