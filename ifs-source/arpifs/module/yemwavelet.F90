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

MODULE YEMWAVELET


!*  Module for ALADIN Wavelet Jb 
! TYPE_LAMWAVELETJB_CONFIG    contains basic configuration
!    %NWAVITER     = Number of iterations in wavelet transform
!    %NWAVSCALES   = Number of wavelet scales used in 3d-Var
!                This can be less than NWAVITER if we truncate the
!                smallest scalesi, reflected in NWAVDIMX.
!    %NGPX         = Dimension of the "full" domain (i.e. extended to a quasi-dyadic number)
!    %NGPY         =

! TYPE_LAMWAVELETJB_INFO     contains Data distribution:
!    %NWAVDIMX     = Reduction of the domain dimensions (=NGPX/2**n, n may be e.g. 1, 2)
!    %NWAVDIMY     =
!    %NWAVBLOCKS   = number of scales & orientations in a wavelet transform
!                    this is NWAVSCALES*3 + 1 (father wav)
!    %NWAVMAXX,NWAVMINX = Max and min "X" value in the vertical stritification for this PE
!    %NWAVXDIST(NPROC,2) contains the min and max X from all processors
!    %NWAVLATS : number of "X" columns locally
!    %SCALESIZEX(1:NWAVBLOCKS) = the number of points in every block
!    %SCALESIZEY
!    %NBLOCKXY     for every wavelet "position", the corresponding block number
!    %NLEVDIST     = for every PE, the corresponding min. and max. level
!    %NFLEVW       = the number of vertical levels in the local PE
!    %LWAV_SP      = Does this PE have SP in horizontal stratification (currently
!                    this is equivalent to MYPROC==NPROC)
!    %NWAV_SP      = the processor that holds SP (should be =NPROC)
!    %NWAVXDIST    = for every PE, the corresponding min. and max. "X"

! Variables for storing Jb parameters:
! ------------------------------------
! WAV_B_EIGVAL    : local eigenvalues of horizontal matrices 
! WAV_B_EIGVEC    : eigenvectors (=projection matrices)
! WAV_B_SIGMAB    : Local standard deviations in wavlet co-ordinates

! Arrays for control vector manipulations
! WAV_CV_VER      : The control vector in vertical stratification
!                   Every PE contains all vertical levels, but not all X
! WAV_CV_HOR      : Horizontal stratification. Every PE contains full horizontal              
!                   levels, but not necessarily all levels.
! ------------------------------------------------------------------
!     Author :
!     --------
!     Alex Deckmyn 
!     Modifications.
!     --------------
!      29-Sep-2008     ORIGINAL
!      18-May-2011     OOPSification
!     M. Fisher   7-March-2012 DEALLOCATES moved out of DEALGES
!-------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

PRIVATE
PUBLIC :: TYPE_LAMWAVELETJB_CONFIG, TYPE_LAMWAVELETJB_INFO, &
        & TYPE_LAMWAVELET_B_EIGVAL_STRUCT, TYPE_LAMWAVELET_B_EIGVEC_STRUCT, &
        & TYPE_LAMWAVELET_B_SIGMAB_STRUCT, TYPE_LAMWAVELET_VCOR_STRUCT, &
        & TYPE_LAMWAVELET_CV_HOR_STRUCT, TYPE_LAMWAVELET_CV_VER_STRUCT, &
        & TYPE_LAMWAVELET_GRID_STRUCT

! ------------------------------------------------------------------

TYPE TYPE_LAMWAVELETJB_CONFIG
  INTEGER(KIND=JPIM)   :: NPERIODIC
  INTEGER(KIND=JPIM)   :: NWAVITER
  INTEGER(KIND=JPIM)   :: NWAVSCALES
  INTEGER(KIND=JPIM)   :: NGPX
  INTEGER(KIND=JPIM)   :: NGPY
END TYPE TYPE_LAMWAVELETJB_CONFIG

! ------------------------------------------------------------------

TYPE TYPE_LAMWAVELETJB_INFO
  INTEGER(KIND=JPIM)   :: NWAVDIMX
  INTEGER(KIND=JPIM)   :: NWAVDIMY
  INTEGER(KIND=JPIM)   :: NWAVBLOCKS
  INTEGER(KIND=JPIM)   :: NWAVMINX
  INTEGER(KIND=JPIM)   :: NWAVMAXX
  INTEGER(KIND=JPIM)   :: NWAVLATS
  INTEGER(KIND=JPIM)   :: NFLEVW
  INTEGER(KIND=JPIM)   :: NMINLEV
  INTEGER(KIND=JPIM)   :: NMAXLEV
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:,:) :: NLEVDIST
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:,:) :: NWAVXDIST
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:)   :: SCALESIZEX
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:)   :: SCALESIZEY
  INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:,:) :: NBLOCKXY
  LOGICAL              :: LWAV_SP
  INTEGER(KIND=JPIM)   :: NWAV_SP 
END TYPE TYPE_LAMWAVELETJB_INFO

! ------------------------------------------------------------------

! Structures for B matrix

TYPE TYPE_LAMWAVELET_B_EIGVAL_STRUCT
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:,:) :: DATA_3D
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:,:) :: DATA_3DSP
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: VOR
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: DIV
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: HUM
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: TMPSP
END TYPE TYPE_LAMWAVELET_B_EIGVAL_STRUCT

! ------------------------------------------------------------------

TYPE TYPE_LAMWAVELET_B_EIGVEC_STRUCT
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:,:) :: DATA_3D
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:,:) :: DATA_3DSP
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: VOR
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: DIV
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: HUM
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: TMPSP
END TYPE TYPE_LAMWAVELET_B_EIGVEC_STRUCT

! ------------------------------------------------------------------

TYPE TYPE_LAMWAVELET_B_SIGMAB_STRUCT
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:,:) :: DATA_3D
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: DATA_2D
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: VOR
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: DIV
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: TMP
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: HUM
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:)     :: SP
END TYPE TYPE_LAMWAVELET_B_SIGMAB_STRUCT

! ------------------------------------------------------------------

TYPE TYPE_LAMWAVELET_VCOR_STRUCT
  TYPE(TYPE_LAMWAVELET_B_EIGVAL_STRUCT)   :: WAV_B_EIGVAL
  TYPE(TYPE_LAMWAVELET_B_EIGVEC_STRUCT)   :: WAV_B_EIGVEC
  TYPE(TYPE_LAMWAVELET_B_SIGMAB_STRUCT)   :: WAV_B_SIGMAB
END TYPE TYPE_LAMWAVELET_VCOR_STRUCT

! ------------------------------------------------------------------

! The Control vector for wavelet part

TYPE TYPE_LAMWAVELET_CV_HOR_STRUCT
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:,:) :: DATA_3D
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: DATA_2D
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: VOR
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: DIV
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: TMP
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: HUM
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:)     :: SP
END TYPE TYPE_LAMWAVELET_CV_HOR_STRUCT

! ------------------------------------------------------------------

TYPE TYPE_LAMWAVELET_CV_VER_STRUCT
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:,:) :: DATA_3D
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:,:) :: DATA_3DSP
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: VOR
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: DIV
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: HUM
  REAL(KIND=JPRB), ALLOCATABLE, DIMENSION (:,:,:,:)   :: TMPSP
END TYPE TYPE_LAMWAVELET_CV_VER_STRUCT

! ------------------------------------------------------------------

TYPE TYPE_LAMWAVELET_GRID_STRUCT
  TYPE(TYPE_LAMWAVELET_CV_HOR_STRUCT) :: WAV_CV_HOR
  TYPE(TYPE_LAMWAVELET_CV_VER_STRUCT) :: WAV_CV_VER
END TYPE TYPE_LAMWAVELET_GRID_STRUCT

!-------------------------------------------------------------------

END MODULE YEMWAVELET
