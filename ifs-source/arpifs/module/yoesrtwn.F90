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

MODULE YOESRTWN

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

INTEGER(KIND=JPIM) , DIMENSION(16:29) :: NG
INTEGER(KIND=JPIM) , DIMENSION(16:29) :: NSPA
INTEGER(KIND=JPIM) , DIMENSION(16:29) :: NSPB
INTEGER(KIND=JPIM) , DIMENSION(14)    :: NMPSRTM

REAL(KIND=JPRB), DIMENSION(59) :: PREF
REAL(KIND=JPRB), DIMENSION(59) :: PREFLOG
REAL(KIND=JPRB), DIMENSION(59) :: TREF

INTEGER(KIND=JPIM), DIMENSION(224) :: NGM
INTEGER(KIND=JPIM), DIMENSION(14)  :: NGC, NGS

REAL(KIND=JPRB), DIMENSION(16)  :: WT, WTSM
REAL(KIND=JPRB), DIMENSION(224) :: RWGT


! Use for 56 g-points
!INTEGER(KIND=JPIM), DIMENSION(56) :: NGBSW, NGN
! Use for 112 g-points
!INTEGER(KIND=JPIM), DIMENSION(112) :: NGBSW, NGN
! Use for 224 g-points
!INTEGER(KIND=JPIM), DIMENSION(224) :: NGBSW, NGN

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      03-03-07
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----    : -------
!  NG     : INTEGER : Number of k-coefficients in spectral intervals
!  NSPA   : INTEGER :
!  NSPB   : INTEGER :
! NMPSRTM : INTEGER : MAPPING INDICES FOR 6-SPECTRAL INT. SURFACE ALBEDO 
! WAVENUM1: REAL    : Lower wavenumber spectral limit
! WAVENUM2: REAL    : Higher wavenumber spectral limit
! DELWAVE : REAL    : Spectral interval width
! PREF    : REAL    : Reference pressure
! PREFLOG : REAL    : Log reference pressure
! TREF    : REAL    : Reference temperature
!  NGC    : INTEGER : The number of new g-points in each band
!  NGS    : INTEGER : The cumulative sum of new g-points for each band
!  NGM    : INTEGER : The index of each new g-point relative to the
!                     original 16 g-points for each band.
!
!  WT     : REAL    : RRTM weights for 16 g-points.
!  WTSUM  : REAL    : Sum of the weights
!  RWGT   : REAL    :
! 
!  NGN    : INTEGER : The number of original g-points that are combined 
!                     to make each new g-point in each band.
!  NGBSW  : INTEGER : The band index for each new g-point.
!     -----------------------------------------------------------------
END MODULE YOESRTWN

