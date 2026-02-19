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

MODULE YOMFPGEO

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

! target grid or interpolation grid parameters
! ============================================

! NFPRGPG     : global number of points in output
! NFPRGPL     : local number of output points
! NFPRGPLX    : maximum number of output points where among all MPI tasks
! NFPROMA     : cache-blocking factor
! NFPBLOCS    : Number of NFPROMA-sized blocks
! NFPEND      : actual size of each NFPROMA-sized block
! NFPRGPNUM   : number of grindpoints for each MPI task
! NFPRGPIND   : if "iloc" is the DM-local index in a NFPRGPL-array and "jroc" the processor,
!                 NFPRGPIND(iloc,jroc) gives the DM-global index "iglo" in a NFPRGPG-array.

! -- Geographical fields actually used on target geometry only :
! RFPGM   : map factor
! RFPNORX : compass, first momentum
! RFPNORY : compass, second momentum
! RFPMASK : binary mask for missing geographical points (1. = present ; 0. = missing)

! -- Geographical fields actually used on interpolation geometry only :
! RFPLA   : latitudes of the output points
! RFPLO   : longitudes of the output points
! RFPGMS  : output geographic resolution = "rfpms"/"rfpgm"
! NFPNUMD : subdomain index of the output points


TYPE TFPGEO

INTEGER(KIND=JPIM) :: NFPRGPG
INTEGER(KIND=JPIM) :: NFPRGPL
INTEGER(KIND=JPIM) :: NFPRGPLX
INTEGER(KIND=JPIM) :: NFPROMA
INTEGER(KIND=JPIM) :: NFPBLOCS
INTEGER(KIND=JPIM),ALLOCATABLE :: NFPEND(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: NFPRGPNUM(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: NFPRGPIND(:,:)

REAL(KIND=JPRB),ALLOCATABLE :: RFPGM(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: RFPNORX(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: RFPNORY(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: RFPMASK(:,:)

REAL(KIND=JPRB),ALLOCATABLE :: RFPLA(:)
REAL(KIND=JPRB),ALLOCATABLE :: RFPLO(:)
REAL(KIND=JPRB),ALLOCATABLE :: RFPGMS(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: NFPNUMD(:)

END TYPE TFPGEO

!     ------------------------------------------------------------------
END MODULE YOMFPGEO
