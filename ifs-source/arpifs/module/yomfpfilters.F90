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

MODULE YOMFPFILTERS

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

PRIVATE
PUBLIC TFPFILTERS

!     ------------------------------------------------------------------

!*    Full-POS spectral filter

! LFPFIL   : .TRUE. if the filter for dilated fields is active (arpege only) for
!            the corresponding (sub)domain
! RFPFIL   : value of the filter for each zonal wavenumber and each subdomain

! LFPMAT   : .TRUE. if filtering matrix in the homogenous resolution space must be computed
! RFPMAT   : Filtering matrix in the homogenous resolution space (ARPEGE only)

TYPE TFPFILTERS

LOGICAL, ALLOCATABLE :: LFPFIL(:)
REAL(KIND=JPRB), ALLOCATABLE :: RFPFIL(:,:)

LOGICAL :: LFPMAT
REAL(KIND=JPRB), ALLOCATABLE :: RFPMAT(:,:)

END TYPE TFPFILTERS

!     ------------------------------------------------------------------
END MODULE YOMFPFILTERS
