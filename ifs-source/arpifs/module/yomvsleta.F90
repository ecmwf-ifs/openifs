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

MODULE YOMVSLETA

USE PARKIND1 , ONLY : JPIM, JPRB
IMPLICIT NONE

SAVE

!*    * Vertical geometry-dependent intermediate quantities used to
!       compute vertical weights in the interpolator.
TYPE TVSLETA
  !   VCUICO: is used to compute denominators of weights
  !           for vertical interpolations
  !           applied to full-level variables.
  REAL(KIND=JPRB), ALLOCATABLE :: VCUICO(:,:)
  !   VCUICOH:is used to compute denominators of weights 
  !           for vertical interpolations
  !           applied to half-level variables.
  REAL(KIND=JPRB), ALLOCATABLE :: VCUICOH(:,:)
  !   VSLD  : auxiliary quantity for SLHD interpolation weights on full levels
  REAL(KIND=JPRB), ALLOCATABLE :: VSLD(:,:)
  !   VSLDH : auxiliary quantity for SLHD interpolation weights on half levels
  REAL(KIND=JPRB), ALLOCATABLE :: VSLDH(:,:)
  !   VSLDW : weights for SLHD vertical Laplacian smoother on full levels
  REAL(KIND=JPRB), ALLOCATABLE :: VSLDW(:,:,:)
  !   VSLDWH: weights for SLHD vertical Laplacian smoother on half levels
  REAL(KIND=JPRB), ALLOCATABLE :: VSLDWH(:,:,:)
  !   VSLVF : weights for SLVF vertical Laplacian smoother on full levels
  REAL(KIND=JPRB), ALLOCATABLE :: VSLVF(:,:,:)
  !   VSLVFH: weights for SLVF vertical Laplacian smoother on half levels
  REAL(KIND=JPRB), ALLOCATABLE :: VSLVFH(:,:,:)
  !   GAMMA_WENO: weights for vertical WENO interpolations determined from
  !               the polynomial identity
  REAL(KIND=JPRB), ALLOCATABLE :: GAMMA_WENO(:,:)
  !   VRDETAR: ratio (eta(lbar)-eta(lbar-1))/(eta(lbar-1)-eta(lbar-2)),
  !            is used in the interpolation routine LAITVSPCQM to
  !            ensure monotonicity and conservation properties.
  REAL(KIND=JPRB), ALLOCATABLE :: VRDETAR(:)
  !   NRLEVX: used to dimension NVAUTF and NVAUTH.
  INTEGER(KIND=JPIM)      :: NRLEVX = 1
  !   VRLEVX: REAL(NRLEVX)
  REAL(KIND=JPRB)         :: VRLEVX = 1.0_JPRB
  !   NVAUTF: NVAUTF(VRLEVX*eta) is the number of the layer (full level)
  !           immediately above "eta", and is bounded by 1 and nflevg-1.
  INTEGER(KIND=JPIM), ALLOCATABLE:: NVAUTF(:)
  !   NVAUTH: NVAUTH(VRLEVX*eta) is the number of the interlayer (half level)
  !           immediately above "eta", and is bounded by 0 and nflevg-1.
  INTEGER(KIND=JPIM), ALLOCATABLE :: NVAUTH(:)
END TYPE TVSLETA

!!TYPE(TVSLETA), POINTER :: YRVSLETA => NULL()

END MODULE YOMVSLETA
