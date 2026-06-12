MODULE YOMGF1S
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE
SAVE

!*    grid point array for large scale forcing of the 1D surface model
!          (current time step)

REAL(KIND=JPRB),ALLOCATABLE,TARGET:: GF(:,:)

! -   Splitting of the array GF

!  === Near surface atmospheric variables ===

!  Variables at current time step

! UNLEV0   -   U-wind at lowest atmospheric model level
! VNLEV0   -   V-wind at lowest atmospheric model level
! TNLEV0   -   Temperature at lowest atmospheric model level
! QNLEV0   -   Sp. humidity at lowest atmospheric model level
! CNLEV0   -   tracer at lowest atmospheric model level
! PNLP0   -   Surface pressure

!  Variables at next time step

! UNLEV1   -   U-wind at lowest atmospheric model level
! VNLEV1   -   V-wind at lowest atmospheric model level
! TNLEV1   -   Temperature at lowest atmospheric model level
! QNLEV1   -   Sp. humidity at lowest atmospheric model level
! CNLEV1   -   tracer at lowest atmospheric model level
! PNLP1   -   Surface pressure 

!  Variables at current time step (Forcing surface fluxes)

! FSSRD    -   Surface solar radiation downwards
! FSTRD    -   Surface thermal radiation downwards
! FLSRF    -   Large scale rainfall
! FCRF     -   Convective rainfall
! FLSSF    -   Large scale snowfall
! FCSF     -   Convective snowfall 

!  === Constants ===

! RALT     -   Observation height  (metre)
! RZUV     -   Wind observation height (m)


REAL(KIND=JPRB),POINTER:: UNLEV0(:)
REAL(KIND=JPRB),POINTER:: VNLEV0(:)
REAL(KIND=JPRB),POINTER:: TNLEV0(:)
REAL(KIND=JPRB),POINTER:: QNLEV0(:)
REAL(KIND=JPRB),POINTER:: CNLEV0(:,:)
REAL(KIND=JPRB),POINTER:: PNLP0(:)
REAL(KIND=JPRB),POINTER:: UNLEV1(:)
REAL(KIND=JPRB),POINTER:: VNLEV1(:)
REAL(KIND=JPRB),POINTER:: TNLEV1(:)
REAL(KIND=JPRB),POINTER:: QNLEV1(:)
REAL(KIND=JPRB),POINTER:: CNLEV1(:,:)
REAL(KIND=JPRB),POINTER:: PNLP1(:)
REAL(KIND=JPRB),POINTER:: FSSRD(:)
REAL(KIND=JPRB),POINTER:: FSTRD(:)
REAL(KIND=JPRB),POINTER:: FLSRF(:)
REAL(KIND=JPRB),POINTER:: FCRF(:)
REAL(KIND=JPRB),POINTER:: FLSSF(:)
REAL(KIND=JPRB),POINTER:: FCSF(:)

REAL(KIND=JPRB) :: RALT , RZUV

END MODULE YOMGF1S
