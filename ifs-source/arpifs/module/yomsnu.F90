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

MODULE YOMSNU

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

! -----------------------------------------------------------------------------

!     SPECTRAL ARRAYS FOR NUDGING

! XPNUDG: temporal coefficient.
! TNUDTE: coefficient for temperature.
! TNUDSH: coefficient for humidity.
! TNUDDI: coefficient for divergence.
! TNUDVO: coefficient for vorticity.
! TNUDSV: coefficient for extra-GFL variables (old passive scalars).
! TNUDLP: coefficient for log(prehyds).
! XWNUDG: time variable weight for nudging

REAL(KIND=JPRB),ALLOCATABLE:: XPNUDG(:)
REAL(KIND=JPRB),ALLOCATABLE:: TNUDTE(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: TNUDSH(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: TNUDDI(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: TNUDVO(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: TNUDSV(:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: TNUDLP(:,:)
REAL(KIND=JPRB) :: XWNUDG

! -----------------------------------------------------------------------------

END MODULE YOMSNU
