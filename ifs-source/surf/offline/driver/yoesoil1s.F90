MODULE YOESOIL1S
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

!       ----------------------------------------------------------------
!*     ** *YOESOIL* SOIL PARAMETERS USED ELSEWHERE
!       ----------------------------------------------------------------

REAL(KIND=JPRB) :: RWSAT
REAL(KIND=JPRB) :: RWCAP
REAL(KIND=JPRB) :: RWPWP
REAL(KIND=JPRB),ALLOCATABLE :: RWSATM(:)
REAL(KIND=JPRB),ALLOCATABLE :: RWCAPM(:)
REAL(KIND=JPRB),ALLOCATABLE :: RWPWPM(:)
REAL(KIND=JPRB),ALLOCATABLE :: RDAW(:)

!*     *YOESOIL* CONTAINS SOIL PARAMETERS
!     USED IN *VDF...* AND *SRF...*.


!     *RWSAT*     REAL     *SOIL WATER CONTENT AT SATURATION
!     *RWCAP*     REAL     *SOIL WATER CONTENT AT FIELD CAPACITY
!     *RWPWP*     REAL     *SOIL WATER CONTENT AT WILTING POINT
!     *RDAW*      REAL     *ARRAY OF LAYER THICKNESSES FOR MOISTURE

!     ------------------------------------------------------------------
END MODULE YOESOIL1S
